export kepler_propagate, propagate_segment, propagate_leg
export compute_mismatch, compute_total_Δv
export compute_sundman_segments

#=
Propagation functions for Sims-Flanagan transcription.

The Sims-Flanagan method uses Kepler propagation between impulsive ΔV maneuvers
applied at segment midpoints. This file provides the core propagation routines
using universal variables (Stumpff functions) for robust propagation.

Supports multiple propulsion types:
- Spacecraft: Constant thrust (chemical/electric)
- SEPSpacecraft: Solar Electric Propulsion (thrust varies with solar distance)
- SolarSail: Solar radiation pressure (no propellant consumption)
=#

# =============================================================================
# Stumpff Functions
# =============================================================================

"""
    stumpff_c2(z)

Stumpff function c2(z) = (1 - cos(√z))/z for z > 0
                       = (cosh(√-z) - 1)/(-z) for z < 0
                       = 1/2 for z ≈ 0
"""
function stumpff_c2(z::T) where {T<:Number}
    if z > T(1e-6)
        sqrt_z = sqrt(z)
        return (one(T) - cos(sqrt_z)) / z
    elseif z < T(-1e-6)
        sqrt_neg_z = sqrt(-z)
        return (cosh(sqrt_neg_z) - one(T)) / (-z)
    else
        # Taylor series for small |z|: c2 = 1/2 - z/24 + z²/720 - ...
        return one(T)/2 - z/24 + z^2/720 - z^3/40320
    end
end

"""
    stumpff_c3(z)

Stumpff function c3(z) = (√z - sin(√z))/(√z)³ for z > 0
                       = (sinh(√-z) - √-z)/(√-z)³ for z < 0
                       = 1/6 for z ≈ 0
"""
function stumpff_c3(z::T) where {T<:Number}
    if z > T(1e-6)
        sqrt_z = sqrt(z)
        return (sqrt_z - sin(sqrt_z)) / (sqrt_z^3)
    elseif z < T(-1e-6)
        sqrt_neg_z = sqrt(-z)
        return (sinh(sqrt_neg_z) - sqrt_neg_z) / (sqrt_neg_z^3)
    else
        # Taylor series for small |z|: c3 = 1/6 - z/120 + z²/5040 - ...
        return one(T)/6 - z/120 + z^2/5040 - z^3/362880
    end
end

# =============================================================================
# Sundman Transformation
# =============================================================================

"""
    compute_sundman_segments(r_start, r_end, tof_leg, n_segments; c=1.0)

Compute segment durations for a single leg using the Sundman transformation.

The Sundman transformation uses: dt/ds = r^c

This causes segments near the central body to have shorter durations (more 
resolution where gravity is stronger and thrust is less effective for SEP).

For a trajectory from r_start to r_end split into n_segments equal steps in the 
Sundman variable s, each segment i has time duration Δt_i ∝ r_i^c where r_i is 
the radius at the midpoint of segment i.

# Arguments
- `r_start::AbstractVector`: Start position of leg [km]
- `r_end::AbstractVector`: End position of leg [km]  
- `tof_leg::Number`: Time of flight for this leg [s]
- `n_segments::Int`: Number of segments in this leg

# Keyword Arguments
- `c::Number`: Sundman exponent (default: 1.0)
  - c = 0: Equal time segments (no Sundman transformation)
  - c = 1: Classic Sundman (dt ∝ r)
  - c = 1.5: Related to eccentric anomaly

# Returns
- `Δt_segments::Vector`: Duration of each segment [s]

# Reference
Based on PyKEP's automated mesh adaptation approach.
"""
function compute_sundman_segments(
    r_start::AbstractVector,
    r_end::AbstractVector,
    tof_leg::Number,
    n_segments::Int;
    c::Number = 1.0,
)
    if c == 0 || n_segments == 0
        # No Sundman transformation: equal time segments
        return fill(tof_leg / max(n_segments, 1), max(n_segments, 1))
    end

    # Interpolate radius at each segment midpoint
    # We use linear interpolation in position space, which is an approximation
    # The actual trajectory is curved, but this gives a reasonable estimate
    r_start_mag = norm(r_start)
    r_end_mag = norm(r_end)

    # Compute weights: w_i = r_i^c (segment at larger r gets more time)
    # This implements dt = r^c ds, where ds is uniform
    T = typeof(tof_leg)
    weights = Vector{T}(undef, n_segments)
    for i = 1:n_segments
        # Midpoint of segment i in normalized leg [0, 1]
        s_mid = (i - T(0.5)) / n_segments
        # Linear interpolation of radius magnitude
        r_mid = r_start_mag + s_mid * (r_end_mag - r_start_mag)
        weights[i] = r_mid^c
    end

    # Normalize weights to sum to tof_leg
    total_weight = sum(weights)
    Δt_segments = weights .* (tof_leg / total_weight)

    return Δt_segments
end

"""
    compute_sundman_leg_times(r0, rf, tof, n_fwd, n_bwd; c=1.0)

Compute segment durations for forward and backward legs using Sundman transformation.

This function:
1. Estimates the match point position based on segment distribution
2. Allocates total TOF between legs proportionally to their Sundman-weighted lengths
3. Computes individual segment times within each leg

# Arguments
- `r0::AbstractVector`: Initial position [km]
- `rf::AbstractVector`: Final position [km]
- `tof::Number`: Total time of flight [s]
- `n_fwd::Int`: Number of forward segments
- `n_bwd::Int`: Number of backward segments
- `c::Number`: Sundman exponent

# Returns
- `(Δt_fwd, Δt_bwd)`: Segment duration vectors for forward and backward legs

# Notes
For c > 0, segments at larger radius get longer durations. This is appropriate for:
- SEP where thrust decreases with distance (more time needed at outer regions)
- Better resolution near perihelion where dynamics are faster

The match point estimation uses linear interpolation between r0 and rf, which is
approximate but avoids the need for iterative refinement.
"""
function compute_sundman_leg_times(
    r0::AbstractVector,
    rf::AbstractVector,
    tof::Number,
    n_fwd::Int,
    n_bwd::Int;
    c::Number = 1.0,
)
    n_total = n_fwd + n_bwd

    if c == 0 || n_total == 0
        # No Sundman: equal time segments
        Δt_seg = tof / max(n_total, 1)
        return fill(Δt_seg, n_fwd), fill(Δt_seg, n_bwd)
    end

    # Compute segment times for the ENTIRE trajectory first (all n_total segments)
    # This ensures consistency between forward and backward legs
    T = typeof(tof)
    r0_mag = norm(r0)
    rf_mag = norm(rf)

    # Compute weights for all segments
    all_weights = Vector{T}(undef, n_total)
    for i = 1:n_total
        # Midpoint of segment i in normalized trajectory [0, 1]
        s_mid = (i - T(0.5)) / n_total
        # Linear interpolation of radius magnitude
        r_mid = r0_mag + s_mid * (rf_mag - r0_mag)
        all_weights[i] = r_mid^c
    end

    # Normalize to get segment times
    total_weight = sum(all_weights)
    all_Δt = all_weights .* (tof / total_weight)

    # Split into forward and backward legs
    Δt_fwd = all_Δt[1:n_fwd]
    Δt_bwd = all_Δt[(n_fwd+1):end]

    return Δt_fwd, Δt_bwd
end

# =============================================================================
# Universal Kepler Propagation
# =============================================================================

#TODO: THIS DEFINITELY DOESNT BELONG HERE
"""
    kepler_propagate(r0, v0, Δt, μ)

Propagate a Keplerian orbit from initial state `(r0, v0)` for time `Δt`.

Uses universal variables (Stumpff functions) which:
- Handle all orbit types (elliptic, parabolic, hyperbolic)
- Have no singularities for circular or equatorial orbits
- Work correctly with automatic differentiation

# Arguments
- `r0::SVector{3}`: Initial position [km]
- `v0::SVector{3}`: Initial velocity [km/s]
- `Δt::Number`: Time of flight [s]
- `μ::Number`: Gravitational parameter [km³/s²]

# Returns
- `(rf, vf)`: Final position [km] and velocity [km/s]
"""
function kepler_propagate(
    r0::AbstractVector{RT},
    v0::AbstractVector{VT},
    Δt::TT,
    μ::MT;
    tol::Float64 = 1e-12,
    maxiter::Int = 50,
) where {RT<:Number,VT<:Number,TT<:Number,MT<:Number}

    T = promote_type(RT, VT, TT, MT)

    # Initial state quantities
    r0_mag = sqrt(r0[1]^2 + r0[2]^2 + r0[3]^2)
    v0_mag_sq = v0[1]^2 + v0[2]^2 + v0[3]^2
    rdotv = r0[1]*v0[1] + r0[2]*v0[2] + r0[3]*v0[3]

    # Semi-major axis (negative for hyperbolic)
    α = 2/r0_mag - v0_mag_sq/μ  # α = 1/a

    # Initial guess for universal anomaly χ
    χ = sqrt(μ) * abs(α) * Δt  # Good initial guess

    # Newton-Raphson iteration to solve universal Kepler's equation
    for _ = 1:maxiter
        χ2 = χ * χ
        z = α * χ2

        c2 = stumpff_c2(z)
        c3 = stumpff_c3(z)

        # Universal Kepler's equation: F(χ) = 0
        r = χ2 * c2 + rdotv / sqrt(μ) * χ * (1 - z * c3) + r0_mag * (1 - z * c2)
        F =
            rdotv / sqrt(μ) * χ2 * c2 + (1 - r0_mag * α) * χ^3 * c3 + r0_mag * χ -
            sqrt(μ) * Δt

        # Derivative dF/dχ = r
        dF = r

        # Newton step
        δχ = F / dF
        χ = χ - δχ

        if abs(δχ) < tol
            break
        end
    end

    # Compute final Stumpff values
    χ2 = χ * χ
    z = α * χ2
    c2 = stumpff_c2(z)
    c3 = stumpff_c3(z)

    # Final radius
    r_mag = χ2 * c2 + rdotv / sqrt(μ) * χ * (1 - z * c3) + r0_mag * (1 - z * c2)

    # Lagrange coefficients
    f = one(T) - χ2 / r0_mag * c2
    g = Δt - χ^3 / sqrt(μ) * c3
    fdot = sqrt(μ) / (r_mag * r0_mag) * χ * (z * c3 - one(T))
    gdot = one(T) - χ2 / r_mag * c2

    # Final state
    rf = SVector{3,T}(f * r0[1] + g * v0[1], f * r0[2] + g * v0[2], f * r0[3] + g * v0[3])
    vf = SVector{3,T}(
        fdot * r0[1] + gdot * v0[1],
        fdot * r0[2] + gdot * v0[2],
        fdot * r0[3] + gdot * v0[3],
    )

    return rf, vf
end

"""
    propagate_segment(r, v, m, throttle, Δt, μ, spacecraft; forward=true)

Propagate a single Sims-Flanagan segment.

The segment is divided into two halves:
1. Coast for Δt/2 from initial state
2. Apply impulsive ΔV at midpoint
3. Coast for Δt/2 to final state

# Arguments
- `r::SVector{3}`: Position at segment start [km]
- `v::SVector{3}`: Velocity at segment start [km/s]
- `m::Number`: Mass at segment start [kg]
- `throttle::SVector{3}`: Throttle vector (magnitude ≤ 1, direction = thrust direction)
- `Δt::Number`: Segment duration [s]
- `μ::Number`: Gravitational parameter [km³/s²]
- `spacecraft::AbstractSpacecraft`: Spacecraft/propulsion parameters
- `forward::Bool`: If true, propagate forward; if false, propagate backward

# Returns
- `(rf, vf, mf, Δv)`: Final state and ΔV magnitude [km/s]
"""
function propagate_segment(
    r::AbstractVector{RT},
    v::AbstractVector{VT},
    m::MT,
    throttle::AbstractVector{ThT},
    Δt::DTT,
    μ::MuT,
    spacecraft::AbstractSpacecraft;
    forward::Val{direction} = Val(true),
) where {RT<:Number,VT<:Number,MT<:Number,ThT<:Number,DTT<:Number,MuT<:Number,direction}

    T = promote_type(RT, VT, MT, ThT, DTT, MuT)

    half_Δt = Δt / 2

    # Throttle magnitude (clamped to [0, 1])
    throttle_mag = clamp(safe_norm(throttle), zero(ThT), one(ThT))

    # Coast to midpoint first to get position for thrust calculation
    if direction
        r_mid, v_mid = kepler_propagate(r, v, half_Δt, μ)
    else
        r_mid, v_mid = kepler_propagate(r, v, -half_Δt, μ)
    end

    # Compute thrust at midpoint position (important for SEP and solar sails)
    thrust = compute_thrust(spacecraft, r_mid, throttle_mag)  # [N]

    # Mass flow rate [kg/s]
    mdot = compute_mass_flow(spacecraft, thrust)

    # Mass consumed during segment
    Δm = mdot * Δt

    # Clamp mass change to avoid consuming dry mass (only for forward propagation)
    if has_propellant_consumption(spacecraft) && direction
        # Forward: can only consume available propellant, not dry mass
        available_propellant = m - spacecraft.dry_mass
        Δm = min(Δm, max(available_propellant, zero(T)))
    elseif !has_propellant_consumption(spacecraft)
        Δm = zero(T)
    end
    # For backward propagation: no clamping needed, we're adding mass back

    # Average mass during burn (use current mass for solar sails)
    # For forward propagation: m is mass at segment START (larger), m_avg = m - Δm/2
    # For backward propagation: m is mass at segment END (smaller), m_avg = m + Δm/2
    if has_propellant_consumption(spacecraft)
        m_avg = direction ? m - Δm / 2 : m + Δm / 2
    else
        m_avg = m
    end

    # ΔV magnitude [m/s -> km/s]
    Δv_mag = (thrust * Δt / m_avg) / 1000  # Convert m/s to km/s

    # ΔV vector (in thrust direction)
    # Use safe division to avoid 0/0 = NaN when throttle is zero
    if throttle_mag > eps(T)
        Δv_vec = throttle * (Δv_mag / throttle_mag)
    else
        Δv_vec = @SVector zeros(T, 3)
    end

    if direction
        # Forward propagation: we already coasted to midpoint, apply impulse, coast rest
        v_mid_plus = v_mid + Δv_vec
        rf, vf = kepler_propagate(r_mid, v_mid_plus, half_Δt, μ)
        mf = m - Δm
    else
        # Backward propagation: we already coasted backward to midpoint, apply negative impulse
        v_mid_minus = v_mid - Δv_vec  # Subtract because we're going backward
        rf, vf = kepler_propagate(r_mid, v_mid_minus, -half_Δt, μ)
        mf = m + Δm  # Mass increases going backward
    end

    return rf, vf, mf, Δv_mag
end

"""
    propagate_leg(r0, v0, m0, throttles, Δt_segments, μ, spacecraft; forward=true)

Propagate multiple segments of a Sims-Flanagan leg.

# Arguments
- `r0::SVector{3}`: Initial position [km]
- `v0::SVector{3}`: Initial velocity [km/s]
- `m0::Number`: Initial mass [kg]
- `throttles::AbstractVector{SVector{3}}`: Throttle vectors for each segment
- `Δt_segments::Union{Number, AbstractVector}`: Duration of each segment [s]
  - If scalar: all segments have equal duration
  - If vector: segment i has duration Δt_segments[i] (for Sundman transformation)
- `μ::Number`: Gravitational parameter [km³/s²]
- `spacecraft::AbstractSpacecraft`: Spacecraft/propulsion parameters
- `forward::Bool`: If true, propagate forward; if false, propagate backward

# Returns
- `(rf, vf, mf, total_Δv)`: Final state and total ΔV [km/s]
"""
function propagate_leg(
    r0::AbstractVector{RT},
    v0::AbstractVector{VT},
    m0::MT,
    throttles::AbstractVector{<:AbstractVector{ThT}},
    Δt_segments::Union{DTT,AbstractVector{DTT}},
    μ::MuT,
    spacecraft::AbstractSpacecraft;
    forward::Val{direction} = Val(true),
) where {RT<:Number,VT<:Number,MT<:Number,ThT<:Number,DTT<:Number,MuT<:Number,direction}

    T = promote_type(RT, VT, MT, ThT, DTT, MuT)

    r, v, m = r0, v0, m0
    total_Δv = zero(T)

    n_seg = length(throttles)
    for i = 1:n_seg
        # Get segment duration (scalar or from array)
        Δt = Δt_segments isa Number ? Δt_segments : Δt_segments[i]

        r, v, m, Δv =
            propagate_segment(r, v, m, throttles[i], Δt, μ, spacecraft; forward = forward)
        total_Δv += Δv
    end

    return r, v, m, total_Δv
end

"""
    compute_mismatch(problem, throttles)

Compute the mismatch at the match point.

The trajectory is propagated forward from the initial state and backward
from the final state. The mismatch is the difference at the match point.

If `sundman_c > 0` in problem options, uses Sundman transformation for 
adaptive segment durations based on distance from central body.

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition
- `throttles::AbstractVector{SVector{3}}`: Throttle vectors for all segments

# Returns
- `mismatch::SVector{7}`: Mismatch vector [Δr; Δv; Δm] at match point
"""
function compute_mismatch(
    problem::SimsFlanaganProblem{T},
    throttles::AbstractVector{<:AbstractVector},
) where {T<:Number}

    n_seg = problem.options.n_segments
    n_fwd = problem.options.n_fwd
    n_bwd = n_seg - n_fwd
    sundman_c = problem.options.sundman_c

    # Compute segment durations for each leg (with optional Sundman transformation)
    # This allocates TOF between legs and computes segment times within each leg
    Δt_fwd, Δt_bwd = compute_sundman_leg_times(
        problem.r0,
        problem.rf,
        problem.tof,
        n_fwd,
        n_bwd;
        c = sundman_c,
    )

    # Forward propagation
    throttles_fwd = throttles[1:n_fwd]
    r_fwd, v_fwd, m_fwd, _ = propagate_leg(
        problem.r0,
        problem.v0,
        mass(problem.spacecraft),
        throttles_fwd,
        Δt_fwd,
        problem.μ,
        problem.spacecraft,
        forward = Val(true),
    )

    # Backward propagation
    throttles_bwd = throttles[(n_fwd+1):end]

    # Estimate final mass using all segment times
    Δt_all = vcat(Δt_fwd, Δt_bwd)
    mf_estimate = estimate_final_mass(problem, throttles, Δt_all)

    # Reverse the backward throttles and times for proper ordering
    r_bwd, v_bwd, m_bwd, _ = propagate_leg(
        problem.rf,
        problem.vf,
        mf_estimate,
        reverse(throttles_bwd),
        reverse(Δt_bwd),
        problem.μ,
        problem.spacecraft,
        forward = Val(false),
    )

    # Compute mismatch
    Δr = r_fwd - r_bwd
    Δv = v_fwd - v_bwd
    Δm = m_fwd - m_bwd

    # Infer element type from computed values (supports ForwardDiff.Dual)
    ET = promote_type(eltype(Δr), eltype(Δv), typeof(Δm))
    return SVector{7,ET}(Δr[1], Δr[2], Δr[3], Δv[1], Δv[2], Δv[3], Δm)
end

"""
    estimate_final_mass(problem, throttles, Δt_segments)

Estimate the final mass for backward propagation initialization.

For propellant-consuming systems, estimates mass consumption.
For solar sails, mass is constant.

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition
- `throttles::AbstractVector`: All throttle vectors
- `Δt_segments::AbstractVector`: Duration of each segment [s] (optional, uses equal segments if not provided)
"""
function estimate_final_mass(
    problem::SimsFlanaganProblem,
    throttles::AbstractVector,
    Δt_segments::Union{Nothing,AbstractVector} = nothing,
)
    spacecraft = problem.spacecraft

    if !has_propellant_consumption(spacecraft)
        # Solar sails have constant mass
        return mass(spacecraft)
    end

    # For thrust systems, estimate mass consumption
    vex = exhaust_velocity(spacecraft)

    # Approximate using reference thrust (for SEP, this is conservative)
    thrust_ref = if spacecraft isa SEPSpacecraft
        spacecraft.thrust_ref
    else
        spacecraft.thrust
    end

    # Use provided segment times or equal segments
    if Δt_segments === nothing
        Δt_seg = problem.tof / problem.options.n_segments
        total_mass_flow = sum(safe_norm(t) * thrust_ref / vex * Δt_seg for t in throttles)
    else
        total_mass_flow = sum(
            safe_norm(throttles[i]) * thrust_ref / vex * Δt_segments[i] for
            i in eachindex(throttles)
        )
    end

    return max(mass(spacecraft) - total_mass_flow, spacecraft.dry_mass)
end

"""
    compute_total_Δv(throttles, Δt_segments, spacecraft)

Compute the total ΔV for a given set of throttles.

# Arguments
- `throttles::AbstractVector{SVector{3}}`: Throttle vectors for all segments
- `Δt_segments::Union{Number, AbstractVector}`: Duration of each segment [s]
  - If scalar: all segments have equal duration
  - If vector: segment i has duration Δt_segments[i]
- `spacecraft::AbstractSpacecraft`: Spacecraft/propulsion parameters

# Returns
- `total_Δv::Number`: Total ΔV [km/s]
"""
function compute_total_Δv(
    throttles::AbstractVector{<:AbstractVector{ThT}},
    Δt_segments::Union{DT,AbstractVector{DT}},
    spacecraft::AbstractSpacecraft;
    r_ref::AbstractVector = SVector{3,Float64}(1.495978707e8, 0.0, 0.0),  # Reference position for thrust calc
) where {ThT<:Number,DT<:Number}

    T = promote_type(ThT, DT, typeof(mass(spacecraft)))
    total_Δv = zero(T)
    m = T(mass(spacecraft))

    for (i, throttle) in enumerate(throttles)
        # Get segment duration (scalar or from array)
        Δt_seg = Δt_segments isa Number ? Δt_segments : Δt_segments[i]

        throttle_mag = clamp(safe_norm(throttle), zero(T), one(T))

        # Compute thrust (use reference position for consistency)
        thrust = compute_thrust(spacecraft, r_ref, throttle_mag)

        # Mass flow
        mdot = compute_mass_flow(spacecraft, thrust)
        Δm = mdot * Δt_seg

        if has_propellant_consumption(spacecraft)
            # Can only consume propellant, not dry mass
            available_propellant = m - spacecraft.dry_mass
            Δm = min(Δm, max(available_propellant, zero(T)))
            m_avg = m - Δm / 2
        else
            Δm = zero(T)
            m_avg = m
        end

        Δv = (thrust * Δt_seg / m_avg) / 1000  # Convert to km/s
        total_Δv += Δv
        m -= Δm
    end

    return total_Δv
end

# Backward-compatible version for constant-thrust spacecraft
function compute_total_Δv(
    throttles::AbstractVector{<:AbstractVector{ThT}},
    Δt_segments::Union{DT,AbstractVector{DT}},
    spacecraft::Spacecraft,
) where {ThT<:Number,DT<:Number}

    T = promote_type(ThT, DT, typeof(mass(spacecraft)))
    total_Δv = zero(T)
    m = mass(spacecraft)
    dry_mass = spacecraft.dry_mass
    vex = exhaust_velocity(spacecraft)

    for (i, throttle) in enumerate(throttles)
        # Get segment duration (scalar or from array)
        Δt_seg = Δt_segments isa Number ? Δt_segments : Δt_segments[i]

        throttle_mag = clamp(safe_norm(throttle), zero(T), one(T))
        thrust = spacecraft.thrust * throttle_mag
        mdot = thrust / vex
        Δm = mdot * Δt_seg
        # Can only consume propellant, not dry mass
        available_propellant = m - dry_mass
        Δm = min(Δm, max(available_propellant, zero(T)))
        m_avg = m - Δm / 2
        Δv = (thrust * Δt_seg / m_avg) / 1000  # Convert to km/s
        total_Δv += Δv
        m -= Δm
    end

    return total_Δv
end
