export kepler_propagate, propagate_segment, propagate_leg
export compute_mismatch, compute_total_Δv

#=
Propagation functions for Sims-Flanagan transcription.

The Sims-Flanagan method uses Kepler propagation between impulsive ΔV maneuvers
applied at segment midpoints. This file provides the core propagation routines
using AstroCoords for coordinate transformations.
=#

"""
    kepler_propagate(r0, v0, Δt, μ)

Propagate a Keplerian orbit from initial state `(r0, v0)` for time `Δt`.

Uses AstroCoords coordinate transformations: Cartesian → Keplerian → propagate 
mean anomaly → Keplerian → Cartesian.

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
    μ::MT,
) where {RT<:Number,VT<:Number,TT<:Number,MT<:Number}

    T = promote_type(RT, VT, TT, MT)

    # Create Cartesian state and convert to Keplerian
    cart = Cartesian(r0[1], r0[2], r0[3], v0[1], v0[2], v0[3])
    kep = Keplerian(cart, μ)

    # Extract orbital elements
    a = kep.a      # Semi-major axis [km]
    e = kep.e      # Eccentricity
    i = kep.i      # Inclination [rad]
    Ω = kep.Ω      # RAAN [rad]
    ω = kep.ω      # Argument of periapsis [rad]
    f = kep.f      # True anomaly [rad]

    # Compute mean motion [rad/s]
    n = sqrt(μ / abs(a)^3)

    # Convert true anomaly to mean anomaly
    M0 = trueAnomaly2MeanAnomaly(f, e)

    # Propagate mean anomaly
    M = M0 + n * Δt

    # Convert back to true anomaly
    f_new = meanAnomaly2TrueAnomaly(M, e)

    # Create new Keplerian state and convert to Cartesian
    kep_new = Keplerian(a, e, i, Ω, ω, f_new)
    cart_new = Cartesian(kep_new, μ)

    rf = SVector{3,T}(cart_new[1], cart_new[2], cart_new[3])
    vf = SVector{3,T}(cart_new[4], cart_new[5], cart_new[6])

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
- `spacecraft::Spacecraft`: Spacecraft parameters
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
    spacecraft::Spacecraft;
    forward::Val{direction} = Val(true),
    cutoff_threshold::Number = 10.0 * eps(promote_type(RT, VT, MT, ThT, DTT, MuT)),
) where {RT<:Number,VT<:Number,MT<:Number,ThT<:Number,DTT<:Number,MuT<:Number,direction}

    T = promote_type(RT, VT, MT, ThT, DTT, MuT)

    half_Δt = Δt / 2

    # Throttle magnitude (clamped to [0, 1])
    throttle_mag = min(norm(throttle), one(ThT))

    # Compute ΔV from rocket equation approximation
    vex = exhaust_velocity(spacecraft)  # [m/s]
    thrust = spacecraft.thrust * throttle_mag  # Actual thrust [N]

    # Mass flow rate: ṁ = thrust / vex
    mdot = thrust / vex  # [kg/s]

    # Mass consumed during segment
    Δm = mdot * Δt

    # Clamp mass change to avoid negative mass
    Δm = min(Δm, m)

    # Average mass during burn
    m_avg = m - Δm / 2

    # ΔV magnitude [m/s -> km/s]
    Δv_mag = (thrust * Δt / m_avg) / 1000  # Convert m/s to km/s

    # ΔV vector (in thrust direction)
    if throttle_mag > cutoff_threshold
        Δv_vec = throttle * (Δv_mag / throttle_mag)
    else
        Δv_vec = @SVector zeros(T, 3)
    end

    if direction
        # Forward propagation: coast half, impulse, coast half
        r1, v1 = kepler_propagate(r, v, half_Δt, μ)
        v1_plus = v1 + Δv_vec
        rf, vf = kepler_propagate(r1, v1_plus, half_Δt, μ)
        mf = m - Δm
    else
        # Backward propagation: coast half backward, impulse (negative), coast half backward
        r1, v1 = kepler_propagate(r, v, -half_Δt, μ)
        v1_minus = v1 - Δv_vec  # Subtract because we're going backward
        rf, vf = kepler_propagate(r1, v1_minus, -half_Δt, μ)
        mf = m + Δm  # Mass increases going backward
    end

    return rf, vf, mf, Δv_mag
end

"""
    propagate_leg(r0, v0, m0, throttles, Δt_seg, μ, spacecraft; forward=true)

Propagate multiple segments of a Sims-Flanagan leg.

# Arguments
- `r0::SVector{3}`: Initial position [km]
- `v0::SVector{3}`: Initial velocity [km/s]
- `m0::Number`: Initial mass [kg]
- `throttles::AbstractVector{SVector{3}}`: Throttle vectors for each segment
- `Δt_seg::Number`: Duration of each segment [s]
- `μ::Number`: Gravitational parameter [km³/s²]
- `spacecraft::Spacecraft`: Spacecraft parameters
- `forward::Bool`: If true, propagate forward; if false, propagate backward

# Returns
- `(rf, vf, mf, total_Δv)`: Final state and total ΔV [km/s]
"""
function propagate_leg(
    r0::AbstractVector{RT},
    v0::AbstractVector{VT},
    m0::MT,
    throttles::AbstractVector{<:AbstractVector{ThT}},
    Δt_seg::DTT,
    μ::MuT,
    spacecraft::Spacecraft;
    forward::Val{direction} = Val(true),
    cutoff_threshold::Number = eps(promote_type(RT, VT, MT, ThT, DTT, MuT)),
) where {RT<:Number,VT<:Number,MT<:Number,ThT<:Number,DTT<:Number,MuT<:Number,direction}

    T = promote_type(RT, VT, MT, ThT, DTT, MuT)

    r, v, m = r0, v0, m0
    total_Δv = zero(T)

    for throttle in throttles
        r, v, m, Δv = propagate_segment(
            r,
            v,
            m,
            throttle,
            Δt_seg,
            μ,
            spacecraft;
            forward = forward,
            cutoff_threshold = cutoff_threshold,
        )
        total_Δv += Δv
    end

    return r, v, m, total_Δv
end

"""
    compute_mismatch(problem, throttles)

Compute the mismatch at the match point.

The trajectory is propagated forward from the initial state and backward
from the final state. The mismatch is the difference at the match point.

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

    Δt_seg = problem.tof / n_seg

    # Forward propagation
    throttles_fwd = throttles[1:n_fwd]
    r_fwd, v_fwd, m_fwd, _ = propagate_leg(
        problem.r0,
        problem.v0,
        problem.spacecraft.mass,
        throttles_fwd,
        Δt_seg,
        problem.μ,
        problem.spacecraft,
        forward = Val(true),
    )

    # Backward propagation
    throttles_bwd = throttles[(n_fwd+1):end]

    # Estimate final mass (simple approximation for backward propagation start)
    vex = exhaust_velocity(problem.spacecraft)
    total_mass_flow =
        sum(norm(t) * problem.spacecraft.thrust / vex * Δt_seg for t in throttles)
    mf_estimate = max(problem.spacecraft.mass - total_mass_flow, problem.spacecraft.mass)

    # Reverse the backward throttles for proper ordering
    r_bwd, v_bwd, m_bwd, _ = propagate_leg(
        problem.rf,
        problem.vf,
        mf_estimate,
        reverse(throttles_bwd),
        Δt_seg,
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
    compute_total_Δv(throttles, Δt_seg, spacecraft)

Compute the total ΔV for a given set of throttles.

# Arguments
- `throttles::AbstractVector{SVector{3}}`: Throttle vectors for all segments
- `Δt_seg::Number`: Duration of each segment [s]
- `spacecraft::Spacecraft`: Spacecraft parameters

# Returns
- `total_Δv::Number`: Total ΔV [km/s]
"""
function compute_total_Δv(
    throttles::AbstractVector{<:AbstractVector{ThT}},
    Δt_seg::DT,
    spacecraft::Spacecraft,
) where {ThT<:Number,DT<:Number}

    T = promote_type(ThT, DT, typeof(spacecraft.mass))
    total_Δv = zero(T)
    m = spacecraft.mass
    vex = exhaust_velocity(spacecraft)

    for throttle in throttles
        throttle_mag = min(norm(throttle), one(T))
        thrust = spacecraft.thrust * throttle_mag
        mdot = thrust / vex
        Δm = mdot * Δt_seg
        Δm = min(Δm, m)
        m_avg = m - Δm / 2
        Δv = (thrust * Δt_seg / m_avg) / 1000  # Convert to km/s
        total_Δv += Δv
        m -= Δm
    end

    return total_Δv
end

