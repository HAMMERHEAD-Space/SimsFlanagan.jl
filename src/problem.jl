export simsflanagan_problem, initial_guess_lambert, initial_guess_zero

"""
    simsflanagan_problem(r0, v0, rf, vf, tof, μ, spacecraft; kwargs...)

Construct a Sims-Flanagan trajectory optimization problem.

# Arguments
- `r0::AbstractVector`: Initial position [km]
- `v0::AbstractVector`: Initial velocity [km/s]
- `rf::AbstractVector`: Final position [km]
- `vf::AbstractVector`: Final velocity [km/s]
- `tof::Number`: Time of flight [s]
- `μ::Number`: Gravitational parameter [km³/s²]
- `spacecraft::AbstractSpacecraft`: Spacecraft/propulsion parameters

# Keyword Arguments
- `n_segments::Int=10`: Number of thrust segments
- `n_fwd::Int=n_segments÷2`: Number of forward propagation segments
- `tol::Float64=1e-8`: Convergence tolerance
- `max_iter::Int=1000`: Maximum iterations
- `verbosity::Int=1`: Output verbosity level

# Returns
- `SimsFlanaganProblem`: The problem definition

# Supported Spacecraft Types
- `Spacecraft`: Constant thrust (chemical/electric propulsion)
- `SEPSpacecraft`: Solar Electric Propulsion (thrust varies with solar distance)
- `SolarSail`: Solar radiation pressure propulsion
"""
function simsflanagan_problem(
    r0::AbstractVector{RT1},
    v0::AbstractVector{VT1},
    rf::AbstractVector{RT2},
    vf::AbstractVector{VT2},
    tof::TT,
    μ::MT,
    spacecraft::AbstractSpacecraft;
    n_segments::Int = 10,
    n_fwd::Int = n_segments ÷ 2,
    tol::Float64 = 1e-8,
    max_iter::Int = 1000,
    verbosity::Int = 1,
    cutoff_threshold::Number = 10.0 * eps(promote_type(RT1, VT1, RT2, VT2, TT, MT)),
) where {RT1<:Number,VT1<:Number,RT2<:Number,VT2<:Number,TT<:Number,MT<:Number}

    options = SimsFlanaganOptions(;
        n_segments = n_segments,
        n_fwd = n_fwd,
        tol = tol,
        max_iter = max_iter,
        verbosity = verbosity,
        cutoff_threshold = cutoff_threshold,
    )

    return SimsFlanaganProblem(
        SVector{3}(r0[1], r0[2], r0[3]),
        SVector{3}(v0[1], v0[2], v0[3]),
        SVector{3}(rf[1], rf[2], rf[3]),
        SVector{3}(vf[1], vf[2], vf[3]),
        tof,
        μ,
        spacecraft,
        options,
    )
end

"""
    initial_guess_lambert(problem)

Generate an initial guess for throttle vectors using Lambert's problem.

Solves a Lambert problem between initial and final positions, then distributes
the required ΔV across segments proportionally.

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition

# Returns
- `throttles::Vector{SVector{3}}`: Initial guess for throttle vectors

# Notes
For solar sails, returns a zero guess since Lambert-based initialization
doesn't apply to radiation pressure propulsion.
"""
function initial_guess_lambert(problem::SimsFlanaganProblem{T}) where {T}
    n_seg = problem.options.n_segments
    spacecraft = problem.spacecraft

    # Solar sails can't use Lambert-based initialization
    if spacecraft isa SolarSail
        return initial_guess_zero(problem)
    end

    # Solve Lambert problem
    lambert_prob = LambertProblem(problem.μ, problem.r0, problem.rf, problem.tof)

    lambert_sol = Lambert.solve(lambert_prob)

    if lambert_sol.retcode != :SUCCESS
        throw(ErrorException("Lambert problem failed to solve"))
    end

    # Compute departure and arrival ΔVs
    Δv_dep =
        SVector{3,T}(lambert_sol.v1[1], lambert_sol.v1[2], lambert_sol.v1[3]) - problem.v0
    Δv_arr =
        problem.vf - SVector{3,T}(lambert_sol.v2[1], lambert_sol.v2[2], lambert_sol.v2[3])

    total_Δv_vec = Δv_dep + Δv_arr
    Δv_per_seg = total_Δv_vec / n_seg

    # Get reference thrust for throttle normalization
    ref_thrust = get_reference_thrust(spacecraft)

    # Convert to throttle (normalize by max ΔV capacity per segment)
    Δt_seg = problem.tof / n_seg
    max_Δv_seg = (ref_thrust * Δt_seg / spacecraft.mass) / 1000  # km/s

    throttle_mag = norm(Δv_per_seg) / max_Δv_seg
    throttle_mag = min(throttle_mag, 1.0)  # Clamp to [0, 1]

    if norm(Δv_per_seg) > eps(T)
        throttle_dir = Δv_per_seg / norm(Δv_per_seg)
        throttle = throttle_mag * throttle_dir
    else
        throttle = SVector{3,T}(0.0, 0.0, 0.0)
    end

    return [throttle for _ = 1:n_seg]
end

"""
    get_reference_thrust(spacecraft)

Get a reference thrust value for initial guess normalization.
"""
function get_reference_thrust(spacecraft::Spacecraft)
    return spacecraft.thrust
end

function get_reference_thrust(spacecraft::SEPSpacecraft)
    return spacecraft.thrust_ref
end

function get_reference_thrust(spacecraft::SolarSail)
    # Return characteristic thrust at reference distance
    a_c = characteristic_acceleration(spacecraft)
    return spacecraft.mass * a_c  # F = ma
end

"""
    initial_guess_zero(problem)

Generate a zero initial guess for throttle vectors.

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition

# Returns
- `throttles::Vector{SVector{3}}`: Zero throttle vectors
"""
function initial_guess_zero(problem::SimsFlanaganProblem{T}) where {T}
    n_seg = problem.options.n_segments
    return [SVector{3,T}(0.0, 0.0, 0.0) for _ = 1:n_seg]
end

"""
    initial_guess_radial(problem)

Generate an initial guess with throttles pointing radially (for solar sails).

For solar sails, a reasonable starting point is to thrust radially outward
or along the velocity direction.

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition

# Returns
- `throttles::Vector{SVector{3}}`: Radial throttle vectors
"""
function initial_guess_radial(
    problem::SimsFlanaganProblem{T};
    magnitude::Number = 0.5,
) where {T}
    n_seg = problem.options.n_segments

    # Use direction from initial position (approximate Sun direction)
    r_dir = problem.r0 / norm(problem.r0)
    throttle =
        SVector{3,T}(magnitude * r_dir[1], magnitude * r_dir[2], magnitude * r_dir[3])

    return [throttle for _ = 1:n_seg]
end
