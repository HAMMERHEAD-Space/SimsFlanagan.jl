export simsflanagan_problem
export initial_guess_lambert, initial_guess_zero, initial_guess_radial
export initial_guess_random, initial_guess_constant
export generate_initial_guess

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
- `tol::Float64=1e-6`: Convergence tolerance
- `max_iter::Int=1000`: Maximum iterations
- `verbosity::Int=1`: Output verbosity level
- `sundman_c::Float64=0.0`: Sundman transformation exponent
  - 0.0: Equal time segments (default, recommended for initial debugging)
  - 1.0: Classic Sundman (dt ∝ r), more segments near central body
  - 1.5: Related to eccentric anomaly

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
    tol::Float64 = 1e-6,
    max_iter::Int = 1000,
    verbosity::Int = 1,
    sundman_c::Float64 = 0.0,  # Default to equal segments
) where {RT1<:Number,VT1<:Number,RT2<:Number,VT2<:Number,TT<:Number,MT<:Number}

    options = SimsFlanaganOptions(;
        n_segments = n_segments,
        n_fwd = n_fwd,
        tol = tol,
        max_iter = max_iter,
        verbosity = verbosity,
        sundman_c = sundman_c,
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
the required ΔV across segments. Departure ΔV is applied to forward segments,
arrival ΔV to backward segments.

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
    n_fwd = problem.options.n_fwd
    n_bwd = n_seg - n_fwd
    spacecraft = problem.spacecraft

    # Solar sails can't use Lambert-based initialization
    if spacecraft isa SolarSail
        return initial_guess_radial(problem; magnitude = 0.5)
    end

    # Solve Lambert problem
    lambert_prob = Lambert.LambertProblem(problem.μ, problem.r0, problem.rf, problem.tof)

    # Check for zero transfer angle and perturb if needed
    transfer_angle = AstroCoords.angle_between_vectors(problem.r0, problem.rf)
    if transfer_angle < T(1e-6) || transfer_angle > T(π - 1e-6)
        # Add small perpendicular perturbation (1 km)
        rf_hat = normalize(problem.rf)
        perp = abs(rf_hat[1]) < T(0.9) ? SVector{3,T}(1, 0, 0) : SVector{3,T}(0, 1, 0)
        perp = normalize(perp - dot(perp, rf_hat) * rf_hat)
        rf_perturbed = problem.rf + perp * one(T)
        lambert_prob =
            Lambert.LambertProblem(problem.μ, problem.r0, rf_perturbed, problem.tof)
    end

    lambert_sol = Lambert.solve(lambert_prob)

    if lambert_sol.retcode != :SUCCESS
        # Fall back to zero guess if Lambert fails
        return initial_guess_zero(problem)
    end

    # Compute departure and arrival ΔVs
    v1_lambert = SVector{3,T}(lambert_sol.v1[1], lambert_sol.v1[2], lambert_sol.v1[3])
    v2_lambert = SVector{3,T}(lambert_sol.v2[1], lambert_sol.v2[2], lambert_sol.v2[3])

    Δv_dep = v1_lambert - problem.v0  # ΔV needed at departure
    Δv_arr = problem.vf - v2_lambert  # ΔV needed at arrival

    # Get reference thrust for throttle normalization
    ref_thrust = get_reference_thrust(spacecraft)
    Δt_seg = problem.tof / n_seg
    max_Δv_seg = (ref_thrust * Δt_seg / mass(spacecraft)) / 1000  # km/s

    throttles = Vector{SVector{3,T}}(undef, n_seg)

    # Distribute departure ΔV over forward segments
    if n_fwd > 0
        Δv_per_fwd = Δv_dep / n_fwd
        throttle_mag_fwd = min(norm(Δv_per_fwd) / max_Δv_seg, one(T))
        if norm(Δv_per_fwd) > eps(T)
            throttle_fwd = throttle_mag_fwd * (Δv_per_fwd / norm(Δv_per_fwd))
        else
            throttle_fwd = SVector{3,T}(0.0, 0.0, 0.0)
        end
        for i = 1:n_fwd
            throttles[i] = throttle_fwd
        end
    end

    # Distribute arrival ΔV over backward segments
    if n_bwd > 0
        Δv_per_bwd = Δv_arr / n_bwd
        throttle_mag_bwd = min(norm(Δv_per_bwd) / max_Δv_seg, one(T))
        if norm(Δv_per_bwd) > eps(T)
            throttle_bwd = throttle_mag_bwd * (Δv_per_bwd / norm(Δv_per_bwd))
        else
            throttle_bwd = SVector{3,T}(0.0, 0.0, 0.0)
        end
        for i = (n_fwd+1):n_seg
            throttles[i] = throttle_bwd
        end
    end

    return throttles
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
    return mass(spacecraft) * a_c  # F = ma
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
    initial_guess_radial(problem; magnitude=0.5)

Generate an initial guess with throttles pointing radially (for solar sails).

For solar sails, a reasonable starting point is to thrust radially outward
or along the velocity direction.

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition
- `magnitude::Number`: Throttle magnitude (default: 0.5)

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

"""
    initial_guess_random(problem; seed=nothing)

Generate a random initial guess for throttle vectors.

Each throttle gets a random direction (uniform on unit sphere) and 
random magnitude in [0, 1].

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition
- `seed::Union{Nothing,Int}`: Random seed for reproducibility (default: nothing)

# Returns
- `throttles::Vector{SVector{3}}`: Random throttle vectors with magnitude in [0, 1]
"""
function initial_guess_random(problem::SimsFlanaganProblem{T}; seed::Int = 1234) where {T}
    n_seg = problem.options.n_segments
    throttles = Vector{SVector{3,T}}(undef, n_seg)

    rng = Xoshiro(seed)

    for i = 1:n_seg
        # Random direction on unit sphere
        θ = 2π * rand(rng, T)
        φ = acos(2 * rand(rng, T) - 1)

        # Random magnitude in [0, 1]
        mag = rand(rng, T)

        throttles[i] =
            SVector{3,T}(mag * sin(φ) * cos(θ), mag * sin(φ) * sin(θ), mag * cos(φ))
    end

    return throttles
end

"""
    initial_guess_constant(problem; direction=[1,1,1], magnitude=0.5)

Generate a constant initial guess with all throttles pointing in the same direction.

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition
- `direction::AbstractVector`: Thrust direction (will be normalized)
- `magnitude::Number`: Throttle magnitude (default: 0.5)

# Returns
- `throttles::Vector{SVector{3}}`: Constant throttle vectors
"""
function initial_guess_constant(
    problem::SimsFlanaganProblem{T};
    direction::AbstractVector = [1.0, 1.0, 1.0],
    magnitude::Number = 0.5,
) where {T}
    n_seg = problem.options.n_segments

    # Normalize direction
    dir_norm = norm(direction)
    if dir_norm < eps(T)
        dir = SVector{3,T}(1, 0, 0)
    else
        dir = SVector{3,T}(
            direction[1]/dir_norm,
            direction[2]/dir_norm,
            direction[3]/dir_norm,
        )
    end

    throttle = T(magnitude) * dir

    return [throttle for _ = 1:n_seg]
end

# =============================================================================
# Strategy-based Initial Guess Generation
# =============================================================================

"""
    generate_initial_guess(problem, strategy::AbstractInitialGuess)

Generate initial throttle guess based on the specified strategy.
"""
function generate_initial_guess(problem::SimsFlanaganProblem, ::LambertGuess)
    return initial_guess_lambert(problem)
end

function generate_initial_guess(problem::SimsFlanaganProblem, ::ZeroGuess)
    return initial_guess_zero(problem)
end

function generate_initial_guess(problem::SimsFlanaganProblem, strategy::RandomGuess)
    return initial_guess_random(problem; seed = strategy.seed)
end

function generate_initial_guess(problem::SimsFlanaganProblem, strategy::ConstantGuess)
    return initial_guess_constant(
        problem;
        direction = strategy.direction,
        magnitude = strategy.magnitude,
    )
end

function generate_initial_guess(problem::SimsFlanaganProblem, strategy::RadialGuess)
    return initial_guess_radial(problem; magnitude = strategy.magnitude)
end
