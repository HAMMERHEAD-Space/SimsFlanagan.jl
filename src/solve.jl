export simsflanagan_solve, scaled_mismatch_constraints

#=
Solver functions for Sims-Flanagan transcription.

Uses Optimization.jl with constrained NLP formulation matching PyKEP's approach.
Default solver is MadNLP (interior-point method).

Supports multiple propulsion types: constant thrust, SEP, and solar sails.

Follows SciMLBase interface:
- solve(problem) - solve with default algorithm
- solve(problem, alg) - solve with specified algorithm
- init(problem, alg) - create iterator for advanced usage
- solve!(iterator) - execute solve on iterator
=#

# =============================================================================
# Canonical Scaling (astrodynamics standard)
# =============================================================================

"""
    get_canonical_scales(problem)

Get canonical units for non-dimensionalizing the problem.

Uses standard astrodynamics canonical units:
- LU (Length Unit): characteristic distance (max of initial/final radius)
- TU (Time Unit): √(LU³/μ)
- VU (Velocity Unit): LU/TU = √(μ/LU)
- MU (Mass Unit): initial spacecraft mass

Returns (r_scale, v_scale, m_scale) for normalizing constraints.
"""
function get_canonical_scales(problem::SimsFlanaganProblem)
    # Length unit: characteristic orbit size
    LU = max(norm(problem.r0), norm(problem.rf))

    # Velocity unit: canonical velocity = √(μ/LU)
    VU = sqrt(problem.μ / LU)

    # Mass unit: initial mass
    MU = mass(problem.spacecraft)

    return LU, VU, MU
end

# =============================================================================
# Helper Functions
# =============================================================================

"""
    throttles_to_vector(throttles)

Flatten throttle SVectors to a plain vector for optimization.
"""
function throttles_to_vector(throttles::AbstractVector{<:AbstractVector{T}}) where {T}
    n = length(throttles)
    x = Vector{T}(undef, 3 * n)
    for i = 1:n
        x[3*(i-1)+1] = throttles[i][1]
        x[3*(i-1)+2] = throttles[i][2]
        x[3*(i-1)+3] = throttles[i][3]
    end
    return x
end

"""
    vector_to_throttles(x, n_segments)

Convert flat vector back to throttle SVectors.
"""
function vector_to_throttles(x::AbstractVector{T}, n_segments::Int) where {T}
    throttles = Vector{SVector{3,T}}(undef, n_segments)
    for i = 1:n_segments
        throttles[i] = SVector{3,T}(x[3*(i-1)+1], x[3*(i-1)+2], x[3*(i-1)+3])
    end
    return throttles
end

# =============================================================================
# Constraint Scaling (PyKEP-style)
# =============================================================================

"""
    scaled_mismatch_constraints(problem, throttles)

Compute mismatch constraints with canonical normalization.

Returns normalized constraints using canonical units:
- Position: Δr / LU (length unit)
- Velocity: Δv / VU (velocity unit = √(μ/LU))
- Mass: Δm / MU (mass unit)
"""
function scaled_mismatch_constraints(
    problem::SimsFlanaganProblem,
    throttles::AbstractVector{<:AbstractVector},
)
    mismatch = compute_mismatch(problem, throttles)

    # Scale using canonical units
    LU, VU, MU = get_canonical_scales(problem)
    T = eltype(mismatch)

    # Return as SVector for type stability
    return SVector{7,T}(
        mismatch[1] / LU,
        mismatch[2] / LU,
        mismatch[3] / LU,
        mismatch[4] / VU,
        mismatch[5] / VU,
        mismatch[6] / VU,
        mismatch[7] / MU,
    )
end

"""
    throttle_magnitude_constraints(throttles)

Compute throttle magnitude inequality constraints.

Returns n_segments values, each being (|throttle|² - 1).
Constraint is satisfied when ≤ 0.
"""
function throttle_magnitude_constraints(
    throttles::AbstractVector{<:AbstractVector{T}},
) where {T}
    n = length(throttles)
    constraints = Vector{T}(undef, n)
    for i = 1:n
        # |throttle|² - 1 ≤ 0
        constraints[i] = dot(throttles[i], throttles[i]) - one(T)
    end
    return constraints
end

# =============================================================================
# SciMLBase Interface: solve
# =============================================================================

"""
    solve(prob::SimsFlanaganProblem; kwargs...)

Solve a Sims-Flanagan trajectory optimization problem using MadNLP.

# Keyword Arguments
- `initial_guess_strategy::AbstractInitialGuess`: Strategy for generating initial guess
  - `RandomGuess(; seed=1234)` (default): Use random throttles
  - `LambertGuess()`: Use Lambert arc solution
  - `ZeroGuess()`: Use zero throttles
  - `ConstantGuess(; direction, magnitude)`: Use constant direction
  - `RadialGuess(; magnitude)`: Use radial direction (for solar sails)
- `initial_guess`: Direct throttle vectors (overrides strategy if provided)
- `max_iter::Int`: Maximum iterations (default: from problem.options.max_iter)
- `tol::Float64`: Convergence tolerance (default: from problem.options.tol)

# Returns
- `SimsFlanaganSolution`: The optimization result with `retcode` from SciMLBase
"""
function SciMLBase.solve(
    prob::SimsFlanaganProblem;
    initial_guess::Union{Nothing,AbstractVector{<:AbstractVector}} = nothing,
    initial_guess_strategy::AbstractInitialGuess = RandomGuess(),
    max_iter::Int = prob.options.max_iter,
    tol::Float64 = prob.options.tol,
)
    return _solve_internal(
        prob;
        initial_guess = initial_guess,
        initial_guess_strategy = initial_guess_strategy,
        max_iter = max_iter,
        tol = tol,
    )
end

# =============================================================================
# Legacy API (for backwards compatibility)
# =============================================================================

"""
    simsflanagan_solve(problem; kwargs...)

Legacy solver interface. Prefer `solve(problem)` for SciMLBase compatibility.

See `solve` for full documentation.
"""
function simsflanagan_solve(
    problem::SimsFlanaganProblem;
    initial_guess::Union{Nothing,AbstractVector{<:AbstractVector}} = nothing,
    initial_guess_strategy::AbstractInitialGuess = RandomGuess(),
    maxiters::Int = problem.options.max_iter,
)
    return SciMLBase.solve(
        problem;
        initial_guess = initial_guess,
        initial_guess_strategy = initial_guess_strategy,
        max_iter = maxiters,
    )
end

# =============================================================================
# Internal Solver Implementation
# =============================================================================

"""
Internal solve implementation using MadNLP.
"""
function _solve_internal(
    problem::SimsFlanaganProblem;
    initial_guess::Union{Nothing,AbstractVector{<:AbstractVector}} = nothing,
    initial_guess_strategy::AbstractInitialGuess = RandomGuess(),
    max_iter::Int = problem.options.max_iter,
    tol::Float64 = problem.options.tol,
)
    n_seg = problem.options.n_segments
    n_fwd = problem.options.n_fwd
    n_bwd = n_seg - n_fwd
    n_vars = 3 * n_seg
    sundman_c = problem.options.sundman_c

    # Get initial guess
    throttles0 = if initial_guess !== nothing
        initial_guess
    else
        generate_initial_guess(problem, initial_guess_strategy)
    end

    x0 = throttles_to_vector(throttles0)

    # Pre-compute segment times (same as used in constraints for consistency)
    Δt_fwd, Δt_bwd = compute_sundman_leg_times(
        problem.r0,
        problem.rf,
        problem.tof,
        n_fwd,
        n_bwd;
        c = sundman_c,
    )
    Δt_segments = vcat(Δt_fwd, Δt_bwd)

    # Objective: minimize total ΔV
    function objective(x, p)
        throttles = vector_to_throttles(x, n_seg)
        return compute_total_Δv(throttles, Δt_segments, problem.spacecraft)
    end

    # Constraints: [7 equality (mismatch), n_seg inequality (throttle)]
    n_eq = 7
    n_ineq = n_seg

    function constraints!(res, x, p)
        throttles = vector_to_throttles(x, n_seg)

        # Equality constraints: scaled mismatch = 0
        mismatch = scaled_mismatch_constraints(problem, throttles)
        for i = 1:7
            res[i] = mismatch[i]
        end

        # Inequality constraints: |throttle|² - 1 ≤ 0
        for i = 1:n_seg
            res[7+i] = dot(throttles[i], throttles[i]) - 1.0
        end

        return nothing
    end

    # Bounds: throttle components in [-1, 1]
    lb = fill(-1.0, n_vars)
    ub = fill(1.0, n_vars)

    # Constraint bounds: equality = 0, inequality ≤ 0
    lcons = vcat(zeros(n_eq), fill(-Inf, n_ineq))
    ucons = vcat(zeros(n_eq), zeros(n_ineq))

    # Build Optimization.jl problem with constraints
    opt_f = Optimization.OptimizationFunction(
        objective,
        Optimization.AutoForwardDiff();
        cons = constraints!,
    )

    opt_prob = Optimization.OptimizationProblem(
        opt_f,
        x0,
        nothing;
        lb = lb,
        ub = ub,
        lcons = lcons,
        ucons = ucons,
    )

    # Create MadNLP optimizer
    optimizer = MadNLP.Optimizer(; linear_solver = MadNLPMumps.MumpsSolver)

    # Solve
    opt_sol = Optimization.solve(
        opt_prob,
        optimizer;
        max_iter = max_iter,
        print_level = problem.options.verbosity > 0 ? MadNLP.INFO : MadNLP.ERROR,
        tol = tol,
        nlp_scaling = false,
    )

    return _build_solution(problem, opt_sol, tol)
end

# =============================================================================
# Solution Building
# =============================================================================

"""
Build SimsFlanaganSolution from optimization result.
"""
function _build_solution(problem::SimsFlanaganProblem, opt_sol, tol::Float64)
    n_seg = problem.options.n_segments
    n_fwd = problem.options.n_fwd
    n_bwd = n_seg - n_fwd
    sundman_c = problem.options.sundman_c

    x_opt = opt_sol.u
    throttles = vector_to_throttles(x_opt, n_seg)

    # Compute segment times (consistent with constraints)
    Δt_fwd, Δt_bwd = compute_sundman_leg_times(
        problem.r0,
        problem.rf,
        problem.tof,
        n_fwd,
        n_bwd;
        c = sundman_c,
    )
    Δt_segments = vcat(Δt_fwd, Δt_bwd)

    # Compute final solution metrics
    mismatch = compute_mismatch(problem, throttles)
    Δv_total = compute_total_Δv(throttles, Δt_segments, problem.spacecraft)
    masses = compute_segment_masses(problem, throttles, Δt_segments)

    # Check convergence based on scaled constraints
    scaled_mismatch = scaled_mismatch_constraints(problem, throttles)
    scaled_norm = norm(scaled_mismatch)

    # Get retcode directly from optimization result
    retcode = opt_sol.retcode

    # Try to extract iteration count
    iterations = 0
    if opt_sol.original !== nothing && hasproperty(opt_sol.original, :iterations)
        iterations = opt_sol.original.iterations
    end

    if problem.options.verbosity > 0
        Δr_norm = norm(SVector{3}(mismatch[1], mismatch[2], mismatch[3]))
        Δv_norm = norm(SVector{3}(mismatch[4], mismatch[5], mismatch[6]))
        Δm = mismatch[7]
        @info "SimsFlanagan solve complete" retcode Δv_total Δr_norm Δv_norm Δm scaled_norm iterations
    end

    return SimsFlanaganSolution(
        problem,
        throttles,
        masses,
        mismatch,
        Δv_total,
        retcode,
        iterations,
    )
end

# =============================================================================
# Mass Computation
# =============================================================================

"""
    compute_segment_masses(problem, throttles, Δt_segments=nothing)

Compute the mass at each segment boundary during forward propagation.

For solar sails, mass remains constant throughout the trajectory.

# Arguments
- `problem::SimsFlanaganProblem`: The problem definition
- `throttles::AbstractVector`: Throttle vectors for all segments
- `Δt_segments::Union{Nothing, AbstractVector}`: Duration of each segment [s] (uses equal if nothing)
"""
function compute_segment_masses(
    problem::SimsFlanaganProblem,
    throttles::AbstractVector{<:AbstractVector},
    Δt_segments::Union{Nothing,AbstractVector} = nothing,
)
    n_seg = problem.options.n_segments
    spacecraft = problem.spacecraft

    T = typeof(mass(spacecraft))
    masses = Vector{T}(undef, n_seg + 1)
    masses[1] = mass(spacecraft)

    # Solar sails have constant mass
    if !has_propellant_consumption(spacecraft)
        for i = 1:n_seg
            masses[i+1] = masses[1]
        end
        return masses
    end

    # For propellant-consuming systems
    vex = exhaust_velocity(spacecraft)
    ref_thrust = get_reference_thrust(spacecraft)

    dry_mass = spacecraft.dry_mass
    for i = 1:n_seg
        # Get segment duration
        Δt_seg = Δt_segments === nothing ? problem.tof / n_seg : Δt_segments[i]

        throttle_mag = clamp(safe_norm(throttles[i]), zero(T), one(T))
        thrust = ref_thrust * throttle_mag
        mdot = thrust / vex
        Δm = mdot * Δt_seg
        # Can only consume propellant, not dry mass
        available_propellant = masses[i] - dry_mass
        Δm = min(Δm, max(available_propellant, zero(T)))
        masses[i+1] = masses[i] - Δm
    end

    return masses
end
