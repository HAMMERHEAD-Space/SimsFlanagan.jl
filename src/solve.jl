export simsflanagan_solve

#=
Solver functions for Sims-Flanagan transcription.

Uses Optimization.jl with OptimizationOptimJL for the NLP solve.
=#
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

"""
    simsflanagan_solve(problem; kwargs...)

Solve a Sims-Flanagan trajectory optimization problem.

# Arguments
- `problem::SimsFlanaganProblem`: The problem to solve

# Keyword Arguments
- `alg`: Optimization algorithm (default: `Optim.LBFGS()` with AutoForwardDiff)
- `initial_guess::Union{Nothing, Vector{SVector{3}}}`: Initial throttle guess
- `use_lambert_guess::Bool=true`: Use Lambert solution for initial guess
- `penalty_weight::Float64=1e6`: Penalty weight for mismatch constraint
- `maxiters::Int`: Maximum iterations (default: from problem options)

# Returns
- `SimsFlanaganSolution`: The optimization result
"""
function simsflanagan_solve(
    problem::SimsFlanaganProblem;
    alg = Optim.LBFGS(),  # Gradient-based default with AD
    initial_guess::Union{Nothing,AbstractVector{<:AbstractVector}} = nothing,
    use_lambert_guess::Val{lambert} = Val(true),
    penalty_weight::Float64 = 1e6,
    maxiters::Int = problem.options.max_iter,
) where {lambert}

    n_seg = problem.options.n_segments

    # Get initial guess
    if initial_guess === nothing
        if lambert
            throttles0 = initial_guess_lambert(problem)
        else
            throttles0 = initial_guess_zero(problem)
        end
    else
        throttles0 = initial_guess
    end

    x0 = throttles_to_vector(throttles0)

    # Objective function: minimize ΔV + penalty for constraint violation
    function objective(x, p)
        throttles = vector_to_throttles(x, n_seg)
        Δt_seg = problem.tof / n_seg

        # Total ΔV (to minimize)
        Δv = compute_total_Δv(throttles, Δt_seg, problem.spacecraft)

        # Mismatch constraint penalty (match point continuity)
        mismatch = compute_mismatch(problem, throttles)
        mismatch_penalty = penalty_weight * dot(mismatch, mismatch)

        # Throttle magnitude penalty (should be ≤ 1)
        throttle_penalty = zero(eltype(x))
        for t in throttles
            mag = norm(t)
            if mag > 1
                throttle_penalty += 100.0 * (mag - 1)^2
            end
        end

        return Δv + mismatch_penalty + throttle_penalty
    end

    # Build and solve Optimization.jl problem with automatic differentiation
    opt_f = OptimizationFunction(objective, AutoForwardDiff())
    opt_prob = OptimizationProblem(opt_f, x0, nothing)
    opt_sol = Optimization.solve(opt_prob, alg; maxiters = maxiters)

    # Extract solution
    x_opt = opt_sol.u
    throttles = vector_to_throttles(x_opt, n_seg)

    # Compute final solution metrics
    Δt_seg = problem.tof / n_seg
    mismatch = compute_mismatch(problem, throttles)
    Δv_total = compute_total_Δv(throttles, Δt_seg, problem.spacecraft)
    masses = compute_segment_masses(problem, throttles)

    mismatch_norm = norm(mismatch)
    converged = mismatch_norm < problem.options.tol

    # Try to extract iteration count
    iterations = 0
    if opt_sol.original !== nothing && hasproperty(opt_sol.original, :iterations)
        iterations = opt_sol.original.iterations
    end

    if problem.options.verbosity > 0
        status = converged ? "converged" : "not converged"
        @info "SimsFlanagan solve complete" status Δv_total mismatch_norm iterations
    end

    return SimsFlanaganSolution(
        problem,
        throttles,
        masses,
        mismatch,
        Δv_total,
        converged,
        iterations,
    )
end

"""
    compute_segment_masses(problem, throttles)

Compute the mass at each segment boundary during forward propagation.
"""
function compute_segment_masses(
    problem::SimsFlanaganProblem,
    throttles::AbstractVector{<:AbstractVector},
)
    n_seg = problem.options.n_segments
    Δt_seg = problem.tof / n_seg

    T = typeof(problem.spacecraft.mass)
    masses = Vector{T}(undef, n_seg + 1)
    masses[1] = problem.spacecraft.mass

    vex = exhaust_velocity(problem.spacecraft)

    for i = 1:n_seg
        throttle_mag = min(norm(throttles[i]), one(T))
        thrust = problem.spacecraft.thrust * throttle_mag
        mdot = thrust / vex
        Δm = mdot * Δt_seg
        Δm = min(Δm, masses[i])
        masses[i+1] = masses[i] - Δm
    end

    return masses
end
