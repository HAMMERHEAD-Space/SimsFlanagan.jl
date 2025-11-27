#=
Core types for Sims-Flanagan low-thrust trajectory optimization.

The Sims-Flanagan method discretizes a low-thrust trajectory into segments,
applying impulsive ΔV at segment midpoints to approximate continuous thrust.
The trajectory is propagated forward from the initial state and backward from
the final state, meeting at a match point where continuity is enforced.
=#

export Spacecraft, SimsFlanaganOptions, SimsFlanaganProblem, SimsFlanaganSolution

"""
    Spacecraft{MT, TT, IT}

Spacecraft parameters for the Sims-Flanagan transcription.

# Fields
- `mass::MT`: Initial spacecraft mass [kg]
- `thrust::TT`: Maximum thrust magnitude [N]
- `isp::IT`: Specific impulse [s]

# Notes
The exhaust velocity is computed as `vex = isp * g0` where `g0 = 9.80665 m/s²`.
"""
struct Spacecraft{MT<:Number,TT<:Number,IT<:Number}
    mass::MT
    thrust::TT
    isp::IT
end

"""
    exhaust_velocity(sc::Spacecraft)

Compute the exhaust velocity [m/s] from specific impulse.
"""
function exhaust_velocity(sc::Spacecraft; g0::Number = 9.80665)
    return sc.isp * g0
end

"""
    SimsFlanaganOptions

Options for the Sims-Flanagan solver.

# Fields
- `n_segments::Int`: Number of thrust segments (default: 10)
- `n_fwd::Int`: Number of forward propagation segments (default: n_segments ÷ 2)
- `tol::Float64`: Convergence tolerance for mismatch constraints (default: 1e-8)
- `max_iter::Int`: Maximum number of optimizer iterations (default: 1000)
- `verbosity::Int`: Output verbosity level (0=silent, 1=summary, 2=detailed)
"""
struct SimsFlanaganOptions
    n_segments::Int
    n_fwd::Int
    tol::Float64
    max_iter::Int
    verbosity::Int
    cutoff_threshold::Number

    function SimsFlanaganOptions(;
        n_segments::Int = 10,
        n_fwd::Int = n_segments ÷ 2,
        tol::Float64 = 1e-8,
        max_iter::Int = 1000,
        verbosity::Int = 1,
        cutoff_threshold::Number = 10.0*eps(Float64),
    )
        n_segments > 0 || throw(ArgumentError("n_segments must be positive"))
        0 < n_fwd <= n_segments || throw(ArgumentError("n_fwd must be in (0, n_segments]"))
        tol > 0 || throw(ArgumentError("tol must be positive"))
        max_iter > 0 || throw(ArgumentError("max_iter must be positive"))
        cutoff_threshold > 0 || throw(ArgumentError("cutoff_threshold must be positive"))
        new(n_segments, n_fwd, tol, max_iter, verbosity, cutoff_threshold)
    end
end

"""
    SimsFlanaganProblem{T}

A Sims-Flanagan low-thrust trajectory optimization problem.

# Fields
- `r0::SVector{3,T}`: Initial position [km]
- `v0::SVector{3,T}`: Initial velocity [km/s]
- `rf::SVector{3,T}`: Final position [km]
- `vf::SVector{3,T}`: Final velocity [km/s]
- `tof::T`: Time of flight [s]
- `μ::T`: Gravitational parameter [km³/s²]
- `spacecraft::Spacecraft{T}`: Spacecraft parameters
- `options::SimsFlanaganOptions`: Solver options

# Notes
The problem is parameterized by throttle vectors at each segment midpoint.
Each throttle vector has magnitude in [0, 1] representing the fraction of 
maximum thrust, and direction representing the thrust pointing.
"""
struct SimsFlanaganProblem{
    RT1<:Number,
    VT1<:Number,
    RT2<:Number,
    VT2<:Number,
    TT<:Number,
    MT<:Number,
}
    r0::AbstractVector{RT1}
    v0::AbstractVector{VT1}
    rf::AbstractVector{RT2}
    vf::AbstractVector{VT2}
    tof::TT
    μ::MT
    spacecraft::Spacecraft
    options::SimsFlanaganOptions
end

"""
    SimsFlanaganSolution

Solution to a Sims-Flanagan problem.

# Fields
- `problem::SimsFlanaganProblem`: The original problem
- `throttles::AbstractVector{AbstractVector}`: Throttle vectors at each segment [dimensionless]
- `masses::AbstractVector`: Mass at each segment boundary [kg]
- `mismatch::AbstractVector`: Final mismatch at match point [Δr; Δv; Δm]
- `Δv_total::Number`: Total ΔV expended [km/s]
- `converged::Bool`: Whether the solver converged
- `iterations::Int`: Number of iterations taken
"""
struct SimsFlanaganSolution{
    P<:SimsFlanaganProblem,
    TH<:AbstractVector,
    M<:AbstractVector,
    MM<:AbstractVector,
    DV<:Number,
}
    problem::P
    throttles::TH
    masses::M
    mismatch::MM
    Δv_total::DV
    converged::Bool
    iterations::Int
end

# Pretty printing
function Base.show(io::IO, sc::Spacecraft)
    print(io, "Spacecraft(mass=$(sc.mass) kg, thrust=$(sc.thrust) N, Isp=$(sc.isp) s)")
end

function Base.show(io::IO, opts::SimsFlanaganOptions)
    print(io, "SimsFlanaganOptions(n_segments=$(opts.n_segments), n_fwd=$(opts.n_fwd))")
end

function Base.show(io::IO, prob::SimsFlanaganProblem)
    print(
        io,
        "SimsFlanaganProblem(tof=$(prob.tof) s, μ=$(prob.μ) km³/s², $(prob.options.n_segments) segments)",
    )
end

function Base.show(io::IO, sol::SimsFlanaganSolution)
    status = sol.converged ? "converged" : "not converged"
    mismatch_norm = norm(sol.mismatch)
    print(
        io,
        "SimsFlanaganSolution($(status), Δv=$(sol.Δv_total) km/s, mismatch=$(mismatch_norm))",
    )
end
