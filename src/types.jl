#=
Core types for Sims-Flanagan low-thrust trajectory optimization.

The Sims-Flanagan method discretizes a low-thrust trajectory into segments,
applying impulsive ΔV at segment midpoints to approximate continuous thrust.
The trajectory is propagated forward from the initial state and backward from
the final state, meeting at a match point where continuity is enforced.
=#

export AbstractSpacecraft, Spacecraft, SEPSpacecraft, SolarSail
export SimsFlanaganOptions, SimsFlanaganProblem, SimsFlanaganSolution
export exhaust_velocity, characteristic_acceleration
export AbstractInitialGuess,
    LambertGuess, ZeroGuess, RandomGuess, ConstantGuess, RadialGuess
export mass

# SciMLBase exports
export solve, remake

# =============================================================================
# Abstract Problem Type (SciMLBase Interface)
# =============================================================================

"""
    AbstractSimsFlanaganProblem

Abstract base type for Sims-Flanagan trajectory optimization problems.
Inherits from SciMLBase.AbstractSciMLProblem for compatibility with SciML ecosystem.
"""
abstract type AbstractSimsFlanaganProblem <: SciMLBase.AbstractSciMLProblem end

# =============================================================================
# Abstract Spacecraft Type
# =============================================================================

"""
    AbstractSpacecraft

Abstract base type for all spacecraft propulsion systems.

Concrete subtypes must implement:
- `compute_thrust(sc, r, throttle_mag)` - Compute thrust magnitude at position r
- `compute_mass_flow(sc, thrust)` - Compute mass flow rate for given thrust
- `has_propellant_consumption(sc)` - Whether the system consumes propellant
"""
abstract type AbstractSpacecraft end

# =============================================================================
# Constant Thrust Spacecraft (Chemical/Electric Propulsion)
# =============================================================================

"""Total initial mass [kg]"""
mass(sc::AbstractSpacecraft) = sc.dry_mass + sc.wet_mass

"""
    Spacecraft{MT, TT, IT}

Constant-thrust spacecraft (chemical or electric propulsion at constant power).

# Fields
- `dry_mass::MT`: Dry mass (without propellant) [kg]
- `wet_mass::MT`: Propellant mass [kg]
- `thrust::TT`: Maximum thrust magnitude [N]
- `isp::IT`: Specific impulse [s]

# Computed Properties
- `mass`: Total initial mass = dry_mass + wet_mass

# Notes
The exhaust velocity is computed as `vex = isp * g0` where `g0 = 9.80665 m/s²`.
Thrust is independent of distance from the Sun.

# Constructors
```julia
# Full specification
Spacecraft(dry_mass, wet_mass, thrust, isp)
```
"""
struct Spacecraft{MT<:Number,TT<:Number,IT<:Number} <: AbstractSpacecraft
    dry_mass::MT
    wet_mass::MT
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
    compute_thrust(sc::Spacecraft, r, throttle_mag)

Compute thrust magnitude [N] for constant-thrust spacecraft.
Thrust is independent of position.
"""
function compute_thrust(sc::Spacecraft, r::AbstractVector, throttle_mag::Number)
    return sc.thrust * min(throttle_mag, one(throttle_mag))
end

"""
    compute_mass_flow(sc::Spacecraft, thrust)

Compute mass flow rate [kg/s] for given thrust.
"""
function compute_mass_flow(sc::Spacecraft, thrust::Number)
    vex = exhaust_velocity(sc)
    return thrust / vex
end

"""
    has_propellant_consumption(::Spacecraft)

Constant-thrust spacecraft consume propellant.
"""
has_propellant_consumption(::Spacecraft) = true

# =============================================================================
# Solar Electric Propulsion (SEP) Spacecraft
# =============================================================================

"""
    SEPSpacecraft{MT, TT, IT, RT}

Solar Electric Propulsion spacecraft where thrust depends on solar distance.

Power available scales as (r_ref/r)² where r is the heliocentric distance.
Thrust is proportional to available power.

# Fields
- `dry_mass::MT`: Dry mass (without propellant) [kg]
- `wet_mass::MT`: Propellant mass [kg]
- `thrust_ref::TT`: Maximum thrust at reference distance [N]
- `isp::IT`: Specific impulse [s]
- `r_ref::RT`: Reference distance [km] (default: 1 AU = 1.495978707e8 km)

# Computed Properties
- `mass`: Total initial mass = dry_mass + wet_mass

# Notes
At distance r from the Sun:
  thrust(r) = thrust_ref × (r_ref / r)²

This models the inverse-square falloff of solar power with distance.
"""
struct SEPSpacecraft{MT<:Number,TT<:Number,IT<:Number,RT<:Number} <: AbstractSpacecraft
    dry_mass::MT
    wet_mass::MT
    thrust_ref::TT
    isp::IT
    r_ref::RT

    function SEPSpacecraft(
        dry_mass::MT,
        wet_mass::MT,
        thrust_ref::TT,
        isp::IT,
        r_ref::RT,
    ) where {MT<:Number,TT<:Number,IT<:Number,RT<:Number}
        dry_mass > 0 || throw(ArgumentError("dry_mass must be positive"))
        wet_mass >= 0 || throw(ArgumentError("wet_mass must be non-negative"))
        thrust_ref >= 0 || throw(ArgumentError("thrust_ref must be non-negative"))
        isp > 0 || throw(ArgumentError("isp must be positive"))
        r_ref > 0 || throw(ArgumentError("r_ref must be positive"))
        new{MT,TT,IT,RT}(dry_mass, wet_mass, thrust_ref, isp, r_ref)
    end
end

"""
    exhaust_velocity(sc::SEPSpacecraft)

Compute the exhaust velocity [m/s] from specific impulse.
"""
function exhaust_velocity(sc::SEPSpacecraft; g0::Number = 9.80665)
    return sc.isp * g0
end

"""
    compute_thrust(sc::SEPSpacecraft, r, throttle_mag)

Compute thrust magnitude [N] for SEP spacecraft at position r.
Thrust scales as (r_ref/|r|)² due to solar power availability.
"""
function compute_thrust(sc::SEPSpacecraft, r::AbstractVector, throttle_mag::Number)
    r_mag = norm(r)
    power_factor = (sc.r_ref / r_mag)^2
    return sc.thrust_ref * power_factor * min(throttle_mag, one(throttle_mag))
end

"""
    compute_mass_flow(sc::SEPSpacecraft, thrust)

Compute mass flow rate [kg/s] for given thrust.
"""
function compute_mass_flow(sc::SEPSpacecraft, thrust::Number)
    vex = exhaust_velocity(sc)
    return thrust / vex
end

"""
    has_propellant_consumption(::SEPSpacecraft)

SEP spacecraft consume propellant.
"""
has_propellant_consumption(::SEPSpacecraft) = true

# =============================================================================
# Solar Sail
# =============================================================================

"""
    SolarSail{MT, AT, CT, RT}

Solar sail spacecraft using solar radiation pressure for propulsion.

# Fields
- `mass::MT`: Spacecraft mass [kg] (constant - no propellant consumption)
- `area::AT`: Sail area [m²]
- `reflectivity::CT`: Sail reflectivity coefficient (0 to 1, typically ~0.9)
- `r_ref::RT`: Reference distance [km] (default: 1 AU = 1.495978707e8 km)

# Notes
The characteristic acceleration at reference distance is:
  a_c = (1 + η) × P_sun × A / (c × m)

where:
- η = reflectivity (0 for absorption, 1 for perfect reflection)
- P_sun ≈ 1361 W/m² (solar constant at 1 AU)
- c = 299792458 m/s (speed of light)

At distance r from the Sun:
  a(r) = a_c × (r_ref / r)²

The thrust direction is controlled by the sail orientation (cone angle).
For this implementation, the throttle vector magnitude represents cos²(α) 
where α is the cone angle, and direction is the thrust direction.
"""
struct SolarSail{MT<:Number,AT<:Number,CT<:Number,RT<:Number} <: AbstractSpacecraft
    dry_mass::MT
    wet_mass::MT
    area::AT
    reflectivity::CT
    r_ref::RT

    function SolarSail(
        dry_mass::MT,
        area::AT,
        reflectivity::CT = 0.9,
        r_ref::RT = 1.495978707e8;  # 1 AU in km
        wet_mass::MT = 0.0,
    ) where {MT<:Number,AT<:Number,CT<:Number,RT<:Number}
        dry_mass > 0 || throw(ArgumentError("dry_mass must be positive"))
        area > 0 || throw(ArgumentError("area must be positive"))
        0 <= reflectivity <= 1 || throw(ArgumentError("reflectivity must be in [0, 1]"))
        r_ref > 0 || throw(ArgumentError("r_ref must be positive"))
        new{MT,AT,CT,RT}(dry_mass, wet_mass, area, reflectivity, r_ref)
    end
end

# Solar constant at 1 AU [W/m²]
const SOLAR_CONSTANT = 1361.0

# Speed of light [m/s]
const SPEED_OF_LIGHT = 299792458.0

"""
    characteristic_acceleration(sail::SolarSail)

Compute the characteristic acceleration [m/s²] at the reference distance.

This is the maximum acceleration achievable with the sail perpendicular to the Sun.
"""
function characteristic_acceleration(sail::SolarSail)
    # a_c = (1 + η) × P / (c × m) × A
    # Units: [W/m² × m²] / [m/s × kg] = [W] / [m/s × kg] = [kg⋅m²/s³] / [m/s × kg] = [m/s²]
    return (1 + sail.reflectivity) * SOLAR_CONSTANT * sail.area /
           (SPEED_OF_LIGHT * mass(sail))
end

"""
    compute_thrust(sail::SolarSail, r, throttle_mag)

Compute effective thrust magnitude [N] for solar sail at position r.

The throttle_mag represents the efficiency factor (0 to 1) which depends on 
sail orientation. For an ideal flat sail: efficiency = cos²(α) where α is 
the cone angle between sail normal and sun direction.

Returns thrust in Newtons (force = mass × acceleration).
"""
function compute_thrust(sail::SolarSail, r::AbstractVector, throttle_mag::Number)
    r_mag = norm(r)
    distance_factor = (sail.r_ref / r_mag)^2
    a_c = characteristic_acceleration(sail)  # [m/s²]
    efficiency = min(throttle_mag, one(throttle_mag))
    acceleration = a_c * distance_factor * efficiency  # [m/s²]
    return mass(sail) * acceleration  # [N]
end

"""
    compute_mass_flow(::SolarSail, thrust)

Solar sails have zero mass flow (no propellant consumption).
"""
function compute_mass_flow(::SolarSail, thrust::Number)
    return zero(thrust)
end

"""
    has_propellant_consumption(::SolarSail)

Solar sails do not consume propellant.
"""
has_propellant_consumption(::SolarSail) = false

"""
    exhaust_velocity(::SolarSail)

Solar sails have no exhaust - return Inf for compatibility.
"""
function exhaust_velocity(::SolarSail; g0::Number = 9.80665)
    return Inf
end

# =============================================================================
# Initial Guess Strategies
# =============================================================================

"""
    AbstractInitialGuess

Abstract type for initial guess strategies in Sims-Flanagan optimization.
"""
abstract type AbstractInitialGuess end

"""
    LambertGuess()

Use Lambert arc solution to generate initial throttle guess.
Falls back to zero guess if Lambert fails (e.g., zero transfer angle).
"""
struct LambertGuess <: AbstractInitialGuess end

"""
    ZeroGuess()

Use zero throttles as initial guess.
"""
struct ZeroGuess <: AbstractInitialGuess end

"""
    RandomGuess(; seed=1234)

Use random throttles as initial guess.
Each throttle gets random direction and magnitude in [0, 1].

# Fields
- `seed::Int`: Random seed for reproducibility (default: 1234)
"""
struct RandomGuess <: AbstractInitialGuess
    seed::Int
    RandomGuess(; seed::Int = 1234) = new(seed)
end

"""
    ConstantGuess(; direction=[1,1,1], magnitude=0.5)

Use constant throttle direction for all segments.

# Fields
- `direction::Vector{Float64}`: Thrust direction (will be normalized)
- `magnitude::Float64`: Throttle magnitude (default: 0.5)
"""
struct ConstantGuess <: AbstractInitialGuess
    direction::Vector{Float64}
    magnitude::Float64
    ConstantGuess(; direction::AbstractVector = [1.0, 1.0, 1.0], magnitude::Float64 = 0.5) =
        new(collect(Float64, direction), magnitude)
end

"""
    RadialGuess(; magnitude=0.5)

Use radial (Sun-pointing) direction for all segments.
Useful for solar sail problems.

# Fields
- `magnitude::Float64`: Throttle magnitude (default: 0.5)
"""
struct RadialGuess <: AbstractInitialGuess
    magnitude::Float64
    RadialGuess(; magnitude::Float64 = 0.5) = new(magnitude)
end

# =============================================================================
# Solver Options
# =============================================================================

"""
    SimsFlanaganOptions

Options for the Sims-Flanagan solver.

# Fields
- `n_segments::Int`: Number of thrust segments (default: 10)
- `n_fwd::Int`: Number of forward propagation segments (default: n_segments ÷ 2)
- `tol::Float64`: Convergence tolerance for mismatch constraints (default: 1e-6)
- `max_iter::Int`: Maximum number of optimizer iterations (default: 1000)
- `verbosity::Int`: Output verbosity level (0=silent, 1=summary, 2=detailed)
- `sundman_c::Float64`: Sundman transformation exponent (default: 0.0)
  - 0.0: Equal time segments (no transformation)
  - 1.0: Classic Sundman (dt ∝ r), more resolution near central body
  - 1.5: Related to eccentric anomaly
"""
struct SimsFlanaganOptions
    n_segments::Int
    n_fwd::Int
    tol::Float64
    max_iter::Int
    verbosity::Int
    sundman_c::Float64

    function SimsFlanaganOptions(;
        n_segments::Int = 10,
        n_fwd::Int = n_segments ÷ 2,
        tol::Float64 = 1e-6,
        max_iter::Int = 1000,
        verbosity::Int = 1,
        sundman_c::Float64 = 0.0,  # Default to no Sundman (equal segments) - enable after debugging``
    )
        n_segments > 0 || throw(ArgumentError("n_segments must be positive"))
        0 < n_fwd <= n_segments || throw(ArgumentError("n_fwd must be in (0, n_segments]"))
        tol > 0 || throw(ArgumentError("tol must be positive"))
        max_iter > 0 || throw(ArgumentError("max_iter must be positive"))
        sundman_c >= 0 || throw(ArgumentError("sundman_c must be non-negative"))
        new(n_segments, n_fwd, tol, max_iter, verbosity, sundman_c)
    end
end

"""
    SimsFlanaganProblem{T}

A Sims-Flanagan low-thrust trajectory optimization problem.
Inherits from AbstractSimsFlanaganProblem for SciMLBase compatibility.

# Fields
- `r0::SVector{3,T}`: Initial position [km]
- `v0::SVector{3,T}`: Initial velocity [km/s]
- `rf::SVector{3,T}`: Final position [km]
- `vf::SVector{3,T}`: Final velocity [km/s]
- `tof::T`: Time of flight [s]
- `μ::T`: Gravitational parameter [km³/s²]
- `spacecraft::AbstractSpacecraft`: Spacecraft/propulsion system
- `options::SimsFlanaganOptions`: Solver options

# Notes
The problem is parameterized by throttle vectors at each segment midpoint.
Each throttle vector has magnitude in [0, 1] representing the fraction of 
maximum thrust (or sail efficiency), and direction representing the thrust pointing.
"""
struct SimsFlanaganProblem{
    RT1<:Number,
    VT1<:Number,
    RT2<:Number,
    VT2<:Number,
    TT<:Number,
    MT<:Number,
    SC<:AbstractSpacecraft,
} <: AbstractSimsFlanaganProblem
    r0::AbstractVector{RT1}
    v0::AbstractVector{VT1}
    rf::AbstractVector{RT2}
    vf::AbstractVector{VT2}
    tof::TT
    μ::MT
    spacecraft::SC
    options::SimsFlanaganOptions
end

"""
    SimsFlanaganSolution

Solution to a Sims-Flanagan problem.
Follows SciMLBase conventions with `retcode` for solution status.

# Fields
- `problem::SimsFlanaganProblem`: The original problem
- `throttles::AbstractVector{AbstractVector}`: Throttle vectors at each segment [dimensionless]
- `masses::AbstractVector`: Mass at each segment boundary [kg]
- `mismatch::AbstractVector`: Full mismatch vector at match point [Δr; Δv; Δm]
- `Δv_total::Number`: Total ΔV expended [km/s] (or effective ΔV for sails)
- `retcode`: Return code from SciMLBase (e.g., `ReturnCode.Success`)
- `iterations::Int`: Number of iterations taken

# Accessors
- `position_mismatch(sol)`: Position mismatch [km] as SVector{3}
- `velocity_mismatch(sol)`: Velocity mismatch [km/s] as SVector{3}
- `mass_mismatch(sol)`: Mass mismatch [kg]
- `position_mismatch_norm(sol)`: ||Δr|| [km]
- `velocity_mismatch_norm(sol)`: ||Δv|| [km/s]
- `converged`: true if `retcode == ReturnCode.Success`
"""
struct SimsFlanaganSolution{
    P<:SimsFlanaganProblem,
    TH<:AbstractVector,
    M<:AbstractVector,
    MM<:AbstractVector,
    DV<:Number,
    RC,
}
    problem::P
    throttles::TH
    masses::M
    mismatch::MM
    Δv_total::DV
    retcode::RC
    iterations::Int
end

# Convenience accessor for backwards compatibility
Base.propertynames(::SimsFlanaganSolution) =
    (:problem, :throttles, :masses, :mismatch, :Δv_total, :retcode, :iterations, :converged)

function Base.getproperty(sol::SimsFlanaganSolution, s::Symbol)
    if s === :converged
        rc = getfield(sol, :retcode)
        return rc == SciMLBase.ReturnCode.Success
    else
        return getfield(sol, s)
    end
end

# Mismatch component accessors
"""Position mismatch at match point [km]"""
position_mismatch(sol::SimsFlanaganSolution) =
    SVector{3}(sol.mismatch[1], sol.mismatch[2], sol.mismatch[3])

"""Velocity mismatch at match point [km/s]"""
velocity_mismatch(sol::SimsFlanaganSolution) =
    SVector{3}(sol.mismatch[4], sol.mismatch[5], sol.mismatch[6])

"""Mass mismatch at match point [kg]"""
mass_mismatch(sol::SimsFlanaganSolution) = sol.mismatch[7]

"""Position mismatch magnitude [km]"""
position_mismatch_norm(sol::SimsFlanaganSolution) = norm(position_mismatch(sol))

"""Velocity mismatch magnitude [km/s]"""
velocity_mismatch_norm(sol::SimsFlanaganSolution) = norm(velocity_mismatch(sol))

export position_mismatch, velocity_mismatch, mass_mismatch
export position_mismatch_norm, velocity_mismatch_norm

# =============================================================================
# Pretty Printing
# =============================================================================

function Base.show(io::IO, sc::Spacecraft)
    print(
        io,
        "Spacecraft(dry=$(sc.dry_mass) kg, wet=$(sc.wet_mass) kg, thrust=$(sc.thrust) N, Isp=$(sc.isp) s)",
    )
end

function Base.show(io::IO, sc::SEPSpacecraft)
    print(
        io,
        "SEPSpacecraft(dry=$(sc.dry_mass) kg, wet=$(sc.wet_mass) kg, thrust@1AU=$(sc.thrust_ref) N, Isp=$(sc.isp) s)",
    )
end

function Base.show(io::IO, sail::SolarSail)
    a_c = characteristic_acceleration(sail)
    print(
        io,
        "SolarSail(dry=$(sail.dry_mass) kg, wet=$(sail.wet_mass) kg, area=$(sail.area) m², a_c=$(round(a_c*1000, digits=4)) mm/s²)",
    )
end

function Base.show(io::IO, opts::SimsFlanaganOptions)
    print(io, "SimsFlanaganOptions(n_segments=$(opts.n_segments), n_fwd=$(opts.n_fwd))")
end

function Base.show(io::IO, prob::SimsFlanaganProblem)
    sc_type = typeof(prob.spacecraft).name.name
    print(
        io,
        "SimsFlanaganProblem($sc_type, tof=$(prob.tof) s, μ=$(prob.μ) km³/s², $(prob.options.n_segments) segments)",
    )
end

function Base.show(io::IO, sol::SimsFlanaganSolution)
    Δr = position_mismatch_norm(sol)
    Δv = velocity_mismatch_norm(sol)
    Δm = mass_mismatch(sol)
    print(
        io,
        "SimsFlanaganSolution(retcode=$(sol.retcode), Δv_total=$(round(sol.Δv_total, digits=4)) km/s, Δr=$(round(Δr, digits=4)) km, Δv=$(round(Δv, digits=6)) km/s, Δm=$(round(Δm, digits=4)) kg)",
    )
end

# =============================================================================
# SciMLBase Interface: remake
# =============================================================================

"""
    remake(prob::SimsFlanaganProblem; kwargs...)

Create a new problem with modified parameters.

# Keyword Arguments
Any field of SimsFlanaganProblem can be overridden:
- `r0`, `v0`, `rf`, `vf`: Boundary conditions
- `tof`: Time of flight
- `μ`: Gravitational parameter
- `spacecraft`: Spacecraft parameters
- `options`: Solver options

# Example
```julia
prob2 = remake(prob; tof=prob.tof * 1.1)  # 10% longer transfer
prob3 = remake(prob; spacecraft=new_spacecraft)
```
"""
function SciMLBase.remake(
    prob::SimsFlanaganProblem;
    r0 = nothing,
    v0 = nothing,
    rf = nothing,
    vf = nothing,
    tof = nothing,
    μ = nothing,
    spacecraft = nothing,
    options = nothing,
)
    return SimsFlanaganProblem(
        r0 === nothing ? prob.r0 : r0,
        v0 === nothing ? prob.v0 : v0,
        rf === nothing ? prob.rf : rf,
        vf === nothing ? prob.vf : vf,
        tof === nothing ? prob.tof : tof,
        μ === nothing ? prob.μ : μ,
        spacecraft === nothing ? prob.spacecraft : spacecraft,
        options === nothing ? prob.options : options,
    )
end
