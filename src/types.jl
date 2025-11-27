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

"""
    Spacecraft{MT, TT, IT}

Constant-thrust spacecraft (chemical or electric propulsion at constant power).

# Fields
- `mass::MT`: Initial spacecraft mass [kg]
- `thrust::TT`: Maximum thrust magnitude [N]
- `isp::IT`: Specific impulse [s]

# Notes
The exhaust velocity is computed as `vex = isp * g0` where `g0 = 9.80665 m/s²`.
Thrust is independent of distance from the Sun.
"""
struct Spacecraft{MT<:Number,TT<:Number,IT<:Number} <: AbstractSpacecraft
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
- `mass::MT`: Initial spacecraft mass [kg]
- `thrust_ref::TT`: Maximum thrust at reference distance [N]
- `isp::IT`: Specific impulse [s]
- `r_ref::RT`: Reference distance [km] (default: 1 AU = 1.495978707e8 km)

# Notes
At distance r from the Sun:
  thrust(r) = thrust_ref × (r_ref / r)²

This models the inverse-square falloff of solar power with distance.
"""
struct SEPSpacecraft{MT<:Number,TT<:Number,IT<:Number,RT<:Number} <: AbstractSpacecraft
    mass::MT
    thrust_ref::TT
    isp::IT
    r_ref::RT

    function SEPSpacecraft(
        mass::MT,
        thrust_ref::TT,
        isp::IT,
        r_ref::RT = 1.495978707e8,  # 1 AU in km
    ) where {MT<:Number,TT<:Number,IT<:Number,RT<:Number}
        mass > 0 || throw(ArgumentError("mass must be positive"))
        thrust_ref > 0 || throw(ArgumentError("thrust_ref must be positive"))
        isp > 0 || throw(ArgumentError("isp must be positive"))
        r_ref > 0 || throw(ArgumentError("r_ref must be positive"))
        new{MT,TT,IT,RT}(mass, thrust_ref, isp, r_ref)
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
    mass::MT
    area::AT
    reflectivity::CT
    r_ref::RT

    function SolarSail(
        mass::MT,
        area::AT,
        reflectivity::CT = 0.9,
        r_ref::RT = 1.495978707e8,  # 1 AU in km
    ) where {MT<:Number,AT<:Number,CT<:Number,RT<:Number}
        mass > 0 || throw(ArgumentError("mass must be positive"))
        area > 0 || throw(ArgumentError("area must be positive"))
        0 <= reflectivity <= 1 || throw(ArgumentError("reflectivity must be in [0, 1]"))
        r_ref > 0 || throw(ArgumentError("r_ref must be positive"))
        new{MT,AT,CT,RT}(mass, area, reflectivity, r_ref)
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
           (SPEED_OF_LIGHT * sail.mass)
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
    return sail.mass * acceleration  # [N]
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
# Solver Options
# =============================================================================

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
}
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

# Fields
- `problem::SimsFlanaganProblem`: The original problem
- `throttles::AbstractVector{AbstractVector}`: Throttle vectors at each segment [dimensionless]
- `masses::AbstractVector`: Mass at each segment boundary [kg]
- `mismatch::AbstractVector`: Final mismatch at match point [Δr; Δv; Δm]
- `Δv_total::Number`: Total ΔV expended [km/s] (or effective ΔV for sails)
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

# =============================================================================
# Pretty Printing
# =============================================================================

function Base.show(io::IO, sc::Spacecraft)
    print(io, "Spacecraft(mass=$(sc.mass) kg, thrust=$(sc.thrust) N, Isp=$(sc.isp) s)")
end

function Base.show(io::IO, sc::SEPSpacecraft)
    print(
        io,
        "SEPSpacecraft(mass=$(sc.mass) kg, thrust@1AU=$(sc.thrust_ref) N, Isp=$(sc.isp) s)",
    )
end

function Base.show(io::IO, sail::SolarSail)
    a_c = characteristic_acceleration(sail)
    print(
        io,
        "SolarSail(mass=$(sail.mass) kg, area=$(sail.area) m², a_c=$(round(a_c*1000, digits=4)) mm/s²)",
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
    status = sol.converged ? "converged" : "not converged"
    mismatch_norm = norm(sol.mismatch)
    print(
        io,
        "SimsFlanaganSolution($(status), Δv=$(sol.Δv_total) km/s, mismatch=$(mismatch_norm))",
    )
end
