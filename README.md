# SimsFlanagan.jl

[![Build Status](https://github.com/HAMMERHEAD-Space/SimsFlanagan.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HAMMERHEAD-Space/SimsFlanagan.jl/actions/workflows/CI.yml?query=branch%3Amaster)

A Julia package for low-thrust trajectory optimization using the Sims-Flanagan transcription method. Supports multiple propulsion types including constant thrust, solar electric propulsion (SEP), and solar sails.

## Overview

The Sims-Flanagan method discretizes a low-thrust trajectory into segments, applying impulsive ΔV at segment midpoints to approximate continuous thrust. The trajectory is propagated forward from the initial state and backward from the final state, meeting at a match point where continuity is enforced.

## Installation

```julia
using Pkg
Pkg.add("SimsFlanagan")
```

## Quick Start

```julia
using SimsFlanagan
using LinearAlgebra  # for norm()

# Define spacecraft parameters
# mass = 1000 kg, thrust = 0.5 N, Isp = 3000 s
spacecraft = Spacecraft(1000.0, 0.5, 3000.0)

# Define boundary conditions
μ = 398600.4418  # Earth gravitational parameter [km³/s²]
r0 = [7000.0, 0.0, 0.0]     # Initial position [km]
v0 = [0.0, 7.546, 0.0]      # Initial velocity [km/s]
rf = [42164.0, 0.0, 1000.0] # Final position [km]
vf = [0.0, 3.075, 0.0]      # Final velocity [km/s]
tof = 86400.0 * 30          # Time of flight [s] (30 days)

# Create problem
prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, spacecraft;
    n_segments = 20,
    verbosity = 1,
)

# Solve
sol = simsflanagan_solve(prob)

# Access solution
println("Converged: ", sol.converged)
println("Total ΔV: ", sol.Δv_total, " km/s")
println("Final mismatch: ", norm(sol.mismatch))
```

## Propulsion Types

SimsFlanagan.jl supports three propulsion models:

### Constant Thrust (Chemical/Electric)

Standard propulsion with constant thrust magnitude, independent of position.

```julia
# mass = 1000 kg, thrust = 0.5 N, Isp = 3000 s
sc = Spacecraft(1000.0, 0.5, 3000.0)
```

### Solar Electric Propulsion (SEP)

Thrust scales with solar distance as `(r_ref/r)²` due to inverse-square law for solar power.

```julia
# mass = 1000 kg, thrust at 1 AU = 0.5 N, Isp = 3000 s
r_1AU = 1.495978707e8  # km
sep = SEPSpacecraft(1000.0, 0.5, 3000.0, r_1AU)

# Thrust at 2 AU will be 1/4 of thrust at 1 AU
```

### Solar Sail

Uses solar radiation pressure for propulsion. No propellant consumption—mass remains constant.

```julia
# mass = 100 kg, sail area = 1000 m², reflectivity = 0.9
sail = SolarSail(100.0, 1000.0, 0.9)

# Characteristic acceleration at 1 AU:
a_c = characteristic_acceleration(sail)  # [m/s²]
```

For solar sails:
- Throttle magnitude represents sail efficiency (0-1), related to cone angle α as cos²(α)
- Throttle direction is the thrust direction
- Mass is constant throughout the trajectory
- Acceleration scales as `(r_ref/r)²` with solar distance

## API Reference

### Types

- `Spacecraft(mass, thrust, isp)` - Constant thrust propulsion
- `SEPSpacecraft(mass, thrust_ref, isp, r_ref)` - Solar electric propulsion
- `SolarSail(mass, area, reflectivity, r_ref)` - Solar radiation pressure
- `SimsFlanaganOptions` - Solver options
- `SimsFlanaganProblem` - Problem definition
- `SimsFlanaganSolution` - Optimization result

### Functions

- `simsflanagan_problem(r0, v0, rf, vf, tof, μ, spacecraft; kwargs...)` - Create a problem
- `simsflanagan_solve(problem; kwargs...)` - Solve the problem
- `characteristic_acceleration(sail)` - Get solar sail characteristic acceleration [m/s²]

### Problem Options

```julia
prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, spacecraft;
    n_segments = 10,      # Number of trajectory segments
    n_fwd = nothing,      # Forward segments (default: n_segments ÷ 2)
    tol = 1e-6,           # Convergence tolerance
    max_iter = 1000,      # Maximum iterations
    verbosity = 0,        # Output verbosity (0 = silent)
)
```

### Solver Options

```julia
using OptimizationOptimJL: Optim

sol = simsflanagan_solve(prob;
    alg = Optim.LBFGS(),           # Optimization algorithm (gradient-based with AD)
    use_lambert_guess = Val(true), # Use Lambert solution for initial guess
    penalty_weight = 1e6,          # Penalty weight for constraints
    maxiters = 1000,               # Maximum iterations
)

# For derivative-free optimization:
sol = simsflanagan_solve(prob; alg = Optim.NelderMead())
```

## Examples

### Heliocentric SEP Transfer

```julia
using SimsFlanagan

# Sun parameters
μ_sun = 1.32712440018e11  # km³/s²
r_1AU = 1.495978707e8     # km

# SEP spacecraft: 1000 kg, 0.5 N at 1 AU, 3000s Isp
sep = SEPSpacecraft(1000.0, 0.5, 3000.0, r_1AU)

# Earth to Mars-like orbit
r0 = [r_1AU, 0.0, 0.0]
v0 = [0.0, 29.78, 0.0]      # Earth orbital velocity
rf = [0.0, 1.524 * r_1AU, 0.0]  # Mars orbit radius
vf = [-24.1, 0.0, 0.0]
tof = 86400.0 * 180         # 180 days

prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ_sun, sep; n_segments = 20)
sol = simsflanagan_solve(prob)
```

### Solar Sail Trajectory

```julia
using SimsFlanagan

μ_sun = 1.32712440018e11
r_1AU = 1.495978707e8

# Lightweight solar sail: 50 kg, 2000 m² sail
sail = SolarSail(50.0, 2000.0, 0.9)

# Characteristic acceleration
a_c = characteristic_acceleration(sail)
println("Characteristic acceleration: $(a_c * 1000) mm/s²")

r0 = [r_1AU, 0.0, 0.0]
v0 = [0.0, 29.78, 0.0]
rf = [r_1AU * 0.7, 0.0, 0.0]  # Toward inner solar system
vf = [0.0, 35.0, 0.0]
tof = 86400.0 * 365          # 1 year

prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ_sun, sail; n_segments = 30)
# Use zero initial guess for solar sails (Lambert doesn't apply)
sol = simsflanagan_solve(prob; use_lambert_guess = Val(false))
```

## References

1. **Sims, J. A., & Flanagan, S. N.** (1999). Preliminary design of low-thrust interplanetary missions. *AAS/AIAA Astrodynamics Specialist Conference*, Paper AAS 99-338, Girdwood, Alaska.

2. **Yam, C. H., Di Lorenzo, D., & Izzo, D.** (2011). Low-thrust trajectory design as a constrained global optimization problem. *Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering*, 225(11), 1243-1251.

3. **Izzo, D.** (2012). PyGMO and PyKEP: Open source tools for massively parallel optimization in astrodynamics. *5th International Conference on Astrodynamics Tools and Techniques (ICATT)*, ESA/ESTEC.

4. **McInnes, C. R.** (1999). *Solar Sailing: Technology, Dynamics and Mission Applications*. Springer-Praxis.
