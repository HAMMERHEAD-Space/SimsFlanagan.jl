# SimsFlanagan.jl

[![Build Status](https://github.com/HAMMERHEAD-Space/SimsFlanagan.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HAMMERHEAD-Space/SimsFlanagan.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/HAMMERHEAD-Space/SimsFlanagan.jl/graph/badge.svg?token=RVSG7F2BNO)](https://codecov.io/gh/HAMMERHEAD-Space/SimsFlanagan.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package for low-thrust trajectory optimization using the Sims-Flanagan transcription method. Supports multiple propulsion types including constant thrust, solar electric propulsion (SEP), and solar sails. Built on the SciML ecosystem with automatic differentiation support.

## Overview

The Sims-Flanagan method discretizes a low-thrust trajectory into segments, applying impulsive ΔV at segment midpoints to approximate continuous thrust. The trajectory is propagated forward from the initial state and backward from the final state, meeting at a match point where continuity is enforced.

### Key Features

- **Multiple propulsion models**: Constant thrust, SEP, and solar sails
- **SciML interface**: Compatible with `solve()` and `remake()` patterns
- **Automatic differentiation**: Uses ForwardDiff for gradient-based optimization
- **Interior-point optimization**: MadNLP solver with MUMPS linear solver
- **Sundman transformation**: Adaptive segment sizing based on orbital distance
- **Lambert-based initialization**: Smart initial guesses for faster convergence

## Installation

```julia
using Pkg
Pkg.add("SimsFlanagan")
```

## Quick Start

```julia
using SimsFlanagan
using LinearAlgebra

# Define spacecraft: dry_mass=500kg, wet_mass=500kg, thrust=0.5N, Isp=3000s
spacecraft = Spacecraft(500.0, 500.0, 0.5, 3000.0)

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
sol = solve(prob)

# Access solution
println("Converged: ", sol.converged)
println("Total ΔV: ", sol.Δv_total, " km/s")
println("Final mismatch: ", norm(sol.mismatch))
```

## Propulsion Types

SimsFlanagan.jl supports three propulsion models through the `AbstractSpacecraft` type hierarchy:

### Constant Thrust (Chemical/Electric)

Standard propulsion with constant thrust magnitude, independent of position.

```julia
# dry_mass=500kg, wet_mass=500kg, thrust=0.5N, Isp=3000s
sc = Spacecraft(500.0, 500.0, 0.5, 3000.0)

# Access properties
total_mass = mass(sc)          # 1000.0 kg
vex = exhaust_velocity(sc)     # ~29420 m/s
```

### Solar Electric Propulsion (SEP)

Thrust scales with solar distance as `(r_ref/r)²` due to inverse-square law for solar power.

```julia
# dry_mass=500kg, wet_mass=500kg, thrust at 1AU=0.5N, Isp=3000s
r_1AU = 1.495978707e8  # km
sep = SEPSpacecraft(500.0, 500.0, 0.5, 3000.0, r_1AU)

# Thrust at 2 AU will be 1/4 of thrust at 1 AU
```

### Solar Sail

Uses solar radiation pressure for propulsion. No propellant consumption—mass remains constant.

```julia
# dry_mass=100kg, sail_area=1000m², reflectivity=0.9
sail = SolarSail(100.0, 1000.0, 0.9)

# Characteristic acceleration at 1 AU:
a_c = characteristic_acceleration(sail)  # [m/s²]
```

For solar sails:
- Throttle magnitude represents sail efficiency (0-1), related to cone angle α as cos²(α)
- Throttle direction is the thrust direction
- Mass is constant throughout the trajectory
- Acceleration scales as `(r_ref/r)²` with solar distance

## SciMLBase Interface

SimsFlanagan.jl follows SciML conventions for problem/solution handling:

### solve

```julia
sol = solve(prob)  # Solve with defaults
sol = solve(prob; max_iter=500, tol=1e-8)  # Custom options
```

### remake

Create a modified problem without rebuilding from scratch:

```julia
prob2 = remake(prob; tof=prob.tof * 1.1)  # 10% longer transfer
prob3 = remake(prob; spacecraft=new_spacecraft)
```

## Initial Guess Strategies

Control how the optimizer starts using `initial_guess_strategy`:

| Strategy | Description | Best For |
|----------|-------------|----------|
| `RandomGuess(seed=1234)` | Random throttle directions and magnitudes | Default, general use |
| `LambertGuess()` | Uses Lambert arc solution | Ballistic-like transfers |
| `ZeroGuess()` | All throttles = 0 | Simple problems |
| `ConstantGuess(direction, magnitude)` | Uniform thrust direction | Known transfer geometry |
| `RadialGuess(magnitude)` | Radial (Sun-pointing) direction | Solar sails |

```julia
# Use Lambert-based initialization
sol = solve(prob; initial_guess_strategy=LambertGuess())

# Use random with specific seed
sol = solve(prob; initial_guess_strategy=RandomGuess(seed=42))

# Custom constant direction
sol = solve(prob; initial_guess_strategy=ConstantGuess(direction=[1,0,0], magnitude=0.7))
```

You can also provide throttles directly:

```julia
my_throttles = [SVector{3}(0.5, 0.0, 0.0) for _ in 1:n_segments]
sol = solve(prob; initial_guess=my_throttles)
```

## Sundman Transformation

Enable adaptive segment durations based on distance from the central body:

```julia
prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, spacecraft;
    n_segments = 20,
    sundman_c = 1.0,  # Classic Sundman: dt ∝ r
)
```

Sundman exponent options:
- `0.0` (default): Equal time segments
- `1.0`: Classic Sundman (dt ∝ r), more resolution near central body
- `1.5`: Related to eccentric anomaly

This is particularly useful for SEP where thrust decreases with distance, requiring more time for maneuvers at outer regions.

## API Reference

### Types

| Type | Description |
|------|-------------|
| `Spacecraft(dry_mass, wet_mass, thrust, isp)` | Constant thrust propulsion |
| `SEPSpacecraft(dry_mass, wet_mass, thrust_ref, isp, r_ref)` | Solar electric propulsion |
| `SolarSail(dry_mass, area, reflectivity, r_ref)` | Solar radiation pressure |
| `SimsFlanaganOptions` | Solver configuration |
| `SimsFlanaganProblem` | Problem definition |
| `SimsFlanaganSolution` | Optimization result |

### Initial Guess Types

| Type | Description |
|------|-------------|
| `LambertGuess()` | Lambert arc-based initialization |
| `ZeroGuess()` | Zero throttles |
| `RandomGuess(seed)` | Random throttles |
| `ConstantGuess(direction, magnitude)` | Uniform direction |
| `RadialGuess(magnitude)` | Radial direction |

### Core Functions

| Function | Description |
|----------|-------------|
| `simsflanagan_problem(r0, v0, rf, vf, tof, μ, spacecraft; kwargs...)` | Create a problem |
| `solve(problem; kwargs...)` | Solve the problem (SciMLBase) |
| `simsflanagan_solve(problem; kwargs...)` | Legacy solver interface |
| `remake(problem; kwargs...)` | Create modified problem |
| `mass(spacecraft)` | Get total initial mass |
| `exhaust_velocity(spacecraft)` | Get exhaust velocity [m/s] |
| `characteristic_acceleration(sail)` | Get solar sail characteristic acceleration [m/s²] |

### Propagation Functions

| Function | Description |
|----------|-------------|
| `kepler_propagate(r0, v0, Δt, μ)` | Universal Kepler propagation |
| `propagate_segment(r, v, m, throttle, Δt, μ, sc)` | Propagate single segment |
| `propagate_leg(r0, v0, m0, throttles, Δt, μ, sc)` | Propagate multiple segments |
| `compute_mismatch(problem, throttles)` | Compute match point mismatch |
| `compute_total_Δv(throttles, Δt, spacecraft)` | Compute total ΔV |

### Solution Accessors

| Function | Description |
|----------|-------------|
| `position_mismatch(sol)` | Position mismatch at match point [km] |
| `velocity_mismatch(sol)` | Velocity mismatch at match point [km/s] |
| `mass_mismatch(sol)` | Mass mismatch at match point [kg] |
| `position_mismatch_norm(sol)` | ‖Δr‖ [km] |
| `velocity_mismatch_norm(sol)` | ‖Δv‖ [km/s] |

### Problem Options

```julia
prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, spacecraft;
    n_segments = 10,      # Number of trajectory segments
    n_fwd = nothing,      # Forward segments (default: n_segments ÷ 2)
    tol = 1e-6,           # Convergence tolerance
    max_iter = 1000,      # Maximum iterations
    verbosity = 0,        # Output verbosity (0 = silent)
    sundman_c = 0.0,      # Sundman transformation exponent
)
```

### Solver Options

```julia
sol = solve(prob;
    initial_guess_strategy = RandomGuess(),  # Initial guess strategy
    initial_guess = nothing,                 # Direct throttle vectors
    max_iter = 1000,                         # Maximum iterations
    tol = 1e-6,                              # Convergence tolerance
)
```

## Examples

### LEO Orbit Raising

```julia
using SimsFlanagan

μ = 398600.4418  # Earth gravitational parameter

# Circular orbit radii
r1, r2 = 7000.0, 7500.0
v1, v2 = sqrt(μ / r1), sqrt(μ / r2)

r0 = [r1, 0.0, 0.0]
v0 = [0.0, v1, 0.0]
rf = [r2, 0.0, 0.0]
vf = [0.0, v2, 0.0]
tof = 86400.0 * 2  # 2 days

sc = Spacecraft(500.0, 500.0, 0.5, 3000.0)
prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ, sc; n_segments=10)
sol = solve(prob)

println("ΔV: $(sol.Δv_total) km/s")
println("Converged: $(sol.converged)")
```

### Heliocentric SEP Transfer

```julia
using SimsFlanagan

# Sun parameters
μ_sun = 1.32712440018e11  # km³/s²
r_1AU = 1.495978707e8     # km

# SEP spacecraft: 500kg dry, 500kg propellant, 0.5N at 1AU, 3000s Isp
sep = SEPSpacecraft(500.0, 500.0, 0.5, 3000.0, r_1AU)

# Earth to Mars-like orbit
r0 = [r_1AU, 0.0, 0.0]
v0 = [0.0, 29.78, 0.0]          # Earth orbital velocity
rf = [0.0, 1.524 * r_1AU, 0.0]  # Mars orbit radius
vf = [-24.1, 0.0, 0.0]
tof = 86400.0 * 180             # 180 days

prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ_sun, sep; 
    n_segments = 20,
    sundman_c = 1.0,  # Use Sundman for better segment distribution
)
sol = solve(prob; initial_guess_strategy=LambertGuess())
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
tof = 86400.0 * 365            # 1 year

prob = simsflanagan_problem(r0, v0, rf, vf, tof, μ_sun, sail; n_segments = 30)
# Use radial guess for solar sails (Lambert doesn't apply)
sol = solve(prob; initial_guess_strategy=RadialGuess())
```

## Utility Functions

### AD-Safe Norm

For use in optimization with automatic differentiation:

```julia
using SimsFlanagan: safe_norm

x = [0.0, 0.0, 0.0]
n = safe_norm(x)  # Returns 0.0 with well-defined gradient
```

Standard `norm(x)` has undefined gradient at x=0, which can cause issues with ForwardDiff.

## Dependencies

- **Optimization.jl**: Unified optimization interface
- **MadNLP**: Interior-point nonlinear optimizer
- **ForwardDiff**: Automatic differentiation
- **Lambert.jl**: Lambert problem solver for initial guesses
- **AstroCoords.jl**: Astrodynamics utilities
- **StaticArrays.jl**: Fast fixed-size arrays
- **SciMLBase.jl**: SciML ecosystem interface

## References

1. **Sims, J. A., & Flanagan, S. N.** (1999). Preliminary design of low-thrust interplanetary missions. *AAS/AIAA Astrodynamics Specialist Conference*, Paper AAS 99-338, Girdwood, Alaska.

2. **Yam, C. H., Di Lorenzo, D., & Izzo, D.** (2011). Low-thrust trajectory design as a constrained global optimization problem. *Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering*, 225(11), 1243-1251.

3. **Izzo, D.** (2012). PyGMO and PyKEP: Open source tools for massively parallel optimization in astrodynamics. *5th International Conference on Astrodynamics Tools and Techniques (ICATT)*, ESA/ESTEC.

4. **McInnes, C. R.** (1999). *Solar Sailing: Technology, Dynamics and Mission Applications*. Springer-Praxis.

## License

MIT License - see [LICENSE](LICENSE) for details.
