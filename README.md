# SimsFlanagan.jl

[![Build Status](https://github.com/HAMMERHEAD-Space/SimsFlanagan.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HAMMERHEAD-Space/SimsFlanagan.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/HAMMERHEAD-Space/SimsFlanagan.jl/graph/badge.svg?token=RVSG7F2BNO)](https://codecov.io/gh/HAMMERHEAD-Space/SimsFlanagan.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia package for low-thrust trajectory optimization using the Sims-Flanagan transcription method. Supports multiple propulsion types including constant thrust, solar electric propulsion (SEP), and solar sails. Built on the SciML ecosystem with automatic differentiation support.

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

## Example

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

## References

1. **Sims, J. A., & Flanagan, S. N.** (1999). Preliminary design of low-thrust interplanetary missions. *AAS/AIAA Astrodynamics Specialist Conference*, Paper AAS 99-338, Girdwood, Alaska.

2. **Yam, C. H., Di Lorenzo, D., & Izzo, D.** (2011). Low-thrust trajectory design as a constrained global optimization problem. *Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering*, 225(11), 1243-1251.

3. **Izzo, D.** (2012). PyGMO and PyKEP: Open source tools for massively parallel optimization in astrodynamics. *5th International Conference on Astrodynamics Tools and Techniques (ICATT)*, ESA/ESTEC.

4. **McInnes, C. R.** (1999). *Solar Sailing: Technology, Dynamics and Mission Applications*. Springer-Praxis.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgements
This implementation is based on the fantastic ESA toolbox [Pykep](https://esa.github.io/pykep/)
