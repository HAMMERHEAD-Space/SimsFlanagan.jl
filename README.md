# SimsFlanagan.jl

[![Build Status](https://github.com/HAMMERHEAD-Space/SimsFlanagan.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/HAMMERHEAD-Space/SimsFlanagan.jl/actions/workflows/CI.yml?query=branch%3Amaster)

A Julia package for low-thrust trajectory optimization using the Sims-Flanagan transcription method.

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
    verbosity = 1
)

# Solve
sol = simsflanagan_solve(prob)

# Access solution
println("Converged: ", sol.converged)
println("Total ΔV: ", sol.Δv_total, " km/s")
println("Final mismatch: ", norm(sol.mismatch))
```

## API Reference

### Types

- `Spacecraft` - Spacecraft parameters (mass, thrust, Isp)
- `SimsFlanaganOptions` - Solver options (n_segments, n_fwd, tol, max_iter, verbosity)
- `SimsFlanaganProblem` - Problem definition (boundary conditions, spacecraft, options)
- `SimsFlanaganSolution` - Solution (throttles, masses, mismatch, Δv_total, converged)

### Functions

- `simsflanagan_problem(r0, v0, rf, vf, tof, μ, spacecraft; kwargs...)` - Create a problem
- `simsflanagan_solve(problem; kwargs...)` - Solve the problem

### Solver Options

```julia
sol = simsflanagan_solve(prob;
    alg = Optim.LBFGS(),           # Optimization algorithm (gradient-based with AD)
    use_lambert_guess = true,       # Use Lambert solution for initial guess
    penalty_weight = 1e6,           # Penalty weight for constraints
    maxiters = 1000                 # Maximum iterations
)
```

## References

1. **Sims, J. A., & Flanagan, S. N.** (1999). Preliminary design of low-thrust interplanetary missions. *AAS/AIAA Astrodynamics Specialist Conference*, Paper AAS 99-338, Girdwood, Alaska.

2. **Yam, C. H., Di Lorenzo, D., & Izzo, D.** (2011). Low-thrust trajectory design as a constrained global optimization problem. *Proceedings of the Institution of Mechanical Engineers, Part G: Journal of Aerospace Engineering*, 225(11), 1243-1251. DOI: 10.1177/0954410011401686

3. **Izzo, D.** (2012). PyGMO and PyKEP: Open source tools for massively parallel optimization in astrodynamics. *5th International Conference on Astrodynamics Tools and Techniques (ICATT)*, ESA/ESTEC.