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
This implementation is a port of the Sims-Flanagan algorithm from the fantastic ESA toolbox [Pykep](https://esa.github.io/pykep/).
