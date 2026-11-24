# Langevin bivariate integrator for LAMMPS

A native C++ extension for **LAMMPS** (originaly implemented with version `22Jul2025` but expected to work with most LAMMPS versions) implementing an exact bivariate stochastic integrator for underdamped and overdamped Langevin dynamics based on Chandrasekhar's joint phase-space covariance matrix [Santamaria2016].

## Installation

1. Copy `fix_lang_chand.cpp` and `fix_lang_chand.h` into the `src/` directory of your LAMMPS source tree:
   ```bash
   cp fix_lang_chand.cpp fix_lang_chand.h /path/to/lammps/src/
   ```
2. Recompile LAMMPS using your preferred build system (CMake or Make):
   ```bash
   cd /path/to/lammps/build
   cmake --build . -j
   ```

## Usage & Syntax

```lammps
fix ID group-ID langevin/chand Tstart Tstop damp seed
```

* `Tstart`, `Tstop`: Target temperature at the beginning and end of the run (temperature units).
* `damp`: Damping parameter $\tau_v = \gamma^{-1}$ (time units).
* `seed`: Positive integer seed for the random number generator.

> **Important Note:** Unlike the standard `fix langevin` (which requires a separate `fix nve`), `fix langevin/chand` performs the full time integration of both positions and velocities (`time_integrate = 1`). Do **not** combine it with `fix nve`.

### Example
```lammps
# Canonical (NVT) integration at T* = 2.50 with damp = 0.10
fix 1 all langevin/chand 2.50 2.50 0.10 482719
```

## Citation

If you use this code in your research, please cite:
* R. Santamaria, A. Alvarez de la Paz, L. Roskop, and L. Adamowicz; *Statistical Contact Model for Confined Molecules*, Stat. Phys. **164**, 1000-1025 (2016).
* O. Gutiérrez-Varela, M. A. Alderete-Lezama, and R. Santamaria, *"Langevin Dynamics via Bivariate Stochastic Integration"*, submitted (2026).
