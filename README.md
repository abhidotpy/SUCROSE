# SUCROSE

SUCROSE is a collection of high-performance modules for molecular dynamics simulation of biomolecular systems written primarily in FORTRAN (and partly in C++) supporting both simple point particles and rigid, anisotropic bodies (ellipsoids, not just spheres). It handles translational and rotational dynamics via quaternions, bond constraints via RATTLE, and thermostatting via Nosé-Hoover chains and Langevin dynamics, with a Monte Carlo barostat available for NPT sampling. It is suitable for simulation of small systems of spheres/rigid bodies.

If you know what a Gay-Berne potential is, this is for you.

## Features

- **Integrators**: Velocity Verlet (with and without holonomic constraints), Nosé-Hoover thermostatted variants, and Langevin dynamics for both translational and rotational degrees of freedom
- **Rotational dynamics**: Full quaternion-based rigid body propagation — no gimbal lock, no Euler angle bookkeeping
- **Pair potentials**:
  - Lennard-Jones (standard, WCA-truncated, and "welded" cutoff variants)
  - Morse potential
  - Gay-Berne anisotropic potential (standard, repulsive-only, welded, and chiral variants)
  - Ellipsoid Contact Potential (ECP) — a Perram-Wertheim-style contact-distance formulation, also with repulsive, welded, and chiral flavors
  - Debye-Hückel and bare Coulombic electrostatics
  - Point-dipole interactions with reaction-field electrostatics
- **Bonded terms**: harmonic and cosine-based bond-angle bending, proper/improper dihedrals
- **Constraints**: RATTLE-style bond length constraints for rigid topologies
- **Neighbor lists**: Verlet lists with skin-distance rebuild triggers
- **Reporters**: state data (thermodynamic time series), LAMMPS-style trajectory dumps, raw configuration snapshots, and full checkpoint save/restore
- **RNG**: a counter-based Philox4x32 generator for reproducible stochastic dynamics

## Module overview

| Module | File | Purpose |
|---|---|---|
| `system` | `system.f90` | Core state: positions, velocities, forces, quaternions, box, thermodynamic accumulators |
| `neighbour_list` | `system.f90` | Verlet neighbor list construction and rebuild logic |
| `constraints` | `system.f90` | RATTLE constraint solver (position and velocity passes) |
| `xmath` | `xmath.f90` | Vector/matrix algebra, quaternion products, Cholesky solver |
| `philox_rng` | `xmath.f90` | Counter-based PRNG + Box-Muller normal variates |
| `verlet`, `quaternion`, `nose_hoover`, `langevin` | `integrators.f90` | Integration schemes for translation, rotation, and thermostatting |
| `monte_carlo_barostat` | `integrators.f90` | MC volume moves for NPT |
| `lennard_jones`, `wca`, `morse`, `gay_berne`, `ecp`, `bonded`, `debye_huckel`, `coulumbic`, `dipole` | `force.f90` | Force field implementations |
| `state_data_reporter`, `trajectory_reporter`, `config_reporter`, `checkpoint_reporter` | `reporters.f90` | Output and restart file handling |

## Requirements

- A Fortran compiler with F2008+ support (gfortran ≥ 9, ifx, or equivalent)
- All codes are tested with Intel ifx LLVM compiler

## Building

The makefile_f90 python files compile the modules in dependency order. Simply run the respective makefile depending on your platform in the working directory of your code. Ensure the files exist in the folder named *fortran*.

```bash
mkdir build
mkdir modules
python makefile_f90_linux.py
```

A Makefile or CMake setup is too complicated for me.

## Basic usage pattern

```fortran
program example
    use system
    use verlet
    use reporters
    implicit none
    integer :: tsteps = 20000
    integer :: IX, JX, time
    real :: dt = 0.0005
    integer :: I, J
    real :: pot, fxij(3), txi(3), txj(3)

    call open_state_data_file( 1, "outputs/state_data.txt" )
    call open_trajectory_file( 2, "outputs/trajectory.txt" )

    call generate_fcc_lattice( 3, 0.1 )
    call set_velocity_to_temperature( 0.722, seed=15 )
    call set_periodic_boundary

    do IX = 1, N
        call set_atom_type( IX, 1 )
    enddo

    call initialize_system
    call report_trajectory( 2, 0 )

    do time = 1, tsteps
        call vv_initial_step( dt )

        call initialize_forces
        do I = 1, N-1
            do J = I+1, N

                call lj_calculate_forces( I, J, 1.0, 1.0, 5.0, pot, fxij )

                pot_eng = pot_eng + pot

                FX(I) = FX(I) + fxij(1)
                FY(I) = FY(I) + fxij(2)
                FZ(I) = FZ(I) + fxij(3)

                FX(J) = FX(J) - fxij(1)
                FY(J) = FY(J) - fxij(2)
                FZ(J) = FZ(J) - fxij(3)

            enddo
        enddo

        call vv_final_step( dt )

        call calculate_thermo
        call report_state_data( 1, time )

        if( mod(time, 100) == 0 ) then
            write(*, "(I0, 2x 2(F12.5))") time, pot_eng, temp
            call report_trajectory( 2, time )
        end if

    enddo

    call close_state_data_file( 1 )
end program example
```

Force evaluation isn't wired up automatically — you loop over neighbor list pairs yourself and call the relevant `*_calculate_forces` subroutine, accumulating into `FX`/`FY`/`FZ`/`TX`/`TY`/`TZ` and `pot_eng`. This gives full control over which potentials apply to which atom types, at the cost of writing that loop yourself.

## License

GNU General Public License v3.0

---
