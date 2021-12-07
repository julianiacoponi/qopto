# qopto
## Quantum Optomechanics scripts and bits

For the time being, this is just a repository for some Fortran scripts (`.f90` files).

These are compiled with [`gfortran`](https://gcc.gnu.org/wiki/GFortranBinaries) and produce, by default, an `a.out` binary executable, which will run the compiled programme.

For example, running the following:
```
➜  qopto git:(master) ✗ gfortran -O2 1DSIMPLE.f90
➜  qopto git:(master) ✗ ./a.out
```
May produce an output like:
```
 Chosen method           3
 Delta, kappa, omega_x, Gamma_M, g_coupling, kB*T/hbar
  -1.25664E+06
   2.75224E+06
   8.04641E+05
   6.48136E-03
   3.33396E+05
   4.00000E+13
 nbar (thermal bath phonon occupancy) =
      4.97E+07
 x-axis back action limited phonons =
      5.19E-01
 Area under S_XX =
      9.90E+00
 x-axis phonons from optomechancial cooling formula:
      3.36E+00
 x-axis phonons from S_XX sideband asymmetry:
      4.45E+00
 Mechanical oscillator`s temperature /K
      2.69E-05

```