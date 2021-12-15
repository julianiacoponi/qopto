# qopto
## Quantum Optomechanics scripts and bits

For the time being, this is just a repository for some Fortran scripts (`.f90` files) I have been sent and have edited to help myself understand. (The original versions of the scripts live in `/originals`.)

These scripts are compiled with [`gfortran`](https://gcc.gnu.org/wiki/GFortranBinaries) and produce, by default, an `a.out` binary executable, which will run the compiled programme.

For example, running the following:
```
gfortran -O2 1DSIMPLE.f90
./a.out
```
May ask for some choices to be made:
```
 Choose input parameters:
 1 = original parameters
 2 = artificial parameters
>>> 1
 Choose expressions to calculate noise vectors with:
 1 = original expressions (in original scripts)
 2 = slightly reworked expressions (using back-action susceptibilities)
 3 = my own manually-derived expressions
>>> 1
```
And then produce an output like:
```
 Delta, kappa, omega_M, Gamma_M, g_coupling, n = nbar * omega_M = (kB * TBATH)/hbar
  -1.25664E+06
   2.75224E+06
   8.04641E+05
   6.48136E-03
   3.33396E+05
   4.00000E+13
 nbar = (kB*TBATH)/(hbar*omega_M) = thermal bath phonon occupancy:
   4.97116E+07
 Back-action limited phonons:
   5.18723E-01
 Area under S_XX:
   9.90407E+00
 Area under S_heterodyne:
   9.90407E+00
 Area under S_homodyne:
   9.65312E+00
 Phonons from optomechancial cooling formula:
   3.35585E+00
 Phonons from S_XX sideband asymmetry:
   4.45204E+00
 Mechanical oscillator`s temperature /K:
   2.68672E-05
```

`TODO: add more explanation of what the output of the scripts is?`