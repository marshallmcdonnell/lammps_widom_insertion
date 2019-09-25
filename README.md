## Widom particle insertion using [LAMMPS](https://lammps.sandia.gov/) as a [shared library](https://lammps.sandia.gov/doc/Python_shlib.html)

This is for single monatomic particles only (no poly-atomic systems).

This is only a proof-of-concept project.

#### DISCLAIMER
**This project is no longer maintained!**

##### Why?
The work was moved to adding a `compute widom ...` command directly into LAMMPS, still a WIP at this branch (*Last update: 7/2018*):

https://github.com/marshallmcdonnell/lammps/tree/widom

The support for polyatomic systems was to be added in the lammps feature branch, not this project.

##### The reason moved to other project and not here

Due to larger effort in this project and isolation of this work to larger LAMMPS community.

## Reference to Theory

 * Daan Frenkel's and Berend Smit's seminal textbook, specifically Chapter 7 Free Energy Calcutions:
 
   [Understanding Molecular Simulation: From Algorithms to Applications, second edition](http://www.acmm.nl/molsim/frenkel_smit/)

* I also have a lecture covering some of the topics in Chapter 7 of the textbook with applications in LAMMPS:

   http://utkstair.org/clausius/docs/mse614/pdf/chemicalpotential_intro_v01.pdf

