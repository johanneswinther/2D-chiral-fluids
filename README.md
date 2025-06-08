# 2D-chiral-fluids

This repository includes LAMMPS computes for finding the center-of-mass coarse-grained kinetic and virial stress for molecules, including the spin "temperature" and translational "temperature" on a molecular level.

Note 1: The virial compute only works, if molecule ID's are shared between processors, which isn't done by default. This can be done by using the "full" atom style, it can also be custom implemented in LAMMPS (there is guidance inside the source code on how to do this) or possibly implemented using the fix with for example: "fix shareMol all property/atom mol ghost yes".

Note 2: The minimum image convention used in the virial compute could give problems if cutoff > boxsize/2. It is a good idea to use large box sizes to avoid this.

Note 3: These computes are currently still being worked on, so changes or corrections might happen.
