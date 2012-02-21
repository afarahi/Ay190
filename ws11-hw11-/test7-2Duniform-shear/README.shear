This is a simple test of a 2D uniform shear flow with SPH, showing the
intrinsic "diffusion" that can be caused by particle rearrangement in
highly anisotropic flows.

Type "make clean" and "make" to re-compile the 2DSPMHD binary.

Run the code using ./2DSPMHD shear.in
Or with output to a file: ./2DSPMHD shear.in >& shear.out &

Plot the results using "nsplash shear_0*.dat"

Things to try
--------------
- Run the code at with the default parameters:

 ./2DSPMHD shear.in >& shear.out &

- Plot the particle positions and see what happens as a function of time

 nsplash -y 2 -x 1 -dev /xw shear_0*.dat

  You will see that after a large number of sound crossing times (about t=20)
 the particles transition to a glass-like arrangement, and after this time
 the kinetic energy in the shear flow starts to decay as a function of time.

- Plot vx as a function of y to see what happens to the velocity profile.
 
  The overall effect of the particle rearrangement that the shear flow diffuses
 in amplitude as a function of time due to the conversion of kinetic energy in 
 the shear flow into small scale transverse motions, and in turn, heat
 as these motions are dissipated by the artificial viscosity.

- You may wish to verify that the artificial viscosity parameters have very
 little influence on the velocity profile, and similarly that turning on the
 "XSPH" smoothing of the velocity field also makes no difference:
 1  0.500                            ! use xsph, parameter

 This problem illustrates the "intrinsic" minimum level of viscosity present  in
SPH simulations of shear flows because of the isotropy of the kernel interaction
(meaning that particles want to form a locally isotropic grid arrangement
instead of the highly anisotropic structure that would be traced by tracer
particles).

Added by Daniel Price, July 2010
