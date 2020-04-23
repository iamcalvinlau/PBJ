# PBJ
PBJ stands for Particles By Julia. Hopefully, this grows into a Julia-based particle simulation suite of functions and subroutines to be used for testing various particle-in-cell methods and techniques. Some benchmarks are shown.

### 1. Sheath Formation with Fixed Potential Boundaries
At the moment, it has been used for 1d3v benchmark simulations. The flattening from 3d to 1d is done by ignoring the other dimensions and only calling a 1d gather/scatter function and a 1d Poisson solver. 

This shows the formation of the sheath using protons and electrons injected from the left and being absorbed at the left and right. There is a thermal bath of particles on the left. The potential is held fixed at both ends.
![Sheath formation using protons and electrons](https://github.com/iamcalvinlau/PBJ/blob/master/figures/sheath_proton-electron_100ippc_1000eppc.gif)

Electrons and ions are injected and absorbed. To accomodate for this, the particle array size is increased as needed in chunks. At steady-state, this should become roughly constant in number of particles.
![Tracking electron injection and absorption](https://github.com/iamcalvinlau/PBJ/blob/master/figures/electron_tracking.png)

![Tracking ion injection and absorption](https://github.com/iamcalvinlau/PBJ/blob/master/figures/ion_tracking.png)


### Two Stream Instability with Periodic Boundaries


