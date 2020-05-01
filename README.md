# PBJ
PBJ stands for Particles By Julia. Hopefully, this grows into a Julia-based particle simulation suite of functions and subroutines to be used for testing various particle-in-cell methods and techniques. 

I have a personal list of interesting algorithms to look at
1. Speed-limited PIC (SLPIC)
2. Particle-in-Fourier (PIF) 
3. Phase-space remapping
4. Asynchronous PIC

So far, SLPIC has been explored. Some PBJ benchmark cases are shown.

### 1. Sheath Formation with Fixed Potential Boundaries
At the moment, it has been used for 1d3v benchmark simulations. The flattening from 3d to 1d is done by ignoring the other dimensions and only calling a 1d gather/scatter function and a 1d Poisson solver. 

#### a. Without speed-limiting effects
This shows the formation of the sheath using protons and electrons injected from the left and being absorbed at the left and right. There is a thermal bath of particles on the left. The potential is held fixed at both ends.
![Sheath formation using protons and electrons](https://github.com/iamcalvinlau/PBJ/blob/master/figures/sheath_proton-electron_100ippc_1000eppc.gif)

Electrons and ions are injected and absorbed. To accomodate for this, the particle array size is increased as needed in chunks. At steady-state, this should become roughly constant in number of particles.
![Tracking electron injection and absorption](https://github.com/iamcalvinlau/PBJ/blob/master/figures/electron_tracking.png)

![Tracking ion injection and absorption](https://github.com/iamcalvinlau/PBJ/blob/master/figures/ion_tracking.png)

#### b. With speed-limiting effects
This shows the same formation of the sheath, almost the same as above but uses argon ions and electrons like the SLPIC paper. However, the time-step is 320x larger, and the simulation incorporates speed-limiting effects (see SLPIC paper).
![SL sheath formation using argon and electrons.](https://github.com/iamcalvinlau/PBJ/blob/master/figures/slpic_sheath_larger-ppc.gif)

### 2. Two Stream Instability with Periodic Boundaries

This shows the two stream instability of electrons streaming with drift-velocity into stationary ions (with 100x electron mass). The boundary is periodic.
![Two-stream instability with mass-ratio=100](https://github.com/iamcalvinlau/PBJ/blob/master/figures/two_stream_mass-ratio%3D100.gif)

The fastest growing mode is when the period of the drift motion is equal to the plasma frequency. In this case, it is dominated by the m=3 (wavelength = 1/3 simulation length) component. The theory is simple for this case and is compared to the simulation results.
![Comparison of theory and simulation](https://github.com/iamcalvinlau/PBJ/blob/master/figures/two-stream_instability.png)
