# Heating Discs and Levitation
**How to run the codes?** Download all the files from https://github.com/PaulMcMillan-Astro/GalPot.git. Then, run the file newpotcalc4.cpp with the following command: g++ -O3 newpotcalc4.cpp -std=c++11 -pthread ./GalPot-master/obj/libPot.a -Wl,-rpath,./GalPot-master/obj -I /usr/include/gsl/1.16-gcc/include/ -L /usr/include/gsl/1.16-gcc/lib/ -lgsl -lgslcblas -lm -I /usr/include/gsl/gsl_spline.h:47:1/ -L /usr/include/gsl/gsl_spline.h:47:1/ -lgsl -lgslcblas -lm -o ./a

**What is the algorithm about?** To investigate the stellar levitation upon a heating disk, the algorithm contains simple orbit simulations of stars in the McMillan (2017) potential, a leapfrog algorithm to track the trajectories of stars within the potential, ramp model, Gaussian velocity distribution, Fourier analysers and so on.  

**Levitation:** The following diagram demonstrates the levitaion process after modelling the velocity of 500 stellar orbits as having a Gaussian distribution.

<img src="https://github.com/JialunSimonLiu/Levitation/blob/main/Picture/Picture1.png" width="800" height="600"/>

**Result**: the colour indicated the initial velocity of the stellar orbit in z direction. Levitation occurs when Rmax is between 12.8(kpc) and 12.9(kpc).

**1:1 resonance:** The following plot containes the resoances in a stellar orbit with a fixed potential.

<img src="https://github.com/JialunSimonLiu/Levitation/blob/main/Picture/Picture2.png" width="700" height="500"/>

**Result**: the colour indicated the initial velocity of the stellar orbit in z direction. Levitation occurs when Rmax is between 12.8(kpc) and 12.9(kpc).
