# Heating Discs and Levitation
**How to run the codes?** Download all the files from https://github.com/PaulMcMillan-Astro/GalPot.git. Then, run the file newpotcalc4.cpp with the following command: g++ -O3 newpotcalc4.cpp -std=c++11 -pthread ./GalPot-master/obj/libPot.a -Wl,-rpath,./GalPot-master/obj -I /usr/include/gsl/1.16-gcc/include/ -L /usr/include/gsl/1.16-gcc/lib/ -lgsl -lgslcblas -lm -I /usr/include/gsl/gsl_spline.h:47:1/ -L /usr/include/gsl/gsl_spline.h:47:1/ -lgsl -lgslcblas -lm -o ./a

![Alt Text](https://github.com/JialunSimonLiu/Levitation/blob/main/Picture.png)
<img src="https://github.com/JialunSimonLiu/Levitation/blob/main/Picture.png" width="800" height="600"/>
