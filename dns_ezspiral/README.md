# EZSpiral files

Initial conditions and driver files are located as compressed archive files for each model. The following changes in the code are needed to run the different models:

* Bar-Eiswirth:
  * ezstep.h: change G(u,v)
* FHN
  * ezstep.h: change U_KINETICS
  * eztip.c: change V_CONT
* Karma:
  * ezgraphGL.c: change colors for Karma model as indicated
  * ezstep.h: change U_KINETICS, G(u,v), U_THRESHOLD(v)
