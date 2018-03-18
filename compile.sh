#!bin/sh
gcc -o inspectBinaryFile inspectBinaryFile.c simulate.c T2toolkit.c -O3 -lm 
gcc -o simulateSystemNoise simulateSystemNoise.c simulate.c T2toolkit.c -O3 -lm 
gcc -o simulateCal simulateCal.c simulate.c T2toolkit.c -O3 -lm 
gcc -o createSearchFile createSearchFile.c simulate.c T2toolkit.c -L/usr/lib/ -I/usr/include/ -w -lcfitsio -O3 -lm 
gcc -o simulateSimplePsr simulateSimplePsr.c simulate.c T2toolkit.c -O3 -lm 
gcc -o simulateRFI simulateRFI.c simulate.c T2toolkit.c -O3 -lm 
gcc -o simulateComplexPsr simulateComplexPsr.c simulate.c T2too lkit.c t1polyco.c tempo2pred.c cheby2d.c -O3 -lm 
gcc -o simulateSimplePsr_mccV1 simulateSimplePsr_mccV1.c simulate.c T2toolkit.c -w -lm -O3 -lgomp -fopenmp -DHAVE_OPENMP 

#gcc -O3 -lm -o simulateSimplePsr_zww simulateSimplePsr_zww.c simulate.c T2toolkit.c -w
