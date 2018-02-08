#!bin/sh
gcc -O3 -lm -o inspectBinaryFile inspectBinaryFile.c simulate.c T2toolkit.c
gcc -O3 -lm -o simulateSystemNoise simulateSystemNoise.c simulate.c T2toolkit.c
gcc -O3 -lm -o simulateCal simulateCal.c simulate.c T2toolkit.c
gcc -lm -o createSearchFile createSearchFile.c simulate.c T2toolkit.c -L/nfshome/mcc/psrsoft/cfitsio/lib -lcfitsio -O3
gcc -O3 -lm -o simulateSimplePsr simulateSimplePsr.c simulate.c T2toolkit.c
gcc -O3 -lm -o simulateRFI simulateRFI.c simulate.c T2toolkit.c
gcc -O3 -lm -o simulateComplexPsr simulateComplexPsr.c simulate.c T2toolkit.c t1polyco.c tempo2pred.c cheby2d.c


gcc -O3 -lm -o simulateSimplePsr_zww simulateSimplePsr_zww.c simulate.c T2toolkit.c -w
