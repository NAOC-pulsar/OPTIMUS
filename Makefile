cc = gcc
objects = simulateSimplePsr_mccV1 inspectBinaryFile simulateSystemNoise createSearchFile simulateRFI simulateComplexPsr


#simulate simple pulsar changed by mcc
simulateSimplePsr_mccV1:simulateSimplePsr_mccV1.c
	$(cc) -o simulateSimplePsr_mccV1 simulateSimplePsr_mccV1.c simulate.c T2toolkit.c -w -lm -O3 -lgomp -fopenmp -DHAVE_OPENMP 
#read file information
inspectBinaryFile:inspectBinaryFile.c
	$(cc) -o inspectBinaryFile inspectBinaryFile.c simulate.c T2toolkit.c -lm -O3
#simulate system noise	
simulateSystemNoise:simulateSystemNoise.c
	$(cc) -o simulateSystemNoise simulateSystemNoise.c simulate.c T2toolkit.c -lm -O3
#
simulateCal:simulateCal.c
	$(cc) -o simulateCal simulateCal.c simulate.c T2toolkit.c -lm -O3
#create search file 
createSearchFile:createSearchFile.c
	$(cc) -o createSearchFile createSearchFile.c simulate.c T2toolkit.c -L/usr/lib/ -I/usr/include/ -lcfitsio -O3 -lm -w
#simulate RFI
simulateRFI:simulateRFI.c
	$(cc) -o simulateRFI simulateRFI.c simulate.c T2toolkit.c -lm -O3
#simulate complex pulsar
simulateComplexPsr:simulateComplexPsr.c
	$(cc) -o simulateComplexPsr simulateComplexPsr.c simulate.c T2toolkit.c t1polyco.c tempo2pred.c cheby2d.c -lm -O3


#rm files compiled

.PHONY : clean
clean :
	rm $(objects)
