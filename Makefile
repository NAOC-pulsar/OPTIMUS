cc = gcc

objects = simulateSimplePsr_mcc inspectBinaryFile simulateSystemNoise createSearchFile simulateRFI simulateComplexPsr

all:$(objects)

#simulate simple pulsar changed by mcc
simulateSimplePsr_mcc:simulateSimplePsr_mcc.c
	$(cc) -o simulateSimplePsr_mcc simulateSimplePsr_mcc.c simulate.c T2toolkit.c -w -lm -O3 -lgomp -fopenmp -DHAVE_OPENMP 
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
	$(cc) -o createSearchFile createSearchFile.c simulate.c T2toolkit.c -L/public/home/mcc/psrsoft/cfitsio/lib -I/public/home/mcc/psrsoft/cfitsio/include -lcfitsio -O3 -lm -w
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
