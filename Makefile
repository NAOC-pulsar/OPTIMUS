simulateSimplePsr_mcc:simulateSimplePsr_mcc.c
	gcc -O3 -o simulateSimplePsr_mcc simulateSimplePsr_mcc.c simulate.c T2toolkit.c -w -lm -lgomp -fopenmp -DHAVE_OPENMP 
