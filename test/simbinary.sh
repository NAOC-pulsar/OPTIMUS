#simulate
echo "Start"
#simulate a simple pulsar
echo "Simulate a simple pulsar. Create file: pulsar.dat" 
#./simulateSimplePsr_zww -p parkes.params -p pulsar1.params -o test.dat
#time ../simulateSimplePsr_mcc -p FAST.params -p pulsar.params -o pulsarbinary.dat

#python ../python/fitsio_combinePolBinary.py FP20180126_1-2GHz_Dec+4642_drifting_0572.fits pulsarbinary.dat new.fits
time python ../python/fitsio_combinePolBinary.py FP20180103_0-1GHz_Dec+39.7_drifting_0305.fits pulsarbinary.dat new.fits
time python ../python/fitsio_cutfreq.py new.fits 300 556 new_cut.fits
