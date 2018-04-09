#simulate
echo "Start"
#simulate a simple pulsar
echo "Simulate a simple pulsar. Create file: pulsar.dat" 
#./simulateSimplePsr_zww -p parkes.params -p pulsar1.params -o test.dat
./../simulateSimplePsr_mcc -p FAST.params -p pulsar.params -o pulsarbinary.dat
python ../python/fitsio_combinePolBinary.py FP20180126_1-2GHz_Dec+4642_drifting_0572.fits pulsarbinary.dat new.fits
#python ../python/fitsio_combinePolBinary.py FP20171027_0-1GHz_Dec-0701_drifting_0455.fits pulsarbinary.dat new.fits
python ../python/fitsio_cutfreq.py new.fits 300 556 new_cut.fits
