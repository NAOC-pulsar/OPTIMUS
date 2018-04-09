#simulate
echo "Start"
##simulate system noise
echo "Simulate system noise. Create file: snoise.dat"
./../simulateSystemNoise -p FAST.params -o snoise.dat
#simulate a simple pulsar
echo "Simulate a simple pulsar. Create file: pulsar.dat" 
#./simulateSimplePsr_zww -p parkes.params -p pulsar1.params -o test.dat
./../simulateSimplePsr -p FAST.params -p pulsar.params -o pulsar.dat
##create search file
echo "Create search file. Create file: search.sf"
./../createSearchFile -o test.sf -p FAST.params -p digitiser.params -f snoise.dat -f pulsar.dat
