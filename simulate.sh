#simulate
echo "Start"
#simulate system noise
echo "Simulate system noise. Create file: snoise.dat"
./simulateSystemNoise -p parkes.params -o snoise.dat
#simulate a simple pulsar
echo "Simulate a simple pulsar. Create file: pulsar.dat" 
#./simulateSimplePsr_zww -p parkes.params -p pulsar1.params -o test.dat
./simulateSimplePsr_mcc -p parkes.params -p pulsar1.params -o test.dat

#create search file
echo "Create search file. Create file: search.sf"
./createSearchFile -o test.sf -p parkes.params -p digitiser.params -f snoise.dat -f test.dat
#createSearchFile -o search.sf -p parkes.params -f snoise1.dat -f pulsar1.dat
#createSearchFile -o search.sf -p digitiser.params -f snoise1.dat -f pulsar1.dat
echo "View the '.sf' file"
#view the .sf file 
pfits_plot -f test.sf -s1 1 -s2 5
#draw the .sf file 
pfits_draw -f test.sf -s1 1 -s2 10
