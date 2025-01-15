#!/bin/bash

#./bash_serial.sh 1 4929 csv/gemat11.ncol > output/gemat11/cilk1_g11.txt
#./bash_serial.sh 1232 4929 csv/gemat11.ncol > output/gemat11/cilk1232_g11.txt
#./bash_serial.sh 2464 4929 csv/gemat11.ncol > output/gemat11/cilk2464_g11.txt
#./bash_serial.sh 4929 4929 csv/gemat11.ncol > output/gemat11/cilk4929_g11.txt
#./bash_serial.sh 6161 4929 csv/gemat11.ncol > output/gemat11/cilk6161_g11.txt

#./bash_serial.sh 1 5300 csv/bcspwr10.ncol > output/bcspwr10/cilk1_b10.txt
#./bash_serial.sh 1325 5300 csv/bcspwr10.ncol > output/bcspwr10/cilk1325_b10.txt
#./bash_serial.sh 2650 5300 csv/bcspwr10.ncol > output/bcspwr10/cilk2650_b10.txt
#./bash_serial.sh 5300 5300 csv/bcspwr10.ncol > output/bcspwr10/cilk5300_b10.txt
#./bash_serial.sh 6625 5300 csv/bcspwr10.ncol > output/bcspwr10/cilk6625_b10.txt


./bash_serial.sh 1 10974 csv/bcsstk17.ncol > output/bcsstk17/cilk1_b17.txt
./bash_serial.sh 2743 10974 csv/bcsstk17.ncol > output/bcsstk17/cilk2743_b17.txt
./bash_serial.sh 5487 10974 csv/bcsstk17.ncol > output/bcsstk17/cilk5487_b17.txt
./bash_serial2.sh 10974 csv/bcsstk17.ncol > output/bcsstk17/cilkX_b17.txt
#./bash_serial.sh 10974 10974 csv/bcsstk17.ncol > output/bcsstk17/cilk10974_b17.txt
#./bash_serial.sh 13717 10974 csv/bcsstk17.ncol > output/bcsstk17/cilk13717_b17.txt

./bash_serial.sh 1 36476 csv/shock-9.ncol > output/shock9/cilk1_s9.txt
./bash_serial.sh 6000 36476 csv/shock-9.ncol > output/shock9/cilk6000_s9.txt
./bash_serial2.sh 36476 csv/shock-9.ncol > output/shock9/cilkX_s9.txt
#./bash_serial.sh 9119 36476 csv/shock-9.ncol > output/shock9/cilk9119_s9.txt
#./bash_serial.sh 18238 36476 csv/shock-9.ncol > output/shock9/cilk18238_s9.txt
#./bash_serial.sh 36476 36476 csv/shock-9.ncol > output/shock9/cilk36476_s9.txt
#./bash_serial.sh 45595 36476 csv/shock-9.ncol > output/shock9/cilk45595_s9.txt

python3 ./read2csv.py 1961 "output/netz4504/cilk1_n4504.txt"
python3 ./read2csv.py 1961 "output/netz4504/cilk490_n4504.txt"
python3 ./read2csv.py 1961 "output/netz4504/cilk980_n4504.txt"
python3 ./read2csv.py 1961 "output/netz4504/cilk1961_n4504.txt"
python3 ./read2csv.py 1961 "output/netz4504/cilk2451_n4504.txt"

python3 ./read2csv.py 4929 "output/gemat11/cilk1_g11.txt"
python3 ./read2csv.py 4929 "output/gemat11/cilk1232_g11.txt"
python3 ./read2csv.py 4929 "output/gemat11/cilk2464_g11.txt"
python3 ./read2csv.py 4929 "output/gemat11/cilk4929_g11.txt"
python3 ./read2csv.py 4929 "output/gemat11/cilk6161_g11.txt"

python3 ./read2csv.py 5300 "output/bcspwr10/cilk1_b10.txt"
python3 ./read2csv.py 5300 "output/bcspwr10/cilk1325_b10.txt"
python3 ./read2csv.py 5300 "output/bcspwr10/cilk2650_b10.txt"
python3 ./read2csv.py 5300 "output/bcspwr10/cilk5300_b10.txt"
python3 ./read2csv.py 5300 "output/bcspwr10/cilk6625_b10.txt"

python3 ./read2csv.py 10974 "output/bcsstk17/cilk1_b17.txt"
python3 ./read2csv.py 10974 "output/bcsstk17/cilk2743_b17.txt"
python3 ./read2csv.py 10974 "output/bcsstk17/cilk5487_b17.txt"
python3 ./read2csv.py 10974 "output/bcsstk17/cilkX_b17.txt"
#python3 ./read2csv.py 10974 "output/bcsstk17/cilk10974_b17.txt"
#python3 ./read2csv.py 10974 "output/bcsstk17/cilk13717_b17.txt"

python3 ./read2csv.py 36476 "output/shock9/cilk1_s9.txt"
python3 ./read2csv.py 36476 "output/shock9/cilk6000_s9.txt"
python3 ./read2csv.py 36476 "output/shock9/cilkX_s9.txt"
#python3 ./read2csv.py 36476 "output/shock9/cilk9119_s9.txt"
#python3 ./read2csv.py 36476 "output/shock9/cilk18238_s9.txt"
#python3 ./read2csv.py 36476 "output/shock9/cilk36476_s9.txt"
#python3 ./read2csv.py 36476 "output/shock9/cilk45595_s9.txt"