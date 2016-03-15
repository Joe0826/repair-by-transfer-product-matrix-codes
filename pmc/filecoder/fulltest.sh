#!/bin/bash
set -e

TESTDATA=./testdata
META=coded_files/alice_metadata
TMPF=coded_files/alice_tmp
cp $TESTDATA/alice_golden ./alice
make


# HELPER FUNCTIONS

# range MIN MAX: echos an inclusive range.
#function range () {
    #for i in $(seq $1 $2); do
         #echo -n "$i "
    #done
#}

rm -f coded_files/*
rm -f $TMPF
./encoder_pmc 20 10 19 8 MSR alice #alpha=10
ls -v -1 coded_files/alice_node* | head -100 | xargs cat >> $TMPF
truncate -r $TESTDATA/alice_golden $TMPF # truncate systematic node data to filesize
cmp $TMPF $TESTDATA/alice_golden
echo "-----------------"
echo "MSR Systematic Test: Passed!"
echo "-----------------"

rm -f $TMPF
./encoder_pmc 20 10 19 8 MSR -rbt CYC alice #alpha=10
ls -v -1 coded_files/alice_node* | head -100 | xargs cat >> $TMPF
truncate -r $TESTDATA/alice_golden $TMPF # truncate systematic node data to filesize
cmp $TMPF $TESTDATA/alice_golden
echo "-----------------"
echo "MSR+RBT(CYC) Systematic Test: Passed!"
echo "-----------------"

# test MSR encoding, collection
./encoder_pmc 20 10 19 8 MSR alice
./collector_pmc $META `seq 10 19`
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MSR Data-Collection Test: Passed!"
echo "-----------------"

# test MSR repair
./encoder_pmc 20 10 19 8 MSR alice
rm coded_files/alice_node0_sym*
./repair_pmc help 0 $META `seq 1 19`
./repair_pmc fix 0 $META `seq 1 19`
./collector_pmc coded_files/alice_metadata
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MSR Repair Test: Passed!"
echo "-----------------"

# test MSR+RBT encoding, collection
./encoder_pmc 20 10 19 8 MSR -rbt CYC alice
./collector_pmc coded_files/alice_metadata `seq 10 19`
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MSR+RBT Data-Collection Test: Passed!"
echo "-----------------"

# test MSR+RBT encoding, collection
./encoder_pmc 20 10 19 8 MSR -rbt SYS alice
./collector_pmc coded_files/alice_metadata 0 1 2 `seq 10 16`
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MSR+RBT(SYS) Data-Collection Test: Passed!"
echo "-----------------"

# test MSR+RBT repair
./encoder_pmc 20 10 19 8 MSR -rbt CYC alice
rm coded_files/alice_node0_sym*
./repair_pmc help 0 $META `seq 1 19`
./repair_pmc fix 0 $META `seq 1 19`
./collector_pmc coded_files/alice_metadata
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MSR+RBT Repair Test: Passed!"
echo "-----------------"


# test MSR+SYSP RBT encoding, collection
./encoder_pmc 20 10 19 8 MSR -rbt SYS18 alice
./collector_pmc coded_files/alice_metadata
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MSR+SYSP RBT Data-Collection Test: Passed!"
echo "-----------------"

# test MSR+SYSP RBT repair
./encoder_pmc 20 10 19 8 MSR -rbt SYS18 alice
rm coded_files/alice_node0_sym*
./repair_pmc help 0 $META `seq 1 19`
./repair_pmc fix 0 $META `seq 1 19`
./collector_pmc coded_files/alice_metadata
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MSR+SYSP RBT Repair Test: Passed!"
echo "-----------------"

#read -p "pause\n"



# test MBR encoding, collection
./encoder_pmc 20 10 15 8 MBR alice
./collector_pmc coded_files/alice_metadata
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MBR Data-Collection Test: Passed!"
echo "-----------------"

# test MBR repair
./encoder_pmc 20 10 15 8 MBR alice
rm coded_files/alice_node0_sym*
./repair_pmc help 0 $META `seq 1 15`
./repair_pmc fix 0 $META `seq 1 15`
./collector_pmc coded_files/alice_metadata
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MBR Repair Test: Passed!"
echo "-----------------"

# test MBR+RBT encoding, collection
./encoder_pmc 20 10 15 8 MBR -rbt CYC alice
./collector_pmc coded_files/alice_metadata
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MBR+RBT Data-Collection Test: Passed!"
echo "-----------------"

# test MBR+RBT repair
./encoder_pmc 20 10 15 8 MBR -rbt CYC alice
rm coded_files/alice_node0_sym*
./repair_pmc help 0 $META `seq 1 15`
./repair_pmc fix 0 $META `seq 1 15`
./collector_pmc coded_files/alice_metadata
cmp ./alice $TESTDATA/alice_golden
echo "-----------------"
echo "MBR+RBT Repair Test: Passed!"
echo "-----------------"


echo "-----------------"
echo "ALL TESTS PASSED!"
echo "-----------------"
