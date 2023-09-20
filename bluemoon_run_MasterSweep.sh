#!/bin/bash
# this is a scrupt for running lots of tests in a batch. Its a bit of a pain but the way to do it is to create either run or test file with a number (run_DRL_{num}.py) corresponding to the numbers in the loop here for each test you want to run
for ((i=130 ; i<154 ; i++ )); do
        echo $i
        sbatch --job-name=Test_Trained_Agent_${i} --output=%x.out  bluemoon_run_ActionSweepSubmit.sh $i
done


