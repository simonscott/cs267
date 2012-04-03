#!/bin/bash

make clean
make knapsack
qsub job-knapsack
watch -n 1 qstat -u simons
cat knapsack.stdout
