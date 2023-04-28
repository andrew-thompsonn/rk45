#!/bin/bash

clearances=(0 10 100 1000 5000 10000 50000 100000)
accuracy=0.5

for clearance in ${clearances[*]}; do
    ./three_body 1 $clearance $accuracy
    ./three_body 2 $clearance $accuracy
done