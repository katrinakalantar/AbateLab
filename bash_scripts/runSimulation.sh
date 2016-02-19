#!/bin/bash

echo !!Initiating Simulation Script!!
date 

output_dir="~/simulation/$1"

nohup ~/simulation/scripts/simulation2_v1.sh $output_dir &

echo !!Finished Simulation Script!!
date


