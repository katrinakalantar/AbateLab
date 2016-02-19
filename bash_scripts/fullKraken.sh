#!/bin/bash

cd $1
for f in *.f*q
do
echo $f
~/scripts/runKraken.sh $f
done

