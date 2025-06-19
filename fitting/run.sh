#!/bin/bash

model=1

for (( model=0; model<=116; model++ ))
do
  R CMD BATCH "--args model=$model" run_fit.R &
  wait
done
