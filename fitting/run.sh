#!/bin/bash

model=1
deep=1

for (( model=12; model<=20; model++ ))
do
  R CMD BATCH "--args model=$model deep=$deep" run_fit.R &
done
