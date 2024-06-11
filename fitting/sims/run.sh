#!/bin/bash

func=1
setting=1

for (( seed=1; seed<=10; seed++ ))
do
  R CMD BATCH "--args func=$func setting=$setting seed=$seed" run_sim.R &
done
