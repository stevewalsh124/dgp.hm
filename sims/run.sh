#!/bin/bash

func=2
setting=2

for (( func=1; func<=2; func++ ))
do
  for ((setting=1; setting<=4; setting++ ))
  do
    for (( seed=1; seed<=20; seed++ ))
    do
    R CMD BATCH "--args func=$func setting=$setting seed=$seed" run_sim.R &
    wait
    done
  done
done
