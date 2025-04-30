#!/bin/bash 

for i in {11..15}; do
  mkdir -p "run_$i"

  cp nvt.in "run_$i/"

  sed -i "s|variable run_idx equal 1|variable run_idx equal $i|" "run_$i/nvt.in"

  sed -i "s|variable rng     equal 389218|variable rng     equal $RANDOM|" "run_$i/nvt.in"
done
