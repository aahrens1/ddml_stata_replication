#!/bin/sh
rm -rf log
rm -rf out
mkdir log
mkdir out
module load gcc/8.2.0 python/3.9.9 stata/16
Seed="7257 9484 7433 7931 4195 722 5534 5504 1739 306 2612 6840 9453"
Dgp="2 3 4 5 7"
for val in $Seed; do
for dgp in $Dgp; do
sbatch -n 4 --time=5-0 -e error.txt -o output.txt --mem-per-cpu=4000 --wrap="stata-se -b do sim_SJ $val 1000 $dgp 5"
sbatch -n 4 --time=5-0 -e error.txt -o output.txt --mem-per-cpu=4000 --wrap="stata-se -b do sim_SJ $val 100 $dgp 20"
done
done
