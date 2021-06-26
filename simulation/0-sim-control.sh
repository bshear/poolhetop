# control file to run simulations for shear & reardon (2020)
# updated 05nov2019 - run trend and equal SD conditions separately
# updated 31may2019
# author: benjamin shear

#

# required files in directory:

## hetop2_lfw.ado
## hetop2_p.ado
## hetop2.ado

## 0-sim-programs.do
## 0-sim-make-counts.do
## 0-sim-base.do
## 0-sim-base-pool.do

## 0-sim-control.sh
## 0-sim-make-dofiles.sh
## 0-sim-make-runfiles.sh
## 0-sim-make-dofiles-pool.sh
## 0-sim-make-runfiles-pool.sh

# update directory
cd "/Users/bshear/Dropbox/GitHub/poolhetop/simulation/"

# simulate count data
statamp -e do 0-sim-make-counts.do

# make the do-files to run main simulation
bash 0-sim-make-dofiles.sh

# make the command files to run main simulation
bash 0-sim-make-runfiles.sh

# make the do-files to run the variable number of pooled datasets
bash 0-sim-make-dofiles-pool.sh

# make the command files to run the variable number of pooled datasets
bash 0-sim-make-runfiles-pool.sh


# PART 1: main simulation results

# run main simulation (warning this takes substantial time to execute; 
# potentially days)

# trend SD conditions

# run 25, 50, 200 conditions in parallel
for n in "25" "50" "200"
do
bash sim-run-n${n}-trend.sh &
done

# run n 10, 100 in parallel
for n in "10" "100"
do
bash sim-run-n${n}-trend.sh &
done

# equal SD conditions

# run 25, 50, 200 conditions in parallel
for n in "10" "25" "50" "100" "200"
do
bash sim-run-n${n}-equal.sh &
done

# run n 10, 100 in parallel
for n in "10" "100"
do
bash sim-run-n${n}-equal.sh &
done


# PART 2: fit pooled models with varying numbers of datasets

# run just pooled hetop model on different numbers of data sets

# n 25, 50, 200conditions
for n in "25" "50" "200"
do
bash sim-run-pool-n${n}.sh &
done

# n 10 100 conditions
for n in "10" "100"
do
bash sim-run-pool-n${n}.sh &
done


