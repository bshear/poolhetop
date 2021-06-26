
#foldername=h2-sim2-zeros
matsize=800

for reps in "1_250" "251_500" "501_750" "751_1000"
do
for nc in "10" "25" "50" "100" "200"
do
for sd in "equal" "trend"
do
for ct in "mixed" "skewed" "mid" "wide"
do

filename=sim-a25-n_$nc-sd_$sd-$ct-rep${reps}

cat <<EOF > $filename.do

* ---------------------------------------------------------------------------- *
* condition parameters

clear all
set more off
macro drop _all
set type double
set matsize $matsize

set seed $RANDOM
global ncond	"$nc"
global cutname	"$ct"
global sdcond	"$sd"
global reps		"$reps"

EOF

cat "0-sim-base.do" >> $filename.do 	# now, add the base script at the end

done
done
done
done




