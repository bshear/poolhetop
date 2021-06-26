## make the run files

for sdc in "trend" "equal"
do
for nc in "10" "25" "50" "100" "200"
do

cat <<EOF > sim-run-n${nc}-${sdc}.sh

n=${nc}
sd=${sdc}

## "1_250" 

for reps in "251_500" "501_750" "751_1000"
do
for ct in "skewed" "mixed" "mid" "wide"
do

filename=sim-a25-n_\${n}-sd_\${sd}-\${ct}-rep\${reps}

statamp -e do \$filename.do

done
done

EOF

done
done
