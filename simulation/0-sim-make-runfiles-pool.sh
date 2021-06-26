## make the run files

for nc in "10" "25" "50" "100" "200"
do

cat <<EOF > sim-run-pool-n${nc}.sh

n=${nc}

for reps in "1_250" "251_500" "501_750" "751_1000"
do
for sd in "equal"
do
for ct in "skewed" "mixed" "mid" "wide"
do

filename=pool-sim-a25-n_\${n}-sd_\${sd}-\${ct}-rep\${reps}

statamp -e do \$filename.do

done
done
done

EOF

done
