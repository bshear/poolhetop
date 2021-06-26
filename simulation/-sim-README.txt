
README file for Stata simulation files to accompany Shear and Reardon (2020) "Using Pooled Heteroskedastic Ordered Probit Models to Improve Small-Sample Estimates of Latent Test Score Distributions."

Required Files:

hetop2_lfw.ado
hetop2_p.ado
hetop2.ado

0-sim-programs.do
0-sim-make-counts.do
0-sim-base.do
0-sim-base-pool.do

0-sim-control.sh
0-sim-make-dofiles.sh
0-sim-make-runfiles.sh
0-sim-make-dofiles-pool.sh ---- pool files are just to do different number of pooled datasets
0-sim-make-runfiles-pool.sh

1-sim-rerun-nonconverged.do
2-sim-analysis-v5.do
3-sim-figs-tables.do

All files begin in a single directory.
Entire simulation can be run from 0-sim-control.sh
This will create necessary data, do-files, and results. 
The process has some parallelization, but still takes multiple days to run.
Run 1-sim-rerun-nonconverged.do to re-fit models that did not converge with higher max iterations.
Run 2-sim-analysis.do to create summary results files.
Run 3-sim-figs-tables.do to create figures and tables for paper.
