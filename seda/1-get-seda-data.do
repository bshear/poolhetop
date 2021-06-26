*26jul2018
*author: benjamin shear
*make seda subset for pooled hetop HLM models

clear all
set more off
version 14.2

cd "/Users/bshear/Dropbox/GitHub/poolhetop/"

* Load SEDA v2.1 data.
* Please complete data use agreement at this link before loading:
* https://edopportunity.org/get-the-data/

use "https://stacks.stanford.edu/file/druid:db586ns4974/SEDA_geodist_long_State_v21.dta"

keep leaidC stateabb grade year subject totgyb_all ///
	sd_all sd_all_se mn_all mn_all_se totgyb_all

* drop districts with missing SD estimates
drop if sd_all == .

* drop dgyb obs with fewer than 50 students
drop if totgyb_all < 50

* drop FL 2009 / TX 2012 / SC 2011 / NM 2015
drop if stateabb == "FL" & year == 2009
drop if stateabb == "NM" & year == 2015
drop if stateabb == "SC" & year == 2011
drop if stateabb == "TX" & year == 2012

* drop states with only one grade-year observation (DC and HI)
drop if inlist(stateabb, "DC", "HI")
bys stateabb grade year subject : g n1 = _N
count if n1==1

egen tagstate = tag(stateabb subject grade year)
table grade year subject if tagstate

g gprime = log(sd_all)
g gprime_se = sqrt(sd_all_se^2*(1/(sd_all^2)))

keep stateabb leaidC grade year subject totgyb_all ///
	mn_all mn_all_se sd_all sd_all_se gprime gprime_se

compress

save "seda/SEDA_geodist_long_State_v21_sub.dta" , replace

