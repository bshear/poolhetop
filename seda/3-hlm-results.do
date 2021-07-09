*10jun2019
*31may2019
*28mar2019
*02aug2018
*author: benjamin shear

* run entire file to produce the data needed to paste into Excel and make
* Table 1 of paper.
* run after step 2 that fits HLM models
* 

clear all
set more off
version 14.2

global ests_folder	"seda"
global outname		"hlm-table.dta"

use "seda/SEDA_geodist_long_State_v21_sub.dta", clear

* ---------------------------------------------------------------------------- *
* load estimates & record sample counts

levelsof stateabb , local(states)

tempname memhold
tempfile results
postfile `memhold' ///
	str2(stateabb) str4(sub) ///
	ndist ncells ny ng n1 n3 nunique meanobs minobs maxobs ///
	omega2 t00 t11 t22 t10 t20 t21 dev pdev df mod ///
	using "`results'"

qui foreach s in `states' {
	
	qui levelsof subject if stateabb == "`s'" , local(subjects)
	
	noi di "." _c
	
	foreach sub in `subjects' {
				
		preserve
			
			unique leaidC if stateabb == "`s'"
			local nunique = r(unique)
			
			keep if stateabb == "`s'" & subject == "`sub'"
			
			sort leaidC subject grade year
			
			by leaidC : g distobs = _N
			
			egen tag = tag(leaidC)
			qui su tag
			local ndist = r(sum)
			
			qui su distobs if tag
			local meanobs = r(mean)
			local minobs = r(min)
			local maxobs = r(max)
			
			drop tag
			
			egen tag = tag(grade year)
			qui su tag
			local ncells = r(sum)
			drop tag

			egen tag = tag(grade)
			qui su tag
			local ngrades = r(sum)
			drop tag

			egen tag = tag(year)
			qui su tag
			local nyears = r(sum)
			drop tag
			
			forv m = 1/3 {

				cap est use "${ests_folder}/hlm`sub'est`s'`m'.ster"
				
				if _rc == 0 {
				
					local tau00		= e(tau00)
					local omega2	= e(s2)
					local df		= e(df)
					local dev		= e(dev)
					local pdev		= e(p_dev)
					local n1		= e(N_1)
					local n3		= e(N_3)
					foreach v in 11 22 10 20 21 {
						local tau`v' = .
					}
					if `m' == 2 {
						foreach v in 11 10 {
							local tau`v' = e(tau`v')
						}
					}
					if `m' == 3 {
						foreach v in 11 10 22 20 21 {
							local tau`v' = e(tau`v')
						}
					}
					
					local posts
					foreach v in ndist ncells nyears ngrades n1 n3 nunique ///
						meanobs minobs maxobs ///
						omega2 tau00 tau11 tau22 tau10 tau20 tau21 dev pdev df m {
						local posts "`posts' (``v'')"
					}
					
					post `memhold' ("`s'") ("`sub'") ///
					`posts'
				}
			
			}
		
		restore
	
	}

}

postclose `memhold'

use `results' , clear

bys stateabb sub : g nmods=_N


* ---------------------------------------------------------------------------- *
* add the null model omega^2

g omega2_null = .
g mn_var = .

qui levelsof stateabb

foreach s in `r(levels)' {
	
	qui levelsof sub
	
	foreach b in `r(levels)' {
		
		est use "${ests_folder}/hlm`b'est`s'0.ster"
		
		qui replace omega2_null = e(s2) if stateabb=="`s'" & sub=="`b'"
		qui replace mn_var = e(mn_evar) if stateabb=="`s'" & sub=="`b'"
		
	}
	
}

preserve
	g ratio = omega2_null/mn_var if mod==1
	g ratio2 = (t00+omega2)/omega2_null if mod==3
	su ratio*
restore


* ---------------------------------------------------------------------------- *
* subset to state-sub with N>=50 districts

keep if ndist > 49
table stateabb sub // all states have all 3 models in both subjects here


* ---------------------------------------------------------------------------- *
* final sample sizes/counts

* confirm that number of obs and ndist remain constant across models in states
bys stateabb sub: egen v = sd(n1)
su v
drop v
bys stateabb sub: egen v = sd(n3)
su v
drop v

* total observations - 620,588
su n1 if mod==1
di r(sum)

* number of unique districts - 9,266
su nunique if mod==1 & sub == "math"
di r(sum)
su nunique if mod==1 & sub == "ela"
di r(sum)

	* slight variation in # districts per subject
	su n3 if sub=="math" & mod==1
	di r(sum)

	su n3 if sub=="ela" & mod==1
	di r(sum)

	* compare to counts and summarize gprime in long data final sample
	preserve
		keep if mod==1
		keep stateabb sub
		rename sub subject
		merge 1:m stateabb subject using "seda/SEDA_geodist_long_State_v21_sub.dta"
		keep if _merge == 3
		unique leaidC
		su
		table subject , c(m gprime sd gprime n gprime)
	restore

* ---------------------------------------------------------------------------- *
* model result calculations

g r = t00/(t00+omega2)
table mod , c(m r)
g lor = exp(-1.96*sqrt(omega2))
g hir = exp(1.96*sqrt(omega2))
g lor_homop = exp(-1.96*sqrt(t00+omega2))
g hir_homop = exp(1.96*sqrt(t00+omega2))
su lo* hi* if mod==3

g deltae = .
g sig2 = 0
g sig3 = 0
g t00_1 = .
g omega2_3 = .

qui levelsof stateabb , local(states)

qui foreach s in `states' {
	
	qui levelsof sub if stateabb == "`s'" , local(subs)
	
	foreach v in `subs' {
		
		* compute delta omega
		qui su omega2 if stateabb == "`s'" & sub == "`v'" & mod == 1
			local a1 = r(mean)
		qui su omega2 if stateabb == "`s'" & sub == "`v'" & mod == 3
			local a2 = r(mean)
		replace deltae = (`a1'-`a2')/`a1' if stateabb == "`s'" & sub == "`v'"
		replace omega2_3 = `a2' if sub == "`v'" & stateabb == "`s'"
		
		* significance tests for slope variance
		replace sig2 = 1 if sub == "`v'" & stateabb == "`s'" & mod == 2 & pdev < 0.01
		replace sig3 = 1 if sub == "`v'" & stateabb == "`s'" & mod == 3 & pdev < 0.01
		
		su t00 if sub == "`v'" & stateabb == "`s'" & mod == 1
		replace t00_1 = r(mean) if sub == "`v'" & stateabb == "`s'"
	
	}

}

g r2=t00/(t00+omega2)

g sqrt_t00=sqrt(t00)
g sqrt_omega2=sqrt(omega2)

* calculate RE covariances from correlations and variances
g rho10=t10
g rho20=t20
g rho21=t21
replace t10=rho10*sqrt(t00*t11)
replace t20=rho20*sqrt(t00*t22)
replace t21=rho21*sqrt(t11*t22)

g rr = omega2/t00

* ---------------------------------------------------------------------------- *
* table for paper

table sig2 if mod == 2
table sig3 if mod == 3

g r3 = omega2/omega2_null
g v00_omega21 = t00+omega2

tempname memhold
tempfile results
postfile `memhold' str12(statistic) str4(sub) min mean max n using `results'

foreach s in math ela {

	qui su omega2 if mod == 1 & sub == "`s'" , d
	post `memhold' ("omega2_1") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')

	qui su t00 if mod == 1 & sub == "`s'" , d
	post `memhold' ("v00") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')
	
	qui su rr if mod == 1 & sub == "`s'" , d
	post `memhold' ("ratio_err") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')
	
	qui su r2 if mod == 1 & sub == "`s'" , d
	post `memhold' ("ratio_bw") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')
	
	qui su omega2 if mod == 3 & sub == "`s'" , d
	post `memhold' ("omega2_2") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')
	
	qui su t00 if mod == 3 & sub == "`s'" , d
	post `memhold' ("t00") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')
	
	qui su t11 if mod == 3 & sub == "`s'" , d
	post `memhold' ("t11") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')

	qui su t22 if mod == 3 & sub == "`s'" , d
	post `memhold' ("t22") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')

	qui su deltae if mod == 3 & sub == "`s'" , d
	post `memhold' ("delta_e") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')	

	qui su omega2_null if mod == 1 & sub == "`s'" , d
	post `memhold' ("omega_null") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')	

	qui su v00_omega21 if mod == 1 & sub == "`s'" , d
	post `memhold' ("v00_omega21") ("`s'") (`r(min)') (`r(mean)') (`r(max)') (`r(N)')	
	
}

postclose `memhold'

preserve
use `results' , clear
list
keep if inlist(statistic, "omega2_1", "omega2_2", "v00", "ratio_bw", "t00", "t11", "t22", "delta_e")
reshape wide min mean max n , i(statistic) j(sub) string
order n* , last
g order = 1
replace order = 2 if statistic=="v00"
replace order = 3 if statistic=="ratio_bw"
replace order = 4 if statistic=="omega2_2"
replace order = 5 if statistic=="t00"
replace order = 6 if statistic=="t11"
replace order = 7 if statistic=="t22"
replace order = 8 if statistic=="delta_e"
order order, last
sort order
list , clean noobs
restore

preserve
sort r2
table mod sub, c(m r2)
count if r2<0.5 & mod==1
*list sub r2 if mod==1 , clean noobs
restore


