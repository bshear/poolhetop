* updated: 07jun2019
* author: benjamin shear
* loop through simulation results and try re-running any models that did not 
* converge, using a greater number of max iterations.
* the code will add new estimates (if models converge) to the results file,
* but will also save out the unconverged results file with the suffix "unconverged"
* in the case there is a need to return to the original results.

set more off
clear all
version 14.2

cd "/Users/bshear/Dropbox/GitHub/poolhetop/simulation/"
cap log close
log using "1-sim-rerun-nonconverged_$S_DATE.txt" , replace

global reg_rerun	= 1
global pool_rerun	= 0

* ---------------------------------------------------------------------------- *
* re-run regular sim results

if ${reg_rerun}==1 {

qui foreach n in 10 25 50 100 200 {
	foreach sd in equal trend {
		foreach c in skewed wide mid mixed {
			foreach r in 1_250 251_500 501_750 751_1000 {
				
				local fn "ests-a25-n_`n'-sd_`sd'-`c'-`r'"
				
				use "fullsim/`fn'.dta" , clear
				
				* check for non-convergence
								
				foreach m in fhet fhom phet trend {
					qui count if conv`m' == . | conv`m' == 0
					local `m'N = r(N)
				}
				
				if (`fhetN'+`fhomN'+`phetN'+`trendN')>0 {
					
					save "fullsim/`fn'-unconverged.dta" , replace
					
					noi di "re-running `fn'"
					
					foreach m in fhet fhom phet trend {
						
						cap noi di "model: `m', reps:"
						cap noi table rep grade if conv`m' == . | conv`m' == 0

					}	
						
					foreach m in fhet fhom {
					
					if ``m'N' > 0 {
						
						egen condrepg = group(sdcond ncond cutcond rep grade)
						
						levelsof condrepg if conv`m'==. , local(condrepgs)

						foreach r in `condrepgs' {

							preserve
							
							constraint drop _all
							
							keep if condrepg == `r'
							
							noi di "ests-a25-n_`n'-sd_`sd'-`c'-`r' , REP=`r', `m'"
							
							noi list condrepg sdcond cutcond ncond rep grade id freq? ///
								idrtg mstar_`m' mstar_`m'_se, clean noobs
							
							if "`m'" == "fhet" {
							noi hetop2 freq , k(3) ///
								grpmn(idrtg) datamn(rt) ///
								grpsd(idrtg) datasd(rt) ///
								datact(rt) iter(500)
							}
							
							if "`m'" == "fhom" {
							noi hetop2 freq , k(3) ///
								grpmn(idrtg) datamn(rt) ///
								grpsd(rt)    datasd(rt) ///
								datact(rt) iter(500)
							}
							
							drop mstar_`m'* sstar_`m'* icchat_`m'* 
							predict `m' , star se
							replace atmp`m' = 5
							replace conv`m' = e(converged)
							g double icchat_`m' = r(icchat1) if id == 1
							g double icchat_`m'_se = r(icchat1_var)^.5 if id == 1
							keep condrepg sdcond cutcond ncond rep grade id idrtg *`m'*
							list , clean noobs
							tempfile res`r'
							save `res`r'' , replace
							
							restore
							
							merge 1:1 condrepg sdcond cutcond ncond rep grade ///
								id idrtg using "`res`r''" , ///
								update replace gen(rerun`r')
								
						}
						
						drop condrepg

					} // loop through reps
					} // fhom v fhet
					
					* ------------------------------------------------------------ *
					* phet and trend re-runs
					
					foreach m in phet trend {
					
						if (``m'N' > 0) {
						
							egen condrep = group(sdcond ncond cutcond rep)
						
							levelsof condrep if (conv`m'== . | conv`m' == 0), local(condreps)
							
							noi di `condreps'
							
							local phetmod hetop2 freq , k(3) ///
										grpmn(idrtg) datamn(rt) ///
										grpsd(id)	 datasd(one) ///
										datact(rt)

							local trendmod hetop2 freq , k(3) ///
										grpmn(idrtg) datamn(rt) ///
										grpsd(id)	 datasd(one) bsd(grade) ///
										datact(rt)					
							
							foreach r in `condreps' {
									
								preserve

									constraint drop _all
									
									keep if condrep == `r'
									
									noi di "ests-a25-n_`n'-sd_`sd'-`c'-`r' , REP=`r', `m'"
																								
									list condrep sdcond cutcond ncond rep grade id freq? ///
										idrtg mstar_`m' mstar_`m'_se, clean noobs
									
									``m'mod' iter(500)
									
									drop mstar_`m'* sstar_`m'* icchat_`m'* 
									
									levelsof rt, local(rt_levs)
									predict `m' , star se
									replace atmp`m' = 5
									replace conv`m' = e(converged)
									g double icchat_`m' = .
									g double icchat_`m'_se = .
									
									foreach lev in `rt_levs' {
										replace icchat_`m' = r(icchat`lev')  ///
											if rt == `lev' & id == 1
										replace icchat_`m'_se = r(icchat`lev'_var)^.5  ///
											if rt == `lev' & id == 1
									}			
									
									list , clean noobs
									
									keep condrep sdcond cutcond ncond rep grade id idrtg ///
										?star_`m' ?star_`m'_se icchat_`m' icchat_`m'_se ///
										atmp`m' conv`m'
									
									g rerun=1
									
									tempfile res`r'
									save `res`r'' , replace

								restore

								merge 1:1 condrep sdcond cutcond ncond rep ///
									grade id idrtg using "`res`r''" , ///
									update replace nogen
							
							}	// reps
							
							drop condrep							
						
						}		// check non-converged
						
					}			// models
										
					compress
					
					save "fullsim/`fn'.dta" , replace
				
				}
				
			}
		}
	}
}

}


* ---------------------------------------------------------------------------- *
* re-run pooled sim results

if ${pool_rerun}==1 {

clear all
constraint drop _all
set more off
foreach n in 10 25 50 100 200 {
	foreach sd in equal {
		foreach c in skewed wide mid mixed {
			foreach r in 1_250 251_500 501_750 751_1000 {
				
				local fn "pool-ests-a25-n_`n'-sd_`sd'-`c'-`r'"
				
				use  "poolsim/`fn'.dta"
				
				local m phet
				
				qui count if conv`m' == . | conv`m' == 0
				
				if r(N) > 0 {
					
					noi di "re-running: `fn'"
					
					save "poolsim/`fn'-unconverged.dta" , replace

					local `m'mod hetop2 freq , k(3) ///
								grpmn(idrtg) datamn(rt) ///
								grpsd(id)	 datasd(one) ///
								datact(rt)
										
					egen cond = group(rep npool)
					
					qui levelsof cond if conv`m' == .
					
					foreach r in `r(levels)' {
						
						preserve
						
							keep if cond == `r'
							
							qui su npool
							
							local p = r(mean)
							
							keep if grade < `p'

							list cond sdcond cutcond ncond rep grade id ///
								freq? idrtg ///
								mstar_`m' mstar_`m'_se, clean noobs
	
							`phetmod' iter(500)
	
							drop mstar_`m'* sstar_`m'* icchat_`m'* 
	
							levelsof rt, local(rt_levs)
							predict `m' , star se
							replace conv`m' = e(converged)
							g double icchat_`m' = .
							g double icchat_`m'_se = .
							foreach lev in `rt_levs' {
								replace icchat_`m' = r(icchat`lev')  ///
									if rt == `lev' & id == 1
								replace icchat_`m'_se = r(icchat`lev'_var)^.5  ///
									if rt == `lev' & id == 1
							}			
							
							keep npool sdcond cutcond ncond rep grade id idrtg ///
								?star_`m' ?star_`m'_se icchat_`m' icchat_`m'_se conv`m'	
							
							g rerun = 1
							
							tempfile f`r'
							save `f`r'', replace
							
						restore

						merge 1:1 sdcond cutcond ncond npool rep grade id ///
							idrtg using "`f`r''" , ///
							update replace nogen
						
					}
					
					compress
					
					save "poolsim/`fn'.dta" , replace
					
				}
			}
		}
	}
}

}

cap log close





























