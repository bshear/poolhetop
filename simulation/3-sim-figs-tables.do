* updated: 12dec2019
* author: benjamin shear

/*
make tables and figures from compiled simulation results

accompanies paper shear & reardon (2020) "Using Pooled Heteroskedastic Ordered Probit 
Models to Improve Small-Sample Estimates of Latent Test Score Distributions"

file will create figures for the paper, plus many additional figures for
more outcomes, by conditions, etc. that can be used to further explore results.
figures will be saved to a "figs" subfolder within current directory.
also summarizes convergence, amount of smoothing, true ICC and CV values,
and creates summary tables for additional insepction of results.
*/


* ---------------------------------------------------------------------------- *
* preliminaries

	clear all
	set more off
	version 14.2

	local date: display %td_CCYY-NN-DD date(c(current_date), "DMY")
	di "`date'"

	global all_graphs = 1	// save all graphs or only paper figures?
	
	cd "/Users/bshear/Dropbox/GitHub/poolhetop/simulation/"
	
	* update accordingly if re-run
	use "sim-summary-mean_12 Dec 2019.dta"

	g sdlab = "Equal SD"
	replace sdlab = "Trend SD" if sdcond == "trend"
	bys sdcond cutcond ncond : egen nnum = mean(nrtg)
	destring grade, replace

	
* ---------------------------------------------------------------------------- *
* report smoothing proportions


	preserve
		destring ncond , replace
		collapse (mean) xhet , by(cutcond sdcond ncond)
		table cutcond sdcond ncond , c(m xhet) format(%9.3f)
		su xhet
		su xhet if ncond == 10 & cutcond == "wide"
		su xhet if ncond == 10 & cutcond == "mixed"
	restore
	

/*

----------------------------------------------------------------------------------------
          |                               ncond and sdcond                              
          | ---- 10 ----    ---- 25 ----    ---- 50 ----    ---- 100 ---    ---- 200 ---
  cutcond | equal  trend    equal  trend    equal  trend    equal  trend    equal  trend
----------+-----------------------------------------------------------------------------
      mid | 0.035  0.057    0.000  0.003    0.000  0.000    0.000  0.000    0.000  0.000
    mixed | 0.182  0.218    0.054  0.097    0.012  0.053    0.003  0.028    0.000  0.022
   skewed | 0.107  0.137    0.013  0.030    0.003  0.005    0.000  0.001    0.000  0.000
     wide | 0.412  0.455    0.139  0.224    0.040  0.116    0.007  0.058    0.000  0.040
----------------------------------------------------------------------------------------

    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
        xhet |         40     .063779    .1057359          0   .4545533

*/


* check convergence rates


	preserve
		destring ncond , replace
		collapse (mean) conv* , by(cutcond sdcond ncond)
		su conv*
	restore

/*
    Variable |        Obs        Mean    Std. Dev.       Min        Max
-------------+---------------------------------------------------------
    convphet |         40           1           0          1          1
   convtrend |         40           1           0          1          1
    convfhet |         40           1           0          1          1
    convfhom |         40           1           0          1          1
   convphet1 |         20           1           0          1          1
-------------+---------------------------------------------------------
   convphet2 |         20           1           0          1          1
   convphet3 |         20           1           0          1          1
   convphet4 |         20           1           0          1          1
   convphet5 |         20           1           0          1          1
*/

* ---------------------------------------------------------------------------- *
* true ICC values and CV values for paper text


	preserve
	destring ncond , replace
	bys sdcond: table ncond grade , c(m icctrue)
	su icctrue
	bys sdcond: table ncond grade , c(m cvtrue)
	su cvtrue
	restore

	
* ---------------------------------------------------------------------------- *
* bias, rmse, and inci rates for means and SDs
* tables of results to look at specific values referenced in paper


	table sdcond ncond cutcond , c(m merrfhom)

	set more off
	
	local s serr	

	foreach v in fhom fhet phet trend {
		
		preserve
			
			qui destring ncond , replace
			
			qui collapse `s'* (count) nobs=`s'fhom , by(sdcond ncond cutcond)			
			
			qui keep ncond cutcond sdcond `s'`v' 
			
			qui reshape wide `s'`v' , i(sdcond ncond) j(cutcond) str
			order sdcond ncond *mid *mixed *skewed *wide
			noi di in red "`s'`v'"
			list , clean noobs

		restore

	}
	
* ---------------------------------------------------------------------------- *
* bias, rmse, and inci rates for means and SDs
* figures 1 and 2 in here


	preserve
		
		collapse (mean) nrtg ssqerr* msqerr* merr* serr* mrelerr* srelerr* ///
			iccerr* minci* sinci*  , ///
			by(cutcond nnum sdcond)
		
		replace sdcond = "Constant SD" if sdcond == "equal"
		replace sdcond = "Trend SD" if sdcond == "trend"
			
		replace sincifhom = .
		
		foreach m in fhet fhom phet trend {
			g srmse`m' = sqrt(ssqerr`m')
			g mrmse`m' = sqrt(msqerr`m')
		}
		
	
		
		foreach v in serr srelerr sinci srmse mrmse merr mrmse minci {
			local yline 
			if inlist("`v'", "minci", "sinci") local yline yline(0.95, lp(dash))
			local ylab "`v'"
			if "`v'" == "serr" 		local ylab "Bias"
			if "`v'" == "ssqerr" 	local ylab "Mean Squared Error"
			if "`v'" == "sinci" 	local ylab "CI Coverage Rate"
			if "`v'" == "srmse" 	local ylab "Root Mean Squared Error"
			
			twoway ///
				(connect `v'fhet nrtg) (connect `v'fhom nrtg) ///
				(connect `v'phet nrtg) (connect `v'trend nrtg) , ///
				by(sdcond cutcond, compact note("") r(2)) ///
				ytitle("`ylab'") xtitle("Group Sample Size") ///
				ylab(,angle(0)) `yline' ///
				legend(order(1 "HETOP" 2 "HOMOP" 3 "Pooled" 4 "Trend") r(1)) ///
				scheme(sj) name(`v'cutsline, replace) xsize(6.5) ysize(4)
			
			if ${all_graphs} ==1 graph export "figs/fig-`v'cutsline_`date'.png", replace
			
			if "`v'" == "serr" graph export "figs/FIGURE1_`date'.png", replace
			if "`v'" == "srmse" graph export "figs/FIGURE2_`date'.png", replace
			
		}
		
	restore



* ---------------------------------------------------------------------------- *
* efficiency ratios for PHET
* figure 3


	preserve

		keep if sdcond == "equal"
		
		local v s
		
		forv i = 1/5 {
			g eff`i' = `v'varfhet/`v'varphet`i'
		}
		g eff6 = `v'varfhet/`v'varphet
		
		collapse (mean) eff? , by(nnum ncond cutcond sdcond)
		
		twoway ///
			(connect eff1 nnum) ///
			(connect eff2 nnum) ///
			(connect eff3 nnum) ///
			(connect eff4 nnum) ///
			(connect eff5 nnum) ///
			(connect eff6 nnum), ///
			by(cutcond, r(1) compact note("")) scheme(sj) xsize(6.5) ysize(3) ///
			ytitle("Var(HETOP)/Var(Pooled)") xtitle("Group Sample Size") ///
			ylab(0(1)10 , angle(0)) ///
			legend(order(1 "p=1" 2 "p=2" 3 "p=3" 4 "p=4" 5 "p=5" 6 "p=6") r(1)) ///
			name(effratio_pool_`v', replace)
			
		if ${all_graphs} ==1 graph export "figs/fig-phet-effratio-`v'-`date'.png", replace
		
		graph export "figs/FIGURE3_`date'.png" , replace
			
	restore


* ---------------------------------------------------------------------------- *
* efficiency ratio plot for SDs at different grades
* figure 4


	preserve

		* look at ratio within grades
		
		keep if sdcond == "trend"
		
		g eff = svarfhet/svartrend
		
		sort cutcond grade nnum
			
		* collapse (mean) eff? , by(nnum ncond cutcond sdcond)

		* 1.90909084 3.387096556 5.526315208 5.526315208 3.387096556 1.90909084

		twoway ///
			(connect eff nnum if grade == 0) ///
			(connect eff nnum if grade == 1) ///
			(connect eff nnum if grade == 2) ///
			(connect eff nnum if grade == 3) ///
			(connect eff nnum if grade == 4) ///
			(connect eff nnum if grade == 5) , ///
			yline(1.91,lp(dash)) yline(3.39,lp(dash)) yline(5.53,lp(dash)) ///
			by(cutcond, r(1) compact note("")) scheme(sj) xsize(6.5) ysize(3) ///
			ytitle("Var(HETOP)/Var(Trend)") xtitle("Group Sample Size") ///
			ylab(0(2)10, angle(0)) ///
			legend(order(1 "g=0" 2 "g=1" 3 "g=2" 4 "g=3" 5 "g=4" 6 "g=5") r(1)) ///
			name(effratio_trend2, replace)
			
			if ${all_graphs} == 1 graph export "figs/fig-effratio_trend_grades-`date'.pdf", replace
			
			graph export "figs/FIGURE4_`date'.png", replace

	restore

