*10jun2020: clean up and remove extra code for distribution
*31dec2019
*ben shear

/*

Fit different HETOP models to example data and create table of summary
results.

*/

* ---------------------------------------------------------------------------- *
* setup

	
	clear all
	version 14.2
	set type double
	set matsize 10000

	local date: display %td_CCYY-NN-DD date(c(current_date), "DMY")
	di "`date'"
	
	global run_ests 	= 1		// estimate models
	global make_table	= 1		// summarize results
	

* ---------------------------------------------------------------------------- *
* pooled hetop functions

	
	do ado/hetop2.ado
	do ado/hetop2_lfw.ado
	do ado/hetop2_p.ado


* ---------------------------------------------------------------------------- *
* estimation

if ${run_ests} == 1 {

	
	use "example/example_data.dta"
	
	
	* ------------------------------------------------------------------------ *
	* use separate HETOP models in each grade (which allows for PHOP)

	
		preserve
		drop _all
		tempfile res
		save `res' , emptyok
		restore

		forv g = 3/8 {
			preserve
				keep if grade==`g'
				g phopvar = num0>1
				hetop numid level , numcats(4) phop(phopvar, mean) ///
					identify(cuts) save(_hetop, star se)
					
				if e(converged) != 1 {
				
					noi di in red "warning model did not converge"
					
					stop
				
				}	
					
				g loglik_hetop = e(ll)
				g icc_hetop = e(icchat)
				g icc_hetop_se = e(icchat_var)
				replace icc_hetop_se = icc_hetop_se^.5
				replace loglik_hetop = . if _n !=1
				est save example/mod_hetop`g' , replace
				est store hetop`g'
				append using `res'
				save `res' , replace
			restore
		}

		use `res' , clear

		
		* save hetop cutscores
		* these are available in star metric
		
		local m hetop
		forv k = 1/3 {	
			g c`k'_`m' = .
		}

		forv g = 3/8 {
			
			est restore hetop`g'
			mat c = e(cstar)
			qui forv k = 1/3 {
				replace c`k'_`m' = c[`k',1] if grade == `g'
			}

		}
		
		
	* ------------------------------------------------------------------------ *	
	* homop model within grade

	
		hetop2 level , k(3) ///
			grpmn(numid_grd) datamn(grdc) ///
			grpsd(one) datasd(one) ///
			datact(grdc)
			
		if e(converged) != 1 {
		
			noi di in red "warning model did not converge"
			
			stop
		
		}
		
		g loglik_homop = e(ll)
		
		est save example/mod_homop , replace
		est store homop
		
		g icc_homop = .
		g icc_homop_se = .
		qui levelsof grdc , local(grdc_levs)
		predict homop , star se
		foreach g in `grdc_levs' {
			replace icc_homop = r(icchat`g') if grdc == `g'
			replace icc_homop_se = `r(icchat`g'_var)'^.5 if grdc == `g'
		}
				
		
		* estimate fixed cuts model
		hetop2 level , k(3) ///
			grpmn(numid_grd) datamn(grdc) ///
			grpsd(one) datasd(one) ///
			datact(one)

		if e(converged) != 1 {
		
			noi di in red "warning model did not converge"
			
			stop
		
		}			
			
		est save example/mod_homop_fixcuts , replace
		
		
	* ------------------------------------------------------------------------ *
	* pooled hetop 

	
		hetop2 level , k(3) ///
			grpmn(numid_grd) datamn(grdc) ///
			grpsd(numid) datasd(one) ///
			datact(grdc) startvals

		if e(converged) != 1 {
		
			noi di in red "warning model did not converge"
			
			stop
		
		}
			
		g loglik_phop = e(ll)
			
		* add estimated intercepts
		g b0_phop = .
		g b0_phop_se = .
		
		qui levelsof numid , local(intlevels)
		mat b = e(b)
		mat v = e(V)
		qui foreach i in `intlevels' {
			local x = b[1,colnumb(b,"lnsigma:`i'.numid")]
			replace b0_phop = `x' if numid == `i'
			local v = v[rownumb(v,"lnsigma:`i'.numid"), ///
				colnumb(v,"lnsigma:`i'.numid")]
			replace b0_phop_se = sqrt(`v') if numid == `i'
		}
		
		est save example/mod_phop , replace
		est store phop
		
		g icc_phop = .
		g icc_phop_se = .
		qui levelsof grdc , local(grdc_levs)
		predict phop , star se
		foreach g in `grdc_levs' {
			replace icc_phop = r(icchat`g') if grdc == `g'
			replace icc_phop_se = `r(icchat`g'_var)'^.5 if grdc == `g'
		}

		
		* estimate constant cuts version
		hetop2 level , k(3) ///
			grpmn(numid_grd) datamn(grdc) ///
			grpsd(numid) datasd(one) ///
			datact(one) startvals difficult iter(100)

		if e(converged) != 1 {
		
			noi di in red "warning model did not converge"
			
			stop
		
		}			
			
		est save example/mod_phop_fixcuts , replace
				
				
	* ------------------------------------------------------------------------ *		
	* trend hetop

	
		hetop2 level , k(3)	///
			grpmn(numid_grd) datamn(grdc)	///
			grpsd(numid) datasd(one) bsd(grdc)	///
			datact(grdc) startvals
			
		if e(converged) != 1 {
		
			noi di in red "warning model did not converge"
			
			stop
		
		}
			
		* add estimated intercepts and slopes
		
		g loglik_trend = e(ll)
		
		g b0_trend = .
		g b0_trend_se = .
		g b1_trend = .
		g b1_trend_se = .
		
		qui levelsof numid , local(intlevels)
		
		mat b = e(b)
		mat v = e(V)
		
		qui foreach i in `intlevels' {

			local x = b[1,colnumb(b,"lnsigma:`i'.numid")]
			replace b0_trend = `x' if numid == `i'
			local v = v[rownumb(v,"lnsigma:`i'.numid"), ///
				colnumb(v,"lnsigma:`i'.numid")]
			replace b0_trend_se = sqrt(`v') if numid == `i'
			
			local x = b[1,colnumb(b,"lnsigma:`i'.numid#c.grdc")]
			replace b1_trend = `x' if numid == `i'
			local v = v[rownumb(v,"lnsigma:`i'.numid#c.grdc"), ///
				colnumb(v,"lnsigma:`i'.numid#c.grdc")]
			replace b1_trend_se = sqrt(`v') if numid == `i'
			
		}
		
		est save example/mod_trend , replace
		est store trend
		
		g icc_trend = .
		g icc_trend_se = .
		qui levelsof grdc , local(grdc_levs)
		predict trend , star se
		foreach g in `grdc_levs' {
			replace icc_trend = r(icchat`g') if grdc == `g'
			replace icc_trend_se = `r(icchat`g'_var)'^.5 if grdc == `g'
		}

		
		* estimate constant cuts version
		
		hetop2 level , k(3)	///
			grpmn(numid_grd) datamn(grdc)	///
			grpsd(numid) datasd(one) bsd(grdc)	///
			datact(one) startvals difficult iter(100)
			
		if e(converged) != 1 {
		
			noi di in red "warning model did not converge"
			
			stop
		
		}	
		
		est save example/mod_trend_fixcuts , replace
		

	* ------------------------------------------------------------------------ *		
	* add estimated cut scores and prime metric estimates 
	* to compute predicted proportions
	
	
		foreach m in phop trend homop {
			
			est restore `m'
			
			mat b = e(b)
			
			qui forv k = 1/3 {
				
				g c`k'_`m' = .
				
				forv g = 0/5 {
					
					local i = 1
					local j = colnumb(b,"cut`k':`g'.grdc")
					local x = b[`i',`j']
					replace c`k'_`m' = `x' if grdc == `g'
					
				}
				
			}
			
			g mprime_`m' = .
			
			qui levelsof numid_grd
			qui foreach r in `r(levels)' {
				local i = 1
				local j = colnumb(b,"means:`r'.numid_grd")
				local x = b[`i',`j']
				replace mprime_`m' = `x' if numid_grd==`r'
			}
			
			
			if "`m'" == "phop" 		g sprime_`m' = exp(b0_phop)
			if "`m'" == "trend" 	g sprime_`m' = exp(b0_trend + b1_trend*grdc)
			if "`m'" == "homop" 	g sprime_`m' = 1
			
		}

	
	* ------------------------------------------------------------------------ *
	* compute fit stats
	
		
		* observed proportions in each cell
		
		
		forv k = 1/4 {
			g p`k'_obs = level`k'/nobs
		}

		
		* predicted proportions in each cell
		
		foreach m in hetop trend phop homop {

			if "`m'" == "hetop" local v star
			if "`m'" != "hetop" local v prime

			g p1_`m' = normal((c1_`m'-m`v'_`m')/s`v'_`m')
			g p2_`m' = normal((c2_`m'-m`v'_`m')/s`v'_`m')-p1_`m'
			g p3_`m' = normal((c3_`m'-m`v'_`m')/s`v'_`m')-normal((c2_`m'-m`v'_`m')/s`v'_`m')
			g p4_`m' = 1-normal((c3_`m'-m`v'_`m')/s`v'_`m')
			
		}
		
		
		* average absolute differences in proportions for each school-grade
		
		foreach m in hetop trend phop homop {

			g pdif_`m' = ///
				(abs(p1_obs-p1_`m') + ///
				abs(p2_obs-p2_`m') + ///
				abs(p3_obs-p3_`m') + ///
				abs(p4_obs-p4_`m')) / 4
			
			
		}
	
		
		* chi-square GOF values
		* compute chi squares as Sum[ ((O_k-E_k)^2) / E_k ]
		
		foreach m in hetop trend phop homop {
					
			* chi-square GoF for each row (school-grade unit)
			g x2_`m' = ///
				((level1-(p1_`m'*nobs))^2)/((p1_`m'*nobs)) + ///
				((level2-(p2_`m'*nobs))^2)/((p2_`m'*nobs)) + ///
				((level3-(p3_`m'*nobs))^2)/((p3_`m'*nobs)) + ///
				((level4-(p4_`m'*nobs))^2)/((p4_`m'*nobs))
							
		}
				
	
	* ------------------------------------------------------------------------ *
	* confirm convergence of models
	
		foreach m in trend phop homop {
			est restore `m'
			if e(converged) != 1 {
				noi di in red "`m' not converged"
				stop
			}
		}

		foreach m in trend phop homop {
			est use example/mod_`m'_fixcuts
			if e(converged) != 1 {
				noi di in red "`m' fixcuts not converged"
				stop
			}
		}			
	
	
	* ------------------------------------------------------------------------ *
	* save
	
		compress
		
		save "example/example_estimates.dta" , replace
	
	
}


* ---------------------------------------------------------------------------- *
* summary table and LRTs

if $make_table == 1 {
	
	
	* load stored model results

	foreach m in homop phop trend {
		est use example/mod_`m'.ster
		est store `m'
		est use example/mod_`m'_fixcuts.ster
		est store `m'_fixcuts
	}	
	
	use "example/example_estimates.dta" 
	
	
	* add number of parameters and constraints
	
	foreach m in phop trend homop {
		est restore `m'
		g npar_`m' = e(k)
		matrix C = e(Cns)
		local c = rowsof(C)
		g ncns_`m' = `c'
	}	
	
	
	* marginal frequencies
	
	preserve
		collapse (sum) level1 level2 level3 level4 , by(grdc)
		egen N = rowtotal(level1 level2 level3 level4)
		forv i = 1/4 {
			replace level`i' = level`i'/N
		}
		sort grdc
		list , clean noobs
	restore
	
	
	* LRT's
		
		* homop model with fixed vs free cuts over grades
		* homop fixcuts df=(3 cuts)*(5 grades)=15
			
			lrtest homop homop_fixcuts
		
		* phop model with fixed vs free cuts over grades
			
			lrtest phop phop_fixcuts , stats dir
			
		* trend model with fixed vs free cuts over grades
		
			lrtest trend trend_fixcuts
		
		* homop v phop model
		
			lrtest homop phop
		
		* phop v trend model
		
			lrtest phop trend
		
	
	* table
	
	preserve
		
		egen tag = tag(numid)
		replace b0_phop = . if tag==0
		replace b0_trend = . if tag==0
		replace b1_trend = . if tag==0
		

		collapse ///
			(mean) ///
			loglik_phop loglik_trend loglik_homop ///
			pdif_phop pdif_trend pdif_homop ///
			npar_phop npar_trend npar_homop ///
			ncns_phop ncns_trend ncns_homop ///
			b0_m_phop=b0_phop ///
			b0_m_trend=b0_trend ///
			b1_m_trend=b1_trend ///
			sstar_m_phop=sstar_phop ///
			sstar_m_trend=sstar_trend ///
			sstar_m_homop=sstar_homop ///
			mstar_m_phop=mstar_phop ///
			mstar_m_trend=mstar_trend ///
			mstar_m_homop=mstar_homop ///
			(sum) ///
			x2_phop x2_trend x2_homop ///
			(sd) ///
			b0_sd_phop=b0_phop ///
			b0_sd_trend=b0_trend ///
			b1_sd_trend=b1_trend ///
			sstar_sd_phop=sstar_phop ///
			sstar_sd_trend=sstar_trend ///
			sstar_sd_homop=sstar_homop ///
			mstar_sd_phop=mstar_phop ///
			mstar_sd_trend=mstar_trend ///
			mstar_sd_homop=mstar_homop ///
			(min) ///
			b0_min_phop=b0_phop ///
			b0_min_trend=b0_trend ///
			b1_min_trend=b1_trend ///
			sstar_min_homop=sstar_homop ///
			sstar_min_phop=sstar_phop ///
			sstar_min_trend=sstar_trend ///
			mstar_min_homop=mstar_homop ///
			mstar_min_phop=mstar_phop ///
			mstar_min_trend=mstar_trend ///
			(max) ///
			b0_max_phop=b0_phop ///
			b0_max_trend=b0_trend ///
			b1_max_trend=b1_trend ///
			sstar_max_homop=sstar_homop ///
			sstar_max_phop=sstar_phop ///
			sstar_max_trend=sstar_trend ///
			mstar_max_homop=mstar_homop ///
			mstar_max_phop=mstar_phop ///
			mstar_max_trend=mstar_trend ///
			(count) ///
			b0_n_phop=b0_phop ///
			b0_n_trend=b0_trend ///
			sstar_n_phop=sstar_phop ///
			sstar_n_trend=sstar_trend ///
			sstar_n_homop=sstar_homop ///
			, ///
			by(one)
		
		reshape long ///
			loglik_ npar_ ncns_ ///
			x2_ pdif_ ///
			b0_m_ b1_m_ b0_sd_ b1_sd_ ///
			b0_min_ b0_max_ b1_max_ ///
			sstar_sd_ sstar_m_ sstar_min_ sstar_max_ ///
			mstar_sd_ mstar_m_ mstar_min_ mstar_max_ ///
			b0_n_ sstar_n_ mstar_n_ , i(one) j(model) str
		
		order model loglik npar_ ncns_ x2 pdif ///
			b0_m_ b0_sd b0_min_ b0_max_ b1_m_ b1_sd_ b1_min_ b1_max_ ///
			sstar_m_ sstar_sd_ sstar_min_ sstar_max_ ///
			mstar_m_ mstar_sd_ mstar_min_ mstar_max_ ///
			sstar_n_ mstar_n_ 
		
		pause on
		
		browse
		
		pause
		
		pause off
	
	restore		
	
	
}


* ---------------------------------------------------------------------------- *


	

