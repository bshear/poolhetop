* updated: 12dec2019
* author: benjamin shear
* compile simulation results and create summary .dta files with bias, rmse, etc.
* accompanies paper shear & reardon (2019) "Using Pooled Heteroskedastic Ordered Probit 
* Models to Improve Small-Sample Estimates of Latent Test Score Distributions"

/*
this do-file reads in the model estimates from simulation, appends together
different conditions, then saves out summary files each with
summary statistics. for example one file will have the mean error,
another the min, max, etc. for each outcome. the outcomes are reported at the 
cell level meaning that each condition of the simulation has multiple rows, 
summarizing the labeled outcome across each "grade" within each condition. 
*/


* ---------------------------------------------------------------------------- *
* preliminaries


	clear all
	set more off
	version 14.2

	cd "/Users/bshear/Dropbox/GitHub/poolhetop/simulation/"
	
	global do_log = 0
	
	cap log close
	
	if $do_log == 1 log using "1-sim-analysis-$S_DATE.txt" , text replace

	
* ---------------------------------------------------------------------------- *
* Load pooled HETOP model results
	
	// These are models fit to 2, 3, 4, 5, or 6 grades at a time rather than
	// the full pooled model.

	foreach n in 10 25 50 100 200 {
		foreach sd in equal {
			foreach c in skewed wide mid mixed {
				foreach r in 1_250 251_500 501_750 751_1000 {
					cap append using "poolsim/pool-ests-a25-n_`n'-sd_`sd'-`c'-`r'.dta"
					if _rc di "pool-ests-a25-n_`n'-sd_`sd'-`c'-`r'.dta"
				}
			}
		}
	}

	replace rerun=0 if rerun==.
	*tostring rerun , gen(rerunp)
	*replace rerunp=rerunp+string(npool)

	keep sdcond ncond cutcond rep grade year id nrtg npool convphet ///
		mstar* sstar* icchat* rerun

	reshape wide ///
		mstar_phet mstar_phet_se sstar_phet sstar_phet_se icchat_phet ///
		icchat_phet_se convphet rerun ///
		, ///
		i(sdcond ncond cutcond rep grade year id) j(npool)

	forv i = 1/5 {
		local mods phet
		if `i' == 1 local mods phet
		foreach m in `mods' {
			rename mstar_`m'_se`i' mstar_`m'`i'_se
			rename sstar_`m'_se`i' sstar_`m'`i'_se
			rename icchat_`m'_se`i' icchat_`m'`i'_se
		}
	}

	tempfile npool
	save `npool' , replace

	drop _all


* ---------------------------------------------------------------------------- *
* Load regular sim results


	foreach n in 10 25 50 100 200 {
		foreach sd in equal trend {
			foreach c in skewed wide mid mixed {
				foreach r in 1_250 251_500 501_750 751_1000 {
					cap append using "fullsim/ests-a25-n_`n'-sd_`sd'-`c'-`r'.dta"
					if _rc di "ests-a25-n_`n'-sd_`sd'-`c'-`r'.dta"
				}
			}
		}
	}


* ---------------------------------------------------------------------------- *
* Merge results


	merge 1:1 sdcond ncond cutcond rep grade year id using "`npool'"


* ---------------------------------------------------------------------------- *
* Confirm convergence


	local mods fhet fhom phet trend phet1 phet2 phet3 phet4 phet5
	foreach m in `mods' {
		replace conv`m' = 0 if conv`m' == .
	}

	egen tagrep = tag(ncond cutcond sdcond rep CC)
	egen tagcell = tag(ncond cutcond sdcond rep CC grade year)
	bys ncond cutcond sdcond CC grade year : egen nnum = mean(nrtg)
	replace nnum = round(nnum)

	forv i = 1/5 {
		replace convphet`i' = . if sdcond == "trend"
		replace convphet`i' = . if grade >= `i'
	}
	
	sum conv*
	
	table ncond sdcond cutcond , c(m convfhet m convphet m convtrend m convfhom)
	table ncond sdcond cutcond if tagcell , ///
		c(sum convfhet sum convphet sum convtrend sum convfhom)
	table ncond sdcond cutcond if tagrep , ///
		c(sum convphet1 sum convphet2 sum convphet3 sum convphet4 sum convphet5)


	* number of non-converged reps
	count if tagcell & convfhom == 0
	count if tagcell & convfhet == 0
	count if tagrep & convphet 	== 0
	count if tagrep & convtrend == 0
	count if tagrep & convphet1 == 0
	count if tagrep & convphet2 == 0
	count if tagrep & convphet3 == 0
	count if tagrep & convphet4 == 0
	count if tagrep & convphet5 == 0
	
	
* ---------------------------------------------------------------------------- *	
* Calculate true icc values and CV values


	g m2 = mstar^2
	g v2 = sstar^2
	bys ncond cutcond sdcond CC rep year grade (id) : egen N = sum(nrtg)
	g prtg = nrtg/N
	bys ncond cutcond sdcond CC rep year grade (id) : egen vb = sum(prtg*m2)
	bys ncond cutcond sdcond CC rep year grade (id) : egen vw = sum(prtg*v2)
	g test = vb+vw
	sum test , d
	g test2 = vb/(vb+vw)
	g icctrue = 1-vw

	drop m2 v2 N prtg vb vw test test2
	
	bys sdcond ncond cutcond rep grade : egen cv1=sd(sstar)
	bys sdcond ncond cutcond rep grade : egen cv2=mean(sstar)
	g cvtrue = cv1/cv2


* ---------------------------------------------------------------------------- *
* Calculate errors, squared errors, etc.


	local mods fhet fhom phet trend phet1 phet2 phet3 phet4 phet5 
	
	qui foreach m in `mods' {
		
		foreach v in m s {
			
			g `v'err`m' = `v'star_`m' - `v'star
			g `v'sqerr`m' = `v'err`m'^2
			g `v'aberr`m' = abs(`v'err`m')
			
			* relative (percent) error
			g `v'relerr`m' = (`v'star_`m' - `v'star)/`v'star
			
			g byte `v'inci`m' = 0
			replace `v'inci`m' = 1 if abs(`v'err`m') < (1.96*`v'star_`m'_se)
			replace `v'inci`m' = . if `v'star_`m' == . | `v'star_`m'_se == .
		
			* sampling variances

			bys ncond cutcond sdcond CC grade year id : egen `v'var`m'=sd(`v'star_`m')
			replace `v'var`m' = `v'var`m'^2
					
			g `v'se2_`m'=`v'star_`m'_se^2
			
		}

		g iccerr`m' = icchat_`m' - icctrue
		g iccsqerr`m' = iccerr`m'^2
		g byte iccinci`m' = 0
		replace iccinci`m' = 1 if abs(iccerr`m') < (1.96*icchat_`m'_se)
		replace iccinci`m' = . if icchat_`m' == . | icchat_`m'_se == .
			
	}

	foreach v in m s {
		forv i = 1/6 {
			if `i' < 6 g `v'ratio`i' = `v'varfhet/`v'varphet`i'
			if `i' == 6 g `v'ratio`i' = `v'varfhet/`v'varphet
		}
	}

	g double b0err = b0hat - b0
	g double b0sqerr = b0err^2
	g double b1err = b1hat - b1
	g double b1sqerr = b1err^2

	g mtrue = mstar
	g strue = sstar

	compress


* ---------------------------------------------------------------------------- *
* Final convergence rates, sample sizes, etc. for paper text


	table sdcond ncond , c(min nrtg mean nrtg max nrtg)
	table ncond grade , by(sdcond) c(min icctrue m icctrue)
	table ncond grade , by(sdcond) c(min cvtrue m cvtrue)
	su icctrue
	su cvtrue

	preserve
	keep icctrue sdcond ncond grade
	tostring grade, replace
	fcollapse (mean) icctrue , by(sdcond ncond grade)
	list , clean noobs
	restore

	* is bias in HOMOP due to weighting?
	su serrfhom [aweight=nrtg]
	table ncond sdcond cutcond [aweight=nrtg] , c(m serrfhom)
	table ncond sdcond cutcond , c(m serrfhom)


* ---------------------------------------------------------------------------- *
* save out summary stat files
	
	
	su conv*

	foreach stat in mean min max { // median min max 
		
		if "`stat'" == "median" {
		preserve
		tostring grade, replace
		fcollapse ///
			(mean) nrtg ///
			(`stat') mvar* svar* mratio* sratio* mse2* sse2* icctrue cvtrue ///
			(sum) tagrep ///
			, ///
			by(ncond cutcond sdcond grade)
		g stat = "`stat'"
		compress
		save "sim-summary-`stat'_$S_DATE.dta" , replace
		restore	
		}
		
		if "`stat'" != "median" {
		preserve
		tostring grade, replace
		fcollapse ///
			(`stat') xhet merr* serr* mrelerr* srelerr* ///
				iccerr* msqerr* ssqerr* iccsqerr* ///
				mvar* svar* mratio* sratio*  ///
				minci* sinci* iccinci* nrtg b?err b?sqerr mse2* sse2* ///
				conv*  icctrue cvtrue ///
			(sum) tagrep ///
			, ///
			by(ncond cutcond sdcond grade)
		g stat = "`stat'"
		compress
		save "sim-summary-`stat'_$S_DATE.dta" , replace
		restore
		}

	}
	

cap log close


