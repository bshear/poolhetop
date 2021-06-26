
* ---------------------------------------------------------------------------- *
* base file:

di "$S_DATE"
di "$S_TIME"

version 14.2

which hetop2

local tempreps : subinstr global reps "_" " ", all 		
global repmin : word 1 of `tempreps'
global repmax : word 2 of `tempreps'	

local outfilename "pool-ests-a25-n_$ncond-sd_$sdcond-$cutname-${repmin}_${repmax}.dta"

use "counts-n_$ncond-sd_$sdcond.dta"	// data already generated

* ------------------------------- *
* programs to define cut constraints and grab correct counts

cap program drop makeC
program define makeC
	
	syntax anything	// but this should be "mid" "skewed" "wide" or "mixed"
	
	local dname "rt"
	
	if "`anything'" == "mixed" {
			
		local midcut1 = invnorm(0.2)
		local midcut2 = invnorm(0.5)
		local midcut3 = invnorm(0.8)

		local skewcut1 = invnorm(0.05)
		local skewcut2 = invnorm(0.3)
		local skewcut3 = invnorm(0.55)

		local widecut1 = invnorm(0.05)
		local widecut2 = invnorm(0.50)
		local widecut3 = invnorm(0.95)
		
		local usecns

		qui levelsof `dname' , local(rts)
		
		foreach r in `rts' {
			forv c = 1/3 {
				constraint free
				if inlist(`r',1,4) constraint `=r(free)' [cut`c']`r'.`dname' = `midcut`c''
				if inlist(`r',2,5) constraint `=r(free)' [cut`c']`r'.`dname' = `skewcut`c''
				if inlist(`r',3,6) constraint `=r(free)' [cut`c']`r'.`dname' = `widecut`c''
				local usecns `usecns' `=r(free)'				
			}
		}	
	}
	
	if "`anything'" != "mixed" {
		
		if "`anything'" == "mid" {
			local cut1 = invnorm(0.2)
			local cut2 = invnorm(0.5)
			local cut3 = invnorm(0.8)
		}
		
		if "`anything'" == "skewed" {
			local cut1 = invnorm(0.05)
			local cut2 = invnorm(0.3)
			local cut3 = invnorm(0.55)
		}
		
		if "`anything'" == "wide" {
			local cut1 = invnorm(0.05)
			local cut2 = invnorm(0.50)
			local cut3 = invnorm(0.95)		
		}
		
		local usecns
	
		qui levelsof `dname' , local(rts)
		
		foreach r in `rts' {
			forv c = 1/3 {
				constraint free
				constraint `=r(free)' [cut`c']`r'.`dname' = `cut`c''
				local usecns `usecns' `=r(free)'				
			}
		}	
	}

	global usecns `usecns'

end

cap program drop makeCounts
program define makeCounts

	syntax anything
	
	forv i = 1/4 {
		g freq`i' = .
	}

	if "`anything'" == "mid" {
		local c1 20
		local c2 50
		local c3 80
		replace freq1 = nb_`c1'
		replace freq2 = nb_`c2'-nb_`c1'
		replace freq3 = nb_`c3'-nb_`c2'
		replace freq4 = nrtg - nb_`c3'
	}
	else if "`anything'" == "skewed" {
		local c1 05
		local c2 30
		local c3 55
		replace freq1 = nb_`c1'
		replace freq2 = nb_`c2'-nb_`c1'
		replace freq3 = nb_`c3'-nb_`c2'
		replace freq4 = nrtg - nb_`c3'
	}
	else if "`anything'" == "wide" {
		local c1 05
		local c2 50
		local c3 95
		replace freq1 = nb_`c1'
		replace freq2 = nb_`c2'-nb_`c1'
		replace freq3 = nb_`c3'-nb_`c2'
		replace freq4 = nrtg - nb_`c3'
	}
	else if "`anything'" == "mixed" {
		
		// grade0 = mid
		// grade1 = skewed
		// grade2 = wide
		// grade3 = mid
		// grade4 = skewed
		// grade5 = wide
		
		local c01 20
		local c02 50
		local c03 80	
		
		local c11 05
		local c12 30
		local c13 55
		
		local c21 05
		local c22 50
		local c23 95
		
		local c31 20
		local c32 50
		local c33 80	
		
		local c41 05
		local c42 30
		local c43 55
		
		local c51 05
		local c52 50
		local c53 95
		
		qui forv g = 0/5 {
			replace freq1 = nb_`c`g'1'				if grade == `g'
			replace freq2 = nb_`c`g'2'-nb_`c`g'1'	if grade == `g'
			replace freq3 = nb_`c`g'3'-nb_`c`g'2'	if grade == `g'
			replace freq4 = nrtg - nb_`c`g'3'		if grade == `g'
		}

	}
	

end


* ------------------------------ *
* set up counts

makeCounts $cutname


* ------------------------------ *
* set up model calls

local fhetmod hetop2 freq , k(3) ///
			grpmn(idrtg) datamn(rt) ///
			grpsd(idrtg) datasd(rt) ///
			datact(rt)

local phetmod hetop2 freq , k(3) ///
			grpmn(idrtg) datamn(rt) ///
			grpsd(id)	 datasd(one) ///
			datact(rt)

local trendmod hetop2 freq , k(3) ///
			grpmn(idrtg) datamn(rt) ///
			grpsd(id)	 datasd(one) bsd(grade) ///
			datact(rt)
			
* use same frequencies for homop for now
			
local fhommod hetop2 freq , k(3) ///
			grpmn(idrtg) datamn(rt) ///
			grpsd(rt)    datasd(rt) ///
			datact(rt)


* ------------------------------ *
preserve
drop _all
g junk = 1
save "`outfilename'", replace
restore


* ------------------------------ *
keep rep grade year id rt idrtg one ///
	mstar sstar mprime b0 b1 nrtg freq*
compress

g sdcond	= "$sdcond"
g ncond		= "$ncond"
g cutcond	= "$cutname"

keep if rep >= $repmin
keep if rep <= $repmax


* ------------------------------ *
* smoothing
* freqfhom`i' will be homop frequency counts
* freq`i' will be counts used for hetop/trend/pooled
egen N = rowtotal(freq1-freq4)
egen num0 = anycount(freq1-freq4) , val(0)
forv i = 1/4 {
	g altfreq`i' = freq`i'
}
g byte xhet = 0
replace xhet = 1 if num0 > 2
replace xhet = 1 if freq1==0 & freq2==0
replace xhet = 1 if freq1==0 & freq4==0
replace xhet = 1 if freq3==0 & freq4==0
replace xhet = 1 if freq2==0 & freq3==0
g byte xhom = 0
replace xhom = 1 if num0==3 & (freq1!=0 | freq4!=0)
forv i = 1/4 {
	replace freq`i'		= (freq`i'+0.25)/(N+1)*N	if xhet==1
	replace altfreq`i' = (altfreq`i'+0.25)/(N+1)*N	if xhom==1
}
drop N num0
* ------------------------------ *



tempfile full
save "`full'"

qui forv r = $repmin / $repmax {
	
	noi di _c "`r'.."
	
	local cutcons "free"
	
	forv ng = 1/5 {
		
		use "`full'" , clear	
		
		keep if rep == `r'
		
		keep if grade<`ng'
		
		g npool = `ng'
		
		tempfile p`ng'
		
		qui levelsof rt , local(rt_levs)

		// only pooled models for now
		
		foreach m in phet {
			
			ereturn clear
			local try = 1
			
			local CC
			constraint drop _all			
			if "`cutcons'" == "fix" makeC $cutname
			if "`cutcons'" == "fix" local CC "allcuts($usecns)"
			
			cap ``m'mod' `CC' iter(200) start
		
			if !_rc & e(converged) == 1 {

				cap predict `m' , star se

				// estimated ICCs
				
				g double icchat_`m' = .
				g double icchat_`m'_se = .
				foreach lev in `rt_levs' {
					replace icchat_`m' = r(icchat`lev')  ///
						if rt == `lev' & id == 1
					replace icchat_`m'_se = r(icchat`lev'_var)^.5  ///
						if rt == `lev' & id == 1
				}			
				
				// save slope/intercept estimates for trend model
				
				if "`m'" == "trend" {
					qui tabulate `=e(sd_groups)' , gen(_x)
					mkmat _x* , mat(X)
					foreach b in 0 1 {
						mat b`b' = X * e(gprime_b`b')
						svmat double b`b' , name(b`b'hat)
						rename b`b'hat1 b`b'hat
					}
					drop _x*
				}
				
				if _rc g pred`m'rc = _rc

			}
				
			g int atmp`m' = `try'
			g byte conv`m' = e(converged)
			
		}
		
		save "`p`ng''" , replace
		
	}
	
	use "`p1'" , clear
	forv i = 2/5 {
		append using `p`i''
	}
	
	compress
	
	append using "`outfilename'"
		
	saveold "`outfilename'" , replace	

}


use "`outfilename'" , clear

describe
nmissing

sum
drop junk

compress

save "`outfilename'" , replace

di "$S_DATE"
di "$S_TIME"

exit
