

cap program drop gen_data
program define gen_data , sclass sortpreserve
	version 13.1
	
	syntax [ , ///
		numcuts(integer 3) ///
		freqname(string) ///
	]
	
	// requires id grade year nrtg mstar sstar
	
	tempname rowid
	g `rowid' = _n	
	
	tempfile orig
	save "`orig'" , replace
	
	qui sum grade
	local R = r(max)
	qui sum year
	local T = r(max)
	qui sum id
	local G = r(max)
		
	if "`freqname'" == "" local freqname "freq"
		
	expand nrtg
	
	g double ystar = rnormal(mstar, sstar)
	
	local cuts "05 20 30 50 55 80 95"
	g tmpc = .
	foreach ct in `cuts' {
		
		local curp_num	= `ct' / 100
		replace tmpc	= 0
		local cur_cut	= invnormal(`curp_num')
		replace tmpc	= 1 if ystar < `cur_cut'
		bysort `rowid' : egen nb_`ct' = sum(tmpc)
		
	}
	
	egen tag = tag(`rowid')
	keep if tag
	
	keep `rowid' nb_*
	
	compress
	
	merge 1:1 `rowid' using "`orig'" , nogen
	
end

cap program drop gen_parms
program define gen_parms , sclass
	version 13.1

	syntax [ , ///
		GROUPS(integer 25) ///
		GRADES(integer 1) ///
		YEARS(integer 1) ///
		NMIN(integer 100) ///
		NMAX(integer 100) ///
		TREND(string) ///
	]

	local R = `grades'
	local T = `years'
	local G = `groups'
	
	/*
	set seed 1234
	local R = 3
	local T = 3
	local nmin = 100
	local nmax = 100
	local trend "trend"
	*/
	
	set obs 5
	
	mat b0 = (0.75 \ 0.85 \ 0.95 \ 1.05 \ 1.15)
	svmat b0
	rename b01 b0
	replace b0 = ln(b0)
	
	if "`trend'" != "" {
		mat b1 = (-`trend' \ -`=`trend'/2' \ 0.0 \ `=`trend'/2' \ `trend')	
		svmat b1
		rename b11 b1
		fillin b0 b1
		drop _fillin
	}
	else {
		mat b1 = (0 \ 0 \ 0 \ 0 \ 0)
		svmat b1
		rename b11 b1
		expand 5
	}
	
	sort b0 b1
	g id = _n
	
	g nrtg = round((runiform()*`=`nmax'-`nmin'') + `nmin') // 10-210
	
	qui sum b0 [fweight = nrtg]
	replace b0 = b0 - r(mean)
	qui sum b1 [fweight = nrtg]
	replace b1 = b1 - r(mean)
	
	expand `R'
	bys id : g grade = _n - 1
*	bys id : g grade = _n - (`R'/2)
	
	expand `T'
	bys id grade : g year = _n
	
	g double gprime = b0 + b1*grade
	g double sprime = exp(gprime)
	
	g random = runiform()
	g mprime = ///
		cond(random < 0.20, -0.6, ///
		cond(random < 0.40, -0.3, ///
		cond(random < 0.60, 0.0, ///
		cond(random < 0.80, 0.3, 0.6 ))))
	
	bys id : egen ng = sum(nrtg)
	egen N = sum(nrtg)
	bys grade year : egen Nrt = sum(nrtg)

	g pg	= ng / N		// proportion of total sample in this group
	g prtg	= nrtg / Nrt	// prop of year/grade sample in this group

	g vprime = sprime^2
	
	// create true star metric parameters
	bys grade year : egen bigm = sum(prtg*mprime)
	replace mprime = mprime - bigm
	drop bigm
	bys grade year : egen bigm = sum(prtg*mprime)	
	g mdevsq = (mprime - bigm)^2
	by grade year : egen vb = sum(prtg*mdevsq)
	by grade year : egen vw = sum(prtg*vprime)
	g sigma = (vb + vw)^.5
	g mstar = (mprime - bigm) / sigma
	g sstar = (vprime^.5) / sigma

	g icc = vb/(vb+vw)

	bys grade year : egen cv1 = sd(sstar)
	bys grade year : egen cv2 = mean(sstar)
	g cv = cv1/cv2
	drop cv1 cv2
	
	noi di "cv:"
	table grade year , c(m cv)
	qui su cv
	di round(`r(mean)',0.001)
	noi di "icc:"
	table grade year , c(m icc)
	qui su icc
	di round(`r(mean)',0.001)
	
	keep grade year id b0 b1 nrtg mstar sstar mprime
	
	compress
	
end	

/*
*test	
		set seed 1234
		drop _all
		qui gen_parms , grades(6) years(1) nmin(50) nmax(50) trend(0.05)
		egen idrtg	= group(grade year id)
		egen rt		= group(grade year)
		g byte one	= 1
		g lnstar = log(sstar)
		tw (scatter lnstar grade) (lfit lnstar grade) , by(id)
		sort grade year id
		su sstar
*/		
