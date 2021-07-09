*! hetop2 v1.1 25jul2018
* uses hetop2_lfw.ado and hetop2_p.ado

* changed start vals to use iweights instead of fweights to allow non-integer
* frequency counts

program define hetop2 , eclass sortpreserve
	
	version 13.1
	
	// setup, e.g., macro drop
	macro drop HET_*
	
	if replay() {
		if "`e(cmd)'" != "hetop2" error 301
		Replay `0'
	}
	else {
		Estimate `0'
	}
	
	// cleanup, e.g. macro or constraint drop
	constraint drop $HET_usecns
	macro drop HET_*
	
end


program Estimate , eclass sortpreserve
	
	syntax 	namelist , 	/// stubname of frequency count variables;
						/// must be ordered 1, 2, ... K+1
		K(integer)	 	/// number of cutscores (same for all years)
		GRPMN(varlist)	/// get an estimated mean for each of these unqiue values
		GRPSD(varlist)	/// get an estimated lnSD for each of these values
		DATAMN(varlist)	/// wgt'd mean estimates within these groups sum to 0
		DATASD(varlist)	/// wgt'd lnsigma estimates within these groups sum to 0
		DATACT(varlist)	/// thresholds estimated uniquely for each of these groups
		[ 					///
			BSD(varlist)	/// linear trend for lnsigma estimates over this continuous var
			GENerate		/// does nothing right now
			USECONSTraints	/// not used now, could be to pass add'l constraints
			ALLCUTS(string)			///
			STARTvals		/// if specified attempt to generate start values
							/// else, just use 0's plus cutscores
			* 				/// ml options
		]

	preserve

	tempfile orig
	qui save "`orig'" , replace
	
	* ------------------------------------ *
	* setup and data checking	
	
	local freq `1'
	mlopts mlopts , `options'

	unab fnames : `freq'*
	local nfnames : word count `fnames'
	if `nfnames' != `=`k'+1' {
		di as error "with `nfnames' frequency count variables need " ///
					"`=`nfnames'-1' cutscores"
		exit 498
	}
	forv i = 1/`nfnames' {
		cap confirm v `freq'`i'
		if _rc != 0 {
			di as error "frequency counts must be named " ///
						"`freq'1,..,`freq'(cuts+1)"
			exit 498
		}
	}

	* ------------------------------------ *	
	* define necessary constraints
	
	* if use has supplied cut score constraints, just use those
	if "`allcuts'" != "" {

		global HET_usecns `allcuts'

	}

	if "`allcuts'" == "" {
	
	global HET_usecns

	foreach p in mn sd {

		tempvar pg`p' ng1 ng N
		egen `ng1' = rowtotal(`freq'*)
		bys `data`p'' : egen `N' = sum(`ng1')
		bys `data`p'' `grp`p'' : egen `ng' = sum(`ng1')
		g double `pg`p'' = `ng'/`N'

		if "`p'" == "mn" local parm "means"
		if "`p'" == "sd" local parm "lnsigma"
		
		qui levelsof `data`p'' , local(datasets)
		qui foreach d in `datasets' {

			constraint free
			local j = r(free)
			
			levelsof `grp`p'' if `data`p'' == `d' , local(groups)
			local words : word count `groups'

			local cns ""
			forv i = 1/`words' {
				local t `: word `i' of `groups''
				qui sum `pg`p'' if `grp`p'' == `t' & `data`p'' == `d'
				if `i' == 1 local cns "`=r(mean)'*[`parm']:`t'.`grp`p''"
				if `i' >  1 local cns "`cns' + `=r(mean)'*[`parm']:`t'.`grp`p''"
			}
			local cns "`cns' = 0"
			constraint define `j' `cns'
			global HET_usecns "$HET_usecns `j'"
			
		}
	
	}
	
	// add'l constraints if we have a linear trend for lnsigma
	
	if "`bsd'" != "" {
		
		qui levelsof `datasd' , local(datasets)
		qui foreach d in `datasets' {
			
			constraint free
			local j = r(free)
			
			levelsof `grpsd' if `datasd' == `d' , local(groups)
			local words : word count `groups'

			local cns ""
			forv i = 1/`words' {
				local t `: word `i' of `groups''
				qui sum `pgsd' if `grpsd' == `t' & `datasd' == `d'
				if `i' == 1 {
					local cns "`=r(mean)'*[lnsigma]:`t'.`grpsd'#c.`bsd'"
				}
				if `i' >  1 {
					local cns "`cns' + `r(mean)'*[lnsigma]:`t'.`grpsd'#c.`bsd'"
				}
			}
			local cns "`cns' = 0"
			constraint define `j' `cns'
			global HET_usecns "$HET_usecns `j'"
			
		}
		
	}
	}	// define constraints

	* ------------------------------------ *
	* these are frequency weights, essentially
	forv i = 1/`=`k'+1' {
		global HET_f`i' "`freq'`i'"
	}
	
	* ------------------------------------ *
	* cutscore equations
	local cuteq
	forv i = 1/`k' {
		local cuteq "`cuteq' (cut`i': ibn.`datact' , nocons)"
	}
	
	if "`bsd'" != "" {
		local lnsigmaslope "ibn.`grpsd'#c.`bsd'"
	}
	
	local meanseqn		"ibn.`grpmn' , nocons"
	local lnsigmaeqn	"ibn.`grpsd' `lnsigmaslope' , nocons"
	
	* ------------------------------------ *
	* TODO: could add code to find better starting values
	* not sure that will help speed, but it might
	
	
	* ------------------------------------ *
	* start values
	
		// if option `startvals' is specified, this will attempt to get
		// very rough starting values
		// otherwise, returns 0s as initial values for all parameters but
		// uses initial values for cutscores
	
	start_vals , mid(`grpmn') md(`datamn') sid(`grpsd') sd(`datasd')		///
			ct(`datact') k(`k') fnames(`fnames') fstub(`freq') bsd(`bsd')	///
			`startvals'
	
	tempname b0
	matrix `b0' = r(b0)

	* ------------------------------------ *
	* call -ml-
	ml model lf hetop2_lfw	(means:   `meanseqn') 	///
							(lnsigma: `lnsigmaeqn')	///
							`cuteq' 				///
							,						///
							maximize 				///
							init(`b0')				///
							`mlopts'				///
							constraints($HET_usecns)
	
	* ------------------------------------ *
	* process estimation results and return
	
	// this puts the parts of e(b) and e(V) into their own matrices
	// not used for anything, but do return them
	
	// NOTE: if there are multiple rows of dataset with same group IDs,
	// these will NOT have multiple rows in the matrices
	// using -predict- to handle that for now.
	
	qui {
	
	tempname b V M MV G GV G0 G1 GbV S SV X B 
	mat `b' = e(b)
	mat `b' = `b''
	mat `V' = e(V)
	
	qui levelsof `grpmn' , local(m)
	local nm : word count `m'
	mat `M'  = J(`nm',1,.)
	mat `MV' = `V'[1..`nm',1..`nm']
	local mnames
	qui forv i = 1/`nm' {
		local g `: word `i' of `m''
		local mnames "`mnames' mprime`g'"
		local j = rownumb(`b', "means:`g'.`grpmn'")
		mat `M'[`i',1] = `b'[`j',1]
	}
	mat rownames `M' = `mnames'
	
	qui levelsof `grpsd' , local(s)
	local ns : word count `s'
	mat `G' = J(`ns',1,.)
	mat `S' = J(`ns',1,.)
	local gnames
	local snames
	qui forv i = 1/`ns' {
		local g `: word `i' of `s''
		local gnames "`gnames' gprime`g'"
		local snames "`snames' sprime`g'"
		local j = rownumb(`b', "lnsigma:`g'.`grpsd'")
		mat `G'[`i',1] = `b'[`j',1]
		mat `S'[`i',1] = exp(`G'[`i',1])
	}

	mat rownames `G' = `gnames'
	mat rownames `S' = `snames'
	
	matrix `GV' = `V'[`=`nm'+1'..`=`nm'+`ns'',`=`nm'+1'..`=`nm'+`ns'']
	matrix `SV' = diag(`S') * `GV' * diag(`S')
	
	if "`bsd'" != "" {
	
		mat `G1' = J(`ns',1,.)		
		qui forv i = 1/`ns' {
			local g `: word `i' of `s''
			local j = rownumb(`b', "lnsigma:`g'.`grpsd'#c.`bsd'")
			mat `G1'[`i',1] = `b'[`j',1]
		}
		matrix rownames `G1' = `gnames'	// b1 estimates
		
		matrix `G0'       = `G'			// b0 estimates
		
		// 
		mat `GbV' = `V'[`=`nm'+1'..`=`nm'+`ns'+`ns'', ///
						`=`nm'+1'..`=`nm'+`ns'+`ns'']
		
		local model ibn.`grpsd' `lnsigmaslope'
		fvrevar `model', stub(_x_)
		local colnames = ""
		local i = 1
		foreach wrd in `r(varlist)' {
			if strpos("`wrd'", "b.") > 0 {
				drop _x_`i'
			}
			else {
				local colnames `colnames' `wrd'
			}
			local i = `i' + 1
		}

		mkmat _x_*, mat(`X')
		mat colnames `X' = `colnames'
		drop _x_*
		
		matrix `B'  = `b'[`=`nm'+1'..`=`nm'+`ns'+`ns'',1]
		matrix `G'  = `X'*`B'
		matrix `GV' = `X'*`GbV'*`X''
		
		local r = rowsof(`G')
		mat `S' = J(`r',1,.)
		forv i = 1/`r' {
			matrix `S'[`i',1] = exp(`G'[`i',1])
		}
		
		matrix `SV' = diag(`S') * `GV' * diag(`S')
		
	}
	
	tempname gprime_se mprime_se sprime_se
	if rowsof(`GV') > 1 {
		// NOTE: this could be improved.
		// now that standardizing is done with -predict- not urgent.
		matrix `gprime_se' = J(rowsof(`GV'), 1,0)
		matrix `sprime_se' = J(rowsof(`GV'), 1,0)
	}
	else {
		matrix `gprime_se' = J(1,1,0)
		matrix `sprime_se' = J(1,1,0)
		local cn : rown `GV'
		noi di "`cn'"
		mat rownames `gprime_se' = "lnsigma:`cn'"
		mat rownames `sprime_se' = "sd:`cn'"
	}
	matrix `mprime_se' = vecdiag(cholesky(diag(vecdiag(`MV'))))'


	* ------------------------------------ *
	* return things
	
	if "`bsd'" != "" {
		ereturn matrix gprime_b_vcov `GbV'
		ereturn matrix gprime_b1 `G1'
		ereturn matrix gprime_b0 `G0'
		ereturn matrix X `X'
		ereturn matrix B `B'
	}
	
	ereturn matrix gprime_se	`gprime_se'
	ereturn matrix gprime_vcov	`GV'
	ereturn matrix gprime		`G'
	
	ereturn matrix sprime_se	`sprime_se'
	ereturn matrix sprime_vcov	`SV'
	ereturn matrix sprime		`S'
	
	ereturn matrix mprime_se	`mprime_se'
	ereturn matrix mprime_vcov	`MV'
	ereturn matrix mprime		`M'
	
	ereturn matrix b0			`b0'
	
	ereturn local freq_stub		"`freq'"
	ereturn local mean_datasets "`datamn'"
	ereturn local sd_datasets	"`datasd'"
	ereturn local mean_groups	"`grpmn'"
	ereturn local sd_groups		"`grpsd'"
	ereturn local cut_datasets	"`datact'"
	
	ereturn local sd_intercept	"`grpsd'"
	ereturn local sd_slope		"`bsd'"		// this can be missing
	
	ereturn local meanseqn		"`meanseqn'"
	ereturn local lnsigmaeqn 	"`lnsigmaeqn'"
	ereturn local predict		"hetop2_p"
	ereturn local cmd			"hetop2"
	
	}
	
	Replay
	
	restore
	
end

program parm_names , sclass
	
	syntax , [ mn(string) sd(string) ct(string) k(integer -1) bsd(string) ]
	
	qui levelsof `mn' , local(mns)
	qui levelsof `sd' , local(sds)
	qui levelsof `ct' , local(cts)
	
	local pnames
	
	foreach n in `mns' {
		local pnames "`pnames' means:`n'.`mn'"
	}
	
	foreach n in `sds' {
		local pnames "`pnames' lnsigma:`n'.`sd'"
	}
	
	if "`bsd'" != "" {
		foreach n in `sds' {
			local pnames "`pnames' lnsigma:`n'.`sd'#c.`bsd'"
		}
	}
	
	foreach n in `cts' {
		forv i = 1/`k' {
			local pnames "`pnames' cut`i':`n'.`ct'"
		}
	}
	
	sreturn local names `pnames'
	
end

program start_vals , rclass
	syntax [ , ///
		mid(varlist max=1)		///
		md(varlist max=1)		///
		sid(varlist max=1)	 	///
		sd(varlist max=1)		///
		ct(varlist max=1) 		///
		k(integer -1)			///
		fnames(varlist) 		///
		fstub(string)			///
		bsd(string) 			///
		STARTvals				///
		]

	tempfile tempstart
	qui save "`tempstart'" , replace
	
	qui levelsof `mid' , local(mn_levs)
	local nlevmn : word count `mn_levs'
	qui levelsof `sid' , local(sd_levs)
	local nlevsd : word count `sd_levs'
	qui levelsof `ct' , local(ct_levs)
	local nlevct : word count `ct_levs'
	
	local totalparms = `nlevmn' + `nlevsd'
	if "`bsd'" != "" local totalparms = `totalparms' + `nlevsd'
	
	if "`startvals'" != "" {
	
		// do custom start vals if specified; else all 0s and cutscores
	
		tempvar rowid y rowm rows rowf rowln dev2 bigm bigln pname
		tempname b0m b0lns b0
		
		quietly {
		
		g `rowid' = _n
		keep `fnames' `rowid' `mid' `md' `sid' `sd'
		reshape long `fstub' , i(`rowid' `mid' `md' `sid' `sd') j(`y')
		sort `rowid' `y'
		by `rowid' : g `rowm' = sum(`y'*`fstub')/sum(`fstub')
		by `rowid' : replace `rowm' = `rowm'[_N]
		g `dev2' = (`y'-`rowm')^2

		by `rowid' : g `rows' = sum(`dev2'*`fstub')/sum(`fstub')
		by `rowid' : replace `rows' = `rows'[_N]^.5
		by `rowid' : egen `rowf' = sum(`fstub')
		collapse `rowm' `rows' `rowf' , by(`rowid' `mid' `md' `sid' `sd')

		sort `md'
		by `md' : g `bigm' = sum(`rowm'*`rowf')/sum(`rowf')
		by `md' : replace `bigm' = `bigm'[_N]
		
		replace `rows' = 1 if `rows' == 0 | `rows' == .
		g `rowln' = log(`rows')
		
		sort `sd'
		by `sd' : g `bigln' = sum(`rowln'*`rowf')/sum(`rowf')
		by `sd' : replace `bigln' = `bigln'[_N]

		replace `rowln' = `rowln'-`bigln'
		replace `rowm' = (`rowm' - `bigm')/exp(`bigln')

		preserve
		collapse `rowm' [iweight=`rowf'], by(`mid')
		tostring `mid' , gen(`pname')
		replace `pname' = "means:"+`pname'+".`mid'"
		sort `pname'
		mkmat `rowm' , matrix(`b0m')
		qui levelsof `pname'
		matrix rownames `b0m' = `r(levels)'
		mat lis `b0m'
		restore

		preserve
		collapse `rowln' [iweight=`rowf'], by(`sid')
		tostring `sid' , gen(`pname')
		replace `pname' = "lnsigma:"+`pname'+".`sid'"
		sort `pname'
		mkmat `rowln' , matrix(`b0lns')
		qui levelsof `pname'
		matrix rownames `b0lns' = `r(levels)'
		mat lis `b0lns'
		restore

		// create matrix for means and lnsigmas
		matrix `b0' = `b0m'' , `b0lns''		
		
		if "`bsd'" != "" {
				
			// TO DO!!
			tempname b0slope
			matrix `b0slope' = J(1,`nlevsd', 0)
			local b0slopenames
			foreach s in `sd_levs' {
				local b0slopenames "`b0slopenames' lnsigma:`s'.`sid'#c.`bsd'"
			}
			
			matrix colnames `b0slope' = `b0slopenames'
			
			matrix `b0' = `b0' , `b0slope'
			
		}
		
		}
		
		qui use "`tempstart'" , clear	
	
	}
	else {
		
		// do this if "`startvals'" == ""
		
		tempname b0
		matrix `b0' = J(1,`totalparms',0)
		
		local b0names
		
		foreach m in `mn_levs' {
			local b0names "`b0names' means:`m'.`mid'"
		}
		
		foreach s in `sd_levs' {
			local b0names "`b0names' lnsigma:`s'.`sid'"
		}
		
		if "`bsd'" != "" {
			foreach s in `sd_levs' {
				local b0names "`b0names' lnsigma:`s'.`sid'#c.`bsd'"
			}
		}
		
		matrix colnames `b0' = `b0names'
	
	}
	
	// now add on cutscores
	
	foreach d in `ct_levs' {
		
		preserve
		qui keep if `ct' == `d'		
		tempname bc
		local bcnames
		mata: get_prop("`fnames'", "`bc'", `k')
		forv i = 1/`k' {
			mat `bc'[1,`i'] = invnormal(`bc'[1,`i'])
			local bcnames "`bcnames' cut`i':`d'.`ct'"
		}
		matrix colnames `bc' = `bcnames'
		mat `b0' = `b0', `bc'		
		restore
		
	}

	return matrix b0 `b0'
	
	
end


mata:
void get_prop(                  ///
          string scalar vnames,	///
          string scalar pmname,
		  real	 scalar numcuts)
{
        real matrix d
		d = st_data(., vnames)
		d = colsum(d)
		d = d:/rowsum(d)
		d = runningsum(d)
		d = d[1,1..numcuts]
		st_matrix(pmname, d)
}
end


program Replay

	ml display

end


exit

