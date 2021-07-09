*! version 1.2 august 12, 2016
* hetop2_p predict means and standard deviations in prime and star metric

* v1.1: introduced use of -matselrc- to create prime matrices for star 
*		calculations to speed up code
* v1.2: no longer fill Gprime with0s for HOMOP; allow to be filled
*		automatically along with others

program define hetop2_p , rclass sortpreserve
	version 13.1
	
	syntax [ namelist(max=1) ] , /// suffix for predicted values
	[ ///
		PRIME	///
		STAR	///
		SE		///
		STATus	///		print status updates?
	]
	
	* ------------------------------------------------------------------------ *
	* data prep and define locals
	
	if "`namelist'" != "" local suffix "_`namelist'"
	if "`namelist'" == "" local suffix ""
	
	local type		"double"
	
	if "`prime'" != "" {
		
		local mprime	"mprime`suffix'"
		local gprime	"gprime`suffix'"
		local sprime	"sprime`suffix'"
		
		foreach v in m s g {
			cap conf var ``v'prime'
			if _rc == 0 {
				di as error "error: variable ``v'prime' already exists."
				exit 498
			}
		}
		
		if "`se'" != "" {
			foreach v in m s g {
				cap conf var ``v'prime'_se
				if _rc == 0 {
					di as error "error: variable ``v'prime'_se already " ///
					"exists."
					exit 498
				}
			}
		}
	
	}

	if "`star'" != "" {
		
		local mstar	"mstar`suffix'"
		local sstar	"sstar`suffix'"
		local frequencies 	"`=e(freq_stub)'"
		local datasets		"`=e(mean_datasets)'"
		
		foreach v in `mstar' `sstar' {
			cap conf var `v'
			if _rc == 0 {
				di as error "error: variable `v' already exists."
				exit 498
			}
		}
		
		if "`se'" != "" {
			foreach v in `mstar'_se `sstar'_se {
				cap conf var `v'
				if _rc == 0 {
					di as error "error: variable `v' already exists."
					exit 498
				}
			}
		}
	
	}	

	* ------------------------------------------------------------------------ *
	* do prime metric predictions
	* NOTE: this is really double-prime, which might or might not be prime
	* it is the metric in which estimates are identified
	
	if "`prime'" != "" {
	
		if "`status'" != "" {
			noi di _n "computing prime estimates and standard errors..."
		}
		_predict `type' `mprime' , xb eq(#1)
		_predict `type' `gprime' , xb eq(#2)
		
		if "`se'" != "" {
			_predict `type' `mprime'_se , stdp eq(#1)
			_predict `type' `gprime'_se , stdp eq(#2)
		}
		
		if "`se'" != "" local sse " , se(`sprime'_se)"
		predictnl `type' `sprime' = exp(xb(#2)) `sse'
		label var `sprime' "exp(`gprime') Prediction"
		
	}
	
	
	* ------------------------------------------------------------------------ *
	* do star metric predictions
	* this is the important part
	
	if "`star'" != "" {
	
	if "`status'" != "" ///
		noi di _n "computing star estimates and standard errors..."
	
	local datamn	"`=e(mean_datasets)'"
	local grpmn		"`=e(mean_groups)'"
	local datasd	"`=e(sd_datasets)'"
	local grpsd		"`=e(sd_groups)'"
	
	tempname Bfull Vfull
	matrix `Bfull' = e(b)
	matrix `Vfull' = e(V)
	
	qui levelsof `datasets' , local(dsets)	// standardize within each of these
	
	tempvar uniqueid
	g `uniqueid' = _n	// used for sorting, merging, etc.

	foreach d in `dsets' {

		// loop through and for each of the datasets in which means are
		// standardized, obtain the matrices needed for star calculations
		
		if "`status'" != "" noi di _n "    `datasets' = `d'..."
		
		preserve
		
		qui keep if `datasets' == `d'
		
		// this sort order is the one that will show up in the matrices
		*tempvar cursort
		*sort `grpmn' 
		*g `cursort' = _n
		sort `uniqueid'		// this shouldn't actually change the order
		
		local G = _N		// number of groups
		
		* -------------------------------------------------------------------- *
		
		if "`status'" != "" noi di "	generate vars..."		
		
		* create these variables for calculations
		tempvar ng pg ninv1 ninv
		
		* create these matrices
		tempname Mprime Gprime Sprime Vprime Lambdaprime ///
				Mprime_se Sprime_se Gprime_se ///
				Omegaprime PI Zprime Wprime P One P2 ninv ninv1 n N

		egen `ng'			= rowtotal(`frequencies'*)
		
		qui sum `ng'
		scalar `N'			= r(sum)
		
		g double `pg'		= `ng'/scalar(`N')
		g double `ninv1'	= 1/(`ng'-1)
		g double `ninv'		= 1/`ng'
		
		mkmat `pg'			, matrix(`P')
		mkmat `ninv1' 		, matrix(`ninv1')
		mkmat `ninv'  		, matrix(`ninv')
		mkmat `ng' 			, matrix(`n')
		
		matrix `n' 			= `n''
		matrix `P'			= `P''
		matrix `One'		= J(`G',1,1)'
		matrix `PI'			= I(`G') - `One''*`P'
		matrix `P2'			= hadamard(`P',`P')

		drop `ninv' `ninv1'

		* -------------------------------------------------------------------- *
		* create matrices
	
		if "`status'" != "" noi di "	generate initial matrices..."		
		
		foreach v in Mprime Gprime Sprime {
			matrix ``v'' = J(`G',1,.)
		}
		foreach v in Vprime Omegaprime Lambdaprime {
			matrix ``v'' = J(`G',`G',.)
		}
		
		local rmnames
		local rlnames
		
		if "`status'" != "" noi di "	step 1..."
		
		// maybe do this with creating a string variable?
		forv r = 1/`G' {
			qui sum `grpmn' if _n == `r'
			local rmnames "`rmnames' means:`=r(mean)'.`grpmn'"
			qui sum `grpsd' if _n == `r'
			local rlnames "`rlnames' lnsigma:`=r(mean)'.`grpsd'"
		}

		
		if "`status'" != "" noi di "	step 2..."
		
		matselrc `Bfull' `Mprime' , r(1) c(`rmnames')
		matselrc `Bfull' `Gprime' , r(1) c(`rlnames')
		
		matrix `Mprime' = `Mprime''
		matrix `Gprime' = `Gprime''
		
		matselrc `Vfull' `Vprime'		, r(`rmnames') c(`rmnames')
		matselrc `Vfull' `Omegaprime'	, r(`rlnames') c(`rlnames')
		matselrc `Vfull' `Lambdaprime'	, r(`rmnames') c(`rlnames')
				
		if "`status'" != "" noi di "	step 3..."		
		
		qui levelsof `grpsd'		
		local grpsd_nlevels : word count `=r(levels)'
				
		
		* -------------------------------------------------------------------- *
		// if linear trend specified, need to modify the Gprime, Omega, 
		// and Lambda matrices
		
		if "`=e(sd_slope)'" != "." {
		
			if "`status'" != "" noi di "	sd slope calculations..."			
			
			/*
			X: design matrix for the current observations
			B: coefficients for int/slope of SDs
			U: Cov(B,B)
			D: Cov(Mprime,B)			
			*/
			
			tempname X B U D
			
			qui levelsof `grpsd' , local(grpsdvals)
			local parms : word count `grpsdvals'
			
			// names of parameters in the B vector
			forv r = 1/`G' {
				qui sum `grpsd' if _n == `r'
				local rlslopenames "`rlslopenames' lnsigma:`=r(mean)'.`grpsd'#c.`=e(sd_slope)'"
			}
			
			if "`status'" != "" noi di "	design matrix..."
			
			// create the design matrix X
			local model "ibn.`grpsd' ibn.`grpsd'#c.`=e(sd_slope)'"
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
			
			// create B matrix

			if "`status'" != "" noi di "	B matrix..."			
			
			matrix `B' = J(`=2*`parms'',1,.)			
			local parmnames
			forv r = 1/`parms' {
				local s `: word `r' of `grpsdvals''
				local parmnames "`parmnames' lnsigma:`s'.`grpsd'"
			}
			forv r = 1/`parms' {
				local s `: word `r' of `grpsdvals''
				local parmnames "`parmnames' lnsigma:`s'.`grpsd'#c.`=e(sd_slope)'"
			}
			
			local r = 1
			foreach p in `parmnames' {
				matrix `B'[`r',1] = ///
					`Bfull'[1,colnumb(`Bfull', ///
							"`p'")]							
				local r = `r'+1
			}
			
			matrix rownames `B' = `parmnames'
			
			// create D matrix

			if "`status'" != "" noi di "	D matrix..."
			
			matselrc `Vfull' `D' , r(`rmnames') c(`parmnames')
			
			// create U matrix
			
			if "`status'" != "" noi di "	U matrix..."
			
			matselrc `Vfull' `U' , r(`parmnames') c(`parmnames')
								
			// transformed matrices
			
			matrix `Gprime'			= `X'*`B'
			matrix `Omegaprime'		= `X'*`U'*`X''
			matrix `Lambdaprime'	= `D'*`X''
			
			*noi di "U"
			*noi mat lis `U'
			
		}
		
		matrix rowna `Vprime'		= `rmnames'
		matrix colna `Vprime'		= `rmnames'
		matrix rowna `Omegaprime'	= `rlnames'
		matrix colna `Omegaprime'	= `rlnames'
		matrix rowna `Lambdaprime'	= `rmnames'
		matrix colna `Lambdaprime'	= `rlnames'

		tempname a
		matrix `a' = `P' * `Gprime'
		local a = `a'[1,1]				// should convert these to scalars...
		local kappa = 1
		
		// ensure prime matrices are in metric for sum-to-zero constraints
		// i.e. convert from double to single prime
		
		if "`status'" != "" noi di "	convert double prime to prime..."		
		
		matrix `Vprime' = ///
			(exp(`a')^-2) * ( ///
				( `kappa'^2 * `PI' * `Vprime' * `PI'' ) - ///
				(`kappa'^(-1)) * ///
				(`PI' * `Mprime' * `P' * `Lambdaprime'' * `PI'' + ///
					`PI' * `Lambdaprime' * `P'' * `Mprime'' * `PI'' ) + ///
				(`kappa'^(-4)) * (`PI' * `Mprime' * `Mprime'' * `PI'' * ///
					( `P' * `Omegaprime' * `P'') ) ///
			)

		matrix `Lambdaprime' = ///
			(exp(`a')^-1) * ( ///
				( `kappa' * `PI' * `Lambdaprime' * `PI'') - ///
				(`kappa'^(-2)) * ///
					`PI' * `Mprime' * `P' * `Omegaprime' * `PI'' ///
			)
		
		matrix `Omegaprime' = ///
			`PI'*`Omegaprime'*`PI''
			
		matrix `Mprime' = (exp(`a')^-1) * ( `PI' * `Mprime' )
		matrix `Gprime' = `PI' * `Gprime'
		
		// compute Sprime, Wprime and Zprime
		
		forv r = 1/`G' {
			matrix `Sprime'[`r',1] = exp(`Gprime'[`r',1])
		}

		matrix `Wprime' = diag(`Sprime') * `Omegaprime' * diag(`Sprime')
		matrix `Zprime' = `Lambdaprime' * diag(`Sprime')		
		
		// prime metric standard errors
		
		matrix `Mprime_se' = vecdiag(cholesky(diag(vecdiag(`Vprime'))))'

		if `grpsd_nlevels' > 1 {
			matrix `Gprime_se' = vecdiag(cholesky(diag(vecdiag(`Omegaprime'))))'
			matrix `Sprime_se' = vecdiag(cholesky(diag(vecdiag(`Wprime'))))'
		}
		if `grpsd_nlevels' == 1 {
			matrix `Gprime_se' = J(`G',1,0)
			matrix `Sprime_se' = J(`G',1,0)
		}
		
		
		// create star matrices and standard errors
		
		if "`status'" != "" noi di "	transform to star metric..."		
		
		tempname sigmaw sigmab sigmaprime ///
			Q Mstar Sstar Mstar_se Sstar_se ///
			icchat R T varsigprime Vstar Wstar Zstar icchatvar ohbg
		
		// create ohbg ... "omega-hat-bar-g"

		if `grpsd_nlevels' == `G' {
			// 
			tempname ntilde
			matrix `ntilde' = invsym(((1/`G') * (`One' * `ninv1')))
			scalar `ohbg' = 1/(2*`ntilde'[1,1])
		}
		if `grpsd_nlevels' > 1 & `grpsd_nlevels' < `G' {
			tempvar nt nt1
			g `nt1' = `ng'-1
			
			bys `grpsd' : egen `nt' = sum(`nt1')
			sort `uniqueid'
			
			qui replace `nt' = 1/(2*`nt')
			qui sum `nt'
			scalar `ohbg' = (1/`G')*r(sum)
		}
		if `grpsd_nlevels' == 1 scalar `ohbg' = 1/(2*(scalar(`N')-`G'))
		
		// standardization
		
		matrix `Q' = ///
			(1/(1+(2*scalar(`ohbg')))) * ///
			hadamard( hadamard( `ninv'' , (`P' + `n' - `One') ) , `P' )

		matrix `sigmaw' = ///
			(1/(1+(2*scalar(`ohbg')))) * (`P' * hadamard( `Sprime' , `Sprime' ))
			
		matrix `sigmab' = ///
			(`P' * hadamard( `Mprime' , `Mprime' )) + ///
			( (1/(1+(2*`ohbg'))) * ( hadamard( `ninv'' , (`P2' - `P') ) * ///
				hadamard( `Sprime' , `Sprime' ) ) )
		
		matrix `sigmaprime' = ///
			cholesky( `P' * hadamard(`Mprime',`Mprime') + ///
			`Q' * hadamard(`Sprime', `Sprime') )

		matrix `varsigprime' = ///
			invsym(hadamard(`sigmaprime', `sigmaprime')) * ///
			(`P' * diag(`Mprime') * `Vprime' * diag(`Mprime') * `P'' + `Q' * ///
			diag(`Sprime') * `Wprime' * diag(`Sprime') * `Q'' + 2 * `P' * ///
			diag(`Mprime') * `Zprime' * diag(`Sprime') * `Q'')
			
		matrix `Mstar' = invsym(`sigmaprime')*`Mprime'
		matrix `Sstar' = invsym(`sigmaprime')*`Sprime'

		matrix `R' = ///
			`P' * diag(`Mstar') * `Vprime' + `Q' * diag(`Sstar') * `Zprime''
		
		matrix `T' = ///
			`P' * diag(`Mstar') * `Zprime' + `Q' * diag(`Sstar') * `Wprime'
				
		matrix `Vstar' = ///
			invsym(hadamard(`sigmaprime', `sigmaprime')) * ///
			( `Vprime' - (`Mstar' * `R' + `R'' * `Mstar'') + ///
				`Mstar' * `Mstar'' * `varsigprime' )
				
		matrix `Wstar' = ///
			invsym(hadamard(`sigmaprime', `sigmaprime')) * ///
			(`Wprime' - (`Sstar' * `T' + `T'' * `Sstar'') + `Sstar' * ///
			`Sstar'' * `varsigprime')				
			
		matrix `Zstar' = ///
			invsym(hadamard(`sigmaprime', `sigmaprime')) * ///
			(`Zprime' - (`Mstar' * `T' + `R'' * `Sstar'') + ///
			`Mstar' * `Sstar'' * `varsigprime')

		matrix `Mstar_se' = vecdiag(cholesky(diag(vecdiag(`Vstar'))))'
		matrix `Sstar_se' = vecdiag(cholesky(diag(vecdiag(`Wstar'))))'
		
		
		// estimate icc and var(icc)
		
		if "`status'" != "" noi di "	estimate icc..."		
		
		matrix `icchat' = ///
			1 - ((1/(1+(2*scalar(`ohbg')))) * `P' * hadamard(`Sstar', `Sstar'))
				
		matrix `icchatvar' = ///
			4*(1/(1+(2*scalar(`ohbg'))))^2 * `P' * ///
			(diag(`Sstar')*`Wstar'*diag(`Sstar')) * `P''
		
		
		// do matrix row and column names again ...
		
		matrix rowna `Mstar' = `rmnames'
		matrix rowna `Mstar' = `d'.mstar:
		matrix rowna `Sstar' = `rlnames'
		matrix rowna `Sstar' = `d'.sstar:

		matrix rowna `Mstar_se' = `rmnames'
		matrix rowna `Mstar_se' = `d'.mstar_se:
		matrix rowna `Sstar_se' = `rlnames'
		matrix rowna `Sstar_se' = `d'.sstar_se:		
				
		// save the mstar, sstar, mstar_se and sstar_se as variables
		// need to enhance options and variable naming here
		
		if "`status'" != "" noi di "	cleanup..."

		* keep `uniqueid' `datasets' `datamn' `grpmn' `datasd' `grpsd'
		keep `uniqueid'
	
		svmat `type' `Mstar' , names(`mstar')
		rename `mstar'1 `mstar'
		svmat `type' `Sstar' , names(`sstar')
		rename `sstar'1 `sstar'
		
		if "`se'" != "" {
			svmat `type' `Mstar_se' , names(`mstar'_se)
			rename `mstar'_se1 `mstar'_se
			svmat `type' `Sstar_se' , names(`sstar'_se)
			rename `sstar'_se1 `sstar'_se
		}
		
		// return matrices in r(`name')
		// could add others
		
		if "`=e(sd_slope)'" != "." {
			return matrix X`d'	`X'
			return matrix B`d'	`B'
		}
		return matrix mstar`d' `Mstar'
		return matrix sstar`d' `Sstar'
		return matrix mstar`d'_se `Mstar_se'
		return matrix sstar`d'_se `Sstar_se'
		
		return local dset_levels "`dsets'"
		
		return scalar icchat`d' = `icchat'[1,1]
		return scalar icchat`d'_var = `icchatvar'[1,1]
		
		sort `uniqueid'
		qui compress

		tempfile temp`d'
		qui save "`temp`d''"

		restore
		
	}

	// combine estimates for each of the datasets
	// odd issues with sorting here...
	qui foreach d in `dsets' {
		merge 1:1 `uniqueid' using "`temp`d''" , nogen update
	}
	
	}
	
end

