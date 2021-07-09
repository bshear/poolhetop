*! hetop2_lfw v1.0 june 9 2016
* hetop2_lfw does the LL calculation in wide format

program hetop2_lfw
	version 13.1
	
	gettoken lnf rest: 0

	gettoken mu rest : rest
	gettoken lnsigma rest : rest
	local sigma exp(`lnsigma')
	
	local j = 0
	foreach kappa in `rest' {
		local j = `j' + 1
		local kappa`j' `kappa'
	}
	local K = `j'
	local M = `j' + 1
	
	forv i = 1/`M' {
		tempvar p`i'
	}
	qui gen double `p1' = normal((`kappa1'-`mu')/`sigma')
	forv i = 2/`K' {
		local j = `i'-1
		qui gen double `p`i'' = normal((`kappa`i''-`mu')/`sigma') - ///
							normal((`kappa`j''-`mu')/`sigma')
	}
	qui gen double `p`M'' = 1 - normal((`kappa`K''-`mu')/`sigma')
	
	local lleqn
	forv i = 1/`M' {
		if `i' == 1 local lleqn "${HET_f1} * ln(`p1')"
		if `i' >  1 local lleqn "`lleqn' + ${HET_f`i'} * ln(`p`i'')"
	}
		
	qui replace `lnf' = `lleqn'
	
end

