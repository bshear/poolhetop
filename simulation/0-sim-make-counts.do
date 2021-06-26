
version 14.2

clear all
set more off
set type double
macro drop _all

global numreps	"1000"
global ncond	"10 25 50 100 200"		// "10 25 50 100 200"
global sdcond	"trend equal"			// "trend equal"

do "0-sim-programs.do"

local sequal10		446692
local sequal25		167014
local sequal100		43158

local strend10		184976
local strend25		233041
local strend100		240021

local sequal200		873816
local strend200		382464

local sequal50		81772
local strend50		911412

qui foreach sd in $sdcond {
	foreach nc in $ncond {
	
		noi di _n "n=`nc', sd=`sd' (seed is `s`sd'`nc''):"
		
		set seed `s`sd'`nc''
		
		local outfilename "counts-n_`nc'-sd_`sd'.dta"
		
		preserve
		drop _all
		g junk=.
		cap save "`outfilename'"
		if _rc != 0 {
			noi di in red "counts already exist!"
			stop
		}
		restore
		
		if "`nc'" == "mixed" {
			local nmin = 10
			local nmax = 400
		}
		else {
			local nmin = `nc'
			local nmax = `nc'
		}
		
		if "`sd'" == "equal" {
			local trend ""
		}
		else {
			local trend "trend(0.1)"
		}
		
		drop _all
		
		qui gen_parms , grades(6) years(1) nmin(`nmin') nmax(`nmax') `trend'

		egen idrtg	= group(grade year id)
		egen rt		= group(grade year)
		g byte one	= 1
		sort grade year id

		compress

		tempfile parms
		save "`parms'" , replace

		qui forv r = 1/$numreps {
			
			noi di _c "..`r'"
			
			use "`parms'" , clear
			g int rep = `r'
			
			qui gen_data
			
			append using "`outfilename'"
			saveold "`outfilename'" , replace
			
		}
	}
}

use "`outfilename'" , clear
compress
drop junk
saveold "`outfilename'" , replace

exit

