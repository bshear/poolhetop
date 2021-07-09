*11jun2020: simplify file for distribution
*22dec2019
*ben shear

/*
Process raw downloaded Excel file and subset to analytic sample used
in the example. The data can be downloaded at:

https://www.cde.state.co.us/assessment/coassess-dataandresults#summarydata
*/

* ---------------------------------------------------------------------------- *
* setup


local date: display %td_CCYY-NN-DD date(c(current_date), "DMY")
di "`date'"

clear all
version 14.2
set type double


* ---------------------------------------------------------------------------- *
* load data


import excel "example/Math-TCAP-District_School_2014.xlsx" , firstrow clear


* ---------------------------------------------------------------------------- *
* process


keep SchoolNo Grade R S U W Y
ren R n 
ren S level1
ren U level2
ren W level3
ren Y level4

replace n = "." if n == "<16"
destring n , replace

destring Grade , replace
drop if Grade > 8

forv i = 1/4 {
	replace level`i' = "." if level`i' == "--"
	destring level`i' , replace
}

drop if SchoolNo == "0000"	// these are state/district totals

ren SchoolNo grpid
ren Grade grade

egen nobs = rowtotal(level1-level4)
drop if nobs == 0

g grdc = grade-3

* potential problematic 0 cells
egen num0 = anycount(level1-level4) , val(0)
table num0

* find schools with data for all grades
bys grpid : g nn = _N

egen tmptag = tag(grpid)
table nn if tmptag
drop tmptag

keep if nn == 6

table num0

drop if inlist(grpid, "1376", "8359")
// these were not matched in other CO datasets we have constructed

su nobs , d

g one = 1
egen numid = group(grpid)
egen numid_grd = group(numid grdc)

drop nn

* ---------------------------------------------------------------------------- *
* save


compress

save "example/example_data.dta" , replace

