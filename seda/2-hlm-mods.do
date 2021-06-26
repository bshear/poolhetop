*28mar2019
*02aug2018
*author: benjamin shear
*fit HLM models for shear & reardon pooled HETOP 
*uses SEDA data subset from step 1

version 14.2
clear all
set more off
macro drop _all
set type double
set matsize 11000

cd "Y:\GitHub\poolhetop\"

// load data
use "seda\SEDA_geodist_long_State_v21_sub.dta", clear

cd "seda\"

// HLM models

	* 0 - HOMOP equivalent: fixed effects for grades and years
	* 1 - fully pooled HETOP model: single constant parameter per district
	* 2 - linear trend pooled HETOP model: district-specific linear grade trends
	* 3 - linear trends pooled HETOP model: add district-specific year trends

// loop through states

levelsof stateabb , local(states)

qui foreach s in `states' {

	qui levelsof subject if stateabb == "`s'" , local(subjects)

	foreach v in `subjects' {

		noi di _n "`s' `v'"

		preserve

			keep if subject == "`v'" & stateabb == "`s'"
			
			g y = gprime
			g y_se = gprime_se	
			
			egen id = group(leaidC)
			
			keep grade year id leaidC y y_se
			
			// center grade and year within states
			qui su grade
			loc gm_grade=r(mean)
			g gradec=grade-`gm_grade'
			qui su year
			loc gm_year=r(mean)
			g yearc=year-`gm_year'
						
			sort id grade year
			
			g l2id		= _n
			g l3x		= 1
			g wgt		= y_se^2
			
			qui su wgt
			local mn_evar=r(mean)
			
			egen x = group(grade year)
			qui tabulate x , gen(gy)
			
			qui su x
			local gymax = `r(max)'
			
			keep  id l2id l3x y wgt gradec yearc gy1-gy`gymax'
			order id l2id l3x y wgt gradec yearc gy1-gy`gymax'
			
			// MODEL 0
			
			// homop equivalent model
			// just a fixed effect for each grade-year
			// want to estimate omega^2, which will be the anticipated variance
			// within grades and years and represent error variance from HOMOP model
			
			noi di _c "..m0"

			* create a 2-level MDM file
			hlm mkmdm using hlm`v'2 , id2(l2id) type(hlm2) ///
				l1(y wgt) ///
				l2(gradec yearc gy1-gy`gymax') ///
				run replace
				
			hlm mdmset hlm`v'2
			
			cap hlm hlm2 y int(gy1-gy`gymax' rand) rand, ///
				pwgt(wgt) cmd(hlm`v'`s'0) run replace robust 

				if _rc == 0 {
				
				mat tau2 = e(tau2)
				scalar s2 = tau2[1,1]
				estadd scalar s2 scalar(s2)
				estadd scalar mn_evar `mn_evar'
				
				est save hlm`v'est`s'0 , replace
				
				local dev = e(dev)
				local df = e(df)
				
				local devtest "dev(`dev' `df')"
				
				}
				else {
					
					local devest ""
					
				}	
				
			// make 3-level MDM file
		
			hlm mkmdm using hlm`v' , id2(l2id) id3(id) type(hlm3) ///
				l1(y wgt) ///
				l2(gradec yearc gy1-gy`gymax') ///
				l3(l3x) ///
				run replace
			
			hlm mdmset hlm`v'
						
			// MODEL 1
			// constant gamma model
			
			noi di _c "..m1"
			
			cap hlm hlm3 y int(int(rand) gy1-gy`gymax' rand) rand, ///
				pwgt(wgt) cmd(hlm`v'`s'1) run replace robust 
				
			if _rc == 0 {
				
				mat tau2 = e(tau2)
				mat tau3 = e(tau3)
				scalar s2 = tau2[1,1]
				estadd scalar s2 scalar(s2)
				scalar tau00 = tau3[1,1]
				estadd scalar tau00 scalar(tau00)
				estadd scalar gm_grade = `gm_grade'
				estadd scalar gm_year = `gm_year'
				
				est save hlm`v'est`s'1 , replace
				
				local dev = e(dev)
				local df = e(df)
		
				local devtest "dev(`dev' `df')"
				
			}
			else {
					
				local devest ""
					
			}
			
			// MODEL 2
			// linear grade trends

			noi di _c "..m2"	
			
			cap hlm hlm3 y int(int(rand) gradec(rand) gy1-gy`gymax' rand) rand, ///
				pwgt(wgt) cmd(hlm`v'`s'2) run replace robust `devtest'  
						
			if _rc == 0 {
				
				mat tau2 = e(tau2)
				mat tau3 = e(tau3)
				scalar s2 = tau2[1,1]
				estadd scalar s2 scalar(s2)
				scalar tau00 = tau3[1,1]
				estadd scalar tau00 scalar(tau00)
				scalar tau11 = tau3[2,2]
				estadd scalar tau11 scalar(tau11)
				scalar tau10 = tau3[1,2]/(sqrt(tau3[1,1])*sqrt(tau3[2,2]))
				estadd scalar tau10 scalar(tau10)		
				estadd scalar gm_grade = `gm_grade'
				estadd scalar gm_year = `gm_year'
				
				est save hlm`v'est`s'2 , replace
				
				local devtest "dev(`dev' `df')"
				
			}
			else {
					
				local devest ""
					
			}

			// MODEL 3
			// add linear year trends to the grade trend model

			noi di _c "..m3"	
			
			cap hlm hlm3 y int(int(rand) gradec(rand) yearc(rand) gy1-gy`gymax' rand) rand, ///
				pwgt(wgt) cmd(hlm`v'`s'3) run replace robust `devtest' 
				
			if _rc == 0 {
				
				mat tau2 = e(tau2)
				mat tau3 = e(tau3)
				scalar s2 = tau2[1,1]
				estadd scalar s2 scalar(s2)
				scalar tau00 = tau3[1,1]
				estadd scalar tau00 scalar(tau00)
				scalar tau11 = tau3[2,2]
				estadd scalar tau11 scalar(tau11)
				scalar tau10 = tau3[1,2]/(sqrt(tau3[1,1])*sqrt(tau3[2,2]))
				estadd scalar tau10 scalar(tau10)		
				scalar tau22 = tau3[3,3]
				estadd scalar tau22 scalar(tau22)
				scalar tau20 = tau3[3,1]/(sqrt(tau3[1,1])*sqrt(tau3[3,3]))
				scalar tau21 = tau3[3,2]/(sqrt(tau3[2,2])*sqrt(tau3[3,3]))
				estadd scalar tau20 scalar(tau20)
				estadd scalar tau21 scalar(tau21)
				estadd scalar gm_grade = `gm_grade'
				estadd scalar gm_year = `gm_year'
				
				est save hlm`v'est`s'3	, replace
				
			}	
		
			scalar drop _all
			
			cap erase tauvc.dat
			cap erase hlm`v'2.mdm
			cap erase hlm`v'2.mdmt
			cap erase hlm`v'2.dta
			cap erase hlm`v'2_mdmvars.dta
			cap erase hlm`v'2.sts
			capture erase hlm`v'.mdm
			capture erase creatmdm.mdmt
			capture erase hlm`v'.mdmt
			capture erase hlm`v'.dta
			capture erase hlm`v'_mdmvars.dta
			capture erase hlm`v'.sts
			capture erase hlm`v'`s'0.hlm
			capture erase hlm`v'`s'1.hlm
			capture erase hlm`v'`s'2.hlm
			capture erase hlm`v'`s'3.hlm
			capture erase newcmd.hlm
	
		restore
	} // subjects
} // states

// additional clean-up
foreach s in math ela {
	foreach m in ri rc rc2 {
		foreach y in hetop scalescore {
			capture erase hlm`s'`m'_res2.dta
			capture erase hlm`s'`m'_res3.dta
		}
	}
}

win man close viewer _all

