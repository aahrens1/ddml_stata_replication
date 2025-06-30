
clear all

if ("`c(os)'"=="MacOSX") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
	cd "/Users/kahrens/MyProjects/ddml_simulations/"
}
else {
	cd "/cluster/home/kahrens/ddml_simulations/"
	cap python set exec /cluster/apps/nss/gcc-8.2.0/python/3.9.9/x86_64/bin/python3
	cap set processors 1
}

global outpath /Users/kahrens/MyProjects/ddml_sjpaper/Simul/sim_SJ
global outpath2 /Users/kahrens/MyProjects/ddml_simulations/sim_SJ/output
cap mkdir $outpath
cd sim_SJ

global folder 

local files : dir out files "*.dta"
foreach f in `files' {
	cap append using out/`f'	
}
 
foreach var of varlist ddml*_b pds_b ols_b pds2_b pds3_b oracle_b {
	local var = subinstr("`var'","_b","",1)
	di "`var'"
	qui gen `var'_p = 2*normal(-abs((`var'_b)/(`var'_se)))
	qui gen `var'_abias = abs(`var'_b-beta)
	qui gen `var'_bias = (`var'_b-beta)
	qui gen `var'_sign=`var'_p<0.05
	qui gen `var'_cov = inrange(beta,`var'_b-1.96*`var'_se,`var'_b+1.96*`var'_se)
} 

drop if missing(dgp)

save sim_temp_all.dta, replace

gen byte count = 1
collapse (mean) *sign *_p *_cov *_b ssw_* stw_* psw_* mspe* rd ry *_bias (median) *abias (sum) count, by(beta dgp0 obs) 

tab dgp, summarize(ry)
tab dgp, summarize(rd)


save sim_temp.dta, replace


********* histogram

use sim_temp_all.dta, clear

keep if dgp0==5 & obs ==1000

label var oracle_bias "Infeasible Oracle"
 
twoway (hist oracle_bias   , width(0.01) start(-0.4) color(midgreen%30)  ) ///
		(hist ddml_x_2_bias  , width(0.01) start(-0.4) color(red%40)) ///
		(hist ddml_x_14_bias  , width(0.01) start(-0.4) color(navy%60))  , ///
		scheme(tab1) ///
		legend(order(1 "Oracle" 2 "CV-Lasso" 3 "Neural net")  size(medsmall)) ///
		xlabel(-.2(0.1).3,labsize(medsmall)) xsc(r(-.4 .2)  ) ///
		ylabel(,labsize(medsmall))  ///
		ytitle("",size(medsmall)) ///
		xtitle("Bias",size(medsmall))  
graph export $outpath/example1.png, replace
graph export $outpath2/example1.png, replace

use sim_temp_all.dta, clear

keep if dgp0==3 & obs ==1000

label var oracle_bias "Infeasible Oracle"
 
twoway (hist oracle_bias   , width(0.01) start(-0.2) color(midgreen%30)  ) ///
		(hist ddml_x_2_bias  , width(0.01) start(-0.2) color(red%40)) ///
		(hist ddml_x_14_bias  , width(0.01) start(-0.2) color(navy%60))  , ///
		scheme(tab1) ///
		legend(order(1 "Oracle" 2 "CV-Lasso" 3 "Neural net")  size(medsmall)) ///
		xlabel(-.2(0.1).3,labsize(medsmall)) xsc(r(-.2 .3)  ) ///
		ylabel(,labsize(medsmall))  ///
		ytitle("",size(medsmall)) ///
		xtitle("Bias",size(medsmall))  
graph export $outpath/example2.png, replace
graph export $outpath2/example2.png, replace

********* bias and coverage 

use sim_temp.dta, clear

keep *_bias *_abias *_cov dgp obs

foreach var of varlist *_bias {
	local var2 = subinstr("`var'","ddml_x_","",1)
	local var2 = subinstr("`var2'","_bias","",1)
	local var2 = subinstr("`var2'","ols_","",1)
	local var2 = subinstr("`var2'","pds_","",1)
	rename `var' bias_`var2'
}
foreach var of varlist *_abias {
	local var2 = subinstr("`var'","ddml_x_","",1)
	local var2 = subinstr("`var2'","_abias","",1)
	local var2 = subinstr("`var2'","ols_","",1)
	local var2 = subinstr("`var2'","pds_","",1)
	rename `var' abias_`var2'
}
foreach var of varlist *cov {
	local var2 = subinstr("`var'","ddml_x_","",1)
	local var2 = subinstr("`var2'","_cov","",1)
	local var2 = subinstr("`var2'","ols_","",1)
	local var2 = subinstr("`var2'","pds_","",1)
	rename `var' cov_`var2'
}
	
reshape long bias_ abias_ cov_, ///
				i(dgp obs) j(estimator ) string
	
*reshape wide bias_ abias_ cov_,	i( estimator) j(dgp )

*graph bar (asis) abias_1 , over(estimator) by(obs)
*graph bar (asis) cov_1 , over(estimator) by(obs)

local v 
foreach i of numlist 2 3 4 5 7 {
	local v `v' abias_`i' bias_`i' cov_`i'
}

reshape wide abias_ bias_ cov_, j(dgp ) i(estimator  obs)
reshape wide `v', j(obs ) i(estimator   )
				
 
gen sortid=1
replace sortid=0.9 if estimator=="oracle"
replace sortid=2 if estimator=="pds"
replace sortid=2.3 if estimator=="pds2"
replace sortid=2.2 if estimator=="pds3"
replace sortid=4 if estimator=="1"
replace sortid=5 if estimator=="2"
replace sortid=6 if estimator=="3"
replace sortid=7 if estimator=="4"
replace sortid=8 if estimator=="5"
replace sortid=9 if estimator=="6"
replace sortid=10 if estimator=="7"
replace sortid=11 if estimator=="8"
replace sortid=12 if estimator=="9"
replace sortid=12.1 if estimator=="10"
replace sortid=12.2 if estimator=="11"
replace sortid=12.3 if estimator=="12"
replace sortid=12.4 if estimator=="13"
replace sortid=12.99999 if estimator=="14"
replace sortid=13.1 if estimator=="0"
replace sortid=13.2 if estimator=="sb"
replace sortid=14 if estimator=="ia"
replace sortid=15.1 if estimator=="ss"
replace sortid=15.2 if estimator=="mse"
replace sortid = 14.1  if estimator=="ddml_ss_nnls1"
replace sortid = 14.2 if estimator=="ddml_ss_singlebest"
replace sortid = 13.1 if estimator=="ddml_st_nnls1"
replace sortid = 13.2 if estimator=="ddml_st_singlebest"
replace sortid = 15.1 if estimator=="ddml_ps_nnls1"
replace sortid = 15.2 if estimator=="ddml_ps_singlebest"
**
replace estimator = "OLS (Base)" if estimator=="ols"
replace estimator = "Oracle" if estimator=="oracle"
replace estimator = "PDS-Lasso (Base)" if estimator=="pds"
replace estimator = "PDS-Lasso (Poly 2 + Inter.)" if estimator=="pds2"
replace estimator = "PDS-Lasso (Poly 5)" if estimator=="pds3"
replace estimator = "Stacking: CLS" if estimator=="0"
replace estimator = "OLS" if estimator=="1"
replace estimator = "Lasso with CV (Base)" if estimator=="2"
replace estimator = "Ridge with CV (Base)" if estimator=="3"
replace estimator = "Lasso with CV (Poly 5)" if estimator=="4"
replace estimator = "Ridge with CV (Poly 5)" if estimator=="5"
replace estimator = "Lasso with CV (Poly 2 + Inter.)" if estimator=="6"
replace estimator = "Ridge with CV (Poly 2 + Inter.)" if estimator=="7"
replace estimator = "Random forest (Low)" if estimator=="8"
replace estimator = "Random forest (Medium)" if estimator=="9"
replace estimator = "Random forest (High)" if estimator=="10"
replace estimator = "Gradient boosting (Low)" if estimator=="11"
replace estimator = "Gradient boosting (Medium)" if estimator=="12"
replace estimator = "Gradient boosting (High)" if estimator=="13"
replace estimator = "Neural net" if estimator=="14"
replace estimator = "Short-stacking: CLS" if estimator=="ddml_ss_nnls1"
replace estimator = "Short-stacking: Single best" if estimator=="ddml_ss_singlebest"
replace estimator = "Stacking: CLS" if estimator=="ddml_st_nnls1"
replace estimator = "Stacking: Single best" if estimator=="ddml_st_singlebest"
replace estimator = "Pooled stacking: CLS" if estimator=="ddml_ps_nnls1"
replace estimator = "Pooled stacking: Single best" if estimator=="ddml_ps_singlebest"

sort sortid
 
tostring cov* abias*, force format(%9.3f) replace	
foreach var of varlist cov* abias* {
	replace `var'="0.\phantom{000}" if `var'=="0.000"
}

gen gap1=""
gen gap2=""
gen gap3=""
sort sortid

texsave gap1 estimator abias_5100 cov_5100 abias_51000 cov_51000 ///
			gap2 abias_2100 cov_2100 abias_21000 cov_21000 ///
			gap3 abias_3100 cov_3100 abias_31000 cov_31000 ///
			using $outpath/bias_1.tex  , dataonly replace nofix
texsave gap1 estimator abias_4100 cov_4100 abias_41000 cov_41000 ///
			gap2 abias_7100 cov_7100 abias_71000 cov_71000 using $outpath/bias_2.tex  , dataonly replace nofix
	
texsave gap1 estimator abias_5100 cov_5100 abias_51000 cov_51000 ///
			gap2 abias_2100 cov_2100 abias_21000 cov_21000 ///
			gap3 abias_3100 cov_3100 abias_31000 cov_31000 ///
			using $outpath2/bias_1.tex  , dataonly replace nofix
texsave gap1 estimator abias_4100 cov_4100 abias_41000 cov_41000 ///
			gap2 abias_7100 cov_7100 abias_71000 cov_71000 using $outpath2/bias_2.tex  , dataonly replace nofix
				
	
********* pystacked weights 

use sim_temp.dta, clear

keep dgp obs stw_nnls1_* ssw_nnls1_* psw_nnls1_*

reshape long stw_nnls1_ ssw_nnls1_  psw_nnls1_, i(obs dgp) j(estimator) string

gen equation = "D" if regexm(estimator,"d")
replace equation = "Y" if regexm(estimator,"y")

gen estimator_desc = ""
replace estimator_desc = "OLS" if regexm(estimator,"1") & !regexm(estimator,"10") 
replace estimator_desc = "Lasso with CV (Base)" if regexm(estimator,"2") 
replace estimator_desc = "Ridge with CV (Base)" if regexm(estimator,"3")
replace estimator_desc = "Lasso with CV (Poly 5)" if regexm(estimator,"4")
replace estimator_desc = "Ridge with CV (Poly 5)" if regexm(estimator,"5")
replace estimator_desc = "Lasso with CV (Poly 2 + Inter.)" if regexm(estimator,"6")
replace estimator_desc = "Ridge with CV (Poly 2 + Inter.)" if regexm(estimator,"7")
replace estimator_desc = "Random forest (Low regularization)" if regexm(estimator,"8")
replace estimator_desc = "Random forest (Medium regularization)" if regexm(estimator,"9")
replace estimator_desc = "Random forest (High regularization)" if regexm(estimator,"10")
replace estimator_desc = "Gradient boosting (Low regularization)" if regexm(estimator,"11")
replace estimator_desc = "Gradient boosting (Medium regularization)" if regexm(estimator,"12")
replace estimator_desc = "Gradient boosting (High regularization)" if regexm(estimator,"13")
replace estimator_desc = "Neural net" if regexm(estimator,"14")

gen sortid = substr(estimator,2,2)
destring sortid, replace
drop estimator 

gen dgp_desc =dgp

tostring stw_nnls1_* , replace format(%9.3f) force
foreach var of varlist stw_nnls1_* {
	replace `var'="0.\phantom{000}" if `var'=="0.000"
}

drop dgp0

reshape wide stw_nnls1_ ssw_nnls1_ psw_nnls1_  , i(sortid obs  dgp_desc estimator_desc) j(equation) string
reshape wide stw_nnls1_Y ssw_nnls1_Y psw_nnls1_Y  stw_nnls1_D ssw_nnls1_D psw_nnls1_D , i(sortid obs  estimator_desc) j(dgp_desc)
ds stw_nnls1_* ssw_nnls1_* psw_nnls1_*
foreach var of varlist `r(varlist)' {
	rename `var' `var'_
}
ds stw_nnls1_* ssw_nnls1_* psw_nnls1_*
reshape wide `r(varlist)'  , i(sortid  estimator_desc) j( obs)

sort sortid

gen gap = ""

foreach stackw in ssw stw psw {
		
	texsave gap estimator_desc `stackw'_nnls1_Y5_100 `stackw'_nnls1_Y5_1000  ///
							`stackw'_nnls1_Y2_100 `stackw'_nnls1_Y2_1000  ///
							`stackw'_nnls1_Y3_100 `stackw'_nnls1_Y3_1000  ///
							`stackw'_nnls1_Y4_100 `stackw'_nnls1_Y4_1000  ///
							`stackw'_nnls1_Y7_100 `stackw'_nnls1_Y7_1000  ///
							using $outpath/`stackw'_Y.tex   , dataonly replace nofix
	texsave gap estimator_desc `stackw'_nnls1_D5_100 `stackw'_nnls1_D5_1000  ///
							`stackw'_nnls1_D2_100 `stackw'_nnls1_D2_1000  ///
							`stackw'_nnls1_D3_100 `stackw'_nnls1_D3_1000  ///
							`stackw'_nnls1_D4_100 `stackw'_nnls1_D4_1000  ///
							`stackw'_nnls1_D7_100 `stackw'_nnls1_D7_1000  ///
							using $outpath/`stackw'_D.tex  , dataonly replace nofix

							
	texsave gap estimator_desc `stackw'_nnls1_Y5_100 `stackw'_nnls1_Y5_1000  ///
							`stackw'_nnls1_Y2_100 `stackw'_nnls1_Y2_1000  ///
							`stackw'_nnls1_Y3_100 `stackw'_nnls1_Y3_1000  ///
							`stackw'_nnls1_Y4_100 `stackw'_nnls1_Y4_1000  ///
							`stackw'_nnls1_Y7_100 `stackw'_nnls1_Y7_1000  ///
							using $outpath2/`stackw'_Y.tex   , dataonly replace nofix
	texsave gap estimator_desc `stackw'_nnls1_D5_100 `stackw'_nnls1_D5_1000  ///
							`stackw'_nnls1_D2_100 `stackw'_nnls1_D2_1000  ///
							`stackw'_nnls1_D3_100 `stackw'_nnls1_D3_1000  ///
							`stackw'_nnls1_D4_100 `stackw'_nnls1_D4_1000  ///
							`stackw'_nnls1_D7_100 `stackw'_nnls1_D7_1000  ///
							using $outpath2/`stackw'_D.tex  , dataonly replace nofix
							
							
}

********* MSPE

use sim_temp.dta, clear

keep dgp obs mspe*

reshape long mspe , i(obs dgp) j(estimator) string

gen equation = "D" if regexm(estimator,"d")
replace equation = "Y" if regexm(estimator,"y")

gen estimator_desc = ""
replace estimator_desc = "OLS" if regexm(estimator,"1") & !regexm(estimator,"10") 
replace estimator_desc = "Lasso with CV (Base)" if regexm(estimator,"2") 
replace estimator_desc = "Ridge with CV (Base)" if regexm(estimator,"3")
replace estimator_desc = "Lasso with CV (Poly 5)" if regexm(estimator,"4")
replace estimator_desc = "Ridge with CV (Poly 5)" if regexm(estimator,"5")
replace estimator_desc = "Lasso with CV (Poly 2 + Inter.)" if regexm(estimator,"6")
replace estimator_desc = "Ridge with CV (Poly 2 + Inter.)" if regexm(estimator,"7")
replace estimator_desc = "Random forest (Low regularization)" if regexm(estimator,"8")
replace estimator_desc = "Random forest (Medium regularization)" if regexm(estimator,"9")
replace estimator_desc = "Random forest (High regularization)" if regexm(estimator,"10")
replace estimator_desc = "Gradient boosting (Low regularization)" if regexm(estimator,"11")
replace estimator_desc = "Gradient boosting (Medium regularization)" if regexm(estimator,"12")
replace estimator_desc = "Gradient boosting (High regularization)" if regexm(estimator,"13")
replace estimator_desc = "Neural net" if regexm(estimator,"14")
replace estimator_desc = "Stacking: CLS" if estimator=="d0" | estimator == "y0"
replace estimator_desc = "Stacking: Single-best" if estimator=="dsb" | estimator == "ysb"
drop if estimator_desc ==""

gen sortid = substr(estimator,2,2)
replace sortid = "16" if sortid =="0"
replace sortid = "17" if sortid =="sb"
destring sortid, replace
drop estimator 

gen dgp_desc =dgp

tostring mspe* , replace format(%9.3f) force
foreach var of varlist mspe* {
	replace `var'="0.\phantom{000}" if `var'=="0.000"
}

drop dgp0

reshape wide mspe  , i(sortid  dgp_desc estimator_desc obs ) j(equation) string
reshape wide mspeY  mspeD  , i(sortid   estimator_desc obs) j(dgp_desc)
ds mspe*
reshape wide `r(varlist)'  , i(sortid  estimator_desc) j( obs)

sort sortid

gen gap = ""

texsave gap estimator_desc mspeY5100 mspeY51000  ///
						mspeY2100 mspeY21000  ///
						mspeY3100 mspeY31000  ///
						mspeY4100 mspeY41000  ///
						mspeY7100 mspeY71000  ///
						using $outpath/mspe_Y.tex , dataonly replace  nofix
texsave gap estimator_desc mspeD5100 mspeD51000  ///
						mspeD2100 mspeD21000  ///
						mspeD3100 mspeD31000  ///
						mspeD4100 mspeD41000  ///
						mspeD7100 mspeD71000  ///
						using $outpath/mspe_D.tex , dataonly replace nofix
 
texsave gap estimator_desc mspeY5100 mspeY51000  ///
						mspeY2100 mspeY21000  ///
						mspeY3100 mspeY31000  ///
						mspeY4100 mspeY41000  ///
						mspeY7100 mspeY71000  ///
						using $outpath2/mspe_Y.tex , dataonly replace  nofix
texsave gap estimator_desc mspeD5100 mspeD51000  ///
						mspeD2100 mspeD21000  ///
						mspeD3100 mspeD31000  ///
						mspeD4100 mspeD41000  ///
						mspeD7100 mspeD71000  ///
						using $outpath2/mspe_D.tex , dataonly replace nofix
 
