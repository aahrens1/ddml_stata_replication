clear all 

if ("`c(os)'"=="MacOSX") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
	adopath + "/Users/kahrens/MyProjects/ddml_simulations/"
	cd "/Users/kahrens/MyProjects/ddml_simulations/"
}
else {
	cd "/cluster/home/kahrens/ddml_simulations/"
	adopath + "/cluster/home/kahrens/ddml_simulations/"
	cap python set exec /cluster/apps/nss/gcc-8.2.0/python/3.9.9/x86_64/bin/python3	
	cap set processors 1
}

cd sim_SJ

cap log close
cap log using log/log_`1'_`2'_`3'.txt, text replace

which ddml
which pystacked

cap set seed `1'
di "`1' `2' `3' `4' `5'"

global reps = 10
global each_reps = 10 

* options
global pyopt njobs(4)  
global rflow max_depth(10) max_features(sqrt) n_estimators(500)
global rfmid max_depth(6) max_features(sqrt) n_estimators(500)
global rfhigh max_depth(2) max_features(sqrt) n_estimators(500)
global gradlow max_depth(3) n_estimators(1000) learning_rate(0.3) validation_fraction(0.2) n_iter_no_change(5) tol(0.01)
global gradmid max_depth(3) n_estimators(1000) learning_rate(0.1) validation_fraction(0.2) n_iter_no_change(5) tol(0.01)
global gradhigh max_depth(3) n_estimators(1000) learning_rate(0.01) validation_fraction(0.2) n_iter_no_change(5) tol(0.01)
global nnetopt hidden_layer_sizes(20 20)

* program 
cap program drop mysim
program define mysim , rclass

	syntax [, 	///
				obs(integer 500) ///
				DGP0(integer 1) ///
				beta(real 0.5) ///
				NOEST ///
				p(integer -1) ///
				cy(real -1) ///
				cd(real -1) ///
				folds(integer 20) ///
			]
				
	clear
	set obs `obs'
	mat Sigma = (1,0\0,1)
	drawnorm e v, cov(Sigma)
	
	** dimension of x
	local p=50
	if `dgp0'==6 {
		local dgp = 2
		local p = 13
	}
	else if `dgp0'==7 {
		local dgp = 4
		local p = 7
	}
	else {
		local dgp = `dgp0'
	}
	
	* correlation matrix
	mat Eps = J(`p',`p',.)
	forvalues i = 1(1)`p' {
		forvalues j = `i'(1)`p' {
			mat Eps[`i',`j']=0.5^(abs(`i'-`j'))
			mat Eps[`j',`i']=0.5^(abs(`i'-`j'))
		}
	}
	drawnorm x1-x`p', cov(Eps)

	* select DGP
	if (`dgp'==1) {
		gen double Xall = 0
		forvalues j = 1(1)`p' {
			replace Xall=Xall +(1/`j')^2*x`j'
		}
		if (`cy'<0) local cy = .515 //.342 
		if (`cd'<0) local cd = .827 //.543
	}	
	if (`dgp'==2) {
		gen double Xall = x1*x2 + x3^2 + x4*x5 + x6*x7 + x8*x9 + x10 + x11^2 + x12*x13
		if (`cy'<0) local cy = .1633 //.11
		if (`cd'<0) local cd = .275 //.18
	}
	if (`dgp'==3) {
		gen double Xall = (x1>.3)*(x2>0)*(x3>(-1))  
		if (`cy'<0) local cy = 1.42 //.94
		if (`cd'<0) local cd = 2.29 //1.48
	}
	if (`dgp'==4) {
		gen double Xall = x1 + sqrt(abs(x2)) + sin(x3) + .3*x4*x5 + x6 + .3*x7^2
		if (`cy'<0) local cy = .335 //.22
		if (`cd'<0) local cd = .547  //.36
	}	
	if (`dgp'==5) {
		gen double Xall = 0
		forvalues j = 1(1)`p' {
			replace Xall=Xall+(.9)^(`j')*x`j'
		}
		if (`cy'<0) local cy = .189 //.125
		if (`cd'<0) local cd = .298 //.195
	}	
	
	* heterosk.
	gen double sigd=(1+`cd'*Xall)^2
	sum sigd, meanonly
	local sigd_m=r(mean)
	replace sigd = sqrt(sigd/`sigd_m')
	
	gen D = `cd'*Xall + v*sigd
	
	gen double sigy=(1+`beta'*D+`cy'*Xall)^2
	sum sigy, meanonly
	local sigy_m=r(mean)
	replace sigy = sqrt(sigy/`sigy_m')	
	
	gen Y = `beta'*D + `cy'*Xall + e*sigy
	
	reg Y Xall
	local ry = e(r2)
	reg D Xall
	local rd = e(r2)
	
	reg Y D Xall, robust
	local oracle_b = _b[D]
	local oracle_se = _se[D]	
	
	drop v e Xall

	** fold var
	gen u = runiform() 
	egen fid = cut(u), group(`folds')
	
	** number of learners
	local nlearners = 14
	
	** 5th order poly
	forvalues i =1(1)5 {
		forvalues j = 1(1)`p' {
			gen poly_x`j'_`i'=(x`j')^(`i')
		}
	}
	
	if "`noest'"=="" { 
 		
		*** ddml partial w/ stacking
		timer on 1
		ddml init partial, foldvar(fid) 
		ddml E[Y|X], l(Y0_py): pystacked Y x*       || ///
						m(ols) || ///
						m(lassocv) || ///
						m(ridgecv) || ///
						m(lassocv) xvars(poly*) || ///
						m(ridgecv) xvars(poly*) || ///
						m(lassocv) xvars(c.(x*)##c.(x*)) || ///
						m(ridgecv) xvars(c.(x*)##c.(x*)) || ///
						m(rf) opt($rflow) || ///
						m(rf) opt($rfmid) || ///
						m(rf) opt($rfhigh) || ///
						m(gradboost) opt($gradlow) || ///
						m(gradboost) opt($gradmid) || ///
						m(gradboost) opt($gradhigh) || ///
						m(nnet) opt($nnetopt) , ///
						$pyopt  
		ddml E[D|X], l(D0_py): pystacked D x*       || ///
						m(ols) || ///
						m(lassocv) || ///
						m(ridgecv) || ///
						m(lassocv) xvars(poly*) || ///
						m(ridgecv) xvars(poly*) || ///
						m(lassocv) xvars(c.(x*)##c.(x*)) || ///
						m(ridgecv) xvars(c.(x*)##c.(x*)) || ///
						m(rf) opt($rflow) || ///
						m(rf) opt($rfmid) || ///
						m(rf) opt($rfhigh) || ///
						m(gradboost) opt($gradlow) || ///
						m(gradboost) opt($gradmid) || ///
						m(gradboost) opt($gradhigh) || ///
						m(nnet) opt($nnetopt) , ///
						$pyopt  
		ddml crossfit , shortstack poolstack
		ddml estimate, robust
		timer off 1
				
		foreach final in nnls1 singlebest {
				
		ddml estimate, robust finalest(`final')
	
			// regular stacking results
			ddml estimate, mname(m0) spec(st) replay  
			local ddml_st_`final'_b = _b[D]
			local ddml_st_`final'_se = _se[D]
			// shortstacking results
			ddml estimate, mname(m0) spec(ss) replay 
			local ddml_ss_`final'_b = _b[D]
			local ddml_ss_`final'_se = _se[D]	
			// poolstacking results
			ddml estimate, mname(m0) spec(ps) replay  
			local ddml_ps_`final'_b = _b[D]
			local ddml_ps_`final'_se = _se[D]	
		
			tempname By Bd
			// pystacked weights
			ddml extract, show(stweights)
			mat `By'=r(Y0_py_w_mn)
			mat `Bd'=r(D0_py_w_mn)
			mat list `Bd'
			mat list `By'
			forvalues i = 1(1)`nlearners' {
				local stw_`final'_d`i' = el(`Bd',`i',2)
				local stw_`final'_y`i' = el(`By',`i',2)
			}
			// shortstacked weights
			ddml extract, show(ssweights)
			mat `By'=r(Y_Y_ss)
			mat `Bd'=r(D_D_ss)
			mat list `Bd'
			mat list `By'
			forvalues i = 1(1)`nlearners' {
				local ssw_`final'_d`i' = el(`Bd',`i',2)
				local ssw_`final'_y`i' = el(`By',`i',2)
			}
			// poolstacked weights
			ddml extract, show(psweights)
			mat `By'=r(Y_Y_ps)
			mat `Bd'=r(D_D_ps)
			mat list `Bd'
			mat list `By'
			forvalues i = 1(1)`nlearners' {
				local psw_`final'_d`i' = el(`Bd',`i',2)
				local psw_`final'_y`i' = el(`By',`i',2)
			}

		}
		
		*** mspe
		foreach i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 {
			cap drop sqpe*
			gen sqpe_y=(Y0_py_L`i'_)^2
			gen sqpe_d=(D0_py_L`i'_)^2
			sum sqpe_y, meanonly
			local mspey`i' = r(mean)
			sum sqpe_d, meanonly
			local msped`i' = r(mean)
		}
		
		*** individual learners
		forvalues l = 1(1)`nlearners' {
				qui ddml estimate, y(Y0_py_L`l'_1) d(D0_py_L`l'_1) robust
				local ddml_x_`l'_b = _b[D]
				local ddml_x_`l'_se = _se[D]
		}	

		
		ddml drop	
				
				
		*** other learners		
		
		local ols_b = .
		local ols_se = .
		if (`p'<`obs') {
			timer on 4
			reg Y D x*, robust
			timer off 4
			local ols_b = _b[D]
			local ols_se = _se[D]
		}		
		
		timer on 5
		pdslasso Y D (x*), robust
		timer off 5
		local pds_b = _b[D]
		local pds_se = _se[D]
	
		timer on 51
		pdslasso Y D (c.(x*)##c.(x*)), robust
		timer off 51
		local pds2_b = _b[D]
		local pds2_se = _se[D]
	
		timer on 52
		pdslasso Y D (poly*), robust
		timer off 52
		local pds3_b = _b[D]
		local pds3_se = _se[D]
	
		*** sd
		foreach var of varlist Y D {
			sum `var'
			local sd_`var' = r(sd)
		}
		
		*** return ************************************
		ereturn clear
		* coefficients per final learner
		foreach final in nnls1 singlebest {
			foreach v in st ss ps {
				return scalar ddml_`v'_`final'_b = `ddml_`v'_`final'_b'
				return scalar ddml_`v'_`final'_se = `ddml_`v'_`final'_se'
			}
		}
		forvalues l = 1(1)`nlearners' {
			return scalar ddml_x_`l'_b = `ddml_x_`l'_b' 
			return scalar ddml_x_`l'_se =`ddml_x_`l'_se'
		}	
		* stacking weights per final learner
		foreach final in nnls1 singlebest {
			foreach type in stw ssw psw {
				forvalues i = 1(1)`nlearners' {
					return scalar `type'_`final'_d`i' = ``type'_`final'_d`i'' 
					return scalar `type'_`final'_y`i' = ``type'_`final'_y`i''
				}
			}
		}
		foreach v in oracle ols pds pds2 pds3 {
			return scalar `v'_b = ``v'_b'
			return scalar `v'_se = ``v'_se'
		}
		* mspe
		foreach i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 {
			return scalar mspey`i' = `mspey`i''
			return scalar msped`i' = `msped`i''
		}
		*** return 
		return scalar obs = `obs'
		return scalar p = `p'
		return scalar dgp = `dgp'
		return scalar dgp0 = `dgp0'
		return scalar beta = `beta'
		return scalar ry = `ry'
		return scalar rd = `rd'
		return scalar folds = `folds'
	} 
	else {
		ereturn clear
		return scalar ry = `ry'
		return scalar rd = `rd'
		return scalar cd = `cd'
		return scalar cy = `cy'
	}
end

timer clear
timer on 90
forvalues i = 1(1)$reps {
	clear
	simulate, reps($each_reps): mysim, obs(`2') dgp(`3') folds(`4')
	gen seed = `1'
	save out/out_`1'_`2'_`3'_`4'_`i'.dta, replace
	cap rm "error.txt"
}
timer off 90
timer list

cap log close
