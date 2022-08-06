*******************************************************
* Construct distance to nearest secondary school using 
* data from Service Availability Roster (SAR)
*******************************************************
clear
use sar.dta

* restrict sample to junior high schools (SMP) and senior high schools (SMU)
    keep if facility == "SMP" | facility == "SMU"
	
* exclude empty community id
    drop if commid00 == ""
	
* exclude invalid facility code
    drop if fascode == "" | fcode00 == "XXXX"
    keep if commid00 == substr(fcode00,1,4)
	
* exclude missing values and extreme outliers of distance
    keep if km2kdx == 1
    drop if km2kd > 999
	
* select the secondary school closest to the office of the head of the community
    collapse (min) dist_sec=km2kd, by (commid00)
	
save distance.dta, replace
	
***********************************************************
* Construct tuition fee of cheapest secondary school using
* data from School Questionnaire, Section E
***********************************************************
clear
use schl_e.dta

* restrict sample to fees for continuing students in school year
    gen form = substr(etype,1,3)
    drop if form != "E1B"
	
* exclude missing values of expenditure
    drop if !(e1ax == 1 | e1bx == 1)
	
* define expenditure as average over 99/00 and 00/01 if both observations are available
    gen exp = e1a if e1bx != 1
    replace exp = e1b if e1ax != 1
    replace exp = (e1a+e1b)/2 if e1ax == 1 & e1bx == 1

* define tuition fee as sum of different type of expenditure
    collapse (sum) exp, by (fcode00)

* generate community id from first 4 digits of facility code
    gen commid00 = substr(fcode00,1,4)

* select the secondary school with lowest tuition fee within the community
    collapse (min) exp, by (commid00)

save tuition.dta, replace


****************************************************
* Construct distance to nearest health post
* using data from Service Availability Roster (SAR)
****************************************************
clear
use sar.dta

* restrict sample to integrated health posts (posyandu)
    keep if facility == "posyandu"

* exclude empty community id
    drop if commid00 == ""

* exclude missing values and extreme outliers of distance
    keep if km2kdx == 1
    drop if km2kd > 999

* select the health post closest to the office of the head of the community
    collapse (min) dist_health=km2kd, by (commid00)

save healthpost.dta, replace


*************************************************************************************
* Construct age and religion using data from IFLS3 Control Book (Book K), Section AR
*************************************************************************************
clear
use bk_ar1.dta

* restrict sample to males aged 25-60
    keep if ar07 == 1
    keep if (ar09>=25) & (ar09<=60)

* generate dummies for religion
    gen protestant = 0
    replace protestant = 1 if ar15 == 2
    gen catholic = 0
    replace catholic = 1 if ar15 == 3
    gen muslim = 0
    replace muslim = 1 if ar15 == 1
    gen religion_other = 0
    replace religion_other = 1 if !(ar15 == 1 | ar15 == 2 | ar15 == 3)

keep pidlink hhid00 ar09 ar15 protestant catholic muslim religion_other


*******************************************************************
* Construct province of residence and indicator of rural residence
* using data from IFLS3 Control Book (Book K), Section SC
*******************************************************************
merge m:1 hhid00 using bk_sc.dta, keep(match) nogen

* generate community id from code of province and code of kabupaten
    tostring sc01, gen(sc01_str)
    tostring sc02, gen(sc02_str)
    replace sc02_str = "0"+sc02_str if strlen(sc02_str)==1
    gen commid00 = sc01_str+sc02_str
	
* generate indicator of rural residence 
	gen rural = 0
replace rural = 1 if sc05 == 2

* generate dummies for province
    gen n_sumatra = 0
    replace n_sumatra = 1 if sc01 == 12
    gen w_sumatra = 0
    replace w_sumatra = 1 if sc01 == 13
    gen s_sumatra = 0
    replace s_sumatra = 1 if sc01 == 16
    gen lampung = 0
    replace lampung = 1 if sc01 == 18
    gen jakarta = 0
    replace jakarta = 1 if sc01 == 31
    gen c_java = 0
    replace c_java = 1 if sc01 == 33
    gen yogyakarta = 0
    replace yogyakarta = 1 if sc01 == 34
    gen e_java = 0
    replace e_java = 1 if sc01 == 35
    gen bali = 0
    replace bali = 1 if sc01 == 51
    gen w_nussa_tengara = 0
    replace w_nussa_tengara = 1 if sc01 == 52
    gen s_kalimanthan = 0
    replace s_kalimanthan = 1 if sc01 == 63
    gen s_sulawesi = 0
    replace s_sulawesi = 1 if sc01 == 73

keep pidlink hhid00 ar09 ar15 protestant catholic muslim religion_other ///
    sc01 sc05 commid00 rural n_sumatra w_sumatra s_sumatra lampung jakarta ///
    c_java yogyakarta e_java bali w_nussa_tengara s_kalimanthan s_sulawesi
	
	
*****************************************************************
* Construct log hourly wages using data from Book 3A, Section TK
*****************************************************************
merge m:1 pidlink using b3a_tk1.dta, keep(match) nogen

* restrict sample to employed
    keep if tk01 == 1
	
keep pidlink hhid00 ar09 ar15 protestant catholic muslim religion_other ///
    sc01 sc05 commid00 rural n_sumatra w_sumatra s_sumatra lampung jakarta ///
    c_java yogyakarta e_java bali w_nussa_tengara s_kalimanthan s_sulawesi

merge m:1 pidlink using b3a_tk2.dta, keep(match) nogen

* exclude invalid working hours and wage
    drop if tk22a == . | tk25a1 == . | tk22a == 0 | tk25a1 == 0

* exclude extreme outliers of working hours (24*7=168)
    drop if tk22a > 168

* calculate hourly wages as monthly wage divided by (4*working hours per week)
	gen lwages = log(tk25a1/(4*tk22a))
	
keep pidlink hhid00 ar09 ar15 protestant catholic muslim religion_other ///
    sc01 sc05 commid00 rural n_sumatra w_sumatra s_sumatra lampung jakarta ///
    c_java yogyakarta e_java bali w_nussa_tengara s_kalimanthan s_sulawesi lwages
	
	
******************************************************************************
* Construct years of education and indicator of attendance of upper secondary 
* school or higher using data from Book 3A, Section DL
******************************************************************************
merge m:1 pidlink using b3a_dl1.dta, keep(match) nogen

* exclude missing values of highest education level and highest grade completed
    drop if dl06 == . | dl07 == .

* exclude other
    drop if dl06 == 10

* exclude adult education
    drop if dl06 == 11 | dl06 == 12
	
* exclude open university
    drop if dl06 == 13

* exclude Islamic school (Pesantren)
    drop if dl06 == 14
	
* generate years of education from highest education level and highest grade completed
    gen edu = 0
    replace edu = 6 if dl06 == 3 | dl06 == 4 | dl06 ==  73
    replace edu = 9 if dl06 == 5 | dl06 == 6 | dl06 ==  74
    replace edu = 12 if dl06 == 60 | dl06 == 61
    replace edu = 16 if dl06 == 62
    replace edu = edu+dl07 if dl07 != 7
    replace edu = edu+6 if dl07 == 7 & (dl06 == 2 | dl06 == 72)
    replace edu = edu+3 if dl07 == 7 & (dl06 == 3 | dl06 == 4 | dl06 ==  73 | dl06 == 5 | dl06 == 6 | dl06 ==  74)
    replace edu = edu+4 if dl07 == 7 & (dl06 == 60 | dl06 == 61)
    replace edu = edu+2 if dl07 == 7 & dl06 == 62
* generate indicator of attendance of upper secondary school or higher
    gen upsec = 0
    replace upsec = 1 if edu>=10

keep pidlink hhid00 ar09 ar15 protestant catholic muslim religion_other ///
sc01 sc05 commid00 rural n_sumatra w_sumatra s_sumatra lampung jakarta ///
c_java yogyakarta e_java bali w_nussa_tengara s_kalimanthan s_sulawesi ///
lwages edu upsec


*******************************************************************
* Construct parental education using data from Book 3B, Section BA
*******************************************************************
merge m:1 pidlink using b3b_ba0.dta, keep(match) nogen

* generate parental education (uneducated, elementary, secondary and higher, missing)
    foreach i in p m{
    gen une_`i' = 0
    replace une_`i' = 1 if ba08`i' == 2 & ba09`i' != 7
    replace une_`i' = 1 if ba08`i' == 72 & ba09`i' != 7
    replace une_`i' = 1 if ba08`i' == 98
    replace une_`i' = 1 if ba08`i' == 7 | ba08`i' == 8
    gen ele_`i' = 0
    replace ele_`i' = 1 if ba08`i' == 2 & ba09`i' == 7
    replace ele_`i' = 1 if ba08`i' == 72 & ba09`i' == 7
    replace ele_`i' = 1 if (ba08`i' == 3 | ba08`i' == 4) & ba09`i' != 7
    replace ele_`i' = 1 if ba08`i' == 73 & ba09`i' != 7
    gen sec_`i' = 0
    replace sec_`i'  = 1 if (ba08`i' == 3 | ba08`i' == 4) & ba09`i' == 7
    replace sec_`i'  = 1 if ba08`i' == 73 & ba09`i' == 7
    replace sec_`i'  = 1 if ba08`i' == 5 | ba08`i' == 6
    replace sec_`i'  = 1 if ba08`i' == 11 | ba08`i' == 13 | ba08`i' == 14
    replace sec_`i'  = 1 if ba08`i' == 60 | ba08`i' == 61 | ba08`i' == 62
    gen missing_`i' = 0
    replace missing_`i' = 1 if ba08`i' == .
    }

keep pidlink hhid00 ar09 ar15 protestant catholic muslim religion_other ///
sc01 sc05 commid00 rural n_sumatra w_sumatra s_sumatra lampung jakarta ///
c_java yogyakarta e_java bali w_nussa_tengara s_kalimanthan s_sulawesi ///
lwages edu upsec une_p ele_p sec_p missing_p une_m ele_m sec_m missing_m


merge m:1 commid00 using distance.dta, keep(match) nogen
merge m:1 commid00 using tuition.dta, keep(match) nogen
merge m:1 commid00 using healthpost.dta, keep(match) nogen

* exclude potential outliers of distance
    drop if dist_sec == 22 

save IFLS2000_main.dta, replace
export delimited using IFLS2000_main.csv, replace
