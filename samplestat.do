import delimited IFLS2000_main.csv, clear

********** Column 1: Upper secondary or higher (treatment group) **********

* Log hourly wages
* Years of education 
* Distance to school (km)
* Distance to health post (km)
* Fees per continuing student (rupiah)
* Age
summarize lwages edu dist_sec dist_health exp ar09 if edu >= 10
* Religion Protestant
*   Catholic
*   Other
*   Muslim
summarize protestant catholic religion_other muslim if edu >= 10
* Father uneducated
*   elementary
*   secondary and higher
*   missing
summarize une_p ele_p sec_p missing_p  if edu >= 10
* Mother uneducated
*   elementary
*   secondary and higher
*   missing
summarize une_m ele_m sec_m missing_m  if edu >= 10
* Rural household
* North Sumatra
* West Sumatra
* South Sumatra
* Lampung
* Jakarta
* Central Java
* Yogyakarta
* East Java
* Bali
* West Nussa Tengara
* South Kalimanthan
* South Sulawesi
summarize rural n_sumatra w_sumatra s_sumatra lampung jakarta c_java yogyakarta ///
 e_java bali w_nussa_tengara s_kalimanthan s_sulawesi  if edu >= 10
 
 ********** Column 2: Less than upper secondary (control group) **********

* Log hourly wages
* Years of education 
* Distance to school (km)
* Distance to health post (km)
* Fees per continuing student (rupiah)
* Age
summarize lwages edu dist_sec dist_health exp ar09 if edu < 10
* Religion Protestant
*   Catholic
*   Other
*   Muslim
summarize protestant catholic religion_other muslim if edu < 10
* Father uneducated
*   elementary
*   secondary and higher
*   missing
summarize une_p ele_p sec_p missing_p  if edu < 10
* Mother uneducated
*   elementary
*   secondary and higher
*   missing
summarize une_m ele_m sec_m missing_m  if edu < 10
* Rural household
* North Sumatra
* West Sumatra
* South Sumatra
* Lampung
* Jakarta
* Central Java
* Yogyakarta
* East Java
* Bali
* West Nussa Tengara
* South Kalimanthan
* South Sulawesi
summarize rural n_sumatra w_sumatra s_sumatra lampung jakarta c_java yogyakarta ///
 e_java bali w_nussa_tengara s_kalimanthan s_sulawesi  if edu < 10
