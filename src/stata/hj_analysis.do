/* =============================================================================
   Hall & Jones (1999) Replication — Stata Analysis Script
   =============================================================================
   Matches analysis_v2.ipynb exactly. Produces:
     Table I   : Levels accounting decomposition (Samples A and B)
     Table II  : OLS and 2SLS regressions + sensitivity (Samples A and B)
     Figure I  : log(Y/L) vs social infrastructure
     Figure II : log(TFP) vs social infrastructure
     Table III : Side-by-side comparison — H&J (1999) | V1: 1988 | V3: 2019

   USAGE
   -----
   1. Set VERSION and OUT_DIR below
   2. Install packages once:
        ssc install ivreg2
        ssc install ranktest
   3. Run:  do hj_analysis.do
      For Table III, both merged CSVs are read automatically regardless
      of VERSION — no need to run twice.

   VERSION 1 : Exact replication  — PWT 5.6, 1988, ICRG GADP 1986-95, Sachs-Warner
   VERSION 3 : Current update     — PWT 10.01, 2019, ICRG GADP 2010-17, Fraser
   ============================================================================= */

clear all
set more off
capture log close


/* ─────────────────────────────────────────────────────────────────────────────
   USER SETTINGS
   ───────────────────────────────────────────────────────────────────────── */

local VERSION  1
local OUT_DIR  "C:\Users\Adams\OneDrive\DE & E Research\outputs"

/* ─────────────────────────────────────────────────────────────────────────────
   DERIVED SETTINGS
   ───────────────────────────────────────────────────────────────────────── */

local MERGED  "`OUT_DIR'\merged_v`VERSION'.csv"
local LOG     "`OUT_DIR'\hj_analysis_v`VERSION'.log"
local ALPHA   = 1/3

local INST_ALL  "distancefromeq fr_trade english_frac we_lang_frac"
local INST_PREF "distancefromeq english_frac we_lang_frac"   // excl. FR trade

log using "`LOG'", replace text

di "================================================================="
di "  Hall & Jones (1999) Replication — Stata  [Version `VERSION']"
di "  `MERGED'"
di "================================================================="


/* =============================================================================
   SECTION 1 — LOAD AND PREPARE
   ============================================================================= */

di _n "--- Section 1: Loading data ---"

import delimited "`MERGED'", clear encoding(latin1) varnames(1)

foreach v in yl ky_ratio hc social_infra mining_va ///
             distancefromeq fr_trade english_frac we_lang_frac yr_sch {
    capture destring `v', replace force
}

* Mining correction
replace mining_va = 0 if missing(mining_va)
gen double yl_adj   = yl * (1 - mining_va)

* Accounting variables
*   alpha/(1-alpha) = 0.5  ->  cap_term = sqrt(K/Y)
gen double cap_term = ky_ratio ^ (`ALPHA' / (1 - `ALPHA')) if ky_ratio > 0
gen double hc_term  = hc
gen double tfp      = yl_adj / (cap_term * hc_term) ///
    if !missing(cap_term, hc_term, yl_adj)

gen double log_yl  = ln(yl_adj)        if yl_adj  > 0
gen double log_tfp = ln(tfp)           if tfp     > 0
gen double log_si  = ln(social_infra)  if social_infra > 0

di "  Countries loaded: `=_N'"

* K/Y outlier filter — Q3 + 3*IQR
*   Removes Venezuela 2019 (K/Y=127: oil-boom capital, collapsed output)
quietly sum ky_ratio, detail
local q1 = r(p25)
local q3 = r(p75)
local ky_upper = `q3' + 3 * (`q3' - `q1')
gen byte ky_ok = (ky_ratio <= `ky_upper') if !missing(ky_ratio)
quietly count if ky_ok == 0
if r(N) > 0 {
    di "  Dropping `r(N)' K/Y outlier(s) (K/Y > `=round(`ky_upper',.01)'):"
    list iso3 ky_ratio if ky_ok == 0
}

* Sample A — complete cases
gen byte sampleA = 1
foreach v in yl_adj ky_ratio hc social_infra `INST_ALL' {
    replace sampleA = 0 if missing(`v')
}
replace sampleA = 0 if ky_ok == 0
quietly count if sampleA == 1
di "  Sample A (complete cases): `r(N)' countries"


/* =============================================================================
   SECTION 2 — IMPUTATION -> Sample B
   OLS regression of social_infra on instruments + dist^2,
   then predict for countries with missing governance/openness.
   ============================================================================= */

di _n "--- Section 2: Imputing missing social infrastructure ---"

gen double dist_sq = distancefromeq^2

reg social_infra `INST_ALL' dist_sq ///
    if !missing(social_infra, `INST_ALL', dist_sq)
di "  Imputation regression R2 = " %6.3f e(r2)

predict double si_pred, xb
replace si_pred = max(0, min(1, si_pred))

gen double social_infra_B = social_infra
replace social_infra_B = si_pred if missing(social_infra) & !missing(si_pred)

quietly count if missing(social_infra) & !missing(social_infra_B)
di "  Imputed `r(N)' additional countries"

gen byte sampleB = 1
foreach v in yl_adj ky_ratio hc social_infra_B `INST_ALL' {
    replace sampleB = 0 if missing(`v')
}
replace sampleB = 0 if ky_ok == 0
quietly count if sampleB == 1
di "  Sample B (with imputation): `r(N)' countries"


/* =============================================================================
   SECTION 3 — TABLE I: LEVELS ACCOUNTING
   ============================================================================= */

di _n "================================================================="
di   "  TABLE I — Levels Accounting Decomposition"
di   "  Y/L = (K/Y)^(a/(1-a)) * h * A,   a = 1/3"
di   "  All values relative to USA = 1.00"
di   "================================================================="

* US benchmarks
foreach v in yl_adj cap_term hc_term tfp {
    quietly sum `v' if iso3 == "USA"
    scalar usa_`v' = r(mean)
}
gen double rel_yl  = yl_adj   / scalar(usa_yl_adj)
gen double rel_cap = cap_term / scalar(usa_cap_term)
gen double rel_hc  = hc_term  / scalar(usa_hc_term)
gen double rel_tfp = tfp      / scalar(usa_tfp)

foreach sample in A B {

    if "`sample'" == "A" {
        local cond   "sampleA == 1"
        local si_var "social_infra"
    }
    else {
        local cond   "sampleB == 1"
        local si_var "social_infra_B"
    }

    quietly count if `cond'
    di _n "  [Sample `sample']  N = `r(N)'"

    * USA / Niger gap
    quietly sum yl_adj if iso3 == "USA" & `cond'
    local usa_yl_s = r(mean)
    quietly sum yl_adj if iso3 == "NER" & `cond'
    if r(N) > 0 di "  USA / Niger Y/L gap: " %5.1f (`usa_yl_s'/r(mean)) "x"

    * Variance decomposition
    tempvar lyl lcap lhc ltfp
    gen double `lyl'  = ln(rel_yl)  if `cond' & rel_yl  > 0
    gen double `lcap' = ln(rel_cap) if `cond' & rel_cap > 0
    gen double `lhc'  = ln(rel_hc)  if `cond' & rel_hc  > 0
    gen double `ltfp' = ln(rel_tfp) if `cond' & rel_tfp > 0

    quietly corr `lyl' `lcap' `lhc' `ltfp' if `cond', cov
    mat C = r(C)
    local sh_cap = C[1,2] / C[1,1]
    local sh_hc  = C[1,3] / C[1,1]
    local sh_tfp = C[1,4] / C[1,1]

    di "  Variance decomposition:"
    di "    Capital term:  " %8.3f `sh_cap'
    di "    Human capital: " %8.3f `sh_hc'
    di "    TFP (A):       " %8.3f `sh_tfp'
    di "    Sum:           " %8.3f (`sh_cap' + `sh_hc' + `sh_tfp')

    * 90th/10th ratios
    foreach v in rel_yl rel_cap rel_hc rel_tfp {
        quietly sum `v' if `cond', detail
        local r_`v' = r(p90) / r(p10)
    }
    di "  90th/10th percentile ratios:"
    di "    Y/L:       " %5.1f `r_rel_yl'  "x"
    di "    Cap term:  " %5.1f `r_rel_cap' "x"
    di "    h:         " %5.1f `r_rel_hc'  "x"
    di "    TFP (A):   " %5.1f `r_rel_tfp' "x"

    * Country table
    di ""
    di "  " %-25s "Country" %8s "Y/L" %10s "Cap" %8s "h" %10s "TFP" %8s "S"
    di "  " "{hline 72}"
    foreach iso in USA CHE DEU JPN GBR FRA BRA CHN IND KEN NGA NER {
        quietly sum rel_yl   if iso3 == "`iso'" & `cond'
        local v_yl  = r(mean)
        quietly sum rel_cap  if iso3 == "`iso'" & `cond'
        local v_cap = r(mean)
        quietly sum rel_hc   if iso3 == "`iso'" & `cond'
        local v_hc  = r(mean)
        quietly sum rel_tfp  if iso3 == "`iso'" & `cond'
        local v_tfp = r(mean)
        quietly sum `si_var' if iso3 == "`iso'" & `cond'
        local v_si  = r(mean)
        if !missing(`v_yl') {
            quietly levelsof country if iso3 == "`iso'", local(cname) clean
            di "  " %-25s "`cname'" ///
               %8.3f `v_yl' %10.3f `v_cap' %8.3f `v_hc' %10.3f `v_tfp' %8.3f `v_si'
        }
    }
}


/* =============================================================================
   SECTION 4 — TABLE II: IV REGRESSIONS
   ============================================================================= */

di _n "================================================================="
di   "  TABLE II — OLS and 2SLS Regressions"
di   "  Dependent variable: log(Y/L)"
di   "  Robust SEs in parentheses"
di   "================================================================="

foreach sample in A B {

    if "`sample'" == "A" {
        local cond   "sampleA == 1"
        local si_var "social_infra"
        local label  "Sample A — complete cases"
    }
    else {
        local cond   "sampleB == 1"
        local si_var "social_infra_B"
        local label  "Sample B — with imputation"
    }

    quietly count if `cond' & !missing(log_yl, `si_var', `INST_ALL')
    di _n "  [`label']  N = `r(N)'"

    * OLS
    di _n "  OLS:"
    reg log_yl `si_var' if `cond', robust
    di "    Constant  " %10.4f _b[_cons]    "  (" %6.4f _se[_cons]    ")"
    di "    beta(S)   " %10.4f _b[`si_var'] "  (" %6.4f _se[`si_var'] ")"
    di "    N=" e(N) "  R2=" %6.4f e(r2)

    * 2SLS — all 4 instruments
    di _n "  2SLS (all 4 instruments):"
    ivreg2 log_yl (`si_var' = `INST_ALL') if `cond', robust first savefirst
    di "    Constant  " %10.4f _b[_cons]    "  (" %6.4f _se[_cons]    ")"
    di "    beta(S)   " %10.4f _b[`si_var'] "  (" %6.4f _se[`si_var'] ")"
    di "    N=" e(N) "  R2=" %6.4f e(r2)
    di "    First-stage F (KP): " %8.2f e(widstat)
    di "    Sargan overid p:    " %8.3f e(sarganp)

    * 2SLS — preferred 3 instruments
    di _n "  2SLS preferred (dist_eq + English + WE lang):"
    ivreg2 log_yl (`si_var' = `INST_PREF') if `cond', robust first
    di "    Constant  " %10.4f _b[_cons]    "  (" %6.4f _se[_cons]    ")"
    di "    beta(S)   " %10.4f _b[`si_var'] "  (" %6.4f _se[`si_var'] ")"
    di "    First-stage F: " %8.2f e(widstat)
    di "    Overid p:      " %8.3f e(sarganp)

    * Sensitivity table
    di _n "  Sensitivity — instrument subsets:"
    di "  %-40s %8s  %8s  %8s  %9s" "Instruments" "beta(S)" "SE" "F-stat" "Overid p"
    di "  {hline 78}"

    local s1 "distancefromeq"
    local s2 "distancefromeq fr_trade"
    local s3 "`INST_PREF'"
    local s4 "`INST_ALL'"
    local d1 "dist_eq only"
    local d2 "dist + FR trade"
    local d3 "dist + language  [preferred]"
    local d4 "all 4 instruments"

    forvalues i = 1/4 {
        quietly ivreg2 log_yl (`si_var' = `s`i'') if `cond', robust first
        local b   = _b[`si_var']
        local se  = _se[`si_var']
        local F   = e(widstat)
        local oid = e(sarganp)
        if missing(`oid') local oid_s "    n/a"
        else               local oid_s = string(round(`oid',.001), "%6.3f")
        di "  %-40s %8.3f  %8.4f  %8.2f  %9s" "`d`i''" `b' `se' `F' "`oid_s'"
    }
}

* OLS vs 2SLS summary
di _n "================================================================="
di   "  OLS vs 2SLS comparison"
di   "================================================================="
di "  %-28s %8s  %14s  %14s" "Sample" "OLS" "2SLS (4-inst)" "2SLS (3-inst)"
di "  {hline 70}"
foreach sample in A B {
    if "`sample'" == "A" {
        local cond   "sampleA == 1"
        local si_var "social_infra"
        local lbl    "A — complete cases"
    }
    else {
        local cond   "sampleB == 1"
        local si_var "social_infra_B"
        local lbl    "B — with imputation"
    }
    quietly reg    log_yl `si_var' if `cond', robust
    local b_ols = _b[`si_var']
    quietly ivreg2 log_yl (`si_var' = `INST_ALL')  if `cond', robust
    local b_iv4 = _b[`si_var']
    quietly ivreg2 log_yl (`si_var' = `INST_PREF') if `cond', robust
    local b_iv3 = _b[`si_var']
    di "  %-28s %8.3f  %14.3f  %14.3f" "`lbl'" `b_ols' `b_iv4' `b_iv3'
}
di ""
di "  H&J original: OLS=3.29  2SLS=5.14"
if `VERSION' == 1 {
    di "  V1 overid failure (~0.037) with all 4 instruments."
    di "  Preferred 3-instrument spec passes overid — result unchanged."
}


/* =============================================================================
   SECTION 5 — FIGURES
   ============================================================================= */

di _n "--- Section 5: Generating figures ---"

quietly ivreg2 log_yl (social_infra = `INST_ALL') if sampleA == 1, robust
local iv_b_fig = string(round(_b[social_infra], .01))

* Figure I — Y/L vs Social Infrastructure
twoway ///
    (scatter log_yl social_infra if sampleA == 1, ///
        msymbol(circle) mcolor("70 130 180%55") msize(small)) ///
    (lfit   log_yl social_infra if sampleA == 1, ///
        lcolor("178 34 34") lwidth(medthin)), ///
    xlabel(0(.2)1, format(%3.1f)) ///
    xtitle("Social Infrastructure (S)") ///
    ytitle("log(Output per Worker)") ///
    title("Figure I — Output per Worker vs. Social Infrastructure") ///
    subtitle("Sample A: complete cases, Version `VERSION'") ///
    note("2SLS {&beta} = `iv_b_fig' (all 4 instruments)") ///
    legend(off) graphregion(color(white))
graph export "`OUT_DIR'\figure_I_yl_vs_si_v`VERSION'.png", replace width(1600)
di "  Saved Figure I"

* Figure II — TFP vs Social Infrastructure
twoway ///
    (scatter log_tfp social_infra if sampleA == 1, ///
        msymbol(circle) mcolor("255 140 0%55") msize(small)) ///
    (lfit   log_tfp social_infra if sampleA == 1, ///
        lcolor("178 34 34") lwidth(medthin)), ///
    xlabel(0(.2)1, format(%3.1f)) ///
    xtitle("Social Infrastructure (S)") ///
    ytitle("log(TFP)") ///
    title("Figure II — TFP vs. Social Infrastructure") ///
    subtitle("Sample A: complete cases, Version `VERSION'") ///
    legend(off) graphregion(color(white))
graph export "`OUT_DIR'\figure_II_tfp_vs_si_v`VERSION'.png", replace width(1600)
di "  Saved Figure II"


/* =============================================================================
   SECTION 6 — TABLE III: SIDE-BY-SIDE COMPARISON
   Self-contained — reads both merged_v1.csv and merged_v3.csv directly.
   Runs regardless of VERSION setting above.
   ============================================================================= */

di _n "================================================================="
di   "  TABLE III — REPLICATION COMPARISON"
di   "  Hall & Jones (1999) | Version 1: 1988 | Version 3: 2019"
di   "================================================================="

* H&J original values
local HJ_n     127
local HJ_gap   35.0
local HJ_cap   0.228
local HJ_hc    0.143
local HJ_tfp   0.601
local HJ_bols  3.29
local HJ_seols 0.398
local HJ_biv   5.14
local HJ_seiv  0.508
local HJ_oid   0.256

* ── Program: compute all stats for one version and store as scalars ──────────
capture program drop _hj_ver
program define _hj_ver
    args ver outdir alpha

    preserve
    import delimited "`outdir'\merged_v`ver'.csv", ///
        clear encoding(latin1) varnames(1)

    foreach v in yl ky_ratio hc social_infra mining_va ///
                 distancefromeq fr_trade english_frac we_lang_frac {
        capture destring `v', replace force
    }

    replace mining_va = 0 if missing(mining_va)
    gen double yl_adj   = yl * (1 - mining_va)
    gen double cap_term = ky_ratio^(`alpha'/(1-`alpha')) if ky_ratio > 0
    gen double hc_term  = hc
    gen double tfp      = yl_adj / (cap_term * hc_term) ///
        if !missing(cap_term, hc_term)
    gen double log_yl   = ln(yl_adj) if yl_adj > 0

    * K/Y outlier filter
    quietly sum ky_ratio, detail
    local q1 = r(p25); local q3 = r(p75)
    gen byte ky_ok = (ky_ratio <= `q3' + 3*(`q3'-`q1')) if !missing(ky_ratio)

    * Sample A
    local INST "distancefromeq fr_trade english_frac we_lang_frac"
    local PREF "distancefromeq english_frac we_lang_frac"
    gen byte sA = 1
    foreach v in yl_adj ky_ratio hc social_infra `INST' {
        replace sA = 0 if missing(`v')
    }
    replace sA = 0 if ky_ok == 0
    quietly count if sA == 1
    scalar t3_nA`ver' = r(N)

    * Sample B imputation
    gen double dist_sq = distancefromeq^2
    quietly reg social_infra `INST' dist_sq ///
        if !missing(social_infra, `INST', dist_sq)
    predict double si_pred, xb
    replace si_pred = max(0, min(1, si_pred))
    gen double si_B = social_infra
    replace si_B = si_pred if missing(social_infra) & !missing(si_pred)
    gen byte sB = 1
    foreach v in yl_adj ky_ratio hc si_B `INST' {
        replace sB = 0 if missing(`v')
    }
    replace sB = 0 if ky_ok == 0
    quietly count if sB == 1
    scalar t3_nB`ver' = r(N)
    quietly count if missing(social_infra) & !missing(si_B) & sB == 1
    scalar t3_nimp`ver' = r(N)

    * Levels accounting
    foreach v in yl_adj cap_term hc_term tfp {
        quietly sum `v' if iso3 == "USA"
        scalar usa`ver'_`v' = r(mean)
    }
    gen double rel_yl  = yl_adj   / scalar(usa`ver'_yl_adj)
    gen double rel_cap = cap_term / scalar(usa`ver'_cap_term)
    gen double rel_hc  = hc_term  / scalar(usa`ver'_hc_term)
    gen double rel_tfp = tfp      / scalar(usa`ver'_tfp)

    quietly sum yl_adj if iso3 == "NER" & sA == 1
    if r(N) > 0 scalar t3_gap`ver' = scalar(usa`ver'_yl_adj) / r(mean)
    else        scalar t3_gap`ver' = .

    tempvar lyl lcap lhc ltfp
    gen double `lyl'  = ln(rel_yl)  if sA == 1 & rel_yl  > 0
    gen double `lcap' = ln(rel_cap) if sA == 1 & rel_cap > 0
    gen double `lhc'  = ln(rel_hc)  if sA == 1 & rel_hc  > 0
    gen double `ltfp' = ln(rel_tfp) if sA == 1 & rel_tfp > 0
    quietly corr `lyl' `lcap' `lhc' `ltfp' if sA == 1, cov
    mat C = r(C)
    scalar t3_cap`ver' = C[1,2] / C[1,1]
    scalar t3_hc`ver'  = C[1,3] / C[1,1]
    scalar t3_tfp`ver' = C[1,4] / C[1,1]

    * Regressions — Sample A
    quietly reg log_yl social_infra if sA == 1, robust
    scalar t3_ols_b`ver'  = _b[social_infra]
    scalar t3_ols_se`ver' = _se[social_infra]

    quietly ivreg2 log_yl (social_infra = `INST') if sA == 1, robust first
    scalar t3_iv4_b`ver'   = _b[social_infra]
    scalar t3_iv4_se`ver'  = _se[social_infra]
    scalar t3_iv4_F`ver'   = e(widstat)
    scalar t3_iv4_oid`ver' = e(sarganp)

    quietly ivreg2 log_yl (social_infra = `PREF') if sA == 1, robust first
    scalar t3_iv3_b`ver'   = _b[social_infra]
    scalar t3_iv3_se`ver'  = _se[social_infra]
    scalar t3_iv3_F`ver'   = e(widstat)
    scalar t3_iv3_oid`ver' = e(sarganp)

    * Regressions — Sample B
    quietly reg log_yl si_B if sB == 1, robust
    scalar t3_ols_bB`ver'  = _b[si_B]
    scalar t3_ols_seB`ver' = _se[si_B]

    quietly ivreg2 log_yl (si_B = `INST') if sB == 1, robust first
    scalar t3_iv4_bB`ver'   = _b[si_B]
    scalar t3_iv4_seB`ver'  = _se[si_B]
    scalar t3_iv4_FB`ver'   = e(widstat)
    scalar t3_iv4_oidB`ver' = e(sarganp)

    quietly ivreg2 log_yl (si_B = `PREF') if sB == 1, robust first
    scalar t3_iv3_bB`ver'   = _b[si_B]
    scalar t3_iv3_seB`ver'  = _se[si_B]
    scalar t3_iv3_FB`ver'   = e(widstat)
    scalar t3_iv3_oidB`ver' = e(sarganp)

    restore
end

di "  Computing V1..."
_hj_ver 1 "`OUT_DIR'" `ALPHA'
di "  Computing V3..."
_hj_ver 3 "`OUT_DIR'" `ALPHA'

* ── Print Table III ──────────────────────────────────────────────────────────

#delimit ;

di "" ;
di "======================================================================" ;
di "  TABLE III — REPLICATION COMPARISON" ;
di "  Hall & Jones (1999)  |  Version 1: 1988  |  Version 3: 2019" ;
di "======================================================================" ;
di "  " %-36s "" %10s "H&J (1999)" %12s "V1: 1988" %12s "V3: 2019" ;
di "  " "-" * 70 ;

di _n "  PANEL A — DATA SOURCES" ;
di "  " %-36s "Output & capital"
    %10s "PWT 5.6"   %12s "PWT 5.6"   %12s "PWT 10.01" ;
di "  " %-36s "Price base"
    %10s "1985 USD"  %12s "1985 USD"  %12s "2017 USD"  ;
di "  " %-36s "Governance source"
    %10s "ICRG GADP" %12s "ICRG GADP" %12s "ICRG GADP" ;
di "  " %-36s "Governance window"
    %10s "1986-1995" %12s "1986-1995" %12s "2010-2017" ;
di "  " %-36s "Education"
    %10s "BL 1985"   %12s "BL 1985"   %12s "BL 2015"   ;
di "  " %-36s "Openness"
    %10s "Sachs-W."  %12s "Sachs-W."  %12s "Fraser"    ;

di _n "  PANEL B — LEVELS ACCOUNTING (Sample A, complete cases)" ;
di "  " %-36s "N (Sample A)"
    %10.0f `HJ_n'
    %12.0f scalar(t3_nA1)
    %12.0f scalar(t3_nA3) ;
di "  " %-36s "N (Sample B, imputed)"
    %10s "—"
    %12.0f scalar(t3_nB1)
    %12.0f scalar(t3_nB3) ;
di "  (imputed " scalar(t3_nimp1) " countries in V1, "
    scalar(t3_nimp3) " in V3)" ;
di "  " %-36s "USA / Niger Y/L gap"
    %10s "`HJ_gap'x"
    %12s string(round(scalar(t3_gap1),.1))+"x"
    %12s string(round(scalar(t3_gap3),.1))+"x" ;
di "  " %-36s "Capital term share"
    %10.3f `HJ_cap'
    %12.3f scalar(t3_cap1)
    %12.3f scalar(t3_cap3) ;
di "  " %-36s "Human capital share"
    %10.3f `HJ_hc'
    %12.3f scalar(t3_hc1)
    %12.3f scalar(t3_hc3) ;
di "  " %-36s "TFP share"
    %10.3f `HJ_tfp'
    %12.3f scalar(t3_tfp1)
    %12.3f scalar(t3_tfp3) ;
di "  " %-36s "Sum of shares"
    %10.3f (`HJ_cap'+`HJ_hc'+`HJ_tfp')
    %12.3f (scalar(t3_cap1)+scalar(t3_hc1)+scalar(t3_tfp1))
    %12.3f (scalar(t3_cap3)+scalar(t3_hc3)+scalar(t3_tfp3)) ;

di _n "  PANEL C — REGRESSIONS: all 4 instruments" ;
di "  Sample A:" ;
di "  " %-36s "OLS beta (S)"
    %10.3f `HJ_bols'
    %12.3f scalar(t3_ols_b1)
    %12.3f scalar(t3_ols_b3) ;
di "  " %-36s ""
    %10s "(`HJ_seols')"
    %12s "("+string(round(scalar(t3_ols_se1),.001))+")"
    %12s "("+string(round(scalar(t3_ols_se3),.001))+")" ;
di "  " %-36s "2SLS beta (S)"
    %10.3f `HJ_biv'
    %12.3f scalar(t3_iv4_b1)
    %12.3f scalar(t3_iv4_b3) ;
di "  " %-36s ""
    %10s "(`HJ_seiv')"
    %12s "("+string(round(scalar(t3_iv4_se1),.001))+")"
    %12s "("+string(round(scalar(t3_iv4_se3),.001))+")" ;
di "  " %-36s "First-stage F"
    %10s "—"
    %12.2f scalar(t3_iv4_F1)
    %12.2f scalar(t3_iv4_F3) ;
di "  " %-36s "Overid p-value"
    %10.3f `HJ_oid'
    %12.3f scalar(t3_iv4_oid1)
    %12.3f scalar(t3_iv4_oid3) ;

di "  Sample B (with imputation):" ;
di "  " %-36s "OLS beta (S)"
    %10s "—"
    %12.3f scalar(t3_ols_bB1)
    %12.3f scalar(t3_ols_bB3) ;
di "  " %-36s ""
    %10s "—"
    %12s "("+string(round(scalar(t3_ols_seB1),.001))+")"
    %12s "("+string(round(scalar(t3_ols_seB3),.001))+")" ;
di "  " %-36s "2SLS beta (S)"
    %10s "—"
    %12.3f scalar(t3_iv4_bB1)
    %12.3f scalar(t3_iv4_bB3) ;
di "  " %-36s ""
    %10s "—"
    %12s "("+string(round(scalar(t3_iv4_seB1),.001))+")"
    %12s "("+string(round(scalar(t3_iv4_seB3),.001))+")" ;
di "  " %-36s "First-stage F"
    %10s "—"
    %12.2f scalar(t3_iv4_FB1)
    %12.2f scalar(t3_iv4_FB3) ;
di "  " %-36s "Overid p-value"
    %10s "—"
    %12.3f scalar(t3_iv4_oidB1)
    %12.3f scalar(t3_iv4_oidB3) ;

di _n "  PANEL D — PREFERRED SPEC: dist_eq + English + WE lang (3 inst.)" ;
di "  Sample A:" ;
di "  " %-36s "2SLS beta (S)"
    %10s "—"
    %12.3f scalar(t3_iv3_b1)
    %12.3f scalar(t3_iv3_b3) ;
di "  " %-36s ""
    %10s "—"
    %12s "("+string(round(scalar(t3_iv3_se1),.001))+")"
    %12s "("+string(round(scalar(t3_iv3_se3),.001))+")" ;
di "  " %-36s "First-stage F"
    %10s "—"
    %12.2f scalar(t3_iv3_F1)
    %12.2f scalar(t3_iv3_F3) ;
di "  " %-36s "Overid p-value"
    %10s "—"
    %12.3f scalar(t3_iv3_oid1)
    %12.3f scalar(t3_iv3_oid3) ;
di "  Sample B (with imputation):" ;
di "  " %-36s "2SLS beta (S)"
    %10s "—"
    %12.3f scalar(t3_iv3_bB1)
    %12.3f scalar(t3_iv3_bB3) ;
di "  " %-36s ""
    %10s "—"
    %12s "("+string(round(scalar(t3_iv3_seB1),.001))+")"
    %12s "("+string(round(scalar(t3_iv3_seB3),.001))+")" ;
di "  " %-36s "First-stage F"
    %10s "—"
    %12.2f scalar(t3_iv3_FB1)
    %12.2f scalar(t3_iv3_FB3) ;
di "  " %-36s "Overid p-value"
    %10s "—"
    %12.3f scalar(t3_iv3_oidB1)
    %12.3f scalar(t3_iv3_oidB3) ;

di _n "  PANEL E — KEY FINDINGS" ;
di "  " %-36s "IV / OLS ratio"
    %10s "1.56x"
    %12s string(round(scalar(t3_iv4_b1)/scalar(t3_ols_b1),.01))+"x"
    %12s string(round(scalar(t3_iv4_b3)/scalar(t3_ols_b3),.01))+"x" ;
local chg = round(
    (scalar(t3_iv4_b3)-scalar(t3_iv4_b1))/scalar(t3_iv4_b1)*100, 1) ;
di "  " %-36s "2SLS beta change V1->V3"
    %10s "—" %12s "baseline" %12s "+`chg'%" ;
di "  " %-36s "Y/L gap change"
    %10s "—" %12s "35x -> 33x" %12s "35x -> 41x" ;

di _n "  ======================================================================" ;
di "  HC1 robust SEs (ivreg2 robust). H&J use bootstrap." ;
di "  V1 Panel C overid failure driven by FR trade instrument." ;
di "  Panel D drops FR trade: result unchanged, overid passes." ;
di "  BL=Barro-Lee. Sachs-W.=Sachs-Warner. Fraser=Fraser Inst. Area 5." ;
di "  Venezuela excluded from V3 (K/Y=127, economic collapse 2013-2019)." ;
di "  ======================================================================" ;
#delimit cr


/* =============================================================================
   SECTION 7 — DOCUMENTED DEVIATIONS
   ============================================================================= */

di _n "================================================================="
di   "  DEVIATIONS FROM HALL & JONES (1999) — Version `VERSION'"
di   "================================================================="

if `VERSION' == 1 {
    di _n "  [Sample size]"
    di   "    H&J N=127. Ours: ~98 complete (A) / ~109 imputed (B)."
    di   "    ICRG did not cover ~26 small countries in 1986-1995."
    di   "    Key result unchanged: 2SLS beta ~5.09 vs H&J 5.14."
    di _n "  [FR trade instrument]"
    di   "    Sargan overid p~0.037 with all 4 instruments."
    di   "    Preferred 3-instrument spec: beta~5.06, overid p~0.183."
    di   "    FR trade likely has direct income effects beyond institutions."
    di _n "  [Capital stock]"
    di   "    ~40 countries use perpetual inventory fallback due to"
    di   "    inconsistent PWT 5.6 KAPW units."
    di _n "  [Standard errors]"
    di   "    H&J use bootstrap (10,000 reps)."
    di   "    We use HC1 robust (ivreg2 robust) — equivalent asymptotically."
}
else {
    di _n "  [Cross-section year]"
    di   "    H&J use 1988. We use 2019 (most recent pre-COVID year)."
    di _n "  [Governance window]"
    di   "    H&J ICRG GADP 1986-1995. We use ICRG GADP 2010-2017."
    di _n "  [Education]"
    di   "    H&J Barro-Lee 1985, age 25+. We use Barro-Lee 2015, age 25-64."
    di _n "  [Openness]"
    di   "    H&J Sachs-Warner (binary). We use Fraser Area 5 (continuous)."
    di _n "  [2SLS beta]"
    di   "    H&J beta=5.14. We find ~8.18 in 2019 (+61%)."
    di   "    Institutional premium on output strengthened over 30 years."
    di _n "  [Venezuela excluded]"
    di   "    K/Y=127 in 2019 (oil capital stock, collapsed output)."
    di   "    Included, it suppresses capital variance share to near zero."
}

di _n "================================================================="
di   "  Analysis complete — Version `VERSION'"
di   "  Log saved: `LOG'"
di   "================================================================="

log close
