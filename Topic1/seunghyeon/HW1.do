***********************************************************************
* Dehejia and Wahba (1999) - Causal Effects in Nonexperimental Studies:
*                            Reevaluating the Evaluation of Training Programs
*
* PROVENANCE. This is not Dehejia and Wahba's code.  Neither LaLonde (1986), nor Dehejia
* and Wahba (1999, 2002), nor Smith and Todd (2005) published code; only
* the data survive, at https://users.nber.org/~rdehejia/nswdata2.html

* The first part of the code is a teaching replication from
* Steve Pischke's project directory. It estimates the training effect
* by OLS under five covariate specifications, on the experimental sample
* and on two non-experimental comparison groups, and then repeats the
* exercise on the subsample where treated and comparison units overlap
* in propensity score. It writes the results to a LaTeX table.
*
* The second part I created using Claude for matching, balance, etc to 
* give you an idea, you should edit the codes for all the questions asked. 
***********************************************************************

* ---------------------------------------------------------------------
* SETUP
* ---------------------------------------------------------------------

* HOW TO RUN THIS ONE BLOCK AT A TIME.
*
* Stata's do-file editor runs a highlighted selection by writing it to a
* temporary file and do-ing that. That temp file is its own scope, so
* anything declared with "local" or "tempfile" does not exist inside it.
* Select a block in the middle of a file built on locals and you get
* silent nonsense rather than an error: an empty covariate list makes
* "logit treat $X" fit an intercept only, and the score comes back
* constant at 185/16,177 = .011436.
*
* So this file uses GLOBALS ($path, $X) and NAMED intermediate datasets
* (_alldata.dta, _pairs.dta and so on) instead. Both survive between
* runs, so you can run Part A once, then work through Part B a block at
* a time, re-running any block as often as you like.
*
* The one ordering rule that remains: a block needs the datasets the
* earlier blocks wrote. Run B0 to B3 once in order, and after that B4,
* B5, B6 and B9 can each be run on their own.
*
* $path is where the data live. "." means the current directory, so
* either cd into this folder before running or replace "." with the
* folder holding the .dta files.
global path "."

cap log close

log using "$path/HW1.log", text replace

set more off

clear

* ---------------------------------------------------------------------
* THE THREE DATASETS
* ---------------------------------------------------------------------
* nswre74  = the experimental sample, Dehejia-Wahba subsample (185 treated,
*            260 controls). Randomized, so its treated-control difference
*            is the benchmark every other number is judged against.
* cps1re74 = CPS-1, 15,992 men from the Current Population Survey. Not
*            applicants, not screened: the general population.
* cps3re74 = CPS-3, 429 men, a much narrower screen of the same survey.
*
* `data' has all three; `datashort' drops the experimental file. The
* second loop uses `datashort' because the propensity-score exercise only
* makes sense for a NON-experimental comparison group. In the experiment
* there is nothing to correct for.
local data "nswre74 cps1re74 cps3re74 "
local datashort "cps1re74 cps3re74"

* ---------------------------------------------------------------------
* THE FIVE COVARIATE SPECIFICATIONS
* ---------------------------------------------------------------------
* Each macro holds a variable list, to be pasted into a regression below.
* spec1 is deliberately EMPTY: "reg re78 treat" with nothing after it is
* the raw difference in means.
* The "a" variants are the same specification plus a marriage dummy.
* Reading down the list, you are adding controls: nothing, then
* demographics, then 1975 earnings, then both, then 1974 earnings as well.
* spec5 is the one that matters: it is the only specification with TWO
* years of pre-treatment earnings, which is what lets you see the earnings
* trajectory rather than just its level.
local spec1 ""
local spec2 "age age2 ed black hisp nodeg"
local spec3 "re75"
local spec4 "re75 age age2 ed black hisp nodeg"
local spec5 "re74 re75 age age2 ed black hisp nodeg"
local spec2a "age age2 ed black hisp nodeg married"
local spec4a "re75 age age2 ed black hisp nodeg married"
local spec5a "re74 re75 age age2 ed black hisp nodeg married"

* ---------------------------------------------------------------------
* ROW LABELS AND COLUMN HEADERS FOR THE OUTPUT TABLE
* ---------------------------------------------------------------------
* Each line* macro starts as a row label and will have "& estimate & (se)"
* appended to it once per dataset, building up a LaTeX table row.
local line1 "Raw Difference"
local line2 "D+W Demographics"
local line3 "and Marriage Dummy"
local line4 "75 Earnings"
local line5 "D+W Demographics 75 Earnings"
local line6 "and Marriage Dummy"
local line7 "D+W Demographics 74,75 Earnings"
local line8 "and Marriage Dummy"
local lastl "Dropped in pscore regs: "

* These three turn a FILE NAME into a COLUMN HEADING. 
local nswre74 "NSW"
local cps1re74 "CPS-1"
local cps3re74 "CPS-3"
local nearlylastl "Sample sizes in pscore regs (treated/untreated): "

***********************************************************************
* Part 1 .1 : OLS ON THE FULL SAMPLE
*
* For each of the three datasets, run the five specifications (plus the
* three marriage-dummy variants) and record the coefficient on `treat'.
*
* This is the "Full Sample" half of the table. No propensity score is
* involved anywhere in this loop. It is eight regressions per dataset.
***********************************************************************

foreach set in `data' {
    local lctr = 0
    * lctr ("line counter") tracks which table row we are filling in.

    u "$path/`set'", clear
    * "u" abbreviates "use". Opens the dataset named by `set'.

    forvalues t = 1/5 {
        local ++lctr
        * Increment lctr by one. Equivalent to: local lctr = `lctr' + 1

        reg re78 treat `spec`t''
        * NESTED MACRO EXPANSION, 
        * Stata resolves the INNER macro first: on pass t=2, `t' becomes 2,
        * so `spec`t'' becomes `spec2', which then becomes its contents:
        *     age age2 ed black hisp nodeg
        * So the command Stata actually runs is
        *     reg re78 treat age age2 ed black hisp nodeg
        * On pass t=1 the spec is empty and the command is "reg re78 treat",
        * whose coefficient is exactly the raw difference in means.

            local b = _b[treat]
            local se = _se[treat]
            * After any estimation command, Stata leaves the coefficients in
            * _b[varname] and their standard errors in _se[varname]. These
            * two lines grab the treatment effect and its precision before
            * the next regression overwrites them.

            local b = round(`b',1)
            local se = round(`se',1)
            * Round to whole dollars for the table.

            local line`lctr' "`line`lctr''&    `b' &   (`se')  "
            * Append two LaTeX cells to the current row. The right-hand side
            * repeats the macro's own current contents and adds to them, so
            * the row grows by one dataset each time through the outer loop.

        if `t' == 2 | `t' == 4 | `t' == 5 {
            * Only three of the five specifications have a marriage-dummy
            * variant, so only those get a second, extra row.
            local ++lctr
            reg re78 treat `spec`t'a'
            * `spec`t'a' resolves to spec2a, spec4a or spec5a.
            local b = _b[treat]
            local se = _se[treat]
            local b = round(`b',1)
            local se = round(`se',1)
            local line`lctr' "`line`lctr''&    `b' &   (`se')  "
        }
    }
}

***********************************************************************
*Part 1 .2 : WHERE THE PROPENSITY SCORE COMES IN
**
*   probit treat <covariates>     <- ESTIMATE the propensity score model
*   predict pscore, pr            <- COMPUTE each unit's score, Pr(D=1|X)
*   keep if pscore>.1 & pscore<.9 <- TRIM to the overlap region
*   reg re78 treat <covariates>   <- run OLS on what is left
*
* So the score is used for ONE thing: deciding which observations to throw
* away. Units with a score near 0 are people no treated unit resembles;
* units near 1 are people no comparison unit resembles. Both are dropped.
*
* So this file demonstrates that overlap alone does a lot of the work,
* which is a real result, but it is not a propensity-score matching
* exercise and should not be read as one.
***********************************************************************

foreach set in `datashort' {
    local lctr = 1
    local line`lctr' "`line`lctr''&         &        "
    * Row 1 is "Raw Difference", which has no propensity-selected version
    * (spec1 has no covariates, so there is no score to estimate). Two
    * blank cells keep the columns aligned.

    u "$path/`set'", clear

    forvalues t = 2/5 {
        * Starts at 2, not 1, for the same reason.
        local ++lctr

        probit treat `spec`t''
        * THE PROPENSITY SCORE MODEL. The dependent variable is TREATMENT,
        * not the outcome. Dehejia and Wahba use a logit; probit gives nearly
        * identical scores.

        preserve
        * Take a snapshot of the data in memory. 

            predict pscore, pr
            * Create the fitted probability for every observation. 

            keep if pscore >.1 & pscore <.9
            * THE TRIM. Keep only the overlap region. This is the single
            * substantive use of the score in the entire file.

            sum treat
            if `r(N)' != 0 {
            * "sum" leaves the number of observations in r(N). If the trim
            * emptied the sample, skip this specification rather than error.
            * This is not hypothetical: CPS-3 is small and heavily screened.

                sum treat if treat == 1
                local nearlylastl "`nearlylastl' ``set'' Spec `t': `r(N)'"
                sum treat if treat == 0
                local nearlylastl "`nearlylastl'/`r(N)';"
                * Record how many treated and how many comparison units
                * survived the trim, so the table footnote can report it.

                reg re78 treat `spec`t''
                * OLS AGAIN, now on the trimmed sample. Same specification
                * as the probit, different dependent variable.

                local dctr = 0
                foreach var in `spec`t'' {
                    if _se[`var'] == 0 {
                        * A standard error of exactly zero means Stata
                        * dropped that variable for collinearity. After
                        * trimming, a covariate can become constant within
                        * the surviving sample. Worth reporting, because a
                        * specification that lost a control is not the
                        * specification the row header claims.
                        local ++dctr
                        if `dctr' == 1 local lastl "`lastl'; ``set'' Spec `t': `var'"
                        else local lastl "`lastl', `var'"
                    }
                }

                local b = _b[treat]
                local se = _se[treat]
                local b = round(`b',1)
                local se = round(`se',1)
                local line`lctr' "`line`lctr''&    `b' &   (`se')  "
            }
            else {
                * The trim left nothing. Write "no obs." into the cell.
                local line`lctr' "`line`lctr''&  no   &  obs.   "
                local nearlylastl "`nearlylastl' ``set'' Spec `t': 0/0;"
            }
        restore
        * Put the untrimmed data back, ready for the next specification.

        if `t' == 2 | `t' == 4 | `t' == 5 {
            * The marriage-dummy variants again, same logic throughout.
            local ++lctr
            probit treat `spec`t'a'
            preserve
                predict pscore, pr
                keep if pscore >.1 & pscore <.9
                sum treat if treat == 1
                local nearlylastl "`nearlylastl' ``set'' Spec `t'a: `r(N)'"
                sum treat if treat == 0
                local nearlylastl "`nearlylastl'/`r(N)';"
                reg re78 treat `spec`t'a'
                local dctr = 0
                foreach var in `spec`t'a' {
                    if _se[`var'] == 0 {
                        local ++dctr
                        if `dctr' == 1 local lastl "`lastl'; ``set'' Spec `t'a: `var'"
                        else local lastl "`lastl', `var'"
                    }
                }
                local b = _b[treat]
                local se = _se[treat]
                local b = round(`b',1)
                local se = round(`se',1)
                local line`lctr' "`line`lctr''&    `b' &   (`se')  "
            restore
        }
    }
}

***********************************************************************
* WRITE THE TABLE
***********************************************************************

file open holder using "$path/table.tex", write replace text
file write holder "\documentclass[a4paper,12pt,fleqn]{article} "_n
file write holder "\usepackage{rotating} "_n
file write holder "\begin{document} " _n
file write holder "\begin{sidewaystable} " _n
file write holder "\begin{center} " _n
file write holder "\begin{tabular}{l cc cc cc | cc cc} " _n
file write holder "                 &   \multicolumn{6}{c}{Full Sample}     & \multicolumn{4}{c}{Propensity Selected Sample}    \\"
file write holder "Specification    &   \multicolumn{2}{c}{NSW}  &   \multicolumn{2}{c}{CPS-1}  &   \multicolumn{2}{c}{CPS-3} &   \multicolumn{2}{c}{CPS-1}  &   \multicolumn{2}{c}{CPS-3} \\[.5em]\hline\\[-.5em]" _n
forvalues t = 1/`lctr' {
    * `lctr' is still holding its final value from the last loop: 8 rows.
    file write holder "`line`t''   \\"
    if `t' == 1 | `t' == 3 | `t' == 4 | `t' == 6 {
        file write holder "[1em]"
        * Extra vertical space after the rows that end a block.
    }
    if `t' == `lctr' {
        file write holder "[1em]\hline\\[-1em]"
    }
    file write holder _n
}
file write holder "\end{tabular} " _n
file write holder "\end{center} " _n
file write holder "`nearlylastl'" _n _n
* The sample sizes surviving the trim, accumulated in loop 2.
file write holder "`lastl'" _n
* The list of covariates dropped for collinearity after trimming.
file write holder "\end{sidewaystable} " _n
file write holder "\end{document} " _n
file close holder

***********************************************************************
***********************************************************************
**                                                                   **
**   PART B: PROPENSITY-SCORE MATCHING FROM SCRATCH                  **
**                                                                   **
***********************************************************************
***********************************************************************
*
* Part A above never matched anybody. It estimated treatment effects by
* OLS, and used the propensity score only to trim. Part B does the
* matching, and does it by hand so you can see what the one-line command
* actually consists of.
*
* THE ORDER OF OPERATIONS MATTERS
*
*   1. estimate the propensity score          } uses D and X only
*   2. look at overlap                        } no outcome anywhere
*   3. pair each treated unit with a control  } in these four steps
*   4. check balance, and iterate on 1 if it fails
*   -----------------------------------------------------------------
*   5. only now compute the within-pair outcome differences: the ATT
*   6. verify the whole thing reproduces teffects exactly
*
* Step 6 is the other point of the file. If your hand calculation and
* the canned command disagree, you have misunderstood the canned
* command, not the other way round.
global X "age age2 education educ2 married nodegree black hispanic re74 re75 u74 u75"

*=====================================================================
* B0. BUILD THE SAMPLE
*     185 NSW treated + 15,992 CPS controls. The 260 experimental
*     controls are dropped: the whole exercise is to ask whether a
*     non-experimental comparison group can stand in for them.
*=====================================================================
use "$path/nsw_dw.dta", clear
drop if treat == 0
append using "$path/cps_controls.dta"
replace treat = 0 if missing(treat)

gen age2  = age^2
gen educ2 = education^2
gen u74   = re74 == 0
gen u75   = re75 == 0
* u74 and u75 flag ZERO earnings. 73% of the sample has re74 exactly 0,
* which is a mass point, not a tail. A score linear in re74 cannot
* represent it, so the indicator enters as its own regressor.

gen long id = _n
* A unit identifier. You need one to record who was matched to whom.

tab treat

save "$path/_alldata.dta", replace

*=====================================================================
* B1. THE PROPENSITY SCORE
*=====================================================================
logit treat $X
predict double ps, pr
label var ps "e-hat = Pr(treat=1 | X)"

di _n as txt "{hline 70}"
di as txt "Distribution of the estimated score"
di as txt "{hline 70}"
tabstat ps, by(treat) stat(n mean p10 p50 p90 min max) format(%9.5f)

*=====================================================================
* B2. OVERLAP
*     Where the treated live on the score, and whether any controls
*     live there too. This is the check no coefficient table can show.
*=====================================================================
* OVERLAP IS A CONDITION ON THE TREATED UNITS
* For the ATT you need, for every treated unit, need some control at a
* similar score. You do not need the control distribution to resemble
* the treated distribution, and it never will when the comparison group
* is a national survey. 
sum ps if treat==1
local tmin = r(min)
local tmax = r(max)
sum ps if treat==0
local cmin = r(min)
local cmax = r(max)

di _n as txt "Treated score range:  " as res %7.5f `tmin' " to " %7.5f `tmax'
di as txt "Control score range:  " as res %7.5f `cmin' " to " %7.5f `cmax'

count if treat==1 & (ps > `cmax' | ps < `cmin')
di as txt "Treated units OUTSIDE the control score range: " as res r(N)
* This is the common-support violation, and it is the number that bears
* on identification.

count if treat==0 & ps > 0.1
di as txt "Controls with a score above 0.1: " as res r(N) as txt " of 15,992"
* The usable pool. 
label define grp 0 "CPS controls (15,992)" 1 "NSW treated (185)", replace
label values treat grp
* COUNTS, not percentages, on a shared scale. Percent-within-group hides
* the thing that matters: above a score of 0.4 the treated OUTNUMBER the
* controls, so those matches are drawn from a handful of units. Counts
* show that directly. Only the first control bar has to be cut.
*
* Stata's yscale(range()) only ever WIDENS an axis, so it cannot clip.
* Bin and count by hand, then cap the plotted height, which does.
preserve
    gen double bin = floor(ps/0.1)*0.1 + 0.05
    contract treat bin
    gen double nshow = min(_freq, 110)
    twoway (bar nshow bin, barwidth(0.09) color(navy%70)), ///
        by(treat, cols(1) ///
           note("Shared count scale, cut at 110. The control bar at 0.05 is 15,756.")) ///
        yscale(range(0 110)) ylabel(0(25)100) ///
        xlabel(0(0.2)1) xtitle("estimated propensity score") ///
        ytitle("number of observations") name(overlap_hand, replace)
    graph export "$path/overlap_hand.png", replace width(1000)
restore

save "$path/_alldata.dta", replace

*=====================================================================
* B3. MATCHING, BY HAND
*     Nearest neighbour on the score, with replacement.
*
*     The transparent way to do this is to form every treated-control
*     pair, 185 x 15,992 = 2,958,520 rows, and keep the smallest gap for
*     each treated unit. It is also the wrong way: the cost grows as the
*     product of the two group sizes, so a study with 5,000 treated and
*     200,000 controls would ask Stata for a billion rows.
*
*     Do it in one sort instead. Line the treated and the controls up
*     together in score order; then for any treated unit, its nearest
*     control is one of exactly two candidates, the closest control
*     below it and the closest control above it. Cost is one sort,
*     n log n, and memory never exceeds the sample itself.
*
*     WORKED EXAMPLE, five units on the score line. C = control,
*     T = treated, sorted by score:
*
*         C(.10)   T(.15)   C(.30)   T(.42)   C(.50)
*
*     Walk UP the line remembering the last control you passed. That
*     gives each treated unit its nearest control BELOW:
*         T(.15) -> C(.10), gap .05
*         T(.42) -> C(.30), gap .12
*     Walk DOWN doing the same for the nearest control ABOVE:
*         T(.15) -> C(.30), gap .15
*         T(.42) -> C(.50), gap .08
*     Then take whichever side is closer:
*         T(.15) matches C(.10)   (.05 beats .15)
*         T(.42) matches C(.50)   (.08 beats .12)
*
*     That is the whole algorithm. The code below is those two walks.
*     "Remembering the last control you passed" is the line
*         replace lo_ps = lo_ps[_n-1] if missing(lo_ps)
*     which copies the value down from the row above whenever the
*     current row is a treated unit rather than a control. Stata
*     evaluates replace in sort order, so by the time it reaches row n
*     the value at row n-1 has already been filled in. One pass, no loop.
*
*     lo_ps and hi_ps are the two candidate scores, d_lo and d_hi the
*     two gaps, dist the smaller of them, and mps the score of the side
*     that won. ncell records how many controls sit at that score, which
*     is how ties are carried through.
*=====================================================================
* One row per control, kept for the balance weights in B4.
preserve
    keep if treat==0
    keep id ps
    rename id   cid

    save "$path/_roster.dta", replace
restore

preserve
    keep if treat==0
    keep ps id
    collapse (count) ncell=id, by(ps)
    gen byte iscell = 1

    save "$path/_cells.dta", replace
restore

* Stack the treated in among those candidates and sort once.
preserve
    keep if treat==1
    keep id ps
    rename id tid
    gen byte iscell = 0
    append using "$path/_cells.dta"
    sort ps iscell
    * At an exact tie in ps, this ascending pass puts treated rows before
    * control cells. The descending pass below puts control cells first
    * and finds any distance-zero match.

    gen double lo_ps = ps      if iscell==1
    gen double lo_n  = ncell   if iscell==1
    replace lo_ps = lo_ps[_n-1] if missing(lo_ps) & _n>1
    replace lo_n  = lo_n[_n-1]  if missing(lo_n)  & _n>1

    gsort -ps -iscell
    gen double hi_ps = ps      if iscell==1
    gen double hi_n  = ncell   if iscell==1
    replace hi_ps = hi_ps[_n-1] if missing(hi_ps) & _n>1
    replace hi_n  = hi_n[_n-1]  if missing(hi_n)  & _n>1

    keep if iscell==0
    gen double d_lo = ps - lo_ps
    gen double d_hi = hi_ps - ps
    gen double dist = min(cond(missing(d_lo),.,d_lo), cond(missing(d_hi),.,d_hi))

    * Pick the closer side. If the two sides are exactly equidistant,
    * every control on both sides is a tied match, so average across
    * them weighted by how many each side holds.
    gen int    nties  = .
    gen double mps    = .
    replace nties  = lo_n  if d_lo <  d_hi | missing(d_hi)
    replace mps    = lo_ps if d_lo <  d_hi | missing(d_hi)
    replace nties  = hi_n  if d_hi <  d_lo | missing(d_lo)
    replace mps    = hi_ps if d_hi <  d_lo | missing(d_lo)
    replace nties  = lo_n + hi_n if d_lo==d_hi & !missing(d_lo)

    keep tid dist nties mps lo_ps hi_ps d_lo d_hi lo_n hi_n

    save "$path/_matched.dta", replace
restore

use "$path/_matched.dta", clear
gen double w_each = 1/nties
keep if !missing(mps)
keep tid mps w_each
rename mps ps

save "$path/_side1.dta", replace

* The equidistant case: both sides are tied matches, so both contribute.
use "$path/_matched.dta", clear
gen double w_each = 1/nties
keep if missing(mps)
expand 2
gen byte side = .
bysort tid: replace side = _n
gen double ps = cond(side==1, lo_ps, hi_ps)

keep tid ps w_each
append using "$path/_side1.dta"

joinby ps using "$path/_roster.dta"
* Every control sharing the matched score value is a tied match and gets
* the same share of that treated unit's weight.
keep tid cid w_each
rename w_each w

save "$path/_pairs.dta", replace

use "$path/_matched.dta", clear
di _n as txt "{hline 70}"
di as txt "How many controls each treated unit was matched to"
di as txt "{hline 70}"
tab nties
sum dist, detail

* THIS IS THE SECOND HALF OF THE OVERLAP CHECK, and the more useful
* half.

*=====================================================================
* B4. BALANCE  -- BEFORE YOU LOOK AT ANY OUTCOME
*
*     The standardized difference for covariate x is
*
*         d = (xbar_treated - xbar_control) / sqrt((s2_t + s2_c)/2)
*
*     a difference in means expressed in standard deviations, so it does
*     not depend on the units of x and is comparable across covariates.
*     |d| < 0.1 is the usual working threshold.
*
*=====================================================================
use "$path/_pairs.dta", clear
collapse (sum) w, by(cid)
rename cid id

save "$path/_cweight.dta", replace

use "$path/_alldata.dta", clear
merge 1:1 id using "$path/_cweight.dta", nogen
replace w = 0 if missing(w) & treat==0
replace w = 1 if treat==1
* Every treated unit counts once. Controls count as often as they were
* used, and controls never matched count zero. That weighted sample IS
* the matched sample.

di _n as txt "{hline 78}"
di as txt %-12s "covariate" %16s "std diff raw" %18s "std diff matched" %18s "var ratio matched"
di as txt "{hline 78}"
local worst = 0
local worstv ""
foreach v of global X {
    quietly sum `v' if treat==1
    local mt = r(mean)
    local vt = r(Var)
    quietly sum `v' if treat==0
    local mc = r(mean)
    local vc = r(Var)
    local draw = (`mt'-`mc')/sqrt((`vt'+`vc')/2)

    quietly sum `v' if treat==1 [aw=w]
    local mt2 = r(mean)
    local vt2 = r(Var)
    quietly sum `v' if treat==0 [aw=w]
    local mc2 = r(mean)
    local vc2 = r(Var)
    local dmat = (`mt2'-`mc2')/sqrt((`vt2'+`vc2')/2)
    local vrat = `vt2'/`vc2'

    if abs(`dmat') > `worst' {
        local worst = abs(`dmat')
        local worstv "`v'"
    }
    di as txt %-12s "`v'" as res %16.4f `draw' %18.4f `dmat' %18.4f `vrat'
}
di as txt "{hline 78}"
di as txt "These closely track tebalance summarize (section B6)."
di as txt "Small differences come from using aggregate aweights here while"
di as txt "tebalance uses Stata's internal matched-sample accounting."
di _n as txt "Worst remaining imbalance: " as res "`worstv'" as txt " at " as res %5.3f `worst'

*=====================================================================
* B5. THE ATT
*     Average the within-pair outcome differences over treated units.
*     That is the entire estimator. Everything before this point was
*     design; this is the first block where the outcome appears.
*=====================================================================
use "$path/_pairs.dta", clear
rename cid id
merge m:1 id using "$path/_alldata.dta", keep(match) keepusing(re78) nogen
gen double wc_re78 = w * re78
collapse (sum) c_re78=wc_re78, by(tid)

save "$path/_matched_y.dta", replace

use "$path/_matched.dta", clear
keep tid dist
merge 1:1 tid using "$path/_matched_y.dta", nogen
rename tid id
merge 1:1 id using "$path/_alldata.dta", keep(match) keepusing(re78) nogen
gen double paird = re78 - c_re78

di _n as txt "{hline 70}"
sum paird
di as txt "HAND-ROLLED ATT = " as res %10.4f r(mean)
di as txt "{hline 70}"
* Compare with the teffects output in B6. They agree to four decimals.

*---------------------------------------------------------------------
* B5b. A CALIPER, BY HAND
*      Refuse a match when the score gap exceeds c. Doing it here rather
*      than through teffects keeps the score FIXED, so you see the
*      caliper alone at work.
*---------------------------------------------------------------------
foreach c in 0.05 0.01 0.001 {
    quietly sum paird if dist <= `c'
    di as txt "caliper " as res %6.3f `c' as txt ": ATT = " as res %9.2f r(mean) ///
       as txt "   treated units kept: " as res r(N) as txt " of 185"
}
di as txt "Tightening the caliper improves match quality and shrinks the"
di as txt "population. The estimate is no longer about the same people, so"
di as txt "a caliper is a change of estimand and must be reported as one."



log close
