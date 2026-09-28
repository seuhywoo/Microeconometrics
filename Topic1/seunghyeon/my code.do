***********************************************************************
* E673 Reproduction: Topic 1 — Propensity score 재현 전용 (new code.do)
*
*   목표: DW(2002, ReStat) Table 2 (NSW-CPS)의
*         "No. of Observations", "Mean Propensity Score",
*         "Treatment Effect (Diff. in Means)" 열을 재현.
*
*   이 파일에 반영한 것:
*     (1) Step 10 매칭 계산 수정
*         - mean propensity score를 treated 기준 가중평균으로 계산
*           (DW 2002 식: tau = 1/|N| sum_i ( Y_i - 1/|J_i| sum_{j in J_i} Y_j ))
*         - NN with replacement에서 동점(tie)이면 control 1개만 선택
*     (2) 1999 score spec 비교: footnote (g) 그대로 vs nodegree 제외
*     (3) 2002 score spec 확정: Table 2 note (A)의 항 + U74/U75를
*         "1 = positive earnings"로 코딩 + constant 없는 logit
*     (+) Figure 캡션의 첫 bin을 [0, 0.05)로 계산
*
*   benchmark, Table 3 stratification, balance, LaLonde sample 등
*   지금 단계에 필요 없는 부분은 모두 뺐다 (my code.do에 그대로 있음).
***********************************************************************

capture log close
log using "new_code_log.log", replace text

clear all
set more off

global path "/Users/seunghyeon/myproject/Microeconometrics/Topic1/seunghyeon"
cd "$path"


* ==========================================================
* STEP 0. Experimental benchmark
*   [과제 PDF Reproduce 1: "The published value is $1,794"]
* ==========================================================

use "nsw_dw.dta", clear
tab treat
* treat==1: 185, treat==0: 260 이어야 DW 표본이다 (LaLonde 표본은 297/425).
regress re78 treat
scalar bench_dw = _b[treat]
di as txt "Experimental benchmark, DW sample = " as res %9.1f bench_dw as txt "   [paper: 1,794]"


* ==========================================================
* STEP 1. 표본 구성: NSW treated (185) + CPS-1 (15,992)
* ==========================================================

use "nsw_dw.dta", clear
keep if treat == 1
append using "cps_controls.dta"
replace treat = 0 if missing(treat)

gen long id = _n
* 이 순서(treated 먼저, 그다음 CPS 원 파일 순서)가 매칭 동점 처리의
* 기준이 된다: 동점이면 id가 가장 작은 control을 고른다.

gen age2      = age^2
gen age3      = age^3
gen educ2     = education^2
gen educ_re74 = education*re74

gen u74 = (re74 == 0)
gen u75 = (re75 == 0)
* 1999 코드 방식: 1 = 소득 0 (unemployed)

gen e74 = 1 - u74
gen e75 = 1 - u75
* 2002 방식: Table 2의 U74/U75 평균이 NSW 0.29/0.40, CPS 0.88/0.89인데,
* 이건 "소득 > 0" 비율과 일치한다 (소득 0 비율은 NSW 0.71/0.60).
* 즉 2002 논문의 U74/U75는 실제로 1 = positive earnings로 코딩되어 있다.
* constant가 있는 logit에서는 u와 e 코딩이 score에 아무 영향이 없지만,
* constant가 없으면 달라진다 (아래 STEP 2 참고).

tab treat
tabstat u74 u75 e74 e75, by(treat) stat(mean) format(%5.2f)
* -> e74/e75 평균이 2002 Table 2의 U74/U75 열과 같아야 한다.

sort id


* ==========================================================
* STEP 2. 세 가지 propensity score
* ==========================================================

* (a) DW 1999 footnote (g)를 적힌 그대로 (my code.do의 ps와 동일)
global Xg  "age age2 age3 education educ2 married nodegree black hispanic re74 re75 u74 u75 educ_re74"
logit treat $Xg
predict double ps_g, pr

* (b) 1999: footnote (g)에서 nodegree 제외
*     -> 1999 Figure 2 캡션(12,611 / 2,969 / C>.8 = 7)과 일치
global X99 "age age2 age3 education educ2 married black hispanic re74 re75 u74 u75 educ_re74"
logit treat $X99
predict double ps_99, pr

* (c) 2002: Table 2 note (A)의 14개 항, U74/U75 = positive earnings,
*     constant 없음 -> 2002 Figure 1 캡션(11,168 / 4,398)과 일치
global X02 "age age2 age3 education educ2 married nodegree black hispanic re74 re75 e74 e75 educ_re74"
logit treat $X02, noconstant
predict double ps_02, pr
* 참고: constant가 없으면 sum_i e-hat_i = N1 성질이 성립하지 않는다.


* ==========================================================
* STEP 3. Figure 캡션 숫자와 비교
*   discarded = treated 최소 score보다 작은 CPS 수
*   bin1      = 버리고 남은 CPS 중 score가 [0, 0.05)인 수
* ==========================================================

capture program drop ps_diag
program define ps_diag
    args v lab
    quietly {
        sum `v' if treat==1
        scalar _tmin = r(min)
        scalar _mT   = r(mean)
        sum `v' if treat==0
        scalar _mC   = r(mean)
        count if treat==0 & `v' <  _tmin
        scalar _disc = r(N)
        count if treat==0 & `v' >= _tmin & `v' < 0.05
        scalar _bin1 = r(N)
        count if treat==1 & `v' > 0.8
        scalar _t8   = r(N)
        count if treat==0 & `v' > 0.8
        scalar _c8   = r(N)
    }
    di as txt %-24s "`lab'" as res %8.4f _mT %8.4f _mC %8.0f _disc %7.0f _bin1 %6.0f _t8 %6.0f _c8
end

di _n as txt "{hline 70}"
di as txt %-24s "score" %8s "mean T" %8s "mean C" %8s "disc" %7s "bin1" %6s "T>.8" %6s "C>.8"
di as txt "{hline 70}"
ps_diag ps_g  "(a) 1999 fn(g) as written"
ps_diag ps_99 "(b) 1999, no nodegree"
di as txt %-24s "    paper 1999 Fig 2" %8s "" %8s "" %8s "12611" %7s "2969" %6s "35" %6s "7"
ps_diag ps_02 "(c) 2002 spec"
di as txt %-24s "    paper 2002 Fig 1" %8s "0.37" %8s "0.01" %8s "11168" %7s "4398" %6s "" %6s ""
di as txt "{hline 70}"


* ==========================================================
* STEP 4. DW 2002 Table 2 (CPS) 재현 — 매칭 (Mata)
*
*   Without replacement: treated를 정해진 순서(random / low-to-high /
*     high-to-low)로 처리하며, 아직 안 쓰인 control 중 가장 가까운 것 1개.
*   Nearest neighbor (with replacement): 가장 가까운 control 1개.
*     동점이면 id가 가장 작은 것 하나만 쓴다.
*   Caliper: |score 차이| <= delta인 control 전부. 하나도 없으면
*     nearest neighbor 1개로 대체 (2002 footnote 10).
*
*   No. of Obs  = 매칭에 쓰인 서로 다른(unique) control 수
*   Mean score  = (1/N1) sum_i (1/|J_i|) sum_{j in J_i} e_j  (treated 가중)
*   TE          = (1/N1) sum_i ( Y_i - (1/|J_i|) sum_{j in J_i} Y_j )
* ==========================================================

mata:
mata clear

void dw_row(string scalar lab, real scalar n, real scalar mps, real scalar te,
            string scalar pn, string scalar pm, string scalar pte)
{
    printf("%-22s %7.0f %8.4f %9.1f    | %6s %5s %6s\n", lab, n, mps, te, pn, pm, pte)
}

void dw_table2(string scalar psv, real scalar seed)
{
    D  = st_data(., "treat")
    P  = st_data(., psv)
    Y  = st_data(., "re78")
    t  = select(P, D:==1)
    c  = select(P, D:==0)
    yt = select(Y, D:==1)
    yc = select(Y, D:==0)
    N1 = rows(t)
    N0 = rows(c)
    idx = .
    w   = .

    printf("\n{hline 76}\n")
    printf("DW 2002 Table 2 (CPS), score = %s\n", psv)
    printf("%-22s %7s %8s %9s    | %6s %5s %6s\n", "", "obs", "mean ps", "TE", "paper:", "ps", "TE")
    printf("{hline 76}\n")
    dw_row("NSW",      N1, mean(t), .,                  "185",   "0.37", "1794")
    dw_row("Full CPS", N0, mean(c), mean(yt)-mean(yc),  "15992", "0.01", "-8498")

    /* ---- without replacement ---- */
    labs  = ("WOR random", "WOR low to high", "WOR high to low")
    pte   = ("1559", "1605", "1559")
    for (r=1; r<=3; r++) {
        if (r==1) {
            rseed(seed)
            o = order((runiform(N1,1), (1::N1)), (1,2))
        }
        else if (r==2) o = order(( t, (1::N1)), (1,2))
        else           o = order((-t, (1::N1)), (1,2))

        used = J(N0, 1, 0)
        m    = J(N1, 1, .)
        for (k=1; k<=N1; k++) {
            i = o[k]
            d = abs(c :- t[i]) + 10*used      /* 쓴 control은 제외 */
            minindex(d, 1, idx, w)
            m[i] = min(idx)                    /* 동점이면 가장 앞 id */
            used[m[i]] = 1
        }
        dw_row(labs[r], N1, mean(c[m]), mean(yt - yc[m]), "185", "0.32", pte[r])
    }

    /* ---- nearest neighbor, with replacement ---- */
    j = J(N1, 1, .)
    ntied = 0
    for (i=1; i<=N1; i++) {
        d = abs(c :- t[i])
        minindex(d, 1, idx, w)
        j[i] = min(idx)
        if (rows(idx) > 1) ntied++
    }
    dw_row("Nearest neighbor", rows(uniqrows(j)), mean(c[j]), mean(yt - yc[j]),
           "119", "0.37", "1360")
    printf("    (treated with tied nearest controls: %4.0f -> kept 1 each)\n", ntied)

    /* ---- caliper ---- */
    deltas = (0.00001, 0.00005, 0.0001)
    clabs  = ("Caliper 0.00001", "Caliper 0.00005", "Caliper 0.0001")
    pn     = ("325", "1043", "1731")
    pte    = ("1119", "1158", "1122")
    for (s=1; s<=3; s++) {
        cnt = J(N0, 1, 0)
        sps = 0
        ste = 0
        nfb = 0
        for (i=1; i<=N1; i++) {
            sel = (abs(c :- t[i]) :<= deltas[s])
            n   = sum(sel)
            if (n == 0) {                      /* footnote 10: NN로 대체 */
                sel = J(N0, 1, 0)
                sel[j[i]] = 1
                n   = 1
                nfb++
            }
            cnt = cnt + sel
            sps = sps + sum(sel :* c)/n
            ste = ste + yt[i] - sum(sel :* yc)/n
        }
        dw_row(clabs[s], sum(cnt :> 0), sps/N1, ste/N1, pn[s], "0.37", pte[s])
        printf("    (pairs = %6.0f, treated with no control in caliper = %4.0f)\n", sum(cnt), nfb)
    }
    printf("{hline 76}\n")
}
end

sort id

* 2002 score로 Table 2 재현 (메인 결과)
mata: dw_table2("ps_02", 20260927)

* 비교용: 1999 footnote (g) score를 그대로 썼을 때 (my code.do와 같은 score)
mata: dw_table2("ps_g", 20260927)

* 참고: 논문의 random 행은 난수 순서에 따라 달라지므로 정확히 같을 수 없다.


* ==========================================================
* STEP 5. DW 1999 Table 3, CPS-1 행 재현 — 1999 score vs 2002 score
*
*   (1) Unadjusted     : re78 on treat (전체 표본)
*   (2) Adjusted       : footnote (a) 공변량 + treat (전체 표본)
*   (3) Quadratic      : re78 on treat, ps, ps^2 (common support 표본)
*   (4) Strat. Unadj.  : block별 평균 차이를 treated 수로 가중평균
*   (5) Strat. Adj.    : block별 (a) 회귀의 treat 계수를 treated 수로 가중평균
*   (6) Observations   : treated 전부 + score가 [treated 최소, 최대]인 control
*   (7) Match. Unadj.  : NN with replacement (동점이면 id가 가장 작은 control)
*   (8) Match. Adj.    : (a) 공변량 WLS, treated 가중치 1,
*                        control 가중치 = 매칭된 횟수 (footnote d)
*
*   Block 규칙 (1999 논문에는 경계가 없어서 2002 appendix 알고리즘을 따름):
*     [0,1]을 nstart개의 같은 폭 구간으로 나눈 뒤, 어떤 구간에서 treated와
*     control의 평균 score가 유의하게 다르면(|t| > 1.96) 그 구간을 반으로
*     쪼갠다. 더 이상 쪼갤 구간이 없을 때까지 반복한다.
* ==========================================================

global XA "age age2 education nodegree black hispanic re74 re75"
* footnote (a): age, age^2, education, no degree, black, Hispanic, RE74, RE75

mata:
void dw_blocks(string scalar psv, string scalar supv, string scalar blkv,
               real scalar nstart, real scalar dosplit)
{
    P = st_data(., psv)
    S = st_data(., supv)
    D = st_data(., "treat")
    edges = rangen(0, 1, nstart+1)
    edges[rows(edges)] = 1 + 1e-7
    changed = 1
    while (changed & dosplit) {
        changed = 0
        newe = edges[1]
        for (k=1; k<rows(edges); k++) {
            lo = edges[k]
            hi = edges[k+1]
            inb = (S:==1) :& (P:>=lo) :& (P:<hi)
            a = select(P, inb :& (D:==1))
            b = select(P, inb :& (D:==0))
            dosp = 0
            if (rows(a)>=2 & rows(b)>=2 & (hi-lo)>1e-4) {
                se = sqrt(variance(a)/rows(a) + variance(b)/rows(b))
                if (se > 0) {
                    if (abs((mean(a)-mean(b))/se) > 1.96) dosp = 1
                }
            }
            if (dosp) {
                newe = newe \ (lo+hi)/2 \ hi
                changed = 1
            }
            else newe = newe \ hi
        }
        edges = newe
    }
    blk = J(rows(P), 1, .)
    for (k=1; k<rows(edges); k++) {
        inb = (S:==1) :& (P:>=edges[k]) :& (P:<edges[k+1])
        if (sum(inb) > 0) blk[selectindex(inb)] = J(sum(inb), 1, k)
    }
    st_store(., blkv, blk)
}

void dw_nn(string scalar psv, string scalar wvar)
{
    P  = st_data(., psv)
    D  = st_data(., "treat")
    Y  = st_data(., "re78")
    ti = selectindex(D:==1)
    ci = selectindex(D:==0)
    t  = P[ti]
    c  = P[ci]
    w  = J(rows(P), 1, 0)
    w[ti] = J(rows(ti), 1, 1)
    att = 0
    idx = .
    ww  = .
    for (i=1; i<=rows(ti); i++) {
        minindex(abs(c :- t[i]), 1, idx, ww)
        k = ci[min(idx)]
        w[k] = w[k] + 1
        att = att + Y[ti[i]] - Y[k]
    }
    st_store(., wvar, w)
    st_numscalar("T3_c7", att/rows(ti))
}
end

capture program drop dw_table3
program define dw_table3
    args ps nstart dosplit
    tempvar sup blk w ps2
    quietly {
        * (1), (2): score와 무관
        regress re78 treat
        scalar T3_c1 = _b[treat]
        regress re78 treat $XA
        scalar T3_c2 = _b[treat]

        * common support, (6)
        sum `ps' if treat==1
        scalar T3_pmin = r(min)
        scalar T3_pmax = r(max)
        gen byte `sup' = (treat==1) | (`ps'>=T3_pmin & `ps'<=T3_pmax)
        count if `sup'
        scalar T3_c6 = r(N)

        * (3)
        gen double `ps2' = `ps'^2
        regress re78 treat `ps' `ps2' if `sup'
        scalar T3_c3 = _b[treat]

        * (4), (5)
        gen int `blk' = .
        mata: dw_blocks("`ps'", "`sup'", "`blk'", `nstart', `dosplit')
        levelsof `blk', local(bl)
        scalar T3_s4 = 0
        scalar T3_s5 = 0
        scalar T3_n1 = 0
        scalar T3_nb = 0
        foreach b of local bl {
            count if `blk'==`b' & treat==1
            scalar T3_nt = r(N)
            count if `blk'==`b' & treat==0
            scalar T3_nc = r(N)
            if T3_nt>0 & T3_nc>0 {
                sum re78 if `blk'==`b' & treat==1
                scalar T3_m1 = r(mean)
                sum re78 if `blk'==`b' & treat==0
                scalar T3_m0 = r(mean)
                regress re78 treat $XA if `blk'==`b'
                scalar T3_s4 = T3_s4 + T3_nt*(T3_m1 - T3_m0)
                scalar T3_s5 = T3_s5 + T3_nt*_b[treat]
                scalar T3_n1 = T3_n1 + T3_nt
                scalar T3_nb = T3_nb + 1
            }
        }
        scalar T3_c4 = T3_s4/T3_n1
        scalar T3_c5 = T3_s5/T3_n1

        * (7), (8)
        gen double `w' = 0
        mata: dw_nn("`ps'", "`w'")
        regress re78 treat $XA [aw=`w'] if `w'>0
        scalar T3_c8 = _b[treat]
    }
end

* ---- 메인: 세 score 각각 (ps_g = 1999 fn(g) 그대로, ps_99 = nodegree 제외, ps_02 = 2002), block 규칙 = 2002 appendix (5개에서 시작, 쪼개기) ----
matrix T3 = J(4, 8, .)
local r = 0
foreach ps in ps_g ps_99 ps_02 {
    local ++r
    dw_table3 `ps' 5 1
    forvalues k = 1/8 {
        local v = T3_c`k'
        matrix T3[`r', `k'] = `v'
    }
    di as txt "`ps': blocks used = " as res T3_nb
}
matrix T3[4, 1] = (-8498, 972, 1117, 1713, 1774, 4117, 1582, 1616)
matrix rownames T3 = "1999 fn(g)" "1999 no nodeg" "2002 score" "Paper 1999"
matrix colnames T3 = "(1)Unadj" "(2)Adj" "(3)Quad" "(4)StrU" "(5)StrA" "(6)Obs" "(7)MatU" "(8)MatA"

di _n as txt "DW 1999 Table 3, CPS-1 row"
matlist T3, format(%9.0f) border(rows) twidth(14)

* ---- 민감도: block 규칙에 따라 (4), (5)가 얼마나 흔들리는지 (1999 score) ----
di _n as txt "Stratification sensitivity (1999 score)   [paper: (4) 1713, (5) 1774]"
di as txt %-26s "rule" %8s "(4)" %8s "(5)" %8s "blocks"
foreach ns in 5 10 20 {
    foreach sp in 0 1 {
        quietly dw_table3 ps_99 `ns' `sp'
        local lab = cond(`sp', "`ns' equal, split", "`ns' equal, no split")
        di as txt %-26s "`lab'" as res %8.0f T3_c4 %8.0f T3_c5 %8.0f T3_nb
    }
}


* ==========================================================
* STEP 6. Audit 1 — Overlap (1999 score, ps_99)
*   [과제 PDF: "Plot the two e-hat densities. Enforce common support,
*    report how many comparison units you dropped, and say what
*    population your ATT now describes."]
*
*   Common-support 규칙 (DW 1999, p.1058):
*     treated의 최소 score보다 작거나 최대 score보다 큰 comparison을 버린다.
*   참고 (Smith & Todd 2005, p.317, 334-335):
*     ST는 density가 0에 가까운 구간을 2% trimming으로 잘라내며,
*     그 결과 treated도 5-10% 버리게 된다. 아래에서는 DW 규칙을 쓰되,
*     comparison의 최대 score보다 큰 treated가 몇 명인지도 함께 보고해서
*     "treated 쪽도 support를 맞추면 ATT가 어떻게 바뀌는지"를 보인다.
* ==========================================================

* ---- (a) 세 score별 common support 요약 ----
di _n as txt "{hline 78}"
di as txt %-8s "score" %10s "T min" %10s "T max" %10s "C max" ///
   %12s "C < T min" %12s "C > T max" %10s "T > C max"
di as txt "{hline 78}"
foreach ps in ps_g ps_99 ps_02 {
    quietly sum `ps' if treat==1
    scalar O_tmin = r(min)
    scalar O_tmax = r(max)
    quietly sum `ps' if treat==0
    scalar O_cmax = r(max)
    quietly count if treat==0 & `ps' < O_tmin
    scalar O_clo = r(N)
    quietly count if treat==0 & `ps' > O_tmax
    scalar O_chi = r(N)
    quietly count if treat==1 & `ps' > O_cmax
    scalar O_thi = r(N)
    di as txt %-8s "`ps'" as res %10.5f O_tmin %10.5f O_tmax %10.5f O_cmax ///
       %12.0f O_clo %12.0f O_chi %10.0f O_thi
}
di as txt "{hline 78}"
di as txt "paper 1999 (Fig. 2): 12,611 comparison units discarded"

* ---- (b) 이후는 메인 score(ps_99)로 ----
quietly sum ps_99 if treat==1
scalar O_tmin = r(min)
scalar O_tmax = r(max)
quietly sum ps_99 if treat==0
scalar O_cmax = r(max)
gen byte cs_dw = (treat==1) | (ps_99 >= O_tmin & ps_99 <= O_tmax)
label var cs_dw "DW common support: all treated + controls within [T min, T max]"

count if treat==0 & cs_dw==0
di as txt "  -> comparison units dropped by the DW rule: " as res r(N)
count if treat==0 & cs_dw==1
di as txt "  -> comparison units kept: " as res r(N)

* ---- (c) bin별 인원: 어디서 support가 얇은지 ----
gen double ps99_bin = min(floor(ps_99/0.1), 9)/10
label var ps99_bin "e-hat bin (lower edge, width 0.1)"
di _n as txt "Treated vs. kept comparison units by e-hat bin (ps_99)"
tab ps99_bin treat if cs_dw

* ---- (d) 그림 1: DW Figure 2 스타일 히스토그램 (support 안의 control만) ----
* histogram 명령은 막대를 자르지 못해서 첫 bin(control 약 3,000명) 때문에
* 나머지가 안 보인다. 그래서 bin별 인원을 직접 세고 200에서 잘라 그린 뒤,
* 잘린 막대에는 실제 인원을 라벨로 단다.
preserve
    keep if cs_dw
    gen double b = floor(ps_99/0.05)*0.05
    contract treat b
    gen double nshow = min(_freq, 200)
    gen double x = b + cond(treat==1, 0.0375, 0.0125)
    gen str8 lab = string(_freq, "%9.0fc") if _freq > 200
    twoway (bar nshow x if treat==0, barwidth(0.022) fcolor(gs12) lcolor(gs7)) ///
           (bar nshow x if treat==1, barwidth(0.022) fcolor(gs2) lcolor(gs2)) ///
           (scatter nshow x if _freq > 200, msymbol(none) mlabel(lab) ///
                mlabposition(12) mlabcolor(black)), ///
           legend(order(2 "NSW treated (185)" 1 "CPS comparison on support (3,381)") ///
                  rows(1) position(6)) ///
           xlabel(0(0.1)1) ylabel(0(50)200) yscale(range(0 215)) ///
           xtitle("estimated propensity score (ps_99), bin width 0.05") ///
           ytitle("number of units") ///
           note("Bars above 200 are cut at 200; the label gives the true count." ///
                "12,611 CPS units below the treated minimum are dropped and not shown.") ///
           name(overlap_hist, replace)
    graph export "overlap_hist_ps99.png", replace width(1200)
restore

* ---- (e) 그림 2: Smith-Todd Figure 1 스타일, log odds의 density ----
gen double lodds_99 = ln(ps_99/(1 - ps_99)) if ps_99 > 0 & ps_99 < 1
label var lodds_99 "log odds of ps_99"
twoway (kdensity lodds_99 if treat==1, lcolor(black)) ///
       (kdensity lodds_99 if treat==0, lcolor(gs8) lpattern(dash)), ///
       legend(order(1 "NSW treated" 2 "CPS comparison (all)")) ///
       xtitle("log odds ratio, ln(e/(1-e))") ytitle("density") ///
       name(overlap_lodds, replace)
graph export "overlap_logodds_ps99.png", replace width(1200)
* score 자체로 그리면 CPS가 0 근처에 몰려서 모양이 안 보이므로,
* Smith-Todd처럼 log odds로 그린다 (순서는 score와 같다).

* ---- (f) treated 쪽 support: 매칭 ATT가 어느 집단의 효과인가 ----
mata:
void dw_nn_y(string scalar psv, string scalar yout)
{
    P  = st_data(., psv)
    D  = st_data(., "treat")
    Y  = st_data(., "re78")
    ti = selectindex(D:==1)
    ci = selectindex(D:==0)
    c  = P[ci]
    out = J(rows(P), 1, .)
    idx = .
    ww  = .
    for (i=1; i<=rows(ti); i++) {
        minindex(abs(c :- P[ti[i]]), 1, idx, ww)
        out[ti[i]] = Y[ci[min(idx)]]
    }
    st_store(., yout, out)
}
end

sort id
gen double y0_match = .
mata: dw_nn_y("ps_99", "y0_match")
gen double nn_diff = re78 - y0_match if treat==1

quietly count if treat==1 & ps_99 > O_cmax
scalar O_thi = r(N)
quietly sum nn_diff
di _n as txt "NN matching ATT, all 185 treated (DW rule)       = " as res %9.1f r(mean) ///
   as txt "   [paper 1999: 1,582]"
quietly sum nn_diff if ps_99 <= O_cmax
di as txt "NN matching ATT, treated with e-hat <= C max (" as res r(N) as txt ") = " ///
   as res %9.1f r(mean)
quietly sum nn_diff if ps_99 > O_cmax
di as txt "  the " as res O_thi as txt " treated above C max: mean matched difference = " ///
   as res %9.1f r(mean)
* 서술 포인트:
*   - DW 규칙은 comparison만 버린다 (12,611명). 모든 treated가 남으므로
*     ATT는 형식상 NSW 참가자 185명 전체에 대한 효과다.
*   - 하지만 treated 중 7명은 가장 높은 comparison score보다 score가 높아서,
*     실제로는 score가 더 낮은 comparison과 매칭된다 (support 밖 외삽).
*   - 이 7명까지 support를 맞추면 ATT는 "comparison과 겹치는 score 구간의
*     참가자"에 대한 효과로 좁혀지고, 값이 크게 바뀐다.



* ==========================================================
* STEP 7. Audit 2 — Balance
*   [과제 PDF: "Within propensity-score strata, compare covariate means
*    across treated and controls. Report balance before and after your
*    adjustment, and say which covariates required refitting the score."]
*
*   세 가지 score를 비교한다.
*     ps_lin : DW appendix의 출발점. 공변량을 선형으로만 넣은 logit
*              (age, education, married, nodegree, black, hispanic, re74, re75)
*     ps_g   : 1999 footnote (g) 그대로
*     ps_99  : footnote (g)에서 nodegree 제외 (메인)
*
*   각 score마다
*     - strata: STEP 5와 같은 block 규칙 (2002 appendix, 5개에서 시작해 쪼개기),
*               DW common support (treated 전부 + [T min, T max] 안의 control)
*     - raw     : 조정 전, treated 185 vs CPS 15,992의 standardized difference
*     - strat   : block 안의 평균 차이를 treated 수로 가중평균한 standardized difference
*     - sig     : block 안 t-test(|t| > 1.96)에서 유의한 block 수 / 검정한 block 수
*                 (DW의 balancing test; 양쪽에 2명 이상 있는 block만 검정)
*     - matched : NN 매칭(with replacement) 표본의 standardized difference
*                 (control은 매칭된 횟수로 가중)
*   standardized difference = (평균 차이) / sqrt((var_T + var_C)/2),
*   분모는 조정 전 전체 표본의 분산으로 고정. |d| < 0.1이면 관례상 balance 양호.
* ==========================================================

logit treat age education married nodegree black hispanic re74 re75
predict double ps_lin, pr

global BALX "age education black hispanic married nodegree re74 re75 u74 u75"

capture program drop dw_balance
program define dw_balance
    args ps lab
    tempvar sup blk w
    quietly {
        sum `ps' if treat==1
        scalar B_tmin = r(min)
        scalar B_tmax = r(max)
        gen byte `sup' = (treat==1) | (`ps'>=B_tmin & `ps'<=B_tmax)
        gen int `blk' = .
        mata: dw_blocks("`ps'", "`sup'", "`blk'", 5, 1)
        levelsof `blk', local(bl)
        gen double `w' = 0
        mata: dw_nn("`ps'", "`w'")
    }
    local nb : word count `bl'
    di _n as txt "{hline 78}"
    di as txt "Balance, score = `lab'   (blocks = `nb')"
    di as txt %-11s "covariate" %11s "raw" %11s "strat" %11s "sig/tested" %11s "matched"
    di as txt "{hline 78}"
    foreach v of global BALX {
        quietly {
            sum `v' if treat==1
            scalar B_m1 = r(mean)
            scalar B_v1 = r(Var)
            sum `v' if treat==0
            scalar B_m0 = r(mean)
            scalar B_v0 = r(Var)
            scalar B_sd = sqrt((B_v1 + B_v0)/2)
            scalar B_raw = (B_m1 - B_m0)/B_sd

            scalar B_acc = 0
            scalar B_n1  = 0
            scalar B_sig = 0
            scalar B_tst = 0
            foreach b of local bl {
                sum `v' if `blk'==`b' & treat==1
                scalar B_a_n = r(N)
                scalar B_a_m = r(mean)
                scalar B_a_v = r(Var)
                sum `v' if `blk'==`b' & treat==0
                scalar B_b_n = r(N)
                scalar B_b_m = r(mean)
                scalar B_b_v = r(Var)
                if B_a_n>0 & B_b_n>0 {
                    scalar B_acc = B_acc + B_a_n*(B_a_m - B_b_m)
                    scalar B_n1  = B_n1 + B_a_n
                }
                if B_a_n>=2 & B_b_n>=2 {
                    scalar B_tst = B_tst + 1
                    scalar B_se = sqrt(B_a_v/B_a_n + B_b_v/B_b_n)
                    if B_se > 0 {
                        if abs((B_a_m - B_b_m)/B_se) > 1.96 scalar B_sig = B_sig + 1
                    }
                }
            }
            scalar B_strat = (B_acc/B_n1)/B_sd

            sum `v' [aw=`w'] if treat==1
            scalar B_mm1 = r(mean)
            sum `v' [aw=`w'] if treat==0 & `w'>0
            scalar B_mm0 = r(mean)
            scalar B_match = (B_mm1 - B_mm0)/B_sd
        }
        local sigtxt = string(B_sig) + "/" + string(B_tst)
        di as txt %-11s "`v'" as res %11.3f B_raw %11.3f B_strat %11s "`sigtxt'" %11.3f B_match
    }
    di as txt "{hline 78}"
end

sort id
dw_balance ps_lin "linear (DW starting point)"
dw_balance ps_g   "1999 fn(g) as written"
dw_balance ps_99  "1999 without nodegree (main)"

* 서술 포인트 (Audit 2):
*   - 조정 전(raw)에는 거의 모든 공변량이 |d| > 0.5로 크게 불균형하다
*     (black 2.4, re74 -1.6, re75 -1.7, married -1.2 등).
*   - 선형 score로 stratify하면 대부분 맞춰지지만 u74, u75가 남는다
*     (strat 약 0.60, 0.23; 매칭 후에도 0.66, 0.20). 그래서 u74, u75와
*     고차항·교차항(age^2, age^3, educ^2, educ*re74)을 넣어 score를 다시
*     추정해야 했다 — 이것이 footnote (g) spec이 된 이유로 볼 수 있다.
*   - footnote (g) score는 stratum 안에서 모든 공변량이 |d| < 0.1이지만,
*     NN 매칭 표본에서는 nodegree(약 0.25), education, hispanic, u74 등이
*     0.1을 넘는다.
*   - nodegree를 뺀 score(메인)는 매칭 표본의 balance가 오히려 더 좋다
*     (nodegree 약 0.14, 나머지는 대부분 0.1 안쪽). 단, stratum 기준으로는
*     nodegree(약 0.13)와 u75(약 0.15)가 조금 남는다.
*   - 매칭 표본 기준으로 가장 큰 잔여 불균형은 nodegree다. 그래서 (8)처럼
*     매칭 후 회귀 조정을 하는 것이 의미가 있다.



* ==========================================================
* STEP 8. Audit 3 — Sample selection (DW sample vs. LaLonde sample)
*   [과제 PDF: "re-estimate your matched ATT on LaLonde's sample, where
*    re74 does not exist and cannot enter the score, and compare the
*    three numbers: benchmark, DW-sample ATT, LaLonde-sample ATT."]
*
*   2 x 2로 나눠 본다.
*     표본: DW (185 treated) vs LaLonde (297 treated), comparison은 둘 다 CPS-1
*     score: RE74 포함 vs RE74 제외
*   DW 표본에서 RE74를 뺀 score도 계산해야, "표본을 바꾼 효과"와
*   "RE74를 score에서 뺀 효과"를 구분할 수 있다.
*
*   RE74 제외 spec = footnote (g)에서 re74, u74, educ*re74 세 항을 뺀 것.
*   메인 spec과 맞추기 위해 nodegree 포함/제외 두 버전을 모두 계산한다.
*   매칭: STEP 5와 같은 NN with replacement (동점이면 id가 작은 control).
* ==========================================================

capture program drop nn_att
program define nn_att, rclass
    args ps
    tempvar w
    quietly {
        gen double `w' = 0
        mata: dw_nn("`ps'", "`w'")
        count if treat==0 & `w' > 0
        return scalar nuc = r(N)
        sum `ps' if treat==1
        scalar S_tmin = r(min)
        count if treat==0 & `ps' < S_tmin
        return scalar disc = r(N)
    }
    return scalar att = T3_c7
end

global XN74  "age age2 age3 education educ2 married nodegree black hispanic re75 u75"
global XN74b "age age2 age3 education educ2 married black hispanic re75 u75"

* ---- (a) DW 표본 (현재 메모리의 데이터) ----
sort id
nn_att ps_g
scalar A_dw_g = r(att)
nn_att ps_99
scalar A_dw_99 = r(att)

quietly logit treat $XN74
predict double ps_dw_n74, pr
nn_att ps_dw_n74
scalar A_dw_n74 = r(att)

quietly logit treat $XN74b
predict double ps_dw_n74b, pr
nn_att ps_dw_n74b
scalar A_dw_n74b = r(att)

* ---- (b) LaLonde 표본: benchmark ----
use "nsw.dta", clear
describe, short
tab treat
* treat==1: 297, treat==0: 425 이어야 LaLonde 표본이다. re74 변수는 없다.
regress re78 treat
scalar bench_la = _b[treat]
di as txt "Experimental benchmark, LaLonde sample = " as res %9.1f bench_la as txt "   [about 886]"

* ---- (c) LaLonde 표본: NSW treated 297 + CPS-1 ----
keep if treat == 1
append using "cps_controls.dta"
replace treat = 0 if missing(treat)
drop re74
* cps_controls.dta에는 re74가 있지만 NSW 쪽에는 없으므로, 실수로 쓰지 않게 지운다.
gen long id  = _n
gen age2     = age^2
gen age3     = age^3
gen educ2    = education^2
gen u75      = (re75 == 0)
tab treat
sort id

quietly logit treat $XN74
predict double ps_la_n74, pr
nn_att ps_la_n74
scalar A_la_n74 = r(att)
scalar D_la_n74 = r(disc)

quietly logit treat $XN74b
predict double ps_la_n74b, pr
nn_att ps_la_n74b
scalar A_la_n74b = r(att)
scalar D_la_n74b = r(disc)

di _n as txt "LaLonde sample: comparison units below the treated minimum = " ///
   as res D_la_n74 as txt " (with nodegree), " as res D_la_n74b as txt " (without nodegree)"

* ---- (d) 비교표 ----
matrix A3 = J(5, 4, .)
local v = bench_dw
matrix A3[1,1] = `v'
local v = bench_la
matrix A3[1,3] = `v'
local r = 1
foreach p in "A_dw_g ." "A_dw_99 ." "A_dw_n74 A_la_n74" "A_dw_n74b A_la_n74b" {
    local ++r
    local a : word 1 of `p'
    local b : word 2 of `p'
    local v = `a'
    matrix A3[`r',1] = `v'
    local v = `a' - bench_dw
    matrix A3[`r',2] = `v'
    if "`b'" != "." {
        local v = `b'
        matrix A3[`r',3] = `v'
        local v = `b' - bench_la
        matrix A3[`r',4] = `v'
    }
}
matrix rownames A3 = "Experimental benchmark" "Matching, RE74, fn(g)" ///
    "Matching, RE74, no nodeg" "Matching, no RE74" "Matching, no RE74 no nodeg"
* 행 이름에 콜론(:)을 쓰면 Stata가 equation 이름으로 읽어서 표가 쪼개지므로 쉼표를 쓴다.
matrix colnames A3 = "DW ATT" "DW bias" "LaLonde ATT" "LaLonde bias"

di _n as txt "Audit 3: benchmark vs. matching ATT, DW sample vs. LaLonde sample (CPS-1)"
di as txt "bias = matching ATT - experimental benchmark of the same sample; . = not available (no RE74)"
matlist A3, format(%12.0f) border(rows) twidth(28)

* 서술 포인트 (Audit 3):
*   - DW가 RE74를 쓰려고 두 해 소득이 있는 사람으로 표본을 줄였다 (DW 1999
*     p.1054: 일찍 등록했거나, 늦게 등록했어도 등록 전에 실업 상태였던 사람).
*     Smith-Todd(Table 2)에 따르면 이 규칙은 1976년 4월 이후 등록자 중
*     소득이 있던 사람을 빼므로, 등록 전 소득이 낮은 사람 쪽으로 표본이
*     기울고 Ashenfelter dip이 약해진다.
*   - 그 결과 실험 benchmark 자체가 886에서 1,794로 바뀐다. 즉 두 표본은
*     서로 다른 모집단의 효과를 추정한다.
*   - 같은 CPS-1 comparison, 같은 매칭 방법인데도 LaLonde 표본의 매칭 ATT는
*     음수로 크게 빗나간다. DW 표본에서 RE74만 뺀 경우와 비교하면, 차이의
*     대부분이 RE74 변수가 아니라 표본 구성에서 온다는 것을 보일 수 있다.
*   - 이것이 Smith and Todd (2005)의 주장이다: DW의 낮은 bias는 특정 표본
*     선택과 score spec에 크게 의존한다.


log close
