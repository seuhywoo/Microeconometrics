***********************************************************************
* E673 Reproduction: Topic 1 (Dehejia & Wahba 1999, reanalyzing LaLonde 1986)
*
* 이 파일은 이전에 따로 있던 두 파일을 하나로 합친 최종 버전이다:
*   - 우리가 처음부터 짠 코드 (benchmark, propensity score, stratification)
*   - HW1.do (교수님 샘플 코드) 중 Part B (matching, overlap, balance)
*
* 각 STEP 옆 주석에 "과제 PDF의 어느 항목에 해당하는지"를 표시해뒀다.
* 과제 PDF 원문 구조:
*   Reproduce 1: experimental benchmark ($1,794)
*   Reproduce 2: DW's propensity-score ATT (NSW treated vs CPS), Table 3와 비교
*   Audit 1 (Overlap), Audit 2 (Balance), Audit 3 (Sample selection)
***********************************************************************

* ---------------------------------------------------------------------
* SETUP
* ---------------------------------------------------------------------
clear all
set more off

global path "/Users/seunghyeon/myproject/Microeconometrics/Topic1/seunghyeon"
* .dta 파일들이 있는 폴더. 이 경로를 고정해서 쓴다.

cd "$path"


* ==========================================================
* STEP 1. 데이터 불러오기 및 확인
*   [과제 PDF: "The two NSW samples" 문단 — 어느 파일을 열었는지 확인]
* ==========================================================

use "nsw_dw.dta", clear
* Dehejia & Wahba의 실험 서브표본 (185 treated + 260 control = 445명).

describe
summarize
tab treat
* treat==1이 185명, treat==0이 260명으로 나와야 정상.
* 다른 숫자가 나오면 nsw.dta(297/425, LaLonde 원본)를 잘못 연 것.


* ==========================================================
* STEP 2. Experimental benchmark 계산
*   [과제 PDF Reproduce 항목 1: "The published value is $1,794"]
* ==========================================================

ttest re78, by(treat)
* treat==1 vs treat==0의 1978년 소득(re78) 평균 차이. $1,794 근처가 나와야 함.

regress re78 treat
* 같은 숫자를 회귀분석으로 재확인 (단순 평균차이 = 단순회귀의 기울기).
scalar benchmark_dw = _b[treat]
* Step 9(Sample selection)에서 이 DW-sample benchmark(약 $1,794)를
* 다시 불러 쓸 것이므로 scalar로 저장해둔다.


* ==========================================================
* STEP 3. Working sample 구성 (NSW treated + CPS-1 comparison group)
*   [과제 PDF Reproduce 항목 2: "the treated units from nsw_dw.dta
*    plus cps_controls.dta"]
*
*   여기서부터는 나중에 Step 4(propensity score), Step 5(stratification),
*   Step 7(matching)에서 전부 재사용할 covariate들(age2, educ2, u74, u75,
*   개인 식별용 id)을 미리 만들어서 저장해둔다.
* ==========================================================

use "nsw_dw.dta", clear
drop if treat == 0
* nsw_dw.dta의 원래 control(260명)은 여기서는 버린다.
* (그건 Step 2의 experimental benchmark 계산에만 쓰는 것.)

append using "cps_controls.dta"
* CPS-1 comparison group(15,992명, 전부 treat==0)을 이어붙인다.
replace treat = 0 if missing(treat)
* append 과정에서 혹시 treat이 비어있는 채로 들어온 관측치가 있으면 0으로 채움
* (cps_controls.dta는 애초에 treat이 다 0으로 채워져 있어야 정상이지만 안전장치).

gen age2  = age^2
gen educ2 = education^2
gen u74   = re74 == 0
gen u75   = re75 == 0
* u74, u75 = "1974/1975년 소득이 정확히 0이었는가" dummy.
* 전체 표본의 상당수가 re74=0인 mass point를 이루므로, re74를 선형으로만
* 넣으면 이 패턴을 놓친다 — 그래서 별도 dummy로 추가.

gen age3      = age^3
gen educ_re74 = education * re74
* DW의 Table 3 각주 (g) CPS-1 propensity score specification에 정확히
* 나온 항들. age3 = age의 세제곱, educ_re74 = 교육연수 x 1974년 소득.
* (u74*black은 PSID-1 전용 항이라 CPS-1에는 안 들어간다 — 이전에 잘못
*  안내했던 부분을 정정함.)

gen long id = _n
* 개인 식별 번호. Step 7(matching)에서 "누가 누구와 매칭됐는지" 기록하려면 필요.

tab treat
* treat==1: 185명, treat==0: 15,992명이 나와야 정상.

save "merged.dta", replace
* 이 파일을 Step 4 이후 전 과정에서 계속 불러 쓴다.


* ==========================================================
* STEP 4. Propensity score 추정
*   [과제 PDF Reproduce 항목 2의 핵심 — Table 3 전체(4/5/6/7/8열)가
*    전부 이 하나의 propensity score 위에서 계산됨]
* ==========================================================

global X "age age2 education educ2 married nodegree black hispanic re74 re75 u74 u75 educ_re74 age3"
* DW의 Table 3 각주 (g) CPS-1, CPS-2, CPS-3 propensity score specification과
* 정확히 일치하는 버전:
*   age, age^2, education, education^2, married, no degree, black,
*   Hispanic, RE74, RE75, u74, u75, education*RE74, age^3
* (참고: u74*black은 PSID-1 전용 항(각주 e)이라 여기엔 포함되지 않는다.)

logit treat $X
predict double ps, pr
label var ps "e-hat = Pr(treat=1 | X)"

di _n as txt "{hline 70}"
di as txt "Distribution of the estimated score"
di as txt "{hline 70}"
tabstat ps, by(treat) stat(n mean p10 p50 p90 min max) format(%9.5f)

save "merged.dta", replace


* ==========================================================
* STEP 4b. Table 3, CPS-1 행의 column (1)(2) 재현
*   [과제 PDF Reproduce 항목 2 — Table 3 전체를 재현하려면 score 기반
*    컬럼(4~8)뿐 아니라 이 두 컬럼도 필요하다]
*   column (1): Unadjusted  — 공변량 없이 re78을 treat에 회귀 (raw diff-in-means)
*   column (2): Adjusted    — 논문 각주 (a)의 스펙으로 최소제곱
*     (a) RE78 on a constant, treat, age, age^2, education, no degree,
*         black, Hispanic, RE74, RE75.
*     주의: 이 스펙은 propensity score용 $X와 다르다 — married, educ2,
*     u74/u75, 상호작용항이 빠지고 "순수 통제변수만" 들어간다.
* ==========================================================

use "merged.dta", clear

* --- Column (1): Unadjusted ---
regress re78 treat
di as txt "Column (1) Unadjusted (raw diff-in-means) = " as res %10.4f _b[treat]
* -> 논문 값 -8,498과 비교. CPS control이 훨씬 부유한 집단이라 부호가
*    음수로 나오는 게 정상 — 이게 바로 selection bias의 크기를 보여준다.

* --- Column (2): Adjusted ---
regress re78 treat age age2 education nodegree black hispanic re74 re75
di as txt "Column (2) Adjusted (footnote a 스펙) = " as res %10.4f _b[treat]
* -> 논문 값 972와 비교. 공변량을 몇 개 통제하는 것만으로 -8,498에서
*    거의 부호가 바뀔 정도로 크게 개선된다 — 그래도 아직 벤치마크 1,794와는
*    거리가 있다는 게 논문이 propensity score 방법을 쓰는 이유다.


* ==========================================================
* STEP 5. Table 3, CPS-1 행의 column (4)(5)(6) 재현 — Stratification
*   [과제 PDF Reproduce 항목 2: "set against the matching estimates
*    in their results table"]
*   column (4): Stratifying on the score, Unadjusted
*   column (5): Stratifying on the score, Adjusted
*   column (6): Observations
*
*   이전에는 xtile로 인원수 기준 5등분(quintile)을 했는데, treated가
*   185명뿐이고 propensity score 상단에 몰려있다 보니 일부 구간에 treated가
*   너무 적게 남아서 추정치가 불안정했다 (Column 4가 -37로 나오는 등
*   논문 값 1,713과 부호까지 달랐던 문제).
*
*   그래서 Becker & Ichino (2002)의 pscore / atts 패키지로 바꾼다.
*   이 패키지는 DW가 설명한 절차 — "구간을 나누고, 안에서 covariate가
*   balance되는지 확인하고, 안 되면 구간을 다시 쪼개는" 반복적 알고리즘을
*   자동으로 수행해준다.
* ==========================================================

cap net install st0026_2, from(http://www.stata-journal.com/software/sj5-3) replace
* pscore는 SSC가 아니라 Stata Journal 소프트웨어 아카이브에 있다
* (Becker & Ichino, SJ 5(3), 최신 버그수정판인 st0026_2).
* 이 한 줄로 pscore, attnd, attk, attr, atts, attnw 명령어가 모두 설치된다.
* cap(capture)를 붙이면 이미 설치돼 있어도 에러 없이 넘어간다.
* 설치가 안 되면 (help pscore로 확인) 아래를 cap 없이 직접 실행해서
* 에러 메시지를 확인할 것:
*   net install st0026_2, from(http://www.stata-journal.com/software/sj5-3) replace

use "merged.dta", clear

pscore treat $X, pscore(ps2) blockid(pblock) comsup logit detail
* treat ~ $X 로 propensity score를 다시 추정하면서(logit 옵션 = DW와 동일한
* 모형), covariate가 balance될 때까지 반복해서 구간을 나눈다.
*   pscore(ps2)  : 새로 추정된 propensity score를 ps2라는 변수로 저장
*   blockid(pblock): balancing test를 통과한 최종 구간 번호를 pblock에 저장
*   comsup       : common support 밖의 관측치는 이후 계산에서 자동 제외
*   detail       : 각 구간에서 balance가 맞는지 확인하는 과정을 화면에 출력
*
* 실행 결과 화면에 "The final number of blocks is ..."라는 메시지가 뜨고,
* 각 covariate별로 balancing test를 통과했는지 나온다 — 이 출력 자체가
* Step 8(Balance audit)에서 "어떤 covariate 때문에 score를 다시 추정해야
* 했는지" 서술할 때 바로 쓸 수 있는 자료다.

describe
* pscore 실행 후 어떤 새 변수들이 생겼는지 확인 (버전에 따라 common-support
* 표시 변수 이름이 다를 수 있으니, 여기서 눈으로 확인하고 아래 count에서
* 그 이름을 맞춰 쓴다. 기존에 Step 4에서 만든 ps, Step 5의 min/max 방식
* common-support 결과와 비교해봐도 좋다).

* --- Column (3): Quadratic in score ---
* 논문 각주 (b): "Least squares regression of RE78 on a quadratic in the
* estimated propensity score and a treatment indicator, for observations
* used under stratification." — 즉 covariate 전체가 아니라 score와 그
* 제곱항만 통제하고, common support 안의 관측치만 사용한다.
gen double ps2sq = ps2^2
regress re78 treat ps2 ps2sq if comsup==1
di as txt "Column (3) Quadratic in score = " as res %10.4f _b[treat]
* -> 논문 값 1,117과 비교. Column (2)의 972보다 벤치마크 1,794에 더
*    가까워지는데, 이게 바로 "score를 비선형으로 쓰는 것"의 가치를
*    보여주는 지점 — 논문 5절(Sensitivity Analysis)의 핵심 메시지이기도 하다.

* --- Column (6): Observations ---
* pscore가 만든 공식 common-support 안에서 몇 명이 남았는지 확인.
* (실행해보니 이 버전에서는 comsup이라는 이름으로 생성됨: 1=common support 안, 0=밖)
count if comsup==1
count if comsup==1 & treat==1
count if comsup==1 & treat==0

* --- Column (4): Stratification, Unadjusted ---
atts re78 treat, pscore(ps2) blockid(pblock) comsup
* Becker-Ichino의 stratification 추정량. 우리가 손으로 짜려던 로직
* (구간별 평균차이를 treated 수로 가중평균)을, balancing test로 제대로
* 만들어진 구간(pblock) 위에서 실행해준다. Bootstrap 표준오차도 같이 나온다.

* --- Column (5): Adjusted (within-block regression) ---
* atts 자체는 unadjusted만 제공하므로, adjusted는 구간을 pblock으로
* 교체해서 이전과 같은 방식(statsby로 구간별 회귀 후 가중평균)으로 계산한다.
*
* 여기서만 쓰는 두 개의 다리(bridge) 파일은 tempfile로 만든다 — 시스템
* 임시 폴더에 저장되고 do-file이 끝나면 자동으로 지워지므로, 프로젝트
* 폴더($path)에는 흔적이 남지 않는다.
tempfile block_n1 block_results

bysort pblock: egen n1_block = total(treat)

preserve
    duplicates drop pblock, force
    keep pblock n1_block
    save `block_n1', replace
restore

statsby delta_k=_b[treat], by(pblock) saving(`block_results', replace): ///
    regress re78 treat $X if comsup==1

use `block_results', clear
describe
* pblock 기준으로 최종 구간 개수만큼 관측치가 나와야 정상.

merge 1:1 pblock using `block_n1'
drop _merge

egen total_n1 = total(n1_block)
gen weight_k = n1_block / total_n1
gen contrib_k = weight_k * delta_k

summarize contrib_k
scalar ATT_strat_adj = r(sum)
display "Column (5) Stratification, Adjusted ATT = " ATT_strat_adj


* ==========================================================
* STEP 6. Overlap 체크
*   [과제 PDF Audit 항목 1: "Plot the two ê densities. Enforce common
*    support, report how many comparison units you dropped, and say
*    what population your ATT now describes."]
* ==========================================================

use "merged.dta", clear

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
* common-support 위반 개수 — 이 숫자를 audit에 "몇 명을 dropped" 했는지로 보고.

count if treat==0 & ps > 0.1
di as txt "Controls with a score above 0.1: " as res r(N) as txt " of 15,992"

label define grp 0 "CPS controls (15,992)" 1 "NSW treated (185)", replace
label values treat grp

* 두 그룹 크기 차이가 커서 단순 겹쳐그리기론 안 보이므로, bin별 개수를
* 세고 y축 상한을 잘라서(capping) 그린다.
preserve
    gen double bin = floor(ps/0.1)*0.1 + 0.05
    contract treat bin
    gen double nshow = min(_freq, 110)
    twoway (bar nshow bin, barwidth(0.09) color(navy%70)), ///
        by(treat, cols(1) ///
           note("Shared count scale, cut at 110.")) ///
        yscale(range(0 110)) ylabel(0(25)100) ///
        xlabel(0(0.2)1) xtitle("estimated propensity score") ///
        ytitle("number of observations") name(overlap_hand, replace)
    graph export "overlap_hand.png", replace width(1000)
restore
* -> 이 그림 + 위의 count 결과 + "이제 ATT는 common support 안에 있는
*    treated 참가자들에 대한 평균효과로 좁혀진다"는 서술이 audit 항목 1의 답.


* ==========================================================
* STEP 7. Table 3, CPS-1 행의 column (7)(8) 재현 — Matching
*   [과제 PDF Reproduce 항목 2 계속]
*   column (7): Matching on the score, Unadjusted
*   column (8): Matching on the score, Adjusted  <- 현재 미구현 (TODO, 아래 참고)
*
*   Nearest neighbor on the score, with replacement.
*   185 x 15,992 쌍을 전부 만드는 대신, score 기준으로 정렬한 뒤 각
*   treated 단위 바로 위/아래 control만 후보로 놓고 더 가까운 쪽을
*   고른다 (동점이면 양쪽 가중평균). 비용은 정렬 한 번 (n log n).
* ==========================================================

use "merged.dta", clear

* 이 Step에서만 쓰는 다리(bridge) 파일들은 전부 tempfile로 만든다.
* tempfile은 시스템 임시 폴더에 저장되고, do-file 실행이 끝나면(또는
* Stata를 닫으면) 자동으로 삭제되므로 프로젝트 폴더에는 남지 않는다.
* pairs, cweight는 Step 8에서도 다시 쓰는데, local macro는 같은 do-file
* 안에서는 끝까지 살아있으므로 여기서 만든 그대로 Step 8에서 재사용한다.
tempfile roster cells matched side1 pairs matched_y paird_outcomes cweight

* B4(Step 8)에서 balance weight 계산에 쓸, control 하나당 한 줄짜리 테이블.
preserve
    keep if treat==0
    keep id ps
    rename id   cid
    save `roster', replace
restore

preserve
    keep if treat==0
    keep ps id
    collapse (count) ncell=id, by(ps)
    gen byte iscell = 1
    save `cells', replace
restore

preserve
    keep if treat==1
    keep id ps
    rename id tid
    gen byte iscell = 0
    append using `cells'
    sort ps iscell

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

    gen int    nties  = .
    gen double mps    = .
    replace nties  = lo_n  if d_lo <  d_hi | missing(d_hi)
    replace mps    = lo_ps if d_lo <  d_hi | missing(d_hi)
    replace nties  = hi_n  if d_hi <  d_lo | missing(d_lo)
    replace mps    = hi_ps if d_hi <  d_lo | missing(d_lo)
    replace nties  = lo_n + hi_n if d_lo==d_hi & !missing(d_lo)

    keep tid dist nties mps lo_ps hi_ps d_lo d_hi lo_n hi_n
    save `matched', replace
restore

use `matched', clear
gen double w_each = 1/nties
keep if !missing(mps)
keep tid mps w_each
rename mps ps
save `side1', replace

* 동점(양쪽이 정확히 같은 거리)인 경우: 양쪽 다 매칭으로 인정.
use `matched', clear
gen double w_each = 1/nties
keep if missing(mps)
expand 2
gen byte side = .
bysort tid: replace side = _n
gen double ps = cond(side==1, lo_ps, hi_ps)

keep tid ps w_each
append using `side1'

joinby ps using `roster'
keep tid cid w_each
rename w_each w

save `pairs', replace

use `matched', clear
di _n as txt "{hline 70}"
di as txt "How many controls each treated unit was matched to"
di as txt "{hline 70}"
tab nties
sum dist, detail
* -> 이 결과가 "matching choices" 서술(몇 명과 매칭됐는지, 동점 처리 방식)에 쓰인다.

* --- Column (7): Matching, Unadjusted ---
use `pairs', clear
rename cid id
merge m:1 id using "merged.dta", keep(match) keepusing(re78) nogen
gen double wc_re78 = w * re78
collapse (sum) c_re78=wc_re78, by(tid)
save `matched_y', replace

use `matched', clear
keep tid dist
merge 1:1 tid using `matched_y', nogen
rename tid id
merge 1:1 id using "merged.dta", keep(match) keepusing(re78) nogen
gen double paird = re78 - c_re78

di _n as txt "{hline 70}"
sum paird
di as txt "Column (7) Matching, Unadjusted ATT = " as res %10.4f r(mean)
di as txt "{hline 70}"
scalar ATT_dw_match = r(mean)
* Step 9(Sample selection)에서 이 DW-sample 매칭 ATT를 다시 불러 쓸 것이므로
* scalar로 저장해둔다 (한 do-file 안에서는 scalar가 끝까지 살아있다).
* teffects psmatch로 나온 결과와 비교해서 소수점 몇 자리까지 일치하는지
* 확인하면, 손으로 짠 알고리즘이 맞게 작동하는지 검증할 수 있다.

save `paird_outcomes', replace

* --- Column (8): Matching, Adjusted ---
* DW는 matching으로 만들어진 표본(treated + 매칭된 control, control은
* 매칭 횟수만큼 가중치)에 대해 re78을 covariate에 회귀시켜서 treat의
* 계수를 조정된 ATT로 쓴다. `pairs'(방금 만든 매칭 결과)에서 control별
* 가중치 w를 만들고, 그 가중치로 가중회귀(aweight)를 돌리면 된다.
* (Step 8에서도 같은 계산을 다시 하는데, 파일을 다시 만들어도 내용은
*  같으니 문제 없다 — 그냥 두 번 계산하는 셈.)
use `pairs', clear
collapse (sum) w, by(cid)
rename cid id
save `cweight', replace

use "merged.dta", clear
merge 1:1 id using `cweight', nogen
replace w = 0 if missing(w) & treat==0
* 매칭에 한 번도 안 쓰인 control은 가중치 0.
replace w = 1 if treat==1
* treated는 전원 가중치 1 (매칭에서 전부 최소 한 번씩은 쓰였으므로).

regress re78 treat $X [aw=w]
* treat의 계수가 Column (8). aweight는 관측치별 신뢰도(=매칭된 횟수)를
* 반영한 가중최소제곱으로, control이 여러 번 매칭됐으면 그만큼 더 비중있게
* 반영되고, 한 번도 안 쓰인 control(w=0)은 회귀에서 빠지는 것과 같다.

di _n as txt "{hline 70}"
di as txt "Column (8) Matching, Adjusted ATT = " as res %10.4f _b[treat]
di as txt "{hline 70}"
* -> 논문의 CPS-1 Column (8) 값 1,616 (표준오차 751)과 비교.

use `paird_outcomes', clear
* Column (8) 회귀 때문에 메모리의 데이터가 merged.dta로 바뀌었으므로,
* caliper 민감도 계산에 필요한 paird/dist 변수가 있는 파일을 다시 불러온다.

* --- Caliper 민감도 (과제의 "common-support rule, including any
*      caliper or trimming rule" 서술에 직접 쓸 수 있는 부분) ---
foreach c in 0.05 0.01 0.001 {
    quietly sum paird if dist <= `c'
    di as txt "caliper " as res %6.3f `c' as txt ": ATT = " as res %9.2f r(mean) ///
       as txt "   treated units kept: " as res r(N) as txt " of 185"
}
di as txt "Caliper를 좁힐수록 매칭 품질은 좋아지지만 표본이 줄어든다."
di as txt "그만큼 estimand 자체가 바뀌는 것이므로 보고서에 명시해야 한다."


* ==========================================================
* STEP 8. Balance 체크
*   [과제 PDF Audit 항목 2: "Within propensity-score strata, compare
*    covariate means across treated and controls. Report balance
*    before and after your adjustment, and say which covariates
*    required refitting the score."]
*
*   표준화된 차이 (standardized difference):
*     d = (xbar_treated - xbar_control) / sqrt((s2_t + s2_c)/2)
*   |d| < 0.1이면 관례적으로 balance가 괜찮다고 판단한다.
* ==========================================================

* pairs, cweight는 Step 7에서 만든 tempfile local macro를 그대로 재사용한다
* (같은 do-file 안에서는 local macro가 끝까지 살아있음).
use `pairs', clear
collapse (sum) w, by(cid)
rename cid id
save `cweight', replace

use "merged.dta", clear
merge 1:1 id using `cweight', nogen
replace w = 0 if missing(w) & treat==0
replace w = 1 if treat==1
* treated는 전원 1번씩, control은 매칭된 횟수만큼, 안 쓰인 control은 0으로 가중.

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
di _n as txt "Worst remaining imbalance: " as res "`worstv'" as txt " at " as res %5.3f `worst'
* -> "raw"열(매칭 전)과 "matched"열(매칭 후)을 비교해서 balance가 개선됐는지
*    보고하고, |d|가 여전히 큰 covariate가 있으면 그게 "score를 다시
*    추정해야 하는 covariate"라고 audit에 서술한다.


* ==========================================================
* STEP 9. Sample selection 비교
*   [과제 PDF Reproduce 항목 1의 두 번째 숫자($886) + Audit 항목 3:
*    "re-estimate your matched ATT on LaLonde's sample, where re74 does
*    not exist and cannot enter the score, and compare the three
*    numbers: benchmark, DW-sample ATT, LaLonde-sample ATT."]
*
*   nsw.dta는 LaLonde 원본 실험 표본(297 treated / 425 control)이고,
*   변수 목록에 re74가 아예 없다 (확인: describe 결과에 re74 없음).
*   그래서 이 표본으로는 $X에서 re74/u74/educ_re74를 뺀 축소판 스펙을 써야 한다.
* ==========================================================

* --- (a) Benchmark, LaLonde sample (Reproduce 항목 1의 두 번째 숫자) ---
use "nsw.dta", clear
describe
* re74 컬럼이 없다는 걸 눈으로 확인.
tab treat
* treat==1이 297명, treat==0이 425명으로 나와야 정상.

ttest re78, by(treat)
regress re78 treat
scalar benchmark_lalonde = _b[treat]
di as txt "Benchmark, LaLonde sample (nsw.dta) = " as res %10.4f benchmark_lalonde
* -> 논문/과제 PDF가 말하는 $886 근처가 나와야 한다.

* --- (b) Working sample: LaLonde의 treated(297명) + CPS-1 comparison ---
use "nsw.dta", clear
drop if treat == 0
* nsw.dta 원래의 실험 control(425명)은 benchmark 계산에만 쓰고 여기선 버린다.
append using "cps_controls.dta"
replace treat = 0 if missing(treat)

gen age2  = age^2
gen educ2 = education^2
gen u75   = re75 == 0
gen age3  = age^3
* re74가 없으니 u74, educ_re74도 만들 수 없다 — $X에서 통째로 뺀다.

gen long id = _n
tab treat
save "_alldata_lalonde.dta", replace
* 이 파일은 merged.dta의 LaLonde-sample 버전이라 핵심 데이터로 남겨둔다
* (스크래치용 다리 파일들과 달리 계속 참고하는 대상이므로 tempfile로 안 만듦).

* --- (c) Propensity score, re74 관련 항 전부 제외 ---
global XL "age age2 education educ2 married nodegree black hispanic re75 u75 age3"
* DW 각주 (g)에서 RE74, u74, education*RE74 세 항만 뺀 버전.
* (re74 자체가 없는 표본이니 당연히 이 세 항은 추정할 수 없다.)

logit treat $XL
predict double psl, pr
label var psl "e-hat (LaLonde sample, no RE74) = Pr(treat=1|X)"

di _n as txt "{hline 70}"
tabstat psl, by(treat) stat(n mean p10 p50 p90 min max) format(%9.5f)
di as txt "{hline 70}"

save "_alldata_lalonde.dta", replace

* --- (d) Overlap 체크 (Step 6과 같은 방식) ---
sum psl if treat==1
local tmin_l = r(min)
local tmax_l = r(max)
sum psl if treat==0
local cmin_l = r(min)
local cmax_l = r(max)

di _n as txt "Treated score range:  " as res %7.5f `tmin_l' " to " %7.5f `tmax_l'
di as txt "Control score range:  " as res %7.5f `cmin_l' " to " %7.5f `cmax_l'

count if treat==1 & (psl > `cmax_l' | psl < `cmin_l')
di as txt "Treated units OUTSIDE the control score range: " as res r(N) as txt " of 297"

* --- (e) Matching: Step 7과 동일한 sort 기반 nearest-neighbor 알고리즘 ---
* (여기서 쓰는 다리 파일들도 전부 tempfile로 만든다 — Step 7의 pairs 등과
*  이름이 겹치지 않도록 _l을 붙인 별도의 local macro를 새로 선언한다.)
tempfile roster_l cells_l matched_l side1_l pairs_l matched_y_l

preserve
    keep if treat==0
    keep id psl
    rename id cid
    save `roster_l', replace
restore

preserve
    keep if treat==0
    keep psl id
    collapse (count) ncell=id, by(psl)
    gen byte iscell = 1
    save `cells_l', replace
restore

preserve
    keep if treat==1
    keep id psl
    rename id tid
    gen byte iscell = 0
    append using `cells_l'
    sort psl iscell

    gen double lo_ps = psl     if iscell==1
    gen double lo_n  = ncell   if iscell==1
    replace lo_ps = lo_ps[_n-1] if missing(lo_ps) & _n>1
    replace lo_n  = lo_n[_n-1]  if missing(lo_n)  & _n>1

    gsort -psl -iscell
    gen double hi_ps = psl     if iscell==1
    gen double hi_n  = ncell   if iscell==1
    replace hi_ps = hi_ps[_n-1] if missing(hi_ps) & _n>1
    replace hi_n  = hi_n[_n-1]  if missing(hi_n)  & _n>1

    keep if iscell==0
    gen double d_lo = psl - lo_ps
    gen double d_hi = hi_ps - psl
    gen double dist = min(cond(missing(d_lo),.,d_lo), cond(missing(d_hi),.,d_hi))

    gen int    nties  = .
    gen double mps    = .
    replace nties  = lo_n  if d_lo <  d_hi | missing(d_hi)
    replace mps    = lo_ps if d_lo <  d_hi | missing(d_hi)
    replace nties  = hi_n  if d_hi <  d_lo | missing(d_lo)
    replace mps    = hi_ps if d_hi <  d_lo | missing(d_lo)
    replace nties  = lo_n + hi_n if d_lo==d_hi & !missing(d_lo)

    keep tid dist nties mps lo_ps hi_ps d_lo d_hi lo_n hi_n
    save `matched_l', replace
restore

use `matched_l', clear
gen double w_each = 1/nties
keep if !missing(mps)
keep tid mps w_each
rename mps psl
save `side1_l', replace

use `matched_l', clear
gen double w_each = 1/nties
keep if missing(mps)
expand 2
gen byte side = .
bysort tid: replace side = _n
gen double psl = cond(side==1, lo_ps, hi_ps)
keep tid psl w_each
append using `side1_l'

joinby psl using `roster_l'
keep tid cid w_each
rename w_each w
save `pairs_l', replace

di _n as txt "How many controls each treated unit was matched to (LaLonde sample)"
use `matched_l', clear
tab nties
sum dist, detail

* --- (f) LaLonde-sample Matching ATT (Unadjusted) ---
use `pairs_l', clear
rename cid id
merge m:1 id using "_alldata_lalonde.dta", keep(match) keepusing(re78) nogen
gen double wc_re78 = w * re78
collapse (sum) c_re78=wc_re78, by(tid)
save `matched_y_l', replace

use `matched_l', clear
keep tid dist
merge 1:1 tid using `matched_y_l', nogen
rename tid id
merge 1:1 id using "_alldata_lalonde.dta", keep(match) keepusing(re78) nogen
gen double paird_l = re78 - c_re78

sum paird_l
scalar ATT_lalonde_match = r(mean)
di as txt "LaLonde-sample Matching ATT (Unadjusted) = " as res %10.4f ATT_lalonde_match

* --- (g) 세 숫자 나란히 비교 — 이게 곧 "discussant's argument, in your own output" ---
di _n as txt "{hline 78}"
di as txt "SAMPLE SELECTION COMPARISON (Audit 항목 3)"
di as txt "{hline 78}"
di as txt "Benchmark, DW sample      (nsw_dw.dta, 185/260) = " as res %9.2f benchmark_dw
di as txt "Benchmark, LaLonde sample (nsw.dta,    297/425) = " as res %9.2f benchmark_lalonde
di as txt "Matching ATT, DW sample      (Step 7, RE74 포함) = " as res %9.2f ATT_dw_match
di as txt "Matching ATT, LaLonde sample (이번 Step, RE74 없음) = " as res %9.2f ATT_lalonde_match
di as txt "{hline 78}"
* 서술 포인트 (audit 문서에 옮겨 쓸 것):
*   - DW가 "2년치 소득 요건"을 부과하면서 benchmark 자체가 886 -> 1,794로
*     크게 바뀐다. 표본이 달라지면 "실험이 답하는 질문" 자체가 달라진다는
*     뜻이고, 이게 두 논문 사이 논쟁(Smith-Todd)의 핵심이다.
*   - RE74를 score에서 뺀 LaLonde-sample 매칭 ATT가 DW-sample 매칭 ATT와
*     어떻게 다른지 비교하면, "선구조적 소득(RE74)을 통제하는 것이
*     selection bias를 얼마나 줄여주는가"를 직접 보여주는 숫자가 된다.
