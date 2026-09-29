# Dehejia & Wahba (1999, *JASA*): Tables 3–4 and covariate balance

Three versions of each result are compared:

| Label | What it is |
|---|---|
| **(i) Paper** | the values printed in the paper |
| **(ii) Reproduction** | our reproduction using the method the paper describes (Table 3 notes and the Appendix) |
| **(iii) SH** | our reproduction using the setup in SH's note (`DW_propensity_score_notes.html`) |

- **Part 1:** Table 3 (estimates) and Table 4 (matched samples), with a covariate-balance check for matching.
- **Part 2:** the covariate-balance check for stratification, showing which covariates fail and in which strata. No adjustments are made here.

Comparison groups: **PSID-1** (`psid_controls`), **CPS-1** (`cps_controls`) and **CPS-3** (the 429 comparison units in `cps3re74`).
Put the notebook in the folder with the `.dta` files. Packages: `numpy`, `scipy`, `pandas`.
**Run time:** a few minutes with 500 bootstrap replications. Set `REPS = 100` for a quick run.

### Reproduction vs. SH

| | (ii) Reproduction (the paper's method) | (iii) SH |
|---|---|---|
| Score | logit with a constant; covariates as in Table 3 notes (e) and (g), **with `nodegree`** | logit with a constant; the same covariates **without `nodegree`** |
| Strata for (3)–(6) | Appendix algorithm: 5 equal score ranges on the common support; split a stratum in half when the **score or any covariate** differs significantly ($|t|>1.96$), provided each half keeps ≥ 3 treated and ≥ 3 comparison units | 5 equal-width blocks on [0, 1]; split a block in half while the **mean score** differs between treated and comparison units ($|t|>1.96$) |
| Covariates in (2), (5), (8) | footnote (a): age, age², education, no degree, black, Hispanic, RE74, RE75 | the same |
| Matching (7)–(8) | nearest neighbor with replacement; ties go to the first comparison unit in data order | the same |

Both setups keep all treated units, plus the comparison units whose score lies between the treated minimum and maximum (note c).

**Our choices where neither source is specific:**

- **Merging strata:** a stratum with treated units but no comparison units is merged into its neighbor.
- **Standard errors:** stratification SEs are analytic within-stratum SEs. Matching SEs come from a bootstrap that re-estimates the score and redoes the matching in each replication. JASA does not say how its SEs were computed.
- **Column (8) under SH:** SH's covariates (footnote a) give 1,565 for CPS-1. The paper's 1,616 is reproduced exactly only with the linear set: age, education, black, Hispanic, no degree, married, RE74, RE75, u74, u75.

## 0. Setup

```python
import numpy as np
import pandas as pd
from scipy.special import expit

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
pd.set_option("display.max_colwidth", 90)

REPS = 500      # bootstrap replications for the matching SEs
SEED = 12345
```

### Data and helpers

- **Score terms** are written as strings: `"age^2"` is a square, `"education*re74"` is an interaction. `add_terms` creates the column when it is needed.
- **Estimation:** `logit_ps` fits the logit by Newton–Raphson; `ols_se` returns a coefficient and its homoskedastic SE (weights optional).

```python
def load(name):
    d = pd.read_stata(f"{name}.dta")
    d = d.drop(columns=[c for c in ["data_id", "age2"] if c in d])
    d = d.rename(columns={"ed": "education", "hisp": "hispanic", "nodeg": "nodegree"})
    d = d.astype(float)
    d["u74"] = (d.re74 == 0).astype(float)          # zero earnings in 1974
    d["u75"] = (d.re75 == 0).astype(float)
    return d.reset_index(drop=True)

nsw = load("nsw_dw")
treated = nsw[nsw.treat == 1].reset_index(drop=True)
c3 = load("cps3re74")
GROUPS = {"PSID-1": load("psid_controls"),
          "CPS-1":  load("cps_controls"),
          "CPS-3":  c3[c3.treat == 0].reset_index(drop=True)}

def add_terms(d, terms):
    for t in terms:
        if t in d:
            continue
        if t.endswith("^2"):
            d[t] = d[t[:-2]] ** 2
        elif t.endswith("^3"):
            d[t] = d[t[:-2]] ** 3
        else:
            a, b = t.split("*")
            d[t] = d[a] * d[b]
    return d

def logit_ps(d, xs):
    X = d[xs].to_numpy(float)
    sd = X.std(0); sd[sd == 0] = 1
    X = np.column_stack([np.ones(len(X)), (X - X.mean(0)) / sd])   # constant included
    y = d.treat.to_numpy()
    b = np.zeros(X.shape[1])
    for _ in range(100):
        p = expit(X @ b)
        step = np.linalg.solve((X * (p * (1 - p))[:, None]).T @ X + 1e-9 * np.eye(len(b)), X.T @ (y - p))
        b += step
        if np.abs(step).max() < 1e-10:
            break
    return expit(X @ b)

def ols_se(d, y, xs, w=None):
    """Coefficient and SE on the first regressor (xs[0]); WLS if w is given."""
    X = d[xs].to_numpy(float).copy()
    sd = X[:, 1:].std(0); sd[sd == 0] = 1
    X[:, 1:] /= sd
    X = np.column_stack([np.ones(len(d)), X])
    w = np.ones(len(d)) if w is None else np.asarray(w, float)
    Y = d[y].to_numpy(float)
    A = (X * w[:, None]).T @ X
    b = np.linalg.lstsq(A, (X * w[:, None]).T @ Y, rcond=None)[0]
    e = Y - X @ b
    s2 = (w * e ** 2).sum() / (len(Y) - np.linalg.matrix_rank(X))
    return b[1], np.sqrt(s2 * np.linalg.pinv(A)[1, 1])

def te(b, se):
    return f"{b:,.0f} ({se:,.0f})"

print("treated:", len(treated), {g: len(c) for g, c in GROUPS.items()})
```

```text
treated: 185 {'PSID-1': 2490, 'CPS-1': 15992, 'CPS-3': 429}
```

### Score specifications and covariate lists

```python
# Table 3, notes (e) and (g), in term notation
PAPER_SCORE = {
 "PSID-1": ["age", "age^2", "education", "education^2", "married", "nodegree", "black", "hispanic",
            "re74", "re75", "re74^2", "re75^2", "black*u74"],
 "CPS-1":  ["age", "age^2", "age^3", "education", "education^2", "married", "nodegree", "black",
            "hispanic", "re74", "re75", "u74", "u75", "education*re74"],
}
PAPER_SCORE["CPS-3"] = PAPER_SCORE["CPS-1"]
SH_SCORE = {g: [x for x in xs if x != "nodegree"] for g, xs in PAPER_SCORE.items()}

REG_A = ["age", "age^2", "education", "nodegree", "black", "hispanic", "re74", "re75"]   # footnote (a)
BAL   = ["age", "education", "black", "hispanic", "nodegree", "married", "re74", "re75", "u74", "u75"]

def pool(comp, xs, t=None):
    """Stack treated + comparison units, create the needed terms, attach the score."""
    t = treated if t is None else t
    d = pd.concat([t, comp], ignore_index=True)
    d = add_terms(d, sorted(set(xs + REG_A)))
    d["ps"] = logit_ps(d, xs)
    return d

def support(d):
    pt = d.ps[d.treat == 1]
    return d[(d.treat == 1) | ((d.ps >= pt.min()) & (d.ps <= pt.max()))]
```

### Strata rules and estimators

- `strata_reproduction` and `strata_sh` are the two splitting rules described above.
- `t_stat` is the two-sample t-statistic used by both rules and by the balance tests.
- `stratify` gives columns (4) and (5).
- `match_nn` gives each comparison unit's weight: the number of times it is used as a nearest neighbor.

```python
def t_stat(g, v):
    a, b = g[g.treat == 1][v], g[g.treat == 0][v]
    if len(a) < 2 or len(b) < 2:
        return 0.0
    se = np.sqrt(a.var() / len(a) + b.var() / len(b))
    return 0.0 if se == 0 else (a.mean() - b.mean()) / se

def in_block(s, a, b):
    return s[(s.ps >= a) & ((s.ps < b) if b < 1 else (s.ps <= b))]

def merge_one_sided(s, blocks):
    """Merge a block with treated but no comparison units into its neighbor."""
    out = []
    for a, b in blocks:
        g = in_block(s, a, b)
        if out and (g.treat == 0).sum() == 0 and (g.treat == 1).sum() > 0:
            out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    while len(out) > 1 and (in_block(s, *out[0]).treat == 0).sum() == 0:
        out[:2] = [(out[0][0], out[1][1])]
    return out

def strata_reproduction(s, min_n=3):
    """Appendix rule: split when the score or any covariate differs; halves need >= min_n of each group."""
    lo, hi = s.ps.min(), 1.0
    edges = sorted({lo, hi, *[e for e in (0.2, 0.4, 0.6, 0.8) if lo < e < hi]})
    todo, done = list(zip(edges[:-1], edges[1:])), []
    while todo:
        a, b = todo.pop()
        g = in_block(s, a, b)
        bad = any(abs(t_stat(g, v)) > 1.96 for v in ["ps"] + BAL)
        if bad and (g.treat == 1).sum() >= 2 and (g.treat == 0).sum() >= 2:
            m = (a + b) / 2
            L, R = in_block(s, a, m), in_block(s, m, b)
            if min((L.treat == 1).sum(), (L.treat == 0).sum(), (R.treat == 1).sum(), (R.treat == 0).sum()) >= min_n:
                todo += [(a, m), (m, b)]
                continue
        done.append((a, b))
    return merge_one_sided(s, sorted(x for x in done if len(in_block(s, *x))))

def strata_sh(s):
    """SH's rule: 5 equal-width blocks on [0,1]; split while the mean score differs."""
    todo, done = [(i / 5, (i + 1) / 5) for i in range(5)], []
    while todo:
        a, b = todo.pop()
        g = in_block(s, a, b)
        if (g.treat == 1).sum() >= 2 and (g.treat == 0).sum() >= 2 and abs(t_stat(g, "ps")) > 1.96 and b - a > 1e-6:
            m = (a + b) / 2
            todo += [(a, m), (m, b)]
            continue
        done.append((a, b))
    return merge_one_sided(s, sorted(x for x in done if len(in_block(s, *x))))

def stratify(s, blocks, regx):
    n1 = (s.treat == 1).sum()
    e4 = v4 = e5 = v5 = 0.0
    for a, b in blocks:
        g = in_block(s, a, b)
        t, c = g[g.treat == 1], g[g.treat == 0]
        if len(t) == 0:
            continue
        w = len(t) / n1
        d4 = t.re78.mean() - c.re78.mean()
        var4 = np.nan_to_num(t.re78.var() / len(t)) + np.nan_to_num(c.re78.var() / len(c))
        if len(g) > len(regx) + 3:
            d5, se5 = ols_se(g, "re78", ["treat"] + regx)
        else:                                         # too few units for the regression
            d5, se5 = d4, np.sqrt(var4)
        e4 += w * d4; v4 += w * w * var4
        e5 += w * d5; v5 += w * w * se5 ** 2
    return (e4, np.sqrt(v4)), (e5, np.sqrt(v5))

def match_nn(d):
    t, c = d[d.treat == 1], d[d.treat == 0]
    pc = c.ps.to_numpy()
    w = np.zeros(len(pc))
    for p in t.ps:
        w[np.abs(pc - p).argmin()] += 1                # argmin: first unit in data order on ties
    return w

def matching(d, regx):
    t, c = d[d.treat == 1], d[d.treat == 0]
    w = match_nn(d)
    att = t.re78.mean() - (w @ c.re78.to_numpy()) / len(t)
    u = w > 0
    dm = pd.concat([t, c[u]])
    reg = ols_se(dm, "re78", ["treat"] + regx, np.r_[np.ones(len(t)), w[u]])[0]
    return att, reg, w

def boot_matching(comp, xs, regx, reps=REPS, seed=SEED):
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(reps):
        tb = treated.iloc[rng.integers(0, len(treated), len(treated))]
        cb = comp.iloc[rng.integers(0, len(comp), len(comp))]
        draws.append(matching(pool(cb, xs, tb), regx)[:2])
    return np.std(draws, axis=0, ddof=1)

def table3_row(g, xs, strata_rule, regx=REG_A):
    comp = GROUPS[g]
    d = pool(comp, xs)
    s = support(d)
    blocks = strata_rule(s)
    (e4, s4), (e5, s5) = stratify(s, blocks, regx)
    b3, se3 = ols_se(s.assign(ps2=s.ps ** 2), "re78", ["treat", "ps", "ps2"])
    att, reg, w = matching(d, regx)
    se = boot_matching(comp, xs, regx)
    row = {"(1)": te(*ols_se(d, "re78", ["treat"])), "(2)": te(*ols_se(d, "re78", ["treat"] + REG_A)),
           "(3)": te(b3, se3), "(4)": te(e4, s4), "(5)": te(e5, s5), "(6)": f"{len(s):,}",
           "(7)": te(att, se[0]), "(8)": te(reg, se[1])}
    return row, {"d": d, "s": s, "blocks": blocks, "w": w}
```

# Part 1. Tables 3 and 4

## Table 3: estimated training effects

For each group: **(i) Paper / (ii) Reproduction / (iii) SH**. The NSW row is the experimental benchmark.

```python
PAPER_T3 = {
 "PSID-1": ["-15,205 (1,154)", "731 (886)", "294 (1,389)", "1,608 (1,571)", "1,494 (1,581)", "1,255", "1,691 (2,209)", "1,473 (809)"],
 "CPS-1":  ["-8,498 (712)", "972 (550)", "1,117 (747)", "1,713 (1,115)", "1,774 (1,152)", "4,117", "1,582 (1,069)", "1,616 (751)"],
 "CPS-3":  ["-635 (657)", "1,326 (798)", "556 (951)", "1,252 (1,617)", "2,219 (2,082)", "514", "587 (1,496)", "662 (776)"],
}
COLS = ["(1)", "(2)", "(3)", "(4)", "(5)", "(6)", "(7)", "(8)"]
SETUPS = {"(ii) Reproduction": (PAPER_SCORE, strata_reproduction),
          "(iii) SH":          (SH_SCORE,    strata_sh)}

nsw_a = add_terms(nsw.copy(), REG_A)
rows = {("NSW", "(i) Paper"): ["1,794 (633)", "1,672 (638)"] + [""] * 6,
        ("NSW", "(ii) Reproduction"): [te(*ols_se(nsw_a, "re78", ["treat"])),
                                       te(*ols_se(nsw_a, "re78", ["treat"] + REG_A))] + [""] * 6}
results = {}
for g in GROUPS:
    rows[(g, "(i) Paper")] = PAPER_T3[g]
    for label, (scores, rule) in SETUPS.items():
        r, results[(g, label)] = table3_row(g, scores[g], rule)
        rows[(g, label)] = [r[c] for c in COLS]
table3 = pd.DataFrame(rows, index=COLS).T
table3
```

<div class="table-wrap">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>(1)</th>
      <th>(2)</th>
      <th>(3)</th>
      <th>(4)</th>
      <th>(5)</th>
      <th>(6)</th>
      <th>(7)</th>
      <th>(8)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">NSW</th>
      <th>(i) Paper</th>
      <td>1,794 (633)</td>
      <td>1,672 (638)</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>1,794 (633)</td>
      <td>1,672 (638)</td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
      <td></td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">PSID-1</th>
      <th>(i) Paper</th>
      <td>-15,205 (1,154)</td>
      <td>731 (886)</td>
      <td>294 (1,389)</td>
      <td>1,608 (1,571)</td>
      <td>1,494 (1,581)</td>
      <td>1,255</td>
      <td>1,691 (2,209)</td>
      <td>1,473 (809)</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>-15,205 (1,155)</td>
      <td>218 (866)</td>
      <td>265 (1,387)</td>
      <td>2,154 (790)</td>
      <td>2,492 (1,894)</td>
      <td>1,331</td>
      <td>552 (1,924)</td>
      <td>495 (1,590)</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>-15,205 (1,155)</td>
      <td>218 (866)</td>
      <td>286 (1,388)</td>
      <td>1,039 (918)</td>
      <td>1,072 (1,336)</td>
      <td>1,342</td>
      <td>1,691 (1,997)</td>
      <td>1,674 (1,316)</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">CPS-1</th>
      <th>(i) Paper</th>
      <td>-8,498 (712)</td>
      <td>972 (550)</td>
      <td>1,117 (747)</td>
      <td>1,713 (1,115)</td>
      <td>1,774 (1,152)</td>
      <td>4,117</td>
      <td>1,582 (1,069)</td>
      <td>1,616 (751)</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>-8,498 (712)</td>
      <td>739 (547)</td>
      <td>1,071 (762)</td>
      <td>1,564 (733)</td>
      <td>1,781 (1,050)</td>
      <td>4,041</td>
      <td>952 (925)</td>
      <td>1,098 (916)</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>-8,498 (712)</td>
      <td>739 (547)</td>
      <td>1,109 (747)</td>
      <td>1,833 (675)</td>
      <td>1,877 (1,033)</td>
      <td>3,566</td>
      <td>1,582 (909)</td>
      <td>1,565 (884)</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">CPS-3</th>
      <th>(i) Paper</th>
      <td>-635 (657)</td>
      <td>1,326 (798)</td>
      <td>556 (951)</td>
      <td>1,252 (1,617)</td>
      <td>2,219 (2,082)</td>
      <td>514</td>
      <td>587 (1,496)</td>
      <td>662 (776)</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>-635 (657)</td>
      <td>1,326 (796)</td>
      <td>625 (950)</td>
      <td>1,587 (884)</td>
      <td>1,434 (1,869)</td>
      <td>506</td>
      <td>-467 (1,182)</td>
      <td>-419 (1,100)</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>-635 (657)</td>
      <td>1,326 (796)</td>
      <td>649 (950)</td>
      <td>1,211 (936)</td>
      <td>1,541 (1,655)</td>
      <td>515</td>
      <td>587 (1,152)</td>
      <td>777 (1,118)</td>
    </tr>
  </tbody>
</table>
</div>

**Reading Table 3:**

- **Matching, column (7):** SH reproduces the paper exactly for all three groups. Reproduction does not.
- **Column (8):** under SH it is close to the paper, and exact only with the linear covariate set (see above).
- **Stratification, columns (3)–(6):** neither version reproduces the paper. The column (6) sample differs, and the paper's strata are not published.
- **Column (2):** the paper's PSID-1 and CPS-1 values (731, 972) do not follow its own footnote (a), which gives 218 and 739.

## Table 4: matched comparison samples

Means of the matched comparison units, **weighted by the number of times each unit is matched**; that is how the paper's means reproduce. `N` = distinct matched comparison units.

```python
T4 = ["age", "education", "black", "hispanic", "nodegree", "married", "re74", "re75"]
T4_LAB = ["Age", "Education", "Black", "Hispanic", "No degree", "Married", "RE74", "RE75"]
PAPER_T4 = {   # means only, transcribed from Table 4
 "PSID-1": [56, 26.39, 10.62, .86, .02, .55, .15, 1794, 1126],
 "CPS-1":  [119, 26.91, 10.52, .86, .04, .64, .19, 2110, 1396],
 "CPS-3":  [63, 25.94, 10.69, .87, .06, .53, .13, 2709, 1587],
}

def fmt(v, x):
    return f"{x:,.0f}" if v.startswith("re") else f"{x:.2f}"

def matched_means(info):
    d, w = info["d"], info["w"]
    c = d[d.treat == 0]
    u = w > 0
    return [int(u.sum())] + [np.average(c[v][u], weights=w[u]) for v in T4]

rows4 = {("NSW treated", ""): [185] + [fmt(v, treated[v].mean()) for v in T4]}
for g in GROUPS:
    p = PAPER_T4[g]
    rows4[("M" + g, "(i) Paper")] = [p[0]] + [fmt(v, x) for v, x in zip(T4, p[1:])]
    for label in SETUPS:
        m = matched_means(results[(g, label)])
        rows4[("M" + g, label)] = [m[0]] + [fmt(v, x) for v, x in zip(T4, m[1:])]
table4 = pd.DataFrame(rows4, index=["N"] + T4_LAB).T
table4
```

<div class="table-wrap">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>N</th>
      <th>Age</th>
      <th>Education</th>
      <th>Black</th>
      <th>Hispanic</th>
      <th>No degree</th>
      <th>Married</th>
      <th>RE74</th>
      <th>RE75</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>NSW treated</th>
      <th></th>
      <td>185</td>
      <td>25.82</td>
      <td>10.35</td>
      <td>0.84</td>
      <td>0.06</td>
      <td>0.71</td>
      <td>0.19</td>
      <td>2,096</td>
      <td>1,532</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">MPSID-1</th>
      <th>(i) Paper</th>
      <td>56</td>
      <td>26.39</td>
      <td>10.62</td>
      <td>0.86</td>
      <td>0.02</td>
      <td>0.55</td>
      <td>0.15</td>
      <td>1,794</td>
      <td>1,126</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>55</td>
      <td>26.13</td>
      <td>10.31</td>
      <td>0.92</td>
      <td>0.01</td>
      <td>0.56</td>
      <td>0.16</td>
      <td>1,873</td>
      <td>1,528</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>52</td>
      <td>26.39</td>
      <td>10.63</td>
      <td>0.86</td>
      <td>0.02</td>
      <td>0.55</td>
      <td>0.15</td>
      <td>1,794</td>
      <td>1,126</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">MCPS-1</th>
      <th>(i) Paper</th>
      <td>119</td>
      <td>26.91</td>
      <td>10.52</td>
      <td>0.86</td>
      <td>0.04</td>
      <td>0.64</td>
      <td>0.19</td>
      <td>2,110</td>
      <td>1,396</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>106</td>
      <td>25.66</td>
      <td>10.80</td>
      <td>0.88</td>
      <td>0.03</td>
      <td>0.59</td>
      <td>0.18</td>
      <td>1,839</td>
      <td>1,475</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>112</td>
      <td>26.91</td>
      <td>10.52</td>
      <td>0.86</td>
      <td>0.04</td>
      <td>0.64</td>
      <td>0.19</td>
      <td>2,110</td>
      <td>1,396</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">MCPS-3</th>
      <th>(i) Paper</th>
      <td>63</td>
      <td>25.94</td>
      <td>10.69</td>
      <td>0.87</td>
      <td>0.06</td>
      <td>0.53</td>
      <td>0.13</td>
      <td>2,709</td>
      <td>1,587</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>63</td>
      <td>25.49</td>
      <td>10.50</td>
      <td>0.89</td>
      <td>0.05</td>
      <td>0.68</td>
      <td>0.16</td>
      <td>1,771</td>
      <td>1,666</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>64</td>
      <td>25.95</td>
      <td>10.69</td>
      <td>0.87</td>
      <td>0.06</td>
      <td>0.53</td>
      <td>0.13</td>
      <td>2,709</td>
      <td>1,587</td>
    </tr>
  </tbody>
</table>
</div>

Under SH, every mean equals the paper's to rounding. The paper's `N` (56 / 119 / 63) does not match the number of distinct matched units (52 / 112 / 64).

### Covariate balance for matching

For each covariate, the standardized difference between the NSW treated and the matched comparison sample:

$$d = 100\cdot\frac{\bar x_T-\bar x_{C,\text{matched}}}{\sqrt{(s_T^2+s_C^2)/2}}$$

Here $s_T^2$ and $s_C^2$ are the variances in the full treated and full comparison samples. **Balanced** means $|d|<10$ for every covariate.

- **(i) Paper:** uses the matched means printed in Table 4. Only the 8 covariates in that table are available; u74 and u75 are not reported.
- **(ii) Reproduction and (iii) SH:** use our matched samples, for all 10 covariates.

```python
def std_diff(xbar_c, v, comp):
    return 100 * (treated[v].mean() - xbar_c) / np.sqrt((treated[v].var() + comp[v].var()) / 2)

def matching_balance(g, label):
    comp = GROUPS[g]
    if label == "(i) Paper":
        return {v: std_diff(x, v, comp) for v, x in zip(T4, PAPER_T4[g][1:])}
    info = results[(g, label)]
    d, w = info["d"], info["w"]
    c = d[d.treat == 0]
    u = w > 0
    return {v: std_diff(np.average(c[v][u], weights=w[u]), v, comp) for v in BAL}

LABELS = ["(i) Paper", "(ii) Reproduction", "(iii) SH"]
mbal = pd.DataFrame({(g, lab): matching_balance(g, lab) for g in GROUPS for lab in LABELS}).reindex(BAL)
print("Standardized differences after matching (%); |d| >= 10 is flagged in the summary below")
mbal.round(1)
```

```text
Standardized differences after matching (%); |d| >= 10 is flagged in the summary below
```

<div class="table-wrap">
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th colspan="3" halign="left">PSID-1</th>
      <th colspan="3" halign="left">CPS-1</th>
      <th colspan="3" halign="left">CPS-3</th>
    </tr>
    <tr>
      <th></th>
      <th>(i) Paper</th>
      <th>(ii) Reproduction</th>
      <th>(iii) SH</th>
      <th>(i) Paper</th>
      <th>(ii) Reproduction</th>
      <th>(iii) SH</th>
      <th>(i) Paper</th>
      <th>(ii) Reproduction</th>
      <th>(iii) SH</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>age</th>
      <td>-6.4</td>
      <td>-3.5</td>
      <td>-6.4</td>
      <td>-11.8</td>
      <td>1.7</td>
      <td>-11.7</td>
      <td>-1.4</td>
      <td>3.6</td>
      <td>-1.4</td>
    </tr>
    <tr>
      <th>education</th>
      <td>-10.5</td>
      <td>1.5</td>
      <td>-10.8</td>
      <td>-7.0</td>
      <td>-18.3</td>
      <td>-7.0</td>
      <td>-13.9</td>
      <td>-6.3</td>
      <td>-14.0</td>
    </tr>
    <tr>
      <th>black</th>
      <td>-4.2</td>
      <td>-20.2</td>
      <td>-5.4</td>
      <td>-5.3</td>
      <td>-10.2</td>
      <td>-6.8</td>
      <td>-7.0</td>
      <td>-11.3</td>
      <td>-7.0</td>
    </tr>
    <tr>
      <th>hispanic</th>
      <td>18.8</td>
      <td>25.8</td>
      <td>20.6</td>
      <td>7.8</td>
      <td>13.1</td>
      <td>8.7</td>
      <td>-0.2</td>
      <td>1.8</td>
      <td>-1.8</td>
    </tr>
    <tr>
      <th>nodegree</th>
      <td>34.5</td>
      <td>31.8</td>
      <td>35.4</td>
      <td>14.9</td>
      <td>24.9</td>
      <td>14.2</td>
      <td>37.6</td>
      <td>5.7</td>
      <td>37.6</td>
    </tr>
    <tr>
      <th>married</th>
      <td>10.7</td>
      <td>7.4</td>
      <td>10.3</td>
      <td>-0.2</td>
      <td>1.3</td>
      <td>-1.3</td>
      <td>13.2</td>
      <td>6.0</td>
      <td>13.2</td>
    </tr>
    <tr>
      <th>re74</th>
      <td>3.0</td>
      <td>2.2</td>
      <td>3.0</td>
      <td>-0.2</td>
      <td>3.4</td>
      <td>-0.2</td>
      <td>-10.4</td>
      <td>5.5</td>
      <td>-10.4</td>
    </tr>
    <tr>
      <th>re75</th>
      <td>4.1</td>
      <td>0.0</td>
      <td>4.1</td>
      <td>2.0</td>
      <td>0.8</td>
      <td>2.0</td>
      <td>-1.7</td>
      <td>-4.1</td>
      <td>-1.7</td>
    </tr>
    <tr>
      <th>u74</th>
      <td>NaN</td>
      <td>21.4</td>
      <td>14.3</td>
      <td>NaN</td>
      <td>-15.0</td>
      <td>-5.5</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>10.9</td>
    </tr>
    <tr>
      <th>u75</th>
      <td>NaN</td>
      <td>1.3</td>
      <td>-18.6</td>
      <td>NaN</td>
      <td>0.0</td>
      <td>-1.3</td>
      <td>NaN</td>
      <td>-14.7</td>
      <td>4.5</td>
    </tr>
  </tbody>
</table>
</div>

```python
summary = {}
for (g, lab), col in mbal.items():
    col = col.dropna()
    bad = col[col.abs() >= 10]
    summary[(g, lab)] = {"max |d|": round(col.abs().max(), 1),
                         "covariates with |d| >= 10": ", ".join(f"{k} ({v:.0f})" for k, v in bad.items()) or "none",
                         "balanced?": "yes" if bad.empty else "no"}
pd.DataFrame(summary).T
```

<div class="table-wrap">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>max |d|</th>
      <th>covariates with |d| &gt;= 10</th>
      <th>balanced?</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="3" valign="top">PSID-1</th>
      <th>(i) Paper</th>
      <td>34.5</td>
      <td>education (-11), hispanic (19), nodegree (35), married (11)</td>
      <td>no</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>31.8</td>
      <td>black (-20), hispanic (26), nodegree (32), u74 (21)</td>
      <td>no</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>35.4</td>
      <td>education (-11), hispanic (21), nodegree (35), married (10), u74 (14), u75 (-19)</td>
      <td>no</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">CPS-1</th>
      <th>(i) Paper</th>
      <td>14.9</td>
      <td>age (-12), nodegree (15)</td>
      <td>no</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>24.9</td>
      <td>education (-18), black (-10), hispanic (13), nodegree (25), u74 (-15)</td>
      <td>no</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>14.2</td>
      <td>age (-12), nodegree (14)</td>
      <td>no</td>
    </tr>
    <tr>
      <th rowspan="3" valign="top">CPS-3</th>
      <th>(i) Paper</th>
      <td>37.6</td>
      <td>education (-14), nodegree (38), married (13), re74 (-10)</td>
      <td>no</td>
    </tr>
    <tr>
      <th>(ii) Reproduction</th>
      <td>14.7</td>
      <td>black (-11), u75 (-15)</td>
      <td>no</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>37.6</td>
      <td>education (-14), nodegree (38), married (13), re74 (-10), u74 (11)</td>
      <td>no</td>
    </tr>
  </tbody>
</table>
</div>

**Reading the matching balance:**

- **No version is balanced** for any group, including the paper's own matched samples (computed from its Table 4 means).
- **`nodegree` is the most frequent problem.** It exceeds 10% for PSID-1 and CPS-1 in all three versions, and for CPS-3 under Paper and SH (38%).
- **CPS-1 under Paper and SH:** only `nodegree` (15% and 14%) and age (−12%) fail.
- **PSID-1:** Hispanic fails in all three versions, and u74 in both reproductions. Reproduction also fails on black; Paper and SH also fail on education and married, and SH on u75.
- **CPS-3:** Paper and SH fail on `nodegree`, education, married and RE74 (SH also on u74). Reproduction fails only on black and u75.

# Part 2. Covariate balance for stratification

The Appendix test: within each stratum that has at least 2 treated and 2 comparison units, a two-sample t-test of each covariate's mean, treated vs. comparison. A covariate **fails** in a stratum if $|t|>1.96$.

- **Strictly balanced:** no test fails.
- **Chance level:** with many tests, about 5% fail by chance even when every stratum is balanced. We report that benchmark next to the count.

The paper publishes no strata, so only (ii) Reproduction and (iii) SH can be checked. No adjustment is made here; the goal is to see which covariates cause the problems.

```python
def strata_table(info):
    s, rows = info["s"], []
    for a, b in info["blocks"]:
        blk = in_block(s, a, b)
        n1, n0 = int((blk.treat == 1).sum()), int((blk.treat == 0).sum())
        tested = n1 >= 2 and n0 >= 2
        bad = [v for v in BAL if tested and abs(t_stat(blk, v)) > 1.96]
        rows.append({"score range": f"[{a:.3f}, {b:.3f})", "treated": n1, "comparison": n0,
                     "tested": "yes" if tested else "no (too few units)",
                     "failing covariates": ", ".join(f"{v} (t={t_stat(blk, v):.1f})" for v in bad) or "-"})
    return pd.DataFrame(rows)

def strata_counts(info):
    s, per_var, tests, n_tested = info["s"], {v: 0 for v in BAL}, 0, 0
    for a, b in info["blocks"]:
        blk = in_block(s, a, b)
        if (blk.treat == 1).sum() < 2 or (blk.treat == 0).sum() < 2:
            continue
        n_tested += 1
        for v in BAL:
            tests += 1
            per_var[v] += abs(t_stat(blk, v)) > 1.96
    return per_var, tests, n_tested

sbal, overview = {}, {}
for g in GROUPS:
    for lab in SETUPS:
        info = results[(g, lab)]
        per_var, tests, n_tested = strata_counts(info)
        sbal[(g, lab)] = per_var
        f = sum(per_var.values())
        overview[(g, lab)] = {"strata": len(info["blocks"]), "strata tested": n_tested,
                              "failing tests / tests": f"{f} / {tests}",
                              "expected by chance (5%)": round(0.05 * tests, 1),
                              "result": "strictly balanced" if f == 0 else
                                        ("within chance level" if f <= 0.05 * tests else "NOT balanced")}
pd.DataFrame(overview).T
```

<div class="table-wrap">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>strata</th>
      <th>strata tested</th>
      <th>failing tests / tests</th>
      <th>expected by chance (5%)</th>
      <th>result</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">PSID-1</th>
      <th>(ii) Reproduction</th>
      <td>8</td>
      <td>8</td>
      <td>11 / 80</td>
      <td>4.0</td>
      <td>NOT balanced</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>8</td>
      <td>8</td>
      <td>10 / 80</td>
      <td>4.0</td>
      <td>NOT balanced</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">CPS-1</th>
      <th>(ii) Reproduction</th>
      <td>18</td>
      <td>18</td>
      <td>27 / 180</td>
      <td>9.0</td>
      <td>NOT balanced</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>9</td>
      <td>9</td>
      <td>8 / 90</td>
      <td>4.5</td>
      <td>NOT balanced</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">CPS-3</th>
      <th>(ii) Reproduction</th>
      <td>10</td>
      <td>10</td>
      <td>8 / 100</td>
      <td>5.0</td>
      <td>NOT balanced</td>
    </tr>
    <tr>
      <th>(iii) SH</th>
      <td>5</td>
      <td>5</td>
      <td>6 / 50</td>
      <td>2.5</td>
      <td>NOT balanced</td>
    </tr>
  </tbody>
</table>
</div>

### Which covariates cause the problems

Number of strata in which each covariate fails the t-test (out of the strata tested):

```python
problems = pd.DataFrame(sbal).reindex(BAL).astype(int)
problems.loc["total"] = problems.sum()
problems
```

<div class="table-wrap">
<table border="1" class="dataframe">
  <thead>
    <tr>
      <th></th>
      <th colspan="2" halign="left">PSID-1</th>
      <th colspan="2" halign="left">CPS-1</th>
      <th colspan="2" halign="left">CPS-3</th>
    </tr>
    <tr>
      <th></th>
      <th>(ii) Reproduction</th>
      <th>(iii) SH</th>
      <th>(ii) Reproduction</th>
      <th>(iii) SH</th>
      <th>(ii) Reproduction</th>
      <th>(iii) SH</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>age</th>
      <td>0</td>
      <td>0</td>
      <td>3</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>education</th>
      <td>0</td>
      <td>0</td>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>1</td>
    </tr>
    <tr>
      <th>black</th>
      <td>3</td>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>2</td>
      <td>1</td>
    </tr>
    <tr>
      <th>hispanic</th>
      <td>3</td>
      <td>2</td>
      <td>3</td>
      <td>1</td>
      <td>1</td>
      <td>2</td>
    </tr>
    <tr>
      <th>nodegree</th>
      <td>0</td>
      <td>0</td>
      <td>4</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>married</th>
      <td>2</td>
      <td>1</td>
      <td>4</td>
      <td>1</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <th>re74</th>
      <td>1</td>
      <td>2</td>
      <td>2</td>
      <td>1</td>
      <td>3</td>
      <td>1</td>
    </tr>
    <tr>
      <th>re75</th>
      <td>0</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>u74</th>
      <td>1</td>
      <td>2</td>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
    </tr>
    <tr>
      <th>u75</th>
      <td>1</td>
      <td>1</td>
      <td>3</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <th>total</th>
      <td>11</td>
      <td>10</td>
      <td>27</td>
      <td>8</td>
      <td>8</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>

### Stratum by stratum

```python
for g in GROUPS:
    for lab in SETUPS:
        print(f"{g} — {lab}")
        print(strata_table(results[(g, lab)]).to_string(index=False), "\n")
```

```text
PSID-1 — (ii) Reproduction
   score range  treated  comparison tested                                                               failing covariates
[0.001, 0.050)        7         918    yes black (t=3.5), hispanic (t=-6.8), married (t=12.7), u74 (t=-11.4), u75 (t=-13.2)
[0.050, 0.100)        4         100    yes                                                                hispanic (t=-2.5)
[0.100, 0.200)        7          53    yes                                                                    re74 (t=-3.5)
[0.200, 0.400)       27          41    yes                                                                                -
[0.400, 0.500)        5           8    yes                                                                   black (t=-2.4)
[0.500, 0.600)       15           6    yes                                                                                -
[0.600, 0.800)       22          13    yes                                                                                -
[0.800, 1.000)       98           7    yes                                black (t=-2.3), hispanic (t=2.3), married (t=2.7) 

PSID-1 — (iii) SH
   score range  treated  comparison tested                                                               failing covariates
[0.000, 0.050)        7         924    yes black (t=3.5), hispanic (t=-6.8), married (t=12.7), u74 (t=-11.2), u75 (t=-13.1)
[0.050, 0.100)        4         102    yes                                                                hispanic (t=-2.5)
[0.100, 0.200)        7          56    yes                                                                    re74 (t=-3.7)
[0.200, 0.300)       10          25    yes                                                                                -
[0.300, 0.400)       18          16    yes                                                                                -
[0.400, 0.600)       21          14    yes                                                                   black (t=-2.4)
[0.600, 0.950)       58          13    yes                                                       re74 (t=-2.5), u74 (t=3.7)
[0.950, 1.000)       60           7    yes                                                                                - 

CPS-1 — (ii) Reproduction
   score range  treated  comparison tested                                                                            failing covariates
[0.001, 0.026)       18        3221    yes                                                                                             -
[0.026, 0.051)        8         242    yes                                                                                  age (t=-2.0)
[0.051, 0.063)        3          64    yes                                  black (t=5.4), hispanic (t=-3.8), u74 (t=-4.8), u75 (t=-5.0)
[0.063, 0.076)        4          44    yes hispanic (t=-2.6), married (t=-6.6), re74 (t=-6.3), re75 (t=-6.4), u74 (t=10.1), u75 (t=10.1)
[0.076, 0.088)        4          35    yes                                             age (t=-3.5), hispanic (t=-2.4), married (t=-2.7)
[0.088, 0.101)        3          21    yes                                                                                 black (t=2.5)
[0.101, 0.125)        5          35    yes                                              nodegree (t=4.5), married (t=-2.9), u75 (t=-3.2)
[0.125, 0.138)        3          11    yes                                         education (t=-2.3), nodegree (t=2.9), married (t=5.2)
[0.138, 0.150)        3          16    yes                                                  nodegree (t=5.0), re74 (t=-2.3), u74 (t=3.4)
[0.150, 0.200)        8          33    yes                                                                                             -
[0.200, 0.225)        4          22    yes                                                          education (t=3.0), nodegree (t=-2.2)
[0.225, 0.250)        3          12    yes                                                                                  age (t=-2.4)
[0.250, 0.300)        6          20    yes                                                                                             -
[0.300, 0.400)       19          32    yes                                                                                             -
[0.400, 0.500)       16          14    yes                                                                                             -
[0.500, 0.600)       22          17    yes                                                                                             -
[0.600, 0.800)       20           9    yes                                                                                             -
[0.800, 1.000)       36           8    yes                                                                                             - 

CPS-1 — (iii) SH
   score range  treated  comparison tested                                                   failing covariates
[0.000, 0.025)       18        2738    yes                                                                    -
[0.025, 0.050)        8         231    yes                                                         age (t=-2.2)
[0.050, 0.100)       15         183    yes hispanic (t=-5.6), nodegree (t=2.1), married (t=-2.1), re74 (t=-3.6)
[0.100, 0.200)       18          99    yes                                                   education (t=-2.1)
[0.200, 0.400)       33          78    yes                                                    education (t=2.2)
[0.400, 0.500)       19          24    yes                                                                    -
[0.500, 0.600)       20           8    yes                                                          u75 (t=2.4)
[0.600, 0.800)       20          13    yes                                                                    -
[0.800, 1.000)       34           7    yes                                                                    - 

CPS-3 — (ii) Reproduction
   score range  treated  comparison tested                          failing covariates
[0.019, 0.110)        9         194    yes               black (t=-2.7), re74 (t=-2.5)
[0.110, 0.200)        4          38    yes                                           -
[0.200, 0.250)        5           8    yes                           hispanic (t=-2.0)
[0.250, 0.300)        5          16    yes                               re74 (t=-2.5)
[0.300, 0.400)        6          16    yes                               black (t=2.6)
[0.400, 0.600)       30          17    yes                                           -
[0.600, 0.700)       21          11    yes                                           -
[0.700, 0.800)       22          14    yes                                           -
[0.800, 0.900)       35           3    yes married (t=4.2), re74 (t=2.0), u74 (t=-2.9)
[0.900, 1.000)       48           4    yes                                           - 

CPS-3 — (iii) SH
   score range  treated  comparison tested               failing covariates
[0.000, 0.200)       13         240    yes                education (t=2.2)
[0.200, 0.400)       19          42    yes hispanic (t=-2.9), re74 (t=-2.9)
[0.400, 0.600)       26          15    yes                                -
[0.600, 0.800)       45          25    yes black (t=-2.3), hispanic (t=2.1)
[0.800, 1.000)       82           8    yes                  married (t=3.7)
```

**Reading the stratification balance:**

- **No group is balanced under either version.** Failing tests are about 1.6–3 times what chance alone would produce.
- **Covariates that fail most often:** Hispanic, black, married and RE74, in every group. `nodegree`, age and u74/u75 fail repeatedly under Reproduction for CPS-1. u74 and u75 also fail in PSID-1's lowest stratum.
- **Where the failures concentrate:**
  - *The lowest-score strata,* where a few treated units face hundreds or thousands of comparison units. PSID-1's first stratum (7 treated vs. about 920 comparison units) fails on five covariates, and most CPS-1 failures are in strata below a score of 0.25.
  - *The highest-score strata,* where many treated units share a handful of comparison units. For example, PSID-1 above 0.6 under SH has 118 treated vs. 20 comparison units, and CPS-3 above 0.8 has 83 treated vs. 7 comparison units under Reproduction.
- **Next step, if wanted:** the Appendix's remedy is to add interaction or higher-order terms for these covariates and re-estimate the score. That is not done here.
