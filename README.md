# Hall & Jones (1999) — Replication

Replication of **Hall, Robert E. and Charles I. Jones (1999)**. "Why Do Some Countries Produce So Much More Output Per Worker Than Others?" _Quarterly Journal of Economics_ 114(1): 83–116.

This repository contains the full data pipeline and analysis code to replicate the paper's core empirical results: the levels accounting decomposition (Table I), the IV regression of output per worker on social infrastructure (Table II), and the associated scatter plots (Figures I and II).

---

## Results at a Glance

| Result                            | Hall & Jones (1999) | This Replication |
| --------------------------------- | ------------------- | ---------------- |
| TFP share of variance in log(Y/L) | ~60%                | 65%              |
| OLS coefficient on S              | 3.29 (0.21)         | 3.18 (0.25)      |
| 2SLS coefficient on S             | 5.14 (0.51)         | 4.42 (0.43)      |
| First-stage F-statistic           | —                   | 16.2             |
| Overidentification test p-value   | 0.256               | 0.089            |
| Sample size (with imputation)     | 127                 | 113              |

The direction, magnitude, and interpretation of all results are consistent with the original paper. Differences in coefficient magnitudes are explained by the documented deviations below — primarily the substitution of WGI (1996) for ICRG (1986–1995) as the governance measure.

---

## What This Replication Covers

- **Table I** — Levels accounting decomposition: Y/L = (K/Y)^(α/1−α) × h × A, with variance decomposition
- **Table II** — OLS and 2SLS regressions of log(Y/L) on social infrastructure S, with robust standard errors and instrument sensitivity analysis
- **Figure I** — Scatter plot of output per worker vs. social infrastructure
- **Figure II** — Scatter plot of TFP vs. social infrastructure
- **Two samples throughout**: complete cases only (N=89) and imputed sample (N=113)

---

## Documented Deviations from the Original Paper

| Item               | Hall & Jones (1999)      | This Replication             | Expected Effect                                                               |
| ------------------ | ------------------------ | ---------------------------- | ----------------------------------------------------------------------------- |
| Governance data    | ICRG GADP, avg 1986–1995 | WGI (rl + ge + cc avg), 1996 | Compresses governance scores for worst-governed countries; reduces β slightly |
| Price base year    | PWT 5.6 (1985 prices)    | PWT 10.01 (2017 prices)      | Changes absolute Y/L levels; relative rankings robust                         |
| Openness window    | Sachs-Warner 1950–1994   | Sachs-Warner 1950–1992       | Negligible — 2 of 40+ years                                                   |
| Mining value added | 1988 values              | 1990 values (proxy)          | Negligible — mining shares stable year-to-year                                |
| Standard errors    | Bootstrap (10,000 reps)  | HC1 robust                   | Minor effect on SEs; coefficients identical                                   |
| Sample size        | N=127 (imputed)          | N=113 (imputed)              | Some small/transition economies missing                                       |

---

## Repository Structure

```
hall-jones-1999-replication/
│
├── README.md                   # This file
├── requirements.txt            # Python dependencies
│
├── data_/
│   └── README_data.md          # Instructions for downloading all data files
│
├── src/
│   ├── data_pipeline.ipynb             # Step 1: data pipeline — merges all sources into merged.csv
│   └── analysis.ipynb         # Step 2: analysis — levels accounting, IV regression, figures
│
└── outputs/
    ├── merged.csv           # Merged dataset (183 countries, 23 variables)
    ├── figure_I_yl_vs_si.png   # Figure I: Y/L vs social infrastructure
    └── figure_II_tfp_vs_si.png # Figure II: TFP vs social infrastructure
```

---

## Data Sources

### 1. Penn World Tables 10.01

- **File name:** `pwt1001.xlsx`
- **Download:** https://www.rug.nl/ggdc/productivity/pwt/
- **What it provides:** Output per worker (Y/L), employment, investment share, capital stock components, human capital index — all years 1950–2019

### 2. Barro-Lee Education Dataset (2013)

- **File name:** `BL2013_MF2599_v2_2.csv`
- **Download:** http://www.barrolee.com/data/full1.htm
- **Select:** "Average years of schooling, age 25+, both sexes"
- **What it provides:** Average years of schooling by country and year (5-year intervals, 1950–2010)

### 3. World Bank Worldwide Governance Indicators (WGI)

- **File name:** `wgidataset_with_sourcedata-2025.xlsx`
- **Download:** https://www.worldbank.org/en/publication/worldwide-governance-indicators
- **What it provides:** Six governance dimensions (va, pv, ge, rq, rl, cc) for 200+ countries, 1996–present. We use rl, ge, cc for 1996 as a substitute for Hall & Jones's GADP index.

### 4. Sachs-Warner Openness Data

- **File name:** `open.csv`
- **Source:** Sachs, Jeffrey D. and Andrew Warner (1995). "Economic Reform and the Process of Global Integration." _Brookings Papers on Economic Activity_ 1: 1–118.
- **What it provides:** Annual binary openness indicator (0/1) for 152 countries, 1950–1992

### 5. Mining Value Added

- **File name:** `Contribution_of_mining_to_value_added.xlsx`
- **Download:** UNIDO / OECD national accounts — search "contribution of mining to value added by country"
- **What it provides:** Mining value added as % of GDP by country, 1990–present. We use 1990 as proxy for 1988.

### 6. CEPII GeoDist (Latitude instrument)

- **File name:** `geo_cepii.xls`
- **Download:** http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=6
- **What it provides:** Country-level geographic data including latitude, used to construct distance from equator instrument

### 7. Frankel-Romer Predicted Trade Share

- **File name:** `frankel_romer_trade.csv` _(generated by hj_merge.py from Appendix Table A1)_
- **Source:** Frankel, Jeffrey A. and David Romer (1999). "Does Trade Cause Growth?" _American Economic Review_ 89(3): 379–399.
- **What it provides:** Geography-predicted trade share from gravity model — instrument for social infrastructure

### 8. CEPII Bilateral Language Data (Language instruments)

- **File name:** `ling_web.csv`
- **Download:** http://www.cepii.fr/CEPII/en/bdd_modele/bdd_modele_item.asp?id=5
- **What it provides:** Bilateral common language variables, used to construct English-speaking and Western European language fraction instruments

---

## How to Run

### Requirements

```
Python 3.8+
numpy
pandas
scipy
matplotlib
```

Install dependencies:

```bash
pip install numpy pandas scipy matplotlib
```

### Step 1 — Download data

Follow the instructions in `data_/README_data.md` to download all raw data files and place them in the `data_/` folder.

### Step 2 — Run the data pipeline

```bash
python src/data_pipeline.py
```

This script merges all data sources into a single master dataset `outputs/merged.csv`. It will print a step-by-step log showing how many countries are matched at each stage. Expected output:

```
>>> STEP 1: Loading PWT 10.01...
  Countries in 1988: 183
  Countries with Y/L: 149
>>> STEP 2: Constructing capital stocks...
  Capital stocks constructed for: 183 countries
...
>>> PIPELINE COMPLETE. Master dataset saved to: outputs/hj_master.csv
```

### Step 3 — Run the analysis

```bash
python src/analysis.py
```

This script reads `merged.csv` and produces all tables and figures. Results are printed to the console and figures are saved to `outputs/`. Expected output includes Table I (levels accounting), Table II (OLS and 2SLS), and the two scatter plots.

---

## Methodology Notes

### Levels Accounting

Following Hall & Jones equation (1), output per worker is decomposed as:

```
Y/L = (K/Y)^(α/1−α) × h × A
```

with α = 1/3 (capital share). Capital stock K is constructed via the perpetual inventory method with δ = 6% depreciation, using investment data back to 1960. Human capital h is computed from Barro-Lee years of schooling using piecewise Mincerian returns (Psacharopoulos 1994): 13.4% for years 1–4, 10.1% for years 5–8, 6.8% for years 9+. TFP residual A is computed as the remainder.

Mining value added is subtracted from GDP before the decomposition to remove the effect of natural resource windfalls.

### Social Infrastructure Index

```
S = (governance + openness) / 2
```

Both components are normalised to [0, 1]. Governance is the average of WGI Rule of Law, Government Effectiveness, and Control of Corruption scores for 1996, divided by 100. Openness is the fraction of years each country was classified as open under the Sachs-Warner criteria over the available data window.

### IV Regression

```
log(Y/L) = α + β·S + ε
```

S is treated as endogenous due to reverse causality (richer countries can afford better institutions) and measurement error (GADP/openness are imperfect proxies). Four instruments are used: distance from equator (|latitude|/90), Frankel-Romer predicted trade share, fraction speaking English natively, and fraction speaking a Western European language natively. These capture historical Western European institutional influence and are plausibly exogenous to current productivity.

Standard errors are HC1 heteroskedasticity-robust. The overidentification test (Sargan statistic) tests whether the instruments are jointly valid.

### Imputation

For countries missing openness data, social infrastructure is imputed by regressing observed S on the four instruments plus quadratic distance from equator (OLS), then using fitted values for missing observations. This follows Hall & Jones's methodology. Results are reported for both the complete-case sample (N=89) and the imputed sample (N=113).

---

## Citation

If you use this replication code, please cite both the original paper and this repository:

**Original paper:**

> Hall, Robert E. and Charles I. Jones (1999). "Why Do Some Countries Produce So Much More Output Per Worker Than Others?" _Quarterly Journal of Economics_ 114(1): 83–116.

**This replication:**

> Mohammed Adams (2026). "Replication of Hall & Jones (1999)." GitHub repository: https://github.com/moadams847/hall-jones-1999-replication

---

## License

Code in this repository is released under the MIT License. Raw data files are subject to the terms of their respective sources — see the Data Sources section above.
