# Replication of Hall and Jones (1999)

**Why Do Some Countries Produce So Much More Output Per Worker Than Others?**  
A Replication and Extension of Hall and Jones (1999), *Quarterly Journal of Economics* 114(1): 83–116

**Group 16 | Design, Execution and Evaluation of Research in Economics and Finance | Summer 2026**  
ABDULAI Yussif · ADAMS Mohammed · ADEOLA Ayankunle Babatunde · KHAYANKHYARVAA Gerelmaa

---

## Repository structure

```
├── 01_data_pipeline.ipynb     # Build master datasets from raw sources
├── 02_analysis.ipynb          # Levels accounting + OLS/2SLS regressions
├── 03_figures.ipynb           # Scatter plots (Figures 1 and 2)
└── README.md
```

Run notebooks **in order**: `01` → `02` → `03`.

The pipeline uses a single `VERSION` switch:
- `VERSION = 1` → Replication (PWT 5.6, 1988, Sachs-Warner) → outputs `merged_v1.csv`
- `VERSION = 3` → Extension (PWT 10.01, 2019, Fraser Area 5) → outputs `merged_v3.csv`

Run `01_data_pipeline.ipynb` twice (once with each VERSION) before running `02` and `03`.

---

## Data sources

Download each file and place it in the path matching `BASE` in the config cell.
The default `BASE` is `C:\Users\Adams\OneDrive\DE & E Research\Data`.

| File | Subfolder | Source | Notes |
|------|-----------|--------|-------|
| `pwt56_forweb.xlsx` | `GDP, Investment, and Labor force\` | [Penn World Tables 5.6](https://pwt.sas.upenn.edu/) | Replication only |
| `pwt1001.xlsx` | `GDP, Investment, and Labor force\` | [Penn World Tables 10.01](https://www.rug.nl/ggdc/productivity/pwt/) | Extension only |
| `BL2013_MF2599_v2.csv` | `Barro-Lee dataset\` | [Barro-Lee](http://www.barrolee.com/) | Replication: 1985 data |
| `OUP_proj_MF2564_v1.csv` | `Barro-Lee dataset\` | [Barro-Lee](http://www.barrolee.com/) | Extension: 2015 data |
| `3BResearchersDataset2017.xlsx` | `For Index\World Bank Governance Indicators\` | [ICRG — PRS Group](https://www.prsgroup.com/) | Proprietary. Contains 5 GADP sub-index sheets |
| `open.csv` | `For Index\SACHS-WARNER OPENNESS INDEX\` | Sachs-Warner | Annual binary openness 1950–1992 |
| `efotw-2025-master-index-data-for-researchers-iso.xlsx` | `For Index\SACHS-WARNER OPENNESS INDEX\` | [Fraser Institute EFW](https://www.fraserinstitute.org/economic-freedom/dataset) | Extension only, Area 5 |
| `Contribution of mining to value added.xlsx` | `MINING VALUE ADDED\` | [UN Statistics](https://unstats.un.org/unsd/snaama/) | Sheet index 1 used |
| `geo_instruments.csv` | `INSTRUMENTAL Variables Data\` | Replication file | Columns: iso3, distancefromeq |
| `frankel_romer_trade.csv` | `INSTRUMENTAL Variables Data\` | Replication file | Columns: iso3, fr_constructed_trade_share |
| `language_instruments.csv` | `INSTRUMENTAL Variables Data\` | Replication file | Columns: iso3, english_frac, we_lang_frac |

> **Note on ICRG:** The `3BResearchersDataset2017.xlsx` file is proprietary. It must contain
> five sheets named exactly: `A-Government Stability`, `F-Corruption`, `I-Law and Order`,
> `K-Democratic Accountability`, `L-Bureaucracy Quality`. Each sheet has a `Country` column
> and year columns. The pipeline auto-detects the header row.

---

## Software requirements

```
Python >= 3.8
numpy
pandas
scipy
matplotlib
openpyxl    # for reading ICRG .xlsx
xlrd        # for reading legacy .xls files (if needed)
```

Install with:
```bash
pip install numpy pandas scipy matplotlib openpyxl xlrd
```

No additional econometrics packages required. OLS and 2SLS are implemented from scratch
in `02_analysis.ipynb` with HC1 heteroskedasticity-robust standard errors.

---

## How to run

**Step 1 — Build replication dataset:**
Open `01_data_pipeline.ipynb`, set `VERSION = 1`, run all cells.
Produces: `merged_v1.csv` in `output_base`.

**Step 2 — Build extension dataset:**
In `01_data_pipeline.ipynb`, set `VERSION = 3`, run all cells.
Produces: `merged_v3.csv` in `output_base`.

**Step 3 — Run analysis:**
Open `02_analysis.ipynb`, run all cells.
Prints Tables 3–6 and robustness checks to the notebook output.

**Step 4 — Generate figures:**
Open `03_figures.ipynb`, run all cells.
Saves `figure1_replication_scatter.png` and `figure2_extension_scatter.png` to `output_base`.

---

## Key results

| | H&J (1999) | Replication (1988) | Extension (2019) |
|---|---|---|---|
| Sample N | 127 | 109 | 111 |
| OLS β(S) | 3.290 | 2.922 | 6.603 |
| 2SLS β(S) [4 instr.] | 5.140 | 5.450 | 8.415 |
| 2SLS β(S) [preferred, 3 instr.] | — | 5.447 | 8.480 |
| First-stage F [preferred] | — | 16.08 | 52.08 |
| Sargan p [preferred] | — | 0.175 | 0.013 |
| TFP share | 0.601 | 0.610 | 0.795 |
| USA/Niger Y/L gap | 35.0× | 34.9× | 43.2× |

**Preferred specification:** 3 instruments (DIST, ENGFRAC, WEFRAC).
FRTRADE dropped: Sargan test rejects 4-instrument spec at p = 0.037 (replication) and p = 0.031 (extension),
suggesting Frankel-Romer predicted trade has a direct income effect beyond the institutional channel.

---

## Documented deviations from Hall and Jones (1999)

| Deviation | Detail | Impact |
|-----------|--------|--------|
| Sample size | N = 109/111 vs 127 — ICRG coverage gaps for small island states and former Soviet republics | Minor |
| Sachs-Warner window | Ends 1992 in our data vs 1994 in original | Negligible |
| Mining value added | 1990 proxy used for 1988 replication | Negligible |
| Standard errors | HC1 robust throughout (vs bootstrap in original) | Asymptotically equivalent |

---

## Citation

Hall, Robert E., and Charles I. Jones, 1999, Why do some countries produce so much more output
per worker than others? *Quarterly Journal of Economics* 114, 83–116.
