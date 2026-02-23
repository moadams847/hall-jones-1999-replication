# Data Download Instructions

Place all files in this `data/` folder before running the pipeline.
Do not commit these files to GitHub — they are subject to the terms of their respective sources.

## File Checklist

| File                                         | Status            | Source                                                                   |
| -------------------------------------------- | ----------------- | ------------------------------------------------------------------------ |
| `pwt1001.xlsx`                               | Download required | https://www.rug.nl/ggdc/productivity/pwt/                                |
| `BL2013_MF2599_v2_2.csv`                     | Download required | http://www.barrolee.com/data/full1.htm                                   |
| `wgidataset_with_sourcedata-2025.xlsx`       | Download required | https://www.worldbank.org/en/publication/worldwide-governance-indicators |
| `open.csv`                                   | Download required | Sachs & Warner (1995) replication files                                  |
| `Contribution_of_mining_to_value_added.xlsx` | Download required | UNIDO / OECD                                                             |
| `geo_cepii.xls`                              | Download required | http://www.cepii.fr — GeoDist dataset                                    |
| `ling_web.csv`                               | Download required | http://www.cepii.fr — Language dataset                                   |

## Notes

- `frankel_romer_trade.csv` and `language_instruments.csv` and `geo_instruments.csv`
  are **generated automatically** by `hj_merge.py` — you do not need to download them separately.
- The Barro-Lee file may download with a numeric prefix in the filename
  (e.g. `1771859530022_BL2013_MF2599_v2_2.csv`). Rename it to `BL2013_MF2599_v2_2.csv`
  or update the `BL_FILE` path in `hj_merge.py` to match the downloaded filename.
- All file paths are set in the CONFIG section at the top of `hj_merge.py`.
  Update them to match your local folder structure if needed.
