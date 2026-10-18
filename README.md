# Microbiome Evidence Atlas

An exploratory microbiome dashboard centered on Sid's five public longitudinal stool profiles, treatment timing, literature context, and comparison with the University of Florida Parkinson's Foundation (UFPF) cohort. The original iHMP/Greengenes explorer is also included.

Live site: https://microbiome-evidence-atlas.amandadilmore.chatgpt.site/

This repository includes the version 12 data export, local research workflows, and a guided Research frontend. It is a standalone snapshot; pushing here does **not** automatically deploy the hosted site.

## Run the website

The website is plain HTML, CSS and JavaScript modules. No npm install, frontend build, API key, or backend is needed to view the included data.

The homepage and `research-story.html` provide a guided interpretation of Sid's five observations. They read the prepared `sid-data.json`, `sid-ufpf-drivers.json`, `sid-ufpf-profiles.tsv`, and `sid-events.json` files alongside reviewed literature and treatment-context records. No raw-data pipeline runs in the browser. The established `sid.html` explorer remains available as the broader view.

```sh
python3 -m http.server 8000 --bind 127.0.0.1 --directory site/dist
```

Open http://localhost:8000/ for Research, http://localhost:8000/sid.html for the full dashboard, or http://localhost:8000/ihmp.html for the original cohort explorer. Use an HTTP server rather than opening HTML directly, because the application fetches data files.

## Repository contents

| Directory | Contents |
| --- | --- |
| `site/dist/` | Editable website source and the currently displayed data; this is a buildless app, despite the directory name. |
| `site/prepare_*.py` | Data preparation, treatment annotations, phylogenetic distances, UF ordination and taxon contributions. |
| `site/test-*` | Existing statistical, longitudinal, phylogenetic and cohort-comparison checks. |
| `site/source-data/sid/` | Public input profiles, event annotations, dates, report and reference tree. |
| `data/`, `count-audit/` | Original Stanford metadata, HMP2Data count inputs and count reconciliation. |
| `ufpf-audit/` | Public UF inputs, provenance and compatibility audit. |
| `oncobiotics-audit/` | Downloaded cancer inputs and audit explaining why a joint comparison was deferred. |
| `cancer-fallback-audit/` | Downloaded CRC/melanoma candidates and compatibility findings; not integrated into the website. |
| `functional-profiling/` | Optional future HUMAnN/HPC workflow. No real Sid functional profiles have been generated. |
| `osteosarc-audit/` | Original downloaded public tables used to inspect Sid's timeline. |
| `docs/` | Reproduction instructions, source inventory and validation results. |

Literature cards are curated JSON records with citations and evidence limitations. The deployed website does not run an autonomous literature-search service. The `literature-*.json` files retain earlier research outputs.

## Validate and regenerate

Python 3.11+ and Node.js 20+ are recommended for development. The website itself only needs an HTTP server. Create a Python environment for analysis:

```sh
python3 -m venv .venv
# Windows: .venv\Scripts\activate
source .venv/bin/activate
python -m pip install -r requirements.txt
python scripts/check.py
```

The check command runs the existing JavaScript and Python tests, including independent weighted UniFrac calculations and UF comparison checks. It does not download sequencing reads or submit jobs.

See [REPRODUCING.md](docs/REPRODUCING.md) for the ordered data preparation commands and [VALIDATION.md](docs/VALIDATION.md) for the checks actually run on this export.

## Interpretation and provenance

- Sid's displayed proportions are normalized reported SGB abundance weights. Their underlying integer units are not established as sequencing counts, cell counts, or absolute microbial load.
- Five visits belong to **one person**. Treatment overlays show temporal overlap and do not establish treatment effects. Some displayed dates are receipt dates or have conflicting source dates.
- UF and Sid use the same reported MetaPhlAn reference database, but complete preprocessing equivalence has not been established. Cross-study distances are descriptive and are not a diagnostic classifier.
- The CRC fallback uses an older reference database. The melanoma deposit has a paper/header discrepancy. Neither has been silently added to Sid's joint ordination.
- The iHMP explorer aggregates repeated visits per participant; it is separate from Sid's shotgun analysis.

Source URLs, checksums, method details and exclusions accompany the data. See [DATA_SOURCES.md](docs/DATA_SOURCES.md). Publicly downloadable data retain their original terms; this export does not assign a new license to external data or research code. No blanket open-source license has been selected for the original project code.

## Hosting and updates

Any ordinary static host can serve `site/dist/`. The original Sites deployment identity, credentials, nested Git history, local software environments, caches and temporary deployment bundles are excluded. To update the existing hosted site, use its Sites workflow; to deploy elsewhere, configure that host separately.

Start future changes in a branch and run `python scripts/check.py` before publishing. The included functional-profiling workflow is intentionally manual and must not be treated as completed analysis.

The research story includes one explorer for five bacterial groups, backed by `immune-evidence.json`. Bacteroides focuses on B. uniformis and B. thetaiotaomicron, with an expandable five-species breakdown and reviewed roles/cautions in `bacteroides-species.json`. The catalog records study design, cancer/treatment context, exact taxonomic scope, reported direction, limitations and primary links. Profile measurements are joined by genus or exact species name, with all matching SGBs summed before calculating Control quantiles. Missing matches/measurements remain unavailable; explicit zero measurements are reported as non-detections. These are research leads, not a complete literature review or an immune-strength/clinical-response score.

One shared treatment section has its own interval selector, separate from bacterial sample selection. Reviewed explanations in `treatment-evidence.json` and procedure notes in `antibiotic-context.json` describe documented exposures, not inferred causes or unrecorded perioperative antibiotics. The largest-change link uses absolute consecutive percentage-point change, not statistical significance. March receipt dates and February date conflicts retain temporal-order uncertainty. Supporting evidence forms one vertical stack of expandable items.
