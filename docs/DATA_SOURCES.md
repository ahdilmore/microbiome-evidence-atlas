# Data sources and attribution

| Component | Primary source | Local provenance |
| --- | --- | --- |
| Sid's public stool profiles and timing | https://osteosarc.com/microbiome/ | `site/dist/sid-provenance.json`, `site/source-data/sid/` |
| Treatment overlays | URLs in the event annotations | `site/source-data/sid/events.tsv`, `site/dist/sid-events.json` |
| MetaPhlAn reference tree | https://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.nwk | `site/dist/sid-phylogeny.json` |
| UF/UFPF cohort | https://github.com/maevekrueger/UFPF_metagenomics and https://zenodo.org/records/10912505 | `ufpf-audit/sources.json`, `site/dist/sid-ufpf.json` |
| Original Stanford abundance table | http://ipop_public.gbsc.os.scg.stanford.edu/HMP/gut_16s_abundance/gut_16s_abundance.txt | `data/metadata_provenance.json` |
| HMP2Data count source | https://github.com/dozmorovlab/HMP2Data | `count-audit/audit-log.json`, `site/dist/data-provenance.json` |
| ONCOBIOTICS deposit | https://doi.org/10.6084/m9.figshare.30336373 | `oncobiotics-audit/sources.json`, `oncobiotics-audit/audit-report.txt` |
| CRC fallback | https://zenodo.org/records/15069069 | `cancer-fallback-audit/crc-audit-results.json` |
| Melanoma fallback and external R scripts | https://github.com/omicsEye/Cancer-Microbiome | `cancer-fallback-audit/melanoma-sources.json`; pinned commit `61bd2c75927a77051dd409e284e02b7f50f8286b` |

Evidence-card citations and limitations are stored in `site/dist/literature.json` and `site/dist/sid-evidence.json`. Use those individual studies for claims; dataset membership alone is not evidence of an association.

The archived public inputs and third-party scripts remain attributable to their original authors and retain their applicable licenses/terms. The melanoma repository's MIT license is retained alongside its copied files. No new license is applied to research articles, clinical reports, datasets, or reference databases by this export.
