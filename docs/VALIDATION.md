# Export validation — 2026-10-17

The exported application source matches the existing site checkout at commit `602ebe9a80996d7303e2016a422f708c74d37186` (hosted version 12). Every original project code file was copied without changes; new repository documentation and a check runner were added. The functional-profiling README was corrected to say that real-data processing was deferred.

`python scripts/check.py` passed:

- Five JavaScript test files: cohort statistics, Sid normalization and evidence matching, abundance/prevalence filters, longitudinal stability and treatment overlap, taxonomic aggregation and scale changes.
- Three Python test files: hand-calculated phylogenetic distances, independently reconstructed UF distances and ordination, and cohort-specific taxon contributions and visit consistency.
- Three functional-workflow unit tests, using synthetic temporary inputs.

Validation used Python 3.12.14, Node.js 24.19.0, NumPy 2.3.5, pandas 2.2.3 and Biopython 1.85. The analysis dependencies are declared in the root requirements file; the check command does not need `rdata`.

The full raw-read/HUMAnN workflow was not run. The R-data import and all source download procedures were not rerun for this export. Existing result JSON and public inputs are retained. No live website deployment was changed.

A credential-pattern scan of the exported text found no credential material; its only match was the ordinary visit-token variable in the Greengenes preparation script. Git metadata, hosting configuration, environment directories and temporary deployment bundles were excluded. No file exceeds GitHub's 100 MiB individual-file limit; the largest archived input is approximately 18 MB.
