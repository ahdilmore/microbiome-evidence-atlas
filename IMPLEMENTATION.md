# Microbiome evidence explorer

## Verified starting data

Source: http://ipop_public.gbsc.os.scg.stanford.edu/HMP/gut_16s_abundance/gut_16s_abundance.txt
Local copy: data/gut_16s_abundance.tsv

The table contains 855 unique SampleID rows and 96 taxonomic features: 6 phyla, 12 classes, 12 orders, 21 families, and 45 genera. No clinical or demographic metadata is included. Values are fractional, not raw sequencing counts. Within-rank totals are sometimes below one (genus totals range approximately 0.000424–1). Do not silently renormalize or treat multiple ranks as independent parts of one composition. Sample identifiers alone do not establish participant or visit relationships.

## First usable workflow

1. Import this abundance table and a CSV/TSV sample metadata table. Explicitly map the sample ID, report unmatched and duplicate IDs, and require one metadata row per sample after resolving duplicates. Preserve original files and provenance.
2. Choose taxonomic rank (default genus), exposure column, categorical reference/comparison or continuous exposure, covariates, and participant ID for repeated measurements. Report exclusions, missing values, group sample sizes and independent participant counts before analysis.
3. Choose a method appropriate to the confirmed abundance representation and study design. Relative abundance inputs require a suitable relative abundance association model; do not feed proportions to count-based methods or reconstruct counts. Raw counts would permit evaluating ANCOM-BC2. Include prevalence filtering, multiplicity correction across all tested taxa for the chosen contrast, uncertainty, and a documented zero-handling sensitivity analysis when transformations need one. Reject unsupported or non-identifiable designs.
4. Display a ranked effect plot with uncertainty, q-values and prevalence. Clicking a taxon shows observed distributions and participant trajectories where applicable. Clearly label the effect scale and reference group. Export all tested taxa, not only significant hits, with method settings and exclusions.
5. Retrieve literature for selected taxa and the exposure context. Use separate retrieval, evidence extraction and synthesis stages. Resolve taxon synonyms and flag ambiguous legacy names. Cache by taxon, context, query and retrieval date.
6. Show evidence cards containing a concise claim, supporting passage, PMID/DOI and link, taxonomic scope, study design, population, sample size if reported, direction, limitations and contradictions. Distinguish abstract-only from full-text evidence. Keep general biological background separate from evidence matching the selected exposure.
7. Produce a plain-language summary grounded only in verified evidence cards. Explain population association, uncertainty and applicability. Never imply that a genus-level result establishes a strain-specific function, causality, an individual diagnosis or a treatment recommendation. Show an explicit no-evidence state when retrieval yields no relevant support.

## Screen layout

Left: dataset quality and analysis controls.
Center: ranked effect plot, sortable results, selected-taxon distribution.
Right: literature evidence with expandable source passages and plain-language explanation.

Results and literature have separate statuses. Literature coverage must not change statistical significance. No synthetic findings should be presented as results from this dataset.

## Delivery dependencies

- Matching sample metadata and explanation of participant/visit/replicate identifiers.
- Abundance preprocessing provenance, particularly incomplete within-rank totals.
- A statistical compute service for fitted models; a web frontend alone does not supply an R analysis runtime.
- A literature retrieval service and configured model access for live literature-agent synthesis. No model credentials belong in browser code.

## Method references

Lin and Peddada, ANCOM-BC2: https://www.nature.com/articles/s41592-023-02092-7
Nearing et al., comparison across 38 datasets: https://www.nature.com/articles/s41467-022-28034-z

Method selection is conditional on data and design; benchmark findings do not establish one universally best differential abundance method.
