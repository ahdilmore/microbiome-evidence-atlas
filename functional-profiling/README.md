# Sid: read-supported functional profiling

Status: workflow prepared; **no Sid reads have been profiled yet**. Processing was deferred at the user’s request. This is an optional future workflow; no job is scheduled or submitted. Do not publish placeholders as pathway results.

## Verified inputs

`reads.tsv` contains 24 public FASTQ.gz objects (12 paired lanes), totaling 25,532,202,533 bytes. Every lane has both mates. Five visits correspond to the existing app: October 22 and December 2, 2025; February 13, March 10, and March 27, 2026. `samples.json` preserves reported/received dates, date conflicts, and alternate dates. Folder dates are not necessarily collection dates.

The archived complete S3 listing is `archive-inventory.xml`. The source is [Sid's microbiome archive](https://osteosarc.com/microbiome/). Three later visits are available but are outside this initial comparison. File sizes and lane pairing were checked against the listing; file contents/checksums are verified during download. S3 multipart ETags are not treated as MD5 checksums.

## Run on Linux HPC

1. Copy this directory to project/scratch storage, excluding `.env`, `work`, and databases. Allow several hundred GB of working storage pending the actual reference database expansion and intermediate sizes.
2. Create and activate the environment using `micromamba create -f environment.yml` or the equivalent conda command. This version specification still requires a target-host smoke test; it is not a platform-specific lockfile. Export the successfully resolved explicit package lock and record `humann --version`, `metaphlan --version`, `bowtie2 --version`, `diamond version`, and `fastp --version`.
3. Install the reference databases once in shared storage. HUMAnN commands:

   ```bash
   humann_databases --download chocophlan full /path/to/databases
   humann_databases --download uniref uniref90_diamond /path/to/databases
   ```

   Use the actual extracted directories as `SID_NUCLEOTIDE_DB` and `SID_PROTEIN_DB`. Use a single compatible release for all samples. Record database download URLs, release names and SHA256 hashes. HUMAnN 3.9's database catalog specifies ChocoPhlAn v201901_v31 and UniRef90 v201901b. Avoid silently upgrading databases midway through the study.

4. Download and extract the [official human Bowtie2 reference](https://huttenhower.sph.harvard.edu/kneadData_databases/Homo_sapiens_hg39_T2T_Bowtie2_v0.1.tar.gz). Set `SID_HOST_INDEX` to the actual prefix shared by its `.bt2` files, not the directory. The official download was reachable during preparation (3,824,531,235 compressed bytes). Record the reference checksum. Raw host reads remain only in the input files; no host alignment file is retained or uploaded.
5. Install the MetaPhlAn 4.1.1 marker database with the fixed `mpa_vJun23_CHOCOPhlAnSGB_202307` index matching the previous taxonomy work. Set `SID_METAPHLAN_OPTIONS` to `--bowtie2db /actual/database/path --index mpa_vJun23_CHOCOPhlAnSGB_202307`. Verify MetaPhlAn/HUMAnN compatibility with a small test before processing all data; do not reuse the public integer-weight species CSV as a native MetaPhlAn profile.
6. Run `python -m unittest test_workflow.py`, then the official HUMAnN installation tests. Profile October first as the end-to-end pilot. Review trimming, host fraction, mapping rates, pathway output, runtime, and memory before running the remaining samples.
7. Submit `sbatch slurm.sh` with the site's required account/partition arguments and exported `SID_*` variables. The template initially requests 8 CPUs, 64 GB RAM and 24 hours and runs visits sequentially. Adjust based on the pilot and cluster policy. It assumes compute nodes have external download access; otherwise stage downloads on an approved transfer node first.

For a single visit outside Slurm:

```bash
python run.py --sample H37004-stool-001 --threads 8 \
  --host-index "$SID_HOST_INDEX" --protein-db "$SID_PROTEIN_DB" \
  --nucleotide-db "$SID_NUCLEOTIDE_DB" --metaphlan-options="$SID_METAPHLAN_OPTIONS"
```

The runner downloads all lanes, checks expected bytes, computes SHA256 hashes, trims paired reads with fastp, gives mates unique identifiers, screens each read against the human reference, and runs HUMAnN. QC uses Q20, at most fastp's default fraction of low-quality bases, and minimum length 60; it is not a claim that every retained base has Q20. Reads are concatenated, not overlap-merged. HUMAnN treats these as individual reads. The runner keeps original outputs and logs, and normalizes gene-family and pathway abundance to copies per million including special features in the denominator.

Interrupted profiling leaves `run.json` marked running and refuses automatic reuse. Archive the partial sample directory or inspect it before retrying; don't label a partial run complete. Downloads before the profiling stage support curl resume. This conservative initial runner does not yet implement a scheduler-aware checkpoint/resume system.

## Comparison outputs

`summarize.py` requires all five complete runs with matching profiling settings. It writes `pathway-summary.json`, preserving community totals separately from stratified taxon rows and unmapped/unintegrated features. It provides consecutive-visit deltas, first-to-last deltas, median abundance, range/median, and consistently detected pathways. Default minimum is an exploratory 10 copies per million; vary it and report sensitivity. This threshold is not a universal detection limit.

Before interpreting stability, compare QC and mapped-read depth across visits and repeat the analysis on equal-depth, randomly subsampled non-host reads if needed. Persistence and low variability are descriptive, not statistical evidence of equivalence. Nondetection is not proven biological absence. Pathway coverage is distinct from abundance; review available reaction support/coverage before interpreting weakly supported pathways. Retain raw gene-family outputs for more targeted functional questions.

Read-supported DNA profiles measure functional potential, not expression, metabolite production, or treatment effects. Use treatment annotations as temporal context only. After QC, add actual pathway trajectories and persistent-pathway summaries to the app, with provenance, mapping rates, and literature tied to the specific pathway and evidence type.

## Fallback hierarchy

If the full nucleotide-plus-protein workflow is too expensive, a consistent protein-only HUMAnN run (`--translated-only`) with its EC-filtered UniRef90 database is a read-supported exploratory alternative. It loses the standard species-stratified nucleotide component and restricts functional breadth; do not mix those results with full-database profiles.

If reads cannot be processed at all, use reference-genome gene content weighted by taxon abundance, preferably multiple genomes per species/SGB. Quantify unmatched abundance and between-strain uncertainty; infer pathway completeness from gene sets rather than summing arbitrary pathway labels. Label results **taxon-inferred functional potential**. This cannot independently establish functional stability because it is derived from the same taxonomic table. No such estimates have been computed in this package.

References: [HUMAnN documentation](https://github.com/biobakery/humann), [bioBakery 3 methods](https://doi.org/10.7554/eLife.65088), [KneadData host-filtering documentation](https://github.com/biobakery/kneaddata).

## Local validation and limitations

Three focused tests pass: FASTQ validation/unique IDs, separation of stratified/special rows, and rejection of invalid numeric abundances. The real-data pipeline is not yet validated end to end. The local Mac packaging/architecture conflicts were repaired; HUMAnN 3.9 starts, native Bowtie2 and DIAMOND start, and the official HUMAnN test suite passed 186 tests. Target-server verification is still required before launch. The site has not been modified and contains no new functional results.
