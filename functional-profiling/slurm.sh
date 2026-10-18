#!/usr/bin/env bash
#SBATCH --job-name=sid-functions
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=sid-functions-%j.log
set -euo pipefail
# Run sbatch from the workflow directory, with the environment already activated.
# These are initial resource requests, not a measured runtime guarantee.
: "${SID_HOST_INDEX:?Set the human Bowtie2 index prefix}"
: "${SID_PROTEIN_DB:?Set the HUMAnN UniRef90 database directory}"
: "${SID_NUCLEOTIDE_DB:?Set the HUMAnN ChocoPhlAn database directory}"
: "${SID_METAPHLAN_OPTIONS:?Set fixed MetaPhlAn --bowtie2db and --index options}"
for sample in H37004-stool-001 H37004-stool-002 H37004-stool-003 H37004-stool-004 H37004-stool-005; do
  if python -c 'import json,sys; from pathlib import Path; p=Path("work")/sys.argv[1]/"run.json"; sys.exit(0 if p.exists() and json.loads(p.read_text()).get("status")=="complete" else 1)' "$sample"; then
    echo "Already complete: $sample (settings checked by final summarizer)"
    continue
  fi
  python run.py --sample "$sample" --threads "${SLURM_CPUS_PER_TASK:-8}" \
    --host-index "$SID_HOST_INDEX" --protein-db "$SID_PROTEIN_DB" \
    --nucleotide-db "$SID_NUCLEOTIDE_DB" --metaphlan-options="$SID_METAPHLAN_OPTIONS"
done
python summarize.py --minimum-cpm 10
