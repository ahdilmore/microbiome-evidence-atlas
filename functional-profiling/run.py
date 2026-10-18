"""Download, QC, remove host reads, and profile one visit. Standard-library only."""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def run(args, log):
    with log.open('a') as f:
        f.write('\n' + json.dumps([str(a) for a in args]) + '\n')
        f.flush()
        subprocess.run([str(a) for a in args], stdout=f, stderr=f, check=True)

def download(row, folder):
    path = folder / row['url'].rsplit('/', 1)[1]
    size = int(row['bytes'])
    if not path.exists():
        partial = path.with_suffix(path.suffix + '.partial')
        subprocess.run(['curl', '--fail', '--location', '--retry', '3',
                        '--continue-at', '-', '--output', str(partial), row['url']], check=True)
        if partial.stat().st_size != size:
            raise ValueError('Download size mismatch: ' + str(partial))
        partial.rename(path)
    if path.stat().st_size != size:
        raise ValueError('Existing file size mismatch: ' + str(path))
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(chunk)
    return path, digest.hexdigest()

def concatenate_fastq(paths, output):
    """Unique IDs distinguish lanes and mates; validates FASTQ and gzip CRC."""
    count = 0
    with gzip.open(output, 'wt', compresslevel=1) as dest:
        for source in paths:
            with gzip.open(source, 'rt') as f:
                while True:
                    header = f.readline()
                    if not header:
                        break
                    seq, plus, qual = f.readline(), f.readline(), f.readline()
                    if not header.startswith('@') or not plus.startswith('+') or len(seq.rstrip()) != len(qual.rstrip()) or not qual:
                        raise ValueError('Malformed FASTQ: ' + str(source))
                    count += 1
                    dest.write(f'@read_{count}\n{seq}+\n{qual}')
    if not count:
        raise ValueError('No QC-passing reads')
    return count

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', required=True)
    p.add_argument('--host-index', required=True, help='Bowtie2 human index prefix')
    p.add_argument('--protein-db', required=True)
    p.add_argument('--nucleotide-db')
    p.add_argument('--metaphlan-options', default='')
    p.add_argument('--translated-only', action='store_true', help='Exploratory protein-only profiling; no species attribution')
    p.add_argument('--threads', type=int, default=4)
    a = p.parse_args()
    if a.threads < 1 or (not a.translated_only and not a.nucleotide_db):
        p.error('Positive threads and nucleotide-db required (unless translated-only)')
    for tool in ['curl', 'fastp', 'bowtie2', 'diamond', 'humann', 'humann_renorm_table']:
        if not shutil.which(tool):
            p.error('Missing executable: ' + tool)
    versions = {}
    for tool, flag in [('fastp','--version'),('bowtie2','--version'),('diamond','version'),('humann','--version')]:
        result = subprocess.run([tool, flag], text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True)
        versions[tool] = result.stdout.strip()
    if not list(Path(a.host_index).parent.glob(Path(a.host_index).name + '*.bt2*')):
        p.error('Human Bowtie2 index not found')
    if not list(Path(a.protein_db).glob('*.dmnd')):
        p.error('Protein database directory must contain DIAMOND files')
    rows = [r for r in csv.DictReader((ROOT/'reads.tsv').open(), delimiter='\t') if r['sample_id'] == a.sample]
    if not rows:
        p.error('Unknown sample')
    out = ROOT/'work'/a.sample
    out.mkdir(parents=True, exist_ok=True)
    # Refuse accidental mixing of settings or incomplete runs. Raw downloads are reusable.
    if (out/'run.json').exists():
        p.error('Run already exists; archive its outputs before rerunning with changed settings')
    raw = out/'raw'; raw.mkdir(exist_ok=True)
    qc = out/'qc'; qc.mkdir(exist_ok=True)
    log = out/'commands.log'
    files = {}; hashes = {}
    for r in rows:
        print('Downloading/verifying', r['lane'], r['mate'], flush=True)
        files[r['lane'], r['mate']], hashes[r['url']] = download(r, raw)
    provenance = dict(settings=vars(a), versions=versions, reads=rows, sha256=hashes, status='running')
    (out/'run.json').write_text(json.dumps(provenance, indent=2))
    clean = []
    for lane in sorted({r['lane'] for r in rows}):
        r1, r2 = qc/f'{lane}_R1.fastq.gz', qc/f'{lane}_R2.fastq.gz'
        run(['fastp', '-i', files[lane,'R1'], '-I', files[lane,'R2'], '-o', r1, '-O', r2,
             '--detect_adapter_for_pe', '--qualified_quality_phred', '20', '--length_required', '60',
             '--thread', min(a.threads, 16), '--json', qc/f'{lane}.json', '--html', qc/f'{lane}.html'], log)
        clean.extend([r1, r2])
    merged = qc/'merged.fastq.gz'
    provenance['qc_passing_reads'] = concatenate_fastq(clean, merged)
    nonhost = qc/'nonhost.fastq.gz'
    # Each mate screened independently; only unmapped reads proceed. No host SAM retained.
    run(['bowtie2', '--very-sensitive-local', '-p', a.threads, '-x', a.host_index,
         '-U', merged, '--un-gz', nonhost, '-S', os.devnull], log)
    prof = out/'humann'
    cmd = ['humann', '--input', nonhost, '--output', prof, '--output-basename', a.sample,
           '--threads', a.threads, '--memory-use', 'minimum', '--protein-database', a.protein_db]
    if a.translated_only:
        cmd += ['--bypass-nucleotide-search']
    else:
        cmd += ['--nucleotide-database', a.nucleotide_db]
        if a.metaphlan_options:
            cmd += ['--metaphlan-options', a.metaphlan_options]
    run(cmd, log)
    for kind in ['pathabundance', 'genefamilies']:
        run(['humann_renorm_table', '--input', prof/f'{a.sample}_{kind}.tsv',
             '--output', prof/f'{a.sample}_{kind}_cpm.tsv', '--units', 'cpm', '--special', 'y'], log)
    provenance['status'] = 'complete'
    (out/'run.json').write_text(json.dumps(provenance, indent=2))
    print('Complete:', out)

if __name__ == '__main__':
    main()
