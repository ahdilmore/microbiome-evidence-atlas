"""Run existing checks without downloading reads, submitting jobs, or deploying."""
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

def main():
    node = shutil.which('node')
    if not node:
        raise SystemExit('Node.js is required for JavaScript tests; install Node.js 20+ and retry.')
    for path in sorted((ROOT / 'site').glob('test-*.mjs')):
        print(f'Running {path.relative_to(ROOT)}', flush=True)
        subprocess.run([node, str(path)], cwd=ROOT / 'site', check=True)
    for name in ('test-phylogeny.py', 'test-ufpf.py', 'test-ufpf-drivers.py'):
        print(f'Running site/{name}', flush=True)
        subprocess.run([sys.executable, name], cwd=ROOT / 'site', check=True)
    subprocess.run([sys.executable, '-m', 'unittest', 'test_workflow.py'],
                   cwd=ROOT / 'functional-profiling', check=True)
    print('All existing checks passed. No raw-read analysis or deployment was run.')

if __name__ == '__main__':
    main()
