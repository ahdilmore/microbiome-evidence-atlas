"""Build clinical context from the downloaded public events-table snapshot."""
from pathlib import Path
import csv, hashlib, json
root = Path(__file__).resolve().parent
source = root / 'source-data/sid/events.tsv'
events = list(csv.DictReader(source.open(), delimiter='\t'))
assert all(e['date'] and (not e['end_date'] or e['end_date'] >= e['date']) for e in events)
(root / 'dist/sid-events.json').write_text(json.dumps({
 'source': 'https://osteosarc.com/data/events.tsv',
 'overview': 'https://osteosarc.com/timeline/',
 'retrieved': '2026-10-17',
 'sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
 'events': events,
}, indent=2))
