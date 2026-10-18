# Reproducing the dashboard

Run these commands from the repository root after installing `requirements.txt`. Input locations are relative to the scripts, not to the original developer's computer. The checked-in website can be used without regeneration.

## Sid and the UF comparison

```sh
python site/prepare_sid.py
python site/prepare_events.py
python site/prepare_phylogeny.py
python site/prepare_ufpf.py
python site/prepare_ufpf_drivers.py
python scripts/check.py
```

These commands read the archived public inputs, overwrite their generated JSON/TSV files in `site/dist/`, and run integrity checks. `prepare_ufpf.py` also supports explicit input/output locations; inspect `python site/prepare_ufpf.py --help`.

Treatment records in `site/source-data/sid/events.tsv` and literature JSON were curated from cited sources. They are not automatically refreshed by these commands. Preserve source links, uncertainty notes and taxonomic applicability when editing them.

## Original iHMP/Greengenes view

The `count-audit/output/` files are included. To rebuild the browser dataset from those outputs:

```sh
python site/prepare_greengenes.py
```

To repeat the upstream R-data import and count reconciliation, install the optional dependencies, then run:

```sh
python -m pip install -r count-audit/requirements.txt
python count-audit/audit.py
python site/prepare_greengenes.py
```

The original metadata extraction was performed during research. Its source URLs and matching rules are retained in `data/metadata_provenance.json`, but a standalone extraction script was not saved. The checked-in metadata tables are the reproduction inputs; this repository does not claim a complete one-command reconstruction from the Internet.

## Compatibility audits

```sh
python ufpf-audit/audit.py
python oncobiotics-audit/audit.py
```

These use the included audit inputs. Cancer fallback results were produced during an interactive audit; the saved JSON, source lists and report are retained, but no standalone fallback-audit script existed at export time. The cancer files are not part of the current dashboard.

## Future functional profiling

See `functional-profiling/README.md`. That separate environment needs bioBakery tools, reference databases and substantial compute/storage. Nothing in the normal website/test commands launches it. The workflow has unit tests; it has not processed Sid's real reads end to end.

## Snapshot provenance

`docs/source-snapshot.json` records hashes of copied project files before repository documentation changes. It identifies the exported Sites commit and excluded paths. It is an export inventory, not a dependency lockfile or a claim that every copied external input has the same redistribution license.
