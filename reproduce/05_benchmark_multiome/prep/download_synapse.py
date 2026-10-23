#!/usr/bin/env python3
"""
Download the bulk-DNA copy-number ground truth and author cell annotations for the
CARE IDH-mutant glioma cohort from Synapse (project syn63600501).

Prerequisites (one-time, see Notes/synapse_setup.md):
    module load conda
    pip install --user synapseclient
    ~/.synapseConfig holding your personal access token, chmod 600

Usage:
    python3 download_synapse.py            # tables + annotation files
    python3 download_synapse.py --tables   # tables only (fast, small)
    python3 download_synapse.py --list     # just print the project tree, download nothing

Note: syn63989096 is a Synapse TABLE, not a file, so it is retrieved with
tableQuery() rather than get(). `synapse get syn63989096` would fail.
"""
import argparse
import os
import sys

OUT = "<MINTCNA_ROOT>/ATAC_RNA/Data/synapse"
PROJECT = "syn63600501"          # CARE IDH-mutant

# Tables worth pulling. syn63989096 is the ground truth; the other two are the
# join keys that map its aliquot_barcode back to a GEO sample.
TABLES = {
    "variants_gatk_seg":         "syn63989096",   # <- the CNA ground truth
    "sample_identifier_mapping": "syn63972901",
    "biospecimen_aliquots":      "syn75099074",
}

# Author cell annotations. Located by name because their synIDs are not published.
# The celltype files matter more than the malignant one: cell type comes from marker
# genes and brain-atlas mapping, whereas the malignant flag was partly defined by
# inferred CNAs, so using it to pick diploid reference cells is circular.
WANTED_FILES = (
    "snatac_celltype_meta_data.RDS",     # syn74567782
    "snrna_celltype_meta_data.RDS",      # syn74567720
    "snrna_malignant_meta_data.RDS",     # syn74567748 -- context only, not truth
    "bulk_dna_methylation_classification.tsv",   # syn64001721, independent subtype check
)


def connect():
    try:
        import synapseclient
    except ImportError:
        sys.exit("synapseclient not installed.  module load conda && pip install --user synapseclient")
    syn = synapseclient.Synapse()
    try:
        syn.login()                      # reads ~/.synapseConfig
    except Exception as e:
        sys.exit(f"Synapse login failed: {e}\n"
                 "Check that ~/.synapseConfig contains:\n"
                 "  [authentication]\n  authtoken = <your token>")
    who = syn.getUserProfile()
    print(f"logged in as {who.get('userName', '?')}")
    return syn


def get_tables(syn):
    """tableQuery already lands a CSV on disk; copy that rather than going through
    asDataFrame(), which is the part most likely to break across pandas versions
    (this env has pandas 3.x)."""
    import csv as _csv
    import shutil
    os.makedirs(OUT, exist_ok=True)
    for name, synid in TABLES.items():
        dest = os.path.join(OUT, f"{name}.csv")
        print(f"  querying {name} ({synid}) ...", end="", flush=True)
        try:
            res = syn.tableQuery(f"SELECT * FROM {synid}")
            shutil.copyfile(res.filepath, dest)
        except Exception as e:
            print(f" FAILED: {e}")
            continue
        with open(dest, newline="") as fh:
            rdr = _csv.reader(fh)
            header = next(rdr, [])
            n = sum(1 for _ in rdr)
        print(f" {n:,} rows -> {dest}")
        if name == "variants_gatk_seg":
            print(f"      columns: {header}")


def walk_project(syn, download=False):
    """List the project tree; optionally download the annotation files."""
    import synapseutils
    os.makedirs(OUT, exist_ok=True)
    found = 0
    for dirpath, _dirnames, filenames in synapseutils.walk(syn, PROJECT):
        folder = dirpath[0]
        for fname, synid in filenames:
            hit = fname in WANTED_FILES
            mark = "  <== wanted" if hit else ""
            print(f"    {folder}/{fname}  ({synid}){mark}")
            if hit and download:
                syn.get(synid, downloadLocation=OUT, ifcollision="overwrite.local")
                found += 1
    if download:
        print(f"  downloaded {found} annotation file(s)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tables", action="store_true", help="tables only")
    ap.add_argument("--list", action="store_true", help="list the project tree, download nothing")
    args = ap.parse_args()

    syn = connect()
    if args.list:
        walk_project(syn, download=False)
        return
    get_tables(syn)
    if not args.tables:
        print("scanning project for annotation files ...")
        walk_project(syn, download=True)
    print(f"\ndone -> {OUT}")


if __name__ == "__main__":
    main()
