#!/usr/bin/env python3
"""Download the CARE cohort's bulk-DNA copy-number segments (the ground truth) and the authors' cell
annotations from Synapse (project syn63600501).

    python3 reproduce/analysis/multiome/download_synapse.py            # tables + annotation files
    python3 reproduce/analysis/multiome/download_synapse.py --tables   # tables only
    python3 reproduce/analysis/multiome/download_synapse.py --list     # print the project tree only
Needs synapseclient and a ~/.synapseConfig with a personal access token.
Output: ATAC_RNA/Data/synapse/  (MINTCNA_DATA must be set, see reproduce/README.md)
"""
import argparse
import csv
import os
import shutil
import sys

OUT = os.path.join(os.environ["MINTCNA_DATA"], "ATAC_RNA/Data/synapse")
PROJECT = "syn63600501"
TABLES = {                                          # Synapse tables, read with tableQuery()
    "variants_gatk_seg":         "syn63989096",     # GATK copy-number segments: the ground truth
    "sample_identifier_mapping": "syn63972901",     # join keys from aliquot barcode ...
    "biospecimen_aliquots":      "syn75099074",     # ... to the GEO sample
}
WANTED_FILES = (                                    # annotation files, located by name
    "snatac_celltype_meta_data.RDS",
    "snrna_celltype_meta_data.RDS",
    "snrna_malignant_meta_data.RDS",
    "bulk_dna_methylation_classification.tsv",
)


def connect():
    try:
        import synapseclient
    except ImportError:
        sys.exit("synapseclient is not installed (pip install synapseclient)")
    syn = synapseclient.Synapse()
    try:
        syn.login()
    except Exception as e:
        sys.exit(f"Synapse login failed: {e}\n~/.synapseConfig needs\n  [authentication]\n  authtoken = <token>")
    print(f"logged in as {syn.getUserProfile().get('userName', '?')}")
    return syn


def get_tables(syn):
    os.makedirs(OUT, exist_ok=True)
    for name, synid in TABLES.items():
        dest = os.path.join(OUT, f"{name}.csv")
        res = syn.tableQuery(f"SELECT * FROM {synid}")
        shutil.copyfile(res.filepath, dest)
        with open(dest, newline="") as fh:
            n = sum(1 for _ in csv.reader(fh)) - 1
        print(f"{name} ({synid}): {n:,} rows -> {dest}")


def walk_project(syn, download):
    import synapseutils
    os.makedirs(OUT, exist_ok=True)
    for dirpath, _dirnames, filenames in synapseutils.walk(syn, PROJECT):
        for fname, synid in filenames:
            wanted = fname in WANTED_FILES
            print(f"  {dirpath[0]}/{fname} ({synid}){'  <- wanted' if wanted else ''}")
            if wanted and download:
                syn.get(synid, downloadLocation=OUT, ifcollision="overwrite.local")


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
        walk_project(syn, download=True)
    print(f"done -> {OUT}")


if __name__ == "__main__":
    main()
