"""Applies the Lipinski filters to the annotated SuperNatural 3.0 table and writes the
two-column input (identifier, Smiles) of the 3D-QSAR screening scripts.

Filters (all inclusive), applied to the four columns of the annotated table:
    Molweight <= 500, cLogP <= 5, H-Acceptors <= 10, H-Donors <= 5

The input is the tab-separated annotated table (plain text or a .zip that contains one text file).
Applied to the table provided in this repository (65,564 compounds, already filtered) the script
returns all the compounds and reproduces supernatural3_cleaned_db.csv; applied to the full
annotated table of SuperNatural 3.0 (124,943 compounds, with the same four columns) it returns
the 65,564 compounds used in the study.

Usage (from any folder):
    python make_cleaned_db.py --input <annotated table> --output supernatural3_cleaned_db.csv
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]          # 1_3D-QSAR/
DEFAULT_INPUT = ROOT / 'supernatural3_database' / 'supernatural3_Lipinski.txt.zip'

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--input', default=str(DEFAULT_INPUT), help='annotated table (.txt or .zip), tab-separated')
ap.add_argument('--output', default='supernatural3_cleaned_db.csv', help='output file (tab-separated)')
args = ap.parse_args()

REQUIRED = ['identifier', 'canonical_smiles', 'Molweight', 'cLogP', 'H-Acceptors', 'H-Donors']
df = pd.read_csv(args.input, sep='\t', low_memory=False)
missing = [c for c in REQUIRED if c not in df.columns]
if missing:
    sys.exit(f'Missing columns in the input table: {missing}')

keep = ((df['Molweight'] <= 500) & (df['cLogP'] <= 5) &
        (df['H-Acceptors'] <= 10) & (df['H-Donors'] <= 5))
print(f'Compounds in the input table: {len(df):,}')
print(f'Compounds passing the Lipinski filters: {int(keep.sum()):,}')

out = df.loc[keep, ['identifier', 'canonical_smiles']].rename(columns={'canonical_smiles': 'Smiles'})
out.to_csv(args.output, sep='\t', index=False)
print(f'Saved: {args.output}')
