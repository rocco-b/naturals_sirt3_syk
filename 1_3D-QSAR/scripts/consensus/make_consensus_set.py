"""Builds the 3D-QSAR consensus set: the compounds predicted active by BOTH models.

The compounds are taken from the list of the predicted actives of the SIRT3 model
(SIRT3/4_screening/svm_top_candidates.csv, sorted by decreasing predicted probability) and kept
when the SYK model (SYK/4_screening/screening_results.csv) also predicted them active. The
position of a compound in the resulting list is its docking ID (the numbers 2314, 3843, 4109
of the manuscript).

Outputs (written in --output-dir):
    consensus_set_with_probabilities.csv   identifier, Smiles, docking_ID, SYK_probability, SIRT3_probability
    consensus_qsar.txt                     tab-separated: identifier, Smiles, ID
    ligands.smi                            SMILES strings, one per line, same order

Usage (from any folder):
    python make_consensus_set.py --output-dir <folder>
"""
import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]          # 1_3D-QSAR/

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--syk-results', default=str(ROOT / 'SYK' / '4_screening' / 'screening_results.csv'))
ap.add_argument('--sirt3-top', default=str(ROOT / 'SIRT3' / '4_screening' / 'svm_top_candidates.csv'))
ap.add_argument('--output-dir', default='.')
args = ap.parse_args()

syk = pd.read_csv(args.syk_results)
sirt3 = pd.read_csv(args.sirt3_top)
syk_active = syk[syk['Predicted_Active'] == 1]
print(f'Predicted active on SYK: {len(syk_active):,}')
print(f'Predicted active on SIRT3: {len(sirt3):,}')

cons = (sirt3[['identifier', 'Smiles', 'Active_Probability']]
        .rename(columns={'Active_Probability': 'SIRT3_probability'})
        .merge(syk_active[['identifier', 'Active_Probability']].rename(columns={'Active_Probability': 'SYK_probability'}),
               on='identifier', how='inner'))
cons.insert(2, 'docking_ID', range(1, len(cons) + 1))
cons = cons[['identifier', 'Smiles', 'docking_ID', 'SYK_probability', 'SIRT3_probability']]
print(f'Consensus set (active on both targets): {len(cons):,}')

out = Path(args.output_dir)
out.mkdir(parents=True, exist_ok=True)
cons.to_csv(out / 'consensus_set_with_probabilities.csv', index=False)
cons[['identifier', 'Smiles', 'docking_ID']].rename(columns={'docking_ID': 'ID'}).to_csv(
    out / 'consensus_qsar.txt', sep='\t', index=False)
cons['Smiles'].to_csv(out / 'ligands.smi', index=False, header=False)
print(f'Saved in: {out.resolve()}')
