"""Extracts the GNINA scores of the consensus set from the two docked.sdf files and selects the
compounds above the CNN_VS cutoff on both targets.

Outputs (written in --output-dir):
    docking_scores_both_targets.csv   one row per compound: CNNscore, CNNaffinity and CNN_VS for SIRT3 and SYK,
                                      and a flag for the compounds with CNN_VS >= cutoff on both targets
    selected_hits.csv                 the selected compounds with the 3D-QSAR probabilities (Table 2 of the manuscript)

Convention of the tables: CNNscore is rounded to 4 decimals and CNNaffinity to 3 decimals (half up), and
CNN_VS is the product of the two rounded values, rounded to 4 decimals. This differs by at most 0.0008
from the CNN_VS stored in docked.sdf (unrounded product), and the selection is the same with both
(checked for the cutoffs 5.9, 6.0, 6.1, and 6.3). Compared with the table provided in 2_Docking/results/,
the output is identical except for 25 of the 8,354 CNN_VS values, which differ by 0.0001 in the last
digit (products that end exactly in 5 at the fifth decimal, rounded differently).

Usage (from any folder):
    python extract_docking_scores.py --output-dir <folder> [--cutoff 6.0]
"""
import argparse
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]          # repository root

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--sirt3-sdf', default=str(ROOT / '2_Docking' / 'SIRT3' / 'docked.sdf'))
ap.add_argument('--syk-sdf', default=str(ROOT / '2_Docking' / 'SYK' / 'docked.sdf'))
ap.add_argument('--consensus', default=str(ROOT / '1_3D-QSAR' / 'consensus_set_with_probabilities.csv'))
ap.add_argument('--cutoff', type=float, default=6.0)
ap.add_argument('--output-dir', default='.')
args = ap.parse_args()


def read_sdf_fields(path):
    """Returns the data fields (> <name> blocks) of every record of an SDF file as a DataFrame."""
    records, cur, key = [], {}, None
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line.startswith('>') and '<' in line:
                key = line[line.index('<') + 1:line.rindex('>')]
            elif line == '$$$$':
                records.append(cur)
                cur, key = {}, None
            elif key is not None:
                if line == '':
                    key = None
                else:
                    cur[key] = line
    return pd.DataFrame(records)


def q(value, places):
    return Decimal(value).quantize(Decimal(1).scaleb(-places), rounding=ROUND_HALF_UP)


def scores(path):
    d = read_sdf_fields(path)
    score = [q(v, 4) for v in d['CNNscore']]
    aff = [q(v, 3) for v in d['CNNaffinity']]
    vs = [q(str(s * a), 4) for s, a in zip(score, aff)]
    return pd.DataFrame({'identifier': d['identifier'], 'docking_ID': d['ID'].astype(int),
                         'CNNscore': [float(x) for x in score], 'CNNaffinity': [float(x) for x in aff],
                         'CNN_VS': [float(x) for x in vs]})


sirt3, syk = scores(args.sirt3_sdf), scores(args.syk_sdf)
if not (sirt3['identifier'].equals(syk['identifier']) and sirt3['docking_ID'].equals(syk['docking_ID'])):
    raise SystemExit('The two docked.sdf files do not list the same compounds in the same order.')
print(f'Compounds in each docked.sdf: {len(sirt3):,}')

flag = f'selected_CNN_VS_ge_{args.cutoff:.1f}_both_targets'
table = pd.DataFrame({'docking_ID': sirt3['docking_ID'], 'identifier': sirt3['identifier']})
for name, d in (('SIRT3', sirt3), ('SYK', syk)):
    for col in ('CNNscore', 'CNNaffinity', 'CNN_VS'):
        table[f'{name}_{col}'] = d[col]
table[flag] = (table['SIRT3_CNN_VS'] >= args.cutoff) & (table['SYK_CNN_VS'] >= args.cutoff)
print(f'Compounds with CNN_VS >= {args.cutoff:.1f} on both targets: {int(table[flag].sum())}')

cons = pd.read_csv(args.consensus)
sel = table[table[flag]].merge(cons[['identifier', 'docking_ID', 'SYK_probability', 'SIRT3_probability']],
                               on=['identifier', 'docking_ID'])
sel = sel.rename(columns={'docking_ID': 'Molecule', 'SYK_probability': 'SYK_QSAR_prob',
                          'SIRT3_probability': 'SIRT3_QSAR_prob'})
sel = sel[['Molecule', 'identifier', 'SYK_QSAR_prob', 'SIRT3_QSAR_prob', 'SIRT3_CNNscore', 'SIRT3_CNNaffinity',
           'SIRT3_CNN_VS', 'SYK_CNNscore', 'SYK_CNNaffinity', 'SYK_CNN_VS']]

out = Path(args.output_dir)
out.mkdir(parents=True, exist_ok=True)
table.to_csv(out / 'docking_scores_both_targets.csv', index=False)
sel.to_csv(out / 'selected_hits.csv', index=False)
print(f'Saved in: {out.resolve()}')
