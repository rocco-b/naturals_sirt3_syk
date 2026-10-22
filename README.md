# Scaffolds from the SuperNatural 3.0 database as potential dual SIRT3–SYK inhibitors: 3D-QSAR and docking data

This repository contains the data and scripts of the 3D-QSAR screening and docking steps described in the manuscript *"Scaffolds from the SuperNatural 3.0 Database as Potential Dual SIRT3–SYK Inhibitors: A 3D-QSAR, Docking, and Molecular Dynamics Study"* (R. Buccheri, L. E. M. Russo, A. Rescifina; Department of Drug and Health Sciences, University of Catania, Italy).

The study screened the 124,943 compounds of SuperNatural 3.0 with two classification models (one for SYK, one for SIRT3), docked the compounds predicted active on both targets with GNINA, and selected three compounds with CNN_VS ≥ 6.0 on both targets for molecular dynamics (MD) simulations.

## What is and is not included

Included:

- the training sets, the saved models, and the screening results of the two 3D-QSAR models;
- the 3D-QSAR consensus set (4,177 compounds) with the predicted probabilities for both targets;
- the receptor files, ligand inputs, GNINA output files and logs, and the score tables of the docking of the consensus set against SIRT3 and SYK;
- the scripts used to train the models and to screen the database, and the scripts that rebuild the input table of the screening, the consensus set, and the docking score tables from the files provided;
- `requirements.txt`, with the versions of the Python packages recorded in the study.

Not included:

- **MD simulation files.** They are large and are available from the corresponding author on request (see Contact).
- **The original SuperNatural 3.0 SDF file** (about 680 MB, above the GitHub file size limit). It can be downloaded from the COCONUT portal (link in `1_3D-QSAR/supernatural3_database/source_info.txt`). The 65,564 compounds that passed the Lipinski filters are provided here.

## Repository layout

```
.
├── README.md
├── LICENSE
├── requirements.txt
├── 1_3D-QSAR/
│   ├── scripts/
│   │   ├── preparation/     make_cleaned_db.py
│   │   ├── training/        run_3d_qsar_RF.py, run_3d_qsar_SVM.py, clean_duplicates_csv.py
│   │   ├── screening/       predict_3d_qsar_RF.py, predict_3d_qsar_SVM.py
│   │   └── consensus/       make_consensus_set.py
│   ├── SYK/                 Random Forest model
│   │   ├── 1_chembl_raw/
│   │   ├── 2_training_dataset/
│   │   ├── 3_model/
│   │   └── 4_screening/
│   ├── SIRT3/               Support Vector Machine model (SMOTE)
│   │   ├── 1_chembl_raw/
│   │   ├── 2_training_dataset/
│   │   ├── 3_model/
│   │   └── 4_screening/
│   ├── supernatural3_database/
│   └── consensus_set_with_probabilities.csv
└── 2_Docking/
    ├── scripts/             extract_docking_scores.py
    ├── consensus_set_for_docking/
    ├── SIRT3/
    ├── SYK/
    └── results/
```

## Workflow and compound counts

1. SuperNatural 3.0: 124,943 compounds; after the Lipinski filters (molecular weight ≤ 500 Da, cLogP ≤ 5, hydrogen-bond donors ≤ 5, hydrogen-bond acceptors ≤ 10, applied to the columns `Molweight`, `cLogP`, `H-Donors`, and `H-Acceptors` of the annotated table; see `scripts/preparation/make_cleaned_db.py`): 65,564 compounds.
2. 3D-QSAR screening (3D autocorrelation descriptors computed from one RDKit conformer per molecule): 64,341 compounds scored by the SYK model (37,026 predicted active, probability ≥ 0.647) and 64,265 by the SIRT3 model (4,576 predicted active, probability ≥ 0.750).
3. Consensus set (predicted active on both targets): 4,177 compounds.
4. Docking of the 4,177 compounds with GNINA 1.3.3 against SIRT3 (PDB 4JSR) and SYK (PDB 3EMG).
5. Selection: CNN_VS ≥ 6.0 on both targets, which gave three compounds.

The three selected compounds (numbers are the docking IDs, see "Compound numbering" below):

- 2314 = CNP0067322.0
- 3843 = CNP0054808.1
- 4109 = CNP0524408.1

## Contents of the folders

### 1_3D-QSAR/

`supernatural3_database/`

- `supernatural3_cleaned_db.csv`: tab-separated, columns `identifier` and `Smiles`; the 65,564 compounds that passed the Lipinski filters. It is the input of the screening scripts.
- `supernatural3_Lipinski.txt.zip`: the same 65,564 compounds with the annotations of the downloaded table (tab-separated; 50 columns, among them `identifier`, `canonical_smiles`, `standard_inchi_key`, `molecular_weight`, `alogp`, `np_likeness`, `np_classifier_pathway`, `np_classifier_superclass`, `np_classifier_class`, `np_classifier_is_glycoside`; the last four columns, `Molweight`, `cLogP`, `H-Acceptors`, and `H-Donors`, are the ones used for the Lipinski filters). It is compressed to keep the file size low. `scripts/preparation/make_cleaned_db.py` applies the filters to this table and rebuilds `supernatural3_cleaned_db.csv`.
- `source_info.txt`: download address and reference of the database.

`SYK/` and `SIRT3/` (same structure; the SYK model is a Random Forest, the SIRT3 model a Support Vector Machine trained after SMOTE oversampling)

- `1_chembl_raw/`: the bioactivity data downloaded from the ChEMBL web portal (`DOWNLOAD-….csv`, semicolon-separated, IC50 data of CHEMBL2599 for SYK and CHEMBL4461 for SIRT3), the query description (`_query.txt`), and `dataset.xlsx`, the spreadsheet version of the labeled data (same columns as `dataset.csv`).
- `2_training_dataset/`: `dataset.csv` (tab-separated; columns `Molecule ChEMBL ID`, `Smiles`, `Standard Type`, `Standard Relation`, `Standard Value`, `Standard Units`, `activity_label`), where `activity_label` is 1 for IC50 < 1000 nM and 0 for IC50 ≥ 1000 nM; and `dataset_unique.csv` (comma-separated), the same table after the removal of duplicate SMILES strings with `scripts/training/clean_duplicates_csv.py`. It is the input of the training scripts. SYK: 5,590 molecules (4,825 active, 765 inactive); SIRT3: 482 molecules (76 active, 406 inactive).
- `3_model/`: the saved model (`3dqsar_rf_model.pkl` for SYK; `3dqsar_svm_model.pkl` for SIRT3). The files are Python pickles (joblib) of a dictionary: the trained model and the classification threshold (`optimal_threshold`), and for SIRT3 also the fitted scaler.
- `4_screening/`: `screening_results.csv` (SYK) and `svm_screening_results.csv` (SIRT3) with the columns `identifier`, `Smiles`, `Active_Probability`, `Predicted_Active` for all compounds that were scored; `top_candidates.csv` (SYK) and `svm_top_candidates.csv` (SIRT3) with the predicted actives only, sorted by probability.

`consensus_set_with_probabilities.csv`: comma-separated; the 4,177 compounds predicted active on both targets, with the columns `identifier`, `Smiles`, `docking_ID`, `SYK_probability`, `SIRT3_probability`. It is built from the predicted actives of the SIRT3 model (sorted by decreasing probability) that were also predicted active by the SYK model, with `scripts/consensus/make_consensus_set.py`.

### 2_Docking/

`consensus_set_for_docking/`

- `consensus_qsar.txt`: tab-separated; columns `identifier`, `Smiles`, `ID`. `ID` is the position of the compound in the consensus set (1–4,177).
- `ligands.smi`: the SMILES strings of the 4,177 compounds, one per line, in the same order.
- `consensus_qsar.sdf`: the consensus set as an SDF file before the Open Babel preparation (V3000 molfile with the 3D coordinates generated with DataWarrior).

`SIRT3/` and `SYK/`

- `protein.pdb`: the receptor structure used for docking.
- `consensus_qsar_pH.sdf`: the ligand input actually docked against the target: the 3D structures of the consensus set (generated with DataWarrior, as in `consensus_qsar.sdf`) with the protonation state set at pH 7.4 with Open Babel. The files of the two targets have the same identifiers and order, but not the same conformers (the headers record the dates of the Open Babel conversion: 2026-07-14 for SIRT3 and 2026-06-23 for SYK).
- `docked.sdf`: the GNINA output, one pose per compound (4,177 records), with the fields `identifier`, `Smiles`, `ID`, `MMFF94 Energy`, `minimizedAffinity`, `CNNscore`, `CNNaffinity`, `CNNaffinity_variance`, `CNN_VS`.
- `log.txt`: the complete GNINA log, which starts with the command line and reports the random seed (SIRT3: −605940204; SYK: 1233175373).
- The docked poses of the three selected compounds are the records of `docked.sdf` with `ID` 2314, 3843, and 4109.
- The GNINA command line used for each target is reported in the first lines of `log.txt` and, with relative paths, in the section "How to reproduce the analyses".

`results/`

- `docking_scores_both_targets.csv`: comma-separated; one row per compound with `docking_ID`, `identifier`, the three GNINA scores for SIRT3 (`SIRT3_CNNscore`, `SIRT3_CNNaffinity`, `SIRT3_CNN_VS`) and for SYK (`SYK_CNNscore`, `SYK_CNNaffinity`, `SYK_CNN_VS`), and a column that is true for the compounds with CNN_VS ≥ 6.0 on both targets.
- `selected_hits.csv`: the three selected compounds with the QSAR probabilities and the docking scores (Table 2 of the manuscript).
- Both tables are produced from the two `docked.sdf` files and `consensus_set_with_probabilities.csv` by `2_Docking/scripts/extract_docking_scores.py`. In the tables, CNNscore is rounded to 4 decimals, CNNaffinity to 3 decimals, and CNN_VS is the product of the two rounded values (4 decimals); the CNN_VS stored in `docked.sdf` is the unrounded product and differs from it by at most 0.0008, without changing the selection.

## Compound numbering

The compound numbers used in the manuscript (2314, 3843, 4109) are the positions of the compounds in the consensus set (`ID` in `consensus_qsar.txt`, `docking_ID` in the other tables). The stable identifiers are the SuperNatural 3.0 identifiers (`CNP…`).

## How to reproduce the analyses

Software used (as reported in the manuscript): RDKit 2025.09.1, scikit-learn 1.8.0, imbalanced-learn 0.14.2 (SMOTE), Open Babel 3.1.1, GNINA 1.3.3 (master:6fe1ce2, built on 2026-06-30). The versions of pandas, numpy, and joblib were not recorded. The Python packages are listed in `requirements.txt` (`pip install -r requirements.txt`); Open Babel and GNINA are separate programs. Python scripts print their results in the terminal.

**Input table of the screening.** The script applies the Lipinski filters to the annotated table of SuperNatural 3.0 and writes `supernatural3_cleaned_db.csv`. With the table provided here (65,564 compounds, already filtered) it returns all of them and reproduces the file of the repository:

```
python <repository>/1_3D-QSAR/scripts/preparation/make_cleaned_db.py --output supernatural3_cleaned_db.csv
```

**Training of a model.** The scripts read `dataset_unique.csv` from the current folder and save the model in the current folder. Run them in a scratch copy of the dataset folder to avoid overwriting `3_model/`:

```
# SYK: run from a folder that contains dataset_unique.csv
python <repository>/1_3D-QSAR/scripts/training/run_3d_qsar_RF.py

# SIRT3: run from a folder that contains dataset_unique.csv
python <repository>/1_3D-QSAR/scripts/training/run_3d_qsar_SVM.py
```

**Screening of the database.** The scripts read the model (`3dqsar_rf_model.pkl` or `3dqsar_svm_model.pkl`) and `supernatural3_cleaned_db.csv` from the current folder and write the two result files in the same folder:

```
# SYK: run from a folder that contains 3dqsar_rf_model.pkl and supernatural3_cleaned_db.csv
python <repository>/1_3D-QSAR/scripts/screening/predict_3d_qsar_RF.py

# SIRT3: run from a folder that contains 3dqsar_svm_model.pkl and supernatural3_cleaned_db.csv
python <repository>/1_3D-QSAR/scripts/screening/predict_3d_qsar_SVM.py
```

**Consensus set.** The compounds predicted active by both models, in the order used for docking (the position is the docking ID):

```
python <repository>/1_3D-QSAR/scripts/consensus/make_consensus_set.py --output-dir <folder>
```

It writes `consensus_set_with_probabilities.csv`, `consensus_qsar.txt`, and `ligands.smi`, identical to the files of the repository. The 3D structures of the consensus set (`consensus_qsar.sdf`) were generated with DataWarrior (MMFF94 force field) and protonated at pH 7.4 with Open Babel (`obabel consensus_qsar.sdf -O consensus_qsar_pH.sdf -p 7.4`).

**Docking.** The commands recorded in the logs (local paths replaced by relative ones) are:

```
# SIRT3
gnina -r protein.pdb -l consensus_qsar_pH.sdf --center_x 26.037 --center_y 46.249 --center_z -1.103 --size_x 18 --size_y 18 --size_z 20 --cnn_scoring rescore --cnn_rotation 15 --num_modes 1 --log log.txt -o docked.sdf

# SYK
gnina -r protein.pdb -l consensus_qsar_pH.sdf --center_x 17.324 --center_y 40.93 --center_z 9.178 --size_x 22 --size_y 30 --size_z 26 --cnn_scoring rescore --cnn_rotation 15 --num_modes 1 --log log.txt -o docked.sdf
```

**Docking scores and selection of the three compounds.** The script reads the two `docked.sdf` files and writes the score table and the selected compounds (CNN_VS ≥ 6.0 on both targets, changeable with `--cutoff`):

```
python <repository>/2_Docking/scripts/extract_docking_scores.py --output-dir <folder>
```

The output is identical to `docking_scores_both_targets.csv` and `selected_hits.csv`, except for 25 of the 8,354 CNN_VS values, which differ by 0.0001 in the last digit (rounding of products that end exactly in 5).

## Notes and caveats

- The 3D embedding of the molecules (RDKit ETKDGv3, with a fallback to random coordinates) was run without a fixed random seed. Re-running the training scripts therefore gives slightly different metrics and thresholds, and a re-run of the screening can score a slightly different set of molecules. The results reported in the manuscript are those of the saved models in `3_model/` and of the files in `4_screening/` (99 compounds were scored only by the SYK model and 23 only by the SIRT3 model). The SIRT3 screening script discards the molecules whose descriptors contain non-finite values; the SYK script does not.
- The Lipinski filters of the study were applied to the columns `Molweight`, `cLogP`, `H-Acceptors`, and `H-Donors` of the annotated table, with inclusive limits (the table contains compounds with 10 acceptors, 5 donors, and a cLogP of 5.0). The columns `molecular_weight`, `alogp`, and `hydrogen_bond_*` of the original SDF are different calculated values: the same limits applied to them select 61,234 compounds, 60,140 of which are among the 65,564 of the study. The original SDF alone therefore does not reproduce the filtered set.
- In the training scripts, the 5-fold cross-validation of the SIRT3 model is run on the SMOTE-resampled training set, and the classification thresholds are derived from the predicted probabilities of the training set (see the Methods of the manuscript).
- The pickled models were saved with scikit-learn 1.8.0. Loading them with another version produces a warning and is not guaranteed to work. Load pickle files only from sources you trust.
- One SYK entry (CHEMBL5782610) has no SMILES string, so 5,589 of the 5,590 SYK molecules were used for modeling.
- Curation of the training sets: `clean_duplicates_csv.py` keeps the first occurrence of each SMILES string in a file ordered by increasing IC50, that is, the most potent measurement. In the SYK set, 56 structures with measurements on both sides of 1000 nM are therefore labeled active, and three records with the relation `<` and a value of 1000 nM (CHEMBL5892954, CHEMBL1078178, CHEMBL5958423) are labeled inactive. The SIRT3 set has no structure with measurements on both sides of 1000 nM. The descriptors of the training sets were not saved: they are computed by the training scripts.
- The structural zinc ion of SIRT3 (coordinated by Cys256, Cys259, Cys280, and Cys283) is not present in `2_Docking/SIRT3/protein.pdb`, and it was not included in the MD systems.
- Some files contain the local paths of the computer on which the calculations were run (the first lines of `log.txt`).
- The SYK receptor used for docking is the crystal structure (PDB 3EMG, retrieved from PDB-REDO); the MD simulations of SYK started from a full-length homology model, which is not part of this repository. The SIRT3 receptor is derived from PDB 4JSR (PDB-REDO).

## Data sources and terms of use

- Bioactivity data: ChEMBL (https://www.ebi.ac.uk/chembl/), targets CHEMBL2599 (SYK) and CHEMBL4461 (SIRT3), standard type IC50, downloaded from the web portal on June 22, 2026 (the downloaded files are in `1_chembl_raw/`).
- Natural products: SuperNatural 3.0 (Gallo et al., Nucleic Acids Res. 2023, 51, D654–D659, https://doi.org/10.1093/nar/gkac1008), obtained through the COCONUT portal (https://coconut.naturalproducts.net/).
- Protein structures: PDB entries 4JSR (SIRT3) and 3EMG (SYK), via PDB-REDO.
- The files derived from these resources are subject to the terms of use of the original databases.

## Citation

The citation of the article will be added when its DOI is available.

## License

This project is released under the MIT License. See [LICENSE](LICENSE).

## Contact

Corresponding author: Antonio Rescifina, Department of Drug and Health Sciences, University of Catania, Viale A. Doria 6, 95125 Catania, Italy; antonio.rescifina@unict.it. The MD simulation files can be requested by email.
