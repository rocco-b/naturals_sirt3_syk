import pandas as pd
import numpy as np
import joblib
import sys
import math
import multiprocessing
from joblib import Parallel, delayed
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit import RDLogger

# Disable RDKit warnings for cleaner output
RDLogger.DisableLog('rdApp.*')

# =====================================================================
# CONFIGURATION
# =====================================================================
MODEL_FILE = '3dqsar_svm_model.pkl'       # Make sure this matches your saved SVM model file
INPUT_CSV = 'supernatural3_cleaned_db.csv'       # Change this to your new dataset file
OUTPUT_ALL_CSV = 'svm_screening_results.csv' # Contains all successfully processed molecules
OUTPUT_HITS_CSV = 'svm_top_candidates.csv'   # Contains ONLY the predicted active hits

# =====================================================================
# 1. LOAD THE VALIDATED SVM MODEL, SCALER, AND THRESHOLD
# =====================================================================
print("Loading the saved 3D QSAR SVM model pipeline...")
try:
    model_package = joblib.load(MODEL_FILE)
    model = model_package['trained_model']
    scaler = model_package['scaler'] # CRITICAL: SVM requires feature scaling
    optimal_threshold = model_package['optimal_threshold']
    
    print(f"✅ Pipeline loaded successfully!")
    print(f"   - Model: Support Vector Classifier (RBF)")
    print(f"   - Scaler: StandardScaler found")
    print(f"   - Optimal Recall Threshold: {optimal_threshold:.3f}\n")
except FileNotFoundError:
    print(f"❌ Error: Model file '{MODEL_FILE}' not found. Please ensure it is in the same directory.")
    sys.exit()
except KeyError as e:
    print(f"❌ Error: The loaded model package is missing a required component: {e}")
    print("Please ensure you are loading the correct SVM .pkl file that contains the scaler.")
    sys.exit()

# =====================================================================
# 2. LOAD THE UNKNOWN MOLECULES DATASET
# =====================================================================
print(f"Loading input dataset: '{INPUT_CSV}'...")
try:
    df_new = pd.read_csv(INPUT_CSV, sep='\t')
    df_new.columns = df_new.columns.str.strip() # Remove accidental whitespaces
    
    if 'Smiles' not in df_new.columns:
        print("❌ Error: The input CSV must contain a column named exactly 'Smiles'.")
        print(f"Columns found: {df_new.columns.tolist()}")
        sys.exit()
        
    print(f"✅ Dataset loaded! Found {len(df_new)} molecules to screen.\n")
except FileNotFoundError:
    print(f"❌ Error: Input file '{INPUT_CSV}' not found.")
    sys.exit()

# =====================================================================
# 3. 3D DESCRIPTOR PIPELINE (PARALLELIZED & CPU-SAFE)
# =====================================================================
def process_single_molecule(i, smiles):
    """Isolated function to process a single molecule for multiprocessing."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (i, None)
        
        # Add explicit hydrogens
        mol = Chem.AddHs(mol)
        
        # 3D Embedding
        status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if status < 0:
            status = AllChem.EmbedMolecule(mol, useRandomCoords=True)
        
        if status >= 0:
            # Force Field Optimization
            AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
            
            # Extract 3D Autocorrelation Descriptors
            features = rdMolDescriptors.CalcAUTOCORR3D(mol)
            
            # --- MODIFICA CRITICA: Controlla la presenza di NaN o Inf ---
            features_array = np.array(features)
            if np.isnan(features_array).any() or np.isinf(features_array).any():
                return (i, None) # Scarta la molecola se i descrittori sono corrotti
            
            return (i, features)
            
    except Exception:
        # Skip molecule if chemistry engine fails
        pass
        
    return (i, None)

print(f"[1/3] Generating 3D coordinates and minimizing energy for {len(df_new)} molecules...")

# Calculate safe CPU usage (70% max to prevent system freeze)
import multiprocessing
from joblib import Parallel, delayed

n_cores_total = multiprocessing.cpu_count()
n_cores_safe = max(1, int(n_cores_total * 0.70))

print(f"      ➔ Found {n_cores_total} total CPU cores.")
print(f"      ➔ Multiprocessing activated: Using {n_cores_safe} cores (70% limit).")
print(f"      ➔ (The remaining 30% is kept free to maintain PC responsiveness)")

smiles_list = df_new['Smiles'].tolist()

# Execute parallel processing
# verbose=10 shows a progress tracker in the terminal
parallel_results = Parallel(n_jobs=n_cores_safe, verbose=10)(
    delayed(process_single_molecule)(i, smiles) for i, smiles in enumerate(smiles_list)
)

# Reconstruct the results maintaining the correct order
total_descriptors = []
valid_indices = []

for index, features in parallel_results:
    if features is not None:
        total_descriptors.append(features)
        valid_indices.append(index)

X_new_3d = np.array(total_descriptors)

print(f"\n➔ 3D processing complete. {len(X_new_3d)} out of {len(df_new)} molecules successfully converted.\n")

if len(X_new_3d) == 0:
    print("❌ Critical Error: No molecules could be converted to 3D. Exiting.")
    sys.exit()

# Filter the original dataframe to keep only the molecules that were successfully processed
df_processed = df_new.iloc[valid_indices].copy()

# =====================================================================
# 4. DATA SCALING, PREDICTION AND THRESHOLD APPLICATION
# =====================================================================
print("[2/3] Scaling descriptors and predicting activity probabilities...")

# CRITICAL STEP FOR SVM: Apply the exact same scaling transformation used during training
print("      Applying Z-score normalization (StandardScaler)...")
X_new_scaled = scaler.transform(X_new_3d)

# Predict the probability of being active (Class 1) using the scaled data
y_pred_proba = model.predict_proba(X_new_scaled)[:, 1]

# Apply the mathematically optimized threshold saved from the training phase
print(f"      Applying the scientific cutoff (Probability >= {optimal_threshold:.3f})...")
y_pred_binary = (y_pred_proba >= optimal_threshold).astype(int)

# Add the results as new columns to the dataframe
df_processed['Active_Probability'] = np.round(y_pred_proba, 4)
df_processed['Predicted_Active'] = y_pred_binary

# =====================================================================
# 5. EXPORT RESULTS TO CSV
# =====================================================================
print("\n[3/3] Exporting results to CSV files...")

# Save the full results (all successfully processed molecules)
df_processed.to_csv(OUTPUT_ALL_CSV, index=False)
print(f"💾 Saved comprehensive results to: '{OUTPUT_ALL_CSV}'")

# Filter to keep ONLY the predicted active candidates (Hits)
df_hits = df_processed[df_processed['Predicted_Active'] == 1].copy()

# Sort the hits by highest probability first (best candidates at the top)
df_hits = df_hits.sort_values(by='Active_Probability', ascending=False)

# Save the hits for docking
df_hits.to_csv(OUTPUT_HITS_CSV, index=False)
print(f"💾 Saved TOP CANDIDATES ({len(df_hits)} molecules) to: '{OUTPUT_HITS_CSV}'")

print("\n" + "="*50)
print(" 🎉 SVM VIRTUAL SCREENING COMPLETE! ")
print("="*50)
print(f" Initial molecules  : {len(df_new)}")
print(f" 3D Processed       : {len(X_new_3d)}")
print(f" Predicted Hits     : {len(df_hits)} (Ready for molecular docking)")
print("="*50)
