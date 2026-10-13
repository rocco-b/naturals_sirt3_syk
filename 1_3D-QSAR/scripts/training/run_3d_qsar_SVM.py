import pandas as pd
import numpy as np
import sys
import joblib
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (roc_auc_score, matthews_corrcoef, f1_score, 
                             accuracy_score, confusion_matrix, roc_curve)
from imblearn.over_sampling import SMOTE
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit import RDLogger

# Disable RDKit warnings for cleaner output
RDLogger.DisableLog('rdApp.*')

# =====================================================================
# 1. LOAD THE CSV DATASET
# =====================================================================
csv_file_name = 'dataset_unique.csv'  # Change this to your actual SIRT3 file

try:
    df = pd.read_csv(csv_file_name, sep=',')
    df.columns = df.columns.str.strip() 
    print(f"✅ Dataset successfully loaded! Found {len(df)} molecules.")
except FileNotFoundError:
    print(f"❌ Error: File '{csv_file_name}' not found.")
    sys.exit()

# =====================================================================
# 2. INTERNAL PIPELINE: CONVERT SMILES TO 3D DESCRIPTORS
# =====================================================================
def generate_3d_descriptors(smiles_list):
    total_descriptors = []
    valid_indices = []
    
    print("\n[1/4] Generating 3D coordinates and minimizing energy...")
    for i, smiles in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None: continue
            
            mol = Chem.AddHs(mol)
            status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if status < 0:
                status = AllChem.EmbedMolecule(mol, useRandomCoords=True)
            
            if status >= 0:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
                features = rdMolDescriptors.CalcAUTOCORR3D(mol)
                total_descriptors.append(features)
                valid_indices.append(i)
        except Exception:
            continue
            
    return np.array(total_descriptors), valid_indices

X_3d, valid_indices = generate_3d_descriptors(df['Smiles'].tolist())
y = df['activity_label'].iloc[valid_indices].values

print(f"➔ Processing complete. {len(X_3d)} molecules successfully converted.")
if len(X_3d) < 10:
    print("❌ Critical Error: Too few molecules converted.")
    sys.exit()

# =====================================================================
# 3. SPLIT & SYNTHETIC OVERSAMPLING (SMOTE)
# =====================================================================
# Split data while maintaining the real-world imbalance
X_train, X_test, y_train, y_test = train_test_split(
    X_3d, y, test_size=0.20, stratify=y, random_state=42
)

print("\n[2/4] Applying SMOTE (Synthetic Minority Over-sampling) to Training Set...")
print(f"      Original Training Set: {sum(y_train==1)} Actives, {sum(y_train==0)} Inactives")

# Apply SMOTE strictly to the training set to prevent data leakage
smote = SMOTE(random_state=42)
X_train_resampled, y_train_resampled = smote.fit_resample(X_train, y_train)

print(f"      Resampled Training Set: {sum(y_train_resampled==1)} Actives, {sum(y_train_resampled==0)} Inactives")

# =====================================================================
# 4. FEATURE SCALING (CRITICAL FOR SVM)
# =====================================================================
# SVM algorithms are highly sensitive to feature magnitudes.
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_resampled)
X_test_scaled = scaler.transform(X_test) # Scale test set using training parameters

# =====================================================================
# 5. TRAINING SUPPORT VECTOR MACHINE (SVM)
# =====================================================================
print("\n[3/4] Training Support Vector Machine (RBF Kernel)...")

# probability=True is required to calculate ROC-AUC and custom thresholds
model = SVC(kernel='rbf', probability=True, class_weight='balanced', random_state=42)

# Cross-Validation on the resampled training set
cv_scores = cross_val_score(model, X_train_scaled, y_train_resampled, cv=5, scoring='roc_auc')
cv_auc_mean = cv_scores.mean()

# Final fit
model.fit(X_train_scaled, y_train_resampled)

# =====================================================================
# 6. RECALL-ORIENTED THRESHOLD OPTIMIZATION
# =====================================================================
print("\n[4/4] Calculating Recall-Oriented scientific threshold...")

# Get probabilities on the resampled training set
y_train_proba = model.predict_proba(X_train_scaled)[:, 1]
fpr, tpr, thresholds = roc_curve(y_train_resampled, y_train_proba)

# GOAL: Find the strictest threshold that still captures at least 85% of active hits
TARGET_SENSITIVITY = 0.85
acceptable_indices = np.where(tpr >= TARGET_SENSITIVITY)[0]

if len(acceptable_indices) > 0:
    # thresholds array is sorted descending. We take the first index (highest possible 
    # threshold) that guarantees >= 85% True Positive Rate.
    best_idx = acceptable_indices[0]
    optimal_threshold = thresholds[best_idx]
else:
    optimal_threshold = 0.50
    print("      Warning: Target sensitivity not reachable. Defaulting to 0.50")

# =====================================================================
# 7. EXTERNAL PREDICTION ON THE HIDDEN TEST SET
# =====================================================================
# Get probabilities for the unseen, un-resampled, real-world Test Set
y_test_pred_proba = model.predict_proba(X_test_scaled)[:, 1]

# Apply our custom Recall-oriented threshold
y_test_pred_optimized = (y_test_pred_proba >= optimal_threshold).astype(int)

# Calculate metrics
test_roc_auc = roc_auc_score(y_test, y_test_pred_proba)
test_mcc = matthews_corrcoef(y_test, y_test_pred_optimized)
test_f1 = f1_score(y_test, y_test_pred_optimized)
test_acc = accuracy_score(y_test, y_test_pred_optimized)
tn, fp, fn, tp = confusion_matrix(y_test, y_test_pred_optimized).ravel()

# =====================================================================
# 8. FINAL VALIDATION REPORT
# =====================================================================
print("\n" + "="*55)
print("         FINAL SVM MODEL VALIDATION REPORT")
print("="*55)
print(f" Total molecules in Training Set : {len(X_train_resampled)} (SMOTE applied)")
print(f" Total molecules in Hidden Test Set : {len(X_test)} (Real distribution)")
print(f" Calculated Recall Threshold (>= 85% TPR): {optimal_threshold:.3f}")
print(f" Confusion Matrix (Test Set) : TP={tp}, TN={tn}, FP={fp}, FN={fn}")
print("-" * 55)
print(f" 📊 Internal CV ROC-AUC    : {cv_auc_mean:.3f}  (Minimum required: > 0.70)")
print(f" 📊 External Test ROC-AUC  : {test_roc_auc:.3f}  (Overall discriminative power)")
print(f" 📊 External Test MCC      : {test_mcc:.3f}  (Minimum required: > 0.40)")
print(f" 📊 External Test F1-Score : {test_f1:.3f}  (Precision/Recall balance)")
print(f" 📊 Overall Accuracy       : {test_acc:.3f}")
print("="*55)

# =====================================================================
# 9. SAVE MODEL & SCALER
# =====================================================================
# We MUST save the scaler along with the model and threshold, 
# otherwise future unknown molecules cannot be processed.
model_packing = {
    'trained_model': model,
    'scaler': scaler,
    'optimal_threshold': optimal_threshold
}
        
model_file_name = '3dqsar_svm_model.pkl'
joblib.dump(model_packing, model_file_name)
      
print(f"\n💾 Saved model, scaler, and threshold parameters to: '{model_file_name}'")