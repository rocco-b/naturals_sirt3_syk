import pandas as pd
import numpy as np
import sys
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (roc_auc_score, matthews_corrcoef, f1_score, 
                             accuracy_score, confusion_matrix, roc_curve)
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

# =====================================================================
# 1. LOAD THE CSV DATASET
# =====================================================================
csv_file_name = 'dataset_unique.csv'  # Change this to your actual file path

try:
    df = pd.read_csv(csv_file_name, sep=',')
    # Strip invisible whitespaces from column names just in case
    df.columns = df.columns.str.strip() 
    print(f"✅ Dataset successfully loaded! Found {len(df)} molecules.")
except FileNotFoundError:
    print(f"❌ Error: File '{csv_file_name}' not found.")
    print("Loading a temporary mock dataset with your structure for demonstration...")
    # Mock dataset based on your provided structure
    mock_data = {
        'Molecule ChEMBL ID': ['CHEMBL2332041', 'CHEMBL0000002', 'CHEMBL0000003', 'CHEMBL0000004'],
        'Smiles': [
            'CC(C)(C)OC(=O)c1ccc(C(=O)NCCC2CCN(c3ncnc4cc(C(N)=O)sc34)CC2)s1',
            'CCN(CC)CCNc1ccc(C(=O)Nc2ccc(C)cc2)cc1',
            'O=C(Nc1ccc(Cl)cc1)c1ccccc1O',
            'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1'
        ],
        'Standard Type': ['IC50', 'IC50', 'IC50', 'IC50'],
        'Standard Relation': ['=', '=', '=', '='],
        'Standard Value': [3.2, 1500.0, 45.0, 8500.0],
        'Standard Units': ['nM', 'nM', 'nM', 'nM'],
        'activity_label': [1, 0, 1, 0]
    }
    # Artificially duplicate the data to allow for proper train/test splitting
    df = pd.DataFrame(mock_data)
    df = pd.concat([df] * 15, ignore_index=True) 

# =====================================================================
# 2. INTERNAL PIPELINE: CONVERT SMILES TO 3D DESCRIPTORS
# =====================================================================
def generate_3d_descriptors(smiles_list):
    total_descriptors = []
    valid_indices = []
    
    print("\n[1/3] Generating 3D coordinates and minimizing energy...")
    for i, smiles in enumerate(smiles_list):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            # Add explicit hydrogens (crucial for accurate 3D geometry)
            mol = Chem.AddHs(mol)
            
            # Generate initial 3D conformation (ETKDG v3 algorithm)
            status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            if status < 0:
                # Fallback to random coordinates if the molecule is too rigid/complex
                status = AllChem.EmbedMolecule(mol, useRandomCoords=True)
            
            if status >= 0:
                # Optimize 3D geometry (MMFF94 Force Field)
                AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
                
                # Extract 3D Autocorrelation Descriptors (geometric 3D vector)
                features = rdMolDescriptors.CalcAUTOCORR3D(mol)
                
                total_descriptors.append(features)
                valid_indices.append(i)
        except Exception:
            continue
            
    return np.array(total_descriptors), valid_indices

# Execute the spatial conversion
X_3d, valid_indices = generate_3d_descriptors(df['Smiles'].tolist())

# Align targets (0/1) by excluding molecules that failed the 3D conversion
y = df['activity_label'].iloc[valid_indices].values
print(f"➔ Processing complete. {len(X_3d)} molecules successfully converted into 3D geometric vectors.")

if len(X_3d) < 10:
    print("❌ Critical Error: Too few molecules converted to 3D to train a model. Check your SMILES codes.")
    sys.exit()

# =====================================================================
# 3. AUTOMATIC STRATIFIED SPLIT (80% Training, 20% Test)
# =====================================================================
# 'stratify=y' ensures the active/inactive ratio remains identical in both sets.
X_train, X_test, y_train, y_test = train_test_split(
    X_3d, y, test_size=0.20, stratify=y, random_state=42
)

# =====================================================================
# 4. TRAINING AND INTERNAL VALIDATION (Cross-Validation)
# =====================================================================
# 'class_weight=balanced' handles imbalanced ChEMBL datasets automatically
model = RandomForestClassifier(n_estimators=150, random_state=42, class_weight='balanced')

print("\n[2/3] Calculating 5-Fold Internal Cross-Validation...")
cv_scores = cross_val_score(model, X_train, y_train, cv=5, scoring='roc_auc')
cv_auc_mean = cv_scores.mean()

# Final training on the entire Training Set
model.fit(X_train, y_train)

# =====================================================================
# 5. MATHEMATICAL THRESHOLD OPTIMIZATION (Youden's J)
# =====================================================================
print("\n[3/3] Calculating optimal scientific threshold and evaluating Test Set...")

# 5a. Get probabilities on the TRAINING set to find the best threshold
y_train_proba = model.predict_proba(X_train)[:, 1]
fpr, tpr, thresholds = roc_curve(y_train, y_train_proba)

# 5b. Calculate Youden's J statistic for all thresholds to find the optimum
youden_index = tpr - fpr
best_threshold_idx = np.argmax(youden_index)
optimal_threshold = thresholds[best_threshold_idx]

# =====================================================================
# 6. EXTERNAL PREDICTION ON TEST SET USING OPTIMAL THRESHOLD
# =====================================================================
# Get probabilities for the unseen Test Set
y_test_pred_proba = model.predict_proba(X_test)[:, 1]

# Apply the mathematically derived threshold instead of the default 0.50
y_test_pred_optimized = (y_test_pred_proba >= optimal_threshold).astype(int)

# Calculate classification metrics
test_roc_auc = roc_auc_score(y_test, y_test_pred_proba) # ROC-AUC is threshold-independent
test_mcc = matthews_corrcoef(y_test, y_test_pred_optimized)
test_f1 = f1_score(y_test, y_test_pred_optimized)
test_acc = accuracy_score(y_test, y_test_pred_optimized)
tn, fp, fn, tp = confusion_matrix(y_test, y_test_pred_optimized).ravel()

# =====================================================================
# 7. FINAL VALIDATION REPORT AND AUTOMATED DIAGNOSIS
# =====================================================================
print("\n" + "="*55)
print("         FINAL MODEL VALIDATION REPORT")
print("="*55)
print(f" Total molecules in Training Set : {len(X_train)}")
print(f" Total molecules in Hidden Test Set : {len(X_test)}")
print(f" Calculated Optimal Threshold (Youden's J) : {optimal_threshold:.3f}")
print(f" Confusion Matrix (Test Set) : TP={tp}, TN={tn}, FP={fp}, FN={fn}")
print("-" * 55)
print(f" 📊 Internal CV ROC-AUC    : {cv_auc_mean:.3f}  (Minimum required: > 0.70)")
print(f" 📊 External Test ROC-AUC  : {test_roc_auc:.3f}  (Overall discriminative power)")
print(f" 📊 External Test MCC      : {test_mcc:.3f}  (Minimum required: > 0.40)")
print(f" 📊 External Test F1-Score : {test_f1:.3f}  (Precision/Recall balance)")
print(f" 📊 Overall Accuracy       : {test_acc:.3f}")
print("="*55)

print("\n🤖 AUTOMATED SCREENING FUNNEL DIAGNOSIS:")
THRESHOLD_AUC = 0.70
THRESHOLD_MCC = 0.40
OVERFITTING_TOLERANCE = 0.15
cv_test_difference = cv_auc_mean - test_roc_auc

if cv_test_difference > OVERFITTING_TOLERANCE:
    print("❌ VERDICT: MODEL IS OVERFITTING.")
    print(f"Reason: Internal performance ({cv_auc_mean:.2f}) is significantly higher than external performance ({test_roc_auc:.2f}).")
    print("⚠️ DO NOT USE: You risk losing real active compounds or selecting too many false positives.")
elif test_roc_auc < THRESHOLD_AUC or test_mcc < THRESHOLD_MCC:
    print("⚠️ VERDICT: INSUFFICIENT PERFORMANCE (Underfitting).")
    print("Reason: Test Set metrics are below the minimum acceptable reliability criteria.")
    print("⚠️ DO NOT USE: The extracted 3D descriptors fail to separate actives from inactives.")
else:
    print("✅ VERDICT: MODEL VALIDATED AND READY FOR USE!")
    print("Reason: Excellent balance between internal and external data. Solid MCC and ROC-AUC values.")
    print("🚀 STATUS: Safe to use as a preliminary filter before refined docking.")

# =====================================================================
# 8. BONUS: VIRTUAL SCREENING EXPORT LOGIC
# =====================================================================
# If the model is valid, you can simulate screening the Test Set 
# and saving only the predicted active compounds to a new CSV file.
if test_roc_auc >= THRESHOLD_AUC and test_mcc >= THRESHOLD_MCC:
    print("\n[Virtual Screening] Exporting predicted active candidates for docking...")
    
    # Create a DataFrame for the test set predictions
    test_df = df.iloc[valid_indices].iloc[y_test.argsort()][:len(y_test)].copy() # Align with test subset
    
    # Just an example of how to filter: we take indices where prediction is 1
    predicted_active_indices = np.where(y_test_pred_optimized == 1)[0]
    
    if len(predicted_active_indices) > 0:
        print(f"Found {len(predicted_active_indices)} potential hits based on the optimal threshold.")
        # In a real scenario, you'd apply this to an entirely new database without known labels.
        # df_new_hits.to_csv("ready_for_docking.csv", index=False)
        print("➔ 'ready_for_docking.csv' would be successfully generated.")
    else:
        print("No candidates met the threshold criteria.")
        
import joblib

# =====================================================================
# 9. SAVE MODEL
# =====================================================================
model_packing = {
    'trained_model': model,
    'optimal_threshold': optimal_threshold
}
        
model_file_name = '3dqsar_model.pkl'
joblib.dump(model_packing, model_file_name)
      
print(f"Saved model: '{model_file_name}'")
