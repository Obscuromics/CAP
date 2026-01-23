#!/usr/bin/env python3
import pandas as pd
import joblib
import os
import numpy as np
import sys
import warnings

# Suppress version compatibility warnings (they're non-fatal but noisy)
warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# ───────────────────────────────────────────────
# Inputs from command line arguments
# ───────────────────────────────────────────────
if len(sys.argv) != 4:
    sys.exit("Usage: predict_centromeric.py <input_csv> <model_pkl> <output_csv>")

input_csv = sys.argv[1]
model_path = sys.argv[2]
output_csv = sys.argv[3]

# ───────────────────────────────────────────────
# Load trained model artifacts
# ───────────────────────────────────────────────
if not os.path.exists(model_path):
    sys.exit(f"❌ Model file not found: {model_path}")

try:
    artifacts = joblib.load(model_path)
    best_model = artifacts['model']
    scaler = artifacts['scaler']
    selector = artifacts['selector']
    selected_features = artifacts['selected_features']
    numerical_cols = artifacts['numerical_cols']
    print(f"✅ Loaded model from {model_path}")
    print(f"Using {len(selected_features)} selected features.")
except Exception as e:
    sys.exit(f"❌ Error loading model: {e}\n"
             f"   The model may be incompatible with current XGBoost/scikit-learn versions.\n"
             f"   Please regenerate the model with current package versions.")

# ───────────────────────────────────────────────
# Load new genomic class data
# ───────────────────────────────────────────────
if not os.path.exists(input_csv):
    sys.exit(f"❌ Input data file not found: {input_csv}")

new_data = pd.read_csv(input_csv, na_values=['NA', ''], keep_default_na=True, low_memory=False)
print(f"📄 Loaded {len(new_data)} rows from {input_csv}")

# Check for empty data
if len(new_data) == 0:
    sys.exit(f"❌ Input CSV is empty: {input_csv}\n"
             f"   Please check that SCORE_CENTROMERIC process produced valid output.")

new_data = new_data.replace(-1, np.nan)

# ───────────────────────────────────────────────
# Preprocess (match training pipeline)
# ───────────────────────────────────────────────
columns_to_drop = ['class', 'X', 'new_class_num_ID', 'chromosome', 'is_cen']
new_data = new_data.drop(columns=[c for c in columns_to_drop if c in new_data.columns])

# Ensure all training columns are present
for col in numerical_cols:
    if col not in new_data.columns:
        new_data[col] = np.nan

# Impute missing numeric values (fallback to 0 if all NaN)
for col in numerical_cols:
    if new_data[col].isna().all():
        new_data[col] = 0
    else:
        new_data[col] = new_data[col].fillna(new_data[col].median())

# Scale safely
new_data[numerical_cols] = scaler.transform(new_data[numerical_cols])

# Align and sanitize (ensure no NaNs remain)
new_data = new_data.reindex(columns=numerical_cols, fill_value=0)
new_data = new_data.replace([np.inf, -np.inf], 0).fillna(0)

# Feature selection
new_selected = selector.transform(new_data)

# ───────────────────────────────────────────────
# Predict
# ───────────────────────────────────────────────
predictions = best_model.predict(new_selected)
probabilities = best_model.predict_proba(new_selected)[:, 1]

results = pd.DataFrame({
    'prediction_centromeric': predictions,
    'probability_centromeric': probabilities
})

# ───────────────────────────────────────────────
# Save results
# ───────────────────────────────────────────────
results.to_csv(output_csv, index=False)
print(f"💾 Saved predictions to {output_csv}")
