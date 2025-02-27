# Import Statements
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from xgboost import XGBClassifier
import shap
import sys
from xgboost import XGBClassifier
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV, StratifiedKFold, cross_validate
from sklearn.metrics import make_scorer, precision_score, recall_score, f1_score, average_precision_score

cwd = Path.cwd()
print(cwd)
datasets = cwd / '../results/tax_classification_out/abundance_matrices'
results = cwd / '../results/ML_out'

g1 = sys.argv[1]
g2 = sys.argv[2]
col_name = sys.argv[3]
threshold = sys.argv[4]

print(f'{g1}|{g2}|{col_name}|{threshold}')

## Data Preprocessing
### Load data
raw_df = pd.read_csv(datasets / f'RA.G.zeroed.decontam.{threshold}.csv')
meta = pd.read_csv(cwd / "../data/metadata/parsed_patient_metadata.filt.csv")

merged_df = raw_df.merge(meta, on='run_id', how='left')
merged_filt = merged_df.loc[merged_df[col_name].isin([g1, g2]), :]

# Response
y = merged_filt.loc[:, col_name].copy()

# Features
X = merged_filt.loc[:, ~merged_filt.keys().isin(meta.keys())].copy()

# Rename features
X.columns = X.columns.str.replace('[^A-Za-z0-9]+', '_')
print(X.shape)
print(y.shape)

# Binary encode y
y.loc[y == g1] = 1
y.loc[y == g2] = 0
y = y.astype('int')

n_splits = 5

# Read optimal model parameters
param_df = pd.read_csv(f'../results/ML_out/results_out/{g1}.{g2}.{threshold}.results.csv')
param_df = param_df.loc[0, ['colsample_bytree', 'gamma', 'max_depth', 'n_estimators']]
best_params = param_df.to_dict()
best_params['n_jobs'] = 1

# Custom metrics
precision = make_scorer(precision_score, average='binary')
recall = make_scorer(recall_score, average='binary')
f1 = make_scorer(f1_score, average='binary')
auprc = make_scorer(average_precision_score, average=None)

scoring = {'precision': precision,
           'recall': recall,
           'AUROC': 'roc_auc',
           'F1': f1}

# Set up cross-validation
cv = StratifiedKFold(n_splits=n_splits, shuffle=True)

# Generate permutations
np.random.seed(66)
result_list = []

for i in range(1000):
    print(i)
    # Permute labels
    y_temp = np.random.choice(y, size=len(y), replace=False)

    # Train model
    temp_model = XGBClassifier(**best_params)

    cv_results = cross_validate(temp_model, X=X, y=y_temp, cv=cv, scoring=scoring, n_jobs=5)
    temp_df = pd.DataFrame(cv_results).mean()[['test_precision', 'test_recall', 'test_F1', 'test_AUROC']]
    temp_df = pd.DataFrame(temp_df).transpose()
    result_list.append(temp_df)

merged_results = pd.concat(result_list, axis=0)

# Export results
merged_results.to_csv(results / f'permutation_out/{g1}.{g2}.{threshold}.permutations.csv', index=False)


