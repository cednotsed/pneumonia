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

pos = len(y[y == 1])
neg = len(y[y == 0])
split_sizes = pd.DataFrame({g1: [pos - int(pos / n_splits), int(pos / n_splits)],
                            g2: [neg - int(neg / n_splits), int(neg / n_splits)]}, index=['Train fold', 'Test fold'])

print(split_sizes)

# Get negative to positive ratio
ratio = sum(y == 0) / sum(y == 1)


# Nested CV function
def optimise_evaluate(X, y):
    np.random.seed(66)
    ratio = sum(y == 0) / sum(y == 1)

    # Hyperparemeter Optimisation using grid search (F1)
    model = XGBClassifier()
    n_estimators = range(100, 1000, 100)
    max_depth = range(1, 5, 1)
    gamma = np.linspace(0.1, 3, 10)
    colsample_bytree = np.linspace(0.1, 1, 10)

    param_grid = dict(max_depth=max_depth,
                      n_estimators=n_estimators,
                      colsample_bytree=colsample_bytree,
                      gamma=gamma,
                      n_jobs=[1])

    inner_cv = StratifiedKFold(n_splits=n_splits, shuffle=True)
    outer_cv = StratifiedKFold(n_splits=n_splits, shuffle=True)

    # Inner CV
    model = GridSearchCV(model,
                         param_grid,
                         scoring='f1',
                         n_jobs=6,
                         cv=inner_cv,
                         verbose=1)

    model.fit(X, y)
    best_params = model.best_params_
    print(best_params)

    # Custom metrics
    precision = make_scorer(precision_score, average='binary')
    recall = make_scorer(recall_score, average='binary')
    f1 = make_scorer(f1_score, average='binary')
    auprc = make_scorer(average_precision_score, average=None)

    scoring = {'precision': precision,
               'recall': recall,
               'AUROC': 'roc_auc',
               'F1': f1}

    # Outer CV
    outer_results = cross_validate(model, X=X, y=y, cv=outer_cv, scoring=scoring, n_jobs=4)
    outer_results = pd.DataFrame(outer_results).mean()[['test_precision', 'test_recall', 'test_F1', 'test_AUROC']]

    return outer_results, best_params


# Run nested CV
raw_results, raw_params = optimise_evaluate(X, y)

# Parse results and model parameters
results_df = pd.DataFrame(raw_results).transpose()
param_df = pd.DataFrame(raw_params, index=[0])
param_df['g1'] = g1
param_df['g2'] = g2
param_df['threshold'] = threshold

final_results = pd.concat([param_df, results_df], axis=1)

# Model interpretation
final_model = XGBClassifier(**raw_params)
final_model.fit(X=X, y=y)

explainer = shap.TreeExplainer(final_model,
                               feature_perturbation='interventional',
                               model_output='probability',
                               data=X)
shap_values = explainer.shap_values(X)

# Parse SHAP values
shap_df = pd.DataFrame(shap_values, columns=X.columns).add_suffix('.shap')
merged = pd.concat([X, shap_df], axis=1)

# Export results
merged.to_csv(results / f'shap_out/{g1}.{g2}.{threshold}.shap.csv', index=False)
final_results.to_csv(results / f'results_out/{g1}.{g2}.{threshold}.results.csv', index=False)

shap.summary_plot(shap_values, X,
                  show=False,
                  plot_size=(9, 9),
                  color_bar_label='Relative abundance',
                  max_display=10)
fig, ax = plt.gcf(), plt.gca()
ax.set_xlabel('SHAP Value')
plt.savefig(results / f'beeswarm_out/{g1}.{g2}.{threshold}.pdf', dpi=600, format='pdf')

