# Disease model
This repository provides four machine learning models used for binary classification.

## Training models
Common imports (used by all models)

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
```



### 1.Random forest classifier

```python
from sklearn.ensemble import RandomForestClassifier

# ==============================
# Load data
# ==============================
train_df = pd.read_csv('./Pos_Neg.csv')
unknown_df = pd.read_csv('./all_ATS.csv')

X_train = train_df.drop(columns=['Label'])
y_train = train_df['Label']
X_unknown = unknown_df.copy()

# ==============================
# Build model pipeline
# ==============================
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', RandomForestClassifier(
        n_estimators=100,
        random_state=42,
        class_weight='balanced'
    ))
])

# ==============================
# 10-fold cross-validation (probability output)
# ==============================
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

y_pred_proba_cv = cross_val_predict(
    pipeline,
    X_train,
    y_train,
    cv=cv,
    method='predict_proba'
)[:, 1]

# ==============================
# ROC and PR metrics
# ==============================
fpr, tpr, _ = roc_curve(y_train, y_pred_proba_cv)
roc_auc = auc(fpr, tpr)

precision, recall, _ = precision_recall_curve(y_train, y_pred_proba_cv)
ap_score = average_precision_score(y_train, y_pred_proba_cv)

# ==============================
# Visualization (Random Forest only)
# ==============================
plt.figure(figsize=(8, 4))

plt.subplot(1, 2, 1)
plt.plot(fpr, tpr, lw=2, label=f'AUC = {roc_auc:.2f}')
plt.plot([0, 1], [0, 1], linestyle='--')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC Curve')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(recall, precision, lw=2, label=f'AP = {ap_score:.2f}')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('Precision-Recall Curve')
plt.legend()

plt.tight_layout()
plt.show()

# ==============================
# Save CV predictions
# ==============================
train_df['predictions'] = y_pred_proba_cv
train_df.to_csv('04.10Fold_predictions_rf.csv', index=False)

# ==============================
# Train final model and feature importance
# ==============================
pipeline.fit(X_train, y_train)

importances = pipeline.named_steps['classifier'].feature_importances_
indices = np.argsort(importances)[::-1]

sorted_features = X_train.columns[indices]
sorted_importances = importances[indices]

pd.DataFrame(
    sorted_importances,
    index=sorted_features,
    columns=['Importance']
).to_csv('RF.importance.csv')

plt.figure()
plt.barh(sorted_features, sorted_importances)
plt.xlabel('Importance')
plt.title('Feature Importance (Random Forest)')
plt.gca().invert_yaxis()
plt.savefig('04.rf_importance.pdf', dpi=300)
plt.show()

# ==============================
# Predict unknown samples
# ==============================
X_unknown['Probability'] = pipeline.predict_proba(X_unknown)[:, 1]
X_unknown.to_csv('04.unknown_predictions_rf.csv', index=False)
```



### 2.Logistic regression

```python
from sklearn.linear_model import LogisticRegression

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', LogisticRegression(
        random_state=42,
        max_iter=10000,
        class_weight='balanced'
    ))
])

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

aucs, aps = [], []

for train_idx, test_idx in cv.split(X_train, y_train):
    pipeline.fit(X_train.iloc[train_idx], y_train.iloc[train_idx])
    y_prob = pipeline.predict_proba(X_train.iloc[test_idx])[:, 1]

    fpr, tpr, _ = roc_curve(y_train.iloc[test_idx], y_prob)
    aucs.append(auc(fpr, tpr))
    aps.append(average_precision_score(y_train.iloc[test_idx], y_prob))

print(f"Average AUC: {np.mean(aucs):.3f} ± {np.std(aucs):.3f}")
print(f"Average AP:  {np.mean(aps):.3f} ± {np.std(aps):.3f}")

pipeline.fit(X_train, y_train)

# Coefficient-based feature importance
coef = pipeline.named_steps['classifier'].coef_[0]
indices = np.argsort(np.abs(coef))[::-1]

pd.DataFrame(
    coef[indices],
    index=X_train.columns[indices],
    columns=['Coefficient']
).to_csv('04.lr.importance.csv')

X_unknown['Probability'] = pipeline.predict_proba(X_unknown)[:, 1]
X_unknown.to_csv('04.unknown_predictions_lr.csv', index=False)
```



### 3.Support vector machine (SVM)s

```python
from sklearn.svm import SVC
from sklearn.inspection import permutation_importance

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', SVC(probability=True, random_state=42))
])

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

aucs, aps = [], []

for train_idx, test_idx in cv.split(X_train, y_train):
    pipeline.fit(X_train.iloc[train_idx], y_train.iloc[train_idx])
    y_prob = pipeline.predict_proba(X_train.iloc[test_idx])[:, 1]

    fpr, tpr, _ = roc_curve(y_train.iloc[test_idx], y_prob)
    aucs.append(auc(fpr, tpr))
    aps.append(average_precision_score(y_train.iloc[test_idx], y_prob))

print(f"Average AUC: {np.mean(aucs):.3f} ± {np.std(aucs):.3f}")
print(f"Average AP:  {np.mean(aps):.3f} ± {np.std(aps):.3f}")

pipeline.fit(X_train, y_train)

# Permutation importance
result = permutation_importance(
    pipeline,
    X_train,
    y_train,
    n_repeats=10,
    random_state=42
)

np.savetxt('04.svm.importance.mean.txt', result.importances_mean)
np.savetxt('04.svm.importance.std.txt', result.importances_std)

X_unknown['Probability'] = pipeline.predict_proba(X_unknown)[:, 1]
X_unknown.to_csv('04.unknown_predictions_svm.csv', index=False)
```



### 4.Neural network (MLP)

```python
from sklearn.neural_network import MLPClassifier

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', MLPClassifier(
        hidden_layer_sizes=(100,),
        max_iter=1000,
        random_state=42
    ))
])

cv = StratifiedKFold(n_splits=15, shuffle=True, random_state=42)

aucs, aps = [], []

for train_idx, test_idx in cv.split(X_train, y_train):
    pipeline.fit(X_train.iloc[train_idx], y_train.iloc[train_idx])
    y_prob = pipeline.predict_proba(X_train.iloc[test_idx])[:, 1]

    fpr, tpr, _ = roc_curve(y_train.iloc[test_idx], y_prob)
    aucs.append(auc(fpr, tpr))
    aps.append(average_precision_score(y_train.iloc[test_idx], y_prob))

print(f"Average AUC: {np.mean(aucs):.3f} ± {np.std(aucs):.3f}")
print(f"Average AP:  {np.mean(aps):.3f} ± {np.std(aps):.3f}")

pipeline.fit(X_train, y_train)

X_unknown['Probability'] = pipeline.predict_proba(X_unknown)[:, 1]
X_unknown.to_csv('04.unknown_predictions_nn.csv', index=False)
```



## Using pretrained models for prediction

This document explains how to use a **trained sklearn `.pkl` model** (saved as a full `Pipeline`) to predict probabilities for new, unseen data.

The workflow is designed to be **reproducible**, **safe**, and **GitHub-ready**.



**Saving a trained model**

**Best practice:** always save the entire `Pipeline` (including preprocessing and the classifier).

This ensures that the same data transformations used during training are automatically applied during prediction.

```python
import joblib

# Train the pipeline on the full training dataset
pipeline.fit(X_train, y_train)

# Save the trained pipeline
joblib.dump(pipeline, 'rf_pipeline.pkl')
```

**Why save the full pipeline?**
- Prevents data leakage
- Avoids manual feature scaling during prediction
- Guarantees consistency between training and inference



**Predicting new data**

```python
import pandas as pd
import joblib

# ==============================
# Load trained pipeline
# ==============================
pipeline = joblib.load('rf_pipeline.pkl')

# ==============================
# Load new (unlabeled) data
# ==============================
# The feature columns must match the training data
new_data = pd.read_csv('./new_samples.csv')

# ==============================
# Predict probabilities
# ==============================
# Returns probability of the positive class (Label = 1)
probabilities = pipeline.predict_proba(new_data)[:, 1]

# ==============================
# Save prediction results
# ==============================
new_data['Probability'] = probabilities
new_data.to_csv('rf_predictions_new_data.csv', index=False)
```

This script:
- Does **not** retrain the model
- Automatically applies preprocessing
- Can be reused for any saved pipeline (RF / LR / SVM / NN)



### Notes and important considerations

**Feature consistency**
- Feature **names**, **order**, and **number** must be identical to the training data
- Missing or additional columns will cause errors or invalid predictions

**Recommended safety check:**
```python
assert list(new_data.columns) == list(X_train.columns)
```


**Do not re-scale new data manually**
- The scaler is already stored inside the pipeline
- Do **not** apply `StandardScaler.fit_transform()` again

Incorrect:
```python
scaler = StandardScaler()
new_data = scaler.fit_transform(new_data)
```

Correct:
```python
pipeline.predict_proba(new_data)
```



`predict()` vs `predict_proba()`

- `predict()` returns class labels (0 or 1)
- `predict_proba()` returns probabilities

For ranking, scoring, or biological risk analysis, **`predict_proba()` is strongly recommended**.



**Model compatibility**

The same prediction script works for all saved models:

```python
pipeline = joblib.load('lr_pipeline.pkl')
pipeline = joblib.load('svm_pipeline.pkl')
pipeline = joblib.load('nn_pipeline.pkl')
```

This is a key advantage of using sklearn Pipelines.