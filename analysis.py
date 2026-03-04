"""
Bacillus subtilis Essential Genes Analysis
-------------------------------------------
Dataset: Small dataset of B. subtilis genes sourced from SubtiWiki
         (https://subtiwiki.uni-goettingen.de/)

This script demonstrates:
  1. Linear Regression  – predicting gene expression level from gene properties
  2. Logistic Regression – predicting whether a gene is essential (1) or not (0)
"""

from pathlib import Path

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    mean_squared_error,
    r2_score,
    accuracy_score,
    classification_report,
    confusion_matrix,
)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── 1. Load data ──────────────────────────────────────────────────────────────

df = pd.read_csv(Path(__file__).parent / "data" / "bacillus_subtilis_genes.csv")
print("Dataset shape:", df.shape)
print(df.head(), "\n")

FEATURES = ["gene_length_bp", "gc_content", "num_interactions"]

# ── 2. Linear Regression: predict expression_level ───────────────────────────

print("=" * 60)
print("LINEAR REGRESSION – predicting expression_level")
print("=" * 60)

X_lin = df[FEATURES]
y_lin = df["expression_level"]

X_train_l, X_test_l, y_train_l, y_test_l = train_test_split(
    X_lin, y_lin, test_size=0.2, random_state=42
)

scaler_lin = StandardScaler()
X_train_l_scaled = scaler_lin.fit_transform(X_train_l)
X_test_l_scaled = scaler_lin.transform(X_test_l)

lin_model = LinearRegression()
lin_model.fit(X_train_l_scaled, y_train_l)

y_pred_lin = lin_model.predict(X_test_l_scaled)
mse = mean_squared_error(y_test_l, y_pred_lin)
r2 = r2_score(y_test_l, y_pred_lin)

print(f"MSE : {mse:.4f}")
print(f"R²  : {r2:.4f}")
print("Coefficients:", dict(zip(FEATURES, lin_model.coef_)))
print()

# Plot actual vs predicted expression levels
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(y_test_l, y_pred_lin, color="steelblue", edgecolors="k", alpha=0.7)
ax.plot(
    [y_test_l.min(), y_test_l.max()],
    [y_test_l.min(), y_test_l.max()],
    "r--",
    lw=1.5,
    label="Perfect fit",
)
ax.set_xlabel("Actual Expression Level")
ax.set_ylabel("Predicted Expression Level")
ax.set_title("Linear Regression – Actual vs Predicted")
ax.legend()
plt.tight_layout()
plt.savefig("linear_regression_results.png", dpi=150)
plt.close()
print("Saved: linear_regression_results.png\n")

# ── 3. Logistic Regression: predict is_essential ─────────────────────────────

print("=" * 60)
print("LOGISTIC REGRESSION – predicting gene essentiality")
print("=" * 60)

X_log = df[FEATURES]
y_log = df["is_essential"]

X_train_g, X_test_g, y_train_g, y_test_g = train_test_split(
    X_log, y_log, test_size=0.2, random_state=42
)

scaler_log = StandardScaler()
X_train_g_scaled = scaler_log.fit_transform(X_train_g)
X_test_g_scaled = scaler_log.transform(X_test_g)

log_model = LogisticRegression(random_state=42, max_iter=200)
log_model.fit(X_train_g_scaled, y_train_g)

y_pred_log = log_model.predict(X_test_g_scaled)
accuracy = accuracy_score(y_test_g, y_pred_log)

print(f"Accuracy : {accuracy:.4f}")
print("\nClassification Report:")
print(classification_report(y_test_g, y_pred_log, target_names=["Non-essential", "Essential"]))

cm = confusion_matrix(y_test_g, y_pred_log)
print("Confusion Matrix:")
print(cm)
print()

# Plot confusion matrix
fig, ax = plt.subplots(figsize=(5, 4))
im = ax.imshow(cm, interpolation="nearest", cmap=plt.cm.Blues)
plt.colorbar(im, ax=ax)
tick_labels = ["Non-essential", "Essential"]
ax.set_xticks([0, 1])
ax.set_yticks([0, 1])
ax.set_xticklabels(tick_labels)
ax.set_yticklabels(tick_labels)
ax.set_xlabel("Predicted Label")
ax.set_ylabel("True Label")
ax.set_title("Logistic Regression – Confusion Matrix")
for i in range(2):
    for j in range(2):
        ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                color="white" if cm[i, j] > cm.max() / 2 else "black")
plt.tight_layout()
plt.savefig("logistic_regression_confusion_matrix.png", dpi=150)
plt.close()
print("Saved: logistic_regression_confusion_matrix.png\n")

print("Analysis complete.")
