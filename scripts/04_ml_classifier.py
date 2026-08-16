# scripts/04_ml_classifier.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix

print("Loading processed proteomics dataset...")
df = pd.read_csv("data/processed/normalized_proteomics.csv", index_col=0)

# Transpose matrix: samples as rows, proteins as features
X = df.T
y = [idx.split('_')[0] for idx in X.index]

# Stratified Train-Test Split (80/20)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# Training Random Forest Classifier
print("Training Random Forest Classifier on proteomic profiles...")
rf = RandomForestClassifier(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)

# Model Evaluation
y_pred = rf.predict(X_test)
print("\n--- Model Performance Report ---")
print(classification_report(y_test, y_pred))

# Saving Confusion Matrix Plot
cm = confusion_matrix(y_test, y_pred, labels=rf.classes_)
plt.figure(figsize=(6, 5))
sns.heatmap(cm, annot=True, fmt='d', xticklabels=rf.classes_, yticklabels=rf.classes_, cmap='Blues')
plt.title("Pan-Cancer Classification Confusion Matrix")
plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.tight_layout()
plt.savefig("results/confusion_matrix.png")
print("Saved confusion matrix plot to results/confusion_matrix.png")

# Extracting Top 20 Feature Importances (Biomarkers)
importances = pd.Series(rf.feature_importances_, index=X.columns).sort_values(ascending=False).head(20)
plt.figure(figsize=(8, 6))
sns.barplot(x=importances.values, y=importances.index, hue=importances.index, legend=False, palette="viridis")
plt.title("Top 20 Proteomic Biomarker Feature Importances")
plt.xlabel("Feature Importance Score")
plt.ylabel("Protein")
plt.tight_layout()
plt.savefig("results/top20_biomarkers.png")
print("Saved biomarker importance plot to results/top20_biomarkers.png")

# Exporting Top Biomarkers CSV
importances.to_csv("results/top20_biomarker_features.csv", header=["Importance_Score"])
print("Saved biomarker CSV to results/top20_biomarker_features.csv")
