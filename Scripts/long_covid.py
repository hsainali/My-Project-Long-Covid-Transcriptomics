import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt
import seaborn as sns
import mygene

## Step 1: My first attempt (missed index_col = 0) 
df = pd.read_csv("/Users/hulee/Documents/GitHub/My Project Long Covid Transcriptomics/data/long_covid_annotated.csv")
print (df.head())

## Step 1.1: Fix by adding index_col = 0 
df = pd.read_csv("/Users/hulee/Documents/GitHub/My Project Long Covid Transcriptomics/data/long_covid_annotated.csv", index_col = 0
            )
print (df.head())

## Step 1.2: Verifying across the entire magnitude of the clinical dataset. 
print (f"df.shape: {df.shape[0]:,} genes × {df.shape[1]} samples")

## Step 2: Group the samples.

all_samples = df.columns.tolist()

healthy = [s for s in all_samples if s.startswith("VV")]  # I did a mistake here putting VV for Validated 
pasc = [s for s in all_samples if s.startswith("AA")]     # I swapped this too assuming AA for Acute/Affected

print(f"Healthy samples: {len(healthy)}")
print(f"PASC samples: {len(pasc)}")

## Step 2.1: Correcting the group prefixes. 
all_samples = df.columns.tolist()

healthy = [s for s in all_samples if s.startswith("AA")]  
pasc = [s for s in all_samples if s.startswith("VV")]     

print(f"Healthy samples: {len(healthy)}")
print(f"PASC samples: {len(pasc)}")


## Step 3: Normalization (I forgot the Pseudocount +1) 
cpm = df.div(df.sum(axis=0),axis=1) * 1e6 
log_data = np.log2(cpm)
print (log_data.head())

## Step 3.1: Pseudocount +1
cpm = df.div(df.sum(axis=0),axis=1) * 1e6 
log_data = np.log2(cpm+1)
print (log_data.head())

Calculating the log2 Fold change and The Welch's T-Test. 
log_data ["log2fc"] = log_data[pasc].mean(axis=1) - log_data[healthy].mean(axis=1) 
_, log_data ["pvalue"] = ttest_ind(log_data[pasc], log_data[healthy], axis=1, equal_var=False)
print (f"Analysis Complete: {len(log_data):,} genes processed.")
print (f"Signficiant (Raw): {(log_data["pvalue"] < 0.05).sum():,}") #I named it "Raw" because it still includes False Positives and require further correction. 

# Step 5: FDR correction 
# removing Zero-variance genes (t-test produced Not A Number (NaN))
log_data = log_data.dropna(subset=["pvalue"])
_, log_data["padj"] = fdrcorrection(log_data["pvalue"])
print (f"Genes after FDR correction: {len(log_data):,}")
print (f"Significant after FDR (padj < 0.05): {(log_data["padj"] < 0.05).sum():,}")
                                    
                           
# Step 7 Volcano Plot Initial axis and color coding
# y-axis (-log10 padj) 
log_data["-log10padj"] = - np.log10(log_data["padj"])
# Defining color conditions 
conditions = [(log_data["padj"] < 0.05) & (log_data["log2fc"] >= 0.5), # This means upregulation (RED)
              (log_data["padj"] < 0.05) & (log_data["log2fc"] <= 0.5)] # This measn downregulated (BLUE)

# Assigning colors (red, blue) everything else that does not meet the conditions is grey. 
colors = ["red", "blue"]
log_data["color"] = np.select(conditions, colors, default="grey")

print(f"Red genes (Upregulated): {(log_data["color"] == "red").sum()}")
print(f"Blue genes (Downregulated): {(log_data["color"] == "blue").sum()}")
print(f"Grey genes (Non-significant): {(log_data["color"] == "grey").sum()}")


# Step 7.1: I unintentionally missed the negative sign while setting the conditions which resulted in a condition that accounted for everything. 
 
log_data["-log10padj"] = - np.log10(log_data["padj"])
# Defining color conditions 
conditions = [(log_data["padj"] < 0.05) & (log_data["log2fc"] >= 0.5), # This means upregulation (RED)
              (log_data["padj"] < 0.05) & (log_data["log2fc"] <= -0.5)] # This measn downregulated (BLUE) (added negative sign)

# Assigning colors (red, blue) everything else that does not meet the conditions is grey. 
colors = ["red", "blue"]
log_data["color"] = np.select(conditions, colors, default="grey")

print(f"Red genes (Upregulated): {(log_data["color"] == "red").sum()}")
print(f"Blue genes (Downregulated): {(log_data["color"] == "blue").sum()}")
print(f"Grey genes (Non-significant): {(log_data["color"] == "grey").sum()}")

# Step 7.2 Volcano Plot 
plt.figure(figsize=(10,8))
plt.scatter(x=log_data["log2fc"], y=log_data["-log10padj"], c=log_data["color"], alpha=0.6,
s=10)
        

# Step 7.3 Adding threshold lines for calrity and Labels and title  
plt.axvline(x=0.5, color="black", linestyle="--", linewidth=1)
plt.axvline(x=-0.5, color="black", linestyle="--", linewidth=1)
plt.axhline(y=-np.log10(0.05), color="black", linestyle="--", linewidth=1)

plt.xlabel("Log2 Fold Change", fontsize=12)
plt.ylabel("-Log10 Adjusted P-value", fontsize=12)
plt.title("Long COVID (PASC) vs Healthy Recovered Blood Transcriptome", fontsize=13, fontweight="bold")

plt.savefig("../figures/Long_COVID_Volcano_Plot.png", dpi=300, bbox_inches="tight")
plt.show()

# Step 7.4 Compiling everything for the Volcano Plot. 
log_data["-log10padj"] = - np.log10(log_data["padj"])
# Defining color conditions 
conditions = [(log_data["padj"] < 0.05) & (log_data["log2fc"] >= 0.5), # This means upregulation (RED)
              (log_data["padj"] < 0.05) & (log_data["log2fc"] <= -0.5)] # This measn downregulated (BLUE) (added negative sign)

# Assigning colors (red, blue) everything else that does not meet the conditions is grey. 
colors = ["red", "blue"]
log_data["color"] = np.select(conditions, colors, default="grey")

print(f"Red genes (Upregulated): {(log_data["color"] == "red").sum()}")
print(f"Blue genes (Downregulated): {(log_data["color"] == "blue").sum()}")
print(f"Grey genes (Non-significant): {(log_data["color"] == "grey").sum()}")

plt.figure(figsize=(10,8))
plt.scatter(x=log_data["log2fc"], y=log_data["-log10padj"], c=log_data["color"], alpha=0.6,
s=10)

plt.axvline(x=0.5, color="black", linestyle="--", linewidth=1)
plt.axvline(x=-0.5, color="black", linestyle="--", linewidth=1)
plt.axhline(y=-np.log10(0.05), color="black", linestyle="--", linewidth=1)

plt.xlabel("Log2 Fold Change", fontsize=12)
plt.ylabel("-Log10 Adjusted P-value", fontsize=12)
plt.title("Long COVID (PASC) vs Healthy Recovered Blood Transcriptome", fontsize=13, fontweight="bold")

plt.savefig("../figures/Long_COVID_Volcano_Plot.png", dpi=300, bbox_inches="tight")
plt.show()


# Step 9: PCA 
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Transpose: rows=samples, columns=genes 

expression_cols = healthy + pasc
pca_matrix = log_data[expression_cols].T

scaler = StandardScaler()
pca_scaled = scaler.fit_transform(pca_matrix)

pca = PCA(n_components=2)
pca_coords = pca.fit_transform(pca_scaled)
explained = pca.explained_variance_ratio_ * 100

print(f"PC1 explains: {explained[0]:.1f}% of variance")
print(f"PC2 explains: {explained[1]:.1f}% of variance")


# Step 9.1: Plottig PCA 

labels = ["Healthy"] * len(healthy) + ["PASC"] * len(pasc)
plt.figure(figsize=(10, 8))
for group, color in [("Healthy", '#3498DB'), ('PASC', '#E74C3C')]:
    mask = [l == group for l in labels]
    plt.scatter(pca_coords[mask, 0],pca_coords[mask, 1], c=color, label=group, alpha=0.7, s=60)

plt.xlabel(f"PC1 ({explained[0]:.1f}% variance)", fontsize=13)
plt.ylabel(f"PC2 ({explained[1]:.1f}% variance)", fontsize=13)
plt.title("PCA of Blood Transcriptome: Healthy vs Long COVID (PASC)", fontsize=14, fontweight="bold")
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.savefig("../figures/Long_COVID_PCA.png", dpi=300, bbox_inches="tight")
plt.show()

