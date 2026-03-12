import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt

# 1. Load the Clean, Annotated Data
path = '/Users/hulee/ngly1_project/long_covid_project/annotated_long_covid.csv'
df = pd.read_csv(path, index_col=0)

# 2. Group the Samples
all_samples = df.columns.tolist()
recovered_cols = [s for s in all_samples if s.startswith('AA')] # Healthy
pasc_cols = [s for s in all_samples if s.startswith('VV')]      # Long COVID

# 3. Normalization (Log2 CPM)
cpm = df.div(df.sum(axis=0), axis=1) * 1e6
log2_cpm = np.log2(cpm + 1)

# 4. Math: Log2 Fold Change & Welch's T-Test
log2_cpm['mean_healthy'] = log2_cpm[recovered_cols].mean(axis=1)
log2_cpm['mean_pasc'] = log2_cpm[pasc_cols].mean(axis=1)
log2_cpm['log2FC'] = log2_cpm['mean_pasc'] - log2_cpm['mean_healthy']

t_stats, p_values = ttest_ind(log2_cpm[pasc_cols], log2_cpm[recovered_cols], axis=1, equal_var=False)
log2_cpm['pvalue'] = p_values
log2_cpm = log2_cpm.dropna(subset=['pvalue'])

# 5. Error Control: FDR Correction
_, log2_cpm['padj'] = fdrcorrection(log2_cpm['pvalue'])
log2_cpm['-log10(padj)'] = -np.log10(log2_cpm['padj'])

# 6. Volcano Plot Setup
p_val_thresh = 0.05
logfc_thresh = 0.5 # Using a tighter threshold for blood transcriptomics

conditions = [
    (log2_cpm['padj'] <= p_val_thresh) & (log2_cpm['log2FC'] >= logfc_thresh),
    (log2_cpm['padj'] <= p_val_thresh) & (log2_cpm['log2FC'] <= -logfc_thresh)
]
choices = ['#E74C3C', '#3498DB'] # Red (Upregulated in PASC), Blue (Downregulated)
log2_cpm['color'] = np.select(conditions, choices, default='#95A5A6')

plt.figure(figsize=(10, 8), dpi=300)
plt.scatter(log2_cpm['log2FC'], log2_cpm['-log10(padj)'], c=log2_cpm['color'], alpha=0.5, s=15)

plt.axvline(x=logfc_thresh, color='black', linestyle='--', linewidth=1)
plt.axvline(x=-logfc_thresh, color='black', linestyle='--', linewidth=1)
plt.axhline(y=-np.log10(p_val_thresh), color='black', linestyle='--', linewidth=1)

# 7. Highlight Immune Exhaustion Targets
target_markers = ['CD8A', 'CD4', 'PDCD1', 'LAG3', 'FYN', 'STAT3']
for gene in target_markers:
    if gene in log2_cpm.index:
        x_val = log2_cpm.loc[gene, 'log2FC']
        y_val = log2_cpm.loc[gene, '-log10(padj)']
        plt.annotate(gene, (x_val, y_val), fontweight='bold', color='black', 
                     xytext=(5,5), textcoords="offset points")

plt.title('Long COVID (PASC) vs Healthy Recovered Blood Transcriptome')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-Log10 Adjusted P-value')

# 8. Export Data
plt.savefig('/Users/hulee/ngly1_project/long_covid_project/Long_COVID_Volcano.png')
top_immune_genes = log2_cpm[log2_cpm['padj'] <= 0.05].sort_values(by='padj').head(100)
top_immune_genes[['log2FC', 'padj']].to_csv('/Users/hulee/ngly1_project/long_covid_project/Top_100_Long_COVID_Genes.csv')

print("Success! Volcano plot and top 100 immune genes exported.")
plt.show()