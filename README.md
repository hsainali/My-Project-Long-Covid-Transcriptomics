# Longitudinal Transcriptomic Mapping of Persistent Immune Dysregulation in Post-Acute Sequelae of SARS-CoV-2, PASC. 
# **Currently under review for using inaccurate statistical tool** - Husain.
## Dataset used 

The **dataset** used in this project is adapted from the National Centre for Biotechnology Information, Gene Expression Omnibus (NCBI, GEO) accession number GSE267625 (Peripheral blood transcriptomes RNA sequencing for 111 patients (Healthy Recovered vs Long COVID) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267625 . 

**Cohorts:** Healthy Recovered (AA) vs Long COVID (PASC) (VV).

**Timepoints:** 3-6 months and 12-15 months post infection. (longitudinal, same individuals)
## Background 
Post-Acute Sequelae of SARS-CoV-2 (PASC), colloquially known as Long COVID, remains one of the most poorly understood and debilitating consequences of the COVID-19 pandemic, affecting approximately [10–20% of all COVID-19 survivors](https://www.medrxiv.org/content/10.1101/2024.01.14.24301293v1). PASC presents as a highly heterogeneous, multisystemic syndrome encompassing severe chronic fatigue, cognitive dysfunction, cardiovascular anomalies, and autonomic nervous system dysregulation. The precise transcriptomic mechanisms driving susceptibility and chronicity remain largely undefined,though prevailing hypotheses centre on [immune exhaustion](https://www.nature.com/articles/s41590-023-01724-6), [microvascular clotting cascades](https://www.cell.com/cell/fulltext/S0092-8674(24)00886-9), and [long-term viral persistence within distinct physiological reservoirs](https://www.science.org/doi/10.1126/scitranslmed.adk3295) (tissue sites where the virus evades immune clearance and continues replicating at low levels).
## Pipeline Overview 

Step 1: I wrote python code to ingest the massive raw dataset (contains RNA levels of 34,000+ genes across 111 different human patients) 

Step 2: I filtered out the noise (filtering low-expression biological noise) and separated the 111 patients into "Healthy Recovered" and "Long COVID" groups. 

Step 3: Ran an unsupervised machine learning technique (Principal Component Analysis, PCA) to aid visualization and be able to reduce the number of variables in the dataset while retaining as much information as possible. This way I can prove that long COVID patients have a distinct, overarching immune system signature that can be further investigated. 

Step 4: I translated the unreadable data base using API (MyGene) into recognizable gene symbols (e.g FYN and PDCD1). 

Step 5: Executed advanced statistical testing (beyond p value) such as (Welch's t-test and FDR correction) to be able to isolate the top 100 most significantly dysregulated genes. 

Step 6: Finally I was able to generate a Volcano Plot to map that data and visually highlight specific markers of T-cell exhaustion and inflammation. 

## Results and Findings

I interpreted the biological meaning from the raw math (output) to prove that Long COVID is driven by trapped, exhausted T-cells and chronic hyper activation. I found that there was significant up regulation of SCIMP and the P2RX7/PANX1 axis, which locks macrophages into a state of chronic inflammasome activation and oxidative stress as well. 
In terms of Cytoskeletal paralysis: I identified a failure in T-cell migration (through GRK3 up regulation that blinds T-cells to viral locations, while the co-upregulation of actin accelerators (TRIO) and brakes (ARHGAP31) eventually leaves the cell structurally paralyzed and unable to form an immunological synapse. Additionally, BCL2 massive up regulation makes these broken T-cells immortal and resistant to apoptosis. 
Overall, the PCA results showed 40.6% variance representing a strong signal that Long COVID and Healthy Recovered patients have a distinct biological difference (e.g chronic hyper activation of macrophages, cytoskeletal paralysis preventing T-cells from migrating, and BCL2 up regulation making broken T-cell immortal.)

## Statistical tools 

1. Log2 CPM Normalization: This is important to overcome the bias of raw sequencing counts in library size data. Total number of reads the machine generates for a sample A that has 2 million reads will be upregulated compared to sample B that has 1 million reads simply because the machine ran longer. Counts Per Million (CPM) will help converting counts into relative proportions of the total library data. 
I then applied log2 transformation because gene expression data is not straightforward as some genes have few reads and other genes can have millions. This can help normalize the distribution of the data and ensure a 2-fold change in a lowly expressed gene carries the same mathematical weight as a 2-fold change in a highly expressed gene. This makes the data more suitable for linear statistical testing. 

2. Welch's T-Test: I used this test because standard t-test by default assume Homoscedasticity (Healthy Recovered & Long Covid groups have equal variance). This was not the case. While Healthy group was biologically stable, the PASC group was extremely heterogenous and results in an unequal variance. Welch's t-test would address this issue and for example it did not potentially hide genes dysregulation like GRK3 or PDCD1. 

3. FDR correction: Relying on standard p-value (0.05) solely would increase the false positive that leads to pathways that don't actually exist in identifying diagnostic biomarkers. Benjamini-Hochberg False Discovery Rate would force the correction of raw p-values based on their rank (Adjusted p-values). Thus, for 16,842 high-confidence genes where p-value by random chance would make 842 genes appear significant, FDR would narrow down the list to high-confidence signatures of immune exhaustion. For example, LINC01876 can be said to be true molecular footprint of long COVID rather than statistical artifact. 




