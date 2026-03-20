# Acknowledgment and Reference : 
As I am a beginner in coding and learning how ot build functions and understand statistical tools, this project was completed using multiple online resources that helps me understand the code and implement it in my project.
(This project is done to develop my skills and interest in computational biology and may not follow the professional data analysis accuracy and production yet). 

I reference: 

**DESeq2 (Original):** [Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.](https://doi.org/10.1186/s13059-014-0550-8)
 
**PyDESeq2:** [Muzellec, B., Teleńczuk, M., Cabeli, V., & Andreux, M. (2023). PyDESeq2: a python package for bulk RNA-seq differential expression analysis. Bioinformatics, 39(9), btad547.](https://doi.org/10.1093/bioinformatics/btad547)

## **For data formatting:**
**Mike Saint-Antoine:** [Intro to Bioinformatics 7: Downloading Data from the Gene Expression Omnibus (GEO).](https://www.youtube.com/watch?v=4am5XF_597A) 
[And Intro to Bioinformatics 4: Gene Expression Data Format.](https://www.youtube.com/watch?v=EcLF-L-HDf0)

## **Normalization and Differences**
**Harvard Chan Bioinformatics Core** [Sample-level QC.](https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/sample_level_QC.html)

**Lei Guo:** [Normalize RNA-seq Gene Expression: TPM, FPKM, CPM, Z-scores.](https://www.youtube.com/watch?v=zqB_qdRagH4)

**StatQuest with Josh Starmer:** [StatQuest: PCA in Python.](https://www.youtube.com/watch?v=Lsue2gEM9D0) 

## **Differential Expression Analysis**
**Sanbomics:** [Differential expression in Python with pyDESeq2.](https://www.youtube.com/watch?v=wIvxFEMQVwg)

**Bioinformatics Coach:** [DESeq2 Tutorial | How I analyze RNA Seq Gene Expression data using DESeq2.](https://www.youtube.com/watch?v=kOlMcZujHHA)

**OMGenomics:** [RNA-seq tutorial with DESeq2: Differential gene expression project.](https://www.youtube.com/watch?v=NGbZmlGLG5w)

**RunCell:** [Differential Expression Analysis.](https://www.runcell.dev/blog/differential-expression-analysis)

## **Statistical Testing and methods:**
**Vincent Stevenson:** [Calculate FDR Adjusted P-values in Python - Benjamini Hochberg Procedure.](https://www.youtube.com/watch?v=MHjblrdzRuU)

**GeeksforGeeks:** [Welch's t-test in Python.](http://geeksforgeeks.org/machine-learning/welchs-t-test-in-python/)

**DSC Data Science Concepts:** [Omitted Variable Bias. Wald Test in Python (Jupyter).](https://www.youtube.com/watch?v=nv27tbea-qI)

**Andrew P. Wheeler:** [Wald Tests via Statsmodels (Python).](https://andrewpwheeler.com/2021/06/18/wald-tests-via-statsmodels-python/)

## **Volcano plot**
**BioInfo Tips:** [How to Make a Volcano Plot for RNA-Seq Data Analysis | Python.](https://www.youtube.com/watch?v=qHeIhEPKU98)

# Longitudinal Transcriptomic Mapping of Persistent Immune Dysregulation in Post-Acute Sequelae of SARS-CoV-2, PASC. 
# **Currently under review for using inaccurate statistical tool** - Husain.
## Dataset used 

The **dataset** used in this project is adapted from the National Centre for Biotechnology Information, Gene Expression Omnibus (NCBI, GEO) accession number GSE267625 (Peripheral blood transcriptomes RNA sequencing for 111 patients (Healthy Recovered vs Long COVID) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267625 . 

**Cohorts:** Healthy Recovered (AA) vs Long COVID (PASC) (VV).

**Timepoints:** 3-6 months and 12-15 months post infection. (longitudinal, same individuals)
## Background 
Post-Acute Sequelae of SARS-CoV-2 (PASC), colloquially known as Long COVID, remains one of the most poorly understood and debilitating consequences of the COVID-19 pandemic, affecting approximately [10–20% of all COVID-19 survivors](https://www.medrxiv.org/content/10.1101/2024.01.14.24301293v1). PASC presents as a highly heterogeneous, multisystemic syndrome encompassing severe chronic fatigue, cognitive dysfunction, cardiovascular anomalies, and autonomic nervous system dysregulation. The precise transcriptomic mechanisms driving susceptibility and chronicity remain largely undefined,though prevailing hypotheses centre on [immune exhaustion](https://www.nature.com/articles/s41590-023-01724-6), [microvascular clotting cascades](https://www.cell.com/cell/fulltext/S0092-8674(24)00886-9), and [long-term viral persistence within distinct physiological reservoirs](https://www.science.org/doi/10.1126/scitranslmed.adk3295) (tissue sites where the virus evades immune clearance and continues replicating at low levels).

## Pipeline Overview 

Step 1, Data pre-processing:  I implemented python code to ingest the massive raw dataset (contains RNA levels of 34,000+ genes across 111 different human patients). Then filtered out the lowly-expressed genes and extracted for analysis the sample patients groups (Healthy Recovered vs Long COVID PASC) from the 111 patients' samples. 

Step 2: Gene Annotation: Translated the database IDs into standard gene symbols using "MyGene" to allow for biological interpretation.

Step 3: Reducing Dimentionality: Executred PCA on the transformed count data to visualise cohort variance. 

Step 4: Differencial Gene Expression Analysis: modeled the raw data using PyDESeq2 to identify signfiicantly dysregulated genes between the two groups. (Wald t-test) 

Step 5: Volcano Plot: constructed to visualize the data distribution and highlight specific biological markers markers 

## Results and Findings

This project outcome identified 31 significantly upregulated genes along with 290 significantly downregulated genes. This can in fact suggest against the idea that PASC is driven by the hyperactive inflammatory cascade usually seen in acute COVID-19. In terms of the upregulated genes, there are genes like VWDE, a gene recently established as a network driver for Long COVID susceptibility (Pinero et al., 2025). The overexpression of VWDE points towards ongoing vascular dysfunction, which is a central hallmark of the PASC phenotype (Pinero et al., 2025). Moreover, the upregulation of TRPC1, a mechanosensitive calcium channel, further compounds this vascular deterioration. Research shows that the interaction of the SARS-CoV-2 spike protein with host receptors promotes the aberrant activation of calcium channels, including TRPC1, triggering an excessive intracellular influx that ultimately leads to pulmonary microvascular endothelial apoptosis (Aljadah et al., 2024). Furthermore, the dataset reveals a marked upregulation of ALPK2 (atypical kinase). ALPK2 has been identified in recent studies as a primary pathophysiological driver and therapeutic target for heart failure with preserved ejection fraction (HFpEF) and age-related cardiac diastolic dysfunction (Yoshida et al., 2024). This demonstrates the presence of a persistent cardiovascular vulnerability in PASC patients. The significant upregulation of pEG10 (retrotransposon-derived gene) aligns with the current literature showing that SARS-CoV-2 infection triggers specific epigenetic alterations. Specifically, the virus induces the hypomethylation and subsequent overexpression of PEG10, which ultimately drives the abnormal cellular proliferation and facilitates viral replication (Li et al., 2021). In parallel, the dataset shows an increase in NCR3LG1, a ligand that activates natural killer (NK) cells. This can lead to indiscriminate immune-mediated tissue destruction that drives continuous damage to host organs long after the initial viral clearance (Wang et al., 2022).

Conversely, among the 290 downregulated genes, mitochondrial regulatory pseudogenes had a significant, statistically significant suppression, such as MTND1P23, and core structural genes of the electron transport chain, such as MT-ND1. This targeted suppression of oxidative phosphorylation pathways could potentially force highly metabolic tissues to rely on inefficient anaerobic glycolysis (Ryan et al., 2023). Moreover, this virus-induced metabolic reprogramming directly results in elevated cellular stress and the chronic physical fatigue characteristic of Long COVID (Ryan et al., 2023). The data also shows a massive suppression of IFI27, a primary interferon-stimulated gene. Although IFI27 is highly active and serves as a predictor of disease severity during the acute phase of COVID-19, its downregulation in the chronic phase indicates a potentially depleted, unresponsive innate antiviral network (Alcalde-Herraiz et al., 2024).


## Statistical tools 

1. Negative Binomial GLM & Wald Test (PyDESeq2) : This is used becasue RNA-seq count data is highly overdispersed and heterogeneous (especially in the PASC cohort). Standard linear tests (like a t-test) would be insufficient. 

2. FDR correction: Relying on standard p-value (0.05) solely would increase the false positive that leads to pathways that don't actually exist in identifying diagnostic biomarkers. Benjamini-Hochberg False Discovery Rate would force the correction of raw p-values based on their rank (Adjusted p-values).




