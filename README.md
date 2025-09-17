# Pipeline Development - Meta-analysis
## Effect size estimation, meta-analysis, and statistical testing for gene expression data.


**Cohen's d** - used for comparing two groups, calculated as the differences between group means divided by the pooled standard deviaton.

**Pearson's r** - measures the strength of a linear relationship between two variables, with an effect size ranging from -1 to + 1.

**Pooled Standard Deviation** - a weighted average of standard deviation from two or more groups, used when assuming the groups have equal variances. It combines individual group varainces to provide a single, more precise estimate of the common variability across the populations.
  - Example:
      Pooled SD = sqrt((SD1^2 + SD2^2) / 2)

**Effect Size** - a statistical measure quantifying the magnitude of a relationship or difference between groups, helping to determine a finding's practical significance beyond statistical significance.
  - Example:
      Cohen's d, where a weight loss intervention group loses 10 kg, and the control group loses 5 kg, resulting in a Cohen's d approximately 0.67, indicating a medium effect size that is a noticeable difference in the weight loss between the groups.
    - Interpreting Cohen's d:
        - **0.2:** Small effect
        - **0.5:** Medium effect
        - **0.8:** Large effect

A bioinformatics toolkit for computing effect sizes, pooling meta-analytic estimates, and performing statistical tests on gene expression data. 

This toolkit calculates summary effect sizes for individual genes, summarizes effect sizes across multiple studies, and determines statistical significance using t-statistics while accounting for variance assumptions. The output includes forest plots and a comprehensive table containing:
Effect sizes (Hedge's g/Cohen's d) for each gene, measuring the magnitude of expression between groups.
Summary effect sizes pooled across multiple datasets.
P-values and adjusted q-values, allowing for statistical significance assessment while controlling for false discovery rates.

These results can be used to identify genes or gene sets with differential expression across conditions, aiding in biomarker discovery and predictive modeling. Currently, this workflow only works for transcriptomic datasets but will be continue to be modified to accept other omic data.


Below is an example of the output:

<img src="BOLA1.jpg" alt="Example 1" width="200">
<img src="CYP4F3.jpg" alt="Example 2" width="200">
<img src="VEGFA.jpg" alt="Example 3" width="200">




Citations: 

K.S Pollard, S. Dudoit, M.J. van der Laan (2005). Multiple Testing Procedures: R multtest Package and
  Applications to Genomics, in Bioinformatics and Computational Biology Solutions Using R and Bioconductor,
  R. Gentleman, V. Carey, W. Huber, R. Irizarry, S. Dudoit (Editors). Springer (Statistics for Biology and
  Health Series), pp. 251-272. [https://doi.org/10.1007/0-387-29362-0_15](https://doi.org/10.1007/0-387-29362-0_15)

T. E. Sweeney, A. Shidham, H. R. Wong, P. Khatri, A comprehensive time-course–
based multicohort analysis of sepsis and sterile inflammation reveals a robust diagnostic
gene set. Sci. Transl. Med. 7, 287ra71 (2015). [DOI: 10.1126/scitranslmed.aaa5993](https://doi.org/10.1126/scitranslmed.aaa5993)

