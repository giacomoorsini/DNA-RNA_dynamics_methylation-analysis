# DNA-RNA_dynamics_methylation-analysis
This repository contains the steps to analyze Illumina Infinium array methylation data with R, as well as the report for the final examination of the DNA-RNA dynamics course of the Bioinformatics Master of Bologna University, evaluated with a grade of 30L.

Written by: Giacomo Orsini

# Index
- Introduction
  - DNA Methylation
  - Illumina methylation assay
  - Aim of the project 
- Workflow
  - Load raw data
  - Create R/G dataframes
  - Check probe info by address
  - Create the object MSet.raw
  - Quality check
  - Beta and M values
  - Functional normalization
  - PCA
  - Differential Methylation Analysis
  - Multiple test correction
  - Volcano and Manhattan plots
  - Heatmap

# Introduction
# DNA methylation
**DNA methylation** is an epigenetic modification where a methyl group (–CH₃) is added to the cytosine base in DNA, typically at cytosine-guanine (CpG) dinucleotides. This modification can regulate gene expression without changing the underlying DNA sequence. Methylation usually occurs in clusters called **CpG islands**, which are regions with a high frequency of CpG sites and are often located near gene promoters. When CpG islands are methylated, they can suppress gene expression by preventing the binding of transcription factors or recruiting proteins that compact DNA, making it less accessible.

# Illumina methylation arrays
Methylation arrays are a common tool for detecting methylation levels across the genome. They work by measuring the amount of methylation at specific CpG sites in the DNA. The **Illumina HumanMethylation450 BeadChip** (often referred to as the 450K array) is a widely used platform that can analyze over 450,000 CpG sites across the human genome.
The 450K array uses a technology called bisulfite conversion followed by microarray hybridization to detect methylation levels:

1. DNA is treated with bisulfite, which converts unmethylated cytosines into uracils while leaving methylated cytosines unchanged. This results in sequence differences between methylated and unmethylated DNA that the array can later detect.
2. The converted DNA is then hybridized to probes on the array, which are designed to bind specifically to either the methylated or unmethylated version of each CpG site. Each CpG site is represented by two probes: one that detects methylated DNA and one for unmethylated DNA.
3. The array uses fluorescent labels to indicate whether a site is methylated or unmethylated. The relative fluorescence intensity from each probe pair reflects the degree of methylation at each site.
This information is then processed to generate a beta value (ranging from 0 to 1) for each CpG site, where 0 indicates no methylation, and 1 indicates full methylation.

**Red and green fluorescence** represent signals from the two different types of probes that detect methylated and unmethylated DNA:
- Red fluorescence: it is emitted when the methylated probe binds to methylated DNA. After bisulfite treatment, the methylated cytosine remains unchanged, allowing the methylation-specific probe to hybridize to the site. A strong red fluorescence indicates that a high proportion of the CpG sites are methylated, as the probe binds well to its target.
- Green fluorescence: it is emitted when the unmethylated probe binds to unmethylated DNA. During bisulfite conversion, unmethylated cytosines are converted to uracils, which are then detected by the unmethylated-specific probe. A strong green fluorescence suggests a high proportion of unmethylated CpG sites at that particular locus.

By comparing the intensity of red and green fluorescence for each probe pair, we can determine the methylation level at a specific CpG site: predominantly red fluorescence indicates high methylation (most CpG sites at that locus are methylated); predominantly green fluorescence indicates low or no methylation (most CpG sites are unmethylated); intermediate levels indicates partial methylation, where a subset of cells in the sample may have methylated CpG sites while others do not. Methylation levels are typically represented as beta values (from 0 to 1), where 0 indicates no methylation (green), and 1 indicates full methylation (red).

![image](https://github.com/user-attachments/assets/33586534-0b8e-4bad-acb3-a65a52e9f422)

# Aim of the project
The aim of the project was to analyze human methylation DNA data obtained with the Illumina HumanMethylation450 BeadChip. The analysis included specific tests and analysis assigned by the professor. The report also aims to answer some specific theory questions. The scripts used in the analysis can be found in the `report_GO.R` file while in this `README` file the logic behind each step is explained.

