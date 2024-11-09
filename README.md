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
- Green fluorescence: it is emitted when the methylated probe binds to methylated DNA. After bisulfite treatment, the methylated cytosine remains unchanged, allowing the methylation-specific probe to hybridize to the site. A strong green fluorescence indicates that a high proportion of the CpG sites are methylated, as the probe binds well to its target.
- Red fluorescence: it is emitted when the unmethylated probe binds to unmethylated DNA. During bisulfite conversion, unmethylated cytosines are converted to uracils, which are then detected by the unmethylated-specific probe. A strong red fluorescence suggests a high proportion of unmethylated CpG sites at that particular locus.

By comparing the intensity of red and green fluorescence for each probe pair, we can determine the methylation level at a specific CpG site: predominantly red fluorescence indicates high methylation (most CpG sites at that locus are methylated); predominantly green fluorescence indicates low or no methylation (most CpG sites are unmethylated); intermediate levels indicates partial methylation, where a subset of cells in the sample may have methylated CpG sites while others do not. Methylation levels are typically represented as beta values (from 0 to 1), where 0 indicates no methylation (green), and 1 indicates full methylation (red).

![image](https://github.com/user-attachments/assets/33586534-0b8e-4bad-acb3-a65a52e9f422)

# Aim of the project
The aim of the project was to analyze human methylation DNA data obtained with the Illumina HumanMethylation450 BeadChip. The analysis included specific tests and analysis assigned by the professor. The report also aims to answer some specific theory questions. The scripts used in the analysis can be found in the `report_GO.R` file while in this `README` file the logic behind each step is explained.

# Workflow
The workflow hereby reported can be better followed in the `report_GO.html` or `report_GO.R` files, as no scripts will be provided here.

## 1) Load Raw Data with `minfi` and Create an RGset Object
To start our analysis, we need to load the raw data from the Illumina methylation array and store it in an appropriate format for further processing. We use the `mini` package, which is designed for handling Illumina methylation data, for creating an **RGChannelSet object** (named RGset in our case). This object stores the red and green fluorescence intensities from the methylation array, representing the raw, unprocessed data.

First, we define the base directory where our input files are located. We then load the sample sheet, a metadata file that contains essential details about each sample in the experiment, including sample names, IDs, and phenotypic information. Using the `read.metharray.exp` function, we create the RGChannelSet object from the raw array data files and save it for future use, ensuring that we can always return to the original data if needed. This step prepares the data for preprocessing and analysis by loading it into a standard structure compatible with downstream steps.

## 2) Create Dataframes for Red and Green Fluorescences
In this step, we extract the red and green fluorescence channels from the RGset object. These fluorescence intensities are crucial because they reflect the methylation status of specific CpG sites: red fluorescence generally corresponds to methylated DNA, while green fluorescence corresponds to unmethylated DNA. By isolating these channels into separate dataframes (named Red and Green), we simplify the data structure and make it easier to inspect, manipulate, and analyze.

Using the `getRed` and `getGreen` functions from the `minfi` package, we extract the respective fluorescence channels. This separation is an essential step before we proceed with normalization or other transformations, as it allows us to work directly with the fluorescence data for quality control and further processing.

## 3) Determine Red and Green Fluorescences for a Specific Address
To understand the fluorescence intensities for a given probe address (specific to our samples), we need to know where each probe maps in the genome and whether it is designed to detect methylated or unmethylated states. We use the `Illumina Manifest` file, a reference file containing information on the location, target region, and type of each probe. This file tells us whether a probe is of `Infinium I` or `Infinium II` design, each of which has a different way of handling methylated and unmethylated signals.

For this particular analysis, we examine the probe addresses to identify which are assigned to Address A (used in our dataset), and we check the `Infinium_Design_Type` to confirm the type of probe. Once we identify the relevant probes, we use the Red and Green objects to retrieve the fluorescence intensities for the specific addresses. This step helps us understand the initial fluorescence data for specific probes, enabling us to monitor methylation levels at each targeted site.

## 4) Create the `MSet.raw` Object
In this step, we convert the raw red and green fluorescence signals into a methylation signal using the `preprocessRaw` function from `minfi`. This function combines the fluorescence data into a raw methylation measurement (without normalization), creating an object called `MSet.raw`.

The `MSet.raw` object provides the beta values and methylated/unmethylated intensities for each CpG site, which we can use to assess the methylation state of each site across our samples. This raw preprocessing step is essential for initial quality checks; it allows us to inspect the data before any normalization or adjustment, helping us identify any technical artifacts, batch effects, or quality issues in the raw data. Quality checks on `MSet.raw` ensure that we are working with reliable data and set the foundation for accurate downstream analysis.

## 5) Perform Quality Checks
Quality control (QC) is essential in methylation analysis to ensure that the data are reliable and free of technical artifacts or biases. Below are the key quality checks we performed, with explanations for each step:

- **QC Plot** :To assess the overall quality of the data, we generate a QC plot using the `getQC` function from the minfi package. This plot examines the median methylation and unmethylation signal intensities for each sample. Samples with high-quality data typically show strong signals in both the methylated and unmethylated channels, well above the background noise level. In our case, the plot shows that all samples have high median signals, indicating good-quality data with no evident technical issues or significant biases.

- **Check the Intensity of Negative Control Probes**: Another important QC step is to evaluate the negative control probes included in the Illumina array. These probes target non-CpG regions and are bisulfite converted, serving as a baseline to detect background noise. The intensities of these controls help us verify that there is no unintended signal from sequences that should not be methylated. According to the Illumina guide, each type of control probe has an expected intensity range. Using the `getProbeInfo` function on the RGset object, we identify the names and counts of each type of control probe, which allows us to check them individually.
We then use the `controlStripPlot` function to visualize the intensity values for each type of control probe in our samples. Ideally, the intensities for negative controls should be below a certain threshold (in our case, below 1000, or log2(1000) = 10). The plot shows that all negative control intensities are within acceptable limits, indicating minimal background noise. However, there is a minor labeling issue in the minfi package where the green and red labels are swapped, so caution is advised in interpreting the colors.

- **Calculate Detection p-Values**: Detection p-values provide an objective measure of the quality of each probe in each sample by calculating the probability that a detected signal is above the background noise level. Low detection p-values (below a threshold, such as 0.01) indicate that the signal is reliable, while high p-values suggest that the probe might be unreliable or indistinguishable from background noise.

We calculate detection p-values using the `detectionP` function, applying a p-value threshold of 0.01 to filter out low-quality probes. We then use the `summary(failed)` function to determine the number of failed (p-value above threshold) and passed probes for each sample. This check helps us identify samples or probes that may require additional filtering, ensuring that only reliable data are used in the downstream analysis.

## 6) Calculate Raw Beta and M Values and Plot the Densities of Mean Methylation Values, Dividing the Samples in WT and MUT
To retrieve the Beta and M values, we use the functions `getBeta` and `getM`. Note: Beta values range from 0 to 1, and M values range from -∞ to +∞. Some positions may result in NA values; these occur when both the methylation and unmethylation values in MSet.raw are zero. In later analyses, we will address these NA values by removing them.
As requested, we divide the "Beta" and "M" sets into two subsets based on Mutant (MUT) and Wild Type (WT) groups. To achieve this, we load the pheno dataframe, which contains the data from the sample sheet with the appropriate group labels. Specifically, pheno$Group shows the distribution of the Group column values across samples.
To calculate the mean of each row across the 8 samples, we use the mean() and apply() functions, ensuring that we strip any NA values. We set the MARGIN argument to 1 to apply the mean() function across rows (set to 2 for columns).
With the mean values calculated, we generate the density distributions and plot them.

The resulting plots indicate that there are no significant differences between the Mutant and Wild Type groups. However, there is a slight increase in Beta value density at both low and high levels for the WT group. This minimal difference is also observed in the density of M values.
