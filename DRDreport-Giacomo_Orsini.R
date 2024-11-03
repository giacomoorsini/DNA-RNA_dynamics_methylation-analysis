#' ---
#' title: "DRDII final-report"
#' author: "Giacomo Orsini"
#' date: " DNA-RNA dynamics, Jun 28th 2023"
#' ---

#'
#' ### Introduction
#'This is the final report for module II of the course DNA-RNA dynamics, master of Bioinformatics. The report aims to answer biological questions analyzing an Illumina experiment output.
#'To get started, first make sure that the environment is empty and set the working directory. For our analysis, we used the minfi library.
rm(list=ls())
getwd()
setwd("C:/Users/gianl/Desktop/working_directory/DRD/module_II/Final Report-20230614")
#+ warning=FALSE, message==FALSE
suppressMessages(library(minfi))

#' ### 1) Load raw data with minfi and create an object called RGset storing the RGChannelSet object 

#' We want to download the Sample sheet document from our working directory. To do this, we have to set our base directory to the "Input" directory with the **baseDir** command. We then load the .csv file with the **read.metharray.sheet** function, storing it in a variable:
baseDir <- ("C:/Users/gianl/Desktop/working_directory/DRD/module_II/Final Report-20230614/Input") 
targets <- read.metharray.sheet(baseDir)
targets

#' We create and save an object of class RGChannelSet using the function **read.metharray.exp**. This class represents raw data from an Illumina methylation array. We save the file for future use.
RGset <- read.metharray.exp(targets = targets)
save(RGset,file="RGset.RData")
RGset

#' ### 2) Create the dataframes Red and Green to store the red and green fluorescences respectively 

#' To store the Green and Red channels we use the **getRed/getGreen** functions, which extract the data from the RGset object.
Red <- data.frame(getRed(RGset))
Green <- data.frame(getGreen(RGset))

#' ### 3) What are the Red and Green fluorescences for the address assigned to you? 

#' For this step, we have to load a file containing information about names and locations of targeted reference regions, the Illumina Manifest file. 
load("C:/Users/gianl/Desktop/working_directory/DRD/module_II/Final Report-20230614/Illumina450Manifest_clean.RData")
head(Illumina450Manifest_clean)

#' To retrieve the type of probe we check both the addresses section of the Manifest file. In our cases the ID was only for the Address A, while the Address B is empty. We retrieve the information form the "Infinium_Design_Type" section.
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressA_ID=="45652402",]
Illumina450Manifest_clean[Illumina450Manifest_clean$AddressB_ID=="45652402",]

#' To retrieve the red and green fluorescences we search for our address in the Red and Green objects.
Red[rownames(Red)=="45652402",]
Green[rownames(Green)=="45652402",]

#' Now fill the table:
#+ echo=FALSE, results= 'asis'
suppressMessages(library(knitr))
tab <- matrix(c('R01C01','9754','1014','I','RED',
                'R02C01','10517','1288','I','RED',
                'R03C01','9276','1410','I','RED',
                'R04C01','10325','1731','I','RED',
                'R01C02','2639','1376','I','RED',
                'R02C02','2550','1399','I','RED',
                'R03C02','6758','1315','I','RED',
                'R04C02','2536','1618','I','RED'), ncol=5, byrow=TRUE)
colnames(tab) <- c('Sample', 'Red fluor', 'Green fluor', 'Type','Color')
rownames(tab) <- c('A', 'B','C','D','E','F','G','H')
tab <- as.table(tab)
kable(tab)

#' ### 4) Create the object MSet.raw 

#' The **preprocessRaw** converts Red/Green channels fluorescenses into methylation signal (whitohout normalization). We will use it to perform quality checks in the following point.
MSet.raw <- preprocessRaw(RGset)

#' ### 5)	Perform the following quality checks and provide a brief comment to each step: 
#' #### QCplot 

#' A QC plot allows to identify biases in the data. we consider the median of methylation and unmethylation channels for each sample, using the function getQC.
qc <- getQC(MSet.raw)
qc
plotQC(qc)
#' The resulting graph indicates that the data have good quality, because the samples have both high methylation and unmethylation signals, way above the line.
#'
#' #### Check the intensity of negative controls using minfi
#' Then, we check the intensities of the negative control probes. Negative controls target bisulfite-converted sequences that do not contain CpG dinucleotides. The Illumina guide tells us the expected intensities for each type of control probes.
#' To know the names and the numbers of the different types of control probes, we apply the **getProbeInfo** function to the RGset object.

getProbeInfo(RGset, type = "Control")
df_TypeControl <- data.frame(getProbeInfo(RGset, type = "Control"))

#' We can use the **controlSripPlot** function to plot the intensity values of each type of controls probes in our samples
controlStripPlot(RGset, controls="NEGATIVE")
#' Negative controls are all fine, as all below 1000 (log2(1000)=10). Note that there is a minor error in the package: green and red labels are swapped.
#' 
#' #### Calculate detection pValues. For each sample, how many probes have a detection p-value higher than the threshold assigned to you?
#' The detection p-value is computed from the background model. It indicates the chance that the target sequence signal is distinguishable from the negative controls. It is as an objective measure of the overall probe performance. 
#' To calculate the p-value, we use the **detectionP** function 
detP <- detectionP(RGset)
save(detP,file="detP.RData")
head(detP)

#' We input a p-value threshold of 0.01. **summary(failed)** is particularly useful, as it returns the number of failed (TRUE) and not failed (FALSE) positions for each sample.
failed <- detP>0.01
table(failed)
summary(failed) 

#' Now we fill in the table:
#+ echo=FALSE, results= 'asis'
suppressMessages(library(knitr))
tab <- matrix(c('R01C01','56',
                'R02C01','39',
                'R03C01','42',
                'R04C01','35',
                'R01C02','242',
                'R02C02','186',
                'R03C02','21',
                'R04C02','550'), ncol=2, byrow=TRUE)
colnames(tab) <- c('Sample', 'Failed positions')
rownames(tab) <- c('A', 'B','C','D','E','F','G','H')
tab <- as.table(tab)
kable(tab)

#' ### 6) Calculate raw beta and M values and plot the densities of mean methylation values, dividing the samples in WT and MUT. Do you see any difference between the two groups?

#' To retrieve the Beta and M values, we use the functions **getBeta** and **getM**.
beta <- getBeta(MSet.raw)
summary(beta)
M <- getM(MSet.raw)
summary(M)
#' **Note**: Beta values range from 0 to 1 and M values from -inf to +inf; there are although some positions which result as NA:these are the positions for which both the Methylation and Unmethylation values in the MSet.raw are equal to 0. In later analysis, we will have to deal with this values by removing them.
#' 
#' As requested, we divide the "beta" and "M" sets into two subsets considering MUTant and Wild Type groups.
#' To do that, we can load the pheno dataframe, that contains the data in the sample sheet ordered in levels. Specifically, pheno$Group shows the levels of the Group column, and how they are distributed along the sample.

pheno <- read.csv("C:/Users/gianl/Desktop/working_directory/DRD/module_II/Final Report-20230614/Input/Samplesheet_report_2023.csv",header=T, stringsAsFactors=T)
pheno$Group
betaWT <- beta[,pheno$Group=="WT"]
betaMUT <- beta[,pheno$Group=="MUT"]
MWT <- M[,pheno$Group=="WT"]
MMUT <- M[,pheno$Group=="MUT"]

#' We have to calculate the mean of each row of each matrix across the 8 samples. To calculate the mean, we will use the **mean()** and the **apply()** functions, stripping the NA values. The “MARGIN” argument in the function allows to specify that the function mean() should be applied to each row (set “MARGIN” to 1) or to each column (set “MARGIN” to 2) of the matrix/dataframe.
mean_of_betaWT <- apply(betaWT,1,mean,na.rm=T)
mean_of_betaMUT <- apply(betaMUT,1,mean,na.rm=T)
mean_of_MWT <- apply(MWT,1,mean,na.rm=T)
mean_of_MMUT <- apply(MMUT,1,mean,na.rm=T)

#' Now we can calculate the density distributions and plot them.
d_mean_of_betaWT <- density(mean_of_betaWT,na.rm=T)
d_mean_of_betaMUT <- density(mean_of_betaMUT,na.rm=T)
d_mean_of_MWT <- density(mean_of_MWT,na.rm=T)
d_mean_of_MMUT <- density(mean_of_MMUT,na.rm=T)

par( mfrow = c(1, 2))

#Density of Beta Values plot
plot(d_mean_of_betaWT,main="Density of WT Beta Values",col="orange")
lines(d_mean_of_betaMUT,main="Density of MUT Beta Values",col="green")

#Density of M Values plot
plot(d_mean_of_MWT,main="Density of WT M Values",col="orange")
lines(d_mean_of_MMUT,main="Density of MUT M Values",col="green")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE) 
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n') 
legend('bottom',legend = c("WT","MUT"), col = c("orange","green"), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 0.8, seg.len=1, bty = 'n')

#' The plots indicate that there are no significant differences between Mutant and wild type groups, except for a minimal higher B values density at both low and high levels for the WT group. The same minimal difference can be observed in the density of M values.
#'
#' ### 7) Normalize the data using the function assigned to each student and compare raw data and normalized data. 
#' #### Produce a plot with 6 panels in which, for both raw and normalized data, you show the density plots of beta mean values according to the chemistry of the probes, the density plot of beta standard deviation values according to the chemistry of the probes and the boxplot of beta values. Provide a short comment about the changes you observe.  

#' My assigned function was the **preprocessQuantile** function. 
#' To show the differences between normalized data and raw data for the chemistry of our probes, we have to subgroup the original dataset in two sub sets, type I and type II. 
#' We subset the Illumina450Manifest_clean in two dataframes, containing only type I (dfI) or type II (dfII) probes.
dfI <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="I",]
dfI <- droplevels(dfI)
dfII <- Illumina450Manifest_clean[Illumina450Manifest_clean$Infinium_Design_Type=="II",]
dfII <- droplevels(dfII)

#' Now we subset the beta matrix in order to retain only the rows whose name is in the first column of dfI or dfII.
beta_I <- beta[rownames(beta) %in% dfI$IlmnID,]
beta_II <- beta[rownames(beta) %in% dfII$IlmnID,]

#' Then we calculate the mean and standard deviation of the beta values of the new subsets and we retrieve the densities.
mean_of_beta_I <- apply(beta_I,1,mean)
mean_of_beta_II <- apply(beta_II,1,mean)
d_mean_of_beta_I <- density(mean_of_beta_I,na.rm=T)
d_mean_of_beta_II <- density(mean_of_beta_II,na.rm=T)

sd_of_beta_I <- apply(beta_I,1,sd,na.rm=T)
sd_of_beta_II <- apply(beta_II,1,sd,na.rm=T)
d_sd_of_beta_I <- density(sd_of_beta_I,na.rm=T)
d_sd_of_beta_II <- density(sd_of_beta_II,na.rm=T)


#' We normalize the data with our assigned function.
preprocessQuantile_results <- preprocessQuantile(RGset)
beta_preprocessQuantile <- getBeta(preprocessQuantile_results)

#' Then we retrieve the normalized data for our 2 subsets...
beta_preprocessQuantile_I <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfI$IlmnID,]
beta_preprocessQuantile_II <- beta_preprocessQuantile[rownames(beta_preprocessQuantile) %in% dfII$IlmnID,]

#' ...and calculate again the mean, SD and densities.
mean_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,mean)
mean_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,mean)
d_mean_of_beta_preprocessQuantile_I <- density(mean_of_beta_preprocessQuantile_I,na.rm=T)
d_mean_of_beta_preprocessQuantile_II <- density(mean_of_beta_preprocessQuantile_II,na.rm=T)

sd_of_beta_preprocessQuantile_I <- apply(beta_preprocessQuantile_I,1,sd)
sd_of_beta_preprocessQuantile_II <- apply(beta_preprocessQuantile_II,1,sd)
d_sd_of_beta_preprocessQuantile_I <- density(sd_of_beta_preprocessQuantile_I,na.rm=T)
d_sd_of_beta_preprocessQuantile_II <- density(sd_of_beta_preprocessQuantile_II,na.rm=T)

#' Now we can construct the requested plots with the usage of the pheno dataframe,specifically the pheno$Group.
pheno$Group
par(oma = c(4,1,1,1), mfrow = c(2, 3), mar = c(3, 3, 3, 3))

# Density plot of raw B values
plot(d_mean_of_beta_I,col="blue",main="raw beta")
lines(d_mean_of_beta_II,col="red")

# Density plot of raw SD values
plot(d_sd_of_beta_I,col="blue",main="raw sd")
lines(d_sd_of_beta_II,col="red")

# Box plot of raw B values
boxplot(beta, ylab="probes mean",xlab="samples", main="boxplot raw beta", col=c("orange","green","orange","orange","green","green","orange","green"),names=levels(pheno$Sample))

# Density plot of normalized B values
plot(d_mean_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile beta")
lines(d_mean_of_beta_preprocessQuantile_II,col="red")

# Density plot of normalized SD values
plot(d_sd_of_beta_preprocessQuantile_I,col="blue",main="preprocessQuantile sd")
lines(d_sd_of_beta_preprocessQuantile_II,col="red")

# Box plot of normalized B values
boxplot(beta_preprocessQuantile, ylab="probes mean",xlab="probes mean",main="boxplot normalized beta",col=c("orange","green","orange","orange","green","green","orange","green"),names=levels(pheno$Sample))
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE) 

plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n') 
legend('bottom',legend = c("Type I","Type II","WT","MUT"), col = c("blue","red","orange","green"), lwd = 5, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=1, bty = 'n')
#' From this representation, we can see that data have been efficiently normalized, although not perfectly. Hence, our assigned function (preprocessQuantile) can be considered as suited. 
#'
#' ###  8) 	Perform a PCA on the matrix of normalized beta values generated in step 7, after normalization. Comment the plot
#' PCA can be used as a diagnostic plot for the detection of outliers and batch effects. We used the function **prcomp()** to calculate the PCA on our matrix of beta values.
pca_results <- prcomp(t(beta_preprocessQuantile),scale=T)
#' We can print and plot the variance accounted for each component.
par(mar=c(8,8,8,8))
plot(pca_results, main="PCA results")

#' The principal components of interest are stored in the element named “x” of the list.
pca_results$x
#' We can plot PC1 and PC2 and check if samples cluster according to some variable (Group, Sex, batch).
par(oma = c(4,1,1,1), mfrow = c(2, 2), mar=c(2, 2, 2, 2))
plot(pca_results$x[,1], pca_results$x[,2],cex=1,pch=2,xlab="PC1",ylab="PC2",main="PCA",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=levels(pheno$Sample),cex=0.8,pos=1)
pheno$Group
palette(c("orange","green"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=17,col=pheno$Group,main="PCA by group",xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=levels(pheno$Sample),cex=0.8,pos=1)
legend("bottomright",legend=levels(pheno$Group),col=c(1:nlevels(pheno$Group)),pch=17)
pheno$Sex
palette(c("pink","aquamarine"))
plot(pca_results$x[,1], pca_results$x[,2],cex=2,pch=17,col=pheno$Sex,main="PCA by sex",xlab="PC1",ylab="PC2",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=levels(pheno$Sample),cex=0.8,pos=1)
legend("bottomright",legend=levels(pheno$Sex),col=c(1:nlevels(pheno$Sex)),pch=17)
pheno$Sentrix_ID
plot(pca_results$x[,1], pca_results$x[,2],main="PCA by batch",xlab="PC1",ylab="PC2",cex=2,pch=17,col="purple",xlim=c(-1000,1000),ylim=c(-1000,1000))
text(pca_results$x[,1], pca_results$x[,2],labels=levels(pheno$Sample),cex=0.8,pos=1)
legend("bottomright",legend="200400320115",col="purple",pch=17)

#' The graphs suggest that data might be clusterized by sex, but not by group, as we can observe 4 clusters (3 of female and 1 of males). All the samples belong to the same batch, so clustering by batch did not have lots of sense.
#'
#' ### 9)	Using the matrix of normalized beta values generated in step 7, identify differentially methylated probes between group WT and group MUT using the function assigned to each student

#' My assigned function was the **t-test**. To identify the differentially methylated probes, we need to run the t-test on each row of the WT and MUT dataframes to extract the p-value. To do this, we will use an ad hoc simple R function (**My_ttest_function**)
#'
#' We can test our t-test on one row first.
t_test <- t.test(beta_preprocessQuantile[1,] ~ pheno$Group)
t_test
t_test$p.value

#' To run it on the whole dataframe, and to calculate the value between group WT and group MUT, we need to create a formula that takes a probe (x) and computes the t test between the two groups.
My_ttest_function <- function(x) {
  t_test <- t.test(x~ pheno$Group)
  return(t_test$p.value)
} 
#' Now, we run the formula over each line of the normalized beta values dataframe.
pValues_ttest <- apply(beta_preprocessQuantile,1, My_ttest_function)

final_ttest <- data.frame(beta_preprocessQuantile, pValues_ttest)
head(final_ttest)

#' We can order the probes on the basis of the pValues column (from the smallest to the largest value).
final_ttest <- final_ttest[order(final_ttest$pValues_ttest),]
head(final_ttest)

#' Finally, we can count how many probes have a p-Value smaller than 0.05.
final_ttest_0.05 <- final_ttest[final_ttest$pValues_ttest<=0.05,]
dim(final_ttest_0.05)

#' ### 10)	Apply multiple test correction and set a significant threshold of 0.05. How many probes do you identify as differentially methylated considering nominal pValues? How many after Bonferroni correction? How many after BH correction?

#' To apply a multiple test correction we will use the **p.adjust** function. We have to apply both the Bonferroni and the Benjamini Hochberg correction, we can store the output in a dataframe.
corrected_pValues_BH <- p.adjust(final_ttest$pValues_ttest,"BH")
corrected_pValues_Bonf <- p.adjust(final_ttest$pValues_ttest,"bonferroni")
final_ttest_corrected <- data.frame(final_ttest, corrected_pValues_BH, corrected_pValues_Bonf)
head(final_ttest_corrected)

#' We can visualize the distributions of the p-values and of the corrected p-values by boxplots; the red line indicates our threshold:
par(mfrow = c(1, 1),mar = c(2, 2, 2, 2))
boxplot(final_ttest_corrected[,9:11])
abline(h=0.05,col="red")

#' Number of differentially expressed probes before correction:
dim(final_ttest_corrected[final_ttest_corrected$pValues_ttest<=0.05,])[1]
#' Number of differentially expressed probes after BH correction:
dim(final_ttest_corrected[final_ttest_corrected$corrected_pValues_BH<=0.05,])[1]
#' Number of differentially expressed probes after Bonferroni correction:
dim(final_ttest_corrected[final_ttest_corrected$corrected_pValues_Bonf<=0.05,])[1]

#' ### 11)	Produce a volcano plot and a Manhattan plot of the results of differential methylation analysis 

#' Hereby all the steps to produce a volcano plot
# Create two matrices for normalized beta values of WT and MUT and calculate the means
beta_volcano <- final_ttest_corrected[,1:8]
beta_volcano_WT <- beta_volcano[,pheno$Group=="WT"]
mean_beta_volcano_WT<- apply(beta_volcano_WT,1,mean)
beta_volcano_MUT <- beta_volcano[,pheno$Group=="MUT"]
mean_beta_volcano_MUT <- apply(beta_volcano_MUT,1,mean)

# Calculate the distance between the means
delta_volcano <- mean_beta_volcano_WT-mean_beta_volcano_MUT

# Create a 2 columns data frame containing delta values and -log(p-value) values
toVolcPlot <- data.frame(delta_volcano, -log10(final_ttest_corrected$pValues_ttest))
head(toVolcPlot)

# Make the Volcano plot
plot(toVolcPlot[,1], toVolcPlot[,2],pch=16,cex=0.5,xlab="Δ mean",ylab="-log10(pvalue)")
abline(h=-log10(0.01),col="red")

toHighlight <- toVolcPlot[toVolcPlot[,1]>0.2 & toVolcPlot[,2]>(-log10(0.01)),]
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="paleturquoise")
toHighlight <- toVolcPlot[toVolcPlot[,1]<(-0.2) & toVolcPlot[,2]>(-log10(0.01)),]
points(toHighlight[,1], toHighlight[,2],pch=16,cex=0.7,col="palegreen")
#' Higlightened in light blue the hypermethylated probes, in green the hypomethilated, while in red the significance threshold.
#' 
#' Hereby all the steps to produce a Manhattan plot with the library *qqman*
# Install the library
suppressMessages(library(qqman))

# Retrieve the t test corrected data frame 
final_ttest_reduced <- data.frame(rownames(final_ttest_corrected),final_ttest_corrected)

# Retrieve the column names
colnames(final_ttest_reduced)
colnames(final_ttest_reduced)[1] <- "IlmnID"
colnames(final_ttest_reduced)

# Merge the reduced corrected file with the illimunia manifest file.
final_ttest_reduced_annotated <- merge(final_ttest_reduced, Illumina450Manifest_clean,by="IlmnID")
dim(final_ttest_reduced)
dim(final_ttest_reduced_annotated)

# The Manhattan plot input has to contain 4 items: probe, chromosome, position on the chromosome and p-value. Select them...
input_Manhattan <- final_ttest_reduced_annotated[colnames(final_ttest_reduced_annotated) %in% c("IlmnID","CHR","MAPINFO","pValues_ttest")]
levels(input_Manhattan$CHR)

#...and order the chromosomes. Convert also the Y and X chromosomes into numbers.
order_chr <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
input_Manhattan$CHR <- factor(input_Manhattan$CHR,levels=order_chr )
input_Manhattan$CHR <- as.numeric(input_Manhattan$CHR)
table(input_Manhattan$CHR)

# Create the plot
par(mfrow = c(1, 1),mar = c(2, 2, 2, 2))
manhattan(input_Manhattan, snp="IlmnID",chr="CHR", bp="MAPINFO", p="pValues_ttest",annotatePval = 0.00001,col=rainbow(24) )
#' In blue our threshold, highlightened the samples with value higher than it.
#' 
#' ### 12) Produce an heatmap of the top 100 differentially methylated probes 

#' To create an heatmap, we will use the **heatmap.2** function of the *gplots* library. It takes as input a matrix.
#' 
#' As asked in the exercise title, we will use just the top 100 most significant CpG probes to create some heatmaps. The heat maps were obtained using the **gplots** library and different hierchical clustering methods.
suppressMessages(library(gplots))
input_heatmap=as.matrix(final_ttest[1:100,1:8])
colorbar <- c("orange","green","orange","orange","green","green","orange","green")
# complete linkage hierarchical clustering
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Complete linkage")
# single linkage hierarchical clustering
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'single'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Single linkage")
# average linkage hierarchical clustering
heatmap.2(input_heatmap,col=terrain.colors(100),Rowv=T,Colv=T,hclustfun = function(x) hclust(x,method = 'average'),dendrogram="both",key=T,ColSideColors=colorbar,density.info="none",trace="none",scale="none",symm=F,main="Average linkage")
#' Looking at the results, we can see that the different linkage methods led to different heatmaps, in particular the complete linkage heat map is different from the other two. Considering our type of analysis, an average linkage method would be better. In all the cases, the group are correctly clustered. We can see that there are some probes that are methylated/unmethylated depending on the groups.