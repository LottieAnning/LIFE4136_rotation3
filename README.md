# LIFE4136 Rotation 3
Welcome to the GitHub repository dedicated to Rotation 3 of LIFE4136, where we delve into the intriguing realm of allopolyploidy following whole-genome duplication (WGD). Building upon the groundwork laid by Marburger et al. in 2019 (accessible via this link: https://www.nature.com/articles/s41467-019-13159-5), this platform houses an exploration of discernible patterns indicating the presence of allopolyploidy within hybrid populations.

Our research aims to provide insights crucial to the ongoing investigations led by the Marbuger et al. team, unraveling the mechanisms and significance of allopolyploidy arising from interspecific hybridization post-WGD.

Arabidopsis arenosa (referred to as 'arenosa' henceforth) stands as the sister species to Arabidopsis lyrata ('lyrata' henceforth). Both species manifest in diploid and tetraploid forms, with arenosa preceding lyrata in undergoing WGD. The post-WGD hybridization of lyrata with arenosa tetraploids, though occurring subsequently, introduces a degree of instability due to the infusion of novel genetic elements. Those hybrids that retain sub-genomes reminiscent of their parent species are recognized as allopolyploids.

Anticipated is the gradual reduction of allotetraploids to a diploid state over time, thus fostering increased genetic diversity within the Arabidopsis taxa.

The below code is to be ran in alternating R, Python and UNIX environments

## Contents
- [Dependencies](#dependencies)
- [Files Required](#files_required)
- [Initial Visualisation of Data](#initial_visualisation_of_data)
  - [Converting a VCF into a genlight object](#vcf_to_genlight)
  - [Sense checking data](#sense_check)
  - [Running an initial PCA](#initial_pca)
- [Filtering Data for Further Analysis of Trends](#filter_data_for_further_analysis_of_trends)
  - [Filtering with GATK](#gatk)
  - [Downloading files from a HPC](#hpc)
  - [Re-running the PCA](#second_pca)
  - [Testing geographical influence](#map)
  - [Running a PCA on individuals in a population](#third_pca)
  - [Further filtering with GATK](#further_gatk)
- [Relatedness Calculations](#relatedness_calculations)
  - [Nei's distances](#nei)
  - [SplitsTree](#splitstree)
- [Fast Structure](#fast_structure)
  - [Structure plotting](#structure_plot)
- [Allele Frequency Spectrum](#allele_frequency_spectrum)
  - [Creating an example allotetraploid](#allotetraploid)
  - [Creating allele frequency histograms](#histogram)

## Dependencies

<a name="dependencies"/>

* To run *Principle Component Analysis*, *Sample Mapping*, *Nei's Distance Calculations* and *Allele Frequency Spectrum Plots* **R Studio version 4.3.3** is needed, this can be downloaded here: https://cran.r-project.org/mirrors.html, simply navigate to your country and select the package compatible for your machine. Secondary to this, the following R packages need to be installed for the code to run:
  * **adegenet** version 2.1.10 or higher, simply install by typing: ```install.pacakges(adegenet)``` into the R command line
  * **adegraphics** version 1.0.21 or higher, simply install by typing: ```install.pacakges(adegraphics)``` into the R command line
  * **dplyr** version 1.1.4 or higher, simply install by typing: ```install.pacakges(dplyr)``` into the R command line
  * **ggplot2** version 3.4.4 or higher, simply install by typing: ```install.pacakges(ggplot2)``` into the R command line
  * **ggrepel** version 0.9.5 or higher, simply install by typing: ```install.pacakges(ggrepel)``` into the R command line
  * **leaflet** version 2.2.1 or higher, simply install by typing: ```install.pacakges(leaflet)``` into the R command line
  * **StAMPP** version 1.6.3 or higher, simply install by typing: ```install.pacakges(StAMPP)``` into the R command line
  * **tidyr** version 1.3.0 or higher, simply install by typing: ```install.pacakges(tidyr)``` into the R command line
  * **vcfR** version 1.15.0 or higher, simply install by typing: ```install.pacakges(vcfR)``` into the R command line
* To run the scripts for *fastStructure* **Python version 2.7.18** is required, to enable the use of this, create a virtual enviornmnet: ```conda create -y -n faststructure python=2.7.18```
  * **faststructure** version 0.0.0, which can be cloned here: git clone https://github.com/rajanil/fastStructure
  * **scipy** version 1.2.1, which can be installed in a python 2 environment by: pip install scipy
  * **numpy** version 1.16.5, which can be installed in a python 2 environment by: pip install numpy
  * **utils** version 0.9.0, which can be installed in a python 2 environment by: pip install utils 
  * **parse_bed**
  * **parse_str**
* To run the *Cochlearia_create_structure_file.py* script for *fastStructure* and *downloading files from a High Power Computer (HPC)* **Python version 3.8.12** is required along with:
  * .
* To run *gatk* **version 4.2.2.0** is required, first create a virtual environment: ```conda create --name /[path_to_virtual_environment]/[virtual_environment_name]```, then download the package: https://github.com/broadinstitute/gatk/releases, however a HPC with an existing version is recommended as this is a large package. Other dependencies are:
  * **Conda version 23.11.0**
  * **SAMtools version 1.19.2** which can be downloaded here: https://www.htslib.org/download/
* To plot the *Neighbour Joining (NJ) Trees* **Splits tree of version 6.2.2-beta** is required which can be downloaded at:	https://github.com/husonlab/splitstree6
* **structure plot?**

## Files Required

<a name="files_required"></a>

* A **vcf** with all your samples
* A reference **fasta file** to which your reads were aligned to

## Initial Visualisation of Data

<a name="initial_visualisation_of_data"></a>

Load packages in R:
```
library(vcfR)
library(adegenet)
library(adegraphics)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(StAMPP)
library(leaflet)
```
Set working directory:
```
setwd("[path_to_working_directory]")
```
Read in vcf:
```
vcf <- read.vcfR("[title_of_vcf].vcf")
```
### Convert your vcf into a genlight object:

<a name="vcf_to_genlight"></a>

First create the function:
```
vcfR2genlight.tetra <- function (x, n.cores = 1) 
{
  bi <- is.biallelic(x)
  if (sum(!bi) > 0) {
    msg <- paste("Found", sum(!bi), "loci with more than two alleles.")
    msg <- c(msg, "\n", paste("Objects of class genlight only support loci with two alleles."))
    msg <- c(msg, "\n", paste(sum(!bi), "loci will be omitted from the genlight object."))
    warning(msg)
    x <- x[bi, ]
  }
  x <- addID(x)
  CHROM <- x@fix[, "CHROM"]
  POS <- x@fix[, "POS"]
  ID <- x@fix[, "ID"]
  x <- extract.gt(x)
  x[x == "0|0"] <- 0
  x[x == "0|1"] <- 1
  x[x == "1|0"] <- 1
  x[x == "1|1"] <- 2
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/0"] <- 1
  x[x == "1/1"] <- 2
  x[x == "1/1/1/1"] <- 4
  x[x == "0/1/1/1"] <- 3
  x[x == "0/0/1/1"] <- 2
  x[x == "0/0/0/1"] <- 1
  x[x == "0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/0"] <- 0
  x[x == "0/0/0/0/0/1"] <- 1
  x[x == "0/0/0/0/1/1"] <- 2
  x[x == "0/0/0/1/1/1"] <- 3
  x[x == "0/0/1/1/1/1"] <- 4
  x[x == "0/1/1/1/1/1"] <- 5
  x[x == "1/1/1/1/1/1"] <- 6
  if (requireNamespace("adegenet")) {
    x <- new("genlight", t(x), n.cores = n.cores)
  }
  else {
    warning("adegenet not installed")
  }
  adegenet::chromosome(x) <- CHROM
  adegenet::position(x) <- POS
  adegenet::locNames(x) <- ID
  return(x)
}
```
Then run the conversion:
```
aa.genlight <- vcfR2genlight.tetra(vcf)
```
Warnings may occur after this line, you can use the **warnings()** function to view the specifics. In my initial run, 31 warnings occured all of which returned: 'In initialize(value, ...) : NAs introduced by coercion' meaning there was missing data in the vcf so NA values are introduced. You can inspect your VCF file to identify any missing or unexpected data that might be causing the issue, however at this stage the code can handle missing values and visualisation will still show trends.

Create individual id's and population names:
```
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3) 
```
### Check your data - this section is optional:

<a name="sense_check"></a>

```
aa.genlight
```
This will tell you basic information about your genlight object, for example how many individuals in sample.

If you want to see the individual names: ```indNames(aa.genlight)```

If you want to see the populations and number of them, run: ```unique(pop(aa.genlight))```

If you want to see the variety of ploidy in sample, run: ```unique(ploidy(aa.genlight))```

### Run an initial PCA:

<a name="initial_pca"></a>

This will allow for the discovery of potential trends in your data and allow for preliminary visualisation.

Define the PCA function:
```
glPcaFast <- function(x,
                      center=TRUE,
                      scale=FALSE,
                      nf=NULL,
                      loadings=TRUE,
                      alleleAsUnit=FALSE,
                      returnDotProd=FALSE){
  
  if(!inherits(x, "genlight")) stop("x is not a genlight object")
  # keep the original mean / var code, as it's used further down
  # and has some NA checks..
  if(center) {
    vecMeans <- glMean(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecMeans))) stop("NAs detected in the vector of means")
  }
  if(scale){
    vecVar <- glVar(x, alleleAsUnit=alleleAsUnit)
    if(any(is.na(vecVar))) stop("NAs detected in the vector of variances")
  }
  # convert to full data, try to keep the NA handling as similar
  # to the original as possible
  # - dividing by ploidy keeps the NAs
  mx <- t(sapply(x$gen, as.integer)) / ploidy(x)
  # handle NAs
  NAidx <- which(is.na(mx), arr.ind = T)
  if (center) {
    mx[NAidx] <- vecMeans[NAidx[,2]]
  } else {
    mx[NAidx] <- 0
  }
  # center and scale
  mx <- scale(mx,
              center = if (center) vecMeans else F,
              scale = if (scale) vecVar else F)
  # all dot products at once using underlying BLAS
  # to support thousands of samples, this could be
  # replaced by 'Truncated SVD', but it would require more changes
  # in the code around
  allProd <- tcrossprod(mx) / nInd(x) # assume uniform weights
  ## PERFORM THE ANALYSIS ##
  ## eigenanalysis
  eigRes <- eigen(allProd, symmetric=TRUE, only.values=FALSE)
  rank <- sum(eigRes$values > 1e-12)
  eigRes$values <- eigRes$values[1:rank]
  eigRes$vectors <- eigRes$vectors[, 1:rank, drop=FALSE]
  ## scan nb of axes retained
  if(is.null(nf)){
    barplot(eigRes$values, main="Eigenvalues", col=heat.colors(rank))
    cat("Select the number of axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  ## rescale PCs
  res <- list()
  res$eig <- eigRes$values
  nf <- min(nf, sum(res$eig>1e-10))
  ##res$matprod <- allProd # for debugging
  ## use: li = XQU = V\Lambda^(1/2)
  eigRes$vectors <- eigRes$vectors * sqrt(nInd(x)) # D-normalize vectors
  res$scores <- sweep(eigRes$vectors[, 1:nf, drop=FALSE],2, sqrt(eigRes$values[1:nf]), FUN="*")
  ## GET LOADINGS ##
  ## need to decompose X^TDV into a sum of n matrices of dim p*r
  ## but only two such matrices are represented at a time
  if(loadings){
    if(scale) {
      vecSd <- sqrt(vecVar)
    }
    res$loadings <- matrix(0, nrow=nLoc(x), ncol=nf) # create empty matrix
    ## use: c1 = X^TDV
    ## and X^TV = A_1 + ... + A_n
    ## with A_k = X_[k-]^T v[k-]
    myPloidy <- ploidy(x)
    for(k in 1:nInd(x)){
      temp <- as.integer(x@gen[[k]]) / myPloidy[k]
      if(center) {
        temp[is.na(temp)] <- vecMeans[is.na(temp)]
        temp <- temp - vecMeans
      } else {
        temp[is.na(temp)] <- 0
      }
      if(scale){
        temp <- temp/vecSd
      }
      res$loadings <- res$loadings + matrix(temp) %*% eigRes$vectors[k, 1:nf, drop=FALSE]
    }
    res$loadings <- res$loadings / nInd(x) # don't forget the /n of X_tDV
    res$loadings <- sweep(res$loadings, 2, sqrt(eigRes$values[1:nf]), FUN="/")
  }
  ## FORMAT OUTPUT ##
  colnames(res$scores) <- paste("PC", 1:nf, sep="")
  if(!is.null(indNames(x))){
    rownames(res$scores) <- indNames(x)
  } else {
    rownames(res$scores) <- 1:nInd(x)
  }
  if(!is.null(res$loadings)){
    colnames(res$loadings) <- paste("Axis", 1:nf, sep="")
    if(!is.null(locNames(x)) & !is.null(alleles(x))){
      rownames(res$loadings) <- paste(locNames(x),alleles(x), sep=".")
    } else {
      rownames(res$loadings) <- 1:nLoc(x)
    }
  }
  if(returnDotProd){
    res$dotProd <- allProd
    rownames(res$dotProd) <- colnames(res$dotProd) <- indNames(x)
  }
  res$call <- match.call()
  class(res) <- "glPca"
  return(res)
}
```
Run the PCA:
```
pca.1 <- glPcaFast(aa.genlight, nf=300)
```
Create a colour palette:
```
col <- funky(10)
```
Plot the PCA:
```
s.class(pca.1$scores, pop(aa.genlight), xax=1, yax=2, col=transp(col,.6), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, xlab = "PC1", ylab = "PC2")
```
This is what my first PCA looked like:

![First PCA](Figures/first_pca.png)

For reference, this figure from the Marburger paper, *linked at the top of the page,* shows some information on ploidy and purity of the new data (not all populations are included in this plot and SWA and HAL are SWB and HAB in the new data respectively):

![Marburger Plot](Figures/marburger_plot.png)

At first it is difficult to see trends.

Try colouring by ploidy:

```
ploidy_labels <- factor(ploidy(aa.genlight))
s.class(pca.1$scores, ploidy_labels, xax=1, yax=2, col=transp(col,.6), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, xlab = "PC1", ylab = "PC2")
```
![Ploidy PCA](Figures/ploidy_pca.png)

As you can see PC1 seems to relatively seperate diploids and tetraploids and PC2 hybrids from pure lyrata. However, this isnt a perfect pattern, many diploids cluster positively with PC1 (with the tetraploids). Furthermore, the tetraploids include hybrids (mixture of lyrata and arenosa) and a population genetically similar to arenosa (KEH). In the diploid population there are only pure lyrata, this means that the PCA could be unbalanced as there are no arenosa diploids to balance out the arenosa tetraploids.

## Filter Data for Further Analysis of Trends

<a name="filter_data_for_further_analysis_of_trends"></a>

### Filter your vcf with gatk:

<a name="gatk"></a>

In a UNIX command line, activate your gatk environment: 
```
conda activate /[path_to_virtual_environment]/[virtual_environment_name]
```
Index your reference fasta file: 
```
samtools faidx [name_of_your_reference_file].fasta
```
Create a dictionary file from your reference fasta: 
```
gatk CreateSequenceDictionary -R [name_of_your_reference_file].fasta
```
Create a filtered .args file with the populations you wish to include (in this example I filter for only tetraploids, the structure of the code is 'BZD-....' for example because the sample id's in the data are in the format: 'BZD-01tl' where 'BZD' is the population name and the following characters are the individual id): 
```
grep -o -E 'BZD-....|PEK-....|SCT-....|TEM-....|GYE-....|JOH-....|KAG-....|LIC-....|LOI-....|MAU-....|MOD-....|PIL-....|SCB-....|SWB-....|HAB-....|ROK-....|FRE-....|OCH-....|KEH-....' [name_of_your_vcf].vcf > tetraploid.args
```
Filter your vcf:
```
gatk SelectVariants -R [name_of_your_reference_file].fasta -V [name_of_your_vcf].vcf -sn tetraploid.args -O tetraploid.vcf 
```
Deactivate your gatk environment: 
```
conda deactivate
```

### If you are using a HPC download the vcf via a web browser: 

<a name="hpc"></a>

Create the webbrowser: 
```
python3 -m http.server 36895
```
Access your web address and download the file: ```http://[your_HPC_ip]:36895```

### Re-run the PCA:

<a name="second_pca"></a>

In R, read in the vcf: 
```
vcf <- read.vcfR("tetraploid.vcf")
```
Convert to a genlight object:
```
aa.genlight <- vcfR2genlight.tetra(vcf)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_")
pop(aa.genlight)<-substr(indNames(aa.genlight),1,3) 
```
Run a PCA and plot:
```
pca.2 <- glPcaFast(aa.genlight, nf=300)
s.class(pca.2$scores, pop(aa.genlight), xax=1, yax=2, col=transp(col,.6), ellipseSize=0, starSize=0, ppoints.cex=4, paxes.draw=T, pgrid.draw =F, xlab = "PC1", ylab = "PC2")
```
![Tetraploid PCA](Figures/tetraploids_pca.png)

Again PC1 seperates pure lyrata (left) from hybrids (centre) and populations closer to arenosa (right) and PC2 is potentially seperating populations via geography.

### To test if PC2 is seperating populations via geographical location you can use the leaflet() package:

<a name="map"></a>

```
mymap <- leaflet() %>%
  setView(lng = 15.2551, lat = 50, zoom = 5) %>%
  addTiles()

mymap <- mymap %>% addMarkers(lng = 14.722019, lat = 50.533611, popup = "BZD", label = "BZD",labelOptions = labelOptions(noHide = TRUE))
mymap <- mymap %>% addMarkers(lng = 10.582812, lat = 51.583498, popup = "SCT", label = "SCT",labelOptions = labelOptions(noHide = TRUE))
mymap <- mymap %>% addMarkers(lng = 16.248345, lat = 49.090359, popup = "TEM", label = "TEM",labelOptions = labelOptions(noHide = TRUE))
mymap <- mymap %>% addMarkers(lng = 14.482233, lat = 50.116151, popup = "PEK", label = "PEK",labelOptions = labelOptions(noHide = TRUE))
mymap <- mymap %>% addMarkers(lng = 15.611308, lat = 47.876394, popup = "OCH")
mymap <- mymap %>% addMarkers(lng = 17.285185, lat = 46.773573, popup = "GYE", label = "GYE",labelOptions = labelOptions(noHide = TRUE))
mymap <- mymap %>% addMarkers(lng = 16.038565, lat = 47.780388, popup = "JOH")
mymap <- mymap %>% addMarkers(lng = 15.426401, lat = 48.294661, popup = "KAG")
mymap <- mymap %>% addMarkers(lng = 16.269984, lat = 48.092981, popup = "LIC")
mymap <- mymap %>% addMarkers(lng = 15.552699, lat = 48.396128, popup = "LOI")
mymap <- mymap %>% addMarkers(lng = 15.560472, lat = 48.381430, popup = "MAU")
mymap <- mymap %>% addMarkers(lng = 16.267199, lat = 48.079485, popup = "MOD")
mymap <- mymap %>% addMarkers(lng = 15.349241, lat = 48.239088, popup = "PIL")
mymap <- mymap %>% addMarkers(lng = 15.393005, lat = 48.2742817, popup = "SCB")
mymap <- mymap %>% addMarkers(lng = 15.40084, lat = 48.3403, popup = "SWB")
mymap <- mymap %>% addMarkers(lng = 15.698451, lat = 47.937742, popup = "HAB")
mymap <- mymap %>% addMarkers(lng = 15.68304, lat = 47.9053, popup = "ROK")
mymap <- mymap %>% addMarkers(lng = 15.57118, lat = 47.99405, popup = "FRE")
mymap <- mymap %>% addMarkers(lng = 15.542146, lat = 47.816197, popup = "KEH")

mymap
```
SetView() is used to set the central coordinates of the map produced (the ones used in the above coordinates are central Europe) and the zoom is set to 5. Add in the latitude and longitude coordinates for your sample populations and their names. Coordinates for the data used in this study were obtained from the sample map provided, some data was missing so BZD, SCT, TEM, PEK and GYE coordinates were estimated using the description of where they were sampled.

![Annotated Map](Figures/annotated_map.png)
I have annotated the populations that were sampled outside the main cluster i.e. SCT, BZD, PEK, TEM and GYE using the labelOptions() command. On the PCA SCT plots most positvely with PC2 and was sampled the furtherst North, similar with BZD, TEM and PEK so this supports the hypothesis that PC2 seperates the populations by geography. However, GYE is further South from the sample cluster, but it doesn't have the most negative correlation with PC2 (this is LIC, MOD and JOH). Therefore, the relationship which PC2 describes is not fully clear.

### For further analysis you can plot the PCA by individuals:

<a name="third_pca"></a>

```
df <- extract.gt(vcf)
df[df == "0|0"] <- 0
df[df == "0|1"] <- 1
df[df == "1|0"] <- 1
df[df == "1|1"] <- 2
df[df == "0/0"] <- 0
df[df == "0/1"] <- 1
df[df == "1/0"] <- 1
df[df == "1/1"] <- 2
df[df == "0/0/0/0"] <- 0
df[df == "0/0/0/1"] <- 1
df[df == "0/0/1/1"] <- 2
df[df == "0/1/1/1"] <- 3
df[df == "1/1/1/1"] <- 4
df[df == "0/0/0/0/0/0"] <- 0
df[df == "0/0/0/0/0/1"] <- 1
df[df == "0/0/0/0/1/1"] <- 2
df[df == "0/0/0/1/1/1"] <- 3
df[df == "0/0/1/1/1/1"] <- 4
df[df == "0/1/1/1/1/1"] <- 5
df[df == "1/1/1/1/1/1"] <- 6
df <- data.frame(apply(df,2,function(x)as.numeric(as.character(x))))
```
Remove samples with > 50% missing data:
```
mis <- apply(df,2,function(x)sum(is.na(x))/length(x))
df <- df[,mis <= 0.5]
```
Calculate allele frequencies:
```
ploidy <- apply(df,2,max,na.rm=T)
p <- apply(df,1,function(x)sum(x,na.rm=T)/sum(ploidy[!is.na(x)]))
```
Removing individuals can change allele frequencies, so we make sure that maf >= 0.05:
```
df <- df[p >= 0.05 & p <= 0.95,]
p <- p[p >= 0.05 & p <= 0.95]
```
Estimate a covariance matrix:
```
n <- ncol(df)
cov <- matrix(nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:i){
    x <- mean(c(ploidy[i],ploidy[j]))
    cov[i,j] <- mean((df[,i]-x*p)*(df[,j]-x*p)/(x*p*(1-p)),na.rm=T)
    cov[j,i] <- cov[i,j]
  }	
}
```
Do a PCA on the matrix:
```
pc <- prcomp(cov,scale=T)
xlab <- paste0("PC1 (",round(summary(pc)$importance[2]*100),"%)")
ylab <- paste0("PC2 (",round(summary(pc)$importance[5]*100),"%)")
pcs <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],id=colnames(df),ploidy=ploidy)
```
Plot:
```
ggplot(pcs, aes(PC1, PC2, color=as.factor(ploidy))) +
  geom_point(size=7) +
  labs(x=xlab, y=ylab, color="Ploidy") +
  geom_text_repel(aes(label=id), size=4, force=20, color="black", max.overlaps = 100) +
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(color="black", linewidth=0.5),
        axis.text=element_text(size=11,color="black"),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks=element_line(color="black", linewidth=0.5),
        axis.title=element_text(size=12, color="black"),
        plot.title=element_text(size=14, color="black", hjust = 0.5),
        legend.text=element_text(size=11, color="black"),
        legend.title=element_text(size=12, color="black"),
        legend.key=element_blank(),
        aspect.ratio=1)
```
![Individuals PCA](Figures/individuals_pca.png)

### Filter out impure individuals

<a name="further_gatk"></a>

From observation alone it is clear to see that some individuals are plotting incorrectly. Individuals such as KAG.03tl plot close to the arenosa end of the scale (left) which is incorrect as we know from Marburger et al., it's pure lyrata. Impurity can occur due to many reasons, including: incorrect labelling of samples, sample mix ups in the lab, contamination etc. The below figure shows the other impure individuals, highlighted in yellow:

![Annotated Individuals PCA](Figures/annotated_individuals_pca.png)

Returning to gatk, filter out the impure individuals:
```
grep -o -E 'BZD-....|PEK-....|SCT-....|TEM-....|GYE-....|JOH-....|KAG-(01|02|04|05|06|07|08)..|LIC-....|LOI-....|MAU-....|MOD-....|PIL-....|SCB-....|SWB-....|HAB-....|ROK-....|FRE-(01|02|03|04|05|07)tl|OCH-(01|02|03|04|06|07|08)tl|KEH-(01|02|03|04|05)tl' Chrom_1_noSnakemake.lyrata.bipassed.dp.m.bt.1pct.ld_pruned.vcf > filtered_tetraploid.args
```
Filter the vcf:
```
gatk SelectVariants -R [name_of_your_reference_file].fasta -V [name_of_your_vcf].vcf -sn filtered_tetraploid.args -O filtered_tetraploid.vcf 
```
Download the vcf and then read it into your local pc in R:
```
vcf <- read.vcfR("filtered_tetraploids.vcf")
```
Run PCA:
```
df <- extract.gt(vcf)
df[df == "0|0"] <- 0
df[df == "0|1"] <- 1
df[df == "1|0"] <- 1
df[df == "1|1"] <- 2
df[df == "0/0"] <- 0
df[df == "0/1"] <- 1
df[df == "1/0"] <- 1
df[df == "1/1"] <- 2
df[df == "0/0/0/0"] <- 0
df[df == "0/0/0/1"] <- 1
df[df == "0/0/1/1"] <- 2
df[df == "0/1/1/1"] <- 3
df[df == "1/1/1/1"] <- 4
df[df == "0/0/0/0/0/0"] <- 0
df[df == "0/0/0/0/0/1"] <- 1
df[df == "0/0/0/0/1/1"] <- 2
df[df == "0/0/0/1/1/1"] <- 3
df[df == "0/0/1/1/1/1"] <- 4
df[df == "0/1/1/1/1/1"] <- 5
df[df == "1/1/1/1/1/1"] <- 6
df <- data.frame(apply(df,2,function(x)as.numeric(as.character(x))))
```
Remove samples with > 50% missing data:
```
mis <- apply(df,2,function(x)sum(is.na(x))/length(x))
df <- df[,mis <= 0.5]
```
Calculate allele frequencies:
```
ploidy <- apply(df,2,max,na.rm=T)
p <- apply(df,1,function(x)sum(x,na.rm=T)/sum(ploidy[!is.na(x)]))
```
Removing individuals can change allele frequencies, so we make sure that maf >= 0.05:
```
df <- df[p >= 0.05 & p <= 0.95,]
p <- p[p >= 0.05 & p <= 0.95]
```
Estimate a covariance matrix:
```
n <- ncol(df)
cov <- matrix(nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:i){
    x <- mean(c(ploidy[i],ploidy[j]))
    cov[i,j] <- mean((df[,i]-x*p)*(df[,j]-x*p)/(x*p*(1-p)),na.rm=T)
    cov[j,i] <- cov[i,j]
  }	
}
```
Do PCA on the matrix:
```
pc <- prcomp(cov,scale=T)
xlab <- paste0("PC1 (",round(summary(pc)$importance[2]*100),"%)")
ylab <- paste0("PC2 (",round(summary(pc)$importance[5]*100),"%)")
pcs <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],id=colnames(df),ploidy=ploidy)
```
Plot:
```
ggplot(pcs, aes(PC1, PC2, color=as.factor(ploidy))) +
  geom_point(size=7) +
  labs(x=xlab, y=ylab, color="Ploidy") +
  geom_text_repel(aes(label=id), size=4, force=20, color="black", max.overlaps = 100) +
  theme(panel.background = element_blank(),
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(),
        axis.line=element_line(color="black", linewidth=0.5),
        axis.text=element_text(size=11,color="black"),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks=element_line(color="black", linewidth=0.5),
        axis.title=element_text(size=12, color="black"),
        plot.title=element_text(size=14, color="black", hjust = 0.5),
        legend.text=element_text(size=11, color="black"),
        legend.title=element_text(size=12, color="black"),
        legend.key=element_blank(),
        aspect.ratio=1)
```
![Filtered PCA](Figures/filtered_pca.png)

It is obvious by now that PC1 splits individuals closest to arenosa (left), hybrids (central) and pure lyrata (right), however it is slightly unclear what PC2 annotates. For example, why does BZD cluster on its own? Investigating relatedness may clear this up.

## Relatedness Calculations

<a name="relatedness_calculations"></a>

### Calculate Nei's distances

<a name="nei"></a>

Create a matrix of pairwise distances for **individuals**:
```
aa.D.ind <- stamppNeisD(aa.genlight, pop = FALSE)
```
Create a matrix of pairwise distances for **populations**:
```
aa.D.pop <- stamppNeisD(aa.genlight, pop = TRUE)
```
Create the distance objects:
```
colnames(aa.D.ind) <- rownames(aa.D.ind)
aa.D.ind.dist <-as.dist(aa.D.ind, diag=T)
attr(aa.D.ind.dist, "Labels") <-rownames(aa.D.ind)

colnames(aa.D.pop) <- rownames(aa.D.pop) 
aa.D.pop.dist <-as.dist(aa.D.pop, diag=T)
attr(aa.D.pop.dist, "Labels") <-rownames(aa.D.pop)
```
Plot and save the Neighbour Joining (NJ) tree:
```
plot(nj(aa.D.ind), typ="unrooted", cex=0.7)
title(expression("Neighbour-joining tree of distance-based analysis of "*italic(Arabidposis)*" "))
write.tree(nj(aa.D.pop),file="NJ.distance_tree_outgroups.tre")
```
This code plots the data as individuals and saves a .tre file in their populations for plotting in SplitsTree later. You can alter this in the code just above where it says 'plot(nj(aa.**ind**)...' and 'write.tree(nj(aa.D.**pop**)' change 'ind' and 'pop' accordingly. It is useful to plot as individuals to see if there are still impure data or any occurences which are different to what is expected.

The following plot is produced for the data I am using:
![NJ](Figures/NJ.png)

This tree is unrooted and follows the same trends as seen with the PCA and purity figure: LIC and MOD branch together as they are the purest lyrata, and share a recent common ancester. JOH (one of the new populations) branches with these two and clusters with them on the PCA meaning it must be a pure lyrata. The SCT, PEK and TEM populations, which we first saw clustering seperately on the PCA and then discovered a geographical similarity, branch as the same trio, showing a genetic liklihood. Other trends show and are easier to see when uploading to the SplitsTree software.

### SplitsTree

<a name="splitstree"></a>

After downloading SplitsTree (instructions at the top of the page) ipload your .tre file to the SplitsTree software by navigating to [File] then [Open]. Once loaded select [Tree] then [NJ] to view your tree as a midpoint rooted neighbour joined tree. The output should look as follows:

![SplitsTree](Figures/SplitsTree.png)

Interestingly, SWB and MAU are outgrouped. From the Marburger et al., purity plot these individuals were annotated as approximate pure lyrata. These two populations likely diverged first after undergoing WGD and havent hybridised. Interestingly the population closest to arenosa (KEH) is nested amongst the lineage, this supports the theory that these lyrata tetraploids hybridised post WGD. The distance from KEH from its shared common ancestor with OCH is extensive, meaning the further aquisition of arenosa material occured after divergence. The allopolyploids form a monophyletic group (KEH to PIL) every other population external to this (except GYE) is relatively pure lyrata.

## FastStructure

<a name="fast_structure"></a>

### Prepare your environment:

make a directory for the omicsspeaks csv output
```
mkdir final_omicsspeaks_output
```
make a directory for population specific .str files we produce when reordering the faststructure input
```
mkdir individual_str_output
```
make a directory for all filtered vcf output data
```
mkdir final_populations
```
make directory to store all the output for the individuals
```
mkdir individual_population_files
```
make directories to store faststructure output
```
mkdir -p faststructure_output/vcf_dir
mkdir faststructure_output/final_svg_files
```
get the 'attributes' of your environment 
```
source $HOME/.bash_profile
```
### Getting a filtered vcf. The populations should be given in the order you want the final faststructure output to be in. 
The retrieve_IDs_updated_FIX.py script is required here, along with an .args file of what populations you want to exclude (same from earlier when removing impure populations for the PCA), the file I used can be found above and is called 'samples_to_exclude.args'
activate conda env
```
conda activate /shared/apps/conda/bio2
```
get sample names to include in the filtered vcf, change the --pop flag accordingly
```
python3 retrieve_IDs_updated_FIX.py \
	-i [name_of_your_vcf].vcf \
	--pop 'KEH','BZD','OCH','FRE','ROK','HAB','KAG','MAU','JOH','MOD','LIC','PEK' \
	--same_file 'yes' \
	-odir final_populations \
	-opre diploids_hybird_1arenosa_1lyrata \
	--fast_struc 'yes' \
	-xcl samples_to_exclude.args 
```
produce a dictionary file for your reference fasta (if you already followed these steps earlier under [Filter Data for Further Analysis of Trends / Filter out impure individuals] skip the next two lines.
```
gatk CreateSequenceDictionary \
	-R [name_of_your_reference_file].fasta
```
create index file for the fasta file
```
samtools faidx [name_of_your_reference_file].fasta
```
index the vcf file
```
gatk IndexFeatureFile \
	-I [name_of_your_vcf].vcf
```
filter the vcf to only contain biallelic variants we are interested in
```
gatk SelectVariants \
	-R [name_of_your_reference_file].fasta \
	-V [name_of_your_vcf].vcf \
	-sn final_populations/gatk_args_output/*.args \
	--output final_populations/diploids_hybird_1arenosa_1lyrata.vcf
```
deactivate conda env
```
conda deactivate
```
### Producing individual vcf files for each of the populations you have narrowed down to 
activate conda environment
```
conda activate /shared/apps/conda/bio2
```
Get sample names to include in the filtered vcf, i.e. do everything again but produce files for the individual populations
```
python3 retrieve_IDs_updated_FIX.py \
        -i [name_of_your_vcf].vcf \
        --pop 'KEH','BZD','OCH','FRE','ROK','HAB','KAG','MAU','JOH','MOD','LIC','PEK' \
        --same_file 'no' \
        -odir individual_population_files \
	-xcl samples_to_exclude.args 
```
go through all the arg files to produce a vcf for each individual population
```
for file in individual_population_files/gatk_args_output/*.args ; do
	
	# get name of just file and not the full path
	individual_file=$(basename "$file")
	
	## filter the vcf to only contain biallelic variants we are interested in
	gatk SelectVariants \
		-R [name_of_your_reference_file].fasta \
		-V [name_of_your_vcf].vcf \
		-sn "$file" \
		--output individual_population_files/"$individual_file.vcf"
done
```
deactivate conda environment
```
conda deactivate
```
### Convert filtered vcf file to a str file
The Cochlearia_create_structure_file.py script is required here.

move file with filtered vcf to the vcf directory
```
cp final_populations/*.vcf faststructure_output/vcf_dir
```
#convert polyploids data to format acceptable to fastSTRUCTURE
```
python3 Cochlearia_create_structure_file.py \
        -v faststructure_output/vcf_dir/ \
        -o filtered_vcf_converted_to_str \
        -s true
```
remove first and last line from .str file
```
sed -i '1d;$d' faststructure_output/vcf_dir/vcf_to_str/filtered_vcf_converted_to_str.StructureInputDiploidized.str
```
### Reorder the .str file based on what order we've done
The reorder_str_file.py, structure.py and chooseK.py scripts are required here.
```
python3 reorder_str_file.py \
	-i faststructure_output/vcf_dir/vcf_to_str/filtered_vcf_converted_to_str.StructureInputDiploidized.str \
	-o faststructure_output/vcf_dir/vcf_to_str/filtered_vcf_converted_to_str.StructureInputDiploidized2.str \
	-p 'KEH','BZD','OCH','FRE','ROK','HAB','KAG','MAU','JOH','MOD','LIC','PEK' \
	--str_dir individual_str_output
```
activate your conda environment for faststructure
```
conda activate /shared/conda/faststructure
```
move into final svg directory so all output gets put into there
```
cd faststructure_output/final_svg_files
```
run fas structure command on our reordered .str file with different K values
```
for i in {1..10}; do 
	echo $i
	python /shared/conda/faststructure/bin/structure.py \
		-K $i \
		--input ~/fast_structure/faststructure_output/vcf_dir/vcf_to_str/filtered_vcf_converted_to_str.StructureInputDiploidized2 \
		--output ~/fast_structure/faststructure_output/final_svg_files/filtered_populations \
		--format str \
		--full
done
```
run choose K
```
python /shared/conda/faststructure/bin/chooseK.py \
	--input ~/faststructure_output/final_svg_files/filtered_populations
```
deactivate conda environment
```
conda deactivate
```
### Run code to create csv file necessary for ommics speaks. 
for K3
```
paste -d '\t' ~/final_populations/faststructure_files/*faststructure_popnames.txt ~/faststructure_output/final_svg_files/filtered_populations.3.meanQ > \
~/final_omicsspeaks_output/omics_speaks_K3.tsv

cat ~/final_omicsspeaks_output/omics_speaks_K3.tsv | tr '\t' ',' | tr -s '[:blank:]' ',' > ~/final_omicsspeaks_output/omics_speaks_K3.csv
```
for K2
```
paste -d '\t' ~/final_populations/faststructure_files/*faststructure_popnames.txt ~/faststructure_output/final_svg_files/filtered_populations.2.meanQ > \
~/final_omicsspeaks_output/omics_speaks_K2.tsv

cat ~/final_omicsspeaks_output/omics_speaks_K2.tsv | tr '\t' ',' | tr -s '[:blank:]' ',' > ~/final_omicsspeaks_output/omics_speaks_K2.csv
```
for K4
```
paste -d '\t' ~/final_populations/faststructure_files/*faststructure_popnames.txt ~/faststructure_output/final_svg_files/filtered_populations.4.meanQ > \
~/final_omicsspeaks_output/omics_speaks_K4.tsv

cat ~/final_omicsspeaks_output/omics_speaks_K4.tsv | tr '\t' ',' | tr -s '[:blank:]' ',' > ~/final_omicsspeaks_output/omics_speaks_K4.csv
```
### Structure Plot
Upload your .csv file to Strucutre Plot, which can be access on the web via: http://omicsspeaks.com/strplot2/

Adjust which k value you are using before clicking submit

## Allele Frequency Spectrum

<a name="allele_frequency_spectrum"></a>

Another way of investigating whether a population is **allo-** or **auto-** polyploid is by creating Allele Frequency Spectrums (AFS). Most commonly this is represented on a histogram. The output visualises the genetic variation within a population by describing the distribution of allele frequencies.

We know from previous research by Yant et al. that allo- and auto- **hexaploids** plot on a histogram like such:
![Example Allo/AutoPolyploids](Figures/auto:allo-ploids.png)
The **autohexaploids** plot with an exponential distribution which is right skewed, this is because the least frequent allele frequencies are the rarely occuring single nucleotide polymorphisms (SNPs). There are many different SNPs which occur close to a frequency of zero, causing this high left peak. **Allohexaploids** plot with a left, central and right peak. Again the data contains a high count of in-frequent SNPs which form the left peak. The central peak represents the alleles which occcur at virtually half frequency becuase on average hybrids contain 50/50 alleles from the two populations that merged and thus the frequently occuring alleles in one population are only half in the new population. And the right sided peak is formed from alleles which are conserved over both populations.

Autopolyploids always follow this exponential distribution, where the increase in sets of chromosomes only broadens the allele frequency distribtution. However, it is uncertain how a allo**tetraploid** would plot as the above histogram plots **hexaploid** data. To avoid potential overlooking of data and to draw reasonable comparisons, it is best to create an example allotetraploid to compare the sample data with.

### Creating an example allotetraploid:

<a name="allotetraploid"></a>

Download the files 'arenosa_632.txt' and 'lyrata_272_with_some_hybrids.txt' to the working directory you set earlier. These files contain mostly tetraploid data of arenosa and lyrata samples which we can merge to create an allotetraploid. If you are using different populations locate files which contain atleast allele frequencies at specific positions on chromosomes.

Read the data in:
```
arenosa_632 <- read.table(file ='arenosa_632.txt' ,header = TRUE,sep = '\t')
lyrata_272 <- read.table(file ='lyrata_272_with_some_hybrids.txt' ,header = TRUE,sep = '\t')
```
Merge the two data frames based on the 'POS' column:
```
merged_df <- merge(arenosa_632, lyrata_272, by = "POS", suffixes = c("_arenosa", "_lyrata"))
```
Calculate the mean of the 'AF' column for matching rows:
```
merged_df$Mean_AF <- (merged_df$AF_arenosa + merged_df$AF_lyrata) / 2
```
Create a new data frame with 'POS' and 'Mean_AF':
```
new_df <- merged_df[, c("POS", "Mean_AF")]
```
This new data frame has the mean allele frequencies at the same sites for lyrata and arenosa.

Plot:
```
ggplot(data = new_df, aes(Mean_AF)) +
  geom_histogram(color='black',fill='white', bins = 100)
```
This plot is shown below as 'Un-Filtered' because we had to filter the data as the uncommon SNPs conceal potential trends in the data (you can see a slight peak in the middle however its masked because theres such a high count of in-frequent alleles), so filter the mean allele frequency to be greater than 0.1:
```
filtered_df <- new_df[new_df$Mean_AF > 0.1, ]
```
Now plot:
```
ggplot(data = filtered_df, aes(Mean_AF)) +
  geom_histogram(color='black',fill='white', bins = 100)
```
![Synthetic Plots](Figures/synthetic.png)

From the filtered plot it is clear to see there is still a central peak in allotetraploids it it just less defined than the allohexaploid plot (note this may be becuase of different sample or bin sizing). Now lets compare our data to this allotetraploid.

### Creating allele frequency histograms

<a name="histogram"></a>

Firstly, in a Unix environment:
This requires the 'pops.txt' file created earlier for fastStructure and the poly_freq.c script. 

Compile the scipt into an executable environment called poly_freq:
```
gcc poly_freq.c -o poly_freq -lm
```
Execute the script:
```
 ./poly_freq -vcf [title_of_your_vcf].vcf -pops pops.txt > info.tsv
```
This will create a file with population-specific allele frequencies

Now read the tsv file into R:
```
df <- read.table(file ='info.tsv', header = TRUE, sep = '\t')
```
Store allele frequencies of the population you want to plot in a variable
```
allele_frequencies <- df$BZD
```
#### plot a histogram of the results
```
ggplot(data=df, aes(x=allele_frequencies)) + geom_histogram(color='black',fill='white')
```
![BZD Histogram](Figures/BZD_histogram.png)

## Creating Selection Scans


