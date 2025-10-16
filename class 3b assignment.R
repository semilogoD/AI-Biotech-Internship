# Assignment 
#-------------


#load Required Libraries
library(GEOquery)             # Download GEO datasets (series, matrix,raw CEL files)
library(affy)                 # pre-processing of Affymetrix microarray data (RMA normalization)
library(arrayQualityMetrics)  # QC reports for microarray data
library(limma)
library(AnnotationDbi)        #interface for annotation databases
library(hgu133plus2.db)       # anootation for Affymetrix Human Genome U133 plus 2.0
library(dplyr)                # data manipualtion

# 1. Perform quality control before and after normalization and 

unzip("Raw_Data/E-MTAB-9990_1.zip", exdir = "Raw_Data/E-MTAB-9990")


# Read E-MATB files

raw_data <- ReadAffy(celfile.path = "Raw_Data/E-MTAB-9990")

raw_data
# check whether any arrays are flagged as outliers. 

### QC before pre-processing ###
arrayQualityMetrics(expressionset = raw_data, 
                    outdir = "Results/QC_Raw_Data_E-MTAB-9990",
                    force = TRUE, do.logtransform = TRUE)

# note down how many you found before and after normalization
# none before normalization

# 2. Normalize the data and then apply filtering to remove low-intensity probes 
# and note how many transcripts remain.

# RMA Normalization #
normalized_data <- rma(raw_data)

# extract expression data
expression_data <- exprs(normalized_data)

# extract feature data 
feature_data <- fData(normalized_data)
nrow(feature_data)

#QC after normalization 
arrayQualityMetrics(expressionset = normalized_data, 
                    outdir = "Results/QC_E-MTAB_Normalized_Data",
                    force = TRUE)

# note down how many you found before and after normalization
#none after normalization

# extract normalized expression data

processed_data <- as.data.frame(exprs(normalized_data))
 

# 3. Use the phenotype information to define your target groups and re-label them (e.g normal vs cancer)

row_median <- rowMedians(as.matrix(processed_data))

row_median

#Plot the distribution of median intensities of probes
hist(row_median, 
     breaks = 100, 
     freq = FALSE,
     main = "Median Intensity Distribution")

threshold <- quantile(row_median, 0.25)

abline(v = threshold, col ="black", lwd = 2)

indx <- row_median > threshold

filtered_data <- processed_data[indx,]

#extracting phenotype data 

# Check if the zip file exists
zip_file <- unzip("Raw_Data/E-MTAB-9990.zip", exdir = "Raw_Data/E-MTAB-9990_SDRF_File")
file.exists(zip_file)


phenotype_data <- read.delim("Raw_Data/E-MTAB-9990_SDRF_File/E-MTAB-9990.sdrf.txt", 
                             stringsAsFactors = FALSE)

colnames(filtered_data) <- rownames(phenotype_data)

processed_data <- filtered_data

class(phenotype_data$Characteristics.sampling.site.)

groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("normal tissue adjacent to tumor", "tumor"),
                 label = c ("normal", "cancer"))


class(groups)
levels(groups) 

