# Tessa Ferrari

#-----------
# Libraries
#-----------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)

#-------------
# RNAseq Data
#-------------

getwd()
setwd("../analysis_data")
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
sum_exp <- GDCprepare(query)

#---------------
# Prepping data
#---------------

# Mask created using boolean indexing in order to get rid of rows with NA CEA levels
NAmask = !is.na(colData(sum_exp)$paper_preoperative_pretreatment_cea_level)&colData(sum_exp)$paper_preoperative_pretreatment_cea_level!="NA"
clinicalData = colData(sum_exp)[NAmask, ]     # Data is saved as clinicalData so as to preserve the original
# Converting CEA level column from a vector to a numeric
clinicalData$paper_preoperative_pretreatment_cea_level = as.numeric(levels(clinicalData$paper_preoperative_pretreatment_cea_level))[clinicalData$paper_preoperative_pretreatment_cea_level] # converts CEAdata from factor to numeric

# Column CEAlevel categorizes patients as either having normal CEA levels or high CEA levels
clinicalData$CEAlevel = ifelse(clinicalData$paper_preoperative_pretreatment_cea_level<=3, "Normal", "High")

# The ensemble_gene_id for FUT2 was found by masking the gene id column with the gene name
FUT2geneID = rowData(sum_exp)$ensembl_gene_id[rowData(sum_exp)$external_gene_name=="FUT2"]
# A vector containing the trascript counts for FUT2 (not including rows with NA CEA levels)
FUT2data = assays(sum_exp)$"HTSeq - Counts"[rownames(assays(sum_exp)$"HTSeq - Counts")==FUT2geneID, NAmask]

#--------------------------------------------
# Scatterplot - CEA Level v. FUT2 Expression
#--------------------------------------------

plot(x = FUT2data[clinicalData$paper_preoperative_pretreatment_cea_level<250], # <250 removes extreme CEA level outliers to make the plot easier to read
     y = clinicalData$paper_preoperative_pretreatment_cea_level[clinicalData$paper_preoperative_pretreatment_cea_level<250],
     main = "CEA Level v. FUT2 RNAseq",
     xlab = "FUT2 RNAseq Counts",
     ylab = "CEA Level (ng/mL)",
     pch = 20)

#----------------------------------------
# Boxplot - CEA Level v. FUT2 Expression
#----------------------------------------

boxplot(FUT2data[clinicalData$paper_preoperative_pretreatment_cea_level<250]~clinicalData$CEAlevel[clinicalData$paper_preoperative_pretreatment_cea_level<250],
        main = "CEA Level v. FUT2 RNAseq",
        xlab = "CEA Level (High: >3ng/mL, Normal: <=3ng/mL)",
        ylab = "FUT2 RNAseq Counts",
        notch = TRUE,
        las = 2, 
        cex.axis = 0.5)

#-----------------------------
# Histogram - FUT2 Expression
#-----------------------------

hist(FUT2data,
     xlab = "FUT2 RNAseq Counts",
     main = "Distribution of FUT2 Expression")

#-----------------------
# Histogram - CEA Level
#-----------------------

hist(clinicalData$paper_preoperative_pretreatment_cea_level[clinicalData$paper_preoperative_pretreatment_cea_level<250],
     xlab = "CEA Level (ng/mL)",
     main = "Distribution of CEA Level")

#-----------------------------------------------------
# Kaplan-Meier plot - Survivorship based on CEA Level
#-----------------------------------------------------

# Replaces NA values in days_to_death with the values from days_to_last_follow_up
isDeathNA = is.na(clinicalData$days_to_death)
clinicalData$days_to_death = ifelse(isDeathNA, clinicalData$days_to_last_follow_up, clinicalData$days_to_death)

# Creates a boolean vector with the TRUEs being dead
clinicalData$death_event = ifelse(clinicalData$vital_status=="Dead", 1, 0)

# Initialize the survival object
surv_object <- Surv(time = clinicalData$days_to_death, 
                    event = clinicalData$death_event)
# Create the fit object
CEA_fit <- surv_fit( surv_object ~ clinicalData$CEAlevel, data = clinicalData )
# Formatting
survplot = ggsurvplot(CEA_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
CEA_survplot = survplot$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
CEA_survplot

# Save the plot as a png
ggsave("../midsemester_project_ferrari/kmplot_by_CEAlevel.png", plot = CEA_survplot, width = 12, height = 9)
# Save the clinical data as a csv
write.csv(clinicalData, "../midsemester_project_ferrari/coad_clinical_data.csv", row.names = F)