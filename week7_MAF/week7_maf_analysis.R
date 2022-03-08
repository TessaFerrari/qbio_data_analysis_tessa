#--------------
# Exercise 1.1
#--------------

#BiocManager::install("maftools")
library(maftools)
library(TCGAbiolinks)

#--------------
# Exercise 1.2
#--------------

# SET YOUR WORKING DIRECTORY TO YOUR analysis_data FOLDER
getwd()
setwd("../analysis_data")  # fill in here

# replace the path!
# fread is a way to read in a data frame, but it's much faster than the built-in read.csv()
# by default, it reads it in as a slightly different data type than a data frame, so we set the data.table flag to false
clinic <- data.table::fread("../week4_clinical/coad_clinical_data.csv",
                            data.table = F)

# rename the patient barcode to make it work with maftools
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

#--------------
# Exercise 1.3
#--------------

length(colnames(clinic))                    #1) 78
colnames(clinic) == "bcr_patient_barcode"   #2) 78 TRUE/FALSE if the colname is "bcr_patient_barcode"
                                            #3) None

#--------------
# Exercise 1.4
#--------------

# fill in!!
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

#--------------
# Exercise 1.5
#--------------

# 1. clear your environment either using the following code or by clicking the broom.

# 2. Locate the appropriate csv file
getwd()
setwd("GDCdata")
list.files()

# 3. Use fread to read in the *mutation* data you just downloaded (refer to the fread syntax from the previous example)
maf_dataframe = data.table::fread("TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.csv",
                                  data.table = F)

# 4. Read in your clinic data frame
# you also need to rename the bcr_patient_barcode again
# refer to the previous example!
clinic <- data.table::fread("../../week4_clinical/coad_clinical_data.csv",
                            data.table = F)
colnames(clinic)[ colnames(clinic) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"

# re-create the maf_object as before
mutation_query <- GDCquery_Maf(tumor = "COAD", 
                               pipeline = "mutect2",
                               save.csv = TRUE)

maf_object <- read.maf(maf = mutation_query, 
                       clinicalData = clinic, 
                       isTCGA = TRUE)

#--------------
# Exercise 2.1
#--------------

# explore here!
maf_object
str(maf_object)
maf_object@data
maf_object@clinical.data

#--------------
# Exercise 3.1
#--------------

# play with the number of genes with the "top" argument
oncoplot(maf = maf_object,
         top = 10)
library(tidyverse)
ggsave("../week7_maf/oncoplot.png")

#--------------
# Exercise 3.2
#--------------

# APC gene
# helps control how often a cell divides, how it attaches to other cells within a tissue, how the cell polarizes and the morphogenesis of the 3D structures

#--------------
# Exercise 3.3
#--------------

# 1. Write to clinic again
clinic = maf_object@clinical.data

# 2. Create the young_patients_ids vector
young_patients_ids = clinic$Tumor_Sample_Barcode[clinic$age_category=="young"]

# 3. Create another 
young_maf = subsetMaf(maf = maf_object,
                      tsb = young_patients_ids)

# 4. Repeat steps 2-3 to create an old_maf! Can you do it in one line?
old_patients_ids = clinic$Tumor_Sample_Barcode[clinic$age_category=="old"]
old_maf = subsetMaf(maf = maf_object,
                      tsb = old_patients_ids)


#--------------
# Exercise 3.4
#--------------

coOncoplot(m1 = young_maf, 
           m2 = old_maf, 
           m1Name = "Young", 
           m2Name = "Old")

ggsave("../week7_maf/coOncoplot.png")

#--------------
# Exercise 3.5
#--------------

# APC, TTN, and TP53 are more mutated in old patients
# KRAS is mutated equally in young and old patients
# MUC16 is more mutated in young patients

#--------------
# Exercise 3.6
#--------------

# pick a gene to look at!
lollipopPlot(maf_object, gene = "APC")

ggsave("../week7_maf/lollipopPlot.png")

#--------------
# Exercise 3.7
#--------------

# pick a gene to look at!
lollipopPlot2(m1 = young_maf, 
              m2 = old_maf, 
              m1_name = "Young",
              m2_name = "Old",
              gene = "APC")

ggsave("../week7_maf/lollipopPlot2.png")

# More commonly mutated in old patients
# Mutations see concentrated in the APC_crr region
# Nonsense mutations are the most common 
# Certain types of mutations seem to be more common in certain areas

#--------------
# Exercise 3.8
#--------------

# Mutation A is present in 3/5ths of the samples
# Mutation B is present in 9/20ths of the samples

#--------------
# Exercise 3.9
#--------------

b = 7
c = 2
d = 35
e = 37
f = 42

#---------------
# Exercise 3.10
#---------------

# geneA is TP53
geneA_maf <- subsetMaf(maf = maf_object,
                       genes = "TP53")

# geneB is KRAS
geneB_maf <- subsetMaf(maf = maf_object, 
                       genes = "KRAS")

#---------------
# Exercise 3.11
#---------------

# subsetMaf() is looking only at the data for the gene specified (nGenes is now 1)
# Each sample can have more than 1 mutation in geneA
# There are more samples in the data section because some patients have more than one mutation

#---------------
# Exercise 3.12
#---------------

# 1. Access the barcodes of the patients with mutations in genes A and B
# bc stands for barcode
mut_bc_geneA = geneA_maf@clinical.data$Tumor_Sample_Barcode
mut_bc_geneB = geneB_maf@clinical.data$Tumor_Sample_Barcode

# 2. Get the lengths of these two vectors
num_mut_geneA = length(mut_bc_geneA)
num_mut_geneB = length(mut_bc_geneB)
# Numbers correspond to the sum mutation in geneA and geneB portion

# 3. Fill in the intersect here! Then get the nubmer of patients
mut_bc_geneAB = intersect(mut_bc_geneA, mut_bc_geneB)
num_mut_geneAB = length(mut_bc_geneAB)

#---------------
# Exercise 3.13
#---------------

num_mut_geneA_only = num_mut_geneA - num_mut_geneAB
num_mut_geneB_only = num_mut_geneB - num_mut_geneAB

#---------------
# Exercise 3.14
#---------------

num_mut_neither = length(maf_object@clinical.data$Tumor_Sample_Barcode) - num_mut_geneB - num_mut_geneA_only

contig_table <- matrix(c(num_mut_geneAB, 
                         num_mut_geneB_only,
                         num_mut_geneA_only,
                         num_mut_neither), 
                       nrow=2)

# view the contingency table
contig_table

#---------------
# Exercise 3.15
#---------------

fe_results <- fisher.test(contig_table)
fe_results
# The p-value is 0.06543, which is greater than 0.05 (i.e. not statistically significant)