#--------------
# Exercise 1.1
#--------------
getwd()
setwd("GitHub/qbio_data_analysis_tessa/analysis_data")
#BiocManager::install("SummarizedExperiment")
library(TCGAbiolinks)
library(SummarizedExperiment)

#--------------
# Exercise 2.1
#--------------
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling", # get the RNA-seq transcriptome
                  data.type = "Gene Expression Quantification", # gets the counts
                  workflow.type = "HTSeq - Counts") # gets the raw counts processed by this method
GDCdownload(query) # only need to download the data once! Comment this out once you have completed it once
sum_exp <- GDCprepare(query)
str(sum_exp)

#--------------
# Exercise 2.2
#--------------
counts = assays(sum_exp)$"HTSeq - Counts"[1:5, 1:5]"
#[1:5, 1:5] means rows 1-5 and columns 1-5
#rows are gene names and columns are patients

head(rowData(sum_exp))
#rows are gene IDs and columns are alternate gene names

colData(sum_exp)[1:5, 25:29]

metadata(sum_exp)

#--------------
# Exercise 2.3
#--------------
dim(colData(sum_exp))
dim(rowData(sum_exp))
dim(assays(sum_exp)$"HTSeq - Counts")
#rows in rowData and rows in assays are the same (genes)
#rows in colData and cols in assays are the same (patients)

#--------------
# Exercise 2.4
#--------------
str(colData(sum_exp))
head(colData(sum_exp))
#rows are patients columns are clinical data
#rows in colData correspond to the cols in assays

#--------------
# Exercise 2.5
#--------------
colnames(colData(sum_exp))
#age column is "age_at_diagnosis"

#--------------
# Exercise 2.6
#--------------
sum_exp$age_at_diagnosis[1:10]
#units are days

#--------------
# Exercise 2.7
#--------------
#colData(sum_exp)$age_at_diagnosis = colData(sum_exp)$age_at_diagnosis/365

#--------------
# Exercise 2.8
#--------------
colData(sum_exp)$age_category = ifelse(colData(sum_exp)$age_at_diagnosis>=50, "Old", "Young")

#--------------
# Exercise 2.9
#--------------
head(rowData(sum_exp))
dim(rowData(sum_exp))
#columns are gene ID, gene name, and original gene ID
#rows are genes

#---------------
# Exercise 2.10
#---------------
ex_list <- c("apple", "orange", "grape", "kiwi", "mango")
"apple" %in% ex_list

"KRAS" %in% rowData(sum_exp)$external_gene_name
"BRAF" %in% rowData(sum_exp)$external_gene_name

#---------------
# Exercise 2.11
#---------------
assays(sum_exp)$"HTSeq - Counts"[20:25, 30:35]
#rows are genes
#columns are patients

#---------------
# Exercise 2.12
#---------------
#STEP 1
geneA_id_mask = (rowData(sum_exp)$external_gene_name == "KRAS") #what should go before the $ in this line?
sum(geneA_id_mask) #BEFORE running this line, guess the result. Hint: How many TRUE's should be in this mask? How many times should your gene name appear?

#STEP 2
ensembl_geneA = rowData(sum_exp)$ensembl_gene_id[geneA_id_mask] #fill in the dotted lines. 

#STEP 3: Repeat for geneB
geneB_id_mask = (rowData(sum_exp)$external_gene_name == "BRAF") #what should go before the $ in this line?
sum(geneB_id_mask) #BEFORE running this line, guess the result. Hint: How many TRUE's should be in this mask? How many times should your gene name appear?

#STEP 4
ensembl_geneB = rowData(sum_exp)$ensembl_gene_id[geneB_id_mask] #fill in the dotted lines. 

#A boolean mask was created by finding the external gene names that were equal to the gene of interest
#then we summed the mask to make sure that only 1 gene is marked as true in the vector
#Finally we used the mask to get rid of all the gene IDs that don't correspond to the gene of interest

#---------------
# Exercise 2.13
#---------------
head(assays(sum_exp)$"HTSeq - Counts")
#ensembl gene IDs are the rows

#---------------
# Exercise 2.14
#---------------
min(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ])  #On which side of the comma does ensembl_geneA go?
max(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ])# find max of geneA

summary(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ]) # Use the same thing as for the min() function, but for your second gene.

#---------------
# Exercise 2.15
#---------------
plot(assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, ],
     assays(sum_exp)$"HTSeq - Counts"[ensembl_geneB, ],
     xlab = "KRAS", # remember to rename axes!
     ylab = "BRAF"
)

#---------------
# Exercise 2.16
#---------------
bool_age_na = is.na(colData(sum_exp)$age_category)
num_na = sum(bool_age_na)


#---------------
# Exercise 2.17
#---------------
age_cat_no_NAs = colData(sum_exp)$age_category[!bool_age_na]

#---------------
# Exercise 2.18
#---------------
#use the length() function here
length(age_cat_no_NAs)

#Notes about dim()
dim( colData(sum_exp) ) #gives number of rows then number of columns
dim( colData(sum_exp) )[1] #gives number of rows
dim( colData(sum_exp) )[2] #gives number of columns

#Use the == instead of = to return a TRUE or FALSE from your math equation.
x = 3
y = 4
z = 7
x + y == z
z - x == y

#create your math equation here. Hint: You can either add two of the values to equal the third, or subtract two values to equal the third. 
dim(colData(sum_exp))[1]-num_na==length(age_cat_no_NAs)

#---------------
# Exercise 2.19
#---------------
dim(assays(sum_exp)$"HTSeq - Counts")
length(assays(sum_exp)$"HTSeq - Counts"[1,])
#this is NOT equal to the number of patients in age_cat_no_NAs 
#because we didn't delete the NAs in the original data, only a copy.

#---------------
# Exercise 2.20
#---------------
identical(rownames(colData(sum_exp)), colnames(assays(sum_exp)$"HTSeq - Counts"))
# Fill in the dotted lines with "row" or "col"

gene_counts = assays(sum_exp)$"HTSeq - Counts"[ensembl_geneA, !bool_age_na]

#---------------
# Exercise 2.21
#---------------
length(age_cat_no_NAs)==length(gene_counts)
#we created gene_counts by using 2 boolean masks, 
#one to filter out the patients with NA in their age 
#category, and one to filter out only our gene of interest

#---------------
# Exercise 2.22
#---------------
boxplot(gene_counts ~ age_cat_no_NAs, 
        xlab = "Age Catagory", 
        ylab = "KRAS Gene Counts")

#--------------
# Exercise 3.1
#--------------
#1
#You access the HTSeq - Counts data frame via 
#assays(sum_exp)$"HTSeq - Counts".
#The rows are genes.
#The cols are patients.

#2
#You access the data frame via rowData(sum_exp).
#The rows of the counts data frame are the ensembl_gene_id
#for each of the genes. This corresponds to the rows on the 
#rowData data frame, and rowData stores alternate gene names and IDs.

#3
#You access the data frame via colData(sum_exp).
#The cols of the counts data frame are patient IDs.
#This corresponds to the rows on the colData data frame
#and colData stores additional clininal data about each patient.