#--------------
# Exercise 1.1
#--------------

# 1) Download and load in packages
getwd()
setwd("GitHub/qbio_data_analysis_tessa/analysis_data")
#BiocManager::install("DESeq2")
library(DESeq2)
library(TCGAbiolinks)
library(SummarizedExperiment)

# 2) Set our working directory
getwd()
setwd("C:/Users/tessa/OneDrive/Documents/GitHub/qbio_data_analysis_tessa/analysis_data")

# 3) Create our sum_exp object again
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
sum_exp <- GDCprepare(query)

#--------------
# Exercise 1.2
#--------------

# 1) ID patients with NA as their age
is.na(colData(sum_exp)$age_at_diagnosis)

# 2) Make a copy of clinical and counts data
patients_data = colData(sum_exp)  # contains the clinical data
counts = assays(sum_exp)$"HTSeq - Counts"
  
# 3) Convert counts into a dataframe (DOESN'T WORK... SKIP)

# 4) Remove patients with NA as their age
mask = !is.na(patients_data$age_at_diagnosis)
patients_data = patients_data[mask,]
counts = counts[,mask]

# 5) Create the age_category column in patient_data
patients_data$age_category = ifelse(patients_data$age_at_diagnosis>=18250, "Old", "Young")

# 6) Turn the age_category column into a factor
patients_data$age_category = factor(patients_data$age_category, levels = c("Young", "Old"))

#--------------
# Exercise 1.3
#--------------

# 1) Convert the rownames of counts to make them readable (DOESN'T WORK... SKIP)

# 2) Compute the sum across the row of the counts data frame
counts_row_sums = rowSums(counts)

# 3) Identify the genes that have fewer than 10 reads
low_counts_mask = ifelse(counts_row_sums>=10, TRUE, FALSE)
sum(!low_counts_mask)

# 4) Remove these lowly expressed genes from counts
counts = counts[low_counts_mask,]

#--------------
# Exercise 2.1
#--------------

# 1) Make DESeqDataSet
dds = DESeqDataSetFromMatrix(countData = counts,
                             colData = patients_data,
                             design = ~age_category)

# 2) Comparisons and statistical analysis
dds_obj = DESeq(dds)

# 3) Tells what comparisons get run
resultsNames(dds_obj)

# 4) Get the young vs. old comparison
results = results(dds_obj, format = "DataFrame", contrast = c("age_category", "Young", "Old"))

#--------------
# Exercise 2.2
#--------------

head(results)
my_df = data.frame(x = c('b', 'd', 'c', 'e', 'a'),
                   y = c(2,4,3,5,1))
order_indices = order(my_df$y)
# we expect c(5, 1, 3, 2, 4) because:
# 1 is index 5
# 2 is index 1
# 3 is index 3
# 4 is index 2
# 5 is index 4
order_indices  # note the order!
# note that we sort the rows
my_df = my_df[order_indices, ]
my_df

#--------------
# Exercise 2.3
#--------------

# 1) Create row_order, a vector sorting rows by padj
row_order = order(results$padj)

# 2) Use row_order to sort results
results = results[row_order, ]

# 3) Look at the 1st 20 rows of results using head()
head(results)

# 4) Pick a differentially expressed gene
  # a) More expressed in young or old patients?
  # Old
  # b) What is the gene's name/function?
  # ENSG00000202198 (7SK RNA)/Inhibits P-TEFb, which regulates transcription elongation

#--------------
# Exercise 2.4
#--------------

# 1) Create thresholds for Log2FoldChange and padj
log2FoldChange_threshold = 1
padj_threshold = 0.05

# 2) ID genes with a greater magnitude Log2FoldChange than the threshold
results_aboveLog2FoldChange = ifelse(abs(results$log2FoldChange)>log2FoldChange_threshold, TRUE, FALSE)

# 3) ID genes with a lower padj than the threshold
results_belowPadj = ifelse(results$padj<padj_threshold & !is.na(results$padj), TRUE, FALSE)

# 4) Subset the results that satisfy both criteria
criteria_mask = ifelse(results_aboveLog2FoldChange & results_belowPadj, TRUE, FALSE)
results_passCrit = results[criteria_mask,]

#--------------
# Exercise 2.5
#--------------

fc_threshold = 2  # set a threshold of at least a 2 fold increase (double)
p_threshold = 0.05  # set a threshold of adjusted p-value being <= 0.05

# fill in your plot code here!
# be sure to relabel the axes!
# note: you can perform the log transformation directly in the plot function
plot(x = results$log2FoldChange,
     y = results$padj,
     xlab = "Log2FoldChange (Young/Old)", # be sure the specify that it's young over old!
     ylab = "Adjusted P-Value",
     pch = 20) # smaller solid circles

# these lines put the lines on the plot
# abline() plots straight lines on an R plot.
# v argument is for a vertical line, h argument is for a horizontal line, col argument is color
abline(v=c(-log2(fc_threshold), log2(fc_threshold)), h= c(-log10(p_threshold)), col="green")

library(ggplot2)

volcano_plot = ggplot(data = data.frame(results), 
                      aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(color = ifelse(log2FoldChange < -1 & padj < 0.05, "lower in young",
                                ifelse(log2FoldChange > 1 & padj < 0.05, "higher in young", "NS"))),
             size = 0.5) + 
  theme_minimal() + # make things pretty +
  theme(legend.title = element_blank()) + 
  # next 2 lines draw lines at the thresholds
  geom_vline(xintercept=c(-log2(fc_threshold), log2(fc_threshold)), color="green") + 
  geom_hline(yintercept=-log10(p_threshold), color="green") + 
  scale_color_discrete(type=c("red", "blue", "black")) +
  labs(x = "log2 Fold Change (Young/Old)",
       y = "-log10 Adjusted p-value")


volcano_plot

#--------------
# Exercise 2.6
#--------------

# write CSV to your week6 folder
# We're not putting it into analysis_data since you actually generated the data
write.csv(x = results,
          file = "../week6_DESeq2/results.csv",
          row.names = FALSE)