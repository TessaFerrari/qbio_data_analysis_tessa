# this will install packages (if necessary) and load them
if(!require(BiocManager)) install.packages("BiocManager")

# the double colon syntax calls a function from a specific package
# this avoids loading the entire package
# in this case, we need to download TCGAbiolinks from Bioconductor using BiocManager
if(!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

# this just loads a package
library(TCGAbiolinks)

setwd("GitHub/qbio_data_analysis_tessa/analysis_data")

clin_query <- GDCquery(project = "TCGA-COAD", data.category = ..., file.type = "xml")
# Only use this line ONCE! Comment out after you have downloaded the data. 
#GDCdownload(clin_query)
clinic <- GDCprepare_clinic(clin_query, clinical.info = "patient")
# Just adding an underscore between follow and up
names(clinic)[names(clinic)=="days_to_last_followup"] <- "days_to_last_follow_up"

# Exercise 1.1
str(clinic)
head(clinic)
dim(clinic)

# Exercise 1.2
colnames(clinic)
clinic$gender

# Exercise 2.1
plot(x = clinic$age_at_initial_pathologic_diagnosis, y = clinic$weight, 
     main = "Weight vs. Age",
     xlab = "Age (Years)",
     ylab = "Weight (kg)")

# Exercise 2.2
unique(clinic$race_list) # This shows you all the UNIQUE entries of the race column
#the mar argument in par() sets the plot margins. 
# As the race names may be long, we want to have a large bottom margin. 
# The syntax is par(mar = c(bottom, left, top, right)), where mar is margins
# How does the plot change if you change these margins?
> par(mar=c(10,1,1,1))
# there's a few extra parameters to make the plot display better
boxplot(clinic$age_at_initial_pathologic_diagnosis~clinic$race_list,  # fill in variables to plot here
        las = 2, 
        cex.axis = 0.5)

# Exercise 2.3
levels(clinic$race_list) = c("No data",
                             "AMERICAN INDIAN OR ALASKA NATIVE",
                             "ASIAN",
                             "BLACK OR AFRICAN AMERICAN",
                             "WHITE")

# Exercise 2.4
min(clinic$age_at_initial_pathologic_diagnosis)
max(clinic$age_at_initial_pathologic_diagnosis)
mean(clinic$age_at_initial_pathologic_diagnosis)
median(clinic$age_at_initial_pathologic_diagnosis)
summary(clinic$age_at_initial_pathologic_diagnosis)

# Exercise 2.5
old = function(x){
  if(x>=50){return(TRUE)}else{return(FALSE)}
}
isOld = mapply(old, clinic$age_at_initial_pathologic_diagnosis)
numOld = sum(isOld)
numYoung = length(isOld)-sum(isOld)

# Exercise 2.6
isYoung = !isOld
young_patient_ids = clinic$bcr_patient_barcode[isYoung]
# This is a good sanity check (making sure you have the same number of patients)
# length(young_patient_ids) == num_young #Sanity checks are always important!
old_patient_ids = clinic$bcr_patient_barcode[isOld]

# Exercise 2.7
# create the new column
clinic$age_category = ifelse(isOld, "old", "young")
# Do some sanity checks!
# look at the first few rows -- how do you do this?
head(clinic$age_at_initial_pathologic_diagnosis)
head(clinic$age_category)
# count again to make sure everything looks right
# refer to the previous exercise!

# Exercise 2.8
clinic[1,1] #This is the top left entry of the dataframe. R has "one-based indexing"
clinic[1,]  #1st row
clinic[2:5,]#2nd to 5th rows
clinic[,3]  #3rd column
# A blank before the comma means all rows will be included
# A blank after the comma means all columns will be included

# Exercise 2.9
young_clinic = clinic[isYoung,]
old_clinic = clinic[isOld,]

# Exercise 2.10
young_clinic_one_line = clinic[clinic$age_at_initial_pathologic_diagnosis < 50, ]
# Check if your dimensions are the same! dim() gives dimensions of a data frame in number_of_rows, number_of_cols
identical(dim(young_clinic), dim(young_clinic_one_line))

#install.packages("survival")
library(survival)
#install.packages("survminer")
library(survminer)

# Exercise 3.1
isDeathNA = is.na(clinic$days_to_death)
clinic$days_to_death = ifelse(isDeathNA, clinic$days_to_last_follow_up, clinic$days_to_death)

# Exercise 3.2
clinic$death_event = ifelse(clinic$vital_status=="Dead", 1, 0)
# We initialize a 'survival' object first, which contains the data we need.
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)
# We then create a fit object
race_fit <- surv_fit( surv_object ~ clinic$race_list, data = clinic )
#the ggtheme and legend arguments are for formatting. 
# Feel free to play around with the margins and legend placement
survplot = ggsurvplot(race_fit, 
                      pval=TRUE, 
                      ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                      legend = "right")
p = survplot$plot + 
  theme_bw() +  # changes the appearance to be a bit prettier
  theme(axis.title = element_text(size=20), # increase font sizes
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
p
# save the plot as a png
# if you want to present the data, change the strata labels!
ggsave("../week4_clinical/kmplot_by_race.png", plot = p, width = 12, height = 9)

# Exercise 3.3
# Key takeaways are that the survivorship for white patients is highest, and survivorship for patients with no race data collected is the lowest
# Questions: Why is the survivorship for asian patients 100%? Do we need more data from asian patients? Why do black patients have worse outcomes? Are they receiving a lower standard of care?
# The plot could be improved by changing the Time units from days to years

# Exercise 4.1
# change the file path! make sure it's in your week4 folder
# we set row.names to false since the rows don't contain any relevant info
write.csv(clinic, "/Users/tessa/OneDrive/Documents/GitHub/qbio_data_analysis_tessa/week4_clinical/coad_clinical_data.csv", row.names = F)
#READ IN FILE > clinic_read_in <- read.csv("file_path/.../file_name.csv")