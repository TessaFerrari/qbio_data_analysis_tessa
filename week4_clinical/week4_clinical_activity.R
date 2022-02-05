# reads in the clinic data
clinic <- read.csv("GitHub/qbio_data_analysis_tessa/week4_clinical/coad_clinical_data.csv")

# determines the column names of the clinic data
colnames(clinic)
# determines the number of NAs in the KRAS mutation column
sum(is.na(clinic$kras_mutation_found))
# determines the number of NAs in the number of relatives with cancer column
sum(is.na(clinic$number_of_first_degree_relatives_with_cancer_diagnosis))


#----------------------------------------
# 1) Relate your variables to each other
#----------------------------------------
clinic_hasKRAS = clinic[clinic$kras_mutation_found=="YES", ] #uses a mask to make clinic data of just patients w/ KRAS
# Histogram of the number of relatives w/ cancer of patients w/ KRAS
hist(clinic_hasKRAS$number_of_first_degree_relatives_with_cancer_diagnosis,
     xlab = "Number of Relatives w/ Cancer",
     main = "Patients w/ KRAS")

clinic_noKRAS = clinic[clinic$kras_mutation_found=="NO", ] #uses a mask to make clinic data of just patients w/o KRAS
# Histogram of the number of relatives w/ cancer of patients w/o KRAS
hist(clinic_noKRAS$number_of_first_degree_relatives_with_cancer_diagnosis,
     xlab = "Number of Relatives w/ Cancer",
     main = "Patients w/o KRAS")


#----------------------------------------------------------------
# 2) Relate your first variable to survival in colorectal cancer
#----------------------------------------------------------------
library(survival)
library(survminer)
# initializes a 'survival' object containing the clinic data
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)
# creates a fit object using the data from the KRAS column
KRAS_fit <- surv_fit( surv_object ~ clinic$kras_mutation_found, data = clinic )
# formatting for the KM plot
survplot = ggsurvplot(KRAS_fit, 
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
# saves the plot in my week 4 folder
ggsave("GitHub/qbio_data_analysis_tessa/week4_clinical/kmplot_by_KRAS.png", plot = p, width = 12, height = 9)


#-----------------------------------------------------------------
# 3) Relate your second variable to survival in colorectal cancer
#-----------------------------------------------------------------
# initializes a 'survival' object containing the clinic data
surv_object <- Surv(time = clinic$days_to_death, 
                    event = clinic$death_event)
# creates a fit object using the data from the family members w/ cancer column
fam_fit <- surv_fit( surv_object ~ clinic$number_of_first_degree_relatives_with_cancer_diagnosis, data = clinic )
# formatting for the KM plot
survplot = ggsurvplot(fam_fit, 
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
# saves the plot in my week 4 folder
ggsave("GitHub/qbio_data_analysis_tessa/week4_clinical/kmplot_by_fam.png", plot = p, width = 12, height = 9)