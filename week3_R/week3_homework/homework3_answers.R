#TESSA FERRARI
# exercise 1.1

is.na(attenu$station)

attenu_cleaned = na.omit(attenu, cols=3, invert=FALSE)

head(attenu_cleaned)
dim(attenu_cleaned)

# exercise 1.2

Theoph_2 = Theoph

median(Theoph_2$Dose)

doseClass = function(x, y){
  if (x>=y){return("high")} else{return("low")}
}
Theoph_2$Dose_Class = mapply(doseClass, Theoph_2$Dose, median(Theoph_2$Dose))

head(Theoph_2)
dim(Theoph_2)

# exercise 1.3

getwd()
setwd("GitHub/qbio_data_analysis_tessa/week3_R/week3_homework")
starbucks = read.csv("starbucks.csv")

is_row_empty = rowSums(is.na(starbucks))/6
nrow(starbucks)
length(is_row_empty)
mask = !is_row_empty
starbucks_cleaned = starbucks[mask, ]

plot(x = starbucks_cleaned$Carb, y = starbucks_cleaned$Calories,
     main = "Calories vs. Carb",
     xlab = "Carb (g)",
     ylab = "Calories")

Func1 = function(x,y){
  if(x==y){return(TRUE)}else{return(FALSE)}
}
isMaxCal = mapply(Func1, starbucks_cleaned$Calories, max(starbucks_cleaned$Calories))
MaxCal = starbucks_cleaned[isMaxCal,]
MaxCal[1,1]

starbucks_cleaned$is_highest_fat = mapply(Func1, starbucks_cleaned$Fat, max(starbucks_cleaned$Fat))
plot(x = starbucks_cleaned$Carb, y = starbucks_cleaned$Calories,
     col = factor(starbucks_cleaned$is_highest_fat),
     main = "Calories vs. Carb",
     xlab = "Carb (g)",
     ylab = "Calories")

##BONUS?

# exercise 1.4

batting = read.csv("Batting.csv")

Func2 = function(x){
  if(x>=3){return(TRUE)}else{return(FALSE)}
}
doesHave3orMoreHRs = mapply(Func2, batting$HR)

plot(x = batting$yearID, y = batting$HR,
     main = "Homeruns vs. Year",
     xlab = "Year",
     ylab = "Homeruns")

Func3 = function(x){
  if(x=="LAA"){return(TRUE)}else{return(FALSE)}
}
isLAAngles = mapply(Func3, batting$teamID)
LAAbatting = batting[isLAAngles,]
doesLAAHave3orMoreHRs = mapply(Func2, LAAbatting$HR)
plot(x = LAAbatting$yearID, y = LAAbatting$HR,
     main = "LA Angels Homeruns vs. Year",
     xlab = "Year",
     ylab = "Homeruns")

Func4 = function(x){
  if(x=="ATL"|x=="PIT"){return(TRUE)}else{return(FALSE)}
}
isATLorPIT = mapply(Func4, batting$teamID)
ATLandPITbatting = batting[isATLorPIT,]
Func5 = function(x){
  if(x=="ATL"){return(TRUE)}else{return(FALSE)}
}
ATLandPITbatting$isATL = mapply(Func5, ATLandPITbatting$teamID)
plot(x = ATLandPITbatting$yearID, y = ATLandPITbatting$HR,
     col = factor(ATLandPITbatting$isATL),
     main = "ATL and PIT Homeruns vs. Year",
     xlab = "Year",
     ylab = "Homeruns")

# exercise 1.5

easy_plot = function(x, y, color_data){
  median(color_data)
  Func6 = function(a,b){if(a>=b){return("high")}else{return("low")}}
  levels = mapply(Func6, color_data, median(color_data))
  levels = factor(levels)
  print(levels)
  
  plot(x = x, y = y,
       col = levels,
       pch = 20,)
  
  print(cor.test(x,y))
}
easy_plot(starbucks_cleaned$Carb,starbucks_cleaned$Calories,starbucks_cleaned$Carb)

easy_plot(starbucks_cleaned$Fat, starbucks_cleaned$Carb,starbucks_cleaned$Sodium)
easy_plot(batting$yearID, batting$G, batting$HR)

# exercise 2.1

iris
glimpse(iris)
dim(iris)
# The data set describes the species, sepal
# length/width, and petal length/width of irises.
# The data set contains 150 observations.
# The data set contains 5 features.

# exercise 2.2

# cols
# 1 - SepalL -  continuous  - double
# 2 - SepalW -  continuous  - double
# 3 - PetalL -  continuous  - double
# 4 - PetalW -  continuous  - double
# 5 - Species - catagorical - factor

# exercise 2.3

hist(iris$Sepal.Length)
hist(iris$Sepal.Width)
hist(iris$Petal.Length)
hist(iris$Petal.Width)
# The sepal lengths/widths have a relatively even 
# distribution but the petal lengths/widths have
# frequent small-petal outliers (maybe due to species?)

# exercise 2.4

meanSepWid = mean(iris$Sepal.Width)
iris_copy = iris

Func7 = function(x,y){
  if (x>=y){return("wide-sepaled")}else{return("narrow-sepaled")}
}
SepWid = mapply(Func7, iris_copy$Sepal.Width, meanSepWid)
iris_copy$SepWid = SepWid

boxplot(iris_copy$Sepal.Width ~ iris_copy$SepWid)

# exercise 2.5

# most unique: iris virginica
# most similar: iris setosa
pairs(iris_copy[,1:4], col=iris_copy$Species)

# exercise 3.1

install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library("BiocManager")
