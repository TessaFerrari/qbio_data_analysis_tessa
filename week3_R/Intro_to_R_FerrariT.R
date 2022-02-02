#-------------
#exercise 2.1
#-------------

name = "Tessa" #string
age = 20 #numeric
birthday = "07/17/2001" #string

#-------------
#exercise 2.2
#-------------

list_of_numbers = c(1, 3, 4, 8, 1, 44, -3, 29, -12)
min(list_of_numbers)# min() = -12
max(list_of_numbers)# max() = 44

#-------------
#example 2.3
#-------------

graph_parabola = function(a, b, c){
  curve(a*x^2 + b*x + c, xlim = c(-10, 10), ylim = c(-10, 10), 
        col = "blue", ylab = "y")
  # there is no return statement -- the output is in the plot window!
}

graph_parabola(2, 4, 3) # much easier than typing the above code!
graph_parabola(1, 0, 0)

#-------------
#exercise 2.3
#-------------

quad_form = function(a, b, c){
  x = (-b+sqrt(b*b - (4*a*c)))/(2*a)
  return(x)
}
quad_form(1, 0, -1) # equal to 1

#-------------
#example 2.4
#-------------
  
if (1 > 2){
  print("The if statement is true.")
} else if(2 > 1){
  print("The else if statement is true.") #this line will run
} else{
  print("Neither the if NOR the else if is true.")
}

#-------------
#exercise 2.4
#-------------

print("Green eggs and ham.")

if(1 > 2 | "red" == "blue"){
  print("I do not like green eggs and ham.")
}

if("red" == "green" | 2 > 1){
  print("I do not like them, Sam-I-Am.")
}

# note: when you have a (NOT) in front of a bunch of things,
# you can evaluate from inside-out
if(!(3 != 3 | "pass" == "fail")){
  print("Would you like them here or there?")
}

x = -1
y = 3

if(x < y & x > y){
  print("I would not like them here or there.")
} else if (x < y){
  print("I would not like them anywhere.")
} else if(x > y){
  print("Would you like them in a house?")
} else{
  print("Would you like them with a mouse?")
}
#output:
#Green eggs and ham.
#I do not like them, Sam-I-Am.
#Would you like them here or there?
#I would not like them anywhere.

#-------------
#exercise 2.5
#-------------

x = 0
while (x<=5){
  x = x + runif(1)
  print(x)
}

#-------------
#exercise 2.6
#-------------

for (i in 1:5){
  print("USC")
}

#-------------
#exercise 2.7
#-------------

print(10:1) # way 1

my_vec = 10:1 # way 2
print(my_vec)

#-------------
#exercise 2.8
#-------------

for(i in 1:6){
  for(j in 1:6){
    print(i+j)
  }
}

#-------------
#exercise 2.9
#-------------

mean_mat = matrix(1:100, nrow = 10, ncol = 10)  # data frame of values

rMeans = rowMeans(mean_mat)
cMeans = colMeans(mean_mat)

#-------------
#exercise 3.1
#-------------

library(tidyverse)

print("(1)")
str(mtcars)

print("(2)")
head(mtcars)

print("(3)")
glimpse(mtcars)

#-------------
#exercise 3.2
#-------------

mtcars_4 = mtcars$hp
mtcars_4 = mtcars[ , 4]

#-------------
#exercise 3.3
#-------------

mtcars$hp[3]
mtcars[3,4]

#-------------
#exercise 3.4
#-------------

head(mtcars)

mtcars$wt[5]
mtcars[5,6]

#-------------
#exercise 3.5
#-------------

mtcars[,1:3]
mtcars[1:2,8:11]
evens = 1:16 *2
mtcars[evens,]

#-------------
#exercise 3.6
#-------------

mtcars_copy = mtcars
mtcars_copy$mpg = mtcars_copy$mpg * 2 
mtcars_copy$super_mpg = mtcars_copy$mpg * 100 
glimpse(mtcars_copy)

mtcars_copy$mpg = mtcars_copy$mpg / 2 
mtcars_copy$super_mpg = mtcars$mpg

#-------------
#exercise 3.7
#-------------
