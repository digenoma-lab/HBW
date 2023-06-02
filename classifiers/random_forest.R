#random forest 

install.packages("dplyr")
install.packages('e1071', dependencies=TRUE)
install.packages("caTools")
install.packages('randomForest')
install.packages('caret')

library("dplyr")
library("caTools")
library(randomForest)
library(caret)
library(e1071)
library(stringr)



data = read.csv("/path/to/data.txt",
                sep=" ")
data2 = read.csv("/path/to/data2.txt",
                 sep=" ")

summary(data) #resumen estadistico de los datos ingresados
summary(data2)


data$p[str_detect(data$id,"K1")] <-1 #1=B
data$p[str_detect(data$id,"K2")] <-0 #0=A

data <- transform(data, p=as.factor(p))

#ignorar el id
data <- data[,-1]

#dividir el dataset
sample <- sample(c(TRUE, FALSE), nrow(data), replace = TRUE, prob = c(0.05,0.95))
train <- data[sample, ]
test <- data[!sample, ]

train <- transform(train, p=as.factor(p))
test <- transform(test, p=as.factor(p))

#random forest
rf <- randomForest(p ~ .,data=train)

pred = predict(rf, newdata = test, type = "class")

confusionMatrix(table(pred,test$p))
