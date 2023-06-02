library(MASS)
library(e1071)
library(stringr)

data2 = read.csv("/Users/cvalenzuela/GIT/HBW/classification/kidONT.R94.kmers.txt",
                 sep=" ")
data = read.csv("/Users/cvalenzuela/GIT/HBW/classification/kidONT.R94.kmers.uniq.txt",
                sep=" ")

data$p[str_detect(data$id,"K1")] <-1 #1=B
data$p[str_detect(data$id,"K2")] <-0 #0=A

data <- transform(data, p=as.factor(p))

data <- data[,-1]
C <- 10
C <- 10^3

svm.lineal <- svm(p~.,data=data, kernel='linear', cost=C, cross=2, scale=FALSE)
summary(svm.lineal)

svm.cuad <- svm(p~., data = data, kernel = 'polynomial',
                degree = 2, gamma = 1, coef0 = 1, cost = C, cross = 2, scale= FALSE)

