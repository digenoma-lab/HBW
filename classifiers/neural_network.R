install.packages("tidyverse")
install.packages("neuralnet")

library(tidyverse)
library(neuralnet)

data = read.csv("/Users/cvalenzuela/GIT/HBW/classification/kidONT.R94.kmers.txt",
                sep=" ")
data2 = read.csv("/Users/cvalenzuela/GIT/HBW/classification/kidONT.R94.kmers.uniq.txt",
                 sep=" ")

data$p[str_detect(data$id,"K1")] <-1 #1=B
data$p[str_detect(data$id,"K2")] <-0 #0=A

data <- transform(data, p=as.factor(p))

data <- data[,-1]
sample <- sample(c(TRUE, FALSE), nrow(data), replace = TRUE, prob = c(0.05,0.95))
train <- data[sample, ]
test <- data[!sample, ]

model = neuralnet(p~ha_k18+hb_k18+ha_k21+hb_k21+ha_k24+hb_k24, data= train, hidden=c(5,3,2), linear.outpu = FALSE)

plot(model,rep = "best")

pred <- predict(model,test)
labels <- c("A","B")
prediction_label <- data.frame(max.col(pred)) %>%
  mutate(pred=labels[max.col.pred.]) %>%
  select(2) %>%
  unlist()

table(test$p,prediction_label)

check = as.numeric(test$p) == max.col(pred)
accuracy = (sum(check)/nrow(test))*100
print(accuracy)
