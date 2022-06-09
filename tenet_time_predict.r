TENET_time_predict <-  read.delim("C:/Users/LSCS/Downloads/aada.txt")

for(i in 1:dim(TENET_time_predict)[2]){
  
  imsi = as.numeric(TENET_time_predict[,i]);   
  
  TENET_time_predict[,i] <- imsi;
  
}


str(TENET_time_predict)

model1 <- lm(formula = user_time ~ gene + thread + cell_select, data =  TENET_time_predict)

summary(model1)

model2 <- lm(formula = user_time ~ gene + cell_select, data =  TENET_time_predict)

summary(model2)

time = -3.431e+06 + 2.115e+02 * gene + 3.272e+03*cell_select



