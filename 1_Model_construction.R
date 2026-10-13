
library(readxl)
require(xgboost)
require(Matrix)
require(data.table)
library(survival)
library(dplyr)
library(caret)
library(tidyverse)
library(MASS)
library(reportROC)
library(shapviz)
library(reshape2)
library(ggprism)
library(rms)
library(rmda)
library(tidyverse)
library(glmnet)
library(caret)
library(ROCR)
library(pROC)
agresti_coull_ci <- function(x, n, conf.level = 0.95) {
  z <- qnorm(1 - (1 - conf.level)/2)
  n_tilde <- n + z^2
  p_tilde <- (x + z^2/2) / n_tilde
  se <- sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  lower <- p_tilde - z * se
  upper <- p_tilde + z * se
  p=x/n
  i=paste0(x,"/",n)
  
  return(c(p,lower, upper))
}
library(PRROC)
library(boot)
# 定义计算AP值95%置信区间的函数,查全、查准、阈值
calculate_apci <- function(predictions, labels, R = 1000) {
  # 计算PR曲线
  pr_curve <- pr.curve(scores.class0 = predictions, weights.class0 = as.numeric(as.character(labels)),
                       curve = FALSE)
  
  # 计算AP值
  ap_value <- pr_curve$auc.integral
  # 自定义boot函数
  ap_boot <- function(data,indices) {
    indices=sample(length(as.numeric(data[,1])), replace = TRUE)
    d <- data[indices, ]
    pr_curve_boot <- pr.curve(scores.class0 = as.numeric(d[,1]), weights.class0 = as.numeric(as.character(d[,2])),
                              curve = FALSE)
    return(pr_curve_boot$auc.integral)
  }
  
  # 将数据放入数据框
  data <- data.frame(predictions, labels)
  
  # 使用boot函数计算AP值的95%置信区间
  set.seed(123)
  results <- boot(data, ap_boot, R = 1000)
  
  # 计算95%置信区间
  ci <- quantile(results$t, c(0.025, 0.975))
  
  return(list(ap = ap_value, ci = ci))
} 
create_calibration_curve <- function(truth, prob, n_bins = 10) {
  # 创建分箱
  bins <- cut(prob, 
              breaks = seq(0, 1, length.out = n_bins + 1),
              include.lowest = TRUE)
  # 计算每个分箱的统计量
  cal_data <- data.frame(truth = truth, prob = prob, bin = bins) %>%
    group_by(bin) %>%
    summarise(mean_prob = mean(prob),
              mean_truth = mean(truth),
              n = n(),
              .groups = 'drop')
  # 计算Brier分数（校准度量）
  brier_score <- mean((prob - truth)^2)
  # 绘制校准曲线
  p <- ggplot(cal_data, aes(x = mean_prob, y = mean_truth)) +
    #geom_point(aes(size = n), color = "blue", alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey") +
    geom_smooth(method = "loess", se = FALSE, color = "red") +
    labs(x = "Predicted Probability",
         y = "Observed Probability",
         title = "Calibration Curve")+
    #subtitle = paste("Brier Score:", round(brier_score, 4))) +
    xlim(0, 1) + ylim(0, 1) +
    theme_classic() +
    theme(legend.position = "bottom")
  
  return(p)
}

###*建立模型####
Sys.setenv(LANGUAGE="en")
cntrl = trainControl(
  method = "cv",
  number = 5,
  verboseIter = F,
  search = "random"
)
options(stringsAsFactors = F)
load("Model_construction.RData")
a1 <- column_to_rownames(a1,var = "id")
gene=c("cg04055528", "cg05374531", "cg06555256", "cg06791979", "cg06965637", 
       "cg07401103", "cg08080122", "cg09372617", "cg12768945", "cg13436782", 
       "cg13924807", "cg18272286", "cg19459049", "cg24206084", "cg24604851", 
       "cg24905336", "cg26932830", "cg27356587")
gtm=c("hp", "ratio")
train=a1[a1$gg=="train",-1]
test=a1[a1$gg=="test",-1]
wai=a1[a1$gg=="wai",-1]
qz=a1[a1$gg=="qz",-1]
train$sex=ifelse(train$sex == 1, 1, 0)
test$sex=ifelse(test$sex == 1, 1, 0)
wai$sex=ifelse(wai$sex == 1, 1, 0)
qz$sex=ifelse(qz$sex == 1, 1, 0)
train$hp=ifelse(train$hp == 1, 1, 0)
test$hp=ifelse(test$hp == 1, 1, 0)
wai$hp=ifelse(wai$hp == 1, 1, 0)
qz$hp=ifelse(qz$hp == 1, 1, 0)
####
trainx=as.matrix(train[,-1])
trainy=train[,1]
testx=as.matrix(test[,-1])
testy=test[,1]
waix=as.matrix(wai[,-1])
waiy=wai[,1]
qzx=as.matrix(qz[,-1])
qzy=qz[,1]
#####cg#####
trainx=as.matrix(train[,gene])
testx=as.matrix(test[,gene])
waix=as.matrix(wai[,gene])
qzx=as.matrix(qz[,gene])

trainy=ifelse(trainy=="GC","T","N")
set.seed(466)
tuned_model = train(x = trainx,y = trainy,trControl = cntrl,method = "xgbTree")
params <- as.list(tuned_model$bestTune[-1]) 
nrounds <- tuned_model$bestTune[[1]]
trainy <- ifelse(trainy == "T", 1, 0)
xgb <- xgboost(data = as.matrix(trainx), label = trainy, nrounds = nrounds, params = params, objective = "binary:logistic")
y_predicted.train <- predict(xgb, as.matrix(trainx))
C1=data.frame(class=trainy,cg=y_predicted.train)
rownames(C1)=rownames(train)

y_predicted.test<- predict(xgb, as.matrix(testx))
C2=data.frame(class=testy,cg=y_predicted.test)
rownames(C2)=rownames(test)

y_predicted.wai<- predict(xgb, as.matrix(waix))
C3=data.frame(class=waiy,cg=y_predicted.wai)
rownames(C3)=rownames(wai)

y_predicted.qz<- predict(xgb, as.matrix(qzx))
C4=data.frame(class=qzy,cg=y_predicted.qz)
rownames(C4)=rownames(qz)

roc0<-roc(class~cg,C1)
trcutoff0=data.frame(sensitivity=roc0[[2]],specificity=roc0[[3]],cutoff=roc0[[4]])
#yd_cutoff <- round(as.numeric(coords(roc0, "best", ret = "threshold")),2)
##2.sex+age
g00=c("sex","age")
trainx=as.matrix(train[,g00])
testx=as.matrix(test[,g00])
waix=as.matrix(wai[,g00])
qzx=as.matrix(qz[,g00])
set.seed(466)
trainy=ifelse(trainy == 1, "T", "N")
tuned_model = train(x = trainx,y = trainy,trControl = cntrl,method = "xgbTree")
params <- as.list(tuned_model$bestTune[-1]) 
nrounds <- tuned_model$bestTune[[1]]
trainy <- ifelse(trainy == "T", 1, 0)
xgb <- xgboost(data = as.matrix(trainx), label = trainy, nrounds = nrounds, params = params, objective = "binary:logistic")
y_predicted.train <- predict(xgb, as.matrix(trainx))
y_predicted.test<- predict(xgb, as.matrix(testx))
y_predicted.wai<- predict(xgb, as.matrix(waix))
y_predicted.qz<- predict(xgb, as.matrix(qzx))
C1$cli=y_predicted.train
C2$cli=y_predicted.test
C3$cli=y_predicted.wai
C4$cli=y_predicted.qz
###3.tm
trainx=as.matrix(train[,gtm])
testx=as.matrix(test[,gtm])
waix=as.matrix(wai[,gtm])
qzx=as.matrix(qz[,gtm])
set.seed(466)
trainy=ifelse(trainy == 1, "T", "N")
tuned_model = train(x = trainx,y = trainy,trControl = cntrl,method = "xgbTree")
params <- as.list(tuned_model$bestTune[-1]) 
nrounds <- tuned_model$bestTune[[1]]
trainy <- ifelse(trainy == "T", 1, 0)
xgb <- xgboost(data = as.matrix(trainx), label = trainy, nrounds = nrounds, params = params, objective = "binary:logistic")
y_predicted.train <- predict(xgb, as.matrix(trainx))
y_predicted.test<- predict(xgb, as.matrix(testx))
y_predicted.wai<- predict(xgb, as.matrix(waix))
y_predicted.qz<- predict(xgb, as.matrix(qzx))
C1$tm=y_predicted.train
C2$tm=y_predicted.test
C3$tm=y_predicted.wai
C4$tm=y_predicted.qz

#####二维(cg+cli)##########
set.seed(1)
trainy=ifelse(trainy == 1, "T", "N")
trainx=as.matrix(C1[,2:3])
colnames(trainx) <- colnames(C1)[2:3]
tuned_model = train(x = trainx,y = trainy,trControl = cntrl,method = "xgbTree")
params <- as.list(tuned_model$bestTune[-1]) 
nrounds <- tuned_model$bestTune[[1]]
trainy <- ifelse(trainy == "T", 1, 0)
xgb <- xgboost(data = as.matrix(trainx), label = trainy, nrounds = nrounds, params = params, objective = "binary:logistic")
testx=as.matrix(C2[,2:3])
waix=as.matrix(C3[,2:3])
qzx=as.matrix(C4[,2:3])
y_predicted.train <- predict(xgb, as.matrix(trainx))
y_predicted.test<- predict(xgb, as.matrix(testx))
y_predicted.wai<- predict(xgb, as.matrix(waix))
y_predicted.qz<- predict(xgb, as.matrix(qzx))

D2_1_yp1=data.frame(class=trainy,pre=y_predicted.train)
rownames(D2_1_yp1)=rownames(train)
y_predicted.test<- predict(xgb, as.matrix(testx))
D2_1_yp2=data.frame(class=testy,pre=y_predicted.test)
rownames(D2_1_yp2)=rownames(test)

y_predicted.wai<- predict(xgb, as.matrix(waix))
D2_1_yp3=data.frame(class=waiy,pre=y_predicted.wai)
rownames(D2_1_yp3)=rownames(wai)

y_predicted.qz<- predict(xgb, as.matrix(qzx))
D2_1_yp4=data.frame(class=qzy,pre=y_predicted.qz)
rownames(D2_1_yp4)=rownames(qz)
roc2_1_yp1<-roc(class~pre,D2_1_yp1)
#yd_cutoff2_1 <- round(as.numeric(coords(roc2_1, "best", ret = "threshold")),2)
trcutoff2_1=data.frame(sensitivity=roc2_1_yp1[[2]],specificity=roc2_1_yp1[[3]],cutoff=roc2_1_yp1[[4]])
#####三维###--------------------
set.seed(214)
trainy=ifelse(trainy == 1, "T", "N")
trainx=as.matrix(C1[,2:4])
#colnames(trainx) <- colnames(C1)[2:3]
tuned_model = train(x = trainx,y = trainy,trControl = cntrl,method = "xgbTree")
params <- as.list(tuned_model$bestTune[-1]) # 假设第一列是nrounds
nrounds <- tuned_model$bestTune[[1]]
trainy <- ifelse(trainy == "T", 1, 0)
xgb <- xgboost(data = as.matrix(trainx), label = trainy, nrounds = nrounds, params = params, objective = "binary:logistic")
testx=as.matrix(C2[,2:4])
waix=as.matrix(C3[,2:4])
qzx=as.matrix(C4[,2:4])
y_predicted.train <- predict(xgb, as.matrix(trainx))
y_predicted.test<- predict(xgb, as.matrix(testx))
y_predicted.wai<- predict(xgb, as.matrix(waix))
y_predicted.qz<- predict(xgb, as.matrix(qzx))

yp1=data.frame(class=trainy,pre=y_predicted.train)
rownames(yp1)=rownames(train)
y_predicted.test<- predict(xgb, as.matrix(testx))
yp2=data.frame(class=testy,pre=y_predicted.test)
rownames(yp2)=rownames(test)
####
y_predicted.wai<- predict(xgb, as.matrix(waix))
yp3=data.frame(class=waiy,pre=y_predicted.wai)
rownames(yp3)=rownames(wai)
####
y_predicted.qz<- predict(xgb, as.matrix(qzx))
yp4=data.frame(class=qzy,pre=y_predicted.qz)
rownames(yp4)=rownames(qz)
roc1<-roc(class~pre,yp1)  
#yd_cutoff3 <- round(as.numeric(coords(roc1, "best", ret = "threshold")),2)
trcutoff=data.frame(sensitivity=roc1[[2]],specificity=roc1[[3]],cutoff=roc1[[4]])

