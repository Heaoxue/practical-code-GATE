library(tidyverse)
library(glmnet)
library(caret)
library(ROCR)
library(reportROC)
library(pROC)
library(readxl)
library(randomForest)
library(bestglm)
library(tidymodels)
library(bonsai)
library(lightgbm)
library(e1071)
library(readxl)
require(xgboost)
require(Matrix)
require(data.table)
library(survival)
library(dplyr)
library(caret)
library(tidyverse)
#library(MASS)
library(ROCR)
library(reportROC)
library(shapviz)
library(reshape2)
library(ggprism)
library(rms)
library(rmda)
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

load("Model_construction.RData")
load("Performance.RData")
###Fig.3B####
#Training,Testing,Validation(yp1,yp2,yp3)
C2$class <- ifelse(C2$class == "GC", 1, 0) 
C3$class <- ifelse(C3$class == "GC", 1, 0)
D2_1_yp2$class <- ifelse(D2_1_yp2$class == "GC", 1, 0)
D2_1_yp3$class <- ifelse(D2_1_yp3$class == "GC", 1, 0)
yp2$class <- ifelse(yp2$class == "GC", 1, 0)
yp3$class <- ifelse(yp3$class == "GC", 1, 0)
roc1<-roc(class~cg,C1)
roc2<-roc(class~pre,D2_1_yp1)
roc3<-roc(class~pre,yp1)
ci1=ci(roc1)
ci2=ci(roc2)
ci3=ci(roc3)
plot(roc1, col = "#006BB9", main = "ML ROC", lwd = 2, lty = 1, 
     cex.lab = 1.2, identity = FALSE)  
abline(a = 1, b = -1, col = "gray", lty = "dashed", lwd = 2)
lines(roc2, col = "#fdae61", lwd = 2, lty = 1) 
lines(roc3, col = "#E3353A", lwd = 2, lty = 1)

legend("bottomright",
       legend = c(paste0("Methy:", round(ci1[2], 2), paste0("; 95% CI:", round(ci1[1],2), "-", round(ci1[3],2))), 
                  paste0("Methy+Clinic:", round(ci2[2], 2), paste0("; 95% CI:", round(ci2[1],2), "-", round(ci2[3],2))),
                  paste0("GCscan:", round(ci3[2], 2), paste0("; 95% CI:", round(ci3[1],2), "-", round(ci3[3],2)))), 
       
       col = c("#006BB9","#fdae61","#E3353A"), 
       lty = 1,lwd=1,cex = 0.85,bty = "n")


###*Shap图####
library(shapviz)
load("Performance.RData")
#训练集
shap_xgboost=shapviz(xgb,X_pred=as.matrix(trainx))
p <- sv_importance(shap_xgboost, kind = "beeswarm", 
                   show_numbers = TRUE,  # 显示特征数量
                   bee_width = 0.2,      # 减小宽度使点更密集
                   alpha = 1,          # 设置透明度
                   size = 1.2)           # 调整点的大小（可以调小使点更密集）

# 调整配色和布局
p + theme_bw() +
  scale_color_gradient(low = "blue", high = "#DE3044") +
  scale_fill_gradient(low = "blue", high = "#DE3044") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)  # 调整边距
  )

#测试集
shap_xgboost=shapviz(xgb,X_pred=as.matrix(testx))
p <- sv_importance(shap_xgboost, kind = "beeswarm", 
                   show_numbers = TRUE,  # 显示特征数量
                   bee_width = 0.15,      # 减小宽度使点更密集
                   alpha = 1,          # 设置透明度
                   size = 1.2)           # 调整点的大小（可以调小使点更密集）

# 调整配色和布局
p + theme_bw() +
  scale_color_gradient(low = "blue", high = "#DE3044") +
  scale_fill_gradient(low = "blue", high = "#DE3044") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)  # 调整边距
  )
#外部
shap_xgboost=shapviz(xgb,X_pred=as.matrix(waix))
p <- sv_importance(shap_xgboost, kind = "beeswarm", 
                   show_numbers = TRUE,  # 显示特征数量
                   bee_width = 0.15,      # 减小宽度使点更密集
                   alpha = 1,          # 设置透明度
                   size = 1.2)           # 调整点的大小（可以调小使点更密集）

# 调整配色和布局
p + theme_bw() +
  scale_color_gradient(low = "blue", high = "#DE3044") +
  scale_fill_gradient(low = "blue", high = "#DE3044") +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)  # 调整边距
  )
###*校准曲线####
# 如果尚未安装，请先运行：install.packages("CalibrationCurves")
library(CalibrationCurves)
# 你的数据预处理保持不变
#一维
C1$class <- ifelse(C1$class == 0, 0, 1)
C2$class <- ifelse(C2$class == "NGC", 0, 1)
C3$class <- ifelse(C3$class == "NGC", 0, 1)

#二维
D2_1_yp1$class <- ifelse(D2_1_yp1$class == 0, 0, 1)
D2_1_yp2$class <- ifelse(D2_1_yp2$class == "NGC", 0, 1)
D2_1_yp3$class <- ifelse(D2_1_yp3$class == "NGC", 0, 1)

#三维
yp1$class <- ifelse(yp1$class == 0, 0, 1)
yp2$class <- ifelse(yp2$class == "NGC", 0, 1)
yp3$class <- ifelse(yp3$class == "NGC", 0, 1)

# 一维
calibration_plot_1 <- valProbggplot(p = C1$cg, y = C1$class)
print(calibration_plot_1)

calibration_plot_2 <- valProbggplot(p = C2$cg, y = C2$class)
print(calibration_plot_2)

calibration_plot_3 <- valProbggplot(p = C3$cg, y = C3$class)
print(calibration_plot_3)

# 二维
calibration_plot_1 <- valProbggplot(p = D2_1_yp1$pre, y = D2_1_yp1$class)
print(calibration_plot_1)

calibration_plot_2 <- valProbggplot(p = D2_1_yp2$pre, y = D2_1_yp2$class)
print(calibration_plot_2)

calibration_plot_3 <- valProbggplot(p = D2_1_yp3$pre, y = D2_1_yp3$class)
print(calibration_plot_3)

# 三维
calibration_plot_1 <- valProbggplot(p = yp1$pre, y = yp1$class)
print(calibration_plot_1)

calibration_plot_2 <- valProbggplot(p = yp2$pre, y = yp2$class)
print(calibration_plot_2)

calibration_plot_3 <- valProbggplot(p = yp3$pre, y = yp3$class)
print(calibration_plot_3)

###Table6####
sen <- data.frame()  
spe <- data.frame()
auc <- data.frame()
accu <- data.frame()

load("Performance.RData")
a <- a[a$cohort %in% c("Training set","Testing set","Validation set"), 
       c(2:5, 11:15, 19)]
a <- a[!grepl("preGC", a$Group), ]
a$Hp <- factor(ifelse(a$Hp == "Yes", "yes", "no"), levels = c("yes", "no"))
yp <- rbind(yp1,yp2,yp3)
yp$id <- rownames(yp)
yp <- left_join(yp, a, by = "id")
yp[yp == "Unknown"] <- NA
yp <- yp[, -c(1, 3)]
yp$Group=factor(ifelse(yp$Group=="GC",1,0),levels = c(1,0))

yp <- yp %>%
  mutate(
    Age= ifelse(Age <= 60, "≤60", ">60"),
    Age = factor(Age, levels = c("≤60", ">60")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    Stage = factor(Stage, levels = c("I", "II", "III" ,"Ⅳ")),
    T_stage = factor(T_stage, levels = c("T1", "T2", "T3", "T4")),
    N_stage = factor(N_stage, levels = c("N0", "N1", "N2", "N3")),
    M_stage = factor(M_stage, levels = c("M0", "M1")),
    Size = factor(Size, levels = c("≤3.5cm", ">3.5cm")),
    Hp = factor(Hp, levels = c("yes", "no"))
  )

write.csv(yp,"YZ3.csv")

data <- yp
cutoff =0.43
calculate_metrics <- function(data_subset) {
  if (nrow(data_subset) == 0) {
    return(c(
      n = 0, 
      sensitivity = NA, specificity = NA, 
      sensitivity_ci_lower = NA, sensitivity_ci_upper = NA,
      specificity_ci_lower = NA, specificity_ci_upper = NA,
      auc = NA, auc_ci_lower = NA, auc_ci_upper = NA
    ))
  }
  
  data_subset <- data_subset[!is.na(data_subset$Group) & !is.na(data_subset$pre), ]
  
  if (nrow(data_subset) == 0 || length(unique(data_subset$Group)) < 2) {
    return(c(
      n = nrow(data_subset), 
      sensitivity = NA, specificity = NA,
      sensitivity_ci_lower = NA, sensitivity_ci_upper = NA,
      specificity_ci_lower = NA, specificity_ci_upper = NA,
      auc = NA, auc_ci_lower = NA, auc_ci_upper = NA
    ))
  }
  
  roc_obj <- roc(response = data_subset$Group, predictor = data_subset$pre, levels = c(0, 1), direction = "<")
  auc_val <- auc(roc_obj)
  ci_auc <- ci.auc(roc_obj)
  
  pred_class <- ifelse(data_subset$pre >= cutoff, 1, 0)
  cm <- table(Predicted = pred_class, Actual = data_subset$Group)
  
  sensitivity <- specificity <- NA
  sensitivity_ci_lower <- sensitivity_ci_upper <- NA
  specificity_ci_lower <- specificity_ci_upper <- NA
  
  if (all(c(0,1) %in% colnames(cm)) && all(c(0,1) %in% rownames(cm))) {
    TP <- cm["1", "1"]
    FN <- cm["0", "1"]
    sensitivity <- TP / (TP + FN)
    
    TN <- cm["0", "0"]
    FP <- cm["1", "0"]
    specificity <- TN / (TN + FP)
    
    if (TP + FN > 0) {
      n <- TP + FN
      p <- sensitivity
      z <- qnorm(0.975)
      sensitivity_ci_lower <- (p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
      sensitivity_ci_upper <- (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
    }
    
    
    if (TN + FP > 0) {
      n <- TN + FP
      p <- specificity
      z <- qnorm(0.975)
      specificity_ci_lower <- (p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
      specificity_ci_upper <- (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)
    }
  }
  

  return(c(
    n = nrow(data_subset),
    sensitivity = round(sensitivity, 4),
    sensitivity_ci_lower = round(sensitivity_ci_lower, 4),
    sensitivity_ci_upper = round(sensitivity_ci_upper, 4),
    specificity = round(specificity, 4),
    specificity_ci_lower = round(specificity_ci_lower, 4),
    specificity_ci_upper = round(specificity_ci_upper, 4),
    auc = round(auc_val, 4),
    auc_ci_lower = round(as.numeric(ci_auc[1]), 4),
    auc_ci_upper = round(as.numeric(ci_auc[3]), 4)
  ))
}

results <- data.frame()

# 1. Total
total_metrics <- calculate_metrics(data)
results <- rbind(results, data.frame(
  class = "Total",
  n = total_metrics["n"],
  sensitivity = paste0(
    sprintf("%.4f", total_metrics["sensitivity"]), " (",
    sprintf("%.4f", total_metrics["sensitivity_ci_lower"]), "-",
    sprintf("%.4f", total_metrics["sensitivity_ci_upper"]), ")"
  ),
  specificity = paste0(
    sprintf("%.4f", total_metrics["specificity"]), " (",
    sprintf("%.4f", total_metrics["specificity_ci_lower"]), "-",
    sprintf("%.4f", total_metrics["specificity_ci_upper"]), ")"
  ),
  auc = paste0(
    sprintf("%.4f", total_metrics["auc"]), " (",
    sprintf("%.4f", total_metrics["auc_ci_lower"]), "-",
    sprintf("%.4f", total_metrics["auc_ci_upper"]), ")"
  )
))

# 2. Age
age_le60 <- data %>% filter(Age == "≤60"| Group == 0)
age_gt60 <- data %>% filter(Age == ">60"| Group == 0)
age_na <- data %>% filter(!Age %in% c("≤60", ">60"))

metrics_le60 <- calculate_metrics(age_le60)
metrics_gt60 <- calculate_metrics(age_gt60)
metrics_age_na <- calculate_metrics(age_na)

add_results_row <- function(metrics, class_name) {
  data.frame(
    class = class_name,
    n = metrics["n"],
    sensitivity = ifelse(is.na(metrics["sensitivity"]), NA, 
                         paste0(
                           sprintf("%.4f", metrics["sensitivity"]), " (",
                           sprintf("%.4f", metrics["sensitivity_ci_lower"]), "-",
                           sprintf("%.4f", metrics["sensitivity_ci_upper"]), ")"
                         )),
    specificity = ifelse(is.na(metrics["specificity"]), NA,
                         paste0(
                           sprintf("%.4f", metrics["specificity"]), " (",
                           sprintf("%.4f", metrics["specificity_ci_lower"]), "-",
                           sprintf("%.4f", metrics["specificity_ci_upper"]), ")"
                         )),
    auc = ifelse(is.na(metrics["auc"]), NA,
                 paste0(
                   sprintf("%.4f", metrics["auc"]), " (",
                   sprintf("%.4f", metrics["auc_ci_lower"]), "-",
                   sprintf("%.4f", metrics["auc_ci_upper"]), ")"
                 ))
  )
}

results <- rbind(results, 
                 data.frame(class = "Age", n = NA, sensitivity = NA, specificity = NA, auc = NA),
                 add_results_row(metrics_le60, "≤ 60"),
                 add_results_row(metrics_gt60, "> 60")
)

# 3. Sex
sex_female <- data %>% filter(Sex == "Female"| Group == 0)
sex_male <- data %>% filter(Sex == "Male"| Group == 0)

metrics_female <- calculate_metrics(sex_female)
metrics_male <- calculate_metrics(sex_male)

results <- rbind(results, 
                 data.frame(class = "Sex", n = NA, sensitivity = NA, specificity = NA, auc = NA),
                 add_results_row(metrics_female, "Female"),
                 add_results_row(metrics_male, "Male")
)

# 4. TNM Stage 
stage_i <- data %>% filter(Stage == "I" | Group == 0)
stage_ii <- data %>% filter(Stage == "II"| Group == 0)
stage_iii <- data %>% filter(Stage == "III"| Group == 0)
stage_iv <- data %>% filter(M_stage == "M1" | Stage == "IV"| Group == 0)
stage_unknown <- data %>% filter(is.na(Stage) | Stage == "NA"| Group == 0)

metrics_i <- calculate_metrics(stage_i)
metrics_ii <- calculate_metrics(stage_ii)
metrics_iii <- calculate_metrics(stage_iii)
metrics_iv <- calculate_metrics(stage_iv)
metrics_stage_unknown <- calculate_metrics(stage_unknown)

results <- rbind(results, 
                 data.frame(class = "TNM stage", n = NA, sensitivity = NA, specificity = NA, auc = NA),
                 add_results_row(metrics_i, "I"),
                 add_results_row(metrics_ii, "II"),
                 add_results_row(metrics_iii, "III"),
                 add_results_row(metrics_iv, "IV"),
                 add_results_row(metrics_stage_unknown, "Unknown")
)

# 5. T stage
t1 <- data %>% filter(T_stage == "T1"| Group == 0)
t2 <- data %>% filter(T_stage == "T2"| Group == 0)
t3 <- data %>% filter(T_stage == "T3"| Group == 0)
t4 <- data %>% filter(T_stage == "T4"| Group == 0)
t_unknown <- data %>% filter(is.na(T_stage) | T_stage == "NA"| Group == 0)

metrics_t1 <- calculate_metrics(t1)
metrics_t2 <- calculate_metrics(t2)
metrics_t3 <- calculate_metrics(t3)
metrics_t4 <- calculate_metrics(t4)
metrics_t_unknown <- calculate_metrics(t_unknown)

results <- rbind(results, 
                 data.frame(class = "T stage", n = NA, sensitivity = NA, specificity = NA, auc = NA),
                 add_results_row(metrics_t1, "T1"),
                 add_results_row(metrics_t2, "T2"),
                 add_results_row(metrics_t3, "T3"),
                 add_results_row(metrics_t4, "T4"),
                 add_results_row(metrics_t_unknown, "Unknown")
)

# 6. N stage
n0 <- data %>% filter(N_stage == "N0"| Group == 0)
n1 <- data %>% filter(N_stage == "N1"| Group == 0)
n2 <- data %>% filter(N_stage == "N2"| Group == 0)
n3 <- data %>% filter(N_stage == "N3"| Group == 0)
n_unknown <- data %>% filter(is.na(N_stage) | N_stage == "NA"| Group == 0)

metrics_n0 <- calculate_metrics(n0)
metrics_n1 <- calculate_metrics(n1)
metrics_n2 <- calculate_metrics(n2)
metrics_n3 <- calculate_metrics(n3)
metrics_n_unknown <- calculate_metrics(n_unknown)

results <- rbind(results, 
                 data.frame(class = "N stage", n = NA, sensitivity = NA, specificity = NA, auc = NA),
                 add_results_row(metrics_n0, "N0"),
                 add_results_row(metrics_n1, "N1"),
                 add_results_row(metrics_n2, "N2"),
                 add_results_row(metrics_n3, "N3"),
                 add_results_row(metrics_n_unknown, "Unknown")
)

# 7. M stage
m0 <- data %>% filter(M_stage == "M0"| Group == 0)
m1 <- data %>% filter(M_stage == "M1"| Group == 0)
m_unknown <- data %>% filter(is.na(M_stage) | M_stage == "NA"| Group == 0)

metrics_m0 <- calculate_metrics(m0)
metrics_m1 <- calculate_metrics(m1)
metrics_m_unknown <- calculate_metrics(m_unknown)

results <- rbind(results, 
                 data.frame(class = "M stage", n = NA, sensitivity = NA, specificity = NA, auc = NA),
                 add_results_row(metrics_m0, "M0"),
                 add_results_row(metrics_m1, "M1"),
                 add_results_row(metrics_m_unknown, "Unknown")
)

# 8. Tumor Size
size_le35 <- data %>% filter(Size == "≤3.5cm"| Group == 0)
size_gt35 <- data %>% filter(Size == ">3.5cm"| Group == 0)
size_unknown <- data %>% filter(is.na(Size) | Size == "NA"| Group == 0)

metrics_le35 <- calculate_metrics(size_le35)
metrics_gt35 <- calculate_metrics(size_gt35)
metrics_size_unknown <- calculate_metrics(size_unknown)

results <- rbind(results, 
                 data.frame(class = "Tumor size", n = NA, sensitivity = NA, specificity = NA, auc = NA),
                 add_results_row(metrics_le35, "≤ 3.5cm"),
                 add_results_row(metrics_gt35, ">3.5cm"),
                 add_results_row(metrics_size_unknown, "Unknown")
)

# 9. H.pylori infection
hp_yes <- data %>% filter(Hp == "yes"| Group == 0)
hp_no <- data %>% filter(Hp == "no"| Group == 0)
hp_unknown <- data %>% filter(is.na(Hp) | Hp == "NA"| Group == 0)

metrics_hp_yes <- calculate_metrics(hp_yes)
metrics_hp_no <- calculate_metrics(hp_no)
metrics_hp_unknown <- calculate_metrics(hp_unknown)

results <- rbind(results, 
                 data.frame(class = "H.pylori infection", n = NA, sensitivity = NA, specificity = NA, auc = NA),
                 add_results_row(metrics_hp_yes, "Yes"),
                 add_results_row(metrics_hp_no, "No")
)


print(results)

convert_to_percent <- function(x) {

  parts <- strsplit(as.character(x), " ")[[1]]

  num <- as.numeric(parts[1])
  percent_val <- sprintf("%.2f%%", num * 100)

  if (length(parts) > 1) {
    return(paste(percent_val, paste(parts[-1], collapse = " ")))
  } else {
    return(percent_val)
  }
}

# 应用函数
results <- results %>%
  mutate(
    sensitivity = sapply(sensitivity, convert_to_percent),
    specificity = sapply(specificity, convert_to_percent)
  )
print(results)

write.csv(results, "subgroup_analysis_results.csv", row.names = FALSE)

convert_to_percent_with_ci <- function(x) {
  if(is.na(x) || x == "NA" || grepl("只包含一个类别", x)) {
    return(x)
  }
  
  pattern <- "([0-9.]+)\\s*\\(([0-9.]+)-([0-9.]+)\\)"
  
  if(grepl(pattern, x)) {

    matches <- regmatches(x, regexec(pattern, x))[[1]]
    main_value <- as.numeric(matches[2])
    ci_lower <- as.numeric(matches[3])
    ci_upper <- as.numeric(matches[4])
    
   
    result <- sprintf("%.2f%% (%.2f%%-%.2f%%)", 
                      main_value * 100, 
                      ci_lower * 100, 
                      ci_upper * 100)
    return(result)
  } else {

    main_value <- as.numeric(x)
    if(!is.na(main_value)) {
      return(sprintf("%.2f%%", main_value * 100))
    } else {
      return(x)
    }
  }
}

results <- results %>%
  mutate(
    sensitivity = sapply(sensitivity, convert_to_percent_with_ci),
    specificity = sapply(specificity, convert_to_percent_with_ci)
  )


###Fig.5A####
library(ggplot2)
library(ggpubr)
library(dplyr)
df <- yp4
cutoff1=0.43
data <- data.frame(group=df$class,score=as.numeric(df$pre))
# 从高到低排序数据
data <- data[order(-data$score),]
data <- data %>% mutate(level = ifelse(score >= cutoff1, "High-risk", "Low-risk"))
colnames(data)=c("status","predicted_risk_score","risk_group")
#data$status <- ifelse(data$status == 1, "GC", 
# ifelse(data$status == 0, "NGC", NA))
data$i=1:nrow(data)
points_data <- data.frame(
  x = c(max(which(data$risk_group=="High-risk"))-0.5),  # 手动设置 x 坐标
  y = c(cutoff1)  # 手动设置 y 坐标
)
p3 <- ggplot(data, aes(x = (i-1), y = predicted_risk_score)) +
  geom_area(aes(y = predicted_risk_score), fill = "#9EADD4", alpha = 0.5) +  # 填充曲线下的区域
  geom_smooth(method = "loess", span = 0.2, color = "#9EADD4", size = 1, se = FALSE) +  # 减小span值使曲线更贴合
  # 添加从点向左和向上的线段
  geom_segment(data = points_data, 
               aes(x = x, xend = 0, 
                   y = y, yend = y), 
               linetype = "dashed", color = "gray", linewidth = 1) +  # 向左的水平线
  geom_segment(data = points_data, 
               aes(x = x, xend = x, 
                   y = y, yend = 1), 
               linetype = "dashed", color = "gray", linewidth = 1) +  # 向上的垂直线
  geom_point(data = points_data, aes(x = x, y = y), 
             color = "black", size = 2.5) +
  labs(y = "Predicted risk score", x = "") +
  theme_classic() +  # 使用经典主题，没有网格线
  theme(
    panel.border = element_blank(),  # 去掉面板边框
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.line = element_line(color = "black", linewidth = 1),  # 设置坐标轴线的粗细
    axis.ticks.y = element_line(color = "black", linewidth = 1),
    panel.grid = element_blank()
  ) +
  scale_x_continuous(expand = c(0, 0)) +  # 设置 x 轴从 0 开始
  scale_y_continuous(expand = c(0, 0))    # 设置 y 轴从 0 开始
p3
# 绘制顶部的高危病变条形图
#GC
data$status=factor(ifelse(data$status=="GC","Yes","No"),levels = c("Yes","No"))
data=pivot_longer(data,cols=c("status","risk_group"),names_to = "status",values_to = "value")
p1 <- ggplot(data, aes(x = (i-1), fill = value)) +
  geom_tile(aes(y = status, height = 0.8, width = 1),size=2.5) +
  scale_fill_manual(values = c( "Yes" ="#50669C","No" = "#B2D0C0","High-risk" = "#FF6B6B", "Low-risk" ="#B8D0E7" )) +#"#B20D35"
  labs(y = "", x = "") +
  theme_classic() + 
  theme(
    panel.border = element_blank(), 
    axis.line = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank()
  )
p1

p13 <- ggarrange(
  p1, p3,
  nrow = 2,
  heights = c(1, 2)  
)
p13


###Fig.5E####
library(ggalluvial)
library(ggplot2)
load("Performance.RData")
q=column_to_rownames(q,var="id")
cutoff2=0.43
q=q[q$gg=="qz",]
score_tbl <- q %>% 
  mutate(
    age_score = case_when(
      between(age, 0, 49) ~ 0,
      between(age, 50, 59) ~ 5,
      between(age, 60, 69) ~ 6,
      age >= 70            ~ 10,
      TRUE                 ~ NA_real_),
    sex_score = 1 * (sex == 1),
    hp_score  = 3 * (hp == 1),
    ratio_score = 4 * (ratio > 3),
    G17_score = case_when(
      G17 < 1.5           ~ 0,
      between(G17, 1.5, 5.7) ~ 3,
      G17 > 5.7           ~ 5,
      TRUE                ~ NA_real_),
    total_score = age_score + sex_score + hp_score + ratio_score + G17_score
  )
score_nABC <- q %>% 
  mutate(
    G17_score = ifelse(G17>1&G17<15,0,1),
    PG_score = ifelse(PG1<=70&ratio<=7,1,0),
    total_score = ifelse(G17_score==1&PG_score==1,"High-risk",ifelse(G17_score==0&PG_score==0,"Low-risk","Medium-risk"))
  )

df = data.frame(
  Histopathology = q$group,
  Model = ifelse(yp4$pre > cutoff2, "High risk", "Low risk"), 
  ABC = ifelse(q$hp == 1, "Medium-risk", 
               ifelse(q$PG1 <= 70 & q$ratio <= 3, "High-risk", "Low-risk")),
  newscore = ifelse(score_tbl$total_score >= 17, "High-risk", 
                    ifelse(score_tbl$total_score <= 11, "Low-risk", "Medium-risk"))
)

df[,3:4] <- lapply(df[,3:4], function(x) {
  factor(x, levels = c("High-risk", "Medium-risk", "Low-risk"))
})

df$Model <- factor(df$Model, levels = c("High risk", "Low risk"))
df$Histopathology <- factor(df$Histopathology, levels = c("GC", "preGC", "NGC"))

colors = c(
  "High-risk" = "#FAB79B",
  "Medium-risk" = "#C8E0B4",
  "Low-risk" = "#9EADD4",
  "High risk" = "#FC6057", 
  "Low risk" = "#B8D0E7",   
  "preGC" = "#CCA1C5",
  "NGC" = "#DDDFE0",
  "GC" = "#8E3389"
)

plot_data <- df %>%
  dplyr::count(ABC, newscore, Model, Histopathology)

ggplot(plot_data,
       aes(y = n, axis1 = Histopathology, axis2 = Model,
           axis3 = ABC, axis4 = newscore)) +
  geom_alluvium(aes(fill = Histopathology), width = 2/5, alpha = 0.8) +  # alpha移到外面
  geom_stratum(
    aes(fill = after_stat(stratum)),
    width = 2/5,
    color = NA,
    size = 0.5
  ) +
  scale_fill_manual(values = colors) +
  geom_label(stat = "stratum", 
             aes(label = paste0(round(after_stat(count)/nrow(df)*100, 2), "%")),
             label.size = 0,
             fill = NA,
             size = 3) +
  scale_x_discrete(limits = c("Histopathology", "Model", "ABC", "newscore"),
                   expand = c(0.05, 0.05)) +
  theme_classic() +
  theme(legend.position = "top") +
  labs(title = "",
       y = "Number of Patients")

###Fig.5C####
cutoff1=0.43
df=data.frame(Histopathology=q$group,Model=ifelse(yp4$pre<=cutoff1,"Low risk",ifelse(yp4$pre>cutoff2,"High risk","Medium risk")),
              ABC=ifelse(q$hp==1,"Medium-risk",ifelse(q$PG1<=70&q$ratio<=3,"High-risk","Low-risk")),
              newscore=ifelse(score_tbl$total_score>=17,"High-risk",ifelse(score_tbl$total_score<=11,"Low-risk","Medium-risk")))
df[,3:4] <- lapply(df[,3:4], function(x) {
  factor(x, levels = c( "High-risk","Medium-risk","Low-risk"))
})
df$Model=factor(df$Model,levels = c( "High risk","Medium risk","Low risk") )
df$Histopathology=factor(df$Histopathology,levels = c("GC","preGC","NGC"))
df_clean <- list('GC+pre'=c("preGC","GC"),GC=c("GC"))
names(df_clean)=c('GC+pre',"GC")

da=data.frame()
for (h in 1:2) {
  i=df_clean[[h]]
  data <- df %>%
    mutate(Histopathology = ifelse(Histopathology %in% i, "Yes", "No"))
  da2=data.frame()
  for (i in 2:4) {
    da1 <- data %>%
      dplyr::count(!!sym(colnames(data)[i]), Histopathology)
    colnames(da1)[1]="group"
    da1$panel=colnames(df)[i]
    da2=rbind(da2,da1)
  }
  da2$name=names(df_clean[h])
  da=rbind(da,da2)
}
da1=data.frame()
for (i in seq(1,nrow(da),2)) {
  data_summary=da[c(i:(i+1)),]
  set.seed(225)
  ci=agresti_coull_ci(data_summary[2,3], sum(data_summary$n))
  prevalence=ci[1]
  lower=ci[2]
  upper=ci[3]
  panel=data_summary[1,4]
  name=data_summary[1,5]
  group=as.character(data_summary[1,1])
  No.Tatal.=paste0(data_summary[2,3],"/",sum(data_summary$n))
  da1=rbind(da1,cbind(name=name,panel=panel,group=group,"No./Total No."= No.Tatal.,prevalence=prevalence,lower=lower,upper=upper))
}
da1$panel=factor(da1$panel,levels = c("Model","ABC","newscore"))
da1$group=factor(da1$group,levels = c("Low risk","Medium risk","High risk","Low-risk","Medium-risk","High-risk"))
da1[,5:7]=apply(da1[,5:7],2,as.numeric)
# 发病率置信区间
library(ggplot2)
da=da1[da1$name=="GC",]
h=df_clean[[2]]
da <- da1[da1$name == "GC", ]

da_filtered <- da %>%
  mutate(group_clean = case_when(
    group %in% c("High-risk", "High risk") ~ "High-risk",
    group %in% c("Low-risk", "Low risk") ~ "Low-risk",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group_clean)) %>%
  
  mutate(panel_risk = paste(panel, group_clean, sep = "_"))

da_filtered$panel_risk <- factor(
  da_filtered$panel_risk,
  levels = c(
    "Model_High-risk", "ABC_High-risk", "newscore_High-risk",
    "Model_Low-risk", "ABC_Low-risk", "newscore_Low-risk"
  )
)

ggplot(da_filtered, aes(x = panel_risk, y = prevalence, color = group_clean)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25, linewidth = 1.05) +
  geom_text(
    aes(label = sprintf("%.2f%%", prevalence * 100), 
        y = upper + 0.01),  
    vjust = 0,  
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("High-risk" = "#E3353A", "Low-risk" = "#006bb9")
  ) +
  geom_hline(
    yintercept = sum(df$Histopathology %in% h) / nrow(df),
    linetype = "dashed", color = "grey"
  ) +
  labs(x = NULL, y = "Prevalence") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))  

###Fig.5D####
library(dplyr)
library(tidyr)
library(ggplot2)

da = data.frame()
for (i in 2:4) {
  da1 <- df %>%
    dplyr::count(!!sym(colnames(df)[i]), Histopathology)
  colnames(da1)[1] = "group"
  da1$panel = colnames(df)[i]
  da = rbind(da, da1)
}

da$panel = factor(da$panel, levels = c( "Model","ABC", "newscore"))

da$Histopathology = factor(da$Histopathology, levels = c("GC", "preGC", "NGC"))

da$group <- as.character(da$group)

da_model <- da %>% 
  filter(panel == "Model") %>%
  mutate(group = factor(group, levels = c("Low risk", "High risk")))

da_other <- da %>% 
  filter(panel != "Model") %>%
  mutate(group = factor(group, levels = c("Low-risk", "Medium-risk", "High-risk")))

da <- bind_rows(da_model, da_other)
da <- da %>% 
  arrange(panel, Histopathology)

label <- da %>% 
  group_by(panel, group) %>% 
  dplyr::summarise(
    labels = paste(Histopathology, n, sep = "=", collapse = "\n"),
    .groups = "drop"
  )

da_total <- da %>% 
  group_by(panel, group) %>% 
  dplyr::summarise(total = sum(n), .groups = "drop")

label <- left_join(label, da_total, by = c("panel", "group"))

ggplot(da, aes(x = group, y = n, fill = Histopathology)) +
  geom_col(width = 0.7, colour = NA)  +
  facet_wrap(~panel, scales = "free_x", strip.position = "bottom") + 
  labs(title = "", y = "Number of patients", x = "") +
  theme_classic() +
  scale_y_continuous(limits = c(0.0, 1800), expand = c(0, 0)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    panel.grid = element_blank(),
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key = element_blank(),
    legend.spacing = unit(0.1, "cm")
  ) +
  scale_fill_manual(values = c("NGC" = "#E2E4E5", "preGC" = "#A88BC1", "GC" = "#830B70")) +
  geom_text(
    data = label,
    aes(x = group, y = total + 20, label = labels),
    inherit.aes = FALSE,
    vjust = 0, size = 3, fontface = "bold"
  )
###Fig.5F####
library(pROC)
load("Performance.RData")
df2=column_to_rownames(df2,var="id")
#GC
data <- df2 %>% 
  filter(Histopathology != "preGC")
data$Actual <- ifelse(data$Histopathology == "GC", 1, 0)
data$ABC1 <- as.numeric(factor(data$ABC, levels = c("Low-risk", "Medium-risk", "High-risk")))
data$newscore1 <- as.numeric(factor(data$newscore, levels = c("Low-risk", "Medium-risk", "High-risk")))

model <- glm(Actual ~ ABC1, data = data, family = binomial)
model2 <- glm(Actual ~ newscore, data = data, family = binomial)

data$predicted_score_ABC <- predict(model, newdata = data, type = "response")
data$predicted_score_new <- predict(model2, newdata = data, type = "response")

roc_obj <- roc(data$Actual, data$predicted_score_ABC)
plot(roc_obj, main = "ROC Curve for ABC Method", col = "#1c61b6")
text(0.6, 0.2, paste("AUC =", round(auc(roc_obj), 2)), col = "red", cex = 1.2)

yp4_clean <- yp4 %>% filter(class != "preGC")
yp4_clean$class <- ifelse(yp4_clean$class == "GC", 1, 0)

roc1<-roc(class~pre,yp4_clean)
roc2<-roc(Actual~predicted_score_ABC,data)
roc3<-roc(Actual~predicted_score_new,data)
ci1=ci(roc1)
ci2=ci(roc2)
ci3=ci(roc3)

plot(roc1, col = "#E3353A", main = "ML ROC", lwd = 2, lty = 1, 
     cex.lab = 1.2, identity = FALSE)  # 添加 identity = FALSE

abline(a = 1, b = -1, col = "gray", lty = "dashed", lwd = 2)
lines(roc2, col = "#1A487B",lwd=2,lty=1) 
lines(roc3, col = "#F3B91E",lwd=2,lty=1) 
legend("bottomright",
       legend = c(paste0("Model:", round(ci1[2], 3), paste0("; 95% CI:", round(ci1[1],3), "-", round(ci1[3],3))), 
                  paste0("ABC:", round(ci2[2], 2), paste0("; 95% CI:", round(ci2[1],2), "-", round(ci2[3],2))),
                  paste0("newscore:", round(ci3[2], 2), paste0("; 95% CI:", round(ci3[1],2), "-", round(ci2[3],2)))), 
       col = c("#DF1F25","#1A487B","#F3B91E"), 
       lty = 1,lwd=2,cex = 0.85,bty = "n")


