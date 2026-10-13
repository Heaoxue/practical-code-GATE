
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(Iso)
library(readxl)
library(stringr)
library(diagram)
library(ggalluvial)
library(cowplot)
compliance = 1
#iso_join_seer <- read_excel("D:/260408GC/China_Data.xlsx",sheet = 2)
load("Interception.RData")
#get my functions
#slip rate computation

#exact: integrate over those missed during a screening interval
#see supplemental information for argument that this function is appropriate
integrate_slip_rate<-function(screen_interval, weibull_shape, dwell,compliance){
  # slip rate
  # yield of escape = integral cumululate distribution function F(t), 0<t<screen_interval
  # mean of weibull = scale*gamma(1+1/shape)
  dwell_scale=dwell/gamma(1+1/weibull_shape)
  tiny_delta<-365
  #trapezoidal integration
  days_low<-seq(0,screen_interval*tiny_delta-1,by=1)/tiny_delta
  days_hi<-days_low+1/tiny_delta
  F_by_day<-0.5*(pweibull(days_low,shape=weibull_shape,scale=dwell_scale)+pweibull(days_hi,shape=weibull_shape,scale=dwell_scale))
  escaped_yield<-sum(F_by_day)*(1/tiny_delta) #day width in years
  #convert to slip rate: how many missed
  total_yield<-screen_interval #total incidence is just duration
  slip_rate<-escaped_yield/total_yield
  slip_rate<-pmin(1,(slip_rate+1-compliance))
  slip_rate
}

#use integrated weibull cumulative distribution functions
#"exact" solution to slip rate
exact_slip_rate_from_dwell<-function(dwell_model_all_df,screen_interval=1,weibull_shape=1,compliance){
  
  #slip is "before clinical"
  #slip_clinical = "at stage of clinical detection"
  #assume expected is half-duration of stage of clinical detection
  #completeness in modeling
  dwell_slip_df<-dwell_model_all_df %>%
    mutate(slip = sapply(dwell,function(z){integrate_slip_rate(screen_interval,weibull_shape,z,compliance)}),
           slip_clinical=sapply(dwell*0.5,function(z){integrate_slip_rate(screen_interval,weibull_shape,z,compliance)}),
           screen_interval=screen_interval)
  
  dwell_ideal_df <- dwell_slip_df %>%
    filter(scenario==1) %>%
    mutate(dwell=10000,slip=1-compliance,slip_clinical=0,scenario=as.integer(0),screen_interval=NA)
  
  #add "scenario 0" perfect interception
  dwell_slip_df<-bind_rows(dwell_slip_df,dwell_ideal_df)
  dwell_slip_df
}

xStage<-c("I","II","III","IV")

#okay: this multiplies a final destination IR: all the people who wind up at some final destination
#so everything scales within itself
#assume stage IV destination unless the cases are stopped earlier
effective_sens<-function(cumulative_sens, dwell_detect_rate){
  detect<-rep(0,length(cumulative_sens))
  miss<-0
  i<-1
  arrive<-cumulative_sens[i]
  live<-arrive+miss
  detect[1]<-cumulative_sens[i]*dwell_detect_rate[i]
  miss<-cumulative_sens[i]*(1-dwell_detect_rate[i])
  if (length(cumulative_sens)>1){
    for (i in 2:length(cumulative_sens))
    {
      #newly detectable cases
      arrive<-cumulative_sens[i]-cumulative_sens[i-1]
      live<-arrive+miss #currently detectable is newly detectable + missed at earlier stages
      detect[i]<-live*dwell_detect_rate[i] #would detect all of them, but miss some of them due to timing
      miss<-live*(1-dwell_detect_rate[i])
    }
  }
  arrive<-1-cumulative_sens[i] #final miss
  live=arrive+miss
  #detect[1-4] is detected at each stage
  #clinical detect= miss[4]
  list(intercept=detect,clinical=live)
}

add_survival_to_stage_shift<-function(incidence_sens_source,incidence_intercepted){
  #add survival: original before shift, and shifted survival
  #using 5 year survival as "crude estimate of statistical cure"
  #unlikely to be affected by 1-2 year lead time significantly
  #extract survival by 'stage_at_detection'
  just_survival<-incidence_sens_source %>%
    mutate(prequel=match(Stage,xStage)) %>%
    select(Cancer,prequel,Survival) %>%
    filter(!is.na(prequel))
  
  intercept_survival<-incidence_intercepted %>%
    left_join(just_survival %>%
                select(Cancer,prequel,s_survival=Survival),by=c("Cancer","prequel")) %>% 
    left_join(just_survival %>% 
                select(Cancer,clinical=prequel,c_survival=Survival),by=c("Cancer","clinical")) 
  
  #compute absolute numbers rather than local rates
  intercept_survival<-intercept_survival %>%
    mutate(original_survivors=c_survival*caught,
           shifted_survivors=s_survival*caught,
           original_deaths=(1-c_survival)*caught,
           shifted_deaths=(1-s_survival)*caught)
  intercept_survival
}

compute_effective_detection_with_slip<-function(incidence_sens_source,dwell_slip_df, active_slip_clinical){
  #just detection rate
  just_detection<-incidence_sens_source %>%
    mutate(number_stage=match(Stage,xStage),
           prequel=number_stage,
           detect=iso_sens) %>%
    select(Cancer,prequel,detect)
  
  #differences - marginal detection rate of remaining cases 
  #given that cases already detectable at earlier stage were removed or treated separately
  just_delta<- just_detection %>%
    group_by(Cancer) %>%
    arrange(prequel,.by_group=TRUE) %>%
    mutate(delta_detect=diff(c(0,detect))) %>%
    ungroup() %>%
    arrange(Cancer)
  
  #modify using slip rate
  #intercept using slip rate
  #slip to next
  
  #include modification of slip rate by clinical stage of detection
  #extra 'parameter'
  just_slip_delta_extra<-just_delta %>%
    left_join(dwell_slip_df %>% 
                #filter(scenario==dw_scenario) %>% 
                select(Cancer,prequel=number_stage,slip,slip_clinical),by=c("Cancer","prequel")) %>%
    filter(!is.na(prequel)) %>%
    mutate(unroll=4) %>%
    uncount(unroll,.id="clinical") %>%
    filter(clinical>=prequel) %>%
    mutate(modified_slip=case_when(prequel<clinical ~ slip,
                                   prequel==clinical & active_slip_clinical ~ slip_clinical,
                                   prequel==clinical & !active_slip_clinical ~ slip,
                                   TRUE ~ 1.0)) %>%
    arrange(Cancer,clinical,prequel) %>%
    group_by(Cancer,clinical) %>%
    mutate(sens_slip=effective_sens(detect,1-modified_slip)$intercept) %>%
    ungroup()
  
  just_slip_delta_extra
}

run_intercept_model<-function(incidence_sens_source, dwell_slip_df, active_slip_clinical=TRUE){
  
  #set up all previous stages where cases could be intercepted given clinical detection
  incidence_set<-incidence_sens_source %>% 
    filter(Stage %in% xStage) %>%
    select(Cancer,Stage,IR) %>%
    mutate(number_stage=match(Stage,xStage),
           clinical=number_stage,
           unroll=number_stage) %>%
    uncount(unroll,.id="prequel")
  
  #compute effective detection by stage conditional on slip rate model
  just_slip_delta_extra<-compute_effective_detection_with_slip(incidence_sens_source,
                                                               dwell_slip_df, 
                                                               active_slip_clinical)
  
  #updated: split effective slip rate in 2 for last stage
  #as the "time spent" should be halved
  #this involves a more elaborate model
  #note that "lives saved" is not affected, because those individuals are not stage shifted
  #this assumes that 'stage 4' is just automatically halved anyway
  incidence_intercepted<-incidence_set %>%
    left_join(just_slip_delta_extra,by=c("Cancer","clinical","prequel")) %>% 
    mutate(unroll=1+(number_stage==prequel)) %>%
    uncount(unroll,.id="found_clinical") %>%
    group_by(Cancer,clinical) %>%
    mutate(c_slip=cumsum(sens_slip),
           delta_detect=case_when(
             found_clinical==2 ~ 1-c_slip+sens_slip, #anyone not caught by new screening must be found clinically
             TRUE ~ sens_slip)) %>%
    mutate(caught=IR*delta_detect) %>%
    ungroup()
  
  intercept_survival<-add_survival_to_stage_shift(incidence_sens_source,incidence_intercepted)
  
  intercept_survival
}


run_excluded_model<-function(excluded_source){
  #fills out individuals not staged as they are not modeled
  excluded_survival<-excluded_source %>%
    mutate(number_stage=0,
           clinical=0,
           prequel=0,
           detect=0.0,
           delta_detect=0.0,
           slip=1.0,
           slip_clinical=1.0,
           modified_slip=1.0,
           sens_slip=0.0,
           found_clinical=2,
           c_slip=1.0,
           caught=IR,
           s_survival=Survival,
           c_survival=Survival) %>%
    mutate(original_survivors=c_survival*caught,
           shifted_survivors=s_survival*caught,
           original_deaths=(1-c_survival)*caught,
           shifted_deaths=(1-s_survival)*caught) %>%
    select(Cancer,Stage,IR,
           number_stage,clinical,prequel,
           detect,delta_detect,
           slip,slip_clinical,modified_slip,sens_slip,
           found_clinical,c_slip,
           caught,
           s_survival,c_survival,
           original_survivors,shifted_survivors,
           original_deaths,shifted_deaths)
  
  excluded_survival
}

#load_seer_dwell_model
load_seer_dwell_model<-function(origin_host_dir=here("data")){
  xStage<-c("I","II","III","IV")
  dwell_model_group_df<-read_tsv(sprintf("%s/20200728_dwell_time_groups.tsv",origin_host_dir))
  dwell_model_timing_df<-read_tsv(sprintf("%s/20200728_dwell_group_timing.tsv",origin_host_dir))
  
  dwell_model_all_df<-dwell_model_group_df %>% 
    rename(dwell_group=group) %>%
    full_join(dwell_model_timing_df %>% rename(dwell_group=group)) %>%
    full_join(dwell_model_timing_df) %>%
    mutate(number_stage=match(Stage,xStage))
  
  dwell_model_all_df
}

#standard parameter set
dwell_standard_model<-load_seer_dwell_model()

#load and clean manuscript sensitivity
iso_sens_joined<- iso_join_seer

#remove everything not staged for initial analysis
incidence_sens_source<-iso_sens_joined %>% 
  filter(Stage!="NotStaged") %>%
  mutate(iso_sens=sens)

#keep not staged for adding back
incidence_excluded_source<-iso_sens_joined %>%
  filter(Stage=="NotStaged")

##dwell time
#annual screening examined here for sensitivity
#generate exact slip rates
dwell_slip_rate<-exact_slip_rate_from_dwell(dwell_standard_model,screen_interval=1,weibull_shape=1,compliance)

#generate prevalent slip rate by clever use of very large interval and multiplying expectation
long_interval<-100
dwell_prevalent_rate<-exact_slip_rate_from_dwell(dwell_standard_model,screen_interval=long_interval,weibull_shape=1,compliance)

#no screening is happening, therefore nothing is ever intercepted
dwell_no_rate<-dwell_slip_rate %>% 
  mutate(slip=1.0,
         slip_clinical=1.0)

# accumulate 4 dwell scenarios
# prevalent and incident results
# plus perfect screening
# plus no screening
my_list<-vector("list",4*2+2)
k<-1
for (dw_scenario in 1:4){
  print(k)
  local_performance<-run_intercept_model(incidence_sens_source,
                                         dwell_slip_rate %>% 
                                           filter(scenario==dw_scenario))
  local_excluded<-run_excluded_model(incidence_excluded_source) #does not depend on scenario
  local_performance<-bind_rows(local_performance,local_excluded)
  
  incident_performance<-local_performance
  local_performance<-local_performance %>%
    mutate(screen_interval=1,
           dw_scenario=dw_scenario,
           scan="incident")
  
  my_list[[k]]<-local_performance
  k<-k+1
  
  #generate prevalent round starting off screening
  #only going to use caught by cfdna and combine with incident
  #because expected rates = average over all years, we can reverse identity
  #to obtain first-year screen by multiplying
  #rather than doing a special integral for prevalent rounds
  prevalent_performance<-run_intercept_model(incidence_sens_source,
                                             dwell_prevalent_rate %>% 
                                               filter(scenario==dw_scenario))
  
  prevalent_performance<-prevalent_performance %>% 
    filter(found_clinical==1) %>%
    mutate(caught=caught*long_interval,
           original_survivors=original_survivors*long_interval,
           shifted_survivors=shifted_survivors*long_interval,
           original_deaths=original_deaths*long_interval,
           shifted_deaths=shifted_deaths*long_interval) %>%
    bind_rows(incident_performance %>% 
                filter(found_clinical==2))
  
  prevalent_performance<-prevalent_performance %>%
    mutate(screen_interval=1,
           dw_scenario=dw_scenario,
           scan="prevalent")
  my_list[[k]]<-prevalent_performance
  k<-k+1
}

#this is the MIS scenario where schedule sensitivity is perfect so slip rates are 0
optimal_performance<-run_intercept_model(incidence_sens_source,
                                         dwell_slip_rate %>% 
                                           filter(scenario==0))
optimal_excluded<-run_excluded_model(incidence_excluded_source) #does not depend on scenario
optimal_performance<-bind_rows(optimal_performance,optimal_excluded)
optimal_performance<-optimal_performance %>% mutate(opt="0",dw_scenario=0,scan="incident")

#no screening so nothing found by cfdna operations
no_screening_performance<-run_intercept_model(incidence_sens_source,
                                              dwell_no_rate %>%
                                                filter(scenario==0))
no_screening_performance<-bind_rows(no_screening_performance, optimal_excluded) %>%
  mutate(opt="NO",dw_scenario=0,scan="no")

all_options_df<-bind_rows(my_list,.id="opt") %>%
  bind_rows(optimal_performance) %>%
  bind_rows(no_screening_performance)

#Now we have the full, detailed data frame
#add some helper text fields to clarify states represented by each line of the file 
text_options_df<-all_options_df %>%
  select(opt,Cancer,clinical,prequel,found_clinical,
         caught,s_survival,c_survival, 
         original_survivors,shifted_survivors,
         original_deaths,shifted_deaths,
         screen_interval,dw_scenario,scan) %>%
  mutate(mode_found=c("cfdna","soc")[found_clinical],
         aggressive=c("MIS","VSlow","Slow","Fast","AggFast")[dw_scenario+1])

#exact tracking of stage shift
#including flow per stage and resulting change in 5 year mortality

#pick my big color set


#translate from numbers to stage
cStage<-c("NS","I","II","III","IV")

#read in basic data set
rich_survival<-text_options_df

text_levels<-c("MIS","VSlow","Slow","Fast","AggFast")
cStage<-c("NS","I","II","III","IV")

##MIS: maximum interception scenario
scenario_code = "MIS"
#for each scenario looped over

data <- rich_survival %>%
  filter(scan%in%c("incident")&aggressive!=scenario_code|scan=="no")
intercept_shifted<-data %>% 
  filter(mode_found=="cfdna") %>%
  group_by(aggressive,Cancer,prequel) %>% 
  summarize(caught=sum(caught),.groups="keep") %>% 
  ungroup() %>%
  mutate(Stage=cStage[prequel+1])

intercept_original<-data %>% 
  filter(mode_found=="cfdna") %>%
  group_by(aggressive,Cancer,clinical) %>% 
  summarize(caught=sum(caught),.groups="keep") %>% 
  ungroup() %>%
  mutate(Stage=cStage[clinical+1])
stage_shift <-intercept_shifted %>% 
  select(aggressive,Cancer,Stage,cfdna=caught) %>%
  left_join(intercept_original %>% 
              select(aggressive,Cancer,Stage,clinical=caught)) %>%
  gather(key="case",value="caught",cfdna,clinical) %>%
  mutate(case=case_when(case=="clinical" ~ "pre-intercept",
                        TRUE ~ "intercepted")) %>%
  group_by(aggressive,case,Stage) %>%
  summarize(caught=sum(caught),.groups="keep") %>%
  ungroup() %>%
  mutate(Stage=factor(Stage,levels=rev(cStage)))
df=stage_shift%>%
  mutate(Stage=factor(Stage,levels=cStage),aggressive=ifelse(aggressive == "MIS", "no screen", aggressive),aggressive=factor(aggressive,levels=c("no screen","VSlow","Slow","Fast","AggFast")))

ggplot(df[df$case=="intercepted",], aes(x = aggressive, y = caught, fill = aggressive)) +
  geom_col(width = 0.9,alpha=0.7) +
  labs(title = "",y = "Incidence per 100K", x ="Stage") +
  #scale_y_continuous(limits = c(0.0, 1.1), breaks = seq(0.0, 1.1, 0.2), labels = c("0.0", "0.2", "0.4", "0.6", "0.8",  "1.0"), expand = c(0, 0)) +
  facet_wrap(~Stage,strip.position="bottom",nrow = 1) +#switch="x","y"
  theme_classic() +  # 使用经典主题，没有网格线
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.placement  = "outside", 
    panel.spacing.x = unit(0.2, "cm"), 
    axis.line = element_line(color = "black"),  # 设置坐标轴线的粗细
    axis.ticks.y  = element_line(color = "black"),
    axis.ticks.x  = element_blank(),
    axis.text.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.direction= "horizontal",
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.spacing = unit(0.1, "cm")
  )+
  scale_fill_manual(values = c("no screen"= "gray","VSlow"= "#006bb9","Slow"= "#CB98B2","Fast"= "#FFE799","AggFast"= "#E3353A"))
library(dplyr)

da = df%>%
  filter(aggressive != "no screen") %>%
  group_by(aggressive, case) %>%
  mutate(total = sum(caught)) %>%
  ungroup() %>%
  pivot_wider(
    id_cols = c(aggressive, Stage),
    names_from = case,
    values_from = caught
  ) %>%
  group_by(aggressive) %>%
  summarise(
    intercepted_sum = sum(`intercepted`[Stage %in% c("I", "II")], na.rm = TRUE),
    preintercept_sum = sum(`pre-intercept`[Stage %in% c("I", "II")], na.rm = TRUE),
    total_intercepted = sum(`intercepted`, na.rm = TRUE),
    shift_prob = (intercepted_sum - preintercept_sum) / total_intercepted
  ) %>%
  ungroup()
ggplot(da, aes(x = aggressive, y = shift_prob, fill = aggressive)) +
  geom_col(width = 0.7, alpha = 1) +
  # 修改数值标签为百分数格式
  geom_text(aes(label = scales::percent(shift_prob, accuracy = 0.01)),  # 显示为百分数，精度0.1%
            vjust = -0.5,  # 垂直位置调整
            size = 3,  # 字体大小
            color = "black") +  # 字体颜色
  labs(title = "", y = "Incidence per 100K", x = "") +
  # 将y轴刻度标签改为百分数格式
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # 精度1%
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.spacing.x = unit(0.08, "cm"), 
    axis.line = element_line(color = "black"),
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.spacing = unit(0.1, "cm")
  ) +
  scale_fill_manual(values = c("VSlow" = "#4472C4", 
                               "Slow" = "#A88BC1", 
                               "Fast" = "#F1B615", 
                               "AggFast" = "#DF5A52"))
#stage-shift
data <- rich_survival %>%
  filter(scan%in%c("incident"))#&aggressive!=scenario_code|scan=="no"


all_a_summary <- data %>%
  group_by(aggressive,mode_found,clinical,prequel) %>%
  summarize(IR=sum(caught),
            Deaths=sum(shifted_deaths),
            pre=sum(original_deaths),
            scan=scan[1],
            aggressive=aggressive[1],
            .groups="keep") %>%
  mutate(Delta=pre-Deaths,
         Survived=IR-Deaths-Delta) %>%
  ungroup()
found_code<-"cfdna"
all_b_summary<- all_a_summary %>% 
  filter(mode_found=="cfdna") #%>%&aggressive!="MIS" |mode_found=="soc"&aggressive=="MIS"
#mutate(aggressive=ifelse(aggressive=="MIS","no screen",aggressive))

da <- all_b_summary %>% 
  select(aggressive, clinical, prequel, pre,IR,Dead = Deaths, Alive = Survived, Saved = Delta) %>%
  mutate(
    clinical = cStage[clinical + 1],
    prequel = cStage[prequel + 1]
  ) %>%
  group_by(aggressive) %>% 
  mutate(
    total = sum(IR),
    survival = sum(Alive + Saved),pre_death=sum(pre),
    prob = survival / total
  ) %>%
  ungroup()
h=sum(da$Alive[da$aggressive=="MIS"])/sum(da$IR[da$aggressive=="MIS"])
da=da[,c("aggressive","prob")]
da <- unique(da)
da=da[!da$aggressive == "MIS", ]  
da=as.data.frame(tibble::add_row(da, aggressive = "no screen", prob = h))
da$aggressive=factor(da$aggressive,levels=c("no screen","VSlow","Slow","Fast","AggFast"))
ggplot(da, aes(x = aggressive, y = prob, fill = aggressive)) +
  geom_col(width = 0.7, alpha = 1) +
  # 添加数值标签（百分数格式）
  geom_text(aes(label = scales::percent(prob, accuracy = 0.01)),  # 精度0.1%
            vjust = -0.5,  # 垂直位置调整（在柱子上方）
            size = 3,  # 字体大小
            color = "black") +  # 字体颜色
  labs(title = "", y = "Incidence per 100K", x = "") +
  # 将y轴刻度标签改为百分数格式
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +  # 精度1%
  theme_classic() +  # 使用经典主题，没有网格线
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    panel.spacing.x = unit(0.08, "cm"), 
    axis.line = element_line(color = "black"),  # 设置坐标轴线的粗细
    axis.ticks.y = element_line(color = "black"),
    axis.ticks.x = element_line(color = "black"),
    panel.grid.major.y = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.spacing = unit(0.1, "cm")
  ) +
  scale_fill_manual(values = c("no screen" = "gray", 
                               "VSlow" = "#4472C4", 
                               "Slow" = "#A88BC1", 
                               "Fast" = "#F1B615", 
                               "AggFast" = "#DF5A52"))
a_intercept<- rich_survival %>%
  filter(scan=="incident",aggressive==scenario_code)

#stage-shift
intercept_shifted<-a_intercept %>% 
  filter(mode_found=="cfdna") %>%
  group_by(Cancer,prequel) %>% 
  summarize(caught=sum(caught),.groups="keep") %>% 
  ungroup() %>%
  mutate(Stage=cStage[prequel+1])

intercept_original<-a_intercept %>% 
  filter(mode_found=="cfdna") %>%
  group_by(Cancer,clinical) %>% 
  summarize(caught=sum(caught),.groups="keep") %>% 
  ungroup() %>%
  mutate(Stage=cStage[clinical+1])

#generate horizontal barplot
stage_shift_in_intercepted <- intercept_shifted %>% 
  select(Cancer, Stage, cfdna = caught) %>%
  left_join(intercept_original %>% 
              select(Cancer, Stage, clinical = caught)) %>%
  gather(key = "case", value = "caught", cfdna, clinical) %>%
  mutate(case = case_when(case == "clinical" ~ "pre-intercept",
                          TRUE ~ "intercepted")) %>%
  group_by(case, Stage) %>%
  summarize(caught = sum(caught), .groups = "keep") %>%
  ungroup() %>%
  mutate(Stage = factor(Stage, levels = cStage))  # 注意：这里去掉了rev()，因为现在是垂直显示

# 绘图部分调整
stage_shift1 <- ggplot(stage_shift_in_intercepted, aes(x = Stage, fill = case, y = caught)) +
  geom_bar(stat = "identity", position = position_dodge(0.85), width = 0.82, alpha = 1) + 
  geom_text(aes(label = round(caught, 2)), 
            position = position_dodge(0.85), 
            vjust = -0.5,  # 将hjust改为vjust，标签放在柱子上方
            color = "black", 
            size = 4.3) +
  scale_fill_manual(values = c("#4472C4", "#E36E67")) +
  scale_x_discrete(breaks = cStage, drop = TRUE) +  # 调整x轴刻度顺序
  scale_y_continuous(position = "left", 
                     expand = expansion(mult = c(0, 0.1))) +  # 为标签留出空间
  theme_bw() +
  theme(
    axis.text = element_text(size = 12.5),
    axis.text.x = element_text(angle = 45, hjust = 1),  # 添加x轴标签倾斜
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 15),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 12)
  ) +
  labs(
    y = "Number of diagnoses per 100K",
    x = "Stage"  # 添加x轴标签
  ) +
  ylim(0, max(stage_shift_in_intercepted$caught) * 1.1)  # 调整y轴范围

stage_shift1
#reduce to stage found, stage at, mode found
all_a_summary <- a_intercept %>%
  group_by(mode_found,clinical,prequel) %>%
  summarize(IR=sum(caught),
            Deaths=sum(shifted_deaths),
            Delta=sum(original_deaths),
            scan=scan[1],
            aggressive=aggressive[1],
            .groups="keep") %>%
  mutate(Delta=Delta-Deaths,
         Survived=IR-Deaths-Delta) %>%
  ungroup()
found_code<-"cfdna"
all_b_summary<- all_a_summary %>% 
  filter(mode_found==found_code)

ee_sankey <- all_b_summary %>% 
  select(clinical,prequel,Dead=Deaths,Alive=Survived,Saved=Delta) %>%
  mutate(clinical=cStage[clinical+1],prequel=cStage[prequel+1]) %>%
  pivot_longer(cols=c("Dead","Alive","Saved"),names_to="status",values_to="value") %>%
  mutate(uid=1:length(status)) %>%
  select(uid,StageInSEER=clinical,StageInScenario=prequel,FiveYear=status,value) %>%
  pivot_longer(cols=c("StageInSEER","StageInScenario","FiveYear"),names_to="event",values_to="Status") %>%
  mutate(event=factor(event,levels=c("StageInSEER","StageInScenario","FiveYear"))) %>%
  arrange(uid,event,Status)

#compute percentages to supply as label
ee_pct<-ee_sankey %>% 
  group_by(event,Status) %>% 
  summarize(tv=sum(value),.groups="keep") %>%
  ungroup() %>% 
  group_by(event) %>% 
  mutate(pct=round(100*tv/sum(tv),2)) %>%
  ungroup() %>% 
  mutate(text_pct=paste0(Status,"(",pct,"%)"))

#percentages added
saved_sankey_plot<- ee_sankey %>%
  filter(value>0) %>% 
  left_join(ee_pct %>% select(Status,event,text_pct,pct)) %>%
  mutate(Status=factor(Status,levels=c(cStage,"Alive","Saved","Dead"))) 
#pick my big color set
color_set_sankey<-c("NS"="#BEBFC0",
                    "I"="#87A4D9",
                    "II"="#F9DB96",
                    "III"="#BEBFC0",
                    "IV"="#E4A6AF",
                    "Alive"="#4472C4",
                    "Saved"="#F1B615",
                    "Dead"="#DF5A52")
a_plot = ggplot(saved_sankey_plot,aes(x=event,stratum=Status,alluvium=uid,y=value,fill=Status,label=text_pct))+
  scale_x_discrete(expand = c(0.05,0.1),position = "bottom",labels=c("Pre-intercept","Intercept","Outcome")) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0)),breaks=seq(0,31,height_ticks))+
  geom_flow(color="transparent",linewidth=0.2,alpha=1,width = 0.25)+
  geom_stratum(alpha=1,color="white",linewidth=0.5,width = 2/7)+#0.4
  geom_text(aes(label=Status),stat="stratum",size=3.4,vjust=-0.2)+
  geom_text(aes(label=paste0("(",pct,"%)")),stat="stratum",size=3.2,vjust=1.2)+
  scale_fill_manual(values=color_set_sankey)+
  guides(fill="none")+
  labs(x="",y="Individuals Diagnosed With GC")+theme_bw()+
  #coord_cartesian(ylim=c(0,31)) +
  theme(strip.text=element_text(size=18),
        axis.text.x = element_text(size=13),
        axis.text.y = element_blank(),
        axis.title=element_text(size=14),
        plot.title = element_text(size = 15),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())
a_plot

#load and clean manuscript sensitivity
##dwell time
#annual screening examined here for sensitivity
#generate exact slip rates
h=seq(0,1,0.1)[-1]
data_l=data.frame()
data_j=data.frame()
for (compliance in h) {
  dwell_slip_rate<-exact_slip_rate_from_dwell(dwell_standard_model,screen_interval=1,weibull_shape=1,compliance)
  
  #generate prevalent slip rate by clever use of very large interval and multiplying expectation
  long_interval<-100
  dwell_prevalent_rate<-exact_slip_rate_from_dwell(dwell_standard_model,screen_interval=long_interval,weibull_shape=1,compliance)
  
  #no screening is happening, therefore nothing is ever intercepted
  dwell_no_rate<-dwell_slip_rate %>% 
    mutate(slip=1.0,
           slip_clinical=1.0)
  optimal_performance<-run_intercept_model(incidence_sens_source,
                                           dwell_slip_rate %>% 
                                             filter(scenario==0))
  optimal_excluded<-run_excluded_model(incidence_excluded_source) #does not depend on scenario
  optimal_performance<-bind_rows(optimal_performance,optimal_excluded)
  optimal_performance<-optimal_performance %>% mutate(opt="0",dw_scenario=0,scan="incident")
  
  #no screening so nothing found by cfdna operations
  no_screening_performance<-run_intercept_model(incidence_sens_source,
                                                dwell_no_rate %>%
                                                  filter(scenario==0))
  no_screening_performance<-bind_rows(no_screening_performance, optimal_excluded) %>%
    mutate(opt="NO",dw_scenario=0,scan="no")
  all_options_df<-bind_rows(optimal_performance) %>%
    bind_rows(no_screening_performance)
  
  #Now we have the full, detailed data frame
  #add some helper text fields to clarify states represented by each line of the file 
  text_options_df<-all_options_df %>%
    select(opt,Cancer,clinical,prequel,found_clinical,
           caught,s_survival,c_survival, 
           original_survivors,shifted_survivors,
           original_deaths,shifted_deaths,dw_scenario,scan) %>%#screen_interval,
    mutate(mode_found=c("cfdna","soc")[found_clinical],
           aggressive=c("MIS"))
  #for each scenario looped over
  a_intercept<- text_options_df %>%
    filter(scan=="incident",aggressive=="MIS")
  
  #stage-shift
  intercept_shifted<-a_intercept %>% 
    filter(mode_found=="cfdna") %>%
    group_by(Cancer,prequel) %>% 
    summarize(caught=sum(caught),.groups="keep") %>% 
    ungroup() %>%
    mutate(Stage=cStage[prequel+1])
  
  intercept_original<-a_intercept %>% 
    filter(mode_found=="cfdna") %>%
    group_by(Cancer,clinical) %>% 
    summarize(caught=sum(caught),.groups="keep") %>% 
    ungroup() %>%
    mutate(Stage=cStage[clinical+1])
  
  #generate horizontal barplot.im
  stage_shift_in_intercepted <-intercept_shifted %>% 
    select(Cancer,Stage,cfdna=caught) %>%
    left_join(intercept_original %>% 
                select(Cancer,Stage,clinical=caught)) %>%
    gather(key="case",value="caught",cfdna,clinical) %>%
    mutate(case=case_when(case=="clinical" ~ "pre-intercept",
                          TRUE ~ "intercepted")) %>%
    group_by(case,Stage) %>%
    summarize(caught=sum(caught),.groups="keep") %>%
    ungroup() %>%
    mutate(Stage=factor(Stage,levels=rev(cStage)),compliance=compliance)
  shift_n=sum(stage_shift_in_intercepted$caught[1:2])-sum(stage_shift_in_intercepted$caught[5:6])
  all_a_summary <- a_intercept %>%
    group_by(mode_found,clinical,prequel) %>%
    summarize(IR=sum(caught),
              Deaths=sum(shifted_deaths),
              Delta=sum(original_deaths),
              scan=scan[1],
              aggressive=aggressive[1],
              .groups="keep") %>%
    mutate(Delta=Delta-Deaths,
           Survived=IR-Deaths-Delta) %>%
    ungroup()
  found_code<-"cfdna"
  all_b_summary<- all_a_summary %>% 
    filter(mode_found==found_code)
  ee_sankey <- all_b_summary %>% 
    select(clinical,prequel,Dead=Deaths,Alive=Survived,Saved=Delta) %>%
    mutate(clinical=cStage[clinical+1],prequel=cStage[prequel+1])
  save_n=sum(ee_sankey$Saved)
  data_l=rbind(data_l,cbind(saved_n=save_n,shift_n=shift_n,compliance=compliance))
}
library(scales)
data=data_l%>%
  mutate(label=percent(compliance,accuracy=1))
ggplot(data, aes(x = compliance, y = saved_n)) +
  geom_line(color="black") +
  geom_point(size = 3, shape = 21, fill = "red",stroke = 0.5) +
  theme_light(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.text.x = element_text(size = 12, hjust = 1),
    axis.title.x = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "black"),
    legend.position = c(0.1, 0.9),
    legend.background = element_rect(fill = alpha("white", 0.8)),
    # 添加以下两行来移除网格线
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  ) +
  labs(
    title = "",
    x = NULL,
    y = "Mortality Reductions per 100K"
  )
ggplot(data, aes(x = compliance, y = shift_n )) +
  geom_line(color="black") +#linewidth = 1.5,
  geom_point(size = 3, shape = 21, fill = "red",stroke = 0.5) +
  #scale_y_continuous(labels = label_number_si()) +  # 智能数字格式
  theme_light(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
    axis.text.y = element_text(size = 12, hjust = 1),  # Y轴标签右对齐
    axis.text.x = element_text(size = 12, hjust = 1),
    axis.title.x = element_text(size = 12),
    plot.subtitle = element_text(hjust = 0.5, color = "gray50"),
    legend.position = c(0.1, 0.9),  # 图例在图表内部
    legend.background = element_rect(fill = alpha("white", 0.8)),
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank()   # 移除次要网格线
  ) +
  labs(
    title = "",
    x = NULL,
    y = "Reduction cancers per 100K"
  )

