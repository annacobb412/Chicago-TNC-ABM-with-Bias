library(tidyverse)
library(dplyr)
library(colorBlindness)
library(latex2exp)
library(extrafont)
library(paletteer)
theme_set(theme_bw())
# Date: 07/09/24
# Code Purpose: Takes in "coordinate files" (processed simulation results files)
# for simulations which have been run for both a bias case and a no-bias case 
# and creates boxplots. Assumes all bias results files were for the same 
# fraction (bias_level) of drivers.

# Inputs to modify
list_of_dates = c("01_12_21","01_13_22")
bias_level = "0.023" 
file_dir = "Coordinate Files" # directory where coordinate files are stored
characteristic = "_most_updated" # leave "_" at front unless you have no characteristic
                                 # in your file name

# instantiating data frames
all_dates_avgs = data.frame(matrix(ncol = 4, nrow = 0))
colnames(all_dates_avgs) = c("bias_status","race","avg_WT","date")
all_dates_diffs = data.frame(matrix(ncol = 5,nrow = 0))
colnames(all_dates_diffs) = c("date","difference","trip_count","bias_diff","unbias_diff")
all_dates_meds = data.frame(matrix(ncol = 4, nrow = 0))
colnames(all_dates_meds) = c("bias_status","race","med_WT","date")
all_dates_diffs_m = data.frame(matrix(ncol = 5,nrow = 0))
colnames(all_dates_diffs_m) = c("date","difference","trip_count","bias_diff","unbias_diff")
all_dates_B = data.frame(matrix(ncol = 4, nrow = 0))
colnames(all_dates_B) = c("bias_status","race","WT_sec","date")
all_dates_UB = data.frame(matrix(ncol = 4,nrow = 0))
colnames(all_dates_UB) = c("bias_status","race","WT_sec","date")
all_dates_all_races = data.frame(matrix(ncol = 4,nrow = 0))
colnames(all_dates_all_races) = c("bias_status","race","WT_sec","date")
all_dates_diffs_race = data.frame(matrix(ncol = 8, nrow = 0))
all_coord = data.frame(matrix(ncol = 10, nrow = 0))

# filling in data frame
for (date in list_of_dates) {
  file_UB = paste0(file_dir,"/",date,characteristic,"_0.0_coord.csv")
  file_B = paste0(file_dir,"/",date,characteristic,"_B_",bias_level,"_coord.csv")
  UB_coord = read.csv(file_UB)
  B_coord = read.csv(file_B)
  if ("cancel" %in% colnames(B_coord)) {
    B_coord = B_coord %>% select(-cancel)
  }
  UB_coord_20 = UB_coord %>% filter(WT_sec <= 20*60)
  
  pc_UB_20 = round(100*nrow(UB_coord_20)/nrow(UB_coord),2)
  pc_UB_30 = round(100*nrow(UB_coord %>% filter(WT_sec <= 30*60))/nrow(UB_coord),2)
  pc_UB_60 = round(100*nrow(UB_coord %>% filter(WT_sec <= 60*60))/nrow(UB_coord),2)
  pc_UB_90 = round(100*nrow(UB_coord %>% filter(WT_sec <= 90*60))/nrow(UB_coord),2)
  pc_UB_120 = round(100*nrow(UB_coord %>% filter(WT_sec <= 120*60))/nrow(UB_coord),2)
  UB_coord_BW = UB_coord %>% 
    filter(race_word == "White" | race_word == "Black") %>% 
    rename(race = race_word) %>% 
    mutate(bias_status = "unbiased", date = date)
  B_coord_BW = B_coord %>% 
    filter(race_word == "White" | race_word == "Black") %>% 
    rename(race = race_word) %>% 
    mutate(bias_status = "biased", date = date)
  all_dates_UB = rbind(all_dates_UB,UB_coord_BW)
  intermediate = UB_coord %>% mutate(date = date)
  all_dates_all_races = rbind(all_dates_all_races,intermediate) 
  all_dates_B = rbind(all_dates_B,B_coord_BW)
  UB_avgs = UB_coord %>% 
    group_by(race_word) %>% 
    summarize(avg_WT = mean(WT_sec)) %>% 
    mutate(bias_status = "unbiased")
  B_avgs = B_coord %>% 
    group_by(race_word) %>% 
    summarize(avg_WT = mean(WT_sec)) %>% 
    mutate(bias_status = "biased")
  UB_meds = UB_coord %>% 
    group_by(race_word) %>% 
    summarize(med_WT = median(WT_sec)) %>% 
    mutate(bias_status = "unbiased")
  B_meds = B_coord %>% 
    group_by(race_word) %>% 
    summarize(med_WT = median(WT_sec)) %>% 
    mutate(bias_status = "biased")
  avgs = rbind(UB_avgs,B_avgs) %>% 
    mutate(date = date) %>% 
    rename(race = race_word)
  meds = rbind(UB_meds,B_meds) %>% 
    mutate(date = date) %>% 
    rename(race = race_word)
  all_dates_avgs = rbind(all_dates_avgs,avgs)
  all_dates_meds = rbind(all_dates_meds,meds)
  
  WT_B_B = (B_avgs %>% filter(race_word == "Black"))$avg_WT
  WT_W_B = (B_avgs %>% filter(race_word == "White"))$avg_WT
  WT_H_B = (B_avgs %>% filter(race_word == "Hispanic"))$avg_WT
  WT_A_B = (B_avgs %>% filter(race_word == "Asian"))$avg_WT
  WT_O_B = (B_avgs %>% filter(race_word == "Other"))$avg_WT
  WT_B_UB = (UB_avgs %>% filter(race_word == "Black"))$avg_WT
  WT_W_UB = (UB_avgs %>% filter(race_word == "White"))$avg_WT
  WT_H_UB = (UB_avgs %>% filter(race_word == "Hispanic"))$avg_WT
  WT_A_UB = (UB_avgs %>% filter(race_word == "Asian"))$avg_WT
  WT_O_UB = (UB_avgs %>% filter(race_word == "Other"))$avg_WT
  diff = (WT_B_B - WT_W_B) - (WT_B_UB - WT_W_UB)
  diff_row = data.frame(date = date, difference = diff, trip_count = nrow(UB_coord),
                        bias_diff = (WT_B_B - WT_W_B), unbias_diff = (WT_B_UB - WT_W_UB))
  all_dates_diffs = rbind(all_dates_diffs,diff_row)
  
  WT_B_B_m = (B_meds %>% filter(race_word == "Black"))$med_WT
  WT_W_B_m = (B_meds %>% filter(race_word == "White"))$med_WT
  WT_H_B_m = (B_meds %>% filter(race_word == "Hispanic"))$med_WT
  WT_A_B_m = (B_meds %>% filter(race_word == "Asian"))$med_WT
  WT_O_B_m = (B_meds %>% filter(race_word == "Other"))$med_WT
  WT_B_UB_m = (UB_meds %>% filter(race_word == "Black"))$med_WT
  WT_W_UB_m = (UB_meds %>% filter(race_word == "White"))$med_WT
  WT_H_UB_m = (UB_meds %>% filter(race_word == "Hispanic"))$med_WT
  WT_A_UB_m = (UB_meds %>% filter(race_word == "Asian"))$med_WT
  WT_O_UB_m = (UB_meds %>% filter(race_word == "Other"))$med_WT
  diff_m = (WT_B_B_m - WT_W_B_m) - (WT_B_UB_m - WT_W_UB_m)
  diff_row_m = data.frame(date = date, difference = diff_m, trip_count = nrow(UB_coord), 
                          bias_diff = (WT_B_B_m - WT_W_B_m), unbias_diff = (WT_B_UB_m - WT_W_UB_m))
  all_dates_diffs_m = rbind(all_dates_diffs_m,diff_row_m)
  
  diff_row_W = data.frame(date = date, race = "White", bias_WT_avg = WT_W_B, unbias_WT_avg = WT_W_UB, 
                          diff_avg = (WT_W_B - WT_W_UB), bias_WT_med = WT_W_B_m, unbias_WT_med = WT_W_UB_m,
                          diff_med = (WT_W_B_m - WT_W_UB_m))
  diff_row_B = data.frame(date = date, race = "Black", bias_WT_avg = WT_B_B, unbias_WT_avg = WT_B_UB, 
                          diff_avg = (WT_B_B - WT_B_UB), bias_WT_med = WT_B_B_m, unbias_WT_med = WT_B_UB_m,
                          diff_med = (WT_B_B_m - WT_B_UB_m))
  diff_row_H = data.frame(date = date, race = "Hispanic", bias_WT_avg = WT_H_B, unbias_WT_avg = WT_H_UB, 
                          diff_avg = (WT_H_B - WT_H_UB), bias_WT_med = WT_H_B_m, unbias_WT_med = WT_H_UB_m,
                          diff_med = (WT_H_B_m - WT_H_UB_m))
  diff_row_A = data.frame(date = date, race = "Asian", bias_WT_avg = WT_A_B, unbias_WT_avg = WT_A_UB, 
                          diff_avg = (WT_A_B - WT_A_UB), bias_WT_med = WT_A_B_m, unbias_WT_med = WT_A_UB_m,
                          diff_med = (WT_A_B_m - WT_A_UB_m))
  diff_row_O = data.frame(date = date, race = "Other", bias_WT_avg = WT_O_B, unbias_WT_avg = WT_O_UB, 
                          diff_avg = (WT_O_B - WT_O_UB), bias_WT_med = WT_O_B_m, unbias_WT_med = WT_O_UB_m,
                          diff_med = (WT_O_B_m - WT_O_UB_m))
  all_dates_diffs_race = rbind(all_dates_diffs_race,diff_row_W,diff_row_B,diff_row_H,diff_row_A,diff_row_O)
  all_coord = all_coord %>% rbind(UB_coord_BW,B_coord_BW)
  print(paste(date,"biased:",nrow(B_coord),"trips"))
  print(paste(date,"unbiased:",nrow(UB_coord),"trips"))
}

# PLOTTING UB BLACK VS. WHITE WT TIMES AGGREGATED ACROSS DATES (ONE BOXPLOT PER RACE)-----
UB_coord_BW = UB_coord_BW %>% 
  mutate(WT_min = round(WT_sec/60,2))
ggplot(data = UB_coord_BW) +
  geom_boxplot(aes(x = race,y = WT_min, fill = race),outlier.shape = NA) +
  coord_cartesian(y = c(0,18)) +
  ggtitle("Black and White Wait Times Aggregated \nAcross Simulated Dates with no Biased Drivers") +
  ylab("Wait Time [min]") +
  xlab("Race") +
  theme(legend.position = "none") +
  scale_fill_manual(values = c(PairedColor12Steps[8],PairedColor12Steps[7])) +
  theme(axis.text=element_text(size=11)) +
  theme(text = element_text(family = "Times New Roman"),legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), axis.title = element_text(size = 13)) 

# PLOTTING B & UB BLACK VS. WHITE WT TIMES AGGREGATED ACROSS DATES (2 BOXPLOTS PER RACE)----
all_coord = all_coord %>% 
  mutate(WT_min = round(WT_sec/60,2))
plot_code_order = c("Black (no bias)","Black (bias)","White (no bias)", "White (bias)")
coord_BW = all_coord %>% 
  mutate(plot_code = case_when((race == "White" & bias_status == "biased") ~ "White (bias)",
                               (race == "White" & bias_status == "unbiased") ~ "White (no bias)",
                               (race == "Black" & bias_status == "biased") ~ "Black (bias)",
                               (race == "Black" & bias_status == "unbiased") ~ "Black (no bias)"))
means_BW = coord_BW %>% 
  group_by(plot_code) %>% 
  summarize(avg_WT = mean(WT_min), sd = sd(WT_min)) %>% 
  cbind(c("Black","Black","White","White"))
means_BW_20 = coord_BW %>% 
  filter(WT_min <= 20) %>% 
  group_by(plot_code) %>% 
  summarize(avg_WT = mean(WT_min), sd = sd(WT_min)) %>% 
  cbind(c("Black","Black","White","White"))
colnames(means_BW) = c("plot_code","avg_WT","sd","race")
colnames(means_BW_20) = c("plot_code","avg_WT","sd","race")
dodge = position_dodge(width = .75)
coord_BW_20 = coord_BW %>% filter(WT_min <= 20)

ggplot(data = coord_BW) +
  geom_boxplot(aes(x = race,y = WT_min, fill = factor(plot_code, level = plot_code_order)),outlier.shape = NA) +
  coord_cartesian(y = c(0,16)) +
  geom_point(data = means_BW, aes(x = factor(race, level = c("Black","White")), y = avg_WT, fill = factor(plot_code, level = plot_code_order)), position = dodge, shape = "square", size = 1) +
  ggtitle("Black and White Wait Times Aggregated \nAcross Simulated Dates With and Without Biased Drivers") +
  ylab("Wait Time [min]") +
  xlab("Race") +
  labs(fill = "Race (bias status)") +
  scale_fill_manual(values = c(PairedColor12Steps[8],PairedColor12Steps[10],PairedColor12Steps[7],PairedColor12Steps[9])) +
  theme(axis.text=element_text(size=11)) +
  theme(text = element_text(family = "Times New Roman"),legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"), axis.title = element_text(size = 13)) 

# -----------------------------PLOTTING IMPACT----------------------------------

all_dates_diffs = all_dates_diffs %>% 
  mutate(date_plot = gsub("_","/",date))

ggplot(data = all_dates_diffs) +
  geom_col(aes(x = factor(date_plot, level = date_order),y = difference), fill = PairedColor12Steps[10]) +
  ylab(TeX("$I = \\bar{\\Delta}_B - \\bar{\\Delta}_{UB}$ [seconds]")) +
  xlab("Simulation Date") +
  ggtitle("Impact of Direct Discrimination on Wait Times") +
  theme_bw() +
  ylim(c(-7.5,30)) +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), axis.title = element_text(size = 14),
        axis.text=element_text(size=12))+
  theme(axis.text.x = element_text(angle = 35, vjust = 0.95, hjust=1))

# ---------------------PLOTTING MEDIAN VERSION OF IMPACT------------------------

all_dates_diffs_m = all_dates_diffs_m %>% 
  mutate(date_plot = gsub("_","/",date))
ggplot(data = all_dates_diffs_m) +
  geom_col(aes(x = date_plot,y = difference), fill = PairedColor12Steps[10]) +
  ylab(TeX("$I = \\Delta_B - \\Delta_{UB}$ [seconds]")) +
  xlab("Simulation Date") +
  ggtitle("Impact (Median) of Direct Discrimination on Wait Times") +
  theme(axis.text=element_text(size=11)) +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 13)) 

# ------------------------PLOTTING AVERAGE WAIT TIMES---------------------------
all_dates_avgs_WB = all_dates_avgs %>% 
  filter(race == "White" | race == "Black") %>% 
  mutate(date_plot = gsub("_","/",date)) %>% 
  mutate(plot_code = case_when((race == "White" & bias_status == "biased") ~ "W | Bias ",
                               (race == "White" & bias_status == "unbiased") ~ "W | No Bias",
                               (race == "Black" & bias_status == "biased") ~ "B | Bias",
                               (race == "Black" & bias_status == "unbiased") ~ "B | No Bias"))
ggplot(data = all_dates_avgs_WB) +
  geom_col(aes(x = date_plot,y = avg_WT, fill = plot_code),position = "dodge") +
  xlab("Simulation Date") +
  scale_fill_manual(values = PairedColor12Steps[7:10]) +
  ylab("Average Wait Time [seconds]") +
  ggtitle("Wait Times across Race with and without Bias") +
  theme(axis.text=element_text(size=11)) +
  labs(fill = "Race | Biased Drivers? ") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 13)) 

# ---PLOTTING DISTRIBUTION OF WAIT TIMES (W/ & W/O BIAS) FOR ALL DATES----------
all_dates_avgs_WB = all_dates_avgs %>% 
  filter(race == "White" | race == "Black") %>% 
  mutate(date_plot = gsub("_","/",date), avg_WT_min = round(avg_WT/60,2)) %>% 
  mutate(plot_code = case_when((race == "White" & bias_status == "biased") ~ "White (3.2%)",
                               (race == "White" & bias_status == "unbiased") ~ "White (0.0%)",
                               (race == "Black" & bias_status == "biased") ~ "Black (3.2%)",
                               (race == "Black" & bias_status == "unbiased") ~ "Black (0.0%)"))
all_dates_B_UB = rbind(all_dates_B,all_dates_UB) %>% 
  mutate(date_plot = gsub("_","/",date),WT_min = round(WT_sec/60,2)) %>% 
  mutate(plot_code = case_when((race == "White" & bias_status == "biased") ~ "White (3.2%)",
                               (race == "White" & bias_status == "unbiased") ~ "White (0.0%)",
                               (race == "Black" & bias_status == "biased") ~ "Black (3.2%)",
                               (race == "Black" & bias_status == "unbiased") ~ "Black (0.0%)"))

dodge = position_dodge(width=0.75) 
ggplot(data = all_dates_B_UB) +
  geom_boxplot(aes(x = factor(date_plot, level = date_order), y = WT_min, fill = plot_code), outlier.shape = NA) +
  geom_point(data = all_dates_avgs_WB, aes(x = factor(date_plot, level = date_order), y = avg_WT_min, fill = plot_code), position = dodge, shape = "square", size = 1) +
  xlab("Simulation Date") +
  scale_fill_manual(values = PairedColor12Steps[7:10]) +
  coord_cartesian(y = c(0,20)) +
  ylab("Wait Times [minutes]") +
  ggtitle("Wait Times across Race with and without Bias") +
  theme(axis.text=element_text(size=11)) +
  labs(fill = "Race (% Bias)") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 13)) 

# ---------PLOTTING NO-BIAS DISTRIBUTION FOR ALL RACES FOR ALL DATES------------
all_races_plot = all_dates_all_races %>% 
  mutate(date_plot = gsub("_","/",date),WT_min = round(WT_sec/60,2)) 
all_dates_avgs_WHB = all_dates_avgs %>% 
  mutate(date_plot = gsub("_","/",date), avg_WT_min = round(avg_WT/60,2)) %>% 
  filter(bias_status == "unbiased")
dodge = position_dodge(width=0.75) 
race_order = c("Black","Hispanic","Other","Asian","White")
ggplot(data = all_races_plot) +
  geom_boxplot(aes(x = factor(date_plot, level = date_order), y = WT_min, fill = factor(race_word, level = race_order)), outlier.shape = NA) +
  geom_point(data = all_dates_avgs_WHB, aes(x = factor(date_plot, level = date_order), y = avg_WT_min, fill = factor(race, level = race_order)), position = dodge, shape = "square", size = 1) +
  xlab("Simulation Date") +
  scale_fill_manual(values = c(PairedColor12Steps[8],
                               PairedColor12Steps[9],
                               PairedColor12Steps[10],
                               PairedColor12Steps[11],
                               PairedColor12Steps[7]
                               )) +
  coord_cartesian(y = c(0,20)) +
  ylab("Wait Times [minutes]") +
  ggtitle("Wait Times across All Races with No Biased Drivers") +
  theme_bw() +
  labs(fill = "Race") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 13)) +
  theme(axis.text=element_text(size=13))
all_dates_avgs_UB = all_dates_avgs %>% filter(bias_status == "unbiased") %>% 
  mutate(avg_WT = round(avg_WT,2), date = gsub("_","/",date)) %>% 
  pivot_wider(names_from = "race", values_from = avg_WT) %>% select(-bias_status) %>% 
  mutate(Diff_B_W = Black - White, Diff_H_W = Hispanic - White)
write.csv(all_dates_avgs_UB,"all_dates_avgs_UB.csv")

# ---------------PLOTTING DISTRIBUTION OF WAIT TIMES (ONLY W/O BIAS)------------
all_dates_UB = all_dates_UB %>% 
  mutate(date_plot = gsub("_","/",date)) 
all_dates_avgs_WB_UB = all_dates_avgs %>% 
  filter(race == "White" | race == "Black") %>% 
  filter(bias_status == "unbiased") %>% 
  mutate(date_plot = gsub("_","/",date)) 
  
dodge = position_dodge(width=0.75) 
ggplot(data = all_dates_UB) +
  geom_boxplot(aes(x = factor(date_plot, level = date_order), y = WT_sec/60, fill = race), outlier.shape = NA) +
  geom_point(data = all_dates_avgs_WB_UB, aes(x = factor(date_plot, level = date_order), y = avg_WT/60, fill = race), position = dodge, shape = "square", size = 1) +
  xlab("Simulation Date") +
  coord_cartesian(y = c(0,20)) +
  scale_fill_manual(values = c(SteppedSequential5Steps[12],SteppedSequential5Steps[14])) +
  ylab("Wait Times [minutes]") +
  ggtitle("Wait Times of Black and White Riders\nwith No Biased Drivers") +
  theme_bw() +
  theme(axis.text=element_text(size=13)) +
  labs(fill = "Race") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 14)) 

# ------BLACK & WHITE RIDER WAIT TIMES & TRIP COUNT THROUGHOUT THE DAY----------
all_dates_avg_hourly = all_dates_UB %>% 
  mutate(to_convert = paste(start_time,start_date)) %>% 
  mutate(real_time = as.POSIXct(to_convert,format = "%H:%M:%OS %Y-%m-%d")) %>%
  mutate(just_hour = as.numeric(hour(real_time))) %>% 
  mutate(real_date = paste(substr(start_date,6,7),"_",substr(start_date,9,10),"_",substr(start_date,3,4),sep="")) %>% 
  filter(real_date %in% list_of_dates) %>% 
  group_by(race,just_hour,date) %>% 
  summarize(avg_WT = mean(WT_sec),num_trips = n())
all_dates_avg_hr_BW = all_dates_avg_hourly %>% 
  filter(race == "White" | race == "Black") %>% 
  mutate(plot_code = paste(race,": ",date,sep=""))

ggplot(data = all_dates_avg_hr_BW, aes(x = just_hour, y = avg_WT, color = plot_code)) +
  geom_line() +
  geom_point() +
  ylab("Average Wait Time [sec]") +
  xlab("Hour of the Day") +
  ggtitle("Average Wait Time Throughout the Day \n by Race and Date") +
  scale_color_manual(values = paletteer_d("ggthemes::Hue_Circle")) +
  theme(axis.text=element_text(size=12)) +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 13)) 
ggplot(data = all_dates_avg_hr_BW, aes(x = just_hour, y = num_trips, color = plot_code)) +
  geom_line() +
  geom_point() +
  ylab("Number of Trips") +
  xlab("Hour of the Day") +
  ggtitle("Number of Trips Throughout the Day \nby Race & Date") +
  theme_bw() +
  scale_color_manual(values = paletteer_d("ggthemes::Hue_Circle")) +
  theme(axis.text=element_text(size=12)) +
  theme(text = element_text(family = "Times New Roman")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), axis.title = element_text(size = 13)) 
