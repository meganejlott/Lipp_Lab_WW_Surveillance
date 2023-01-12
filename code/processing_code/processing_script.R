#Weekly Data Update Checklist

#1. Update cfx_n1.csv, cfx_n2.csv, recovery_data.csv, and calfguard.csv with qPCR data.
#2. Update plant_data.csv with influent flow and TSS data from ACC WRFs (Check Slack).
#3. Update the sample_data dataframe in this script with the most recent collection numbers. 
#4. Run this Script. 
#5. Send output XXX.rds to Billy for weekly website update.
#6. Send output XX.csv to Cristina (NWSS).


#Load Libraries
library(tidyverse)
library(plyr)
library(dplyr)
#install.packages("downloader")
library(downloader)


#Load Data 

##Generate Sample Collection Data Set
###Sample Collection + Sample Date
###This line needs to get updated week. 
####Generate Sampling Data Set
#Sampling Data 
####Generate Sampling Data Set
#Sampling Data 
sample_data = data.frame("collection_num" = 5:243, "date" = c("2020-06-16", "2020-06-23", "2020-06-30", "2020-07-07", "2020-07-14", "2020-07-21", "2020-07-28", 
                                                              "2020-08-04", "2020-08-11", "2020-08-18", "2020-08-25", "2020-09-01", 
                                                              "2020-09-08", "2020-09-15", "2020-09-22", "2020-09-29", 
                                                              "2020-10-06", "2020-10-13", "2020-10-20", "2020-10-27",
                                                              "2020-11-02", "2020-11-04", "2020-11-09", "2020-11-11", 
                                                              "2020-11-16", "2020-11-18", "2020-11-23", "2020-11-25", 
                                                              "2020-11-30", "2020-12-02", "2020-12-07", "2020-12-09", 
                                                              "2020-12-14", "2020-12-16", "2020-12-21", "2020-12-23", 
                                                              "2020-12-28", "2021-01-04", "2021-01-11", "2021-01-13", 
                                                              "2021-01-19", "2021-01-20", "2021-01-25", "2021-01-27", 
                                                              "2021-02-01", "2021-02-03", "2021-02-08", "2021-02-10", 
                                                              "2021-02-15", "2021-02-17", "2021-02-22", "2021-02-24", 
                                                              "2021-03-01", "2021-03-03", "2021-03-08", "2021-03-10", 
                                                              "2021-03-15", "2021-03-17", "2021-03-22", "2021-03-24", 
                                                              "2021-03-29", "2021-03-31", "2021-04-05", "2021-04-07", 
                                                              "2021-04-12", "2021-04-14", "2021-04-19", "2021-04-21", 
                                                              "2021-04-26", "2021-04-28", "2021-05-03", "2021-05-05", 
                                                              "2021-05-10", "2021-05-12", "2021-05-17", "2021-05-19", 
                                                              "2021-05-24", "2021-05-26", "2021-06-01", "2021-06-02", 
                                                              "2021-06-07", "2021-06-09", "2021-06-14", "2021-06-16", 
                                                              "2021-06-21", "2021-06-23", "2021-06-28", "2021-06-30", 
                                                              "2021-07-06", "2021-07-06", "2021-07-12", "2021-07-14", 
                                                              "2021-07-19", "2021-07-21", "2021-07-26", "2021-07-28", 
                                                              "2021-08-02", "2021-08-04", "2021-08-09", "2021-08-11", 
                                                              "2021-08-16", "2021-08-18", "2021-08-23", "2021-08-25", 
                                                              "2021-08-30", "2021-09-01", "2021-09-07", "2021-09-08", 
                                                              "2021-09-13", "2021-09-15", "2021-09-20", "2021-09-22", 
                                                              "2021-09-27", "2021-09-29", "2021-10-04", "2021-10-06", 
                                                              "2021-10-11", "2021-10-13", "2021-10-18", "2021-10-20", 
                                                              "2021-10-25", "2021-10-27", "2021-11-01", "2021-11-03", 
                                                              "2021-11-08", "2021-11-10", "2021-11-15", "2021-11-17", 
                                                              "2021-11-22", "2021-11-29", "2021-12-01", "2021-12-06", 
                                                              "2021-12-08", "2021-12-13", "2021-12-15", "2021-12-20", 
                                                              "2021-12-27", "2022-01-03", "2022-01-05", "2022-01-10", 
                                                              "2022-01-12", "2022-01-17", "2022-01-19", "2022-01-24", 
                                                              "2022-01-26", "2022-01-31", "2022-02-02", "2022-02-07", 
                                                              "2022-02-09", "2022-02-14", "2022-02-16", "2022-02-21", 
                                                              "2022-02-23", "2022-02-28", "2022-03-01", "2022-03-07", 
                                                              "2022-03-09", "2022-03-14", "2022-03-16", "2022-03-21", 
                                                              "2022-03-23", "2022-03-28", "2022-03-30", "2022-04-04", 
                                                              "2022-04-06", "2022-04-11", "2022-04-13", "2022-04-18", 
                                                              "2022-04-20", "2022-04-25", "2022-04-27", "2022-05-02", 
                                                              "2022-05-04", "2022-05-09", "2022-05-11", "2022-05-16", 
                                                              "2022-05-18", "2022-05-23", "2022-05-25", "2022-05-31", 
                                                              "2022-06-01", "2022-06-06", "2022-06-08", "2022-06-13", 
                                                              "2022-06-15", "2022-06-20", "2022-06-22", "2022-06-27", 
                                                              "2022-06-29", "2022-07-04", "2022-07-06", "2022-07-11", 
                                                              "2022-07-13", "2022-07-18", "2022-07-20", "2022-07-25", 
                                                              "2022-07-27", "2022-08-01", "2022-08-03", "2022-08-08", 
                                                              "2022-08-10", "2022-08-15", "2022-08-17", "2022-08-22", 
                                                              "2022-08-24", "2022-08-29", "2022-08-31", "2022-09-05", 
                                                              "2022-09-07", "2022-09-13", "2022-09-14", "2022-09-19", 
                                                              "2022-09-21", "2022-09-26", "2022-09-28", "2022-10-03", 
                                                              "2022-10-05", "2022-10-10", "2022-10-12", "2022-10-17", 
                                                              "2022-10-19", "2022-10-24", "2022-10-26", "2022-10-31", 
                                                              "2022-11-02", "2022-11-07", "2022-11-09", "2022-11-14", 
                                                              "2022-11-16", "2022-11-21", "2022-11-28", "2022-11-30",
                                                              "2022-12-05", "2022-12-07", "2022-12-12", "2022-12-14",
                                                              "2022-12-19", "2022-12-21", "2023-01-04"), stringsAsFactors = FALSE)
sample_data$date = as.Date(sample_data$date)
sample_data$collection_num = as.numeric(sample_data$collection_num)


##Load StepOne Data (Year 1)
stepone_n1 = read.csv("./data/raw_data/stepone_n1.csv")
stepone_n2 = read.csv("./data/raw_data/stepone_n2.csv")

##Load CFX Data (Year 2)
cfx_n1 = read_csv("./data/raw_data/cfx_n1.csv")
cfx_n2 = read_csv("./data/raw_data/cfx_n2.csv")

##Load WRF Plant Data
plant_data = read_csv("./data/raw_data/plant_data.csv")
plant_data$date = as.Date(plant_data$date, "%m/%d/%Y")

##Load BCoV Recovery Data

#Load Data 
recovery_data = read_csv("./data/raw_data/recovery_data.csv")
calf_guard = read_csv("./data/raw_data/calfguard.csv")



##Generate County Data Set
###Download Case Data from the GA DPH website
download("https://ga-covid19.ondemand.sas.com/docs/ga_covid_data.zip", dest="dataset.zip", mode="wb") 

###Unzip the file
unzip("dataset.zip", exdir = "./data/raw_data/ga_covid_data") 

###Filter for Athens-Clarke County & Select only relevant columns
case_data = read.csv("./data/raw_data/ga_covid_data/epicurve_rpt_date.csv") %>% 
  filter(county == "Clarke") %>% 
  select(report_date, cases, cases_cum, moving_avg_cases) 

###Format for downstream script
case_data$report_date = as.Date(case_data$report_date, "%Y-%m-%d")
case_data = case_data %>% filter(report_date > "2020-3-14")
case_data$moving_avg_cases = round(case_data$moving_avg_cases, digits = 0)
names(case_data) = c("date", "new_cases_clarke", "cases_cum_clarke", "X7_day_ave_clarke")



#Process the Data

##StepOne Data (Year 1)
stepone_n1$ct = as.numeric(stepone_n1$ct)
###Use the slope and intercept from standard curve, here
stepone_n1 = stepone_n1 %>% mutate(copy_num_uL_rxn = as.numeric(10^((ct-34.008)/-3.389)))
###use the LOD value here
stepone_n1$copy_num_uL_rxn = stepone_n1$copy_num_uL_rxn %>% replace_na(0.0002)


stepone_n2$ct = as.numeric(stepone_n2$ct)
###Use the slope and intercept from standard curve, here
stepone_n2 = stepone_n2 %>% mutate(copy_num_uL_rxn = as.numeric(10^((ct-32.416)/-3.3084)))
###use the LOD value here
stepone_n2$copy_num_uL_rxn = stepone_n2$copy_num_uL_rxn %>% replace_na(0.0002)


###Bind N1 and N2 data together
stepone_n1_n2 = bind_rows(stepone_n1, stepone_n2)
###Average together the qPCR replicates 
stepone_n1_n2_ave = plyr::ddply(stepone_n1_n2,.(sample_id, target, run_num),plyr::summarize, copy_num_uL_rxn = mean(copy_num_uL_rxn)) 
###Now, let's estimate copy num per L, based on the replicates
stepone_n1_n2_ave = mutate(stepone_n1_n2_ave, copy_num_L = copy_num_uL_rxn *20/2*25/3*60/280*1000*1000)



##CFX Data (Year 2)
cfx_n1$ct = as.numeric(cfx_n1$ct)
###Use the slope and intercept from standard curve, here
cfx_n1 = cfx_n1 %>% mutate(copy_num_uL_rxn = as.numeric(10^((ct-36.046)/-3.5293)))
###use the LOD value here to replace NAs
cfx_n1$copy_num_uL_rxn = cfx_n1$copy_num_uL_rxn %>% replace_na(0.004)

cfx_n2$ct = as.numeric(cfx_n2$ct)
###Use the slope and intercept from standard curve, here
cfx_n2 = cfx_n2 %>% mutate(copy_num_uL_rxn = as.numeric(10^((ct-37.731)/-3.2505)))
###use the LOD value here to replace NAs
cfx_n2$copy_num_uL_rxn = cfx_n2$copy_num_uL_rxn %>% replace_na(0.004)

###Bind N1 and N2 data together
cfx_n1_n2 = bind_rows(cfx_n1, cfx_n2)
###Average together the qPCR replicates
cfx_n1_n2_ave = plyr::ddply(cfx_n1_n2,.(sample_id, target, run_num),plyr::summarize, copy_num_uL_rxn = mean(copy_num_uL_rxn)) 
###Estimate copy num per L, based on the replicates
cfx_n1_n2_ave = mutate(cfx_n1_n2_ave, copy_num_L = copy_num_uL_rxn *20/5*60/280*1000*1000)

##Bind StepOne (Year 1) and CFX (Year 2) data together
n1_n2_ave = rbind(stepone_n1_n2_ave, cfx_n1_n2_ave)


##Separate out the WRF, sample week, and Rep ID. 
n1_n2_ave = n1_n2_ave %>% separate(col = sample_id, into = c("wrf","collection_num", "rep_id"), sep = "_") %>% 
  mutate(collection_num = as.numeric(collection_num))

#Add in the sample collection data
n1_n2_ave = left_join(n1_n2_ave, sample_data, by = ("collection_num"))




##Calculate total copies

###First, calculate influent flow in L 
####Conversion 1 gallon = 3.78541 L 
plant_data = mutate(plant_data, influent_flow_L = influent_flow_mg *1000000*3.78541)

#Join the data frames
n1_n2_plant = dplyr::left_join(n1_n2_ave, plant_data, by = c("date", "wrf"))

#Make a new column, where you calculate the total number of copies of the target per day.
n1_n2_plant = n1_n2_plant %>% mutate(total_copies = copy_num_L * influent_flow_L)


#Now, average the extraction replicates 
n1_n2_cleaned = plyr::ddply(n1_n2_plant, c("wrf", "collection_num","date", "target"), 
                            summarize, mean_copy_num_uL_rxn = mean(copy_num_uL_rxn), 
                            mean_copy_num_L = mean(copy_num_L), sd_L = sd(copy_num_L), 
                            se_L = sd(copy_num_L)/sqrt((length(copy_num_L))), mean_total_copies = mean(total_copies), 
                            sd_total_copies = sd(total_copies), lo_95 = mean(copy_num_L) - 2*se_L, 
                            up_95 = mean(copy_num_L) + 2*se_L)


#Add in case data
n1_n2_cleaned_cases = left_join(case_data, n1_n2_cleaned, by = c("date")) 


#Quick look at the viral load
ave_total_copies = plyr::ddply(n1_n2_cleaned, c("date", "wrf"), summarize, total_copies = mean(mean_total_copies))
ave_total_copies = plyr::ddply(ave_total_copies, c("date"), summarize, total_copies = sum(total_copies))
ave_total_copies %>% ggplot(aes(x = date, y = log10(total_copies))) + geom_point()
ggsave("./data/processed_data/Rplot.png")




#Recovery Data

recovery_data$ct = as.numeric(recovery_data$ct)
recovery_data = plyr::ddply(recovery_data,.(sample_id, cg_num),plyr::summarize, avg_ct = mean(ct))
recovery_data = recovery_data %>% drop_na() 

recovery_data = recovery_data %>% separate(sample_id, into = c("wrf", "collection_num", "rep"), sep = "_") %>% drop_na()
recovery_data$collection_num = as.numeric(recovery_data$collection_num)

recovery_data = recovery_data %>% mutate(copies_ul_rxn = 10^((avg_ct-30.7)/-3.238))
recovery_data = recovery_data %>% mutate(copies_ul_sample = copies_ul_rxn *20/2*25/3*60/280)

calf_guard = plyr::ddply(calf_guard,.(sample_id),plyr::summarize, avg_ct = mean(ct))
calf_guard = calf_guard %>% separate(sample_id, into = c("CG", "cg_num"), sep = "_") %>% drop_na()
calf_guard$cg_num = as.numeric(calf_guard$cg_num)
calf_guard = calf_guard %>% drop_na()

calf_guard = calf_guard %>% mutate(copies_ul_rxn = 10^((avg_ct-30.7)/-3.238))
calf_guard = calf_guard %>% mutate(copies_ul_input = copies_ul_rxn *20/2*25/3*60/50*40/40000)


output = recovery_data %>% select(wrf, collection_num, cg_num, copies_ul_sample, avg_ct)
input = calf_guard %>% select(cg_num, copies_ul_input, avg_ct)

recovery_calc = left_join(output, input, by = "cg_num")
recovery_calc = recovery_calc %>% mutate(percent_recovery = 100*(copies_ul_sample/copies_ul_input))



#Write data frames into "processed_data" folder
saveRDS(n1_n2_cleaned_cases, "./data/processed_data/n1_n2_cleaned_cases.RDS") #Send this to Slack (for Billy)
write.csv(n1_n2_cleaned_cases, "./data/processed_data/n1_n2_cleaned_cases.csv") #Send this to Slack (for Erin) 
write.csv(recovery_calc, "./data/processed_data/bcov_recovery.csv") #For NWSS Update


  