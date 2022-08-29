# Lipp Lab Weekly Waster Surveillance Updates

About:
This Repository includes the processing code to translate the raw RT-qPCR data into data ready for reporting on our dashboard and through NWSS. 

# Steps: 

Update Raw Data Files:
1. Go to the raw_data folder (data > raw_data)
2. Update cfx_n1.csv, cfx_n2.csv, recovery_data.csv, and calfguard.csv with qPCR output
3. Update plant_data.csv with influent flow and TSS data from ACC WRFs (Check Slack)

Run Processing Script
1. Open ww_monitoring R Project
2. Open processing_script.R in the processing_code folder (code > pcocessing_code) 
3. Within the script, be sure to update the sample_data dataframe with the most recent collection numbers and dates. 
4. Run the Script. 

Share the Data
1. Processed data will be written into the processed_data folder (data > processed_data) 
2. In Slack, post the weekly viral load chart (copy and paste from Plot Viewer in R). 
3. In Slack, post BOTH the n1_n2_cleaned_cases.RDS and n1_n2_cleaned_cases.csv files.
4. Pull recovery data from the bcov_recovery.csv file for NWSS reporting.
