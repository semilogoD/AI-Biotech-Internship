# creating sub folders
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

# loading data set 
Patient_data <- read.csv(file.choose())

# inspecting the patient_data file
str(Patient_data)
summary(Patient_data)

# identify variables with incorrect or inconsistent data types
# gender
# diagnosis
# smoker

# Convert variables to appropriate data types where needed
Patient_data$gender_fac <- as.factor(Patient_data$gender)
Patient_data$diagnosis_fac <- as.factor(Patient_data$diagnosis)
Patient_data$smoker_fac <-as.factor(Patient_data$smoker)

# Create a new variable for smoking status as a binary factor
Patient_data$smoker_num <- ifelse(Patient_data$smoker_fac == "Yes",1,0)

# Save the cleaned dataset in your clean_data folder with the name patient_info_clean.csv
write.csv(Patient_data, file = "clean_data/paient_info_clean.csv")
