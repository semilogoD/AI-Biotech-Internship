# Practice Exercises 

# ----------------------------------------------------------------------------------------------------------------

# 1. Check Cholesterol level (using if) 
# Write an If statement to check cholesterol level is greater than 240, 
# if true, it will prints “High Cholesterol”

cholesterol <- 230

if(cholesterol > 240){
  print("High Cholesterol")
}

# 2. Blood Pressure Status (using if...else)
# Write an if…else statement to check if blood pressure is normal.
# If it’s less than 120, print: “Blood Pressure is normal”
# If false then print: “Blood Pressure is high”

Systolic_bp <- 130

if(Systolic_bp < 120){
  print("Blood Pressure is normal")
}else{
  print("Blood Pressure is high")
}

# 3. Automating Data Type Conversion with for loop

# Use patient_info.csv data and metadata.csv

# Using  patient_info.csv data and metadata.csv
patient <- read.csv(file.choose())
metadata <-  read.csv(file.choose())

# Creating  a copy of each dataset to work on.
mdata <- metadata
pdata <- patient

#checking the copied dataset 
str(mdata)
str(pdata)
# Identify all columns that should be converted to factor type.

# In mdata,columns height and gender should be converted to factor type 
factor_cols <- c("height","gender")

# In pdata,columns gender, diagnosis and smoker  should be converted to factor type 
factor_cols2 <- c("gender","diagnosis","smoker")


# Using a for loop to convert all the columns in factor_cols and factot_cols2 to factor type.

for (col in factor_cols) {
    mdata[[col]] <- as.factor(mdata[[col]])
}

for (col in factor_cols2) {
  pdata[[col]] <- as.factor(pdata[[col]])
}

# ----------------------------------------------------------------------------------------------------------------

# 4. Converting Factors to Numeric Codes

# Choose one or more factor columns (e.g., smoking_status).
# Convert "Yes" to 1 and "No" to 0 using a for loop.

#converting factor column to numeric 
binary_cols <- c("smoker")
for (col in binary_cols){
  pdata[[col]] <- ifelse(pdata[[col]] == "No", 0, 1)
}
   

# ----------------------------------------------------------------------------------------------------------------

#  Verification:
# Comparing the original and modified datasets to confirm changes.
str(patient)
str(pdata)

str(metadata)
str(mdata)
# ----------------------------------------------------------------------------------------------------------------


  