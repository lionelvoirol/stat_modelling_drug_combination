# clean ws
rm(list=ls())

# pkg
library(magrittr)
library(dplyr)
library(progress)

# define drugs
nbr_drugs = 11
drugs = LETTERS[1:nbr_drugs]
dosages = c(1,2)
nbr_dosages = length(dosages)
df_drugs_dosages=expand.grid(drugs, dosages)
all_drugs_and_dosages = sort(paste(df_drugs_dosages[,1], df_drugs_dosages[,2], sep=""))

# compute all treatments
nbr_of_drug_per_treatment = 4
all_treatments = combn(all_drugs_and_dosages, m = nbr_of_drug_per_treatment)

# we need to remove treatments which contain the same drugs at different dosage as they are impossible in practice.
# define fct to identify same drugs in matrix of combination
identify_same_drug <- function(column) {
  first_letters <- substr(column, 1, 1)  # Extract the first letter of each element in the column
  length(unique(first_letters)) < length(column)  # Check if all first letters are the same
}

# remove treatment where the same drug is at the same dosage
id_treatment_same_drug_diff_dosage = which(apply(all_treatments, MARGIN = 2, FUN = identify_same_drug))
id_treatment_same_drug_diff_dosage
all_treatments_no_duplicate = all_treatments[, -id_treatment_same_drug_diff_dosage]


dim(all_treatments_no_duplicate)[2]
dim(combn(length(drugs), nbr_of_drug_per_treatment))[2] * length(dosages)^nbr_of_drug_per_treatment

# create all interactions from example
vec_interactions = c()
for(i in 2:nbr_of_drug_per_treatment){
  all_comb_i = combn(all_drugs_and_dosages, m = i)
  # remove treatment where the same drug is at the same dosage
  id_treatment_same_drug_diff_dosage = which(apply(all_comb_i, MARGIN = 2, FUN = identify_same_drug))
  all_comb_i_no_duplicate = all_comb_i[, -id_treatment_same_drug_diff_dosage]
  vec_to_add= apply(all_comb_i_no_duplicate, MARGIN = 2, FUN = function(x){paste(x, collapse ="-")})
  # print(length(vec_to_add))
  # print(tail(vec_to_add))
  vec_all_comb_i_no_duplicate = vec_to_add
  # append
  vec_interactions = c(vec_interactions, vec_all_comb_i_no_duplicate)
}


length(vec_interactions)
head(vec_interactions)
tail(vec_interactions)

# Compute nbr of treatment from formula
total_int = 0
for(i in 2:nbr_of_drug_per_treatment){
  nbr_int_i = dim(combn(1:(nbr_drugs), i))[2] * nbr_dosages^i
  total_int = total_int + nbr_int_i
}
total_int
# check if equal
total_int == length(vec_interactions)

