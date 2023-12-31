# clean ws
rm(list=ls())

# pkg
library(magrittr)
library(dplyr)
library(progress)

# define drugs
drugs = c("A", "B", "C")
dosages = c(1,2)
df_drugs_dosages=expand.grid(drugs, dosages)
all_drugs_and_dosages = sort(paste(df_drugs_dosages[,1], df_drugs_dosages[,2], sep=""))

# compute all treatments
nbr_of_drug_per_treatment = 2
all_treatments = combn(all_drugs_and_dosages, m = nbr_of_drug_per_treatment )

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

# create vector for variable name of interactions of 1st order
all_int_1 = paste(all_treatments_no_duplicate[1,], all_treatments_no_duplicate[2,], sep="")

# create design matrix 
# matrix X1 design the marginal effect 
nbr_replicates = 2
Z  = matrix(data = 0, nrow = nbr_replicates * dim(all_treatments_no_duplicate)[2], ncol=length(all_drugs_and_dosages))
colnames(Z) = all_drugs_and_dosages

# fill matrix encoding marginal treatment
for(i in seq(dim(all_treatments_no_duplicate)[2])){
  idx_column_to_fill = match(all_treatments_no_duplicate[,i], colnames(Z))
  idx_row_to_fill = (nbr_replicates*i-(nbr_replicates-1)):(nbr_replicates*i)
  Z[idx_row_to_fill, idx_column_to_fill] = 1
}

# create matrix encoding first order interactions
n = nbr_replicates * dim(all_treatments_no_duplicate)[2]
X1 = matrix(0,
            nrow = n,
            ncol=length(all_int_1))
colnames(X1) = all_int_1

# fill matrix X1
for(i in seq(all_int_1)){
  idx = match(all_int_1[i], colnames(X1))
  X1[(nbr_replicates*i-(nbr_replicates-1)):(nbr_replicates*i), idx] = 1
}

# create dataset and merge
df = data.frame("id" = seq(n),
                "combination" = rep(all_int_1, each=nbr_replicates))
df = dplyr::bind_cols(df, Z, X1)

# create vector of beta
beta = vector(mode = "numeric",length = length(all_drugs_and_dosages)+ length(all_int_1))
names(beta) = c(all_drugs_and_dosages, all_int_1)

# assign marginal effect value
beta[1:length(all_drugs_and_dosages)] = sample(c(-10:-1, 1:10), length(all_drugs_and_dosages))

# assign value to some of the interactions, say 40% of them
percentage_significant_int_1 = 40
nbr_significant_int_1 = round(percentage_significant_int_1/100 * length(all_int_1))
seed_value = 123456
set.seed(seed_value)
id_significant_int_1 = sample((length(all_drugs_and_dosages)+1):length(beta), size = nbr_significant_int_1)
set.seed(seed_value)
beta[id_significant_int_1] = sample(c(-10:-1, 1:10), nbr_significant_int_1)
beta

# generate data from the model
sigma_2 = 2
X_mat = as.matrix(df[,3:dim(df)[2]])
X_mat = matrix(as.numeric(X_mat), nrow = nrow(X_mat))

# create signal
dim(X_mat)
length(beta)
set.seed(seed_value)
eps = rnorm(n = n, sd = sqrt(sigma_2))
y = X_mat %*% beta + eps

# add y to df
df$y = y

# check colnames y
colnames(df)

# we now want to list all possible model that consider all marginal effect and all possible combinations of interactions

# Create an empty list to store model formulas
model_formulas <- list()

# Loop over all possible dimensions (1 to length(variables))
for (i in 1:length(all_int_1)) {

  # Generate all combinations of variables of length i
  variable_combinations <- combn(all_int_1, i, simplify = FALSE)
  
  # Create model formulas for each combination
  formulas <- lapply(variable_combinations, function(combination) {
    # Combine the variables into a formula string
    formula_string <- paste("y ~ A1 + A2 + B1 + B2 + C1 + C2 +", 
                            paste(combination, collapse = " + "))
    formula_string
  })
  
  # Append the formulas for this dimension to the list
  model_formulas[[i]] <- formulas
}

# Flatten the list of formulas
all_model_formulas <- unlist(model_formulas, recursive = FALSE)

# count number of model and verify if equal to 2^p
length(all_model_formulas) == 2^length(all_int_1)

# we only miss the model with only the marginal effect
all_model_formulas = append(all_model_formulas, "y ~ A1 + A2 + B1 + B2 + C1 + C2")

# construct dataframe with all model to record AIC and BIC
df_all_model = data.frame("model" = unlist(all_model_formulas),
                          "AIC"=NA,
                          "BIC"=NA)

# estimate all possible model and compute a measure of fit based on the likelihood
pb <- progress_bar$new(total = length(all_model_formulas))
for(i in 1:length(all_model_formulas)){

  fit = lm(formula = as.formula(all_model_formulas[[i]]), data = df)
  df_all_model[i, "AIC"]= AIC(fit)
  df_all_model[i, "BIC"]= BIC(fit)
  # print(i)
  pb$tick()
}

# sort models per AIC value
df_all_model = df_all_model %>% arrange(AIC)
min_aic = min(df_all_model$AIC)
upper_bound = min_aic +2
equivalent_model = which(df_all_model$AIC < upper_bound)

# construct df with only selected model based on min AIC + 2 as an upper bound
df_equivalent_model = df_all_model[equivalent_model, ]

# define true model based on beta
true_model = paste("y ~ ", paste(names(beta[which(beta != 0)]), collapse = " + "), sep="")

# identify rank of true model
which(df_all_model$model == true_model)

# check if true model is included in selected model
true_model %in% df_equivalent_model$model

# check number of selected model
dim(df_equivalent_model)[1]
