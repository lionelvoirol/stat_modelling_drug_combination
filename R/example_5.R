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

# compute all possible treatments
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
dim(all_treatments_no_duplicate)
vec_all_treatments_no_duplicate = apply(all_treatments_no_duplicate, 
                                        MARGIN = 2,
                                        FUN = function(x){paste(x, collapse="")})
length(vec_all_treatments_no_duplicate)

# create all interactions from example
lst_interactions = list()
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
  lst_interactions[[i-1]] = vec_to_add
  # vec_interactions = c(vec_interactions, vec_all_comb_i_no_duplicate)
}

# create matrix Z encoding the marginal effects
Z_all_treatments  = matrix(data = 0, nrow = dim(all_treatments_no_duplicate)[2], ncol=length(all_drugs_and_dosages))
colnames(Z_all_treatments) = all_drugs_and_dosages


# fill matrix encoding marginal treatment
for(i in seq(dim(all_treatments_no_duplicate)[2])){
  idx_column_to_fill = match(all_treatments_no_duplicate[,i], colnames(Z_all_treatments))
  idx_row_to_fill = (i):(1*i)
  Z_all_treatments[idx_row_to_fill, idx_column_to_fill] = 1
}

# create now all interactions of order 1,2 and 3 
all_interactions <- model.matrix(~ .^4, data.frame(Z_all_treatments))
colnames(all_interactions)

# remove interactons that concerns the same drugs expressed at different dosage
colnames(all_interactions) = gsub(":", "", colnames(all_interactions))

# Function to check for duplicate letters
has_duplicate_letter <- function(string) {

  letters_and_dosage <- strsplit(string, "")[[1]]
  # Identify numeric elements
  numeric_elements <- grepl("[0-9]", letters_and_dosage)
  # Keep only the non-numeric elements (letters)
  letters_only_vector <- letters_and_dosage[!numeric_elements]
  any(duplicated(letters_only_vector))
}

# find column to remove
id_col_to_remove = which(sapply(colnames(all_interactions) , has_duplicate_letter))
X_mat_all_treatments = all_interactions[,-id_col_to_remove]
dim(X_mat_all_treatments)
colnames(X_mat_all_treatments)

# define proportion of treatment screened
prop_treatment_screen = .1
id_treatment_sampled = sample(1:dim(X_mat_all_treatments)[1], size = prop_treatment_screen * dim(X_mat_all_treatments)[1])
X_mat_treatment_experimented = X_mat_all_treatments[id_treatment_sampled, ]

# define number of replicate
nbr_replicate = 20
duplicate_rows <- function(matrix, n) {
  duplicated_indices <- rep(1:nrow(matrix), each = n)
  duplicated_matrix <- matrix[duplicated_indices, , drop = FALSE]
  return(duplicated_matrix)
}

X_mat = duplicate_rows(X_mat_treatment_experimented, nbr_replicate)
dim(X_mat)

# create beta vector
beta = vector(mode = "numeric", length = dim(X_mat)[2])
names(beta) = colnames(X_mat)

# assign effect value to marginal effect and 
nbr_marginal = nbr_drugs * nbr_dosages
prop_signif_1st_order_int = .2
nbr_1st_order_int = length(lst_interactions[[1]])
prop_signif_2nd_order_int = .1
nbr_2nd_order_int = length(lst_interactions[[2]])
prop_signif_3rd_order_int = .05
nbr_3rd_order_int = length(lst_interactions[[3]])

# assign value to beta
seed_value = 123456789
set.seed(seed_value)
beta[1:nbr_marginal] = rnorm(mean = 15, sd = 5, n = nbr_marginal)
set.seed(seed_value)
# assign first order interaction effects
id_signif_1st_int = sample((nbr_marginal+1):(nbr_1st_order_int+ nbr_marginal),
                           size = round(prop_signif_1st_order_int * nbr_1st_order_int) )
beta[id_signif_1st_int] = rnorm(mean = 15, sd = 5, n = round(prop_signif_1st_order_int * nbr_1st_order_int) )
# assign second order interaction effects
id_signif_2nd_int = sample((nbr_marginal+nbr_1st_order_int+1):(nbr_marginal +nbr_1st_order_int+nbr_2nd_order_int ),
                           size = round(prop_signif_2nd_order_int * nbr_2nd_order_int) )
beta[id_signif_2nd_int] = rnorm(mean = 15, sd = 5, n = round(prop_signif_2nd_order_int * nbr_2nd_order_int)  )

# assign third order interaction effects
id_signif_3rd_int = sample((nbr_marginal+nbr_1st_order_int+nbr_2nd_order_int+1):(nbr_marginal +nbr_1st_order_int+nbr_2nd_order_int+nbr_3rd_order_int ),
                           size = round(prop_signif_3rd_order_int * nbr_3rd_order_int) )
beta[id_signif_3rd_int] = rnorm(mean = 15, sd = 5, n = round(prop_signif_3rd_order_int * nbr_3rd_order_int)  )


# identify non zero beta
which(beta!=0)

# generate data
sigma2 = 10
eps=rnorm(n = dim(X_mat)[1], mean = 0, sd = sqrt(sigma2))
y = X_mat%*% beta + eps

# compute true best drugs
true_treatment_effect_all_treatment = X_mat_all_treatments %*% beta
df_all_treatments_true_effect = data.frame("treatment" = vec_all_treatments_no_duplicate,
                                           "treatment_effect" = true_treatment_effect_all_treatment)
df_all_treatments_true_effect = df_all_treatments_true_effect %>% arrange(true_treatment_effect_all_treatment)
head(df_all_treatments_true_effect)

# perform lasso as initial estimator
fit_lasso_init = glmnet::glmnet(x = X_mat, y = y, intercept = F, family ="gaussian")

# define fct for AIC
aic_lasso <- function(lasso_fit, Xi, y, k = 2, intercept=F) {
  
  # AIC based on RSS
  n = nrow(Xi)
  if(intercept){
    pred <- Xi %*% coef(lasso_fit)
  }else{
    pred <- Xi %*% coef(lasso_fit)[-1,]
  }
  rss <- apply((y - pred)^2, 2, sum) 
  non_zero_coef = apply(coef(lasso_fit), MARGIN = 2, FUN = function(x){sum(x!=0)})
  aic = k * non_zero_coef + n*log(rss) 
  sel = which.min(aic)
  return(sel)
}


# select by AIC
select_lambda <- aic_lasso(fit_lasso_init, as.matrix(X_mat), as.vector(y), k = 2, intercept = F)

# extract vector of estimated beta from the initial lasso
beta_lasso_aic <- coef(fit_lasso_init)[, select_lambda]

# get beta
beta_marginal = beta_lasso_aic[2:(nbr_marginal+1)]
beta_marginal
beta_1st_interaction =  beta_lasso_aic[(nbr_marginal+1+1):(nbr_marginal + nbr_1st_order_int+1)]
head(beta_1st_interaction)
tail(beta_1st_interaction)
beta_2nd_interaction =  beta_lasso_aic[(nbr_marginal+nbr_1st_order_int+1+1):(nbr_marginal + nbr_1st_order_int+nbr_2nd_order_int+1)]
head(beta_2nd_interaction)
tail(beta_2nd_interaction)
beta_3rd_interaction =  beta_lasso_aic[(nbr_marginal+nbr_1st_order_int+nbr_2nd_order_int+1+1):(nbr_marginal + nbr_1st_order_int+nbr_2nd_order_int+nbr_3rd_order_int+1)]
head(beta_3rd_interaction)
tail(beta_3rd_interaction)


# create vector of weights
c1 = mean(beta_marginal != 0 ) * mean(abs(beta_marginal[which(beta_marginal != 0)]))
c1
c2 = mean(beta_1st_interaction != 0 ) * mean(abs(beta_1st_interaction[which(beta_1st_interaction != 0)]))
c2
c3 = mean(beta_2nd_interaction != 0 ) * mean(abs(beta_2nd_interaction[which(beta_2nd_interaction != 0)]))
c3
c4 = mean(beta_3rd_interaction != 0 ) * mean(abs(beta_3rd_interaction[which(beta_3rd_interaction != 0)]))
c4

# create vector of weight
pf_adaptive_lasso = c(rep(1/c1, length(beta_marginal)),
                      rep(1/c2, length(beta_1st_interaction)),
                      rep(1/c3, length(beta_2nd_interaction)),
                      rep(1/c4, length(beta_3rd_interaction)))

# perform adaptive lasso
fit_adalasso= glmnet::glmnet(x = X_mat, y = y, 
                             intercept = F, family ="gaussian",
                             penalty.factor = pf_adaptive_lasso
                             )


# select by AIC
select_lambda_adalasso <- aic_lasso(fit_adalasso, as.matrix(X_mat), as.vector(y), k = 2, intercept = F)

# extract vector of estimated beta from the initial lasso
beta_adalasso_aic <- coef(fit_adalasso)[, select_lambda_adalasso]

# predict expected treatment effect for all treatments
predicted_treatment_effect_all_treatments = X_mat_all_treatments %*% coef(fit_adalasso)[, select_lambda_adalasso][-1]

# add to treatment df
df_all_treatments_true_effect$predicted_treatment_effect = predicted_treatment_effect_all_treatments
