# create dataframe
df = data.frame("id"= seq(9),
                "measure" = c(70, 80, 85, 40, 45, 50, 25,30 , 35),
                "treatment" =c(rep("A", 3), rep("B", 3), rep("C", 3)),
                "replicate" = rep(c(1,2,3), times = 3),
                "drug_A" = c(1,1,1,0,0,0,0,0,0),
                "drug_B"=c(0,0,0,1,1,1,0,0,0),
                "drug_C" = c(0,0,0,0,0,0,1,1,1))
# print df
df

# estimate marginal effect with means
mean(df[which(df$treatment == "A"), "measure"])
mean(df[which(df$treatment == "B"), "measure"])
mean(df[which(df$treatment == "C"), "measure"])

# estimate marginal effect with linear regression
fit = lm(measure ~ treatment-1, data = df) # the -1 indicates no intercept
summary(fit)
coef(fit)
