# clean ws
rm(list = ls())

# transparent color
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  t.col
}

# generate design matrix
n = 200
p = 2
# X = dplyr::bind_cols(rep(1, n),
#                      matrix(data = rnorm(n*(p-1)),
#                             ncol=p-1, nrow=n)
#                      )
X = matrix(data = rnorm(n*(p)),
           ncol=p, 
           nrow=n)
colnames(X) = paste("X", 1:p, sep="_")
sigma2 = 5
beta = c(1.5, 10)

# construct signal
seed_value=1234567
set.seed(seed_value)
eps = rnorm(n = n , sd = sqrt(sigma2))
y = X%*% beta + eps

# construct df
df = data.frame("id" = seq(n), "y"=y, X)

# create matrix for all possible value fo parameters
delta_contour = 5
x1 = seq(beta[1]-delta_contour, beta[2]+delta_contour, by=.05)
x2 = seq(beta[1]-delta_contour, beta[2]+delta_contour, by=.05)
Z_mat = matrix(NA, ncol=length(x1), nrow=length(x2))
for(i in seq(length(x1))){
  for(j in seq(length(x2))){

    resid = y - X %*% c(x1[i], x2[j])
    Z_mat[i,j] = norm(resid, type = "2")^2
  }
}

# plot ellipse
beta_hat = solve(t(X)%*%X )%*%t(X)%*%y







#------------ matriy define plot

# tikz
library(tikzDevice)
tikz(file = "img/plot_1.tex", width = 4, height = 5, standAlone = T)

# Create a matrix to define the layout
layout_matrix <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)

# Set the layout
layout(layout_matrix, heights = c(1,.5))

par(mar=c(1,3,.5,.5))
#------------------------- Plot
# Create a blank plot
plot(NA, 
     ylim=c(-10,18),
     xlim=c(-11, 11),
     # xaxs="i",yaxs="i",
    las = 1,
     asp = 1,
     xlab = "",
     ylab = ""
    # type = "n"
    )
grid(col="grey80", lty=1)

# contour
contour(x1, x2, Z_mat,  levels = c(seq(800, 5000, by = 450)),
        drawlabels = F, las = 1, add = T)


# add lozenge
base_y_vec = c(0,-1,0,1)
base_x_vec =c(-1,0,1,0)
scaling_ols = 11.65
polygon(x = base_x_vec*scaling_ols, y=base_y_vec*scaling_ols, col =NA, lwd = 2, border = "#4f9867", lty=2, lwd = .6)


val_scale_1 = 9.7
delta_losange = 1.8
polygon(x = base_x_vec*val_scale_1, y=base_y_vec*val_scale_1, col = t_col("#a7d5ee", percent = 70), lwd = 2, border = NA)

val_scale_2 = val_scale_1-delta_losange
polygon(x = base_x_vec*val_scale_2, y=base_y_vec*val_scale_2, col = t_col("#decf83", percent = 70), lwd = 2, border = NA)

val_scale_3 = val_scale_2-delta_losange
polygon(x = base_x_vec*val_scale_3, y=base_y_vec*val_scale_3, col = t_col("#a1567c", percent = 70), lwd = 2, border = NA)

# contour(x1, x2, Z_mat,  las = 1, add = T)
abline(h = 0, lwd = 1.5)
abline(v = 0, lwd = 1.5)

points(x = beta_hat[1], y=beta_hat[2], col="#4f9867", pch=16, cex=.8)
points(x = .81, y=9, col="#a7d5ee", pch=16, cex=.8)

points(x=0, y = 7.9, col="#decf83", pch=16, cex=.8)

points(x=0, y = 6.1, col="#a1567c", pch=16, cex=.8)


# add the lambda
cex_lambda = .9
text(x = 2.5, y = 2, "$\\lambda_1$" , col= "#a1567c",cex=cex_lambda)
text(x = 5.15, y = 2, "$\\lambda_2$" , col= "#c6b975",cex=cex_lambda)
text(x = 7.1, y = 2, "$\\lambda_3$" , col= "#a7d5ee",cex=cex_lambda)
text(x = 12, y = 2, "$\\lambda=0$" , col= "#4f9867",cex=cex_lambda)



# annotations
text(x = 14.5, y = 1, "$\\beta_1$")
text(x = -1.3, y = 17, "$\\beta_2$")
# lasso path
library(glmnet)
fit = glmnet(x = X, y=y,intercept = F, data=df)


# Extract coefficients for different lambda values
beta_path <- coef(fit)

# Create a plot of the beta path
par(mar=c(3.5,3,1.5,.5))
plot(0, 0, type = "n", xlim =range(fit$lambda), 
     ylim = range(beta_path), 
     xlab = "", ylab = "", main = "" ,las=1)

# Add lines for each coefficient
col_vec = c("black", "black")
# for (i in 2:3) {
#   lines(fit$lambda, beta_path[i,], col = col_vec[i-1])
# }
col_vec = c("black", "black")


grid(col = "grey80", lty=1)


delta = .2
val_1 = -delta/2
rect(xleft = val_1,xright = val_1+delta, ybottom = -1, ytop = 50,
     col = t_col("#4f9867", percent = 80), border = NA)
abline(v=0, col="grey20", lty=2)
val_2 = .7
rect(xleft = val_2,xright = val_2+delta, ybottom = -1, ytop = 50, 
     col = t_col("#a7d5ee", percent = 80),
       # t_col("#decf83", percent = 80), 
     border = NA)

abline(v=val_2+delta/2, col="grey20", lty=2)
val_3 = 1.5
rect(xleft = val_3,xright = val_3+delta, ybottom = -1, ytop = 50, 
     col = t_col("#decf83", percent = 80), border = NA)
abline(v=val_3+delta/2, col="grey20", lty=2)
val_4 = 3.9
rect(xleft = val_4,xright = val_4+delta, ybottom = -1, ytop = 50, 
     col = t_col("#a1567c", percent = 80), border = NA)

abline(v=val_4+delta/2, col="grey20", lty=2)

for (i in 2:3) {
  
  lines(fit$lambda, beta_path[i,], col = col_vec[i-1])
}
points(x = 0, y = beta_hat[1],  col="#4f9867", pch=16, cex=.8)
points(x=0,y=beta_hat[2],  col="#4f9867", pch=16, cex=.8)

points(x = val_2+delta/2, y=9.1, col="#a7d5ee", pch=16, cex=.8)
points(x = val_2+delta/2, y=1.1, col="#a7d5ee", pch=16, cex=.8)

points(x=val_3+delta/2, y = 8.2, col="#decf83", pch=16, cex=.8)
points(x = val_3+delta/2, y=0.1, col="#decf83", pch=16, cex=.8)

points(x=val_4+delta/2, y = 0, col="#a1567c", pch=16, cex=.8)
points(x=val_4+delta/2, y = 5.99, col="#a1567c", pch=16, cex=.8)



# points(x=0, y = 7.9, col="#decf83", pch=16, cex=.8)

# points(x=0, y = 6.2, col="#a1567c", pch=16, cex=.8)
# rect(xleft = -3.67,xright = -3.55, ybottom = -1, ytop = 6, col = t_col("purple", percent = 80), border = NA)

# Add lines for each coefficient



mtext(side = 2, text = "$\\hat{\\beta}$", line = 1.8)
mtext(side = 1, text = "$\\lambda$", line = 2.2)
box()
dev.off()
system("pdflatex -output-directory=img img/plot_1.tex")
