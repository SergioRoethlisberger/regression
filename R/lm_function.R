lm_function <-
function(Y, X, data){
  # (1) estimate the betas
  beta.hat = solve(t(X)%*%X)%*%t(X)%*%Y
  
  # (2) calculate the standard errors
  y.hat = X%*%beta.hat
  
  # (3) calculate residuals
  epsilon = Y-y.hat
  
  # (4) define n and k
  n = nrow(data) # sample size
  k = ncol(X) # number of betas
  
  # (5) calculate variance-covariance-matrix
  vcov = 1/(n-k) * as.numeric(t(epsilon)%*%epsilon) * solve(t(X)%*%X)
  
  # (6) calculate the standard errors
  s.e. = sqrt(diag(vcov))
  
  # (7) calculate t value
  t = beta.hat/s.e.
  
  # (8) calculate p value
  p = 2*pt(abs(t), df=n-k, lower.tail= FALSE)
  p = ifelse(p < 2*10^(-16), 
             c("< 2e-16"), format(signif(p, digits = 3))) 
  
  # (8a) add stars
  stars = ifelse(p<0.001,"***",
                 ifelse(p<0.01,"**", 
                        ifelse(p<0.05,"*",
                               ifelse(p<0.1,"."," "))))
  
  # (9) calculate residual quartils
  res.distr = quantile(epsilon, p=seq(0,1,length.out=5))
  
  residual.error = sqrt(sum((y.hat-Y)^2)/(n-k))
  
  # (10) calculate multiple rsquared
  multiple.r.squared = sum((y.hat-mean(Y))^2)/sum((Y-mean(Y))^2)
  
  # (11) calculate adjusted rsquared
  adjusted.r.squared = 1-(1-multiple.r.squared)*(n-1)/(n-k)
  
  # (12) calculate F-statistics and its p-value
  f =  multiple.r.squared/(k-1)/((1-multiple.r.squared)/(n-k))
  f.p = pf(f, k-1, n-k, lower.tail=F)
  f.p = replace(f.p, which(f.p<(2*10^(-16))), "< 2.2e-16")
  
  # (13) generate a data frame
  res = as.data.frame(rbind(res.distr))
  colnames(res) = c("Min", "1Q", "Median", "3Q", "Max")
  rownames(res) = c(" ")
  
  result = as.data.frame(cbind(c("(Intercept)", colnames(data[2:k])),
                               sprintf("%.6f", round(beta.hat, digits=6)),
                               sprintf("%.6f", round(s.e., digits=6)),
                               sprintf("%.3f", round(t, digits=3)),
                               p))
  colnames(result) = c("Coefficients","Estimate", "Std. Error", "t value", "Pr(>|t|)")
  result$' ' = stars
  
  # (14) print the output
  cat(paste0("Residuals:", "\n"))
  print(round(res, digits = 4))
  cat(paste0("", "\n", ""))
  print(result)
  cat(paste0("---",
             "\n",
             "Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1",
             "\n",
             "---",
             "\n",
             "Residual standard error: ", round(residual.error, digits = 3), " on ", n-k, " degrees of freedom",
             "\n", 
             "Multiple R-squared: ", round(multiple.r.squared, digits = 5),
             ", Adjusted R-squared: ", round(adjusted.r.squared, digits = 5),
             "\n",
             "F-Statistics: ", round(f, digits = 2), " on ", k-1, " and ", n-k, " DF",
             ", p-value: ", f.p))
}
