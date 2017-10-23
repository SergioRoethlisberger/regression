glm_function <-
function(Y,X,stval, data){
  logit.ll <- function(X,Y,theta){
    mu <- theta[1] + X %*% theta[2:length(theta)]
    
    l.link <- function(x){
      exp(x)/(1+exp(x))
    }
    logl <- sum(Y*log(l.link(mu)) + (1 - Y)*log(1-l.link(mu)))
    return(-logl)
  }
  n = nrow(data) # sample size
  k = ncol(X) # number of betas
  
  #optimize
  res <- optim(par=stval, fn=logit.ll, Y=Y, X=X, hessian=TRUE)
  
  # (6) calculate the standard errors
  sd <- sqrt(diag(solve(res$hessian))) #gives us standarddeviation
  
  # (7) calculate z value
  z = res$par/sd
  
  # (8) calculate p value
  p = 2*pt(abs(z), df=n-(k+1), lower.tail= FALSE)
  
  # (8a) add stars
  stars = ifelse(p<0.001,"***",
                 ifelse(p<0.01,"**", 
                        ifelse(p<0.05,"*",
                               ifelse(p<0.1,"."," "))))
  p = ifelse(p < 2*10^(-16), 
             c("< 2e-16"), format(signif(p, digits = 3)))
  
  result = as.data.frame(cbind(c("(Intercept)", colnames(data[2:(k+1)])),
                               sprintf("%.6f", round(res$par, digits=6)),
                               sprintf("%.6f", round(sd, digits=6)),
                               sprintf("%.3f", round(z, digits=3)),
                               p))
  
  colnames(result) = c("Coefficients","Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  result$' ' = stars
  return(result)
}
