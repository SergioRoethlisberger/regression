\name{glm_function}
\alias{glm_function}

\title{
  Fitting Logit Models
}
\description{
  The function "glm_function" is used to fit logit models, specified by giving a symbolic description of the linear predictor and a description of the error distribution.}
  
\usage{
  glm_function(Y, X, stval, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  the dependent variable is a dichotomous dummy, taking the values of 0 and 1.
  }
  \item{X}{
  the independent variable must include the X-values (X = cbind(X1, X2, ..., Xk)).
  }
  \item{stval}{
  a named list of starting values for the parameters in the model (stval = c()).
  }
  \item{data}{
  an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
  }
}
\details{
The "logit model" is a regression model where the dependent variable is categorical. This function covers the case of a binary dependent variable - that is, where the output can take only two values, 0 and 1, which represent outcomes such as pass/fail or win/lose.
How can we interpret the output? The coefficients are on the linked scale - so a direct interpretation is not possible. The only thing we can interpret is the direction. A positive coefficient means an increase of the probability, a negative means a decrease.
}

\author{
  Catherine Ammann: \email{(catherine.ammann@uzh.ch)} and Sergio Roethlisberger: 
    \email{(sergio.roethlisberger@uzh.ch)}
  }
      
\examples{
## model = glm_function(my.data$Y, cbind(my.data$X1, my.data$X2), my.data, c(z0, z1, ..., zk))
## model
}

\keyword{ models }