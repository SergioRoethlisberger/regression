\name{lm_function}
\alias{lm_function}

\title{
  Fitting Linear Models
}
\description{
  This function can be used to create a simple linear regression model.
}
\usage{
  lm_function(Y, X, data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
  the dependent variable is requested as a matrix (Y = as.matrix(Y)).
  }
  \item{X}{
  the independent variable must include the intercept and the X-values (X = cbind(1, X1, X2, ..., Xk)).
  }
  \item{data}{
  an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
 }
 }
\details{
The "linear regression" attempts to model the relationship between two or more variables by fitting a linear equation to observed data.
The most common method for fitting a regression line is the method of least-squares. This method calculates the best-fitting line for the observed data by minimizing the sum of the squares of the vertical deviations from each data point to the line. 
}

\author{
  Catherine Ammann: \email{(catherine.ammann@uzh.ch)} and Sergio Roethlisberger: 
    \email{(sergio.roethlisberger@uzh.ch)}
  }

\examples{
## model = lm_function(Y = as.matrix(my.data$Y), X=cbind(1, my.data$X1, my.data$X2), data = my.data)
## model
}

\keyword{ models }