\name{predictMMSE}
\alias{predictMMSE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction of MMSE scores in their natural scale.
}
\description{
From the estimates of a linear mixed model applied on normalized MMSE scores, the function computes the predicted values of MMSE in its natural 0-30 scale. The function also provides 95\% confidence intervals computed by a Monte Carlo method. 
}
\usage{
predictMMSE(model, VarTime, Timelim, nTime, Xprofile, 
methInteg = "GH", nsim = 20, draws = FALSE, ndraws = 2000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
a hlme object representing a latent class linear mixed model 
}
%  \item{newdata}{
%si pas VarTime, Timelim,Xprofile?
%}
  \item{VarTime}{
  a character string containing the name of the time variable in the model
}
  \item{Timelim}{
  a numeric vector indicating the time limits of the prediction
}
  \item{nTime}{
  an integer indicating the number of prediction times
}
  \item{Xprofile}{
  a named vector containing the values of the model's covariate
}
  \item{methInteg}{
integration method used to compute the numerical integration. If 0 or 'GH' (the default) a Gauss-Hermite numerical approximation is employed. Value 1 or 'MC' specifies a Monte Carlo method.
}
  \item{nsim}{
number of integration points used for the integration method specified in argument methInteg
}
  \item{draws}{
optional logical indicating if 95\% confidence intervals of the predictions should be calculated. If TRUE, the posterior distribution of the predicted values is approximated by a Monte Carlo method. The predicted value is given by the median value, and the confidence interval is given by the 2.5\% and 97.5\% percentiles of this distribution. If FALSE (the default), a point prediction is provided.
}
  \item{ndraws}{
optional integer indicating the number of draws in the Monte Carlo approximation if draws=TRUE.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
a predMMSE object consisting in a matrix containing the prediction time, the predicted values and, if draws=TRUE, the lower and upper limits of the confidence interval.
}

\author{
Cecile Proust-Lima, Viviane Philipps
}

\seealso{
\code{\link{plot.predMMSE}},\code{\link{normMMSE}}
}
\examples{

#import data paquid from lcmm package
library(lcmm)
data(paquid)

# computation of the normalized MMSE
paquid$MMSEnorm <- normMMSE(paquid$MMSE)

# estimation of a model on the normalized MMSE scores
m <- hlme(MMSEnorm~I(age-65)*CEP, random=~I(age-65), subject="ID", data=paquid)

# prediction from the model in the natural scale of MMSE
predictMMSE(m,VarTime="age",Timelim=c(65,85),nTime=30,
Xprofile=c(CEP=1),methInteg='MC',nsim=200,draws=FALSE)

}

