\name{plot.predMMSE}
\alias{plot.predMMSE}
\title{
Plot of predicted scores in the natural scale of the psychometric test
}
\description{
This function plots the predicted trajectories
obtained with \code{predMMSE} function.
}
\usage{
\method{plot}{predMMSE}(x, legend.loc = "topright", legend, add = FALSE, ...)
}

\arguments{
  \item{x}{
a predMMSE object
}
  \item{legend.loc}{
keyword for the position of the legend from the list "bottomright", "bottom",
"bottomleft", "left", "topleft","top", "topright", "right" and "center".
By default, the legend is located in the top right of the plot.
}
  \item{legend}{
character or expression to appear in the legend. If no legend should be added,
legend should be NULL.
}
  \item{add}{
optional logical indicating if the curves should be added to an existing
plot (add=TRUE) or if a new plot should be created (add=FALSE). 
By default, add is FALSE.
}
  \item{\dots}{
other parameters to be passed through to plotting functions or to legend
}
}
\author{
Cecile Proust-Lima, Viviane Philipps
}
\seealso{
\code{\link{predictMMSE}},\code{\link{normMMSE}}
}
\examples{
#import data paquid from lcmm package
library(lcmm)
data(paquid)

#normalization of MMSE scores
paquid$MMSEnorm <- normMMSE(paquid$MMSE)

#estimation of a linear mixed model on the normalized data
m <- hlme(MMSEnorm~I(age-65)*CEP, random=~I(age-65), subject="ID", data=paquid)

#prediction of MMSE scores in the 0-30 scale
pred <- predictMMSE(m,VarTime="age",Timelim=c(65,85),nTime=30,
Xprofile=c(CEP=1),methInteg='MC',nsim=200,draws=FALSE)     

#plot of the predictions
plot(pred)
}

