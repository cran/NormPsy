\name{normMMSE}
\alias{normMMSE}
\title{
Normalized scores for the Mini Mental State Examination (MMSE)
}
\description{
This function transforms crude MMSE scores (range 0 - 30) into normalized scores. The normalized scale ranges from 0 to 100, with MMSE minimum score 0 being transformed into 0 and MMSE maximum score 30 into 100. The normalization preserves the ranking of the test but transforms the gaps between two consecutive values in order to correct the curvilinearity of MMSE. The normalization is specifically adapted for heterogeneous elderly propulations mixing normal and pathological aging. 
}
\usage{
normMMSE(x)
}
\arguments{
  \item{x}{
a numeric vector containing MMSE scores (integer between 0 and 30)
}
}

\value{
a vector containing the transformed MMSE scores
}

\author{
Cecile Proust-Lima, Viviane Philipps
}


\examples{

#import data paquid from lcmm package
library(lcmm)
data(paquid)

# computation of the normalized MMSE
paquid$MMSEnorm <- normMMSE(paquid$MMSE)

# histogram of these data
par(mfrow=c(1,2))
hist(paquid$MMSE,breaks=seq(-0.5,30.5,1),col=2,main="crude MMSE")
hist(paquid$MMSEnorm,breaks=seq(0,100,10),col=3,main="normalized MMSE")
}




%\seealso{
%
%}
