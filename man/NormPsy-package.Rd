\name{NormPsy-package}
\alias{NormPsy-package}
\alias{NormPsy}
\docType{package}
\title{
  \Sexpr[stage=build,results=hide]{lb <- library(help="NormPsy", character.only=TRUE)}
  \Sexpr[stage=build,results=hide]{descr <- packageDescription("NormPsy")}
  \Sexpr[stage=build]{descr$Title}
}
\description{
  \Sexpr[stage=build]{descr$Description}
}
\details{
  \tabular{ll}{
    Package:  \tab \Sexpr[stage=build]{descr$Package}\cr
    Type:     \tab \Sexpr[stage=build]{descr$Type}\cr
    Version:  \tab \Sexpr[stage=build]{descr$Version} \cr
    Date:     \tab \Sexpr[stage=build]{descr$Date}\cr
    License:  \tab \Sexpr[stage=build]{descr$License}\cr
    Built:    \tab \Sexpr[stage=build]{descr$Built}\cr
  }
  
  Index:
  
  \Sexpr[stage=build,results=rd]{paste("\\\\preformatted{", paste(if(!is.null(lb$info[[2]])) lb$info[[2]] else "No entries",collapse="\\n"), "}", sep="")}

}
\author{
  \Sexpr[stage=build]{descr$Author}
  
  Maintainer: \Sexpr[stage=build]{descr$Maintainer}
}
