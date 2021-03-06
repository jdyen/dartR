% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gl2faststructure.r
\name{gl2faststructure}
\alias{gl2faststructure}
\title{Export DArT genlight object \{adegenet\} to faststructure format (to run faststructure elsewhere)}
\usage{
gl2faststructure(gl, outfile = "gl.str", probar = TRUE)
}
\arguments{
\item{gl}{-- genlight object}

\item{outfile}{-- filename of the output fastA file [default genlight.fasta]}

\item{probar}{switch to show/hide progress bar}
}
\description{
Recodes in the quite specific faststructure format (e.g first six columns need to be there, but are ignored...check faststructure documentation (if you find any :-( )))
}
\details{
The script writes out the a file in faststructure format.
}
\author{
Bernd Gruber (glbugs@aerg.canberra.edu.au)
}
