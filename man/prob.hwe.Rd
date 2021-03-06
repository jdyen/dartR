% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prob.hwe.r
\name{prob.hwe}
\alias{prob.hwe}
\title{Exact SNP test of Hardy-Weinberg Equilibrium}
\usage{
prob.hwe(obs_hets, obs_hom1, obs_hom2)
}
\arguments{
\item{obs_hets}{-- count of heterozygotes by locus}

\item{obs_hom1}{-- count of homozygotes, reference state}

\item{obs_hom2}{-- count of homozygotes, alternate state}
}
\value{
Exact probability of agreement with HWE
}
\description{
This code calculates an exact probability of departure from Hardy-Weinberg Equilibrium 
as described in Wigginton, JE, Cutler, DJ, and Abecasis, GR (2005) A Note on Exact Tests of 
Hardy-Weinberg Equilibrium. American Journal of Human Genetics. 76:887-893.
}
\details{
Note: return code of -1.0 signals an error condition; return code of NA signals that 
all alleles are NA for a locus
}
\examples{
hets <- 20
hom_1 <- 5
hom_2 <- 30
p_value <- prob.hwe(hets, hom_1, hom_2)
}
\author{
Arthur Georges (glbugs@aerg.canberra.edu.au)
}
