% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scCNMF.R
\name{scCNMFmodel}
\alias{scCNMFmodel}
\title{Integration operations are performed on the two matrices X1,X2 obtained after initialization, and the rcpp program is called in the middle.
Function 'scCNMF_cpp_function' calls c++ program 'scCNMF.cpp'}
\usage{
scCNMFmodel(
  X1,
  X2,
  s = 0.25,
  K = 20,
  alpha = 1,
  beta = 1,
  gamma = 10000,
  maxIter = 500,
  stop_rule = 1,
  GeneName,
  PeakName,
  GeneLoc,
  PeakLoc,
  feature_cut_perc = 0.01,
  corr_cut_k = 1e+05,
  core = 8
)
}
\arguments{
\item{X1}{scRNA-seq data matrix}

\item{X2}{scATAC-seq data matrix}

\item{s}{Probability of Bernoulli distribution}

\item{K}{The rank of the latent factor}

\item{alpha}{Model parameter}

\item{beta}{Model parameter}

\item{gamma}{Model parameter}

\item{maxIter}{Maximum number of iterations}

\item{stop_rule}{Iteration stopping conditions, taking a fixed number of iterations when equal to 1 and adaptive iterations when equal to 2}

\item{GeneName}{Name of genes}

\item{PeakName}{Name of the ATAC peaks}

\item{GeneLoc}{Location of genes on chromosomes}

\item{PeakLoc}{Location of peaks on chromosomes}

\item{feature_cut_perc}{Feature cut percent}

\item{corr_cut_k}{Top peak-gene pairs}

\item{core}{Core}
}
\description{
Integration operations are performed on the two matrices X1,X2 obtained after initialization, and the rcpp program is called in the middle.
Function 'scCNMF_cpp_function' calls c++ program 'scCNMF.cpp'
}
