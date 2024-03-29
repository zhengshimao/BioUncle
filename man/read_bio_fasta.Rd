% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_bio_fasta.R
\name{read_bio_fasta}
\alias{read_bio_fasta}
\title{Create a data frame from a fasta file.}
\usage{
read_bio_fasta(fasta, simplified_ids = FALSE)
}
\arguments{
\item{fasta}{A fasta file.}

\item{simplified_ids}{A logical value(default:\code{FALSE}). \code{TRUE} for simplified ids of fasta sequence with characters between ‘>’ and first space, and \code{FALSE} for not simplified its ids.}
}
\value{
A data frame.
}
\description{
Reads a file in fasta format and creates a data frame from it.
}
