% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/script.R
\name{best_information_gain}
\alias{best_information_gain}
\title{Find the attribute that provides the best information Information Gain as a splitting attribute}
\usage{
best_information_gain(x, f)
}
\arguments{
\item{x}{the target vector i.e. the entire set of labels (target variable values)}

\item{f}{a data frame with the candidate attributes as columns}
}
\value{
A list with two members: the value of the best information gain and the column index in \code{f} of the 'winning' attribute
}
\description{
\code{best_information_gain} investigates a set of attributes of a dataset to find the one that provides the best Information Gain as a splitter for a given target variable.
}
\seealso{
\code{\link{information_gain}}
}
