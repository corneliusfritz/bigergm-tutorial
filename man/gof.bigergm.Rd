% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gof.bigergm.R
\name{gof.bigergm}
\alias{gof.bigergm}
\title{Conduct Goodness-of-Fit Diagnostics on a Exponential Family Random Graph
Model for big networks}
\usage{
\method{gof}{bigergm}(
  object,
  ...,
  type = "full",
  control_within = ergm::control.simulate.formula(),
  seed = NULL,
  nsim = 100,
  compute_geodesic_distance = TRUE,
  start_from_observed = TRUE,
  simulate_sbm = FALSE
)
}
\arguments{
\item{object}{An \code{\link{bigergm}} object.}

\item{...}{Additional arguments, to be passed to \code{\link{simulate_bigergm}},
which, in turn, passes the information to \code{\link[ergm]{simulate_formula}}.
See documentation for \code{\link{bigergm}}.}

\item{type}{the type of evaluation to perform. Can take the values \code{full} or \code{within}. \code{full} performs the evaluation on all edges, and \code{within} only considers within-block edges.}

\item{control_within}{MCMC parameters as an instance of \code{control.simulate.formula} to be used for the within-block simulations.}

\item{seed}{the seed to be passed to simulate_bigergm. If \code{NULL}, a random seed is used.}

\item{nsim}{the number of simulations to employ for calculating goodness of fit, default is 100.}

\item{compute_geodesic_distance}{if \code{TRUE}, the distribution of geodesic distances is also computed (considerably increases computation time on large networks. \code{FALSE} by default.)}

\item{start_from_observed}{if \code{TRUE}, MCMC uses the observed network as a starting point. If \code{FALSE}, MCMC starts from a random network.}

\item{simulate_sbm}{if \code{TRUE}, the between-block connections are simulated from the estimated stochastic block model from the first stage not the estimated ERGM.}
}
\value{
\code{\link{gof.bigergm}} returns a list with two entries.
The first entry 'original' is another list of the network stats, degree distribution, edgewise-shared partner distribution, and geodesic distance distribution (if \code{compute_geodesic_distance = TRUE}) of the observed network.
The second entry is called 'simulated' is also list compiling the network stats, degree distribution, edgewise-shared partner distribution, and geodesic distance distribution (if \code{compute_geodesic_distance = TRUE}) of all simulated networks.
}
\description{
A sample of graphs is randomly drawn from the specified model.  The first
argument is typically the output of a call to \code{\link{bigergm}} and the
model used for that call is the one fit.

By default, the sample consists of 100 simulated networks, but this sample
size (and many other settings) can be changed using the \code{ergm_control}
argument described above.
}
\examples{
data(toyNet)
\donttest{
# Specify the model that you would like to estimate.
data(toyNet)
# Specify the model that you would like to estimate.
model_formula <- toyNet ~ edges + nodematch("x") + nodematch("y") + triangle
estimate <- bigergm(model_formula,n_blocks = 4)
gof_res <- gof(estimate,
nsim = 100
)
plot(gof_res)
}
}
