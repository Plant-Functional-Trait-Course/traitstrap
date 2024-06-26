% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bootstrap_traits_multivariate.R
\name{trait_multivariate_bootstrap}
\alias{trait_multivariate_bootstrap}
\title{Bootstrap traits}
\usage{
trait_multivariate_bootstrap(
  filled_traits,
  nrep = 100,
  sample_size = 200,
  raw = FALSE,
  id = "ID",
  fun
)
}
\arguments{
\item{filled_traits}{output from the trait_fill function.}

\item{nrep}{number of bootstrap replicates}

\item{sample_size}{bootstrap size}

\item{raw}{logical; argument to extract the raw data of the trait
distributions.
The default is \code{raw = FALSE}.
If \code{raw = TRUE}, \code{nrep} is restricted to 1 to avoid memory issues.}

\item{id}{column name of unique identifiers of each leaf}

\item{fun}{bivariate or multivariate function to apply}
}
\value{
a tibble with columns for the elements of the scale_hierarchy,
and a list column result which includes the output of \code{fun}.
}
\description{
Function for nonparametric bootstrap resampling to
calculate community
weighted trait correlations, other bivariate or multivariate statistics
}
\details{
The observed and filled leaves are re-sampled in proportion to
their weights, e.g. the abundance of a species or the biomass.
Values across all individuals in a community are
resampled \code{sample_size} times to incorporate the full
spectrum of trait variation, generating \code{nrep} trait distributions.
The function \code{fun} is applied to the trait distribution at the finest level
of the filled trait hierarchy.

Unexpected columns in \code{filled_traits} are deleted with a warning.

Note that due to the flexibility of this function,
the output CAN NOT be summarized using
\code{trait_summarise_boot_moments}.
}
\examples{
require(dplyr)
require(tidyr)
require(ggplot2)
require(purrr)

data(community)
data(trait)

filled_traits <- trait_fill(
  comm = community |>
    filter(
      PlotID \%in\% c("A", "B"),
      Site == 1
    ),
  traits = trait,
  scale_hierarchy = c("Site", "PlotID"),
  taxon_col = "Taxon", value_col = "Value",
  trait_col = "Trait", abundance_col = "Cover",
  complete_only = TRUE, leaf_id = "ID"
)

# Note that more replicates and a greater sample size are advisable
# Here we set them low to make the example run quickly
boot_traits <- trait_multivariate_bootstrap(filled_traits,
  fun = cor,
  nrep = 10,
  sample_size = 100
)

boot_traits_long <- boot_traits |>
  mutate(correlations = map(result, ~ cor_to_df(.x))) |>
  select(-result) |>
  unnest(correlations)

boot_traits_long |>
  ggplot(aes(x = paste(row, "v", col), y = value)) +
  geom_violin() +
  facet_grid(Site ~ PlotID) +
  coord_flip() +
  labs(y = "Correlation", x = "")
}
