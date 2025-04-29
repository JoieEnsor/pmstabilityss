---
output:
  pdf_document: default
  html_document: default
---
\name{pmstabilityss}
\alias{pmstabilityss}
\title{Sample Size for Prediction Model Development Based on Prediction Stability}
\description{
  Calculates the minimum sample size required for developing a prediction model targeting precise individual risk estimates.
}
\usage{
pmstabilityss(data, varlist, prevalence, lp = NULL, subgroup = NULL,
              nodraw = FALSE, nocompare = FALSE, pmss = FALSE, n = NULL,
              cstatistic = 0, simobs = 500000, seed = 1234,
              pciwidth = NULL, pcutpoints = NULL, threshold = 0,
              tolerance = 0.005, logitpvarincrement = 0.001,
              color = "darkgrey")
}
\arguments{
  \item{data}{A data frame containing the variables for analysis.}
  \item{varlist}{Character vector of variable names in the linear predictor.}
  \item{prevalence}{Numeric value of prevalence for the outcome.}
  \item{lp}{Character name of linear predictor variable (optional, used instead of prevalence and cstatistic).}
  \item{subgroup}{Character name of subgroup variable (optional).}
  \item{nodraw}{Logical, if TRUE, suppresses drawing plots.}
  \item{nocompare}{Logical, if TRUE, suppresses comparison plots.}
  \item{pmss}{Logical, if TRUE, uses pmsampsize package to calculate sample size.}
  \item{n}{Numeric vector of sample sizes to evaluate.}
  \item{cstatistic}{Numeric C-statistic value (used for pmsampsize).}
  \item{simobs}{Number of observations for simulation (default 500000).}
  \item{seed}{Random seed value (default 1234).}
  \item{pciwidth}{Numeric vector of target prediction interval widths.}
  \item{pcutpoints}{Numeric vector of probability cutpoints (must match pciwidth).}
  \item{threshold}{Numeric classification threshold (default 0).}
  \item{tolerance}{Numeric convergence tolerance (default 0.005).}
  \item{logitpvarincrement}{Numeric increment for logit p variance (default 0.001).}
  \item{color}{Plot color (default "darkgrey").}
}
\details{
  \code{pmstabilityss} computes the minimum sample size required for the development of a new clinical prediction model targeting precise individual risk estimates as proposed by Riley et al. 2024. 
  
  The function calculates prediction intervals for individual risk estimates at different sample sizes, allowing researchers to identify the sample size needed to achieve desired precision in risk predictions. This approach aims to minimize overfitting and ensure precise estimation of key parameters in the prediction model.
  
  When \code{pciwidth} and \code{pcutpoints} are specified, the function calculates the sample size needed to achieve the target prediction interval widths at different probability levels. For example, you might want narrower intervals for low-risk patients and can accept wider intervals for high-risk patients.
  
  If \code{threshold} is specified, the function also calculates the probability of misclassification at the given threshold.
  
  If \code{pmss = TRUE}, the function uses the pmsampsize package to calculate a sample size based on established criteria and compares the uncertainty in predictions at that sample size.
}
\value{
  A list of class "pmstabilityss" containing:
  \item{overall_stats}{Data frame with summary statistics of interval widths for each sample size.}
  \item{input_data}{Data frame with all calculated variables.}
  \item{plots}{List of ggplot objects.}
  \item{threshold_stats}{Data frame with summary statistics of misclassification probabilities (if threshold specified).}
  \item{cuts_stats}{Data frame with summary statistics by probability categories (if pciwidth and pcutpoints specified).}
}
\author{
  Joie Ensor \email{j.ensor@bham.ac.uk}
}
\references{
  Riley, R. D., Collins, G. S., Whittle, R., Archer, L., Snell, K. I., Dhiman, P., ... & Ensor, J. (2024). Sample size for developing a prediction model with a binary outcome: targeting precise individual risk estimates to improve clinical decisions and fairness. \emph{arXiv preprint arXiv:2407.09293}.
}
\examples{
\dontrun{
# Load example data
data(mtcars)

# Create a binary outcome for demonstration
mtcars$outcome <- ifelse(mtcars$mpg > median(mtcars$mpg), 1, 0)

# Variables for model
vars <- c("wt", "hp", "disp")

# Calculate prevalence
prevalence <- mean(mtcars$outcome)

# Run a logistic regression to get lp
model <- glm(outcome ~ wt + hp + disp, data = mtcars, family = binomial())
mtcars$lp <- predict(model, type = "link")

# Basic usage
result <- pmstabilityss(data = mtcars, 
                        varlist = vars, 
                        prevalence = prevalence,
                        lp = "lp")

# Print results
print(result)

# Plot results
plot(result)

# Using pmsampsize and specifying additional sample sizes
result2 <- pmstabilityss(data = mtcars, 
                         varlist = vars, 
                         prevalence = prevalence,
                         lp = "lp", 
                         pmss = TRUE, 
                         n = c(100, 200, 500))

# With target prediction interval widths
result3 <- pmstabilityss(data = mtcars, 
                         varlist = vars, 
                         prevalence = prevalence,
                         lp = "lp",
                         pcutpoints = c(0.3, 0.6, 1),
                         pciwidth = c(0.1, 0.15, 0.2))

# Looking at classification instability
result4 <- pmstabilityss(data = mtcars, 
                         varlist = vars, 
                         prevalence = prevalence,
                         lp = "lp",
                         threshold = 0.5)

# Plot classification instability
plot(result4, type = "classification")
}
}
\seealso{
  \code{\link[pmsampsize]{pmsampsize}} for sample size calculations based on other criteria.
}
