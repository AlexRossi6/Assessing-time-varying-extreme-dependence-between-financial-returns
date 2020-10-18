# Assessing-time-varying-extreme-dependence-between-financial-returns (Bachelor's thesis)
This project provides a method to assess time-varying extreme dependence between financial returns. The main purpose is to state whether a change in the extremal dependence of international stock markets has - or has not- been influenced by important events, which can be political or economic.

## Description

In this project, the focus is on extreme losses since they are the main factor analyzed in risk management. The two extreme cases of extremal dependence structure are asymptotic independence and asymptotic dependence; the last-mentioned situation is the framework of our modeling; indeed, no models for the spectral density of extremes will provide useful information on the extremal dependence structure in a situation of asymptotic independence.
To compare different stock markets, we necessarily deal with bivariate extremes. Therefore, the method proposed is applied and tested in a bivariate space, but can also be extended in a multivariate space. The results are based on Extreme value theory, that is the theory in support of the study of asymptotic distributions of rare events and the basis of risk management, not only in the financial field.
To make inference on extremes, a non-parametric model based on Bernstein polynomials has been performed, in particular, to model the so-called spectral density, which under specific conditions, can be interpreted as the limit distribution of the pseudo-angles of the extremes. Extreme dependence structure is well described by the spectral density and to assess changes in this structure, a change-point analysis has been implemented in the Bernstein polynomial based model. The algorithm is able to detect a single change-point in a set of time-varying observations by computing the likelihood for different change-points and by taking the most likely as hypothesized change-point location.
The efficiency and the efficacy of the proposed model have been analyzed with simulations of random bivariate extreme values. The model is then applied to real stock markets, and three different results are explored. Indeed, the change-point analysis estimated a situation of stationary extremal dependence structure (no change-point), and two situations of nonstationary extremal dependence structure; the first estimated one change-point, and the second one two change-points.

## Simulation

To state the efficacy of the change-point analysis, which works simultaneously with the Bernstein polynomial based model, two types of analyses were performed. In the first analysis, the observations were simulated from a single distribution to see how the code works in a situation in which the distribution doesn’t change. The second one consists in simulating random bivariate extreme values from two distributions with different parameters of dependence, compute their pseudo-angles, and, after having combined them by row, the code was run to state the performance of the analysis in finding the change-point.
