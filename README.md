# Bayesian Density Based Clustering

This folder contains R scripts necessary to reproduce all results featured in the Bayesian Density Based Clustering article by David Buch, Miheer Dewaskar, and David Dunson. Results can be organized as follows:

1. An illustrative univariate toy data example, where data are generated from a 
mixture of a uniform distribution and a Gaussian distribution, centered sufficiently
far apart so that there is a 'low density' valley between two modes. These data are 
fit with a Dirichlet Process mixture of Gaussians model, and it is shown that 
traditional model-based clustering estimates which assign observations to mixture
components lead to 'splitting' observations which ought to be attributed to the
uniform component. Using the same flawed, fitted density model and applying our
density-based clustering approach, we obtain a clustering of the observations which more closely aligns with the intuitive notion of clusters as regions of high density separated by regions of low density.
2. A collection of three toy datasets which exhibit certain pathological characteristics of model-based clustering. For two of them, while the standard Dirichlet Process Mixture of Gaussians approach fails to identify an appropriate clustering structure, our density based clustering method extracts reasonable estimates of the appropriate clustering structures based on the density model. We find we obtain essential the same inferred clustering structure when we replace the Dirichlet Process Mixture of Gaussians with other density models. For the third example (concentric rings) the Dirichlet Process Mixture of Gaussian model is simply too poorly adapted to the dataset to represent the data generating process well, so density-based clustering based on those density estimates fail to perform as expected. However, the method does perform as expected on that dataset when alternative density models are used, such as the Nearest Neighbor Dirichlet Mixture of Chattopadhyay et al. 2023.
3. A real data analysis of Sky Survey data, as described in the article.
4. A simulation study analyzing data meant to emulate the Sky Survey data.


