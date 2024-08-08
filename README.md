# Spatial SIRUS: an explainability algorithm for Spatial Random Forest
Authors: L. Patelli, N. Golini, R. Ignaccolo, M. Cameletti

This repository contains the functions required to use Spatial SIRUS (S-SIRUS) as described in the paper titled "Spatial SIRUS: an explainability algorithm for spatial regression Random Forest" (Patelli et al. 2024+). S-SIRUS is a regression rule extraction algorithm that can be used to explain the random forest (RF) algorithm. In particular, it takes inspiration from the Stable and Interpratble RUle Set (SIRUS) algorithm (Benard et al. 2021ab) and extend it algorithm to the case when data are spatially correlated by using RF-GLS (Saha et al. 2023) instead of the classical RF (Breiman 2001).

Here below ([Guided Example](s.sirus_guided_example.R)) we propose a simplified guided example with some simulated data to run S-SIRUS. All the required functions are available in separated R scripts available in the XXXX folder. Note that currently S-SIRUS can be used only for regression applications in which training and test observations come from the same region.

The following packages are required to run the code: RandomForestsGLS, tidyverse, glmnet, BRISC, sirus.

**References**
- Patelli, L., Cameletti, M., Golini, N., Ignaccolo, R., 2024+. Spatial SIRUS: an explainability algorithm for spatial
regression Random Forest. **(DA COMPLETARE)** [dove si trova](url)
- Saha, A., Basu, S., Datta, A., 2023. Random forests for spatially dependent data. Journal of the American Statistical Association 118, 665–683.
- Bénard, C., Biau, G., da Veiga, S., Scornet, E., 2021a. Interpretable Random Forests via Rule Extraction, in: Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR. pp. 937–945.
- Bénard, C., Biau, G., da Veiga, S., Scornet, E., 2021b. SIRUS: Stable and Interpretable RUle Set for classification. Electronic Journal of Statistics 15, 427–505
- Breiman, L., 2001. Random Forests. Machine Learning 45, 5-32.
