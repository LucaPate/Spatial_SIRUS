# Spatial SIRUS: an explainability algorithm for Spatial Random Forest
Authors: L. Patelli, N. Golini, R. Ignaccolo, M. Cameletti

This repository contains the functions required to use Spatial SIRUS (S-SIRUS) as described in the paper titled "Spatial SIRUS: an explainability algorithm for spatial regression Random Forest" available here (mettere LINK ARXIV). S-SIRUS is a regression rule extraction algorithm that can be used to explain the random forest (RF) algorithm. In particular, it takes inspiration from the Stable and Interpratble RUle Set (SIRUS) algorithm () and extend it algorithm to the case when data are spatially correlated by using RF-GLS (link?) instead of the classical RF.

*** METTERE DIRETTAMENTE I LINK AI RELATIVI PAPER CHE VORREMMO CITARE?

Here below (SE SI Può COPIARE SOTTO OPPURE MENZIONARE IL FILE) we propose a simplified guided example with some simulated data to run S-SIRUS. All the required functions are available in separated R scripts available in the XXXX folder. Note that currently S-SIRUS can be used only for regression applications in which training and test observations come from the same region.

The following packages are required to run the code: RandomForestsGLS, tidyverse, glmnet, BRISC, sirus.
