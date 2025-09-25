# PTCMGH R package

## Promotion Time Cure Models with a General Hazard structure

The `PTCMGH` R package implements promotion time cure models of the type:

$S_c(t) = \exp [ -\theta \tilde{F}(t) ] , \quad t \geq 0,$

If a vector of covariates ${\bf x} \in {\mathbb R}^p$ are available, these can be incorporated into the two different components of the PTCM.

1. Let ${\bf w} \subseteq {\bf x}$, and consider the log-link $\theta({\bf w}) = \exp({\bf w}^{\top}{\boldsymbol{\alpha}})$.

2. Let ${\bf z} \subseteq {\bf x}$, and consider the hazard structure $\tilde{H}(t \mid {\bf z})$. Then, we can define

$\tilde{F}(t \mid {\bf z}) = 1 - \exp[-\tilde{H}(t \mid {\bf z})].$

The `PTCMGH` R package provides functions for simulating survival times from PTCMs with proportional hazards, accelerated failure time, and general hazard structures under a log-link for $\theta$, as well as for performing maximum likelihood estimation of these models. The hazard structure used for modelling $\tilde{F}$ can include the proportional hazards, the accelerated failure time, or more [general hazard structures](https://doi.org/10.1177/0962280218782293).

- General Hazard (GH) model.

- Accelerated Failure Time (AFT) model.

- Proportional Hazards (PH) model.

- Accelerated Hazards (AH) model.


These models are fitted using the R commands `nlminb` and `optim`. Thus, the user needs to specify the initial points and to check the convergence of the optimisation step, as usual.


The current version of the `PTCMGH` R package implements the following parametric baseline hazards for the models discussed in the previous section, using the command `PTCMMLE`.

- [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) (PGW) distribution. 
 
- [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) (EW) distribution. 
 
- [Generalised Gamma](http://rpubs.com/FJRubio/GG) (GenGamma) distribuiton. 

- [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) (Gamma) distribution. 

- [Lognormal](https://en.wikipedia.org/wiki/Log-normal_distribution) (LogNormal) distribution. 

- [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution) (LogLogistic) distribution. 

- [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) (Weibull) distribution. (only for AFT, PH, and AH models) 


All positive parameters are transformed into the real line using a `log` link (reparameterisation).

The package depends on functions from the [`HazReg` R package](https://github.com/FJRubio67/HazReg), which must also be installed and loaded.

```
library(devtools)
install_github("FJRubio67/HazReg")
install_github("FJRubio67/PTCMGH")

library(HazReg)
library(PTCMGH)
?simPTCMGH
?PTCMMLE
```


### See also: 
- [`PTCMGH` R package: Promotion Time Cure Models with a General Hazard structure](https://rpubs.com/FJRubio/PTCMGH)
