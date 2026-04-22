# PTCMGH: Promotion Time Cure Models with a General Hazard Structure

[![R package](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

`PTCMGH` is an R package for fitting and simulating from **Promotion Time Cure
Models (PTCMs)** with flexible hazard structures. PTCMs are a class of survival
models that explicitly account for a cured fraction of the population — individuals
who will never experience the event of interest. The package is built around the
**General Hazard (GH)** structure, which nests the most widely used hazard
regression models as special cases, giving users flexibility in how covariates
enter both the cure and the survival components of the model.

The population survival function takes the form:

$$S_c(t) = \exp[-\theta  \tilde{F}(t)], \quad t \geq 0,$$

where $\theta > 0$ governs the cured fraction and $\tilde{F}(t)$ is a proper
cumulative distribution function. Covariates $\mathbf{x} \in \mathbb{R}^p$ can
enter through either or both components:

1. **Cure component:** via the log-link $\theta(\mathbf{w}) = \exp(\mathbf{w}^\top \boldsymbol{\alpha})$, where $\mathbf{w} \subseteq \mathbf{x}$.
2. **Survival component:** via the hazard structure $\tilde{H}(t \mid \mathbf{z})$, giving $\tilde{F}(t \mid \mathbf{z}) = 1 - \exp[-\tilde{H}(t \mid \mathbf{z})]$, where $\mathbf{z} \subseteq \mathbf{x}$.

## Installation

`PTCMGH` depends on [`HazReg`](https://github.com/FJRubio67/HazReg), which must
be installed first:

```r
# install.packages("devtools")
devtools::install_github("FJRubio67/HazReg")
devtools::install_github("FJRubio67/PTCMGH")

library(HazReg)
library(PTCMGH)
```

## Main functions

| Function | Description |
|---|---|
| `PTCMMLE` | Maximum likelihood estimation for PTCMs |
| `simPTCMGH` | Simulation of survival times from PTCMs |

For full documentation: `?PTCMMLE`, `?simPTCMGH`

A PDF manual is also available in the repository: [`PTCMGH.pdf`](PTCMGH.pdf).

## Supported hazard structures

The following hazard models are available for the survival component $\tilde{F}$:

| Model | Abbreviation |
|---|---|
| General Hazard | GH |
| Proportional Hazards | PH |
| Accelerated Failure Time | AFT |
| Accelerated Hazards | AH |

Models are fitted by maximum likelihood via `nlminb` and `optim`. Users should
specify initial values and verify convergence of the optimisation, as is standard
practice for these routines.

## Supported baseline hazard distributions

All positive parameters are log-transformed for unconstrained optimisation.

| Distribution | Key | GH | PH | AFT | AH |
|---|---|:---:|:---:|:---:|:---:|
| [Power Generalised Weibull](http://rpubs.com/FJRubio/PGW) | PGW | ✓ | ✓ | ✓ | ✓ |
| [Exponentiated Weibull](http://rpubs.com/FJRubio/EWD) | EW | ✓ | ✓ | ✓ | ✓ |
| [Generalised Gamma](http://rpubs.com/FJRubio/GG) | GenGamma | ✓ | ✓ | ✓ | ✓ |
| [Gamma](https://en.wikipedia.org/wiki/Gamma_distribution) | Gamma | ✓ | ✓ | ✓ | ✓ |
| [Log-normal](https://en.wikipedia.org/wiki/Log-normal_distribution) | LogNormal | ✓ | ✓ | ✓ | ✓ |
| [Log-logistic](https://en.wikipedia.org/wiki/Log-logistic_distribution) | LogLogistic | ✓ | ✓ | ✓ | ✓ |
| [Weibull](https://en.wikipedia.org/wiki/Weibull_distribution) | Weibull | — | ✓ | ✓ | ✓ |

## Tutorial

- [`PTCMGH`: Promotion Time Cure Models with a General Hazard structure](https://rpubs.com/FJRubio/PTCMGH) — illustrative example on RPubs
- A rendered tutorial is also available in the [`Tutorial/`](Tutorial/) folder of this repository.

## Related resources

- [`HazReg`](https://github.com/FJRubio67/HazReg) — R package for parametric
  hazard-based regression models (required dependency)
- [`HazReg.jl`](https://github.com/FJRubio67/HazReg.jl) — Julia implementation
- [Short course on Parametric Survival Analysis](https://github.com/FJRubio67/ShortCourseParamSurvival)

## License

This package is licensed under the [MIT License](LICENSE).
