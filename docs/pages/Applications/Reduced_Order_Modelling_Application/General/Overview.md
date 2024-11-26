---
title: Overview
keywords:
tags: [Overview.md]
sidebar: rom_application
summary:
---

# Theory

The Kratos RomApplication implements Projection-based Reduced Order Models (PROMs)[1,2]. PROMs are a family of reduced models that aim to accelerate the evaluation of parametric models, by incurring a fraction of the costs associated to the high-dimensional FOMs, while still taking into account the physics underlying the models at hand. PROMs are comprised of two different stages:

- Offline stage : In this stage, a set of simulations is performed using the computationally expensive FOM, and the resulting solutions are stored in a so-called snapshots matrix. This matrix is then processed to obtain a reduced-space where the discrete equations are projected (therefore the name of the method). Moreover,  we accomplish the decoupling of the ROMs from full dimensional variables through a hyper-reduction mesh sampling and weighting hyper-reduction technique.

- Online stage: With the basis and additional hyper-reduction data available, the hyper-reduced order models (HROMs) can be efficiently launched for unexplored parameters at a fraction of the cost associated with the FOMs.


[1] Hesthaven, J. S., Rozza, G., & Stamm, B. (2016). Certified reduced basis methods for parametrized partial differential equations (Vol. 590, pp. 1-131). Berlin: Springer.

[2] Rozza, G., Stabile, G., & Ballarin, F. (Eds.). (2022). Advanced reduced order methods and applications in computational fluid dynamics. Society for Industrial and Applied Mathematics.




# License

The RomApplication is OPEN SOURCE. The main code and program structure is available and aimed to grow with the need of any user willing to expand it. The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

# Contact

* **Riccardo Rossi** - *Group Leader* - [rrossi@cimne.upc.edu](mailto:rrossi@cimne.upc.edu)
* **Raul Bravo** - *Developer* - [jrbravo@cimne.upc.edu](mailto:jrbravo@cimne.upc.edu)
* **Sebastian Ares de Parga** - *Developer* - [sares@cimne.upc.edu](mailto:sares@cimne.upc.edu)
* **Nicolas Sibuet** - *Developer* - [nsibuet@cimne.upc.edu](mailto:nsibuet@cimne.upc.edu)