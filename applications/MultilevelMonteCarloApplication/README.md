# Multilevel Monte Carlo Application

MultilevelMonteCarloApplication provides different algorithms, belonging to the Monte Carlo (MC) family, and tools to perform statistical analysis of stochastic problems.
The application contains several interfaces to both Kratos and third part libraries.

## Getting started

This application is part of the Kratos Multiphysics Platform. Instructions on how to download, install and run the software in your local machine for development and testing purposes are available for both [Linux](https://github.com/KratosMultiphysics/Kratos/wiki/Linux-Build) and [Windows](https://github.com/KratosMultiphysics/Kratos/wiki/Windows-Install) systems.

### Prerequisites

Build [Kratos](https://github.com/KratosMultiphysics/Kratos/wiki) and make sure to have

``` cmake
-DMULTILEVEL_MONTE_CARLO_APPLICATION=ON
-DMESHING_APPLICATION=ON
-DEXAQUTE_SANDBOX_APPLICATION=ON
```

in the compilation configuration, in order to compile the MultilevelMonteCarloApplication along with other required auxiliary Kratos applications.

Additionally, you need to add
``` cmake
export PYTHONPATH=$PYTHONPATH:/path/to/Kratos/bin/Release/KratosMultiphysics/MultilevelMonteCarloApplication
```
to use the local (empty) PyCOMPSs libraries.

## Algorithms

### Monte Carlo

* Monte Carlo (MC) is the reference method in the stochastic analysis of multiphysics problems with uncertainties in the data parameters.
* The idea is to repeatedly generate the random input and to solve numerically the associated deterministic problem, in order to perform a statistical analysis.
* Convergence to the exact statistics as the number of samples grows.
* Level of parallelism:
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Deterministic,
    * Adaptive.
* Advantages:
    * Problem under consideration considered as a black-box.
    * Does not suffer of the curse of dimensionality.
* Disadvantages:
    * Convergence rate of the mean square error <img src="https://latex.codecogs.com/svg.latex?\sim~O(N^{-0.5})" alt="Solution"/>, thus may be too high for complex problems.
    This leads to the development of other algorithms, such as Multilevel Monte Carlo (MLMC).
* MC estimator of the expectation of a Quantity of Interest (QoI): <p align="center">
  <img src="https://latex.codecogs.com/svg.latex?\mathbb{E}^{MC}[QoI]:=\frac{\sum_{i=1}^{N}QoI(w^{(i)})}{N}" alt="Solution" />
</p>

### Multilevel Monte Carlo

* Multilevel Monte Carlo (MLMC) consists in the simultaneous computation of MC QoI samples on a hierarchy of levels with increasing accuracy.
* Computation of a large number of cheap and lower accuracy QoI samples with few expensive high accuracy samples.
    * Low accuracy levels: capture satistical variability,
    * High acuracy levels: capture bias.
* Example of hierarchy of computational grids, showing increasing accuracy levels (by decreasing mesh size):
<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Readme_files/MultilevelMonteCarloApplication/mesh012.PNG" alt="Solution" style="width: 600px;"/>
</p>

* Level of parallelism:
    * Between levels,
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Deterministic,
    * Adaptive.

### Continuation Multilevel Monte Carlo

* Evolution of the MLMC algorithm: set of decreasing tolerances to reach gradually the final one. Updating on the fly the tolerance hierarchy.
* Hierarchy update:
    * Adaptive.

### Asynchronous Monte Carlo

* This algorithm fills the machine when running the problem in distributed environments, avoiding idle times and keeping at maximum the computational efficiency.
* Level of parallelism:
    * Between batches,
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Deterministic.

### Asynchronous Multilevel Monte Carlo

* This algorithm fills the machine when running the problem in distributed environments, avoiding idle times and keeping at maximum the computational efficiency.
* Level of parallelism:
    * Between batches,
    * Betweenn levels,
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Deterministic.


## Statistical tools

### Power sums

* Update on the fly of the power sums.
* A power sum of order p is defined as: <p align="center">
  <img src="https://latex.codecogs.com/svg.latex?S_p^N:=\sum_{i=1}^{N}(QoI(w^{(i)}))^p" alt="Solution" />
</p>

### h-statistics

* The h-statistic of order p is the unbiased estimator with minimal variance of the central moment of order p.
* h-statistic is computed as: <p align="center">
  <img src="https://latex.codecogs.com/svg.latex?h_P:=h_P(S_p,N),p\in[1,P]" alt="Solution"/>
</p>


## Convergence criteria

* Convergence is achieved if the estimator of the QoI reaches a desired tolerance with respect ot the true estimator with a given confidence.
* The failure probability to satisfy is
<p align="center">
  <img src="https://latex.codecogs.com/svg.latex?\mathbb{P}[\abs{\mathbb{E}[QoI_M]-\mathbb{E}[QoI]}>\varepsilon]<\phi" alt="Solution"/>
</p>

* Available convergence criteria:
    * Sample variance criteria,
    * Higher order (up to the fourth) moments criteria,
    * Total error criteria,
    * Relative total error criteria.


## Adaptive refinement

The choice is to build the hierarchy of levels for MLMC refining in space, performing solution-oriented adaptive space refinement.
* Metric strategies:
    * geometric error estimate: computation of the hessian of the numerical solution, which gives
information about where the mesh needs to be refined more,
    * divergence-free error estimate: the analysis of the mass conservation controls the mesh
refinement (suitable only for CFD cases).
* Refinement approaches:
    * concurrent adaptive refinement: mesh generation and simulation running one after the other for each simulation,
    * single refinement: storage of the mesh for each accuracy level.

The [MeshingApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MeshingApplication) and the [ExaquteSandboxApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/ExaquteSandboxApplication) are exploited to achieve the first and the second desired metric, respectively.


## External Libraries

MultilevelMonteCarloApplication makes use of third part libraries.
Information about these libraries can be found in their respective pages, which are listed below:

### MMG

[MMG](https://www.mmgtools.org/) is an open source software for simplicial remeshing. It provides 3 applications and 4 libraries.
Informations for installing MMG can be found in the [Kratos wiki](https://github.com/KratosMultiphysics/Kratos/wiki/%5BUtilities%5D-MMG-Process).

### PyCOMPSs

PyCOMPSs is the python library required in order to use [COMPSs](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar) in a python environment.
By default PyCOMPSs is not required in order to run the application.
In case you want to run using this library, you will need to remove
``` cmake
export PYTHONPATH=$PYTHONPATH:/path/to/Kratos/bin/Release/KratosMultiphysics/MultilevelMonteCarloApplication
```
since you need to use the path given by the installation.

The instructions for the installation can be found in the [Kratos wiki](https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs). The current version is able to run several thousands of samples at once exploiting PyCOMPSs.

Finally, in the files [mc_utilities.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MultilevelMonteCarloApplication/python_scripts/mc_utilities.py), [mlmc_utilities.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MultilevelMonteCarloApplication/python_scripts/mlmc_utilities.py) and [statistical_variable_utilities.py](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MultilevelMonteCarloApplication/python_scripts/statistical_variable_utilities.py) you need to switch to:
``` cmake
# Import PyCOMPSs
from exaqute.ExaquteTaskPyCOMPSs import *   # to execute with runcompss
# from exaqute.ExaquteTaskHyperLoom import *  # to execute with the IT4 scheduler
# from exaqute.ExaquteTaskLocal import *      # to execute with python3
```
to use the distributed computing capabilities.

## License

The MultilelvelMonteCarloApplication is OPEN SOURCE. The main code and program structure is available and aimed to grow with the need of any user willing to expand it. The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## References
- Dadvand, P., Rossi, R., & Oñate, E. (2010). An object-oriented environment for developing finite element codes for multi-disciplinary applications. *Archives of Computational Methods in Engineering*, 17(3), 253–297.
- Amela, R., Ramon-Cortes, C., Ejarque, J., Conejero, J., & Badia, R. M. (2018). Executing linear algebra kernels in heterogeneous distributed infrastructures with PyCOMPSs. *Oil & Gas Science and Technology--Revue d’IFP Energies Nouvelles*, 73, 47.
- Pisaroni, M., Nobile, F., & Leyland, P. (2017). A Continuation Multi Level Monte Carlo (C-MLMC) method for uncertainty quantification in compressible inviscid aerodynamics. *Computer Methods in Applied Mechanics and Engineering*, 326, 20–50.
- Pisaroni, M., Krumscheid, S., & Nobile, F. (2017). Quantifying uncertain system outputs via the multilevel Monte Carlo method - Part I: Central moment estimation. *Retrieved from MATHICSE Technical report 23*.2017.
- C. Bayer, H. Hoel, E. von Schwerin, R. Tempone; On NonAsymptotyc optimal stopping criteria in Monte Carlo simulations; *SIAM Journal on Scientific Computing*, 2014, Vol. 36, No. 2 : pp. A869-A885
- P. Pébay, T. B. Terriberry, H. Kolla, J. Bennett; Stable, Scalable Formulas for Parallel and Online Computation of Higher-order Multivariate Central Moments with Arbitrary Weights; *Computational Statistics*, 2016, 31:1305-1325
- Dapogny, C., Dobrzynski, C., & Frey, P. (2014). Three-dimensional adaptive domain remeshing, implicit domain meshing, and applications to free and moving boundary problems. *Journal of Computational Physics*, 262, 358–378.

## Contact

* **Riccardo Rossi** - *Group Leader* - [rrossi@cimne.upc.edu](mailto:rrossi@cimne.upc.edu)
* **Riccardo Tosi** - *Developer* - [rtosi@cimne.upc.edu](mailto:rtosi@cimne.upc.edu)
* **Ramon Amela** - *Developer* - [ramon.amela@bsc.es](mailto:ramon.amela@bsc.es)
