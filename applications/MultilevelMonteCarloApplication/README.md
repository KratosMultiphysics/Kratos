# Multilevel Monte Carlo Application

MultilevelMonteCarloApplication provides different algorithms, belonging to the Monte Carlo (MC) family, to estimate statistics of scalar and field quantities of interest.
The application is designed for running on distributed and high performance computing systems, exploiting both OpenMP and MPI parallel strategies.
The application contains several interfaces with external libraries.

## Getting started

This application is part of the Kratos Multiphysics Platform. Instructions on how to download, install and run the software in your local machine for development and testing purposes are available for both [Linux](https://github.com/KratosMultiphysics/Kratos/wiki/Linux-Build) and [Windows](https://github.com/KratosMultiphysics/Kratos/wiki/Windows-Install) systems.

### Prerequisites

Build [Kratos](https://github.com/KratosMultiphysics/Kratos/wiki) and make sure to have
``` cmake
add_app ${KRATOS_APP_DIR}/MultilevelMonteCarloApplication
```
in the compilation configuration, in order to compile the `MultilevelMonteCarloApplication` application.

## Hierarchical Monte Carlo methods
* Repeatedly generate the random input and solve the associated deterministic problem.
* Convergence to the exact statistics as the number of realizations grows.
* Problem under consideration considered as a black-box.
* Convergence rate independent from stochastic space dimension.

### Monte Carlo

* Monte Carlo (MC) is the reference method in the stochastic analysis of multiphysics problems with uncertainties in the data parameters.
* Levels of parallelism:
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Deterministic,
    * Adaptive.

### Multilevel Monte Carlo

* Multilevel Monte Carlo (MLMC) requires a hierarchy of levels with increasing accuracy to solve the statistical problem. Convergence rate is faster with respect to MC if MLMC hypotheses are satisfied.
* Computation of a large number of cheap and lower accuracy realizations, while only few expensive high accuracy realizations are run.
    * Low accuracy levels: capture statistical variability,
    * High accuracy levels: capture discretization error.
* Example of hierarchy of computational grids, showing increasing accuracy levels (by decreasing mesh size):
<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Readme_files/MultilevelMonteCarloApplication/mesh012.PNG" alt="Solution" style="width: 600px;"/>
</p>

* Levels of parallelism:
    * Between levels,
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Deterministic,
    * Adaptive.

### Continuation Multilevel Monte Carlo

* A set of decreasing tolerances is used and updated on the fly to adaptively estimate the hierarchy and run MLMC.
* Levels of parallelism:
    * Between levels,
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Adaptive.

### Asynchronous Monte Carlo

* This algorithm is equivalent to MC, but designed for running in distributed environments. It avoids idle times and keeps at maximum the computational efficiency.
* Levels of parallelism:
    * Between batches,
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Deterministic.

### Asynchronous Multilevel Monte Carlo

* This algorithm is equivalent to MLMC, but designed for running in distributed environments. It avoids idle times and keeps at maximum the computational efficiency.
* Levels of parallelism:
    * Between batches,
    * Between  levels,
    * Between samples,
    * On each sample at solver level.
* Hierarchy update:
    * Deterministic.


## Statistical tools

### Power sums

* Update on the fly of power sums.
* A power sum of order p is defined as: <p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=S_p^N:=\sum_{i=1}^{N}(Q(w^{(i)}))^p"> .
</p>

### h-statistics

* The h-statistic of order p is the unbiased estimator with minimal variance of the central moment of order p.
* h-statistic dependencies are <p align="center">
  <img src="https://render.githubusercontent.com/render/math?math=h_P:=f(S_p,N),p\in[1,P]"> .
</p>


## Convergence criteria

* Convergence is achieved if the estimator of interest reaches a desired tolerance with respect to the true estimator with a given confidence.
* The failure probability to satisfy (for expected value and MC) is
<p align="center">
    <img src="https://render.githubusercontent.com/render/math?math=\mathbb{P}[\left|\mathbb{E}^{MC}[Q_H]-\mathbb{E}[Q]\right|>\varepsilon]<\phi"> .
</p>

* Other convergence criteria available:
    * Mean square error,
    * Sample variance criteria (MC only),
    * Higher order (up to the fourth) moments criteria (MC only).



## Hierarchy

* Hierarchy strategies:
   * stochastic adaptive refinement: hierarchy of levels built refining in space, performing solution-oriented adaptive space refinement. The coarsest mesh is shared between all realizations, and for each realization different meshes are generated, accordingly to the random variable. Requires compiling `MESHING_APPLICATION`.
   * deterministic adaptive refinement: hierarchy of levels built refining in space, performing solution-oriented adaptive space refinement. All meshes are shared between all realizations, and adaptive refinement is done at pre-process, exploiting a user-defined random variable. Requires compiling `MESHING_APPLICATION`.
   * standard: the user takes care of building the hierarchy, using the strategy he prefers (such as uniform refinement).

* Metric strategies:
    * geometric error estimate: the analysis of the hessian of the numerical solution controls the mesh refinement.
    * divergence-free error estimate: the analysis of the mass conservation controls the mesh refinement (suitable only for CFD cases). Requires compiling `EXAQUTE_SANDBOX_APPLICATION`. *In progress*.


## External Libraries

MultilevelMonteCarloApplication makes use of third part libraries.
Information about these libraries can be found in their respective pages, which are listed below.

### XMC

[XMC](https://gitlab.com/RiccardoRossi/exaqute-xmc) is a Python library, with BSD 4 license, designed for hierarchical Monte Carlo methods. The library develops the above-mentioned algorithms, statistical tools and convergence criteria. The library presents a natural integration with Kratos, which is XMC default solver. By default, an internal version of the library is used. If one wants to use an external version of the library, the flag `-DUSING_INTERNAL_XMC=OFF \ ` should be added when compiling Kratos.

### PyCOMPSs

PyCOMPSs is the Python library required in order to use task-based programming software [COMPSs](https://www.bsc.es/research-and-development/software-and-apps/software-list/comp-superscalar) in a Python environment.
By default PyCOMPSs is not required in order to run the application.
In case one wants to run using this library, Kratos needs to be compiled adding the flag `-DUSING_PYCOMPSS=ON \ `.
The current version is able to run several thousands of samples at once exploiting PyCOMPSs in distributed systems, maximizing parallelism and computational efficiency. Optimal strong scalability up to 128 working nodes (6144 CPUs) has been demonstrated.

Instructions for the installation can be found in the [Kratos wiki](https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs).
To run with `runcompss`, the environment variable `EXAQUTE_BACKEND=pycompss` must be set to use the distributed computing capabilities.
Additionally, running with `runcompss` requires to add to the `PYTHONPATH` the path of the XMC library, that is `/path/to/Kratos/applications/MultilevelMonteCarloApplication/external_libraries/XMC`. You can add the library to the `PYTHONPATH` either in the `.bashrc` file, or directly when running the code using the `runcompss` key `--pythonpath`. We refer to the [Kratos wiki](https://github.com/KratosMultiphysics/Kratos/wiki/How-to-run-multiple-cases-using-PyCOMPSs#running-with-pycompss) for details.

### Mmg and ParMmg

[Mmg](https://www.mmgtools.org/) is an open source software for simplicial remeshing. It provides 3 applications and 4 libraries.
Instructions for installing Mmg can be found in the [Kratos wiki](https://github.com/KratosMultiphysics/Kratos/wiki/%5BUtilities%5D-MMG-Process).
ParMmg is the MPI parallel version of the remeshing library Mmg. Instructions for installing ParMmg can be found in the [Kratos wiki](https://github.com/KratosMultiphysics/Kratos/wiki/%5BUtilities%5D-ParMmg-Process).

## License

The MultilevelMonteCarloApplication is OPEN SOURCE. The main code and program structure are available and aimed to grow with the need of any user willing to expand it. The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction, but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## Examples

Many examples can be found in the [Kratos Multiphysics Examples repository](https://github.com/KratosMultiphysics/Examples/tree/master/multilevel_monte_carlo).

## Main References
- Amela, R., Ayoul-Guilmard, Q., Badia, R. M., Ganesh, S., Nobile, F., Rossi, R., & Tosi, R. (2019). ExaQUte XMC. https://doi.org/10.5281/zenodo.3235833
- Krumscheid, S., Nobile, F., & Pisaroni, M. (2020). Quantifying uncertain system outputs via the multilevel Monte Carlo method — Part I: Central moment estimation. Journal of Computational Physics. https://doi.org/10.1016/j.jcp.2020.109466

## Contact

* **Riccardo Rossi** - *Group Leader* - [rrossi@cimne.upc.edu](mailto:rrossi@cimne.upc.edu)
* **Riccardo Tosi** - *Developer* - [rtosi@cimne.upc.edu](mailto:rtosi@cimne.upc.edu)
* **Marc Núñez** - *Developer* - [mnunez@cimne.upc.edu](mailto:mnunez@cimne.upc.edu)
* **Ramon Amela** - *Developer* - [ramon.amela@bsc.es](mailto:ramon.amela@bsc.es)
