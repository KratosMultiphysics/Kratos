# Contribution Calculators
This folder contains contribution calculators that compute matrix and vector contributions of several physical phenomena. These re-usable classes can be used in elements to compute specific contributions to the local element left-hand side and right-hand side.

Currently, the elements that support the use of these contribution calculators are the `PwElement` and the `UPwInterfaceElement`.

The contributions which have calculators available are:
- Compressibility (generic as well as specifically for filter elements)
- Fluid Body Flow (only computes a right hand side contribution)
- Permeability
- Water $\leftrightarrow$ Displacement Coupling
- Stiffness

## Calculator Interface
Each contribution calculator implements the interface defined in `contribution_calculator.h`, which defines the following methods:
```cpp
virtual std::optional<LHSMatrixType>                           LHSContribution()         = 0;
virtual RHSVectorType                                          RHSContribution()         = 0;
virtual std::pair<std::optional<LHSMatrixType>, RHSVectorType> LocalSystemContribution() = 0;
```
Where the `LHSMatrixType` and `RHSVectorType` are defined as bounded matrices and vectors:
```cpp
using LHSMatrixType = BoundedMatrix<double, NumberOfRows, NumberOfColumns>;
using RHSVectorType = BoundedVector<double, NumberOfRows>;
```
Each calculator computes right-hand side contributions, while the left-hand side contributions are optional, as not all contributions require a left-hand side. The calling side is responsible for checking whether the left-hand side contribution is available (i.e. whether the `std::optional` has a value).

## General Usage
All calculators define a calculator-specific `InputProvider`, which is used to provide functions for retrieving the required input data. The calling side is responsible for creating this provider, to supply the required data to the calculator. The provider is passed to the calculator at construction time.

## Available Calculators
The available contributions are listed in the introduction. The following sections provide more details on each of them. The documentation is still work in progress, meaning not all calculators have a detailed description yet.

### Water $\leftrightarrow$ Displacement Coupling
The [UPCouplingCalculator](up_coupling_calculator.hpp) and the [PUCouplingCalculator](pu_coupling_calculator.hpp) compute the coupling matrix and terms ($Q$) to account for the between the water pressure and the displacement degrees of freedom. The mathematical definition is:
$$Q = \int_\Omega B^T \alpha \xi m N_p d\Omega$$
where $B$ is the B-matrix, $\alpha$ is the Biot-alpha (relation between pressure and displacements, material parameter), $\xi$ is the Bishop coefficient, $m$ is the Voigt-vector (for example $[1,1,1,0]$ for plane strain and $[1,1,1,0,0,0]$ for 3D cases) and $N_p$ is the pressure shape function.

The `UPCouplingCalculator` computes the influence of the pressure on the displacement degrees of freedom and needs the following input data:
- The B-matrices for all integration points
- The Biot-alpha values for all integration points
- The Bishop coefficients for all integration points
- The pressure shape functions
- The integration coefficients for all integration points
- The Voigt-vector
- The nodal water pressures

This data is provided via the `InputProvider` and enables the calculator to compute left-hand side contribution ($Q$ as defined before) and right-hand side contributions ($Q p$, in which $p$ is the vector of nodal water pressures).

The `PUCouplingCalculator` computes the influence of the displacements on the pressure degrees of freedom and needs the following input data:
- The B-matrices for all integration points
- The Biot-alpha values for all integration points
- The degrees of saturation for all integration points
- The pressure shape functions
- The integration coefficients for all integration points
- The Voigt-vector
- The nodal velocity values

This data is provided via the `InputProvider` and enables the calculator to compute left-hand side contribution. The PU coupling matrix is related to Q as follows: $$Q_{pu} = Q^{T}$$ This results in left-hand side contribution. **_Note: next to the transpose $Q_{pu}$ is calculated using degrees of saturation instead of Bishop coefficients._**

($Q_{pu}$ as defined before) and right-hand side contributions ($Q_{pu} v$, in which $v$ is the vector of nodal velocities).

For more code-specific information on how to call/use these calculators in your code, please refer to the related [unit tests](../../tests/cpp_tests/custom_elements/contribution_calculators/test_coupling_calculators.cpp)
