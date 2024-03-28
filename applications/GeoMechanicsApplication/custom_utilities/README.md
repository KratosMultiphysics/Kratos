# Utilities


## Transport Equation utilities

Utilities are developed to calculate matrices and vectors in transport equations

![Image](https://github.com/KratosMultiphysics/Kratos/assets/56549273/296486b0-9e5e-408f-9839-aef8d8c7e720)


### Permeability matrix (H)

The mathematical definition of the permeability matrix is:
$$H = \int_\Omega (\nabla N_p)^T \frac{1}{\mu} k \nabla N_p d\Omega$$
where $\nabla N_p$ is the gradient of the pressure shape function, $\mu$ is the dynamic viscosity (material parameter) and $k$ is material permeability matrix. The k matrix allows one to take into account, for example, an anisotropic permeability. 

### Compressibility matrix (C)

The mathematical definition is:
$$C = \int_\Omega N_{p}^T \frac{1}{Q} N_p d\Omega$$

Where $\Omega$ is the domain, $N_p$ is the pressure shape function and $1/Q$ is the inverse biot modulus.



File transport_equation_utilities.hpp includes 

-  CalculatePermeabilityMatrix function
-  CalculateCompressibilityMatrix function







