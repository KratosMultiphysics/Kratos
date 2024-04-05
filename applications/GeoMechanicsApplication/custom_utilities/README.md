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

Where $\Omega$ is the domain, $N_p$ is the pressure shape function and $1/Q$ is the inverse Biot modulus.



File transport_equation_utilities.hpp includes 

-  CalculatePermeabilityMatrix function
-  CalculateCompressibilityMatrix function

## Stress strain utilities

For convenience functions that compute invariants, equivalents and strain definitions.
Given a stress tensor $$\sigma$$ or a strain tensor $$\epsilon$$. The eigenvalues of the stress tensor are $$\sigma_1 \le \sigma_2 \le \sigma_3$$

$$\sigma = \begin{bmatrix} \sigma_xx & \sigma_xy & sigma_xz \\
                           \sigma_xy & \sigma_yy & sigma_yz \\
                           \sigma_xz & \sigma_yz & sigma_zz  \end{bmatrix$$

$$\epsilon = \begin{bmatrix} \epsilon_xx & \epsilon_xy & epsilon_xz \\
                             \epsilon_xy & \epsilon_yy & epsilon_yz \\
                             \epsilon_xz & \epsilon_yz & epsilon_zz  \end{bmatrix$$

### Trace

The first tensor invariant:

$$ I_1 = trace(\sigma) = \Sigma_i \sigma_{i,i} $$

### Mean stress

$$ p = I_1 = \frac{1}{3} trace(\sigma) $$

### Von Mises stress

With $$J_2$$ the second invariant of the tensor:

$$ \overline\sigma = \sqrt{3 J_2} = \sqrt{0.5((\sigma_{xx}-\sigma_{yy})^2 +
                                 (\sigma_{yy}-\sigma_{zz})^2 +
                                 (\sigma_{zz}-\sigma_{xx})^2 ) +
                            3.0(\sigma_{xy}^2 + \sigma_{yz}^2 + \sigma_{xz}^2) } $$

### Von Mises strain

$$ \overline(\epsilon) = \frac{2}{3} \sqrt{3 J_2} $$

### Green Lagrange strain tensor

With current configuration $$x$$ and reference configuration $$X$$, the deformation gradient $$F = \frac(x)(X)$$ and unit tensor I:

$$ \epsilon = 0.5 ( F^T \cdot F - I ) $$

### Hencky strain tensor

$$\epsilon = 0.5 ln( F^T \cdot F )$$

### Lode angle

The negative sine definition for Lode angle is adapted here [Lode coordinates Wikipedia](https://en.wikipedia.org/wiki/Lode_coordinates):

$$ - sin( 3 \bar{\theta}_s) = \frac{J_3}{2} (\frac{3}{J_2})^{\frac{3}{2}} $$

which brings:

$$ \theta_s = \frac{1}{3} asin( - \frac{27}{2} \frac{(\sigma_1 - p)(\sigma_2 - p)(\sigma_3 - p)}{q^3}) $$

### Mohr Coulomb shear capacity

Assesment of how the current stress utilizes the capacity as defined by the Mohr Coulomb yield surface.

$$ \frac{q}{q_mc} $$

where 

$$q_mc = frac{3}{\sqrt{3}cos \theta_s - sin \theta_s sin \phi}(p sin \phi + c cos \phi) $$

## Mohr Coulomb pressure capacity

Assesment of how the current stress utilizes the capacity as defined by the Mohr Coulomb yield surface.

$$ frac{3 sin \phi}{\sqrt{3}cos \theta_s - sin \theta_s sin \phi} (q_mc - q) $$

