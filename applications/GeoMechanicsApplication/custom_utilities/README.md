# Utilities


## Transport Equation utilities

Utilities are developed to calculate matrices and vectors in transport equations

$$ \begin{bmatrix} M & 0 \\
                   0 & 0 \end{bmatrix} \begin{bmatrix} \ddot{u} \\
                                                       \ddot{p} \end{bmatrix}  +
   \begin{bmatrix} D & 0 \\
                   Q^T & C \end{bmatrix} \begin{bmatrix} \dot{u} \\
                                                          \dot{p} \end{bmatrix}  +
   \begin{bmatrix} K & -Q \\
                   0 & H \end{bmatrix} \begin{bmatrix} u \\
                                                       p \end{bmatrix} =
   \begin{bmatrix} f_u \\
                   f_p \end{bmatrix} $$

### Permeability matrix (H)

The mathematical definition of the permeability matrix is:
$$H = \int_\Omega (\nabla N_p)^T \frac{1}{\mu} k \nabla N_p d\Omega$$
where $\Omega$ is the domain, $\nabla N_p$ is the gradient of the pressure shape function, $\mu$ is the dynamic viscosity (material parameter), and $k$ is the material permeability matrix. The k matrix allows one to take into account, for example, an anisotropic permeability. 

### Compressibility matrix (C)

The mathematical definition is:
$$C = \int_\Omega N_{p}^T \frac{1}{Q} N_p d\Omega$$
where $N_p$ is the pressure shape function and $1/Q$ is the inverse Biot modulus.

### Coupling Matrix (Q)

The mathematical definition is:
$$Q = \int_\Omega B^T \alpha \xi m N_p d\Omega$$
where $B$ is the B-matrix, $\alpha$ is the Biot-alpha (relation between pressure and displacements, material parameter), $\xi$ is the Bishop coefficient, $m$ is the Voigt-vector ($[1,1,1,0,0,0]$) and $N_p$ is the pressure shape function.

### Soil density

The soil density is calculated as
$$\rho= S_r n \rho_w + (1 - n ) \rho_s$$
where $S_r$ is the degree of saturation, $n$ is the porosity, $\rho_w$ is the water density, $\rho_s$ is the solid density. 

File transport_equation_utilities.hpp includes 

-  CalculatePermeabilityMatrix function
-  CalculateCompressibilityMatrix function
-  CalculateCouplingMatrix function
-  CalculateSoilDensity function
-  CalculateSoilDensities function that calculates solid density for all integration points

## Equation of motion utilities

### Mass Matrix (M)

The mathematical definition is:
$$M = \int_\Omega N_{u}^T \rho N_u d\Omega$$

Where $\Omega$ is the domain, $N_u$ is the displacement shape function and $\rho$ is the density matrix that holds density for all directions.

### Damping Matrix (D)

The mathematical definition is:
$$D = \alpha_R M + \beta_R K$$

Where $M$ and $K$ are the mass and stiffness  matrices respectively and $\alpha_R$ and $\beta_R$ are the coefficients from the Rayleigh Method.

File equation_of_motion_utilities.hpp includes 
-  CalculateMassMatrix function
-  CalculateDampingMatrix function
-  CalculateIntegrationCoefficientsInitialConfiguration function that calculates integration coefficient for all integration points

## Stress strain utilities

For convenience functions that compute invariants, equivalents and strain definitions.
Given a stress tensor $\sigma$ or a strain tensor $\epsilon$. The eigenvalues of the stress tensor are $\sigma_1 \le \sigma_2 \le \sigma_3$

$$\sigma = \begin{bmatrix} \sigma_{xx} & \sigma_{xy} & \sigma_{xz} \\
                           \sigma_{xy} & \sigma_{yy} & \sigma_{yz} \\
                           \sigma_{xz} & \sigma_{yz} & \sigma_{zz}  \end{bmatrix}$$

$$\epsilon = \begin{bmatrix} \epsilon_{xx} & \epsilon_{xy} & \epsilon_{xz} \\
                             \epsilon_{xy} & \epsilon_{yy} & \epsilon_{yz} \\
                             \epsilon_{xz} & \epsilon_{yz} & \epsilon_{zz}  \end{bmatrix}$$

### Trace

The first tensor invariant:

$$I_1 = trace(\sigma) = \Sigma_i \sigma_{i,i}$$

### Mean stress

$$p = \frac{1}{3} I_1 = \frac{1}{3} trace ( \sigma )$$

### Von Mises stress

With $J_2$ the second invariant of the tensor:

$$\overline\sigma = \sqrt{3 J_2} = \sqrt{0.5((\sigma_{xx}-\sigma_{yy})^2 +
                                 (\sigma_{yy}-\sigma_{zz})^2 +
                                 (\sigma_{zz}-\sigma_{xx})^2 ) +
                            3.0(\sigma_{xy}^2 + \sigma_{yz}^2 + \sigma_{xz}^2) }$$

### Von Mises strain

$$\overline{\epsilon} = \frac{2}{3} \sqrt{3 J_2}$$

### Green Lagrange strain tensor

With current configuration $x$ and reference configuration $X$, the deformation gradient $F = \frac{x}{X}$ and unit tensor $I$:

$$\epsilon = 0.5 ( F^T \cdot F - I )$$

### Hencky strain tensor

$$\epsilon = 0.5 \ln ( F^T \cdot F )$$

### Lode angle

The negative sine definition for Lode angle is adapted here [Lode coordinates Wikipedia](https://en.wikipedia.org/wiki/Lode_coordinates):

$$- \sin( 3 \bar{\theta}_s) = \frac{J_3}{2} (\frac{3}{J_2})^{\frac{3}{2}}$$

which brings:

$$\bar{\theta}_s = \frac{1}{3} \arcsin( - \frac{27}{2} \frac{(\sigma_1 - p)(\sigma_2 - p)(\sigma_3 - p)}{q^3})$$

### Mohr Coulomb shear capacity

Assessment of how the current stress utilizes the capacity as defined by the Mohr Coulomb yield surface.

$$\frac{q}{q_{mc}}$$

where

$$q_{mc} = \frac{3}{\sqrt{3}\cos \bar{\theta}_s - \sin \bar{\theta}_s \sin \phi }(p \sin \phi + c \cos \phi)$$

### Mohr Coulomb pressure capacity

Assessment of how the current stress utilizes the capacity as defined by the Mohr Coulomb yield surface.

$$(q_{mc} - q)\frac{3 \sin{\phi}}{\sqrt{3} \cos{\bar{\theta}_s} - \sin{\bar{\theta}_s} \sin{\phi}}$$
