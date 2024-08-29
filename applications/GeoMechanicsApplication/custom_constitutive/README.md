# Constitutive laws


## Incremental linear elastic interface law
The constitutive law for an interface element relates tractions $\tau$ to relative displacement $\Delta u$.
Relative displacement for interface element is the differential motion between the two sides of the interface. As a
consequence the relative displacement has unit of length [L] and the stiffness has unit of force over cubic length [F/L^3].

### Relative displacement and traction
In 2D plane strain y is the opening/closing direction of the interface, while differential motion in the tangential direction
gives shear. Like for a continuum, normal behaviour is placed first in the relative displacement and traction vector. Shear 
is placed after the normal motion or traction.

$$ \Delta u = \begin{bmatrix} \Delta u_y \\ \Delta u_x \end{bmatrix} $$

$$ \tau = \begin{bmatrix} \tau_{yy} \\ \tau_{xy} \end{bmatrix} $$

### 2D Elastic constitutive tensor

$$ C = \begin{bmatrix} C_{yy} & 0     \\
                       0     & C_{xy} \end{bmatrix}$$

Where $C_{yy}$ is input as `INTERFACE_NORMAL_STIFFNESS` and $C_{xy}$ is input as `INTERFACE_SHEAR_STIFFNESS`.

### Incremental relation

$$ \tau_{t + \Delta t} = \tau_t + C \cdot \Delta u $$

