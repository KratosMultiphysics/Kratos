# Constitutive laws


## Incremental linear elastic interface law

This constitutive law for an interface element linearly relates increments of tractions $\Delta \tau$ to increments of relative displacement $\Delta \Delta u$.
Relative displacement for interface element is the differential motion between the two sides of the interface. As a
consequence the relative displacement has unit of length $[\mathrm{L}]$ and the stiffness has unit of force over cubic length $[\mathrm{F/L^3}]$.

### Relative displacement and traction
In 2D plane strain $y$ is the opening/closing direction of the interface, while differential motion in the tangential direction
gives shear. Similar to the continuum, the normal behavior is placed first in the relative displacement and traction vectors. The shear behavior
follows after the normal motion or traction in these vectors.

$$ \Delta u = \begin{bmatrix} \Delta u_y \\ \Delta u_x \end{bmatrix} $$

$$ \tau = \begin{bmatrix} \tau_{yy} \\ \tau_{xy} \end{bmatrix} $$

Where:
* $\Delta u_y$: Relative displacement in the $y$-direction (normal to the interface).
* $\Delta u_x$: Relative displacement in the $x$-direction (tangential to the interface).
* $\tau_{yy}$: Traction in the $y$-direction (normal to the interface).
* $\tau_{xy}$: Shear traction in the $x$-direction (tangential to the interface).

### 2D Elastic constitutive tensor

The elastic behavior of the interface is characterized by the 2D elastic constitutive tensor $C$, which relates the traction and relative displacement vectors. The constitutive tensor in 2D is expressed as:

$$ C = \begin{bmatrix} C_{yy} & 0     \\
                       0     & C_{xy} \end{bmatrix}$$

Where:
* $C$: Represents the 2D elastic constitutive tensor.
* $C_{yy}$: Represents the `INTERFACE_NORMAL_STIFFNESS`, which characterizes the stiffness in the normal direction.
* $C_{xy}$: Represents the `INTERFACE_SHEAR_STIFFNESS`, which characterizes the stiffness in the tangential direction.

Both stiffness values have dimension $\mathrm{F/L^3}$.

### Incremental relation

The relation between the traction and the relative displacement for a given time increment is given by the following incremental equation:

$$ \tau_{t + \Delta t} = \tau_t + C \cdot Delta \Delta u $$

Where:
* $\tau_{t + \Delta t}$: The traction at the updated time $t + \Delta t$.
* $\tau_t$: The traction at the current time $t$.
* $C$: The elastic constitutive tensor.
* $\Delta \Delta u$: The incremental relative displacement vector.
