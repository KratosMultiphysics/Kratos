# Constitutive laws


## Incremental linear elastic interface law

This constitutive law for an interface element linearly relates increments of tractions $\Delta \tau$ to increments of relative displacement $\Delta \Delta u$.
Relative displacement for interface element is the differential motion between the two sides of the interface. As a
consequence the relative displacement has unit of length $[\mathrm{L}]$ and the stiffness has unit of force over cubic length $[\mathrm{F/L^3}]$.
Currently, this law is implemented for 2D and 3D cases. The constitutive matrix is
$$K = \begin{bmatrix} k_n \quad 0 \quad 0 \\ 0 \quad k_t \quad 0 \\ 0 \quad 0 \quad k_s \end{bmatrix}$$  
where $k_n$ is normal, $k_t$ and $k_s$ are tangential components. The current implementation uses $k_s=k_t$.  

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

$$ \tau_{t + \Delta t} = \tau_t + C \cdot \Delta \Delta u $$

Where:
* $\tau_{t + \Delta t}$: The traction at the updated time $t + \Delta t$.
* $\tau_t$: The traction at the current time $t$.
* $C$: The elastic constitutive tensor.
* $\Delta \Delta u$: The incremental relative displacement vector.


## Mohr-Coulomb with tensile cutoff

The Mohr-Coulomb failure criterion is defined as:

```math
    F_{MC}(\sigma) = \frac{\sigma_1 - \sigma_3}{2} + \frac{\sigma_1 + \sigma_3}{2} \sin⁡{\phi} - c \cos⁡{\phi} = 0
```

where:

- $`\sigma_1`$ = maximum principal stress component
- $`\sigma_3`$ = minimum principal stress component
- $`c`$ = cohesion
- $`\phi`$ = Internal friction angle.

This criterion represents a linear envelope in the Mohr stress space, approximating the shear strength of a material under different stress states.

Since the Mohr-Coulomb criterion primarily accounts for shear failure, it does not limit tensile stresses. In geomechanical applications, materials such as rocks and soils have a limited tensile strength . A tensile cutoff is imposed as:

```math
    F_{tc}(\sigma) = \sigma_1 - t_c = 0
```

where $t_c$ is the tensile strength. If $`\sigma_1`$ exceeds $`t_c`$, failure occurs regardless of the shear strength condition.

Combination of these two, it yiels to the following figure:

<img src="documentation_data/mohr-coulomb-with-tension-cutoff-zones.svg" alt="Mohr-Coulomb with tension cutoff" title="Mohr-Coulomb with tension cutoff" width="800">


### Implementation

To incorporate the Mohr-Coulomb model with tensile cutoff in numerical simulations, the following steps are followed:

1. Calculate the trial stress by: 

```math
    \boldsymbol{\sigma}^{trial} = \boldsymbol{\sigma}^0 + \boldsymbol{\mathrm{C}} \Delta \boldsymbol{\epsilon}
```

2. Extract the principal stresses $`\sigma_1 \ge \sigma_2 \ge \sigma_3`$ for the calculated trial stress, and calculate the rotation matrix.

3. Calculate the values of $`F_{MC}`$ and $`F_{tC}`$ for the calculated principal stresses.

4. Evaluate the condition and mapping
  - If the trial stress falls in the elastic zone, it stays unchanged. No mapping is applied.
  - If the trial stress falls in the tensile apex return zone. The trial stress then needs to be mapped back to the apex.
  - If the trial stress falls in the tensile cutoff zone. The trial stress then needs to be mapped back to the tension cutoff line.
  - If it falls in the tensile corner return zone, then it needs to be mapped to the corner point.
  - In the case of regular failure zone, then it is mapped back to the Mohr-Coulomb curve along the normal direction of flow function. The flow function is defined by
  
```math
    G(\sigma) = \frac{\sigma_1 - \sigma_3}{2} + \frac{\sigma_1 + \sigma_3}{2} \sin⁡{\psi}
```
  where $`\psi`$ is the dilatancy angle.

5. If after mapping, the condidition $`\sigma_1 \ge \sigma_2 \ge \sigma_3`$ is not valid, average the principal stresses of stage 2 and the direction of the mapping and map the principal stresses again.
  - if $`\sigma_1 \le \sigma_2`$ set:
```math
       \sigma_1 = \sigma_2 = \frac{\sigma_1 + \sigma_2}{2}
```
```math
       \frac{\partial G}{\partial \sigma_1} = \frac{\partial G}{\partial \sigma_2} = \frac{1}{2} \left( \frac{\partial G}{\partial \sigma_1} + \frac{\partial G}{\partial \sigma_2} \right)
```
  - if $`\sigma_2 \le \sigma_3`$ set:
```math
       \sigma_3 = \sigma_2 = \frac{\sigma_3 + \sigma_2}{2}
```
```math
       \frac{\partial G}{\partial \sigma_3} = \frac{\partial G}{\partial \sigma_2} = \frac{1}{2} \left( \frac{\partial G}{\partial \sigma_3} + \frac{\partial G}{\partial \sigma_2} \right)
```
This mapping is based on a new Mohr-Coulomb diagram with modified zones, based on the averaged of the derivatives of flow functions $`\frac{\partial G}{\partial \boldsymbol{\sigma}}`$.

6. Rotate the mapped stress vector back, by appying the rotation matrix.

### Detailed formulations
We define the Coulomb yield surface $F_{MC}$ and the tensile cutoff surface $F_{tc}$ based on $\sigma-\tau$ coordinates. as:

```math
    F_{MC} = \tau + \sigma \sin{\phi} - C \cos{\phi} = 0
```
```math
    F_{tc} = \sigma +\tau - t_c = 0
```

Where
```math
    \sigma = \frac{\sigma_1 + \sigma_3}{2}
```
```math
    \tau = \frac{\sigma_1 - \sigma_3}{2}
```

#### Elastic region: 
This condition occurs when $F_{MC} \le 0$ and $F_{tc} \le 0$. Here, the stress vector stays unchanged and no need any return mapping.


#### Tensile apex return zone

First need to find whether there is intersection between the Yield and cutoff functions at the region of $\tau \ge 0$.

Find the root of $F_{MC}$ (apex).
```math
    \sigma_{MC} = \frac{C}{\tan{\phi}}
```

If $t_c < \sigma_{MC}$ we find a line perpendicular to the tensile-cutoff curve passing through the apex point. The equation for the tension-cutoff is:
```math
    \tau = - \sigma + t_c
```

Then the perpendicular line passing from the apex point is:
```math
    \tau - \sigma + t_c = 0
```
Any trial stress which falls below this line is then belong to this region. It is namely:
```math
    \tau^{trial} - \sigma^{trial} + t_c < 0
```
If a point falls in this zone, it will be mapped back to the root point of the tension-cutoff line, namely to point $\sigma = t_c$ and $\tau = 0$. Then update the principal stresses based on the mapper values.
```math
    \sigma_1 = \sigma + \tau
```
```math
    \sigma_3 = \sigma - \tau
```

They are the corrected principal stresses. They need to be rotated back to the element axes system. We need to use the rotation matrix.
```math
   \boldsymbol{\sigma} = \boldsymbol{R \sigma_p R^{-1}}
```

The rotation matrix is orthogonal
```math
    \boldsymbol{R^T} = \boldsymbol{R^{-1}} \Rightarrow \boldsymbol{\sigma} = \boldsymbol{R \sigma_p R^{-1}}
```

Where
```math
    \boldsymbol{\sigma_p} =
    \begin{bmatrix}
    \sigma_1 & 0 & 0 \\
    0 & \sigma_2 & 0 \\
    0 & 0 & \sigma_3
    \end{bmatrix}
```

#### Tensile cutoff return zone

We need to find the shear at the intersection point. Setting  $\sigma_1 = t_c$ in the yield function (if $t_c < \sigma_{MC}$):
```math
    -\sigma + t + \sigma \sin{\phi} - C \cos{\phi} = 0
```

Then 
```math
    \tau_{corner} = \frac{C \cos⁡{\phi} - t_c \sin{\phi}⁡}{1-\sin{\phi}}
```
```math
    \sigma_{corner} = \frac{t_c - C \cos⁡{\phi}}{1 - \sin⁡{ϕ}}
```

A perpendicular to the tension-cutoff curve, passing from the corner point is:
```math
    \tau - \tau_{corner} = \sigma - \sigma_{corner}
```

Then the condition is:
```math
    \left(\tau^{trial} - \tau_{corner}) - (\sigma^{trial} - \sigma_{corner} \right) < 0
```

Each point which is outside the previous region (tensile-apex region) and follows this condition, then we need to map the trial stresses to the tension-cutoff curve. It means:
```math
    \sigma + \tau = t_c
```

We move perpendicular to the tension-cutoff surface. We use here the derivative of the flow function related to the tension cutoff surface.
```math
    \dot{\lambda} = \frac{\sigma + \tau - t_c}{\partial G_t / \partial \boldsymbol{\sigma}}
```
```math
    \frac{\partial G_t}{\partial \boldsymbol{\sigma}} = 
    \begin{bmatrix} 1 \\
    1
    \end{bmatrix}
```
```math
    \sigma^{map} = \sigma^{trial} + \dot{\lambda} \frac{\partial G_t}{\partial \sigma}
```
They are the corrected principal stresses. They need to be rotated back to the element axes system.


#### Zone of tensile corner return

This zone is located under the line which is perpendicular to the flow function and passes through the intersection point of yield function and tension cutoff. 

First we find the intersection point. The intersection point is related whether the tension-cutoff crosses the yield function. 

- The tension-cut-off crosses the yield function:
In this case of crossing the shear can be found like the previous region:
```math
    \tau_{corner} = \frac{C \cos{\phi} - t_c \sin{\phi}}{1 - \sin{\phi}}
```

And the normal stress can be found by
```math
    \sigma_{corner} = \frac{t_c - C \cos{\phi}}{1 - \sin{\phi}}
```

- The tension-cutoff is located outside:
Then this point is equal to the root, which was derived in the previous section:
```math
    \tau_{corner} = 0
```
```math
    \sigma_{corner} = \frac{C}{\tan{\phi}}
```

- Flow function:
The flow function is:
```math
    G_{MC}(\boldsymbol{\sigma}) = \tau + \sigma \sin{\psi} = 0
```
```math
    \tau =- \sigma \sin{\psi}
```

The slope of this line is $-\sin⁡{\psi}$. The slope of a line perpendicular is then
```math
    m = \frac{1}{\sin{\psi}}
```

If we consider a line as $\tau = m \sigma + B$, and substituting the intersection point
```math
    \tau_{corner} =m \sigma_{corner} + B
```
```math
    B = \tau_{corner} -m \sigma_{corner}
```

Then
```math
    g = (\tau - \tau_{corner}) - \frac{1}{\sin{\psi}} (\sigma - \sigma_{corner}) = 0
```

This zone is defined in the region below this function and above the axial region (above the shear at intersection).
```math
    g \ge 0 , \tau_{trial} > \tau_{corner}
```

Then
```math
    \tau = \tau_{corner}
```
```math
    \sigma = \sigma_{corner}
```

They are the corrected principal stresses. They need to be rotated back to the element axes system. We need to use the rotation matrix, similar as done above for other regions. 


#### Zone of regular failure
This zone is associated with the region above the yield function and above the function g derived in the previous section.

$$ F_{MC} > 0, g < 0 $$

We can use the derivative of the flow function to define the direction and find the return point on the yield surface.

$$ \frac{\partial G_{MC}}{\partial \boldsymbol{\sigma}} = 
\begin{bmatrix}
\sin{\psi} \\
1
\end{bmatrix} =
\begin{bmatrix}
n_1 \\
n_2
\end{bmatrix} $$

Then a parametrized line can be defined by:

$$ \sigma = \sigma^{trial} + \dot{\lambda} n_1 $$
$$ \tau = \tau^{trial} + \dot{\lambda} n_2 $$

At yield function,

$$ F_{MC} = \tau + \sigma \sin{\phi} - C \cos⁡{\phi} = 0 $$
$$ \sigma \sin{\phi} + \tau = C \cos{\phi} $$

Solving this 3 equations:

$$ \dot{\lambda} = \frac{C_2 - \sigma^{trial} C_1 - \tau^{trial}}{n_1 C_1 + n_3} $$

Where $C_1 = \sin{\phi}$ and $C_2 = C \cos{\phi} $. Then

$$ \sigma_1 = \dot{\lambda} \frac{\partial G_{MC}}{\partial \boldsymbol{\sigma}} + \sigma_1^{trial} $$

They are the corrected principal stresses. They need to be rotated back to the element axes system.

### Rotation matrix
The rotation matrix is derived from the eigenvectors of the Cauchy stresses. Having three eigenvectors related to the principal stresses

$$ \begin{bmatrix} v_1 & v_2 & v_3 \end{bmatrix} $$

Normalizing the vectors, it results in rotation matrix

$$ \boldsymbol{R} = \begin{bmatrix} \frac{v_1}{\lVert v_1 \rVert} & \frac{v_2}{\lVert v_2 \rVert} & \frac{v_3}{\lVert v_3 \rVert} \end{bmatrix} $$


### Reordering
It can happen that, after mapping, the role of the principal stress change, and the condition $\sigma_1 \ge \sigma_2 \ge \sigma_3$ is no longer valid. In such a case, we apply averaging to the principal stresses and their associated mapping directions.

- Case $\sigma_1 < \sigma_2$: Then we use averaging on the initial principle trial stresses (principal trial stresses before mapping) as follows:
	
$$ \sigma_1 = \sigma_2 = \frac{\sigma_1 + \sigma_2}{2} $$
$$ \frac{\partial G}{\partial \sigma_1} = \frac{\partial G}{\partial \sigma_2} = \frac{1}{2} \left(\frac{\partial G}{\partial \sigma_1} + \frac{\partial G}{\partial \sigma_2} \right) $$

Where $G$ is the flow function. For Mohr-Coulomb, the derivative of flow function is:

$$ \frac{\partial G_{MC}}{\partial \boldsymbol{\sigma}} =
\begin{bmatrix}
\frac{1}{2} \left( 1 + \sin{\psi} \right) \\
0 \\
\frac{1}{2} \left( -1 + \sin{\psi} \right)
\end{bmatrix}$$

Then the averaging leads to:

$$ \frac{\partial G}{\partial \sigma_1} = \frac{\partial G}{\partial \sigma_2} = \frac{1}{4} \left( 1 + \sin⁡{\psi} \right) $$
$$ \frac{\partial G}{\partial \sigma_3} = - \frac{1}{2} \left( 1 - \sin⁡{\psi} \right) $$

As we solve our mapping based on $\sigma$ and $\tau$, we need to convert this to:

$$ \frac{\partial G}{\partial \sigma} = - \frac{1}{4} \left( 1 - 3 \sin⁡{\psi} \right) $$
$$ \frac{\partial G}{\partial \tau} = \frac{1}{4} \left( 3 - \sin{\psi} \right) $$

- Case $\sigma_1 < \sigma_2$: Then we use averaging on the initial principle trial stresses (principal trial stresses before mapping) as follows:
	
$$ \sigma_1 = \sigma_2 = \frac{\sigma_1 + \sigma_2}{2} $$
$$ \frac{\partial G}{\partial \sigma_1} = \frac{\partial G}{\partial \sigma_2} = \frac{1}{2} \left( \frac{\partial G}{\partial \sigma_1} + \frac{\partial G}{\partial \sigma_2} \right) $$

Then the averaging of the mapping direction leads to:

$$ \frac{\partial G}{\partial \sigma_3} = \frac{\partial G}{\partial \sigma_2} = - \frac{1}{4} \left( 1 - \sin⁡{\psi} \right) $$
$$ \frac{\partial G}{\partial \sigma_1} = \frac{1}{2} \left( 1 + \sin⁡{\psi} \right) $$

And based on $\sigma$ and $\tau$:

$$ \frac{\partial G}{\partial \sigma} = \frac{1}{4} \left( 1 + 3 \sin⁡{\psi} \right) $$
$$ \frac{\partial G}{\partial \tau} = \frac{1}{4} \left( 3 + \sin⁡{\psi} \right) $$

Note that after averaging the mapping direction, we modify the Mohr-Coulomb curve to account for the modified mapping direction. 
The mapping direction for tension cutoff stays unchanged because applying such averaging leads to the same form of mapping. After averaging the mapping for tension cutoff stays unchanged.


### Hardening and softening
In the hardening/softening process, the yield parameters are a not constant anymore, but they will be a function of shear plastic strain $\kappa$. The increment of equivalent shear plastic strain is defined by:

$$ \Delta \kappa = \sqrt{2/3} \lVert \Delta \epsilon^p \rVert$$
$$ \kappa_{n+1} = \kappa_n + \Delta \kappa $$

Where:

$$ \Delta \epsilon^p = \dot{\lambda} \frac{\partial G}{\partial \boldsymbol{\sigma}} $$

As the current implementation is based on $\sigma-\tau$, we must map that 2-vector back to a 3×3 (or Voigt) flow tensor. We use chain rule to get derivatives of $G$ to the principal stresses:

$$ \frac{\partial G}{\partial \sigma_1} = 
\frac{\partial G}{\partial \sigma} \frac{\partial \sigma}{\partial \sigma_1}
+ \frac{\partial G}{\partial \tau} \frac{\partial \tau}{\partial \sigma_1} 
= \frac{1}{2} \frac{\partial G}{\partial \sigma} + \frac{1}{2} \frac{\partial G}{\partial \tau} $$

$$ \frac{\partial G}{\partial \sigma_3} = \frac{1}{2} \frac{\partial G}{\partial \sigma} - \frac{1}{2} \frac{\partial G}{\partial \tau}$$

$$ \frac{\partial G}{\partial \sigma_2} = 0 $$

Denote these principal derivatives $g_1$, $g_2$ and $g_3$:

$$ g_1 = \frac{\partial G}{\partial \sigma_1},\qquad g_2 = 0,\qquad  g_3 = \frac{\partial G}{\partial \sigma_3} $$

We build the flow tensor $\boldsymbol{m}$ in principal frame.

$$ \boldsymbol{m} = 
\begin{bmatrix}
g_1 & 0 & 0 \\
0 & g_2 & 0 \\
0 & 0 & g_3
\end{bmatrix} $$

The deviatoric is then calculated:

$$ \boldsymbol{m_{dev}} = \boldsymbol{m} - \frac{1}{3} \left( \mathrm{tr} \, \boldsymbol{m} \right) $$

In principal components (diagonal), the mean value is

$$\bar{m} = \frac{1}{3} \left( g_1 + g_2 + g_3 \right) $$

Then

$$ \lVert \boldsymbol{m_dev} \rVert = \sqrt{(g_1 - \bar{m})^2 + (g_2 - \bar{m})^2 + (g_3 - \bar{m})^2} $$

If the accumulated hardening variable is the usual equivalent plastic strain

$$ \kappa = \epsilon_p, \qquad \dot{\epsilon_p} = \sqrt{\frac{2}{3} \dot{\epsilon^p}:\dot{\epsilon^p}} $$

As $\Delta \epsilon^p = \dot{\lambda} \boldsymbol{m}$

$$ \Delta \kappa = \dot{\lambda} \sqrt{\frac{2}{3} \boldsymbol{m:m}} $$


The Cohesion, friction angle, dilation angle and tangential cutoff become all functions of $\kappa$.

We start from the most simple formulations for hardening and softening. It is namely linear hardening.

#### Linear hardening:
As it is mentioned above, in hardening process, the material properties for Coulomb yield surface are functions of $\kappa$. Here we use the most simple hardening model, namely linear. 

$$ \phi(\kappa) = \phi_0 + H_\phi \kappa $$
$$ C(\kappa) = C_0 + H_C \kappa $$
$$ \psi(\kappa) = \psi_0 + H_\psi \kappa $$

Where $H_\phi$, $H_C$ and $H_\psi$ are hardening modulus for the friction angle, cohession and dilatation angle, respectively.

Note: These formulations will be extended for more physics-based form.


#### Iterative process

1. Compute the yield surface and map the trial stresses:

 $$ F_{MC}(\boldsymbol{\sigma}, \kappa_n) = \tau + \sigma \sin{\phi(\kappa_n)} - C(\kappa_n) \cos{\phi(\kappa_n)} = 0 $$

2. Compute plastic multiplier increment $\dot{\lambda}$

3. Update plastic strain:

$$ \Delta \boldsymbol{\epsilon^p} = \dot{λ} \frac{\partial G}{\partial \boldsymbol{\sigma}} $$

4. Update the hardening variable:

$$ \kappa_{n+1} = \kappa_n + \sqrt{2/3} \lVert \Delta \boldsymbol{\epsilon^p} \rVert $$

5. Update the material parameters:
C_{n+1} = C(\kappa_{n+1})
\phi_{n+1} = \phi(\kappa_{n+1})
\psi_{n+1} = \psi(\kappa_{n+1})

6. Recompute yield surface with updated parameters.
 $$ F_{MC}(\boldsymbol{\sigma}, \kappa_{n+1}) = \tau + \sigma \sin{\phi(\kappa_{n+1})} - C(\kappa_{n+1}) \cos{\phi(\kappa_{n+1})} = 0 $$

7.	Go to 1 until convergence.

The convergence criterion is defined as:

$$ \lvert F_{MC}(\boldsymbol{\sigma}, \kappa_{n+1}) \rvert < \epsilon $$

Where $\epsilon$ is a tolerance.
