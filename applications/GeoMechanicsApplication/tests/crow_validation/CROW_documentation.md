# CROW Validation documentation

This document describes the validation of the CROW case. The case is based on the following work:

Directory: CROW\CROW-case\CC3-nieuwbouw-50jaar_beta=4,3_met_modelonzekerheid_CoV-Qlast=0,13\Case 2 - CUR166 case
![D-Sheet_Piling.png](D-Sheet_Piling.png)

## Soil parameters

This section documents the derivation of the Kratos input values for the **clay** and **sand** layers.

### Given data

From D-Sheet Piling, the following soil properties are given for the clay and sand layers:

| Property                               | Clay                 | Sand                 | Unit                       |
|:---------------------------------------|:---------------------|:---------------------|:---------------------------|
| Unsaturated total unit weight          | 18.0                 | 20.0                 | $\mathrm{kN}/\mathrm{m}^3$ |
| Saturated total unit weight            | 18.0                 | 20.0                 | $\mathrm{kN}/\mathrm{m}^3$ |
| Cohesion                               | 3.0                  | 0.0                  | $\mathrm{kN}/\mathrm{m}^2$ |
| Friction angle                         | 22.5                 | 32.5                 | $`^{\circ}`$               |
| Delta friction angle                   | 11.25                | 20.0                 | $`^{\circ}`$               |
| OCR                                    | 1.0                  | 1.0                  | $-$                        |
| Earth pressure coefficient             | 0.62                 | 0.46                 | $-$                        |
| Modulus of subgrade reaction-Secant k1 | $1.0 \times 10^{3}$  | $1.0 \times 10^{4}$  | $\mathrm{kN}/\mathrm{m}^3$ |
| Modulus of subgrade reaction-Secant k2 | $1.0 \times 10^{3}$  | $1.0 \times 10^{4}$  | $\mathrm{kN}/\mathrm{m}^3$ |
| Modulus of subgrade reaction-Secant k3 | $1.0 \times 10^{3}$  | $1.0 \times 10^{4}$  | $\mathrm{kN}/\mathrm{m}^3$ |
| Horizontal permeability                | $1.0 \times 10^{-4}$ | $1.0 \times 10^{-4}$ | $\mathrm{m}/\mathrm{s}$    |


### Elastic properties

Since D-Sheet Piling does not work with Young's modulus $`E`$ directly, we will calculate it based on the modulus of subgrade reaction $`k`$, the width $`b`$ of the soil in $z$ direction, and Poisson's ratio $`\nu`$.  The latter is assumed to be equal to $`0.2`$ (for clay) and $`0.3`$ (for sand).  $`b = 1.0\ \mathrm{m}`$, since we assume plane strain conditions.  The Young's modulus is now calculated using the following formula (based on the elastic tensor):

```math
E = k \cdot b \cdot \frac{(1 + \nu) (1 - 2 \nu)}{1 - \nu}
```

For clay, the Young's modulus then becomes $`E_{\mathrm{clay}} = 1.0 \times 10^6\ \mathrm{N} / \mathrm{m}^3 \cdot 1.0\ \mathrm{m} \cdot \frac{(1 + 0.2) (1 - 2 \cdot 0.2)}{1 - 0.2} = 9.0 \times 10^5\ \mathrm{N} / \mathrm{m}^2`$.  And for sand, the Young's modulus becomes $`E_{\mathrm{sand}} = 1.0 \times 10^7\ \mathrm{N} / \mathrm{m}^3 \cdot 1.0\ \mathrm{m} \cdot \frac{(1 + 0.3) (1 - 2 \cdot 0.3)}{1 - 0.3} = 7.4286 \times 10^6\ \mathrm{N} / \mathrm{m}^2`$. 


### Conversion to intrinsic permeability

From D-Sheet Piling, the permeability is given as **hydraulic conductivity** for both soils:

```math
K = 1.0 \times 10^{-4}\ \mathrm{m}/\mathrm{s}
```

Kratos requires **intrinsic permeability**:

```math
k\ [\mathrm{m}^2]
```

The relation between hydraulic conductivity and intrinsic permeability is:

$$ k = \frac{K \mu}{\rho_{\mathrm{w}} g} $$

where:

* $K$ is the hydraulic conductivity
* $\mu$ is the dynamic viscosity
* $\rho_{\mathrm{w}}$ is the water density
* $g$ is the gravitational acceleration

The adopted values are:

```math
\mu = 1.0 \times 10^{-3}\ \mathrm{Pa}\cdot\mathrm{s}
```

```math
\rho_{\mathrm{w}} = 1019.37\ \mathrm{kg}/\mathrm{m}^3
```

```math
g = 9.81\ \mathrm{m}/\mathrm{s}^2
```

Substituting:

$$ k = \frac{1.0 \times 10^{-4} \cdot 1.0 \times 10^{-3}}{1019.37 \cdot 9.81} $$

```math
k = 9.999 \times 10^{-12}\ \mathrm{m}^2
```

This is rounded to:

```math
k \approx 1.0 \times 10^{-11}\ \mathrm{m}^2
```

**Note:** Since the water pressure field will be completely prescribed, the above permeability values as well as the dynamic viscosity will not have any influence on the simulation results.


### Final values
Because the Kratos model requires additional parameters, the table below summarizes the final material properties for the soil layers, including both the calculated values and the assumed values adopted for use in the model.

| Property                     | Kratos input parameter | Clay                     | Sand                     | Unit                         |
|:-----------------------------|:-----------------------|:-------------------------|:-------------------------|:-----------------------------|
| Solid density                | `DENSITY_SOLID`        | 1834.86                  | 2038.74                  | $`\mathrm{kg}/\mathrm{m}^3`$ |
| Water density                | `DENSITY_WATER`        | 1019.37                  | 1019.37                  | $`\mathrm{kg}/\mathrm{m}^3`$ |
| Porosity                     | `POROSITY`             | 0.0                      | 0.0                      | $`-`$                        |
| Young's modulus              | `YOUNG_MODULUS`        | $`9.0 \times 10^5`$      | $`7.4286 \times 10^6`$   | $`\mathrm{Pa}`$              |
| Poisson's ratio              | `POISSON_RATIO`        | 0.2                      | 0.3                      | $`-`$                        |
| Saturated saturation         | `SATURATED_SATURATION` | 1.0                      | 1.0                      | $`-`$                        |
| Residual saturation          | `RESIDUAL_SATURATION`  | $`1.0 \times 10^{-10}`$  | $`1.0 \times 10^{-10}`$  | $`-`$                        |
| Earth pressure coefficient   | `K0_NC`                | 0.62                     | 0.46                     | $`-`$                        |
| Cohesion                     | `GEO_COHESION`         | 3000.0                   | 0.0                      | $`\mathrm{Pa}`$              |
| Friction angle               | `GEO_FRICTION_ANGLE`   | 22.5                     | 32.5                     | $`^{\circ}`$                 |
| Dilatancy angle              | `GEO_DILATANCY_ANGLE`  | 0.0                      | 0.0                      | $`^{\circ}`$                 |
| Tensile strength             | `GEO_TENSILE_STRENGTH` | 7242.64                  | 0.0                      | $`\mathrm{Pa}`$              |

In the above table, the value of the tensile strength $`f_{\mathrm{t}}`$ (which fixes the tension cut-off) has been chosen such that it passes through the apex of the Coulomb yield surface.  In other words, it has been calculated as follows: $`f_{\mathrm{t}} = \frac{c}{\mathrm{tan}(\phi)}`$, with $`c`$ the cohesion and $`\phi`$ the friction angle.  Note that this choice effectively disables the tension cut-off. 


## Interface parameters

This section documents the derivation of the Kratos input values for the interfaces adjacent to **clay** and **sand**.


### Given data

The interface stiffness values are based on the adjacent soil shear modulus and a characteristic element size normal to the interface. The friction angles applied to the interfaces have been taken from the D-Sheet Piling analysis, where they are known as the "Delta friction angle". The cohesion values of the interfaces are related to the cohesion values of the adjacent clay and sand layers using a reduction factor that is derived from the tangent of the friction angles, as follows:

```math
c_{\mathrm{interface}} = \frac{\mathrm{tan}(\phi_{\mathrm{interface}})}{\mathrm{tan}(\phi_{\mathrm{soil}})}
```

By substituting the applicable values for the interface that is connected to the clay, the reduced cohesion is calculated as

```math
c_{\mathrm{interface, clay}} = \frac{\mathrm{tan}(11.25^{\circ}})}{\mathrm{tan}(22.5^{\circ})} = 1440.65\ \mathrm{Pa}
```

For the sand-side interface, there is no need to calculate the applicable reduction factor, since the cohesion of the sand equals 0.0 Pa. Consequently, the cohesion applied to the sand-side interface also equals 0.0 Pa.

The following table lists the adopted values.

| Property                         | Clay-side interface | Sand-side interface    | Unit            |
|:---------------------------------|:--------------------|:-----------------------|:----------------|
| Young's modulus of adjacent soil | $`9.0 \times 10^5`$ | $`7.4286 \times 10^6`$ | $`\mathrm{Pa}`$ |
| Poisson's ratio of adjacent soil | 0.2                 | 0.3                    | $`-`$           |
| Friction angle of adjacent soil  | 22.5                | 32.5                   | $`^{\circ}`$    |
| Friction angle of interface      | 11.25               | 20.0                   | $`^{\circ}`$    |
| Cohesion of adjacent soil        | 3000.0              | 0.0                    | $`\mathrm{Pa}`$ |
| Cohesion of interface            | 1440.65             | 0.0                    | $`\mathrm{Pa}`$ |
| Element size normal to interface | 1.0                 | 1.0                    | $`\mathrm{m}`$  |


### Shear modulus of the adjacent soil

The soil shear modulus is calculated as:

```math
G = \frac{E}{2 (1 + \nu)}
```


#### Clay

```math
G_{\mathrm{clay}} = \frac{9.0 \times 10^5}{2 (1 + 0.2)} = 3.75 \times 10^5\ \mathrm{Pa}
```


#### Sand

```math
G_{\mathrm{sand}} = \frac{7.4286 \times 10^6}{2 (1 + 0.3)} = 2.8571 \times 10^6\ \mathrm{Pa}
```


### Interface stiffness formulation

The interface shear stiffness is taken as:

```math
k_{\mathrm{s}} = \frac{G}{h}
```

The interface normal stiffness is taken as:

```math
k_{\mathrm{n}} = 10\ k_{\mathrm{s}}
```

with:

* $`G`$ = shear modulus of the adjacent soil
* $`h`$ = element size normal to the interface


### Clay-side interface

Using:

```math
G_{\mathrm{clay}} = 3.75 \times 10^5\ \mathrm{Pa}
```

and

```math
h = 1.0\ \mathrm{m}
```

the interface shear stiffness becomes:

```math
k_{\mathrm{s, clay}} = \frac{3.75 \times 10^5\ \mathrm{Pa}}{1.0\ \mathrm{m}} = 3.75 \times 10^5\ \mathrm{N} / \mathrm{m}^3
```

and the interface normal stiffness becomes:

```math
k_{\mathrm{n, clay}} = 10 \cdot 3.75 \times 10^5\ \mathrm{N} / \mathrm{m}^3 = 3.75 \times 10^6\ \mathrm{N} / \mathrm{m}^3
```


### Sand-side interface

Using:

```math
G_{\mathrm{sand}} = 2.8571 \times 10^6\ \mathrm{Pa}
```

and

```math
h = 1.0\ \mathrm{m}
```

the interface shear stiffness becomes:

```math
k_{\mathrm{s, sand}} = \frac{2.8571 \times 10^6\ \mathrm{Pa}}{1.0\ \mathrm{m}} = 2.8571 \times 10^6\ \mathrm{N} / \mathrm{m}^3
```

and the interface normal stiffness becomes:

```math
k_{\mathrm{n, sand}} = 10 \cdot 2.8571 \times 10^6\ \mathrm{N} / \mathrm{m}^3 = 2.8571 \times 10^7\ \mathrm{N} / \mathrm{m}^3
```


### Final values

The following table lists the adopted properties for the two types of interfaces. The values of the friction angle have been taken from the D-Sheet Piling input, where they are listed under "Delta friction angle".  Just like for the adjacent soil, the dilatancy angle has been assumed to be equal to $`0.0^{\circ}`$

| Property                            | Kratos input parameter       | Clay-side interface     | Sand-side interface    | Unit                          |
|:------------------------------------|:-----------------------------|:------------------------|:-----------------------|:------------------------------|
| Normal stiffness $`k_{\mathrm{n}}`$ | `INTERFACE_NORMAL_STIFFNESS` | $`3.75 \times 10^6`$    | $`2.8571 \times 10^7`$ | $`\mathrm{N} / \mathrm{m}^3`$ |
| Shear stiffness $`k_{\mathrm{s}}`$  | `INTERFACE_SHEAR_STIFFNESS`  | $`3.75 \times 10^5`$    | $`2.8571 \times 10^6`$ | $`\mathrm{N} / \mathrm{m}^3`$ |
| Cohesion $`c`$                      | `GEO_COHESION`               | $`1.44065 \times 10^3`$ | 0.0                    | $`\mathrm{Pa}`$               | 
| Friction angle $`\phi`$             | `GEO_FRICTION_ANGLE`         | 11.25                   | 20.0                   | $`^{\circ}`$                  |
| Dilatancy angle $`\psi`$            | `GEO_DILATANCY_ANGLE`        | 0.0                     | 0.0                    | $`^{\circ}`$                  |
| Tensile strength $`f_{\mathrm{t}}`$ | `GEO_TENSILE_STRENGTH`       | $`7.24263 \times 10^3`$ | 0.0                    | $`\mathrm{Pa}`$               |


## Sheet pile parameters

This section documents the parameters of the sheet pile, which is represented as a Timoshenko beam in the Kratos model.


### Given data

From the section data for **AZ26** in D-Sheet Piling, the following properties are given for the sheet pile (except for the weight $`G`$):

| Property                          | Value                | Unit                              |
|:----------------------------------|:---------------------|:----------------------------------|
| Bending stiffness $`EI`$          | $`8.40 \times 10^4`$ | $`\mathrm{kNm}^2 / \mathrm{m}^1`$ |
| Section area per meter wall $`A`$ | 198                  | $`\mathrm{cm}^2 / \mathrm{m}^1`$  |
| Weight $`G`$                      | 146.9                | $`\mathrm{kg} / \mathrm{m}^1`$    |

The weight $`G`$ of the sheet pile wall has been taken from the following [manufacturer's information sheet](https://sheetpiling.arcelormittal.com/products/az-sections/az-700-and-az-770/az-26-700).

The Young’s modulus of steel sheet piles is generally considered to be $`210\ \mathrm{GPa}`$.  This is the standard modulus of elasticity for structural steel.


### Adopted Kratos beam representation

In the current Kratos model, the sheet pile is modeled using a Timoshenko beam with a rectangular cross-section.  To ensure an equivalent bending stiffness $`(EI)_{\mathrm{beam}}`$ and an equivalent extensional stiffness $`(EA)_{\mathrm{beam}}`$, the Young's modulus and the thickness of the cross-section have been calculated such that these stiffness values match the ones taken from D-Sheet Piling.  The equivalent bending stiffness is calculated as follows:

```math
(EI)_{\mathrm{beam}} = E_{\mathrm{beam}} \cdot \frac{1}{12} \cdot b \cdot t^3 = EI 
```

The extensional stiffness is calculated as follows:

```math
(EA)_{\mathrm{beam}} = E_{\mathrm{beam}} \cdot b \cdot t = EA
```

In both equations, the width $`b`$ equals $`1.0\ \mathrm{m}`$, since we assume plane strain conditions.  From these two equations, we can solve for the equivalent Young's modulus of the beam $`E_{\mathrm{beam}}`$ and the equivalent wall thickness $`t`$.  We arrive at the following closed-form solution for the equivalent thickness $`t`$:

```math
t = \sqrt{12 \cdot \frac{EI}{EA}}
```

Through back-substitution, we can calculate the equivalent Young's modulus of the beam $`E_{\mathrm{beam}}`$:

```math
E_{\mathrm{beam}} = \frac{EA}{b \cdot t}
```

Finally, since the adopted cross-section has a rectangular shape, also the density of steel needs to be modified such that the weight matches the original value.

```math
\rho_{\mathrm{s}} = \frac{W}{b \cdot t} = \frac{146.9\ \mathrm{kg} / \mathrm{m}}{1.0\ \mathrm{m} \cdot 0.4924\ \mathrm{m}} = 298.36\ \mathrm{kg} / \mathrm{m}^3
```

Note that $`\rho_{\mathrm{s}}`$ represents an _apparent_ density of the steel rather than the actual one.


### Final values

| Property                           | Kratos input parameter  | Value                    | Unit                           |
|:-----------------------------------|:------------------------|:-------------------------|:-------------------------------|
| Young's modulus                    | `YOUNG_MODULUS`         | $`8.4449 \times 10^9`$   | $`\mathrm{Pa}`$                |
| Poisson's ratio                    | `POISSON_RATIO`         | 0.0                      | $`-`$                          |
| Thickness                          | `THICKNESS`             | 0.4924                   | $`\mathrm{m}`$                 |
| Effective thickness in y-direction | `THICKNESS_EFFECTIVE_Y` | 10.0                     | $`\mathrm{m}`$                 |
| Density                            | `DENSITY`               | 298.36                   | $`\mathrm{kg} / \mathrm{m}^3`$ |

The effective thickness in y-direction accounts for the beam's shear deformation.  It has been given a relatively large value to have a negligible shear deformation.  In this way, the Kratos simulation shows similar beam behavior compared to the one adopted by the comparison FE software package.


## Spring support / anchor parameters

This section documents the derivation of the Kratos input values for the **spring support (anchor)** represented with a **truss element**.

### Given data

From **D-Sheet Piling**:

| Property                | Value             | Unit                                   |
|:------------------------|:------------------|:---------------------------------------|
| Spring support level    | -1.50             | $\mathrm{m}$                           |
| Rotational stiffness    | 0                 | $\mathrm{kNm}/\mathrm{rad}/\mathrm{m}$ |
| Translational stiffness | $1.0 \times 10^4$ | $\mathrm{kN}/\mathrm{m}/\mathrm{m}$    |

Adopted for the truss material:

| Property        | Kratos input parameter | Value                | Unit                       |
|:----------------|:-----------------------|:---------------------|:---------------------------|
| Young's modulus | `YOUNG_MODULUS`        | $2.1 \times 10^{11}$ | $\mathrm{Pa}$              |
| Density         | `DENSITY`              | 0.0                  | $\mathrm{kg}/\mathrm{m}^3$ |
| Prestress       | `TRUSS_PRESTRESS_PK2`  | 0.0                  | $\mathrm{Pa}$              |

### Axial stiffness relation

For a truss element, the axial stiffness is:

$$ k = \frac{EA}{L} $$

where:

* $`k`$ is the axial spring stiffness per meter of wall width
* $`E`$ is the Young's modulus of steel
* $`A`$ is the cross-sectional area per meter of wall width
* $`L`$ is the truss length

### Reference translational stiffness from D-Sheet

Converting the D-Sheet Piling translational distributed spring stiffness to SI units:

```math
k = 1.0 \times 10^4\ \mathrm{kN} / \mathrm{m} / \mathrm{m} = 1.0 \times 10^7\ \mathrm{N} / \mathrm{m} / \mathrm{m}
```

The cross-sectional area of the anchor per meter of wall width can now be calculated as follows:

```math
A = \frac{k \cdot L}{E} = \frac{1.0 \times 10^7\ \mathrm{N} / \mathrm{m} / \mathrm{m} \cdot 1.0\ \mathrm{m}}{2.1 \times 10^{11}\ \mathrm{N} / \mathrm{m}^2} = 4.7619 \times 10^{-5}\ \mathrm{m}^2 / \mathrm{m}
```


### Final values

| Property             | Kratos input parameter | Value                   | Unit                        |
|:---------------------|:-----------------------|:------------------------|:----------------------------|
| Young's modulus      | `YOUNG_MODULUS`        | $2.1 \times 10^{11}$    | $\mathrm{Pa}$               |
| Density              | `DENSITY`              | 0.0                     | $\mathrm{kg}/\mathrm{m}^3$  |
| Cross-sectional area | `CROSS_AREA`           | $4.7619 \times 10^{-5}$ | $\mathrm{m}^2 / \mathrm{m}$ |
| Prestress            | `TRUSS_PRESTRESS_PK2`  | 0.0                     | $\mathrm{Pa}$               |


## Staged analysis

The staged construction analysis consists of seven stages:

1. **Initial stage:**

- The entire soil domain is active. But the anchor, the sheet pile as well as the interfaces at both sides of the soils are inactive.
- Master-slave constraints are applied where the soil will later be separated by the sheet pile. This ensures continuity of the displacement field in this early stage of analysis.
- The only load that is being applied is self-weight.
- At the end of the stage, a $`K_0`$ procedure is performed to initialize the horizontal stress field.  **Note:** the $`K_0`$ procedure requires the use of linear elastic materials for all soil parts.


2. **Null step:**

- This stage is reserved for future test cases, when material models change from linear elastic to, for example, Mohr-Coulomb model. This triggers a stiffness redistribution, and hence a stress redistribution. Currently, the linear elastic model is used for all stages.


3. **Sheet pile installation and applying the external load:**

- Activation of the sheet pile and the interfaces that are attached to its left and right sides.
- Deactivate the master-slave constraints, the interface elements will represent the discontinuity in the displacement at the sheet pile wall location.
- Apply a surface load to a part of the top of the soil on the right-hand side.


4. **First excavation stage:**

- Excavate the top part of the clay to the left of the sheet pile (deactivate corresponding model part).
- Deactivate the interface elements that connect the now excavated clay layer to the sheet pile.
- Apply a normal contact stress to the first exposed part of the sheet pile as well as the bottom of the excavation pit, representing the water in the pit.


5. **Installation of an anchor (spring support):**

- Activate the anchor.


6. **Second excavation stage:**

- Excavate the next top portion of clay to the left of the sheet pile.
- Deactivate corresponding model parts and interface elements.
- Apply a normal contact stress to the second exposed part of the sheet pile as well as the bottom of the excavation pit, representing the water in the pit.


7. **Third excavation stage:**

- Excavate the final top portion of clay to the left of the sheet pile.
- Deactivate corresponding model parts and interface elements.
- Apply a normal contact stress to the exposed part of the sheet pile as well as the bottom of the excavation pit, representing the water in the pit.
