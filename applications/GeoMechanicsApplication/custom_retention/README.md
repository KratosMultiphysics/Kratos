# Retention laws

## Table of Contents

- [Van Genuchten Retention Law](#Van-Genuchten-Retention-Law)
- [Chapter 2](#chapter-2)
- [Chapter 3](#chapter-3)
  
## Van Genuchten Retention Law 

`VanGenuchtenLaw` is a C++ class in the **GeoMechanicsApplication** that implements the
Mualem–Van Genuchten (MVG) unsaturated soil model (Van Genuchten, 1980; Mualem, 1976).
It derives from the `RetentionLaw` base class and computes:

- Degree of saturation $S$
- Effective saturation $S_e$
- Derivative of saturation with respect to fluid pressure $\partial S / \partial p$
- Relative permeability $k_r$
- Bishop's effective stress coefficient $\chi$

### Material Parameters

| Kratos Variable                    | Symbol      | Description                                               | Constraints       |
|------------------------------------|-------------|-----------------------------------------------------------|-------------------|
| `SATURATED_SATURATION`             | $S_s$       | Maximum (fully saturated) degree of saturation            | $[0, 1]$          |
| `RESIDUAL_SATURATION`              | $S_r$       | Residual degree of saturation                             | $[0, S_s)$        |
| `VAN_GENUCHTEN_AIR_ENTRY_PRESSURE` | $p_b$       | Air-entry (bubbling) pressure — inverse of VG scale $\alpha$ | $(0, \infty)$  |
| `VAN_GENUCHTEN_GN`                 | $n$         | Pore-size distribution index                              | $(0, \infty)$     |
| `VAN_GENUCHTEN_GL`                 | $l$         | Mualem pore-connectivity exponent                         | $(-\infty, \infty)$ — see note below |
| `MINIMUM_RELATIVE_PERMEABILITY`    | $k_{r,min}$ | Lower bound for $k_r$ (prevents zero permeability)        | $[0, 1]$          |

> **Note on `VAN_GENUCHTEN_GL` ($l$):** Although Mualem (1976) proposed $l = 0.5$ as a  general value, $l$ is a purely empirical pore-connectivity/tortuosity fitting parameter and **can be negative**. Schaap & Leij (2000) fitted the MVG model to hundreds of soils from the UNSODA database and found optimised values of $l$ ranging from approximately $-4$ to $+4$, with many fine-textured soils yielding $l < 0$. The HYDRUS software (Šimůnek, van Genuchten & Šejna, 2022) imposes **no lower bound** on $l$ and treats it as an unconstrained fitting parameter.

### Sign Convention for Fluid Pressure

Kratos uses the **soil mechanics sign convention** for fluid pressure $p$:

- $p > 0$  → suction (above phreatic level, unsaturated regime)  
- $p \le 0$ → positive pore water pressure (at or below phreatic level, fully saturated)

This is equivalent to the standard hydraulic head formulation where $h = -p/(\rho_w g)$ is the matric suction head in metres.

### Mathematical Formulation (as implemented)

#### Internal variable

$$m = \frac{n - 1}{n}, \qquad gc = \frac{1-n}{n} = -m$$

The constraint $m = 1 - 1/n$ (Mualem–Van Genuchten closure) is **hardcoded**; there is no independent $m$ parameter.

#### Degree of Saturation

$$S = \begin{cases}
S_s & \text{if } p \le 0 \\
S_r + (S_s - S_r)\bigl(1 + (p/p_b)^n\bigr)^{-m} & \text{if } p > 0
\end{cases}$$

#### Effective (Normalised) Saturation

$$S_e = \frac{S - S_r}{S_s - S_r} \in [0,\,1]$$

#### Derivative of Saturation with Respect to Fluid Pressure

$$\frac{\partial S}{\partial p} = \begin{cases}
0 & \text{if } p \le 0 \\[6pt]
(S_s - S_r)\,(-m)\bigl(1+(p/p_b)^n\bigr)^{-m-1} \cdot n\,p_b^{-n}\,p^{n-1} & \text{if } p > 0
\end{cases}$$

> Note: The derivative is **negative** for $p > 0$ since saturation decreases as suction increases.

#### Relative Permeability (Mualem–Van Genuchten)

$$k_r = \max\!\Bigl(k_{r,min},\;
  S_e^{\,l} \Bigl[1 - \bigl(1 - S_e^{1/m}\bigr)^m\Bigr]^2\Bigr)$$

where $1/m = n/(n-1)$.

#### Bishop Effective Stress Coefficient

$$\chi = S_e$$

This is the simplest (and most common) choice for Bishop's parameter.

### Comparison with the Standard Van Genuchten (1980) Model

#### Standard formulation

Van Genuchten (1980) defined the soil water characteristic curve (SWCC) as:

$$\Theta \equiv \frac{\theta - \theta_r}{\theta_s - \theta_r}
  = \Bigl[1 + \bigl(\alpha\,|h|\bigr)^n\Bigr]^{-m}$$

with the Mualem constraint $m = 1 - 1/n$, $\alpha$ [1/length], $h$ = matric suction head (positive).

The Mualem (1976) relative permeability:

$$K_r = \Theta^l \Bigl[1 - \bigl(1-\Theta^{1/m}\bigr)^m\Bigr]^2$$

### Mapping to Kratos variables

| Standard symbol | Kratos variable                    | Relationship            |
|-----------------|------------------------------------|-------------------------|
| $\theta_s$      | `SATURATED_SATURATION`             | $S_s = \theta_s$        |
| $\theta_r$      | `RESIDUAL_SATURATION`              | $S_r = \theta_r$        |
| $\alpha$        | `VAN_GENUCHTEN_AIR_ENTRY_PRESSURE` | $\alpha = 1/p_b$        |
| $n$             | `VAN_GENUCHTEN_GN`                 | $n = \mathtt{gn}$       |
| $m$ (derived)   | —                                  | $m = 1 - 1/n$ (hardcoded) |
| $l$             | `VAN_GENUCHTEN_GL`                 | $l = \mathtt{gl}$       |

#### Equivalence

The Kratos formula is **mathematically identical** to the standard MVG model:

$$S = S_r + (S_s - S_r)\cdot\Theta, \qquad
  k_r = S_e^l\bigl[1-(1-S_e^{1/m})^m\bigr]^2$$

with $\Theta = S_e$ evaluated at $p = \alpha^{-1}\,|h|$.

#### Differences from the standard publication

| Aspect | Standard VG (1980) | Kratos `VanGenuchtenLaw` |
|--------|--------------------|--------------------------|
| Input variable | Matric suction head $h$ [m] | Fluid pressure $p$ [Pa or kPa] |
| Scale parameter | $\alpha = 1/p_b$ | $p_b$ directly |
| $m$ independent? | Optional (some software) | No — $m = 1 - 1/n$ is hardcoded |
| Saturated branch | Not defined separately | $S = S_s$ for $p \le 0$ |
| $S_s$ value | Usually $\theta_s = 1$ | Any value $\in (0,1]$ |
| Min permeability clamp | Not in original | Yes — `MINIMUM_RELATIVE_PERMEABILITY` |

### Known Limitations and Potential Issues

1. **Discontinuous derivative at $p = 0$**: `CalculateDerivativeOfSaturation` returns 0    for $p \le 0$ but is non-zero just above $p = 0$. This C⁰ discontinuity may hinder quadratic convergence of Newton iterations near the phreatic surface.

2. **No independent $m$ parameter**: Soils whose measured SWCC cannot be captured by the MVG constraint $m = 1 - 1/n$ cannot be represented. This is a restriction compared to some literature variants.

3. **$n > 1$ required for finite permeability**: If `VAN_GENUCHTEN_GN` $\le 1$, the exponent $1/m = n/(n-1)$ is undefined or negative. The `Check()` method only validates $n > 0$, not $n > 1$. Passing $n \le 1$ will produce `NaN` or erroneous permeability values.

4. **Minimum relative permeability clamp**: The `MINIMUM_RELATIVE_PERMEABILITY` floor prevents zero permeability but also means $k_r$ is not differentiable at the clamped region. This is consistent with FEM practice.

5. **Bishop coefficient**: $\chi = S_e$ is a simplification. More refined models (e.g., Khalili & Khabbaz, 1998) use power-law corrections. Only the simplest linear form is    implemented.

### References

- Van Genuchten, M. Th. (1980). *A closed-form equation for predicting the hydraulic conductivity of unsaturated soils.* Soil Sci. Soc. Am. J., 44(5):892–898.
- Mualem, Y. (1976). *A new model for predicting the hydraulic conductivity of unsaturated porous media.* Water Resour. Res., 12(3):513–522.
- Brinkgreve, R. B. J. et al. (2006). *PLAXIS Geotechnical Code for Soil and Rock Analyses.* Delft, Netherlands.
- Khalili, N. & Khabbaz, M. H. (1998). *A unique relationship for χ for the determination of the shear strength of unsaturated soils.* Géotechnique, 48(5):681–687.
- Bishop, A. W. (1959). *The principle of effective stress.* Teknisk Ukeblad, 106(39):859–863.
- Thom, C. S., Springman, S. M., Rouainia, M., van Esch, J. M. et al. (2017). *Numerical modelling of slope–vegetation–atmosphere interaction: an overview.* Quarterly Journal of Engineering Geology and Hydrogeology, 50(3):249–270. DOI: 10.1144/qjegh2016-079
- van Esch, J. M. (2012). *Modeling groundwater flow through dikes for real time stability assessment.* Proc. Computational Water Resources (CMWR), Urbana-Champaign, Illinois.
- Schaap, M. G. & Leij, F. J. (2000). *Improved prediction of unsaturated hydraulic conductivity with the Mualem–van Genuchten model.* Soil Sci. Soc. Am. J., 64(3):843–851.
  DOI: 10.2136/sssaj2000.643843x — reports optimised $l$ values ranging from $\approx -4$ to $+4$ across the UNSODA database; many fine-textured soils have $l < 0$.
- Šimůnek, J., van Genuchten, M. Th., & Šejna, M. (2022). *The HYDRUS Software Package for Simulating the One-, Two-, and Three-Dimensional Movement of Water, Heat, and Multiple Solutes in Variably-Saturated Media, Technical Manual I, Version 5.0.* PC Progress, Prague — treats $l$ as an unconstrained real-valued fitting parameter.