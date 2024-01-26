# Transient Thermal Element

## Introduction
The extraction of geothermal energy relies on a profound understanding of heat transport mechanisms within the Earth's subsurface. As heat traverses geological formations, its efficient transfer is governed by mathematical equations encapsulating principles of thermal dynamics. Central to this understanding is the heat transport equation, a fundamental tool that elucidates how thermal energy moves through the Earth's crust, influencing the feasibility and optimization of geothermal applications. In this context, the heat transport equation becomes a linchpin for designing sustainable and efficient geothermal energy systems.

## Governing Equations
The primary governing equation for heat transport in geothermal applications is the heat convection-conduction equation, which describes how heat flows through a porous medium. 

$$ \left(n S \rho^w c^w + (1- n) \rho^s c^s \right) \frac{\partial T}{\partial t} = -\rho^w c^w q_i \frac{\partial T}{\partial x_i} + \frac{\partial}{\partial x_i} \left( D_{ij} \frac{\partial T}{\partial x_j} \right) \qquad \text{on} \quad \Omega $$

where,

- $c^w$       = specific heat capacity liquid phase $\mathrm{[J/kg ^{\circ}C]}$
- $c^s$        = specific heat capacity solid phase  $\mathrm{[J/kg ^{\circ}C]}$
- $D_{ij}$     = hydrodynamic thermal dispersion  $\mathrm{[W/m ^{\circ}C]}$
- $T$           = temperature  $\mathrm{[ ^{\circ}C]}$
- $\rho^s$ = density solid phase $\mathrm{[kg/m^3]}$ 

The hydrodynamic thermal dispersion is defined as:

$$ D_{ij}= nS \lambda^w \delta_{ij} + \left(1-n\right) \lambda_{ij}^s + c^w \rho^w \left( (\alpha_l - \alpha_t) \frac{q_i q_j}{q} + \delta_{ij} \alpha_t q \right) $$

where

- $\alpha_l$    = longitudinal dispersivity $\mathrm{[m]}$
- $\alpha_t$    = transverse dispersivity $\mathrm{[m]}$
- $\delta_{ij}$ = Kronecker delta $\mathrm{[-]}$
- $\lambda^w$   = thermal conductivity water $\mathrm{[W/m ^{\circ}C]}$
- $\lambda^s$   = thermal conductivity solid matrix $\mathrm{[W/m ^{\circ}C]}$

In the absence of ground water flow, these equations are simplified to,

$$ \left(n S \rho^w c^w + (1- n) \rho^s c^s \right) \frac{\partial T}{\partial t} = \frac{\partial}{\partial x_i} \left( D_{ij} \frac{\partial T}{\partial x_j} \right) \qquad \text{on} \quad \Omega $$

$$ D_{ij}= nS \lambda^w \delta_{ij} + \left(1-n\right) \lambda_{ij}^s $$

## Boundary Conditions

### Dirichlet boundary condition
$$ T = \overline T \qquad \text{on} \quad \Gamma_{1}^T $$

where $\overline T$ $\mathrm{\left[ ^{\circ}C \right]}$ is a prescribed temperature 

### Neumann Boundary Condition
$$ D_{ij} \frac{\partial T}{\partial x_j} n_i = \overline{f} \qquad \text{on} \quad \Gamma_{2}^T $$

where $\overline f$ $\mathrm{\left[ W/m^2 \right]}$ is a prescribed conductive heat flux.

### Robin bounday condition
$$ D_{ij} \frac{\partial T}{\partial x_j} n_i = \overline{g} -  \rho^w c^w q_n T \qquad \text{on} \quad \Gamma_{3}^T $$

where $\overline g$ $\mathrm{\left[ W/m^2 \right]}$ is a prescribed convective-conductive heat flux.


## Derived Properties

The density of the bulk material $\rho$ $\mathrm{\left[ kg/m^3 \right]}$ is calculated as,

$$ \rho = n S \rho^w + \left( 1 - n \right) \rho^s $$

And the heat capacity of the bulk material $C$ $\mathrm{\left[ J/m{^3 \circ}C \right]}$ is:

$$ C = n S \rho^w c^w + \left( 1 - n \right) \rho^s c^s $$

The thermal conductivity of the bulk material $\lambda$ $\mathrm{\left[ W/m ^{\circ}C \right]}$

$$ \lambda = n S \lambda^w + \left( 1 - n \right) \lambda^s $$

dynamic viscosity  $\mu$ $[\mathrm {Pas}$] of pure water as a function of temperature is given by:

$$ \mu = 2.4318 \cdot 10^{-5} \cdot 10^{{247.8} / {\left(T+133.0\right]}} $$

density of water $\rho^w$ $[\mathrm {kg/m^3}]$ relates to temperature and Diersch proposes a six order Taylor expansion which is approximated here by:

$$ \rho^w = 9.998396 \cdot 10^2 + 6.764771  \cdot 10^{-2} \cdot T - 8.993699  \cdot 10^{-3} \cdot T^2 + 9.143518 \cdot 10^{-5} \cdot T^3 - 8.907391 \cdot 10^{-7} \cdot T^4 + 5.291959  \cdot 10^{-9} \cdot T^5 - 1.359813  \cdot 10^{-11} \cdot T^6 $$



## Finite Element Formulation

Kratos solves the equations based on incremental method. In The fram of Generelized Newmark method (GN11), the fully implicit incremental temperature formulation read as,

$$ \left(\frac{1}{\theta \Delta t} \boldsymbol{S} + \boldsymbol{A} + \boldsymbol{H} + \boldsymbol{W}^l  \right) \boldsymbol{\Delta T} = \left( \frac{1}{\theta} - 1 \right) \boldsymbol{S} \frac{dT^n}{dt} - \left(\boldsymbol{A} + \boldsymbol{H} + \boldsymbol{W}^l \right) \boldsymbol{T}^{n} + \left( \boldsymbol{V} + \boldsymbol{W}^r \right) $$

### Compressibility matrix 

$$ \boldsymbol{S} = \sum_e \int_{\Omega^e} \left( n S \rho^w c^w + \left(1-n\right) \rho^s c^s \right)^{n+1} \boldsymbol{N}^T  \boldsymbol{N} d \Omega $$

### Convectivity matrix

$$ \boldsymbol{A} = \sum_e \int_{\Omega^e} \left(\rho^w c^w\right)^{n+1}  \boldsymbol{N}^T \boldsymbol{q}^{T,n+1} \boldsymbol{\nabla N}   d \Omega $$

### Conductivity matrix

$$ \boldsymbol{H} = \sum_e \int_{\Omega^e} \boldsymbol{\nabla N}^T \boldsymbol{D}^{n+1} \boldsymbol{\nabla N} d \Omega $$

### Neumann condition (dispersive boundary)

$$ \boldsymbol{V} = \sum_e \int_{\Gamma_2^{ep}}  f^{n+1} \boldsymbol{N}^T  d \Gamma $$

### Robin condition (convective boundary)

$$ \boldsymbol{W^r} = \sum_e \int_{\Gamma_2^{ep}}  g^{n+1} \boldsymbol{N}^T  d \Gamma $$

$$ \boldsymbol{W^l} = \sum_e \int_{\Gamma_3^{ep}}  \left( \rho^w c^w q_n \right)^{n+1} \boldsymbol{N}^T \boldsymbol{I} d \Gamma $$