---
title: Glossary
keywords: glossary, acronyms, symbols, ALM, mortar, dual Lagrange multiplier, weighted gap, augmented pressure, active set, KKT, NCP
tags: [reference, glossary, acronyms, notation]
sidebar: contact_structural_mechanics_application
summary: Acronyms, mathematical symbols and code-specific terms used throughout the documentation, taken from the glossaries of the thesis (pp. 359–371) and from the source code.
---

> **Sources.** Thesis glossaries "Acronyms" (pp. 359–363), "Mathematical symbols", "Constraint enforcement and optimization" (p. 365), "Mortar formulation" (p. 368), "Contact mechanics" (pp. 369–371), "Adaptive remeshing" (pp. 371–373); the variables and names of the code are described in [Variables and flags](../Implementation/Variables_And_Flags_Reference.html).

## Acronyms

| Acronym | Meaning | Where it appears |
|---|---|---|
| AABB | Axis-Aligned Bounding Box | search broad phase (`in_box`) |
| AALM | Adapted Augmented Lagrangian Method (Bussetta et al.) | `AALMAdaptPenaltyValueProcess`, `adapt_penalty` |
| AD | Automatic Differentiation | `automatic_differentiation/` |
| ADLM | Augmented Dual Lagrange Multiplier (the thesis' name for the ALM + dual-multiplier formulation) | ALM conditions |
| ALM | Augmented Lagrangian Method | `ALMContact*` formulations |
| AMG | Algebraic MultiGrid | AMGCL inner solver of the `MixedULMLinearSolver` |
| APM | Adapted Penalty Method | thesis App. D.2.3.1 |
| BC | Boundary Conditions | |
| BS | Bounding Spheres | `in_radius` search |
| BVH | Bounding Volume Hierarchies | KD-tree / octree search |
| BVP / IBVP | (Initial) Boundary Value Problem | strong formulation |
| CAS | Computer Algebra System | sympy generators |
| CCM | Computational Contact Mechanics | |
| CD | Collision Detection | |
| CDM | Contact Domain Method | thesis §4.2.2.3 (not implemented) |
| CDT | Constrained Delaunay Triangulation | mortar segmentation alternatives |
| CL | Constitutive Law | |
| COF | Coefficient Of Friction, $$\mu$$ | `FRICTION_COEFFICIENT` |
| CPPM | Closest Point Projection Method | mortar projection |
| DDM | Domain Decomposition Method | origin of mortar methods |
| DLMM | Double (dual) Lagrange Multiplier Method | dual shape functions $$\Phi_j$$ |
| DOF | Degree Of Freedom | |
| DOP / k-DOP | Discretised Orientation Polytopes | bounding volumes (`kdop` not implemented) |
| FAD | Forward Automatic Differentiation | |
| FE / FEM / FEA | Finite Element (Method / Analysis) | |
| GP | Gauss Point | `INTEGRATION_ORDER_CONTACT` |
| HSM | Hertz–Signorini–Moreau conditions (= KKT conditions of contact) | thesis eq. 4.3 |
| IGA | Isogeometric Analysis | thesis §4.2.2.5.1 (not implemented) |
| KKT | Karush–Kuhn–Tucker conditions | $$g_n \ge 0,\ p_n \le 0,\ p_n g_n = 0$$ |
| LM / LMM | Lagrange Multiplier (Method) | |
| MFC / MPC | MultiFreedom / MultiPoint Constraint | `MPCMortarContactCondition`, `ContactMasterSlaveConstraint` |
| MMG | The Mmg remeshing library | `contact_remesh_mmg_process.py` |
| NCP | Non-linear Complementarity Problem (function) | semi-smooth Newton, thesis eqs. 4.44, 4.79 |
| NL | Non-Linear | |
| NR | Newton–Raphson | |
| NTN / NTS / STS | Node-To-Node / Node-To-Segment / Segment-To-Segment discretisations | thesis §4.2.2; mortar = STS |
| OBB | Oriented Bounding Box | `OrientedBoundingBox`, `*_with_obb` |
| PDASS | Primal–Dual Active Set Strategy | `ActiveSetUtilities` |
| PDE | Partial Differential Equation | |
| PM | Penalty Method | `PenaltyContact*` formulations |
| SAT | Separating Axis Theorem | `OBB_intersection_type: SeparatingAxisTheorem` |
| SPR | Superconvergent Patch Recovery | `ContactSPRErrorProcess` |
| SVD | Singular Value Decomposition | condition number study (thesis §4.3.3.3) |
| TL / UL | Total / Updated Lagrangian | structural elements |
| VM | Von Mises (stress) | remeshing metric |

## Symbols of the contact formulation

| Symbol | Meaning | Code counterpart |
|---|---|---|
| $$\Omega^{(i)}$$, $$\Gamma_u^{(i)}$$, $$\Gamma_\sigma^{(i)}$$, $$\Gamma_c^{(i)}$$ | Body $$i$$ ($$1$$ = slave, $$2$$ = master) and its Dirichlet, Neumann and contact boundaries in the reference configuration; $$\gamma$$ denotes the current configuration | sub-model-parts `Contact`, `SlaveSubModelPart<k>`, `MasterSubModelPart<k>` |
| $$\mathbf{n}$$, $$\boldsymbol\tau_1$$, $$\boldsymbol\tau_2$$ | Slave unit normal and tangents of the local contact frame | `NORMAL`, `MortarUtilities::ComputeTangentMatrix` |
| $$\chi$$, $$\chi_h$$ | (Discrete) interface mapping: projection of slave points onto the master along $$\mathbf{n}$$ | `ExactMortarIntegrationUtility`, `MasterShapeFunctionValue` |
| $$g_n$$ | Normal gap $$g_n = \mathbf{n}\cdot(\mathbf{x}^{(1)} - \hat{\mathbf{x}}^{(2)})$$ (positive when open) | `NORMAL_GAP` (mapped) |
| $$\tilde g_n$$ | Weighted (mortar-integrated) nodal gap $$\tilde g_n = \mathbf{n}\cdot(\mathbf{D}\mathbf{x}^{(1)} - \mathbf{M}\mathbf{x}^{(2)})$$ | `WEIGHTED_GAP` |
| $$\tilde{\mathbf{g}}_\tau$$, $$\tilde{\mathbf{u}}_\tau$$, $$\tilde{\mathbf{v}}_\tau$$ | Weighted tangential slip increment and relative velocity (objective / non-objective measures) | `WEIGHTED_SLIP`, `OPERATOR_THRESHOLD` |
| $$p_n$$, $$\mathbf{t}_c$$, $$t_{co}^n$$, $$t_{co}^\tau$$ | Contact pressure and interface traction (normal / tangential) | |
| $$\lambda_n$$, $$\boldsymbol\lambda$$, $$\boldsymbol\lambda_\tau$$ | Scalar normal multiplier (= $$-p_n$$), vector multiplier and its tangential part | `LAGRANGE_MULTIPLIER_CONTACT_PRESSURE`, `VECTOR_LAGRANGE_MULTIPLIER` |
| $$\bar\lambda_n = k\lambda_n + \varepsilon\tilde g_n$$ | Augmented normal pressure (the frictionless NCP function up to the $$\max$$) | `AUGMENTED_NORMAL_CONTACT_PRESSURE` |
| $$\bar{\boldsymbol\lambda}_\tau$$ | Augmented tangential traction | `AUGMENTED_TANGENT_CONTACT_PRESSURE` |
| $$\varepsilon$$, $$\varepsilon_n$$, $$\varepsilon_\tau$$ | Penalty parameter (normal, tangential $$\varepsilon_\tau = $$ `TANGENT_FACTOR` $$\cdot\varepsilon$$) | `INITIAL_PENALTY`, `TANGENT_FACTOR` |
| $$k$$ | Scale factor of the ALM functional (conditioning only) | `SCALE_FACTOR` |
| $$\mu$$ | Coulomb friction coefficient | `FRICTION_COEFFICIENT` |
| $$\mathcal{F}$$, $$g$$ (Tresca) | Friction threshold ($$\mu\vert p_n\vert$$ for Coulomb, constant for Tresca) | `FrictionalLaw::GetThresholdValue`, `TRESCA_FRICTION_THRESHOLD` |
| $$\beta$$ | Velocity–traction ratio of Coulomb's law (stick $$\beta = 0$$, slip $$\beta \gt 0$$) | `SLIP` flag |
| $$\langle\cdot\rangle$$ | Macaulay bracket | ALM functional (thesis eq. 4.10) |
| $$\mathcal{L}$$, $$\mathcal{L}_\lambda$$, $$\mathcal{L}_{\bar\lambda}$$, $$f_p$$ | Lagrangian functionals of the LMM / ALM and penalised function of the PM | generator scripts (`rv_galerkin`) |
| $$\mathcal{U}$$, $$\mathcal{V}$$, $$\mathcal{M}$$, $$\mathcal{M}_h$$ | Solution and weighting spaces of the displacements and (discrete) multipliers | |
| $$C_{\lambda_n}$$, $$C_\tau$$ | NCP functions of the normal and tangential problems (thesis eqs. 4.44, 4.79) | `ActiveSetUtilities` |
| $$\mathcal{A}$$, $$\mathcal{I}$$, $$\mathcal{S}$$, $$\mathcal{M}$$, $$\mathcal{N}$$ | Active / inactive slave sets, slave / master / remaining DoF sets | `ACTIVE` flag, `MixedULMLinearSolver::BlockType` |
| $$\mathrm{sl}$$, $$\mathrm{st}$$ | Slip / stick subsets of the active set | `SLIP` flag |

## Symbols of the mortar formulation

| Symbol | Meaning | Code counterpart |
|---|---|---|
| $$N_k^{(1)}$$, $$N_l^{(2)}$$ | Standard shape functions of slave and master | geometry `ShapeFunctionsValues` |
| $$\Phi_j$$ | Dual (biorthogonal) Lagrange-multiplier shape functions | `DualLagrangeMultiplierOperators` |
| $$\mathbf{A}_e = \mathbf{D}_e\mathbf{M}_e^{-1}$$ | Coefficient matrix of the dual functions ($$\Phi_j = a_{jk}N_k$$) | `CalculateAe`, `CalculateDe` |
| $$\mathbf{D}$$, $$\mathbf{M}$$ | Mortar operators (slave–slave, diagonal with dual multipliers; slave–master) | `MortarOperator::DOperator`, `MOperator` |
| $$\Delta\mathbf{D}$$, $$\Delta\mathbf{M}$$, $$\Delta\mathbf{A}_e$$, $$\Delta\mathbf{n}$$ | Directional derivatives of the operators, dual coefficients and normals | `MortarOperatorWithDerivatives`, `DerivativesUtilities` |
| $$\mathbf{P} = \mathbf{D}^{-1}\mathbf{M}$$ | Mortar projection operator (relation matrix of the MPC route) | `CalculatePOperator`, `UpdateConstraint*` |
| $$\mathbf{B}_{co}$$, $$\mathbf{B}_{mt}$$ | Discrete mortar contact / tying operators $$[\mathbf{0}, -\mathbf{M}^T, \mathbf{D}^T]$$ | assembled by the conditions |
| $$\xi_a^{1}, \xi_b^{1}, \xi_a^{2}, \xi_b^{2}$$ | Local coordinates of the integration segment ends on slave and master (2D) | `ExactMortarIntegrationUtility` |
| $$\mathbf{x}_{clip}$$, $$\hat{\mathbf{x}}^{i}_{j}$$, $$\mathbf{n}_{plane}$$, $$J_{clip}$$, $$\bar N$$ | Clipping points, projected nodes, auxiliary plane normal, Jacobian and shape functions of the 3D integration cells | idem, `CalculateDeltaCellVertex` |
| $$w_g$$, $$J_g$$ | Gauss weights and Jacobians of the integration cells | `INTEGRATION_ORDER_CONTACT` |
| $$m^{(1)}$$, $$n^{(1)}$$, $$n^{(2)}$$ | Number of multiplier nodes, slave nodes and master nodes | `TNumNodes`, `TNumNodesMaster` |

## Symbols of the search

| Symbol | Meaning | Code counterpart |
|---|---|---|
| $$\mathbf{C}$$, $$\mathbf{A}_i$$, $$a_i$$ | Centre, axes and half-extents of an oriented bounding box (thesis eq. 4.80) | `OrientedBoundingBox` |
| $$r$$ | Search radius (multiple of the condition size) | `search_factor` × `NODAL_H` |
| $$h$$, $$h_{mean}$$ | Element / mean element size | `NODAL_H` |
| $$E_{mean}$$ | Mean Young modulus of the interface | `ALMVariablesCalculationProcess` |

## Code-specific terms

| Term | Meaning |
|---|---|
| **Slave / master** | The slave side is where the mortar integration and the multipliers live (Popp's "non-mortar" side); the master side is projected onto it. Chosen with `assume_master_slave`. Note that `PairedCondition::GetParentGeometry()` returns the slave and `GetPairedGeometry()` the master. |
| **Pair / paired condition** | One slave condition coupled with one master condition (`PairedCondition` with a `CouplingGeometry`), created by the search in `ComputingContact`. |
| **`mortar_type`** | The solver key that selects the formulation (`ALMContactFrictionless`, `ALMContactFrictionlessComponents`, `ALMContactFrictional[PureSlip]`, `PenaltyContactFrictionless`, `PenaltyContactFrictional[PureSlip]`, `ScalarMeshTying`, `ComponentsMeshTying`). |
| **Components formulation** | Frictionless ALM with a vector multiplier whose tangential part is penalised to zero; allows the static condensation of the multipliers. |
| **NV** | "Normal Variation" suffix of the conditions whose tangent includes the derivatives of the slave normals (`normal_variation: nodal_elemental_derivatives`). |
| **Active set** | The slave nodes currently in contact (`ACTIVE` flag), updated after every Newton iteration from the sign of $$\bar\lambda_n$$. |
| **Semi-smooth Newton** | Single Newton loop that treats the active-set change as one more non-linearity (thesis Algorithms 2–3); the *simplified* variant freezes the sets inside an inner loop (`simplified_semi_smooth_newton`, `INTERACTION` flag). |
| **Weighted gap / weighted slip** | Mortar-integrated (nodal) gap and slip, `WEIGHTED_GAP` / `WEIGHTED_SLIP`; not lengths but integrals over the nodal support. |
| **Augmented pressure** | $$\bar\lambda_n$$, the effective contact pressure and active-set indicator (`AUGMENTED_NORMAL_CONTACT_PRESSURE`). |
| **Isolated node** | A slave node whose pairs have all been removed; its multiplier DoFs are fixed by the block builder (`ISOLATED` flag). |
| **Explicit contribution** | Residual-only evaluation of a pair (`AddExplicitContributionOfMortarCondition`) used to refresh the weighted gap in `Predict()` and in the criteria, and by the explicit solver. |
| **Dynamic factor** | `DYNAMIC_FACTOR`, a nodal scaling of the contact contribution derived from the gap history (dynamic problems). |
| **Consistent gap** | Gap computed by mapping the master surface onto the slave with the mortar mapper (`check_gap: check_mapping`, thesis Algorithm 8). |
| **AD exception** | A quantity whose derivative is supplied externally to the symbolic differentiation (the mortar operators and normals); see [Automatic differentiation](../Theory/Automatic_Differentiation.html). |
| **`ComputingContact`** | Sub-model-part holding the pair conditions that are assembled; `Contact` holds the interface conditions and nodes. |

See also the [Bibliography](Bibliography.html).
