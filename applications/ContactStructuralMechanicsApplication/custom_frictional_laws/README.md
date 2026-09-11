# Frictional laws

Small class hierarchy that returns the friction threshold used by the frictional mortar conditions.

```
FrictionalLaw                                   frictional_law.h/.cpp
└── FrictionalLawWithDerivative<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>
                                                frictional_law_with_derivative.h/.cpp
    ├── CoulombFrictionalLaw<…>                 coulomb_frictional_law.h/.cpp
    └── TrescaFrictionalLaw<…>                  tresca_frictional_law.h/.cpp
```

| Method | Meaning |
|---|---|
| `double GetFrictionCoefficient(const Node&, const PairedCondition&, const ProcessInfo&)` | $\mu$ looked up in the properties, then the node, then the process info (`FRICTION_COEFFICIENT`), 0 otherwise. |
| `double GetThresholdValue(...)` | Maximum admissible tangential traction. Coulomb: $-\mu\,\bar\lambda_n$ with $\bar\lambda_n$ read from `AUGMENTED_NORMAL_CONTACT_PRESSURE`; Tresca: the constant `TRESCA_FRICTION_THRESHOLD`. |
| `double GetDerivativeThresholdValue(..., const DerivativeDataType&, const MortarConditionMatrices&, IndexDerivative, IndexNode)` | Linearisation of the threshold with respect to the displacement and multiplier DoFs (Coulomb uses $\mathbf{D}$, $\mathbf{M}$ and their derivatives; Tresca returns 0). |

Both laws are instantiated for the five geometry pairs (`2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N`) with and without normal variation and exposed to Python as `CoulombFrictionalLaw<suffix>`, `TrescaFrictionalLaw<suffix>` and `FrictionalLaw<suffix>` (suffix `…NV` for normal variation).

**Status (work in progress).** The frictional conditions currently implement Coulomb's law directly in their generated code and read the friction coefficient nodally (`GetFrictionCoefficient()` in the condition headers, with a `TODO` to delegate to a law). The variable `FRICTIONAL_LAW` is registered but not consumed anywhere, and the JSON key `frictional_law` of the contact processes is accepted but only the Coulomb behaviour is wired. Tresca friction is therefore not available to users yet.

The theory of the frictional formulation (Coulomb cone, stick/slip states, thesis §4.3.4) is in [frictional contact](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Frictional_Contact.md); the implementation reference is [Frictional laws and MPC constraint](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Frictional_Laws_And_MPC_Constraint.html) ([source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Frictional_Laws_And_MPC_Constraint.md)).
