# Automatic differentiation of the mortar contact conditions

The `CalculateLocalLHS` / `CalculateLocalRHS` methods of the ALM and penalty mortar contact conditions are **generated** by the sympy scripts of this folder: each generator writes the Galerkin functional of one contact state, differentiates it symbolically with respect to the test functions and the degrees of freedom, and prints C++ code into a template. The derivatives of the mortar operators, dual shape functions and normals are *not* differentiated symbolically — they are declared as functions of the DoFs ("AD exceptions") and replaced by the arrays that `DerivativesUtilities` fills at run time (thesis Appendix C).

![Automatic differentiation pipeline](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Theory/images/csma_ad_pipeline.svg)

## Folders

| Folder | Generator | Template | Generated file (in `../custom_conditions/`) | Notes |
|---|---|---|---|---|
| `ALM_frictionless_mortar_condition/` | `generate_frictionless_mortar_condition.py`, `generate_frictionless_mortar_condition_non_zero.py` | `ALM_frictionless_mortar_contact_condition_template.cpp` | `ALM_frictionless_mortar_contact_condition.cpp` | scalar multiplier; theory note `alm_frictionless_mortar_contact_condition.tex` |
| `ALM_frictionless_components_mortar_condition/` | `generate_frictionless_components_mortar_condition.py`, `..._non_zero.py` | `ALM_frictionless_components_mortar_contact_condition_template.cpp` | `ALM_frictionless_components_mortar_contact_condition.cpp` | vector multiplier, tangential part penalised |
| `ALM_frictional_mortar_condition/` | `generate_frictional_mortar_condition.py` | `ALM_frictional_mortar_contact_condition_template.cpp` | `ALM_frictional_mortar_contact_condition.cpp` (~170 k lines) | five branches per node (inactive / slip / stick × objective / non-objective) |
| `penalty_frictionless_mortar_condition/` | `generate_penalty_frictionless_mortar_condition.py`, `..._non_zero.py` | `penalty_frictionless_mortar_contact_condition_template.cpp` | `penalty_frictionless_mortar_contact_condition.cpp` | no multiplier DoFs |
| `penalty_frictional_mortar_condition/` | `generate_penalty_frictional_mortar_condition.py` | `penalty_frictional_mortar_contact_condition_template.cpp` | `penalty_frictional_mortar_contact_condition.cpp` | three branches per node |
| `mesh_tying_mortar_condition/` | `generate_mesh_tying_mortar_condition.py` | – | – | **legacy**: mesh tying is hand-written now; the theory note `mesh_tying_mortar_condition.tex` is still valid |

The `_non_zero` variants emit only the non-zero entries (and zero the local matrices first). The shared helpers are in `../python_scripts/custom_sympy_fe_utilities.py`, on top of the core `KratosMultiphysics.sympy_fe_utilities`.

## How the generation works

1. For every geometry pair (`2D2N`, `3D3N`, `3D4N`, `3D3N4N`, `3D4N3N`), every normal-variation mode (`TNormalVariation` false/true) and every active/inactive (frictional: active/stick/slip) pattern of the slave nodes, the script defines the symbols `u1, u2, X1, X2` (displacements and reference coordinates), the test functions `w1, w2, wLM…`, the multipliers, `NormalSlave`, `DOperator`, `MOperator` (+ `DOperatorold`, `MOperatorold` for friction) and the parameters `DynamicFactor`, `PenaltyParameter`, `ScaleFactor`, `TangentFactor`, `mu`.
2. `DefineDofDependencyMatrix` makes `DOperator`, `MOperator` (and `NormalSlave` when the normal variation is on) depend on the DoFs, so that their derivatives appear as symbols.
3. The functional of the state is assembled (e.g. active node: `DynamicFactor·(ScaleFactor·λn + ε·gn)·n·(D w1 − M w2) + ScaleFactor·gn·wλ`; inactive: `−ScaleFactor²/ε·λn·wλ`) and `Compute_RHS_and_LHS` differentiates it.
4. `OutputMatrix_CollectingFactors` / `OutputVector_CollectingFactors` print C++ with common-sub-expression elimination; `DefineVariableLists` and `SubstituteIndex` rename the derivative symbols to `DeltaDOperator[i]`, `DeltaMOperator[i]`, `DeltaNormalSlave[i]`.
5. The code is inserted at the markers `// replace_lhs` / `// replace_rhs` (between the `BEGIN AD REPLACEMENT` / `END AD REPLACEMENT` banners) of the template; each further combination re-reads the produced file, so the result holds one explicit specialisation per combination, dispatched at run time by `rActiveInactive` (`2^i` encoding of the `ACTIVE` flags, `3^i` for the frictional active/stick/slip states).

## Regenerating a condition

Define dependencies into symbols, is not longer available since Sympy 1.3:

https://github.com/sympy/sympy/wiki/Release-Notes-for-1.3

> Symbols no longer automatically convert to functions when called, e.g., if f = Symbol('f'), f(t) is now a TypeError. To create a function, use f = Function('f') or f = symbols('f', cls=Function).

To solve that temporally you can install the 1.2 version of Sympy: (in order to run the AD scripts contained on this folder)

~~~sh
pip install sympy==1.2
~~~

or

~~~sh
python3 -m pip install sympy==1.2
~~~

Then, with a Kratos installation on the `PYTHONPATH` (the scripts import the application's `custom_sympy_fe_utilities`):

~~~sh
cd <family>/
python3 generate_<family>_mortar_condition.py
cp <family>_mortar_contact_condition.cpp ../../custom_conditions/
~~~

and rebuild the application. The headers in `custom_conditions/` are not generated and must stay consistent with the DoF ordering `[master displacements, slave displacements, multipliers]` and the symbol names used by the scripts. Run the small test suite afterwards: a wrong tangent shows up as a loss of the quadratic convergence in the convergence table.

Note: the `README.md` inside `penalty_frictionless_mortar_condition/` still names the frictionless-ALM script and output; the correct files are the ones in the table above.

## Full documentation

- [Automatic differentiation](https://kratosmultiphysics.github.io/Kratos/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Automatic_Differentiation.html) · [source](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Automatic_Differentiation.md)
- [Linearisation and derivatives](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Theory/Linearisation_And_Derivatives.md) (the externally computed derivatives) · [Conditions](../../../docs/pages/Applications/Contact_Structural_Mechanics_Application/Implementation/Conditions.md)
