# ALM frictionless (components) contact condition

## Description

Symbolic generation of `AugmentedLagrangianMethodFrictionlessComponentsMortarContactCondition<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>`, the augmented Lagrangian frictionless mortar contact condition with a vector Lagrange multiplier whose tangential part is penalised to zero (thesis §4.3.3.2.2). The generated file is `../../custom_conditions/ALM_frictionless_components_mortar_contact_condition.cpp`; the header `ALM_frictionless_components_mortar_contact_condition.h` is hand-written and already stored there. Branches per slave node: inactive / active.

## Files

* `ALM_frictionless_components_mortar_condition.ipynb`: **documented generator** (formulation with the thesis equation numbers, sign conventions, symbol table, and the code that generates the condition).
* `generate_frictionless_components_mortar_condition.py`: the same generator as a command-line script (same functional, same generation call). Keep both in sync.
* `ALM_frictionless_components_mortar_contact_condition_template.cpp`: C++ template with the `// replace_lhs` / `// replace_rhs` markers.
* `alm_frictionlesscomponents_mortar_contact_condition.tex`: theory note (LaTeX).

## Instructions

Requirements: Python 3 and sympy (any modern version). No compiled Kratos is needed.

~~~sh
python3 generate_frictionless_components_mortar_condition.py                       # writes ../../custom_conditions/ALM_frictionless_components_mortar_contact_condition.cpp
python3 generate_frictionless_components_mortar_condition.py --help                # options: --combinations, --normal-variation, --output-dir, --simplify
python3 ../run_notebook.py ALM_frictionless_components_mortar_condition.ipynb   # the same through the notebook, headless
~~~

Then rebuild the application and run the tests. See `../README.md` for the shared machinery (`../mortar_condition_generator.py`) and `../compare_generated_conditions.py` to verify a regeneration numerically.
