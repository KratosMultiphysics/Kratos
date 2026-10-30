# Penalty frictional contact condition

## Description

Symbolic generation of `PenaltyMethodFrictionalMortarContactCondition<TDim, TNumNodes, TNormalVariation, TNumNodesMaster>`, the penalty frictional (Coulomb) mortar contact condition, without multiplier DoFs (thesis §4.3.4, eq. 4.71). The generated file is `../../custom_conditions/penalty_frictional_mortar_contact_condition.cpp`; the header `penalty_frictional_mortar_contact_condition.h` is hand-written and already stored there. Branches per slave node: inactive / slip / stick (objective slip only).

## Files

* `penalty_frictional_mortar_condition.ipynb`: **documented generator** (formulation with the thesis equation numbers, sign conventions, symbol table, and the code that generates the condition).
* `generate_penalty_frictional_mortar_condition.py`: the same generator as a command-line script (same functional, same generation call). Keep both in sync.
* `penalty_frictional_mortar_contact_condition_template.cpp`: C++ template with the `// replace_lhs` / `// replace_rhs` markers.

## Instructions

Requirements: Python 3 and sympy (any modern version). No compiled Kratos is needed.

~~~sh
python3 generate_penalty_frictional_mortar_condition.py                       # writes ../../custom_conditions/penalty_frictional_mortar_contact_condition.cpp
python3 generate_penalty_frictional_mortar_condition.py --help                # options: --combinations, --normal-variation, --output-dir, --simplify
python3 ../run_notebook.py penalty_frictional_mortar_condition.ipynb   # the same through the notebook, headless
~~~

Then rebuild the application and run the tests. See `../README.md` for the shared machinery (`../mortar_condition_generator.py`) and `../compare_generated_conditions.py` to verify a regeneration numerically.
