# Mesh-tying condition (LEGACY FILES)

## LEGACY NOTE:

Here only the generation files are preserved. AD has been removed from MeshTying and now is manually constructed for more generality.

## ELEMENT DESCRIPTION:
Current directory contains the documentation for the symbolic derivation of the _"mesh_tying"_ condition. This element includes a formulation of a mesh tying condition using mortar formulation.

## SYMBOLIC GENERATOR SETTINGS:
* Nothing to add

## INSTRUCTIONS:
Run:
~~~py
python generate_mesh_tying_mortar_condition.py
~~~
Then the file "_mesh_tying_mortar_condition.cpp_" would be generated; the template it expects (`mesh_tying_mortar_condition_template.cpp`) is no longer kept, so this script cannot run any more and is preserved only as a record of the original derivation. The hand-written condition lives in `custom_conditions/mesh_tying_mortar_condition.h/.cpp`.
