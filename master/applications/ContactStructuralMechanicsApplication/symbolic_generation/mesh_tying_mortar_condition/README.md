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
Then  file "_mesh_tying_mortar_condition.cpp_" is generated automatically. Such file should be copied within the "_custom_conditions_" folder of the
**ContactStructuralMechanicsApplication**. The corresponding header file "_mesh_tying_mortar_condition.h_", which implements the element is already stored in the custom_elements folder.
