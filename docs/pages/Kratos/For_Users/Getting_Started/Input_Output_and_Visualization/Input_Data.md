---
title: Input data
keywords: 
tags: [Input Data]
sidebar: kratos_for_users
summary: 
---

## Overview 

The current input data consist of one archive with .mdpa extension. The input file is free format and the reading is not depend to the spaces, tabs, endlines etc.

The mesh ordering can be consulted in [here](../../../For_Developers/Data_Structures/Mesh_Node_Ordering). Remember that the *Kratos* mesh does not contain directly any geometry, but elements and conditions defined by that geometry.

## MDPA file structure format

This format contains all ModelPart's data in a minimalistic block form:

```cpp
Begin ModelPartData
  //  VARIABLE_NAME value
End ModelPartData

Begin Table table_id variable1 variable2 
  // table_x table_y
End Table
 
Begin Properties  properties_id 
  //  VARIABLE_NAME value
End Properties
 
Begin Nodes
  // id	  X	Y	Z
End Nodes
 
Begin Elements element_name
  // id prop_id	 n1	n2	n3	...
End Elements
 
Begin Conditions condition_name 
  // id prop_id	 n1	n2	n3	...
End Conditions
 
Begin NodalData VARIABLE_NAME
  //  id is_fixed value                 //data accessible by GetSolutionStepValue
End NodalData
 
Begin ElementalData VARIABLE_NAME
  //  id value                          //data accessible by GetValue
End ElementalData
 
Begin ConditionalData VARIABLE_NAME
  //  id value                          //data accessible by GetValue
End ConditionalData 
 
Begin Mesh mesh_id                      // mesh_id cannot be zero!
  Begin MeshData
    //VARIABLE_NAME value               //data accessible by GetValue or operator []
  End MeshData
  Begin MeshNodes
    // node_id
  End MeshNodes
 
  Begin MeshElements
    // element_id
  End MeshElements
 
  Begin MeshConditions
    // condition_id
  End MeshConditions
End Mesh

Begin SubModelPart SubModelPartName
  Begin SubModelPartData
    // VARIABLE_NAME value 
  End SubModelPartData
  Begin SubModelPartTables
    // Table_id
  End SubModelPartTables

  Begin SubModelPartNodes
    // node_id
  End SubModelPartNodes

  Begin SubModelPartElements
    // element_id
  End SubModelPartElements

  Begin SubModelPartConditions
    // condition_id
  End SubModelPartConditions

  Begin SubModelPart SubModelPartName   // Note that this would be a sub sub modelpart
    Begin SubModelPartTables
      // Table_id
    End SubModelPartTables
  
    Begin SubModelPartNodes
      // node_id
    End SubModelPartNodes

    Begin SubModelPartElements
      // element_id
    End SubModelPartElements

    Begin SubModelPartConditions
      // condition_id
    End SubModelPartConditions
  End SubModelPart
End SubModelPart
```

Each block starts with a Begin statement following by the block name and ends with the End statement again following by the block name. Some block may have some additional parameter like id or variable in their definitions.

## ModelPartData Block

### Example
Here is an example of mdpa file:

```cpp
Begin ModelPartData
  //  VARIABLE_NAME value
  AMBIENT_TEMPERATURE 250.00
End ModelPartData
 
Begin Table 1 TEMPERATURE VISCOSITY
  200. 2e-6
  300. 3e-6
  400. 4e-6
End Table

Begin Properties 1
  DENSITY 3.4E-5  //scalar
  THICKNESS 19.5
  VOLUME_ACCELERATION [3] (0.00,0.00,9.8) //vector
  LOCAL_INERTIA [3,3] ((0, 0.27,0.27),(0.087,0,0.27),(0.075,0.23,0)) // matrix

  Begin Table TEMPERATURE VISCOSITY
    200. 2e-6
    300. 3e-6
    400. 4e-6
  End Table

End Properties

Begin Nodes
  1                  16                   0                   0
  2                  16                 0.4                   0
  3                15.6                   0                   0
  972                   0                 7.2                   0
  973                   0                 7.6                   0
  974                   0                   8                   0
End Nodes

Begin Elements Element2D3N
  1 1        1        2        3
  2 1        2        3        972
  3 1        3        972      973
  1796 1     972      973      974
End Elements
 
Begin Conditions Condition2D
  1    1        1          2
  1800 1        2          3
  1801 1        3          972
  1947 1        972        973
  1948 1        973        974
End Conditions

Begin NodalData DISPLACEMENT_X
  1 1 0.100000
  2 1 0.200000
  973 1 0.000000
  974 1 0.000000
End NodalData

Begin NodalData DISPLACEMENT_Y
  1 1 0.000000
  2 1 0.000000
  973 1 0.000973
  974 1 0.000974
End NodalData

Begin NodalData DISPLACEMENT_Z
  1 1 0.000000
  2 1 0.000000
  973 1 0.000000
  974 1 0.000000
End NodalData

Begin NodalData VISCOSITY
  1 0 0.010000
  2 0 0.010000
  973 0 0.010000
  974 0 0.010000
End NodalData

Begin ElementalData TEMPERATURE
  1 3.6
End ElementalData

Begin ElementalData VELOCITY
  1 [3] (1.0, 2.0, 3.0)
End ElementalData

Begin ElementalData CAUCHY_STRESS_TENSOR
  1 [3, 3] ((1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0))
End ElementalData

Begin SubModelPart Inlets
  Begin SubModelPartTables
    1
  End SubModelPartTables

  Begin SubModelPartNodes
    1
    2
  End SubModelPartNodes

  Begin SubModelPartElements
    1
  End SubModelPartElements

  Begin SubModelPartConditions
    1
    1800
  End SubModelPartConditions

  Begin SubModelPart Inlet1
    Begin SubModelPartNodes
      1
      3
    End SubModelPartNodes

    Begin SubModelPartConditions
      1
      1800
    End SubModelPartConditions

  End SubModelPart // Inlet1

  Begin SubModelPart Inlet2
    Begin SubModelPartConditions
      1800
      1801
    End SubModelPartConditions
  End SubModelPart // Inlet2

End SubModelPart // Inlets

Begin SubModelPart Outlet
  Begin SubModelPartProperties
    1
  End SubModelPartProperties

  Begin SubModelPartConditions
    1948
  End SubModelPartConditions

End SubModelPart // Outlet
```

The old format for data file is still supported due to the backward compatibility but without further improvement. The old format description can be found here

## Importing your meshes from a different code

If you want to import your mesh from a different code, there is an initial support for *.mdpa* format file in the utility [meshio](https://github.com/nschloe/meshio). The list of compatible meshes are:

 * [Abaqus](http://abaqus.software.polimi.it/v6.14/index.html)
 * [ANSYS msh](http://www.afs.enea.it/fluent/Public/Fluent-Doc/PDF/chp03.pdf)
 * [DOLFIN XML](http://manpages.ubuntu.com/manpages/wily/man1/dolfin-convert.1.html)
 * [Exodus](https://cubit.sandia.gov/public/13.2/help_manual/WebHelp/finite_element_model/exodus/block_specification.htm)
 * [H5M](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html)
 * [Medit](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html)
 * [MED/Salome](http://docs.salome-platform.org/latest/dev/MEDCoupling/med-file.html)
 * [Gmsh](http://gmsh.info/doc/texinfo/gmsh.html#File-formats)
 * [OFF](http://segeval.cs.princeton.edu/public/off_format.html)
 * [PERMAS](http://www.intes.de)
 * [STL](https://en.wikipedia.org/wiki/STL_(file_format))
 * [VTK](https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf)
 * [VTU](https://www.vtk.org/Wiki/VTK_XML_Formats)
 * [XDMF](http://www.xdmf.org/index.php/XDMF_Model_and_Format)

### Example

Let's say we want to import an VTK mesh.

```python
import meshio

mesh = meshio.read("test.vtk")
meshio.write("test.mdpa", mesh)
```

Additionally, due to the support of *meshio* it is possible to create a mesh using [pygmsh](https://github.com/nschloe/pygmsh). Check the link for examples.
