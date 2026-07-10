---
title: MDPA file structure format
keywords: mesh  io  input  mdpa
tags: [MDPA file structure format]
sidebar: kratos_for_developers
summary: The current input data consist of one archive with `.mdpa` extension. The input file is free format and the reading is not depend to the spaces, tabs, endlines etc.
---

# Overview 

The current input data consist of one archive with `.mdpa` extension. The input file is free format and the reading is not depend to the spaces, tabs, endlines etc.

The mesh ordering can be consulted in [here](../Data_Structures/Mesh_Node_Ordering.md). Remember that the *Kratos* mesh does not contain directly any geometry, but elements and conditions defined by that geometry.

# MDPA file structure format

This format contains all ModelPart's data in a minimalistic block form:

```c++
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

Begin Geometries	geometry_name
	// id	 n1	n2	n3	...
End Geometries

Begin Elements element_name
  // id prop_id	 n1	n2	n3	...
End Elements
 
Begin Conditions condition_name 
  // id prop_id	 n1	n2	n3	...
End Conditions
 
Begin Constraints	constraint_name	dependent_variable independent_variables
	// id const_value vector_of_relation id_dependent_node id_indepents_nodes
End Constraints

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

  Begin SubModelPartGeometries
    // geometry_ids
  End SubModelPartGeometries
  
  Begin SubModelPartConstraints
    // constraint_id
  End SubModelPartConstraints

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

    Begin SubModelPartGeometries
      // geometry_ids
    End SubModelPartGeometries
    
    Begin SubModelPartConstraints
      // constraint_id
    End SubModelPartConstraints
  End SubModelPart
End SubModelPart
```

## Block details

The **MDPA (Model Part Data)** file format as shown employs a block-based structure to organize all the necessary information for a simulation, including mesh data, material properties, boundary conditions, and solver settings. Each block is delineated by `Begin` and `End` statements, ensuring a clear and organized file structure. Some block may have some additional parameter like id or variable in their definitions. In detal:

### **ModelPartData**

This block is used to define global variables and parameters for the entire model part. These can include settings for the simulation, such as time step, or flags to control solver behavior.

```c++
Begin ModelPartData
  // VARIABLE_NAME value
  TIME_STEPS 1000
  SOLVER_TOLERANCE 1e-6
End ModelPartData
```

* **VARIABLE_NAME**: The name of the global parameter.
* **value**: The value assigned to the parameter.

---

### **Table**

Tables are used to define time-dependent data or other functional relationships. For instance, a table can specify how a load varies over time.

```c++
Begin Table table_id variable1 variable2 
  // table_x table_y
  0.0 0.0
  1.0 10.0
  2.0 5.0
End Table
```

* **table_id**: A unique integer identifier for the table.
* **variable1, variable2**: The names of the variables represented by the columns in the table (e.g., TIME LOAD).
* **table_x, table_y**: The data points defining the relationship.

---

### **Properties**

This block defines the material properties and other physical characteristics that will be assigned to elements. Each `Properties` block is given a unique ID so that it can be referenced by elements.

```c++
Begin Properties 1
  // VARIABLE_NAME value
  DENSITY 7850.0
  YOUNG_MODULUS 2.1e11
  POISSON_RATIO 0.3
End Properties
```

* **properties_id**: A unique integer identifier for this set of properties.
* **VARIABLE_NAME**: The name of the material property (e.g., DENSITY, YOUNG_MODULUS).
* **value**: The value of the material property.

---

### **Nodes**

The `Nodes` block defines the spatial coordinates of the points that make up the mesh.

```c++
Begin Nodes
  // id   X Y Z
  1   0.0 0.0 0.0
  2   1.0 0.0 0.0
  3   1.0 1.0 0.0
End Nodes
```

* **id**: A unique integer identifier for the node.
* **X, Y, Z**: The coordinates of the node in 3D space. For 2D simulations, the Z coordinate is often set to 0.

---

### **Geometries**

Geometries define geometric entities (like lines, surfaces, or volumes) based on a collection of nodes. This block is often used for pre-processing and is not always directly used by the solver.

```c++
Begin Geometries geometry_name
  // id  n1 n2  n3  ...
  1   1 2 3
End Geometries
```

* **geometry_name**: The name of the geometry type (e.g., `Triangle2D3N`).
* **id**: A unique integer identifier for the geometry.
* **n1, n2, n3, ...**: The IDs of the nodes that form the geometry.

---

### **Elements**

Elements are the fundamental building blocks of the computational domain. They connect nodes to form the mesh and have assigned material properties.

```c++
Begin Elements element_name
  // id prop_id  n1 n2  n3  ...
  1 1 1 2 3
End Elements
```

* **element_name**: The name of the element type (e.g., `ShellThin`, `Solid3D4N`).
* **id**: A unique integer identifier for the element.
* **prop_id**: The ID of the `Properties` block that defines the material for this element.
* **n1, n2, n3, ...**: The IDs of the nodes that form the element's connectivity.

---

### **Conditions**

Conditions are used to apply loads, boundary conditions, and other constraints to the model. They are similar to elements but typically represent surfaces or points where external interactions occur.

```c++
Begin Conditions condition_name 
  // id prop_id  n1 n2  n3  ...
  1 1 4 5 6
End Conditions
```

* **condition_name**: The name of the condition type (e.g., `PointLoad`, `SurfacePressure`).
* **id**: A unique integer identifier for the condition.
* **prop_id**: The ID of the `Properties` block associated with this condition (if any).
* **n1, n2, n3, ...**: The IDs of the nodes to which the condition is applied.

---

### **Constraints**

This block defines kinematic constraints between degrees of freedom in the model, such as rigid links or prescribed relationships between nodal displacements.

```c++
Begin Constraints constraint_name dependent_variable independent_variables
  // id const_value vector_of_relation id_dependent_node id_indepents_nodes
  1 0.0 [2,1.0,-0.5] 10 12 15
End Constraints
```

* **constraint_name**: A name for the constraint type.
* **dependent_variable**: The degree of freedom of the dependent node.
* **independent_variables**: The degrees of freedom of the independent nodes.
* **id**: A unique integer identifier for the constraint.
* **const_value**: A constant value in the constraint equation.
* **vector_of_relation**: The coefficients that define the linear relationship between the dependent and independent variables.
* **id_dependent_node**: The ID of the node whose degree of freedom is constrained.
* **id_indepents_nodes**: The IDs of the nodes that constrain the dependent node.

---

### **NodalData**

This block assigns initial conditions or prescribed values to specific degrees of freedom of the nodes.

```c++
Begin NodalData VARIABLE_NAME
  //  id is_fixed value
  1   1   0.0      // This node's variable is fixed to 0.0
  5   0   10.5     // This node's variable has an initial value of 10.5
End NodalData
```

* **VARIABLE_NAME**: The name of the nodal variable (e.g., `DISPLACEMENT_X`, `TEMPERATURE`).
* **id**: The ID of the node.
* **is_fixed**: A boolean (1 for true, 0 for false) indicating whether the degree of freedom is fixed to the specified value. If `is_fixed` is 1, the solver will enforce this value throughout the analysis.
* **value**: The initial or prescribed value for the variable.

---

### **ElementalData**

Used to assign specific data or initial values to elements, which are not part of the material properties.

```c++
Begin ElementalData VARIABLE_NAME
  //  id value
  101  1.25
End ElementalData
```

* **VARIABLE_NAME**: The name of the variable to be assigned to the elements (e.g., `INITIAL_STRESS`).
* **id**: The ID of the element.
* **value**: The value to be assigned to the element's variable.

---

### **ConditionalData**

Similar to `ElementalData`, this block assigns specific data to conditions.

```c++
Begin ConditionalData VARIABLE_NAME
  //  id value
  201  -9.81
End ConditionalData
```

* **VARIABLE_NAME**: The name of the variable associated with the conditions (e.g., `PRESSURE`).
* **id**: The ID of the condition.
* **value**: The value to be assigned to the condition's variable.

---

### **Mesh**

The `Mesh` block allows for the organization of the model into different sub-regions or meshes. This is particularly useful for complex models or multi-physics simulations.

```c++
Begin Mesh 1
  Begin MeshData
    //VARIABLE_NAME value
  End MeshData
  Begin MeshNodes
    1
    2
  End MeshNodes
  Begin MeshElements
    1
  End MeshElements
  Begin MeshConditions
    1
  End MeshConditions
End Mesh
```

* **mesh_id**: A non-zero integer identifier for the mesh.
* **MeshData**: Contains data specific to this mesh.
* **MeshNodes**, **MeshElements**, **MeshConditions**: Lists the IDs of the nodes, elements, and conditions that belong to this mesh.

---

### **SubModelPart**

This provides a powerful way to create a hierarchical structure within the model. A `SubModelPart` can contain its own data, tables, nodes, elements, and even other `SubModelPart`s, allowing for a nested organization of the simulation domain. This is useful for defining complex boundary conditions, contact pairs, or different material regions.

```c++
Begin SubModelPart PartName
  Begin SubModelPartData
    // VARIABLE_NAME value 
  End SubModelPartData
  Begin SubModelPartNodes
    // node_id
    10
    11
  End SubModelPartNodes
  Begin SubModelPartElements
    // element_id
    101
  End SubModelPartElements
  Begin SubModelPart SubSubPartName
    Begin SubModelPartNodes
      // node_id
      10
    End SubModelPartNodes
  End SubModelPart
End SubModelPart
```

* **SubModelPartName**: The name of the sub-model part.
* The blocks within `SubModelPart` serve to define the components (nodes, elements, etc.) that belong to this specific sub-region of the model. The hierarchical nature allows for a very detailed and organized model definition.

## Example
Here is an example of `mdpa` file:

```c++
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
  VOLUME_ACCELERATION  (0.00,0.00,9.8) //vector
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

Begin Geometries	Triangle2D3
  1        1        2        3
  2        2        3        972
  3        3        972      973
  1796     972      973      974
End Geometries

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

Begin Constraints	LinearMasterSlaveConstraint	DISPLACEMENT_X	DISPLACEMENT_X
	1	0.0	[5.0e-01]	2	1	
	2	0.0	[5.0e-01]	3	1	
End Constraints

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
  1  (1.0, 2.0, 3.0)
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

  Begin SubModelPartGeometries
    // No geometry added
  End SubModelPartGeometries
  
  Begin SubModelPartConstraints
    1
  End SubModelPartConstraints

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

The old format for data file is still supported due to the backward compatibility but without further improvement.

# Importing your meshes from a different code

If you want to import your mesh from a different code, there is an initial support for *.mdpa* format file in the utility [meshio](https://github.com/nschloe/meshio).

## File Formats

The list of compatible meshes are:

### **Abaqus (.inp)**

  * **Description:** The Abaqus input file (`.inp`) is a text-based file containing all instructions, data, and parameters to define and run a simulation in the Abaqus Finite Element Analysis (FEA) software. It uses a keyword-based syntax to define model geometry, materials, properties, analysis steps, loads, and output requests. Some advanced simulations might only be executable via `.inp` files.
  * **Link:** While a single official specification URL is not provided, general information can be found through resources like [CAE Assistant](https://caeassistant.com/blog/how-to-open-run-input-file-abaqus-video/) and [Engssoft](https://www.engssoft.com/abaqus-inp-files/).

### **ANSYS msh (.msh)**

  * **Description:** The ANSYS `.msh` file format is primarily used by ANSYS Fluent and other ANSYS meshing tools. It stores mesh data, including node coordinates, connectivity, and zone information (e.g., wall, fluid) for Computational Fluid Dynamics (CFD) and other analyses. It's essentially a subset of an ANSYS Fluent case file (`.cas`).
  * **Link:** Specific details are typically within ANSYS documentation. General information can be found at [ANSYS Help](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/wb_msh/msh_export_fluent.html).

### **AVS-UCD (.avs,.inp)**

  * **Description:** The AVS UCD (Unstructured Cell Data) format stores unstructured cell data (points, lines, triangles, quads, tetrahedra, etc.) and associated scalar or vector data for visualization, commonly in structural analysis and CFD. It can be ASCII or binary. Note: the `.inp` extension can conflict with Abaqus; `.ucd` is also common.
  * **Link:** Documentation snippets can be found at [LBL Archive](https://dav.lbl.gov/archive/NERSC/Software/express/help6.1/help/reference/dvmac/UCD_Form.htm) and [Debian Manpages](https://debian.man.ac.uk/.f/pub/cgu/iavsc/express/source/data_io/wr_ucd/wr_ucd.html).

### **CGNS (.cgns)**

  * **Description:** The CFD General Notation System (CGNS) is a standard for storing and retrieving Computational Fluid Dynamics (CFD) analysis data. It's designed to be general, portable, extensible, and self-descriptive, facilitating data exchange and archiving. It uses either ADF or HDF5 as its low-level data format.
  * **Link:** The official website is [CGNS.org](http://www.cgns.org). Further documentation can be found at [CGNS GitHub Docs](http://cgns.github.io/CGNS_docs_current/sids/index.html).

### **DOLFIN XML (.xml)**

  * **Description:** This XML-based format is used by DOLFIN, part of the FEniCS Project, for solving partial differential equations (PDEs) with the finite element method (FEM). It stores mesh geometry (vertices, cells) and function data (solutions). PyVista can read these files via `meshio`.
  * **Link:** The FEniCS Project documentation is the primary source. An example of the structure can be seen in discussions like [FEniCS Project Discourse](https://fenicsproject.discourse.group/t/getting-information-from-xml-file-with-dof-index/13925).

### **Exodus (.e,.exo)**

  * **Description:** Exodus is a model developed at Sandia National Laboratories for storing finite element analysis data (mesh, material properties, loads, results). It is built on the NetCDF file format, making it self-describing and machine-independent.
  * **Link:** The primary source is Sandia National Laboratories. The [GitHub repository](https://github.com/sandialabs/exodusii) provides information.

### **FLAC3D (.f3grid)**

  * **Description:** This is the native grid file format for FLAC3D, a geotechnical analysis software. It can be ASCII or binary and is capable of preserving metadata like group assignments and extra variable assignments to zones, faces, and gridpoints.
  * **Link:** Documentation is available from Itasca Consulting Group Inc., for example, [Itasca FLAC3D Docs](https://docs.itascacg.com/flac3d700/flac3d/docproject/source/modeling/problemsolving/gridgeneration.html).

### **H5M (.h5m)**

  * **Description:** H5M is MOAB's (Mesh-Oriented datABase) native file format, built on HDF5. It stores mesh entities (vertices, elements, sets) and arbitrary metadata (tags) using a unique entity ID space.
  * **Link:** [MOAB H5M Documentation](https://www.mcs.anl.gov/~fathom/moab-docs/h5mmain.html).

### **Medit (.mesh,.meshb)**

  * **Description:** File format used by the Medit program (by *Pascal Frey*) to define 2D or 3D meshes (triangles, quadrilaterals, tetrahedra, hexahedra) for FEM. It can be ASCII (`.mesh`) or binary (`.meshb`) and includes sections for vertices, elements, and optional reference markers. Note: There's also a "meditmesh" format specific to Medit (the company) software, which is different.
  * **Link:** [Medit File Format Description](https://people.sc.fsu.edu/~jburkardt/data/medit/medit.html). The original specification is in Pascal Frey's technical report RT-0253, INRIA.

### **MED/Salome (.med)**

  * **Description:** MED (Modèle d'Échange de Données) is the standard file format for meshes and fields within the SALOME platform, based on HDF5. It stores meshes (nodes, various element types including polygons/polyhedra) and fields (data associated with mesh entities).
  * **Link:** Information is available at [SALOME MED GitHub](https://github.com/SalomePlatform/med) and [Spack Packages](https://packages.spack.io/package.html?name=salome-med).

### **Nastran (bulk data,.bdf,.fem,.nas)**

  * **Description:** The Nastran Bulk Data File (BDF) is a widely used ASCII input format for Nastran FEA software. It defines the model using "cards" or entries for geometry, elements, materials, loads, and constraints. It supports fixed-format and free-format fields.
  * **Link:** Documentation from MSC Software, NEi Nastran, or Autodesk Nastran are primary sources. Examples include [Autodesk Nastran Help](https://help.autodesk.com/view/NINCAD/2025/ENU/?guid=GUID-42B54ACB-FBE3-47CA-B8FE-475E7AD91A00) and [Strand7 Documentation](https://www.strand7.com/strand7r3help/Content/Topics/FileFormats/FileFormatsNASTRANFileBulkDataEntries.htm?TocPath=File%20Formats|NASTRAN%20File|_____4).

### **Netgen (.vol,.vol.gz)**

  * **Description:** Native mesh file format for Netgen, an automatic 3D tetrahedral mesh generator. `.vol` stores volumetric mesh data, including nodes, elements (typically tetrahedra), and surface patch information crucial for boundary definitions. `.vol.gz` is a gzip-compressed version.
  * **Link:** Netgen/NGSolve documentation is the primary source. Some information is available at [CalculiX Documentation](https://web.mit.edu/calculix_v2.7/CalculiX/cgx_2.7/doc/cgx/node204.html) and [FFEA Readthedocs](https://ffea.readthedocs.io/en/stable/surftovoltut.html).

### **Neuroglancer precomputed format**

  * **Description:** A directory-based format for visualizing large-scale volumetric image data (e.g., electron microscopy) in the Neuroglancer web viewer. It uses an `info` JSON file for metadata and stores image data in chunked, multi-resolution subdirectories.
  * **Link:** Specification details are on [Neuroglancer GitHub](https://github.com/google/neuroglancer/blob/master/src/neuroglancer/datasource/precomputed/volume.md#info-json-file-specification).

### **Gmsh (.msh)**

  * **Description:** Native ASCII mesh file format for Gmsh, an open-source 3D FEM mesh generator. It's organized into sections (e.g., `$MeshFormat`, `$Nodes`, `$Elements`) and supports various element types and post-processing data.
  * **Link:** [Gmsh File Formats](http://gmsh.info/doc/texinfo/gmsh.html#File-formats).

### **OBJ (.obj)**

  * **Description:** A plain text (ASCII) file format for 3D geometry, storing vertices, texture coordinates, vertex normals, and polygonal faces. Material properties are typically defined in a separate `.mtl` file.
  * **Link:** Originally by Wavefront Technologies. Widely documented, e.g., [Adobe's OBJ Page](https://www.adobe.com/products/substance3d/discover/what-are-obj-files.html).

### **OFF (.off)**

  * **Description:** Object File Format (OFF) is a simple ASCII format for 2D or 3D polygonal models, listing vertices and faces. Often used in computational geometry.
  * **Link:** [Princeton University OFF Format](http://segeval.cs.princeton.edu/public/off_format.html).

### **PERMAS (.post,.post.gz,.dato,.dato.gz)**

  * **Description:** File extensions associated with PERMAS FEA software. `.dato` (or `.dat`) likely refers to model input data files, and `.post` to post-processing/results files. These are proprietary to PERMAS.
  * **Link:** [Intes.de](http://www.intes.de).

### **PLY (.ply)**

  * **Description:** Polygon File Format (or Stanford Triangle Format) designed for 3D scanned data and general polygonal models. It has an ASCII header defining elements (vertices, faces) and their properties (coordinates, color, normals), followed by ASCII or binary data.
  * **Link:** Originated at Stanford University Graphics Lab. [Wikipedia](https://en.wikipedia.org/wiki/PLY_\(file_format\)) and [Gatech PLY Information](https://sites.cc.gatech.edu/projects/large_models/ply.html) are key resources.

### **STL (.stl)**

  * **Description:** STereoLithography (or Standard Triangle Language) format describes 3D surface geometry as a triangular mesh. It exists in ASCII and more common binary forms. Widely used for 3D printing, it does not natively support color or materials.
  * **Link:** Original specification "StereoLithography Interface Specification, 3D Systems, Inc., October 1989". Information available at [Library of Congress](https://www.loc.gov/preservation/digital/formats/fdd/fdd000505.shtml) and [Adobe STL File](https://www.adobe.com/creativecloud/file-types/image/vector/stl-file.html).

### **Tecplot (.dat)**

  * **Description:** Tecplot ASCII data file (`.dat`) for the Tecplot visualization software. It organizes data into zones and can store geometry and associated scalar/vector field data. Tecplot also supports binary (`.plt`) and Sub-Zone Load-on-demand (`.szplt`) formats.
  * **Link:** Tecplot's official website ([tecplot.com](https://www.tecplot.com/)) is the primary source. Details from [Burkardt Tecplot](https://people.math.sc.edu/Burkardt/data/tec/tec.html) and [Tecplot Data File Types](https://tecplot.com/2016/09/16/tecplot-data-file-types-dat-plt-szplt/).

### **TetGen (.node/.ele)**

  * **Description:** ASCII file formats used by TetGen for tetrahedral mesh generation. `.node` files list 3D nodal coordinates with attributes and boundary markers. `.ele` files list tetrahedral elements with node connectivity and optional region attributes.
  * **Link:** TetGen website (e.g., [WIAS-Berlin Tetgen](https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual006.html) or [Tetgen.org](https://www.google.com/search?q=https://www.tetgen.org/)).

### **SVG (2D output only) (.svg)**

  * **Description:** Scalable Vector Graphics (SVG) is an XML-based vector image format for 2D graphics. It's a W3C standard, scalable without quality loss, and widely used on the Web.
  * **Link:** The [W3C website](https://www.w3.org/Graphics/SVG/). General information from [Adobe SVG File](https://www.adobe.com/creativecloud/file-types/image/vector/svg-file.html).

### **SU2 (.su2)**

  * **Description:** Native ASCII mesh file format for the SU2 open-source CFD suite. It defines dimensionality, node coordinates, element connectivity (using VTK element type IDs), and boundary markers (tags).
  * **Link:** [SU2 GitHub Documentation](https://su2code.github.io/docs/Mesh-File/).

### **UGRID (.ugrid)**

  * **Description:** Stores 3D grid data, including boundary surface grids (triangles/quads) and optional polyhedral volume grids. Can be FORTRAN unformatted, C binary, or ASCII. Includes header, node coordinates, face/element connectivity, and optional records for IDs and boundary information.
  * **Link:** Documentation from SimCenter at Mississippi State University: [UGRID File Type](https://www.simcenter.msstate.edu/software/documentation/ug_io/3d_grid_file_type_ugrid.html).

### **VTK (.vtk)**

  * **Description:** The simple legacy VTK format for the Visualization Toolkit. It's ASCII or binary and consists of a version identifier, header, format type, dataset structure (e.g., `STRUCTURED_POINTS`, `UNSTRUCTURED_GRID`, `POLYDATA`), and dataset attributes (scalars, vectors, etc.).
  * **Link:** [VTK File Formats](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html).

### **VTU (.vtu)**

  * **Description:** XML-based VTK format for `vtkUnstructuredGrid` data. It defines points, cells (connectivity, offsets, types), and associated point/cell data using `<DataArray>` elements. Supports ASCII, base64 binary, or appended binary data, with optional compression.
  * **Link:** [VTK XML File Formats](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#xml-file-formats).

### **WKT (TIN) (.wkt)**

  * **Description:** Well-Known Text is an OGC and ISO standard markup language for vector geometry objects and spatial reference systems. TIN (Triangulated Irregular Network) is one such geometry, representing a surface as a collection of non-overlapping triangles.
  * **Link:** OGC and ISO specifications (e.g., ISO/IEC 13249-3:2016). General information at [Wikipedia Well-Known Text](https://en.wikipedia.org/wiki/Well-known_text_representation_of_geometry) and [LibGEOS WKT](https://libgeos.org/specifications/wkt/).

### **XDMF (.xdmf,.xmf)**

  * **Description:** The eXtensible Data Model and Format separates metadata (in XML) from bulk numerical data (often in HDF5 or binary files). The XML describes the data model, structure, and pointers to the heavy data, facilitating exchange between HPC codes and visualization tools.
  * **Link:** [XDMF Model and Format](http://www.xdmf.org/index.php/XDMF_Model_and_Format).

## Example

Let's say we want to import an VTK mesh.

~~~py
import meshio
mesh = meshio.read("test.vtk")
meshio.write("test.mdpa", mesh)
~~~

Additionally, due to the support of *meshio* it is possible to create a mesh using [pygmsh](https://github.com/nschloe/pygmsh). Check the link for examples.
