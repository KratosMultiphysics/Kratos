## Mapping Application

The Mapping Application contains the core developments in mapping data between non matching grids. It works both in shared and distributed (**MPI**) memory environments as well as in 1D, 2D and 3D domains.

### Overview

- [List of features](#list-of-features)
- [Dependencies](#dependencies)
- [Mapping in CoSimulation](#Mapping-in-CoSimulation)
- [Basic Usage](#basic-usage)
- [Advanced Usage](#advanced-usage)
- [Available Mappers](#available-mappers)
- [When to use which Mapper?](#when-to-use-which-mapper)
- [Using the Mapper for ModelParts that are not part of all ranks](#using-the-mapper-for-modelparts-that-are-not-part-of-all-ranks)
- [Miscellaneous functionalities](#miscellaneous-functionalities)
- [FAQ](#faq)

### List of features

- Parallelism:
  - Serial (no parallelism)
  - Shared memory (OpenMP)
  - Distributed memory (MPI)
- Domain sizes: 1D / 2D / 3D
- Matching and non matching grids
- Different mapping technologies (see [here](#available-mappers)):
  - Nearest Neighbor
  - Nearest Element
  - Barycentric
  - Radial Basis Function Mapper
  - Beam Mapper
- Metamappers
  - 3D/2D metamapper (metamapper which obtains the solution for the 3D destination model part from the original 2D solution)
- Mapping operations (see [here](#customizing-the-behavior-of-the-mapping-with-flags))

### Dependencies

The serial / shared memory parallel compilation of the Mapping Application doesn't have any dependencies (except the `KratosCore`).

The distributed compilation of the Mapping Application depends on the [Trilinos library](https://trilinos.github.io/). Also most of the MPI-solvers in Kratos depend on Trilinos, see the [Trilinos Application](../TrilinosApplication).

### Mapping in CoSimulation

The Mapping Application can be used for mapping within the [CoSimulation Application](../CoSimulationApplication). This can be done by using the  [KratosMappingDataTransferOperator](../CoSimulationApplication/python_scripts/data_transfer_operators/kratos_mapping.py).

### Basic Usage

The _Mapper_ maps nodal data from one _ModelPart_ to another. This means that the input for the _Mapper_ is two _ModelParts_, the **Origin** and the **Destination**. Furthermore settings in the form of _Kratos::Parameters_ are passed.

The _Mapper_ is constructed using the _MapperFactory_. See the following basic example.

```py
# import the Kratos Core
import KratosMultiphysics as KM
# import the MappingApplication to load the mappers
import KratosMultiphysics.MappingApplication as KratosMapping

# create ModelParts
# ...

mapper_settings = KM.Parameters("""{
    "mapper_type": "nearest_neighbor",
    "echo_level" : 0
}""")

# creating a mapper for shared memory
mapper = KM.MapperFactory.CreateMapper(
    model_part_origin,
    model_part_destination,
    mapper_settings
)
```

For constructing an _MPI-Mapper_ use the `MPIExtension` instead:

```py
# creating a mapper for distributed memory
from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension
mpi_mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(
    model_part_origin,
    model_part_destination,
    mapper_settings
)
```

After constructing the _Mapper_ / _MPI-Mapper_ it can be used immediately to map any scalar and vector quantities, no further initialization is necessary.\
The **Map** function is used to map values from the **Origin** to the **Destination**. For this the _Variables_ have to be specified. See the following example for mapping scalar quantities.

```py
# mapping scalar quantities
# this maps the nodal quantities of TEMPERATURE on the origin-ModelPart
# to the nodal quantities of AMBIENT_TEMPERATURE on the destination-ModelPart

mapper.Map(KM.TEMPERATURE, KM.AMBIENT_TEMPERATURE)
```

The **Map** function is overloaded, this means that mapping vector quantities works in the same way as mapping scalar quantities.

```py
# mapping vector quantities
# this maps the nodal quantities of VELOCITY on the origin-ModelPart
# to the nodal quantities of MESH_VELOCITY on the destination-ModelPart.

mapper.Map(KM.VELOCITY, KM.MESH_VELOCITY)
```

Mapping from **Destination** to **Origin** can be done using the **InverseMap** function which works in the same way as the **Map** function.

```py
# inverse mapping scalar quantities
# this maps the nodal quantities of AMBIENT_TEMPERATURE on the destination-ModelPart
# to the nodal quantities of TEMPERATURE on the origin-ModelPart

mapper.InverseMap(KM.TEMPERATURE, KM.AMBIENT_TEMPERATURE)

# inverse mapping vector quantities
# this maps the nodal quantities of MESH_VELOCITY on the destination-ModelPart
# to the nodal quantities of VELOCITY on the origin-ModelPart

mapper.InverseMap(KM.VELOCITY, KM.MESH_VELOCITY)
```

For the 3D/2D metamapper the settings to consider are the following, where `base_mapper` is the backend mapper to be considered.

```json
mapper_params = KM.Parameters("""{
    "mapper_type"     : "projection_3D_2D",
    "base_mapper"     : "nearest_neighbor",
    "search_settings" : {},
    "echo_level"      : 0
}""")
```

### Advanced Usage

The previous section introduced the basics of using the _MappingApplication_. The more advanced usage is explained in this section.

#### Customizing the behavior of the mapping with Flags

By default the mapping functions **Map** and **InverseMap** will overwrite the values where they map to. In order to add instead of overwrite the values the behavior can be customized by using _Kratos::Flags_. Consider in the following example that several forces are acting on a surface. Overwriting the values would cancel the previously applied forces.

```py
# Instead of overwriting, this will add the values to the existing ones

mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.ADD_VALUES)
```

Sometimes it can be necessary to swap the signs of quantities that are to be mapped. This can be done with the following:

```py
# Swapping the sign, i.e. multiplying the values with (-1)

mapper.Map(KM.DISPLACEMENT, KM.MESH_DISPLACEMENT, KM.Mapper.SWAP_SIGN)
```

The flags can also be combined:

```py
mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.ADD_VALUES | KM.Mapper.SWAP_SIGN)
```

Historical nodal values are used by default. Mapping to an from nonhistorical nodal values is also supported, the following examples show the usage:

This maps the values from the origin (`REACTION`) as historical values to the destination (`FORCE`) as nonhistorical values:

```py
mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.TO_NON_HISTORICAL)
```

This maps the values from the origin (`REACTION`) as nonhistorical values to the destination (`FORCE`) as historical values:

```py
mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.FROM_NON_HISTORICAL)
```

This maps the values from the destination (`FORCE`) as historical values to the origin (`REACTION`) as nonhistorical values:

```py
mapper.InverseMap(KM.REACTION, KM.FORCE, KM.Mapper.TO_NON_HISTORICAL)
```

This maps the values from the destination (`FORCE`) as nonhistorical values to the origin (`REACTION`) as historical values:

```py
mapper.InverseMap(KM.REACTION, KM.FORCE, KM.Mapper.FROM_NON_HISTORICAL)
```

Of course it is possible to use both origin and destination nonhistorical. This maps the values from the origin (`REACTION`) as nonhistorical values to the destination (`FORCE`) as nonhistorical values:

```py
mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.FROM_NON_HISTORICAL | KM.Mapper.TO_NON_HISTORICAL)
```

Many _Mappers_ internally construct a mapping matrix. It is possible to use the transpose of this matrix for mapping with `USE_TRANSPOSE`. This is often used for conservative mapping of forces in FSI, when the virtual work on both interfaces should be preserved.

```py
mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.USE_TRANSPOSE)
```

#### Updating the Interface

In case of moving interfaces (e.g. in a problem involving Contact between bodies) it can become necessary to update the _Mapper_ to take the new geometrical positions of the interfaces into account.\
One way of doing this would be to construct a new _Mapper_, but this is not efficient and sometimes not even possible.

Hence the _Mapper_ provides the **UpdateInterface** function for updating itseld with respect to the new geometrical positions of the interfaces.\
Note that this is potentially an expensive operation due to searching the new geometrical neighbors on the interface.

```py
mapper.UpdateInterface()
```

#### Checking which mappers are available

The following can be used to see which _Mappers_ are available:

```py
# available mappers for shared memory
KM.MapperFactory.GetRegisteredMapperNames()

# available mappers for distributed memory
MappingMPIExtension.MPIMapperFactory.GetRegisteredMapperNames()

# check if mapper for shared memory exists
KM.MapperFactory.HasMapper("mapper_name")

# check if mapper for distributed memory exists
MappingMPIExtension.MPIMapperFactory.HasMapper("mapper_name")
```

#### Search settings
The search of neighbors / partners on the other side of the interface is a crucial task when creating the mapper. Especially in distributed computations (MPI) this can be very expensive and time consuming. Hence the search of the mapper is very optimized to provide robust and fast results. For this the search works in several iterations where the search radius is increased in each iteration.
The default settings of the search are working fine in most cases, but in some special cases it might still be necessary to tweak and optimize the settings. The following settings are available (as sub-parameter `search_settings` of the settings that are given to the mapper):

| name | type | default| description |
|---|---|---|---|
| `search_radius`| `double` | computed | The search radius to start with in the first iteration. In each next iteration it will be increased by multiplying with `search_radius_increase_factor` (`search_radius *= search_radius_increase_factor`) |
| `max_search_radius` | `double` | computed | The max search radius to use. |
| `search_radius_increase_factor`| `double` | `2.0` | factor by which the search radius is increasing in each search iteration (see above). **Tuning this parameter is usually the best way to achieve a faster search**. In many cases decreasing it will speed up the search, especially for volumetric mapping, but it is case dependent. |
| `max_num_search_iterations` | `int` | computed (min 3) | max number of search iterations that is conducted. If the search is successful before then it will terminate earlier. The more heterogeneous the mesh the larger this will be.

It is recommended to set the `echo_level` to 2 or higher for getting useful information from the search. This will help to debug the search in case of problems.

### Available Mappers

This section explains the theory behind the mappers.

#### Nearest Neighbor

The _NearestNeighborMapper_ is a very simple/basic _Mapper_. Searches its closest neighbor (node) on the other interface. During mapping it gets/sets its value to the value of its closest neighbor.

This mapper is best suited for problems where both interfaces have a similar discretization. Furthermore it is very robust and can be used for setting up problems when one does not (yet) want to deal with mapping.

Internally it constructs the mapping matrix, hence it offers the usage of the transposed mapping matrix. When using this, for very inhomogeneous interface discretizations it can come to oscillations in the mapped quantities.

**Supported mesh topologies**: This mapper only works with nodes and hence supports any mesh topology

#### Nearest Neighbor for IGA scenarios

The _NearestNeighborMapperIGA_ is a simple and robust Mapper for IGA/FEM partitioned simulations. For each node on the FEM side, it finds the closest integration point on the IGA interface.
During mapping, it evaluates the IGA shape functions at that location to assemble the mapping matrix.

This mapper is suited for cases where the origin domain is discretized with IGA elements and the destination with any node-based discretization technique (e.g., FEM or FCV).

Internally, it constructs the mapping matrix, which also allows the use of its transpose for conservative mapping (e.g., mapping forces from FEM to IGA). In cases of highly inhomogeneous interface discretizations, using the transpose may introduce oscillations in the mapped values.

**Supported mesh topologies**: This mapper operates on nodes and supports any mesh topology. The only requirement is that the origin domain must be the IGA domain; otherwise, the mapping problem is not well defined.

#### Nearest Element

The _NearestElementMapper_ projects nodes to the elements( or conditions) on other side of the interface. Mapping is then done by interpolating the values of the nodes of the elements by using the shape functions at the projected position. The NearestElementMapper supports IGA/FEM partitioned simulations where the origin must be the IGA domain. Each FEM node is projected onto the IGA surface or its boundary curves, and the shape functions are evaluated at that point to assemble the mapping matrix.

This mapper is best suited for problems where the _NearestNeighborMapper_ cannot be used, i.e. for cases where the discretization on the interfaces is different. Note that it is less robust than the _NearestNeighborMapper_ due to the projections it performs. In case a projection fails, it uses an approximation that is similar to the approach of the _NearestNeighborMapper_. This can be disabled by setting `use_approximation` to `false` in the mapper-settings.

Internally it constructs the mapping matrix, hence it offers the usage of the transposed mapping matrix. When using this, for very inhomogeneous interface discretizations it can come to oscillations in the mapped quantities.

**Supported mesh topologies**: Any mesh topology available in Kratos, which includes the most common linear and quadratic geometries, see [here](../../kratos/geometries).

#### Barycentric

The _BarycentricMapper_ uses the closest nodes to reconstructs a geometry. This geometry is used in the same way as the _NearestElementMapper_ for interpolating the values of the nodes using the shape functions.

This mapper can be used when no geometries are available and interpolative properties of the mapper are required. E.g. for particle methods when only nodes or point-based entities are available. Overall it can be seen as combining the advantages of the _NearestNeighborMapper_ (which only requires points as input) with the advantages of the _NearestElementMapper_ (which has interpolative properties). The disadvantage is that the reconstruction of the geometry can cause problems in complex situations, hence it should only be used if the _NearestElementMapper_ cannot be used.

Furthermore, the geometry type for the reconstruction/interpolation has to be chosen with the `interpolation_type` setting. The following types are available: `line`, `triangle` and `tetrahedra`

Internally it constructs the mapping matrix, hence it offers the usage of the transposed mapping matrix. When using this, for very inhomogeneous interface discretizations it can come to oscillations in the mapped quantities.

**Supported mesh topologies**: This mapper only works with nodes and hence supports any mesh topology

#### Radial Basis Function (RBF) Mapper 
The _RadialBasisFunctionMapper_ is a global, mesh-independent mapper that constructs a smooth interpolation field based on Radial Basis Functions (RBFs). In contrast to purely local methods, this mapper uses all (or a user-defined subset of) points from the origin interface to build an RBF system that is then evaluated at the destination points.

This allows for smooth, high-quality transfer of field quantities between arbitrarily discretized, non-matching, or strongly non-uniform interfaces. It is therefore particularly suitable for multi-physics problems where interface meshes can differ substantially.

The default configuration looks as follows:

```json
"mapper_settings" : {
    "echo_level"                     : 0,
    "radial_basis_function_type"     : "thin_plate_spline",
    "additional_polynomial_degree"   : 0,
    "origin_is_iga"                  : false,
    "destination_is_iga"             : false,
    "max_support_points"  : 0,
    "use_all_rbf_support_points": true,
    "precompute_mapping_matrix"      : true,
    "search_settings"                : {},
    "linear_solver_settings"         : {}
}
```

**Important notes:**
- Only the origin domain can be IGA.
The mapper can extract coordinates from IGA gauss points on the origin side (```json"origin_is_iga": true```). The destination side must currently be a standard finite element mesh (nodes as mapping coordinates).
-Global RBF System
Internally, the mapper assembles and solves a global RBF interpolation system.
If ```json"precompute_mapping_matrix": true```, the resulting mapping matrix is stored and reused, allowing efficient repeated mapping calls.
- Support Points
  - ```json"use_all_rbf_support_points": true``` → all origin points contribute to the RBF system.
  - ```json"max_support_points" > 0``` → restricts support to a local neighborhood for each destination point.


#### Beam Mapper 
The _BeamMapper_ provides support for mapping between 1D beam elements and 2D/3D surface meshes. It follows the formulation of Wang (2019) and is intended for cases where beam DOFs (displacements and rotations) must be transferred consistently to a surrounding surface, e.g. in FSI or beam–solid coupling.

The mapper projects each surface node onto the undeformed beam centerline and assigns a local rigid cross section. The motion of every projected surface point is obtained through rigid body motion of this cross section, using the beam’s axial and rotational DOFs. Hermitian shape functions are used along the beam to interpolate both displacements and rotations.

A typical configuration json looks as follows:

```json
"mapper" : {
    "type" : "kratos_beam_mapping",
    "model_part_name_beam"         : "Structure.Parts_Beam_beam",
    "model_part_name_surface"      : "Structure.Parts_Shell_wet_surface",
    "solver_name_beam": "beam_structure",
    "solver_name_surface" : "dummy_fluid",
    "echo_level": 3,
    "mapper_settings" : {
        "mapper_type" : "beam_mapper",
        "use_corotation" : false,
        "search_settings": {
            "max_num_search_iterations"     : 30,
            "search_radius": 3.0
        },
        "echo_level": 0
    }
}
```

**Explanation of the main entries:**
- ```json"type"```: selects the CoSimulation mapper. "kratos_beam_mapping" activates the _BeamMapper_.
- ```json"model_part_name_beam" / "model_part_name_surface"```: names of the origin beam model part and the target surface model part used for the mapping.
- ```json"solver_name_beam" / "solver_name_surface"```: identifiers of the solvers that own these model parts. They are used internally by the CoSimulation framework to retrieve nodal values.
- ```json"echo_level"```: controls the amount of printed information for debugging.
- ```json"mapper_settings"```
  - ```json"mapper_type"```: must be "beam_mapper" to use the _BeamMapper_.
  - ```json"use_corotation"```: enables the co-rotational formulation (false → linear mapping, true → large-rotation mapping)
  - ```json"search_settings"```: parameters controlling the projection of surface nodes onto the beam centerline.

### When to use which Mapper?

- **Matching Interface**\
  For a matching interface the _NearestNeighborMapper_ is the best / fastes choice. Note that the ordering / numbering of the nodes doesn't matter.

- **Interfaces with almost matching discretizations**\
  In this case both the _NearestNeighborMapper_ and the _NearestElementMapper_ can yield good results.

- **Interfaces with non matching discretizations**\
  The _NearestElementMapper_ is recommended because it results in smoother mapping results due to the interpolation using the shape functions.

- **Interfaces with non matching discretizations when no geometries are available for interpolation**\
  The _NearestElementMapper_ cannot be used as it requires geometries for the ionterpolation. Here the _BarycentricMapper_ is recommended because it reconstructs geometries from the surrounding nodes and then uses it to interpolate.

### Using the Mapper for ModelParts that are not part of all ranks

In MPI parallel simulations usually all `ModelParts` are distributed across all ranks. However in some cases this does not hold, for example in FSI when the fluid runs on all ranks but the structure runs serial on one rank. In this case it is necessary to do the following:

- Create a dummy-`ModelPart` on the ranks that do not have the original ModelPart.
- **IMPORTANT**: This `ModelPart` must have a `DataCommunicator` that is not defined on the ranks that are not part of the original `ModelPart`.
- Create and MPI-mapper as explained [above](#basic-usage), using the original and the dummy `ModelPart`s on the respective ranks.

Check [this test](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MappingApplication/tests/blade_mapping_test.py) for more details and usage example.

For an example the following assumptions are made:

- Overall 4 MPI processes are used
- `model_part_fluid` is distributed across all 4 ranks
- `model_part_structure` is not distributed and exists only on rank 0

```py
import KratosMultiphysics as KM
import KratosMultiphysics.mpi as KratosMPI

# "model_part_fluid" was already read and exists on all ranks
# "model_part_structure" was already read and exists only on rank 0


# getting the DataCommunicator that wraps `MPI_COMM_WORLD` i.e. contains all ranks
world_data_comm = KM.ParallelEnvironment.GetDataCommunicator("World)

# define the ranks on which the structure ModelPart exists
# structure can also be distributed across several (but not all) ranks
structure_ranks = [0]

# create a DataCommunicator containing only the structure ranks
structure_ranks_data_comm_name = "structure_ranks"
data_comm_all_structure_ranks = KratosMPI.DataCommunicatorFactory.CreateFromRanksAndRegister(
    world_data_comm,
    structure_ranks,
    structure_ranks_data_comm_name)

# create a dummy ModelPart on the ranks where the original ModelPart does not exist
if world_data_comm.Rank() not in structure_ranks:
    dummy_model = KM.Model()
    model_part_structure = dummy_model.CreateModelPart("structure_dummy")

    # Important: set the DataCommunicator so that the Mapper knows on which ranks the ModelPart is only a dummy
    KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(model_part_structure, data_comm_all_structure_ranks)

# now the Mapper can be created with the original and the dummy ModelParts
mpi_mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(
    model_part_fluid,
    model_part_structure,
    mapper_settings
)
```

### Miscellaneous functionalities
- [serial_output_process](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/MappingApplication/python_scripts/serial_output_process.py): This process can be used to map results to one rank and then do postprocessing on this rank. This has two advantages:
  - Some output formats write one file per rank in distributed simulations, which leads to many files when running with many cores. This process collects the results on one rank and can hence reduce the number of files significantly
  - Different meshes can be used to do the postprocessing. This is in particular useful when the computational mesh is very fine, but a coarser mesh would be sufficient for postprocessing.

  <ins>The following input parameters are used:</ins>
  - `model_part_name_origin`: name of the origin ModelPart where the data comes from (is being mapped from)
  - `model_part_name_destination`: name of destination ModelPart where the data is mapped to. This ModelPart is being read.
  - `mdpa_file_name_destination`: name of the mdpa file containing the mesh that is used for the destination
  - `historical_variables_destination` list of historical variables that are allocated on the destination ModelPart
  - `destination_rank` rank on which the processing of the destination happens (i.e. the rank on which the destination ModelPart is read). Note that this increases the memory usage significantly, especially for large destination meshes. The default is rank 0, which in most distributed simulations acts as the master rank with already increased computational effort. Hence it can make sense to use another rank, preferably on another compute node, to optimize the memory and computational load balance
  - `mapper_settings`: setting that are passed to the mapper, as explained above
  - `mapping_settings`: list of mapping steps to be executed before the postprocessing is done. `variable_origin` and `variable_destination` must be specified, while `mapping_options` is optional and can contain the flags as explained above.
  - `output_process_settings`: The settings for the output process (which will be only executed on the destination rank). **Important**: For mapping onto a serial ModelPart, the DataCommunicator is set as explained [here](#using-the-mapper-for-modelparts-that-are-not-part-of-all-ranks). This means that the destination ModelPart is not valid on other ranks and can hence not be used in the regular postprocessing (which happens also on the ranks where it is not valid and hence some MPI-functionalities would fail)
  Example input:
  ~~~js
  "python_module" : "serial_output_process",
  "kratos_module" : "KratosMultiphysics.MappingApplication",
  "Parameters"    : {
      "model_part_name_origin"      : "FluidModelPart",
      "model_part_name_destination" : "PostProcessing",
      "mdpa_file_name_destination"  : "coarse_mesh",
      "historical_variables_destination" : ["REACTION", "DISPLACEMENT"],
      "mapper_settings" :  {"mapper_type" : "nearest_neighbor"},
      "mapping_settings" : [{
          "variable_origin" : "REACTION",
          "variable_destination" : "REACTION"
      },{
          "variable_origin" : "REACTION",
          "variable_destination" : "REACTION",
          "mapping_options" : ["add_values"]
      },{
          "variable_origin" : "MESH_DISPLACEMENT",
          "variable_destination" : "DISPLACEMENT"
      }],
      "output_process_settings" : {
          "python_module" : "vtk_output_process",
          "kratos_module" : "KratosMultiphysics",
          "Parameters"    : {
              // ...
          }
      }
  }
  ~~~

### FAQ

- **Is mapping of elemental / conditional data or gauss-point values possible?**\
  The mapper only supports mapping of nodal data. In order to map other quantities, those have to first be inter- / extrapolated to the nodes.

- **Something is not working with the mapping. What should I do?**\
  Problems with mapping can have many sources. The first thing in debugging what is happening is to increase the `echo_level` of the _Mapper_. Then in many times warnings are shown in case of some problems.

- **I get oscillatory solutions when mapping with `USE_TRANSPOSE`**\
  Research has shown that "simple" mappers like _NearestNeighbor_ and _NearestElement_ can have problems with mapping with the transpose (i.e. when using `USE_TRANSPOSE`) if the meshes are very different. Using the _MortarMapper_ technology can improve this situation. This _Mapper_ is currently under development.

- **Projections find the wrong result**\
  For complex geometries the projections can fail to find the correct result if many lines or surfaces are close. In those situations it helps to partition the mapping interface and construct multiple mappers with the smaller interfaces.

- **Creation of the mapper takes very long**\
  Often this is because of of unfit search settings. If the settings are not suitable for the problem then the mapper creation time can increase several magnitudes! Check [here](#search-settings) for an explanation of how to set the search settings in case the defaults are not working well.
