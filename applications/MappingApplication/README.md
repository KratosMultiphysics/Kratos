## Mapping Application
The Mapping Application contains the core developments in mapping data between non matching grids. It works both in shared and distributed (**MPI**) memory environments as well as in 1D, 2D and 3D domains.


### Overview
- [List of features](#list-of-features)
- [Dependencies](#dependencies)
- [Mapping in CoSimulation](#Mapping-in-CoSimulation)
- [Basic Usage](#basic-usage)
- [Advanced Usage](#advanced-usage)
- [Available Mappers](#Available-Mappers)
- [When to use which Mapper?](#When-to-use-which-Mapper)
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
- Mapping operations (see [here](#customizing-the-behavior-of-the-mapping-with-flags))

### Dependencies
The serial / shared memory parallel compilation of the Mapping Application doesn't have any dependencies (except the `KratosCore`).

The distributed compilation of the Mapping Application depends on the [Trilinos library](https://trilinos.github.io/). Also most of the MPI-solvers in Kratos depend on Trilinos, see the [Trilinos Application](../TrilinosApplication).

### Mapping in CoSimulation
The Mapping Application can be used for mapping within the [CoSimulation Application](../CoSimulationApplication). This can be done by using the  [KratosMappingDataTransferOperator](../CoSimulationApplication/python_scripts/data_transfer_operators/kratos_mapping.py).

### Basic Usage
The _Mapper_ maps nodal data from one _ModelPart_ to another. This means that the input for the _Mapper_ is two _ModelParts_, the **Origin** and the **Destination**. Furthermore settings in the form of _Kratos::Parameters_ are passed.

The _Mapper_ is constructed using the _MapperFactory_. See the following basic example.

~~~py
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
~~~

For constructing an _MPI-Mapper_ use the `MPIExtension` instead:

~~~py
# creating a mapper for distributed memory
from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension
mpi_mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(
    model_part_origin,
    model_part_destination,
    mapper_settings
)
~~~

After constructing the _Mapper_ / _MPI-Mapper_ it can be used immediately to map any scalar and vector quantities, no further initialization is necessary.\
The **Map** function is used to map values from the **Origin** to the **Destination**. For this the _Variables_ have to be specified. See the following example for mapping scalar quantities.

~~~py
# mapping scalar quantities
# this maps the nodal quantities of TEMPERATURE on the origin-ModelPart
# to the nodal quantities of AMBIENT_TEMPERATURE on the destination-ModelPart

mapper.Map(KM.TEMPERATURE, KM.AMBIENT_TEMPERATURE)
~~~

The **Map** function is overloaded, this means that mapping vector quantities works in the same way as mapping scalar quantites.

~~~py
# mapping vector quantities
# this maps the nodal quantities of VELOCITY on the origin-ModelPart
# to the nodal quantities of MESH_VELOCITY on the destination-ModelPart.

mapper.Map(KM.VELOCITY, KM.MESH_VELOCITY)
~~~

Mapping from **Destination** to **Origin** can be done using the **InverseMap** function which works in the same way as the **Map** function.

~~~py
# inverse mapping scalar quantities
# this maps the nodal quantities of AMBIENT_TEMPERATURE on the destination-ModelPart
# to the nodal quantities of TEMPERATURE on the origin-ModelPart

mapper.InverseMap(KM.TEMPERATURE, KM.AMBIENT_TEMPERATURE)

# inverse mapping vector quantities
# this maps the nodal quantities of MESH_VELOCITY on the destination-ModelPart
# to the nodal quantities of VELOCITY on the origin-ModelPart

mapper.InverseMap(KM.VELOCITY, KM.MESH_VELOCITY)
~~~


### Advanced Usage
The previous section introduced the basics of using the _MappingApplication_. The more advanced usage is explained in this section.

#### Customizing the behavior of the mapping with Flags
By default the mapping functions **Map** and **InverseMap** will overwrite the values where they map to. In order to add instead of overwrite the values the behavior can be customized by using _Kratos::Flags_. Consider in the following example that several forces are acting on a surface. Overwritting the values would cancel the previously applied forces.
~~~py
# Instead of overwriting, this will add the values to the existing ones

mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.ADD_VALUES)
~~~

Sometimes it can be necessary to swap the signs of quantites that are to be mapped. This can be done with the following:

~~~py
# Swapping the sign, i.e. multiplying the values with (-1)

mapper.Map(KM.DISPLACEMENT, KM.MESH_DISPLACEMENT, KM.Mapper.SWAP_SIGN)
~~~

The flags can also be combined:
~~~py
mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.ADD_VALUES | KM.Mapper.SWAP_SIGN)
~~~

Many _Mappers_ internally construct a mapping matrix. It is possible to use the transpose of this matrix for mapping with `USE_TRANSPOSE`. This is often used for conservative mapping of forces in FSI, when the virtual work on both interfaces should be preserved.

~~~py
mapper.Map(KM.REACTION, KM.FORCE, KM.Mapper.USE_TRANSPOSE)
~~~

#### Updating the Interface
In case of moving interfaces (e.g. in a problem involving Contact between bodies) it can become necessary to update the _Mapper_ to take the new geometrical positions of the interfaces into account.\
One way of doing this would be to construct a new _Mapper_, but this is not efficient and sometimes not even possible.

Hence the _Mapper_ provides the **UpdateInterface** function for updating itseld with respect to the new geometrical positions of the interfaces.\
Note that this is potentially an expensive operation due to searching the new geometrical neighbors on the interface.

~~~py
mapper.UpdateInterface()
~~~

#### Checking which mappers are available
The following can be used to see which _Mappers_ are available:
~~~py
# available mappers for shared memory
KM.MapperFactory.GetRegisteredMapperNames()

# available mappers for distributed memory
MappingMPIExtension.MPIMapperFactory.GetRegisteredMapperNames()

# check if mapper for shared memory exists
KM.MapperFactory.HasMapper("mapper_name")

# check if mapper for distributed memory exists
MappingMPIExtension.MPIMapperFactory.HasMapper("mapper_name")
~~~


### Available Mappers
This section explains the theory behind the mappers.

#### Nearest Neighbor
The _NearestNeighborMapper_ is a very simple/basic _Mapper_. Searches its closest neighbor (node) on the other interface. During mapping it gets/sets its value to the value of its closest neighbor.

This mapper is best suited for problems where both interfaces have a similar discretization. Furthermore it is very robust and can be used for setting up problems when one does not (yet) want to deal with mapping.

Internally it constructs the mapping matrix, hence it offers the usage of the transposed mapping matrix. When using this, for very inhomogenous interface discretizations it can come to oscillations in the mapped quantities.

**Supported mesh topologies**: This mapper only works with nodes and hence supports any mesh topology

#### Nearest Element
The _NearestElementMapper_ projects nodes to the elements( or conditions) on other side of the inteface. Mapping is then done by interpolating the values of the nodes of the elements by using the shape functions at the projected position.

This mapper is best suited for problems where the _NearestNeighborMapper_ cannot be used, i.e. for cases where the discretization on the interfaces is different. Note that it is less robust than the _NearestNeighborMapper_ due to the projections it performs. In case a projection fails it uses an approximation that is similar to the approach of the _NearestNeighborMapper_.

Internally it constructs the mapping matrix, hence it offers the usage of the transposed mapping matrix. When using this, for very inhomogenous interface discretizations it can come to oscillations in the mapped quantities.

**Supported mesh topologies**: Any mesh topology available in Kratos, which includes the most common linear and quadratic geometries, see [here](../../kratos/geometries).


### When to use which Mapper?
- **Matching Interface**\
  For a matching interface the _NearestNeighborMapper_ is the best / fastes choice. Note that the ordering / numbering of the nodes doesn't matter.

- **Interfaces with almost matching discretizations**\
  In this case both the _NearestNeighborMapper_ and the _NearestElementMapper_ can yield good results.

- **Interfaces with non matching discretizations**\
  The _NearestElementMapper_ is recommended because it results in smoother mapping results due to the interpolation using the shape functions.


### FAQ

- **Is mapping of elemental / conditional data or gauss-point values possible?**\
  The mapper only supports mapping of nodal data. In order to map other quantities, those have to first be inter- / extrapolated to the nodes.

- **Something is not working with the mapping. What should I do?**\
  Problems with mapping can have many sources. The first thing in debugging what is happening is to increase the `echo_level` of the _Mapper_. Then in many times warnings are shown in case of some problems.

- **I get oscillatory solutions when mapping with `USE_TRANSPOSE`**\
  Research has shown that "simple" mappers like _NearestNeighbor_ and _NearestElement_ can have problems with mapping with the transpose (i.e. when using `USE_TRANSPOSE`) if the meshes are very different. Using the _MortarMapper_ technology can improve this situation. This _Mapper_ is currently under development.

- **Projections find the wrong result**\
  For complex geometries the projections can fail to find the correct result if many lines or surfaces are close. In those situations it helps to partition the mapping interface and construct multiple mappers with the smaller interfaces.
