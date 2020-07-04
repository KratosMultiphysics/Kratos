## Mapping Application

The Mapping Application contains the core developments in mapping data between non matching grids. It works both in shared and distributed (**MPI**) memory environments as well as in 1D, 2D and 3D domains.

### Basic Usage
The _Mapper_ maps nodal data from one _ModelPart_ to another. This means that the input for the _Mapper_ is two _ModelParts_, the **Origin** and the **Destination**. Furthermore settings in the form of _Kratos::Parameters_ are passed.

The _Mapper_ is constructed using the _MapperFactory_. See the following basic example.

~~~py
# import the Kratos Core
import KratosMultiphysics as KM
# import the MappingApplication
import KratosMultiphysics.MappingApplication as
KratosMapping

# create ModelParts
# ...

mapper_settings = KM.Parameters("""{
    "mapper_type": "nearest_neighbor",
    "echo_level" : 0
}""")

mapper = KratosMapping.MapperFactory.CreateMapper(
    model_part_origin,
    model_part_destination,
    mapper_settings
)
~~~

After constructing the _Mapper_ it can be used to map any scalar and vector quantities. The **Map** function is used to map values from the **Origin** to the **Destination**. For this the _Variables_ have to be specified. See the following example for mapping scalar quantities.

**Note**: In order to demonstrate that the _Mapper_ is not tied to any particular physics, arbitrary _Variables_ are chosen in the following examples.

~~~py
# mapping scalar quantities
# this maps the nodal quantities of TEMPERATURE on the origin-ModelPart
# to the nodal quantities of PRESSURE on the destination-ModelPart

mapper.Map(KM.TEMPERATURE, KM.PRESSURE)
~~~

The **Map** function is overloaded, this means that mapping vector quantities works in the same way as mapping scalar quantites.

~~~py
# mapping vector quantities
# this maps the nodal quantities of VELOCITY on the origin-ModelPart
# to the nodal quantities of DISPLACEMENT on the destination-ModelPart.

mapper.Map(KM.VELOCITY, KM.DISPLACEMENT)
~~~

### Advanced Usage


### When to use which Mapper?

### FAQ

- **Is mapping of elemental / conditional data or gauss-point values possible?**\
  The mapper only supports mapping of nodal data. In order to map other quantities, those have to first be inter- / extrapolated to the nodes.

- **Something is not working with the mapping. What should I do?**\
  Problems with mapping can have many sources. The first thing in debugging what is happening is to increase the `echo_level` of the _Mapper_. Then in many times warnings are shown in case of some problems.


### Available Mappers
This section explains the theory behind the different mappers.

#### Nearest Neighbor

#### Nearest Element