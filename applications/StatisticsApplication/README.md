# Statistics Application

Statistics application consist of widely used methods to calculate statistics in various containers of KratosMultiphysics. There are mainly two groups of statistical methods namely **Spatial** and **Temporal**. **Spatial** methods calculate statistics on spatial containers and output the values whenever they are called. **Temporal** methods calculate statistics on the fly in a transient simulation. All the temporal methods gurantee that, the resultant statistics will be same as the ones if one calculates accumulating all the data upto that time instance and calculate the same statistics. All of these methods in each group is `OpenMP` and `MPI` compatible, and tested.

Following table summarize capabilities of Statistics Application.

| Statistical Methods                   | Norm Types              | Spatial Domain                                                | Temporal Domain                                                              | Data types |
|---------------------------------------|-------------------------|---------------------------------------------------------------|------------------------------------------------------------------------------|------------|
| [Sum](#sum)                           | [Value](#value)         | [Spatial methods](#spatial-methods)                           | [Temporal methods](#temporal-methods)                                        | Double     |
| [Mean](#mean)                         | [Magnitude](#magnitude) | [Spatial containers](#spatial-containers)                     | [Temporal containers](#temporal-containers)                                  | Array 3D   |
| [Root mean square](#root-mean-square) | [Euclidean](#euclidean) | [nodal_historical](#spatial-nodal-historical)                 | [nodal_historical_historical](#temporal-nodal-historical-historical)         | Vector     |
| [Variance](#variance)                 | [Infinity](#infinity)   | [nodal_non_historical](#spatial-nodal-non-historical)         | [nodal_historical_non_historical](#temporal-nodal-historical-non-historical) | Matrix     |
| [Min](#min)                           | [P-Norm](#p-norm)       | [condition_non_historical](#spatial-condition-non-historical) | [nodal_non_historical](#temporal-nodal-non-historical)                       |            |
| [Max](#max)                           | [Lpq-Norm](#lpq-norm)   | [element_non_historical](#spatial-element-non-historical)     | [element_non_historical](#temporal-element-non-historical)                   |            |
| [Median](#median)                     | [Frobenius](#frobenius) | [Spatial statistics process](#spatial-statistics-process)     | [condition_non_historical](#temporal-condition-non-historical)               |            |
| [Distribution](#distribution)         | [Trace](#trace)         |                                                               | [Temporal statistics process](#temporal-statistics-process)                  |            |
| [Norm methods](#norm-methods)         | [Index](#index-based)         |                                                               |                                                                              |            |
|                                       | [Component](#component-based) |                                                               |                                                                              |            |

## JSON Examples

If you prefer to use statistics of a simulation using StatisticsApplication, there is `spatial_statistics_process` for spatial statistics calculations and `temporal_statistics_process` for temporal statistics. These processes can be included via JSON settings under `auxiliary_processes`.

### Spatial statistics process examples

Following example illustrates different methods used in different containers with different norms. `input_variable_settings` holds an array of methods for specified containers, specified norm and specified variables. They can be customized for your requirement. `output_settings` holds information about how the output should be handled.

```json
{
                "kratos_module" : "KratosMultiphysics.StatisticsApplication",
                "python_module" : "spatial_statistics_process",
                "Parameters" : {
                    "model_part_name" : "test_model_part",
                    "input_variable_settings" : [
                        {
                            "method_name"    : "sum",
                            "norm_type"      : "none",
                            "container"      : "nodal_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "mean",
                            "norm_type"      : "none",
                            "container"      : "nodal_non_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "variance",
                            "norm_type"      : "none",
                            "container"      : "element_non_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "rootmeansquare",
                            "norm_type"      : "none",
                            "container"      : "condition_non_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "sum",
                            "norm_type"      : "magnitude",
                            "container"      : "nodal_historical",
                            "variable_names" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "mean",
                            "norm_type"      : "pnorm_2.5",
                            "container"      : "nodal_non_historical",
                            "variable_names" : ["VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "variance",
                            "norm_type"      : "component_x",
                            "container"      : "condition_non_historical",
                            "variable_names" : ["VELOCITY"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "rootmeansquare",
                            "norm_type"      : "index_3",
                            "container"      : "nodal_non_historical",
                            "variable_names" : ["LOAD_MESHES"],
                            "method_settings": {}
                        },
                        {
                            "method_name"    : "min",
                            "norm_type"      : "frobenius",
                            "container"      : "nodal_non_historical",
                            "variable_names" : ["GREEN_LAGRANGE_STRAIN_TENSOR"],
                            "method_settings": {}
                        }
                    ],
                    "output_settings" : {
                        "output_control_variable": "STEP",
                        "output_time_interval"   : 1,
                        "write_kratos_version"   : false,
                        "write_time_stamp"       : false,
                        "output_file_settings"   : {
                            "file_name"  : "<model_part_name>_<container>_<norm_type>_<method_name>.dat",
                            "output_path": "spatial_statistics_process",
                            "write_buffer_size" : -1
                        }
                    }
                }
            }
```

### Temporal statistics process example

Following is an example on how to use `temporal_statistics_process` in json. The output variables data type should be matched with the input variables order and the data type for `"norm_type" = "none"`. If `"norm_type" != "none"`, then output variables should be `scalars`. If `"container" = "nodal_historical_historical"` is used as the container type, then output variables should be added to `NodalSolutionStepVariables` list in Kratos since this container type outputs temporal statistics variables to nodal historical container. This json settings also can be added to `auxiliary_processes` list.

For details about all the available statistical methods, norm_types, etc, please refer to rest of the `README.md` file.

```json
        {
            "kratos_module": "KratosMultiphysics.StatisticsApplication",
            "python_module": "temporal_statistics_process",
            "Parameters": {
                "model_part_name": "FluidModelPart.fluid_computational_model_part",
                "input_variable_settings": [
                    {
                        "method_name": "variance",
                        "norm_type": "none",
                        "container": "nodal_historical_non_historical",
                        "echo_level": 1,
                        "method_settings": {
                            "input_variables": [
                                "VELOCITY",
                                "PRESSURE"
                            ],
                            "output_mean_variables": [
                                "VECTOR_3D_MEAN",
                                "SCALAR_MEAN"
                            ],
                            "output_variance_variables": [
                                "VECTOR_3D_VARIANCE",
                                "SCALAR_VARIANCE"
                            ]
                        }
                    }
                ],
                "statistics_start_point_control_variable_name": "TIME",
                "statistics_start_point_control_value": 2.5
            }
        }
```

## Method definitions

There are two types of methods under each **Spatial** and **Temporal** method groups. They are namely **Value** and **Norm** methods. In these methods, `i` index refers to spatial domain element index, and `k` index refers to time step.

### Value methods

#### Sum

In the case of spatial domain, it adds up all the variable values for a given container and returns summed up value as shown in following equation. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{r}&space;=&space;\sum_{i=1}^N{\underline{x}_i}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{r}&space;=&space;\sum_{i=1}^N{\underline{x}_i}}" title="\color{Black}{\underline{r} = \sum_{i=1}^N{\underline{x}_i}}" /></a>

Following is an example of summation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Sum(model_part, Kratos.VELOCITY)
```

In the case of temporal domain, **Sum** methods is the time integrated quantity for a specific variable. It will be stored each element under user specified variable and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{r}&space;=&space;\sum_{k=1}^{P}{\underline{x}_k\Delta&space;t_k}&space;\quad&space;where&space;\quad&space;\Delta&space;t_k&space;=&space;T_{k}&space;-&space;T_{k-1}&space;\quad&space;\forall&space;T_k&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{r}&space;=&space;\sum_{k=1}^{P}{\underline{x}_k\Delta&space;t_k}&space;\quad&space;where&space;\quad&space;\Delta&space;t_k&space;=&space;T_{k}&space;-&space;T_{k-1}&space;\quad&space;\forall&space;T_k&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{\underline{r} = \sum_{k=1}^{P}{\underline{x}_k\Delta t_k} \quad where \quad \Delta t_k = T_{k} - T_{k-1} \quad \forall T_k \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

Following is an example of integration calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `DISPLACEMENT` where integrated value will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
integration_starting_time = 2.0
sum_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    sum_method.CalculateStatistics()
```

#### Mean

In the case of spatial domain, it calculates mean of a given variable for a given container and returns it as shown in following equation. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user. (If it has higher dimension than a scalar, mean of each dimension will be calculated seperately resulting with a mean having same dimension as the input dimension)

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{r}&space;=&space;\frac{1}{N}\sum_{i=1}^{N}\underline{x}_i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{r}&space;=&space;\frac{1}{N}\sum_{i=1}^{N}\underline{x}_i}" title="\color{Black}{\underline{r} = \frac{1}{N}\sum_{i=1}^{N}\underline{x}_i}" /></a>

Following is an example of mean calculation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
mean = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Mean(model_part, Kratos.VELOCITY)
```

In the case of temporal domain, **Mean** methods is the time integrated quantity's mean for a specific variable. It will be stored each element under user specified variable and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{\bar{x}}&space;=&space;\frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k\Delta&space;t_k}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_k&space;=&space;T_{k}&space;-&space;T_{k-1}&space;\quad&space;\forall&space;T_k&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{\bar{x}}&space;=&space;\frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k\Delta&space;t_k}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_k&space;=&space;T_{k}&space;-&space;T_{k-1}&space;\quad&space;\forall&space;T_k&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{\underline{\bar{x}} = \frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k\Delta t_k} \quad where \quad T_{total} = T_{end} - T_{initial} \quad and \quad \Delta t_k = T_{k} - T_{k-1} \quad \forall T_k \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

Following is an example of mean calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `VECTOR_3D_MEAN` where mean will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
mean_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Mean.Array(model_part, "", Kratos.VELOCITY, 0, KratosStats.VECTOR_3D_MEAN)
integration_starting_time = 2.0
mean_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    mean_method.CalculateStatistics()
```

#### Root mean square

In the case of spatial domain, it calculates root mean square of a given variable for a given container and returns it as shown in following equation. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user. (If it has higher dimension than a scalar, root mean square of each dimension will be calculated seperately resulting with a mean having same dimension as the input dimension)

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{r}&space;=&space;\sqrt{\frac{1}{N}\sum_{i=1}^{N}{\underline{x}^2_i}&space;}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{r}&space;=&space;\sqrt{\frac{1}{N}\sum_{i=1}^{N}{\underline{x}^2_i}&space;}}" title="\color{Black}{\underline{r} = \sqrt{\frac{1}{N}\sum_{i=1}^{N}{\underline{x}^2_i} }}" /></a>

Following is an example of root mean square calculation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
rms = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.RootMeanSquare(model_part, Kratos.VELOCITY)
```

In the case of temporal domain, **Root Mean Square** methods is the time integrated quantity's root mean square for a specific variable. It will be stored each element under user specified variable and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{r}&space;=&space;\sqrt{\frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k^2\Delta&space;t_k}}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_k&space;=&space;T_{k}&space;-&space;T_{k-1}&space;\quad&space;\forall&space;T_k&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{r}&space;=&space;\sqrt{\frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k^2\Delta&space;t_k}}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_k&space;=&space;T_{k}&space;-&space;T_{k-1}&space;\quad&space;\forall&space;T_k&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{\underline{r} = \sqrt{\frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k^2\Delta t_k}} \quad where \quad T_{total} = T_{end} - T_{initial} \quad and \quad \Delta t_k = T_{k} - T_{k-1} \quad \forall T_k \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

Following is an example of root mean square calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `VECTOR_3D_MEAN` where root mean square value will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
rms_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.RootMeanSquare.Array(model_part, "", Kratos.VELOCITY, 0, KratosStats.VECTOR_3D_MEAN)
integration_starting_time = 2.0
rms_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    rms_method.CalculateStatistics()
```

#### Variance

In the case of spatial domain, it calculates variance of a given variable for a given container and returns mean and variance as shown in following equations. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will have the same type as the type of the variable specified by the user. (If it has higher dimension than a scalar, mean of each dimension will be calculated seperately resulting with a mean having same dimension as the input dimension)

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{\bar{x}}&space;=&space;\frac{1}{N}\sum_{i=1}^N{\underline{x}_i}}&space;\\&space;\color{Black}{\underline{v}&space;=&space;\frac{1}{N}\sum_{i=1}^N{\left(\underline{x}_i&space;-&space;\underline{\bar{x}}&space;\right&space;)^2}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{\bar{x}}&space;=&space;\frac{1}{N}\sum_{i=1}^N{\underline{x}_i}}&space;\\&space;\color{Black}{\underline{v}&space;=&space;\frac{1}{N}\sum_{i=1}^N{\left(\underline{x}_i&space;-&space;\underline{\bar{x}}&space;\right&space;)^2}}" title="\color{Black}{\underline{\bar{x}} = \frac{1}{N}\sum_{i=1}^N{\underline{x}_i}} \\ \color{Black}{\underline{v} = \frac{1}{N}\sum_{i=1}^N{\left(\underline{x}_i - \underline{\bar{x}} \right )^2}}" /></a>

Following is an example of variance calculation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
mean, variance = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Variance(model_part, Kratos.VELOCITY)
```

In the case of temporal domain, **Variance** method is the time integrated quantity's variance for a specific variable. Mean and variance will be stored each element under user specified variables and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{\bar{x}}&space;=&space;\frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k\Delta&space;t_k}}&space;\\&space;\color{Black}{Var\left(\underline{x}&space;\right&space;)&space;=&space;\frac{1}{T_{total}}\sum_{k=1}^{P}{\left(\underline{x}_k&space;-&space;\underline{\bar{x}}&space;\right&space;)^2\Delta&space;t_k}}&space;\\&space;\\&space;\color{Black}{&space;\quad&space;where}&space;\\&space;\\&space;\color{Black}{&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_k&space;=&space;T_{k}&space;-&space;T_{k-1}&space;\quad&space;\forall&space;T_k&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{\bar{x}}&space;=&space;\frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k\Delta&space;t_k}}&space;\\&space;\color{Black}{Var\left(\underline{x}&space;\right&space;)&space;=&space;\frac{1}{T_{total}}\sum_{k=1}^{P}{\left(\underline{x}_k&space;-&space;\underline{\bar{x}}&space;\right&space;)^2\Delta&space;t_k}}&space;\\&space;\\&space;\color{Black}{&space;\quad&space;where}&space;\\&space;\\&space;\color{Black}{&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_k&space;=&space;T_{k}&space;-&space;T_{k-1}&space;\quad&space;\forall&space;T_k&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{\underline{\bar{x}} = \frac{1}{T_{total}}\sum_{k=1}^{P}{\underline{x}_k\Delta t_k}} \\ \color{Black}{Var\left(\underline{x} \right ) = \frac{1}{T_{total}}\sum_{k=1}^{P}{\left(\underline{x}_k - \underline{\bar{x}} \right )^2\Delta t_k}} \\ \\ \color{Black}{ \quad where} \\ \\ \color{Black}{ T_{total} = T_{end} - T_{initial} \quad and \quad \Delta t_k = T_{k} - T_{k-1} \quad \forall T_k \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

Following is an example of root mean square calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable `VECTOR_3D_MEAN` will store mean and `VECTOR_3D_VARIANCE` will store variance in same container for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
variance_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Variance.Array(model_part, "", Kratos.VELOCITY, 0, KratosStats.VECTOR_3D_MEAN, KratosStats.VECTOR_3D_VARIANCE)
integration_starting_time = 2.0
variance_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    variance_method.CalculateStatistics()
```

#### Min

In the case of spatial domain, it returns minimum of a given variable's norm for a given container, and the corresponding items id. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will be double and integer, irrespective of the input type, since higher dimensional variable types will be reduced to scalars by the use of norms.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{v&space;=&space;\min_{\underline{x}_i&space;\in&space;\mathbf{T}}&space;|\underline{x}_i|}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{v&space;=&space;\min_{\underline{x}_i&space;\in&space;\mathbf{T}}&space;|\underline{x}_i|}" title="\color{Black}{v = \min_{\underline{x}_i \in \mathbf{T}} |\underline{x}_i|}" /></a>

Following is an example of min method of non historical `VELOCITY`'s magnitude over the whole model part's nodes. It returns a tuple, first argument being the minimum, and the second argument being the id of the node where the minimum is found.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
min_value, min_id = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Min(model_part, Kratos.VELOCITY, "magnitude")
```

In the case of temporal domain, **Min** method returns minimum value in the temporal domain, and the time minimum is found. Minimum and its occurring time will be stored each element under user specified variables and a user specified container. x<sub>k</sub> is the k<sup>th</sup> time step element's variable value of the corresponding container. Results will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{v&space;=&space;\min_{\underline{x}_k&space;\in&space;\mathbf{T}}&space;|\underline{x}_k|}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{v&space;=&space;\min_{\underline{x}_k&space;\in&space;\mathbf{T}}&space;|\underline{x}_k|}" title="\color{Black}{v = \min_{\underline{x}_k \in \mathbf{T}} |\underline{x}_k|}" /></a>

Following is an example of min method in non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable `VECTOR_3D_NORM` will store minimum and `TIME` will store the time minimum occured for each node. The `0` represents echo level for this method object. "magnitude" indicates that magnitude norm is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
min_method = KratosStats.TemporalMethods.NonHistorical.Nodes.NormMethods.Min.Array(model_part, "magnitude", Kratos.VELOCITY, 0, KratosStats.VECTOR_3D_NORM, Kratos.TIME)
integration_starting_time = 2.0
min_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    min_method.CalculateStatistics()
```

#### Max

In the case of spatial domain, it returns maximum of a given variable's norm for a given container, and the corresponding items id. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will be double and integer, irrespective of the input type, since higher dimensional variable types will be reduced to scalars by the use of norms.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{v&space;=&space;\max_{\underline{x}_i&space;\in&space;\Omega}&space;|\underline{x}_i|}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{v&space;=&space;\max_{\underline{x}_i&space;\in&space;\Omega}&space;|\underline{x}_i|}" title="\color{Black}{v = \max_{\underline{x}_i \in \Omega} |\underline{x}_i|}" /></a>

Following is an example of max method of non historical `VELOCITY`'s magnitude over the whole model part's nodes. It returns a tuple, first argument being the maximum, and the second argument being the id of the node where the maximum is found.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
max_value, max_id = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Max(model_part, Kratos.VELOCITY, "magnitude")
```

In the case of temporal domain, **Max** method returns maximum value in the temporal domain, and the time maximum is found. Maximum and its occurring time will be stored each element under user specified variables and a user specified container. x<sub>k</sub> is the k<sup>th</sup> time step element's variable value of the corresponding container. Results will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{v&space;=&space;\max_{\underline{x}_k&space;\in&space;\mathbf{T}}&space;|\underline{x}_k|}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{v&space;=&space;\max_{\underline{x}_k&space;\in&space;\mathbf{T}}&space;|\underline{x}_k|}" title="\color{Black}{v = \max_{\underline{x}_k \in \mathbf{T}} |\underline{x}_k|}" /></a>

Following is an example of max method in non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable `VECTOR_3D_NORM` will store maximum and `TIME` will store the time maximum occured for each node. The `0` represents echo level for this method object. "magnitude" indicates that magnitude norm is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
max_method = KratosStats.TemporalMethods.NonHistorical.Nodes.NormMethods.Max.Array(model_part, "magnitude", Kratos.VELOCITY, 0, KratosStats.VECTOR_3D_NORM, Kratos.TIME)
integration_starting_time = 2.0
max_method.InitializeStatisticsMethod(integration_starting_time)
for t in range(3, 6):
    max_method.CalculateStatistics()
```

#### Median

Median returns the median value in spatial domain of a given variable's norm for a given container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will be double value type, irrespective of the input type, since higher dimensional variable types will be reduced to scalars by the use of norms.

Following is an example of median method of non historical `VELOCITY`'s magnitude over the whole model part's nodes. It returns a double which is the median.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
median_value = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Median(model_part, Kratos.VELOCITY, "magnitude")
```

#### Distribution

Distribution methods calculates distribution of a given variable with respect to given norm in a given container in spatial domain. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will be a tuple with followings in the same order:

1. min in the domain or user specified min value
2. max in the domain or user specified min value
3. group limits of the distribution (only upper limit) (double array)
4. number of occurences of items within each group (int array)
5. percentage distribution of number of occurences of items within each group (double array)
6. Each group's seperate means
7. Eash group's seperate variances

This method requires followings as the parameters as a object of Kratos.Parameters, if nothing is provided, then followings are assumed as defaults.

```json
{
    "number_of_value_groups" : 10,
    "min_value"              : "min",
    "max_value"              : "max"
}
```

In here, `"min_value"` can either be `"min"` or a `double` value. In the case of a double, it will be used as the minimum when creating groups to identify distribution. `"max_value"` also can either be `"max"` of `double` value, which will determine the maximum value when creating the groups. This will create 2 additional groups to represent values below the specified `"min_value"` and `"max_value"` apart from the `"number_of_value_groups"` specified by the user.

Following is an example of median method of non historical `VELOCITY`'s magnitude over the whole model part's nodes.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
min_value, max_value, group_upper_values, group_histogram, group_percentage_distribution, group_means, group_variances = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Distribution(model_part, Kratos.VELOCITY, "magnitude")
```

### Norm methods

All of the value methods mentioned before supports norm version of it, in these methods, higher dimensional values are transformed in to scalar values by a use specified norm, and then statistics are calculated based on the chosen method. Supported norm types may differ based on the variable data type being used. There are few methods, which only supports norm methods. Following table summarize availability of value and norm methods.

| Statistical Methods                   | Value Method | Norm Method |
|---------------------------------------|--------------|-------------|
| [Sum](#sum)                           | [x]          | [x]         |
| [Mean](#mean)                         | [x]          | [x]         |
| [Root mean square](#root-mean-square) | [x]          | [x]         |
| [Variance](#variance)                 | [x]          | [x]         |
| [Min](#min)                           |              | [x]         |
| [Max](#max)                           |              | [x]         |
| [Median](#median)                     |              | [x]         |
| [Distribution](#distribution)         |              | [x]         |

Following example shows variance method, under **value** category and **norm** category for nodal non historical velocity. Norm method uses `"magnitude"` as the norm to reduce `VELOCITY` to a scalar value. `value_mean` and `value_variance` will be of `Array 3D` type, whereas `norm_mean` and `norm_variance` will be of `Double` type.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
value_mean, value_variance = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Variance(model_part, Kratos.VELOCITY)
norm_mean, norm_variance = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Variance(model_part, Kratos.VELOCITY, "magnitude")
```

## Norm definitions

Few different norms are predefined in this application. Following table summarize

|                               | Double | Array 3D | Vector | Matrix |
|-------------------------------|--------|----------|--------|--------|
| [Value](#value)               | [x]    |          |        |        |
| [Magnitude](#magnitude)       | [x]    | [x]      | [x]    | [x]    |
| [Euclidean](#euclidean)       |        | [x]      | [x]    |        |
| [Infinity](#infinity)         |        | [x]      | [x]    | [x]    |
| [P-Norm](#p-norm)             |        | [x]      | [x]    | [x]    |
| [Lpq-Norm](#lpq-norm)         |        |          |        | [x]    |
| [Frobenius](#frobenius)       |        |          |        | [x]    |
| [Trace](#trace)               |        |          |        | [x]    |
| [Index](#index-based)         |        |          | [x]    | [x]    |
| [Component](#component-based) |        | [x]      |        |        |

### Value

This returns the exact value. Only available for `Double` type variables.

### Magnitude

This returns the second norm of the variable.

1. For a `Double` type, it returns the absolute value.
2. For a `Array 3D` type, it returns the magnitude of the 3D vector.
3. For a `Vector` type, it returns root square sum of the vector.
4. For a `Matrix` type, it returns the [frobenius](#frobenius) norm.

### Euclidean

This again returns the second norm of the variable for following variable types.

1. For a `Array 3D` type, it returns the magnitude of the 3D vector.
2. For a `Vector` type, it returns root square sum of the vector.

### Frobenius

This is only available for `Matrix` data type variables only. Following equation illustrates the norm, where `A` is the matrix and a<sub>ij</sub> is the i<sup>th</sup> row, j<sup>th</sup> column value.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||A||_F&space;=&space;\sqrt{\sum_i&space;{\sum_j{||a_{ij}||^2}}}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||A||_F&space;=&space;\sqrt{\sum_i&space;{\sum_j{||a_{ij}||^2}}}}" title="{\color{Black} ||A||_F = \sqrt{\sum_i {\sum_j{||a_{ij}||^2}}}}" /></a>

### Infinity

This is infinity (i.e. max norm) which is available for `Array 3D`, `Vector` and `Matrix` type.
For an `Array 3D` variable following equation is used.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{V}||_\infty&space;=&space;\max\left&space;\lbrace&space;{|V_X|,&space;|V_Y|,&space;|V_Z|}\right&space;\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{V}||_\infty&space;=&space;\max\left&space;\lbrace&space;{|V_X|,&space;|V_Y|,&space;|V_Z|}\right&space;\rbrace}" title="{\color{Black} ||\underline{V}||_\infty = \max\left \lbrace {|V_X|, |V_Y|, |V_Z|}\right \rbrace}" /></a>

For a `Vector` variable following equation is used.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{V}||_\infty&space;=&space;\max\left&space;\lbrace&space;{|V_0|,&space;|V_1|,&space;|V_2|,...,{V_{n-1}}}\right&space;\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{V}||_\infty&space;=&space;\max\left&space;\lbrace&space;{|V_0|,&space;|V_1|,&space;|V_2|,...,{V_{n-1}}}\right&space;\rbrace}" title="{\color{Black} ||\underline{V}||_\infty = \max\left \lbrace {|V_0|, |V_1|, |V_2|,...,{V_{n-1}}}\right \rbrace}" /></a>

For a `Matrix` variable following equation is used.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{A}||_\infty&space;=&space;\max_{0&space;\leq&space;i&space;<&space;n}\left(&space;\sum_{j=0}^{n-1}&space;|a_{ij}|\right&space;)}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{A}||_\infty&space;=&space;\max_{0&space;\leq&space;i&space;<&space;n}\left(&space;\sum_{j=0}^{n-1}&space;|a_{ij}|\right&space;)}" title="{\color{Black} ||\underline{A}||_\infty = \max_{0 \leq i < n}\left( \sum_{j=0}^{n-1} |a_{ij}|\right )}" /></a>

### Trace

Trace is only for `Matrix` type variables. If the given matrix is not squre, an error is thrown.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{A}||_{trace}&space;=&space;\sum_{i=0}^{n-1}a_{ij}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{A}||_{trace}&space;=&space;\sum_{i=0}^{n-1}a_{ij}}" title="{\color{Black} ||\underline{A}||_{trace} = \sum_{i=0}^{n-1}a_{ij}}" /></a>

### P norm

P norm is applicable for `Array 3D`, `Vector` and `Matrix` variables. The p-norm is used as `pnorm_p` where, `p` represents the value to be used as `p` in the p-norm equation given below. `p` should be greater than or equal to 1.0.

For `Array 3D` and `Vector`:

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{V}||_{p}&space;=&space;\left(\sum_{i=0}^{n-1}|v_{i}|^p&space;\right&space;)^{1/p}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{V}||_{p}&space;=&space;\left(\sum_{i=0}^{n-1}|v_{i}|^p&space;\right&space;)^{1/p}}" title="{\color{Black} ||\underline{V}||_{p} = \left(\sum_{i=0}^{n-1}|v_{i}|^p \right )^{1/p}}" /></a>

For `Matrix`:

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{A}||_{p}&space;=&space;\left(\sum_{i=0}^{n-1}\sum_{j=0}^{n-1}|a_{ij}|^p&space;\right&space;)^{1/p}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{A}||_{p}&space;=&space;\left(\sum_{i=0}^{n-1}\sum_{j=0}^{n-1}|a_{ij}|^p&space;\right&space;)^{1/p}}" title="{\color{Black} ||\underline{A}||_{p} = \left(\sum_{i=0}^{n-1}\sum_{j=0}^{n-1}|a_{ij}|^p \right )^{1/p}}" /></a>

### Lpq norm

This norm is only applicable to `Matrix` type variables.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{A}||_{lpq}&space;=&space;\left(\sum_{j=0}^{n-1}\left(\sum_{i=0}^{n-1}|a_{ij}|^p&space;\right&space;)^{q/p}&space;\right&space;)^{1/q}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{A}||_{lpq}&space;=&space;\left(\sum_{j=0}^{n-1}\left(\sum_{i=0}^{n-1}|a_{ij}|^p&space;\right&space;)^{q/p}&space;\right&space;)^{1/q}}" title="{\color{Black} ||\underline{A}||_{lpq} = \left(\sum_{j=0}^{n-1}\left(\sum_{i=0}^{n-1}|a_{ij}|^p \right )^{q/p} \right )^{1/q}}" /></a>

### Index based

Index based norms are available for both `Vector` and `Matrix` based variable types. `index_i` is used for `Vector` variables, where `i` represents the index of the vector to be used as the scalar representing a corresponding vector. `index_(i,j)` is used for `Matrix` variables, where `i` represents the row index and `j` represents the column index of the matrix value which to be used as the scalar representation of the given matrix.

### Component based

Component based norms are available only for `Array 3D` type variables. `component_x` represents the `x` component of 3D variable, `component_y` represents `y` component, `component_z` represents the `z` component.


## Spatial methods

All the spatial statistical methods can be found in the `SpatialMethods` submodule. `ValueMethods` container will contain all available spatial value methods, and `NormMethods` container will contain all available norm methods.

Following example illustrates all the available value methods executed on nodal non historical variable velocity.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Sum(model_part, Kratos.VELOCITY)
mean_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Mean(model_part, Kratos.VELOCITY)
rms_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.RootMeanSquare(model_part, Kratos.VELOCITY)
mean_value, variance_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Variance(model_part, Kratos.VELOCITY)
```

Following example illustrates all the available norm methods executed on nodal non historical variable velocity's magnitude.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
mean_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Mean(model_part, Kratos.VELOCITY, "magnitude")
rms_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.RootMeanSquare(model_part, Kratos.VELOCITY, "magnitude")
mean_norm, variance_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Variance(model_part, Kratos.VELOCITY, "magnitude")
min_norm, min_norm_id = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Min(model_part, Kratos.VELOCITY, "magnitude")
max_norm, max_norm_id = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Max(model_part, Kratos.VELOCITY, "magnitude")
min_norm, max_norm, group_upper_norms, group_histogram, group_percentage_distribution = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Distribution(model_part, Kratos.VELOCITY, "magnitude")
```

### Spatial containers

Four different types of containers are supported for spatial methods.

1. [nodal_historical](#spatial-nodal-historical)
2. [nodal_non_historical](#spatial-nodal-non-historical)
3. [condition_non_historical](#spatial-condition-non-historical)
4. [element_non_historical](#spatial-element-non-historical)

#### Spatial nodal historical

Nodal historical containers methods can be accessed through `Historical` submodule in `SpatialMethods`. This calculates chosen statistics of the variable in the nodal historical data. Following example illustrates use of value and norm methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
sum_value = KratosStats.SpatialMethods.Historical.ValueMethods.Sum(model_part, Kratos.VELOCITY)
sum_norm = KratosStats.SpatialMethods.Historical.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
```

#### Spatial nodal non historical

Nodal non historical containers methods can be accessed through `NonHistorical.Nodes` submodule in `SpatialMethods`. This calculates chosen statistics of the variable in the nodal non historical data. Following example illustrates use of value and norm methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Sum(model_part, Kratos.VELOCITY)
sum_norm = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
```

#### Spatial condition non historical

Condition non historical containers methods can be accessed through `NonHistorical.Conditions` submodule in `SpatialMethods`. This calculates chosen statistics of the variable in the condition non historical data. Following example illustrates use of value and norm methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.SpatialMethods.NonHistorical.Conditions.ValueMethods.Sum(model_part, Kratos.VELOCITY)
sum_norm = KratosStats.SpatialMethods.NonHistorical.Conditions.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
```

#### Spatial element non historical

Element non historical containers methods can be accessed through `NonHistorical.Elements` submodule in `SpatialMethods`. This calculates chosen statistics of the variable in the element non historical data. Following example illustrates use of value and norm methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.SpatialMethods.NonHistorical.Elements.ValueMethods.Sum(model_part, Kratos.VELOCITY)
sum_norm = KratosStats.SpatialMethods.NonHistorical.Elements.NormMethods.Sum(model_part, Kratos.VELOCITY, "magnitude")
```

### Spatial statistics process

This is a seperate process, which can be included in the json file as an auxiliary process. Followings are the defaults used in this process.

```json
{
    "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
    "input_variable_settings" : [
        {
            "method_name"    : "sum",
            "norm_type"      : "none",
            "container"      : "nodal_historical",
            "variable_names" : [],
            "method_settings": {}
        }
    ],
    "output_settings" : {
        "output_control_variable": "STEP",
        "output_time_interval"   : 1,
        "write_kratos_version"   : true,
        "write_time_stamp"       : true,
        "output_file_settings"   : {
            "file_name"  : "<model_part_name>_<container>_<norm_type>_<method_name>.dat",
            "output_path": "spatial_statistics_output",
            "write_buffer_size" : -1
        }
    }
}
```

`model_part_name` indicates which model part to be operated on for statistics calculation. `input_variable_settings` contains list of statistics and/or norm methods along with their acting containers of the model part. `method_name` is the statistics method name (always everything is in lower case). `norm_type` is the applied norm, `none` means, value methods are used, otherwise specified norm is used. `container` is the [container](#spatial-method-containers) to be specified for statistical value calculation. `method_settings` is used to pass additional information required by statistics method such as in [Distribution](#distribution) method. Output is written to a file in the format given at `file_name`, under the folder `folder_name`. This process can be used, if someone prefers not to modify their python scripts to calculate statistics as mentioned in the previous sections

## Temporal methods

All the temporal methods are available through `TemporalMethods` submodule under `StatisticsApplication`. In the case of temporal methods, seperate objects needs to be created for each input variable for different type of variables as shown in following example for nodal non historical data container value methods.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_double_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Double(model_part, "", Kratos.PRESSURE, 0, Kratos.DENSITY)
sum_array_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_vector_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Vector(model_part, "", Kratos.LOAD_MESHES, 0, Kratos.MATERIAL_PARAMETERS)
sum_matrix_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Matrix(model_part, "", Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, 0, Kratos.CAUCHY_STRESS_TENSOR)
```

In order to remove the hassle of using different names (such as `Double`, `Array`...), one can get the advantage of common methods available in each of the containers. Following example will create the same objects as in the previous example in one go using simplified sum method.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
params = Kratos.Parameters(r"""
{
    "input_variables" : ["PRESSURE", "VELOCITY", "LOAD_MESHES", "GREEN_LAGRANGE_STRAIN_TENSOR"],
    "output_variables" : ["DENSITY", "DISPLACEMENT", "MATERIAL_PARAMETERS", "CAUCHY_STRESS_TENSOR"]
}""")
sum_methods = KratosStats.TemporalMethods.NonHistorical.Nodes.Sum(model_part, "none", 0, params)
```

For **Sum**, **Mean**, **RootMeanSquare** methods, above `Kratos.Parameters` object keys will remain the same. But for **Variance** method following `json` parameters are used.

```json
{
    "input_variables"           : [],
    "output_mean_variables"     : [],
    "output_variance_variables" : []
}
```

When using these simplified, one needs to take care of compatibility of input and output variables depending on the method and norm type, otherwise runtime errors will be thrown.

### Temporal containers

In the temporal domain, there are five different containers available for user to calculate temporal statistics. They are,

1. [nodal_historical_historical](#temporal-nodal-historical-historical)
2. [nodal_historical_non_historical](#temporal-nodal-historical-non-historical)
3. [nodal_non_historical](#temporal-nodal-non-historical)
4. [condition_non_historical](#temporal-condition-non-historical)
5. [element_non_historical](#temporal-element-non-historical)

#### Temporal nodal historical historical

In this container, statistics are calculated on user specified variables in the nodal historical container and, output is also written to the nodal historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.Historical.HistoricalOutput.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.Historical.HistoricalOutput.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

#### Temporal nodal historical non historical

In this container, statistics are calculated on user specified variables in the nodal historical container and, output is also written to the nodal non historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.Historical.NonHistoricalOutput.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.Historical.NonHistoricalOutput.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

#### Temporal nodal non historical

In this container, statistics are calculated on user specified variables in the nodal non historical container and, output is also written to the nodal non historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.NonHistorical.Nodes.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

#### Temporal condition non historical

In this container, statistics are calculated on user specified variables in the condition non historical container and, output is also written to the condition non historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.NonHistorical.Conditions.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.NonHistorical.Conditions.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

#### Temporal element non historical

In this container, statistics are calculated on user specified variables in the element non historical container and, output is also written to the element non historical container. Example is given below.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_value = KratosStats.TemporalMethods.NonHistorical.Elements.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_norm = KratosStats.TemporalMethods.NonHistorical.Elements.NormMethods.Sum.Array(model_part, "magnitude", Kratos.VELOCITY, 0, Kratos.DENSITY)
```

### Temporal statistics process

This is a seperate process, which can be included in the json file as an auxiliary process. Followings are the defaults used in this process.

```json
    "model_part_name" : "PLEASE_SPECIFY_MODEL_PART_NAME",
    "input_variable_settings" : [
        {
            "method_name"     : "sum",
            "norm_type"       : "none",
            "container"       : "nodal_historical_non_historical",
            "echo_level"      : 0,
            "method_settings" : {}
        }
    ],
    "echo_level" : 0,
    "statistics_start_point_control_variable_name" : "TIME",
    "statistics_start_point_control_value"         : 0.0
```

`model_part_name` is the model part which all of the mentioned statistics will be calculated. `input_variable_settings` contains list of statistical methods, their norms and container combinations where statistics will be calculated. under each `input_variable_settings_block`, `method_name` is the statistics method, `norm_type` is the norm type (`"none"` for value methods, other wise norms are used), `container` to tell from which container to read/write, `method_settings` to specifiy input/output variables for the method. Apart from that, there are few global settings. `statistics_start_point_control_variable_name` is the temporal statistics calculation start control variable name. This variable should be present in the model_part specified in order to calculate temporal statistics. Statistics calculation will start when this variable's value passes user specified value in `statistics_start_point_control_value`.
