# Statistics Application

Statistics application consist of widely used methods to calculate statistics in various containers of KratosMultiphysics. There are mainly two groups of statistical methods namely **Spatial** and **Temporal**. **Spatial** methods calculate statistics on spatial containers and output the values whenever they are called. **Temporal** methods calculate statistics on the fly in a transient simulation. All the temporal methods gurantee that, the resultant statistics will be same as the ones if one calculates accumulating all the data upto that time instance and calculate the same statistics. All of these methods in each group is `OpenMP` and `MPI` compatible, and tested.

Following table summarize capabilities of Statistics Application.

| Statistical Methods                   | Norm Types | Spatial Containers       | Temporal Containers             | Data types |
|---------------------------------------|------------|--------------------------|---------------------------------|------------|
| [Sum](#sum)                           | Value      | nodal_historical         | nodal_historical_historical     | Double     |
| [Mean](#mean)                         | Magnitude  | nodal_non_historical     | nodal_historical_non_historical | Array 3D   |
| [Root mean square](#root-mean-square) | Euclidean  | element_non_historical   | nodal_non_historical            | Vector     |
| [Variance](#variance)                 | Infinity   | condition_non_historical | element_non_historical          | Matrix     |
| [Min](#min)                           | P-Norm     |                          | condition_non_historical        |            |
| [Max](#max)                           | Lpq-Norm   |                          |                                 |            |
| [Median](#median)                     | Frobenius  |                          |                                 |            |
| [Distribution](#distribution)         | Trace      |                          |                                 |            |
| [Norm methods](#norm-methods)         | Index      |                          |                                 |            |

## Method definitions

There are two types of methods under each **Spatial** and **Temporal** method groups. They are namely **Value** and **Norm** methods.

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

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{r}&space;=&space;\sum_{i=1}^{n}{\underline{x}\Delta&space;t_i}&space;\quad&space;where&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{r}&space;=&space;\sum_{i=1}^{n}{\underline{x}\Delta&space;t_i}&space;\quad&space;where&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{\underline{r} = \sum_{i=1}^{n}{\underline{x}\Delta t_i} \quad where \quad \Delta t_i = T_{i} - T_{i-1} \quad \forall T_i \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

Following is an example of integration calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `DISPLACEMENT` where integrated value will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Sum.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_method.InitializeStatisticsMethod()
previous_time = 2
for t in range(3, 6):
    delta_time = t - previous_time
    sum_method.CalculateStatistics(delta_time)
    previous_time = t
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

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{r}&space;=&space;\frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}\Delta&space;t_i}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{r}&space;=&space;\frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}\Delta&space;t_i}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{\underline{r} = \frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}\Delta t_i} \quad where \quad T_{total} = T_{end} - T_{initial} \quad and \quad \Delta t_i = T_{i} - T_{i-1} \quad \forall T_i \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

Following is an example of mean calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `VELOCITY_MEAN` where mean will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
mean_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Mean.Array(model_part, "", Kratos.VELOCITY, 0, KratosStats.VELOCITY_MEAN)
mean_method.InitializeStatisticsMethod()
previous_time = 2
for t in range(3, 6):
    delta_time = t - previous_time
    mean_method.CalculateStatistics(delta_time)
    previous_time = t
```

#### Root mean square

In the case of spatial domain, it calculates root mean square of a given variable for a given container and returns it as shown in following equation. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user. (If it has higher dimension than a scalar, mean of each dimension will be calculated seperately resulting with a mean having same dimension as the input dimension)

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

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{\underline{r}&space;=&space;\sqrt{\frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}^2\Delta&space;t_i}}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{\underline{r}&space;=&space;\sqrt{\frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}^2\Delta&space;t_i}}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{\underline{r} = \sqrt{\frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}^2\Delta t_i}} \quad where \quad T_{total} = T_{end} - T_{initial} \quad and \quad \Delta t_i = T_{i} - T_{i-1} \quad \forall T_i \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

Following is an example of root mean square calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `VELOCITY_MEAN` where root mean square value will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
rms_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.RootMeanSquare.Array(model_part, "", Kratos.VELOCITY, 0, KratosStats.VELOCITY_MEAN)
rms_method.InitializeStatisticsMethod()
previous_time = 2
for t in range(3, 6):
    delta_time = t - previous_time
    rms_method.CalculateStatistics(delta_time)
    previous_time = t
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

Following is an example of root mean square calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable `VELOCITY_MEAN` will store mean and `VELOCITY_VARIANCE` will store variance in same container for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
variance_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Variance.Array(model_part, "", Kratos.VELOCITY, 0, KratosStats.VELOCITY_MEAN, KratosStats.VELOCITY_VARIANCE)
variance_method.InitializeStatisticsMethod()
previous_time = 2
for t in range(3, 6):
    delta_time = t - previous_time
    variance_method.CalculateStatistics(delta_time)
    previous_time = t
```

#### Min

In the case of spatial domain, it returns minimum of a given variable's norm for a given container, and the corresponding items id. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Results will be double and integer, irrespective of the input type, since higher dimensional variable types will be reduced to scalars by the use of norms.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{v&space;=&space;\min_{\underline{x}_i&space;\in&space;\Omega}&space;|\underline{x}_i|}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{v&space;=&space;\min_{\underline{x}_i&space;\in&space;\Omega}&space;|\underline{x}_i|}" title="\color{Black}{v = \min_{\underline{x}_i \in \Omega} |\underline{x}_i|}" /></a>

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

Following is an example of min method in non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable `VELOCITY_NORM` will store minimum and `TIME` will store the time minimum occured for each node. The `0` represents echo level for this method object. "magnitude" indicates that magnitude norm is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
min_method = KratosStats.TemporalMethods.NonHistorical.Nodes.NormMethods.Min.Array(model_part, "magnitude", Kratos.VELOCITY, 0, KratosStats.VELOCITY_NORM, Kratos.TIME)
min_method.InitializeStatisticsMethod()
previous_time = 2
for t in range(3, 6):
    delta_time = t - previous_time
    min_method.CalculateStatistics(delta_time)
    previous_time = t
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

Following is an example of max method in non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable `VELOCITY_NORM` will store maximum and `TIME` will store the time maximum occured for each node. The `0` represents echo level for this method object. "magnitude" indicates that magnitude norm is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
max_method = KratosStats.TemporalMethods.NonHistorical.Nodes.NormMethods.Max.Array(model_part, "magnitude", Kratos.VELOCITY, 0, KratosStats.VELOCITY_NORM, Kratos.TIME)
max_method.InitializeStatisticsMethod()
previous_time = 2
for t in range(3, 6):
    delta_time = t - previous_time
    max_method.CalculateStatistics(delta_time)
    previous_time = t
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

This method requires followings as the parameters as a object of Kratos.Parameters, if nothing is provided, then followings are assumed as defaults.

```json
{
    "number_of_value_groups" : 10,
    "min_value"              : "min",
    "max_value"              : "max"
}
```

In here, `"min_value"` can either be `"min"` or a `double` value. In the case of a double, it will be used as the minimum when creating groups to identify distribution. `"max_value"` also can either be `"max"` of `double` value, which will determine the maximum value when creating the groups. This will create 2 additional groups to represent values below the specified `"min_value"` and `"max_value"` apart from the `"number_of_value_groups"` specified by the user.

Following is an example of median method of non historical `VELOCITY`'s magnitude over the whole model part's nodes. It returns a double which is the median.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
min_value, max_value, group_upper_values, group_histogram, group_percentage_distribution = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Distribution(model_part, Kratos.VELOCITY, "magnitude")
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

|           | Double | Array 3D | Vector | Matrix |
|-----------|--------|----------|--------|--------|
| Value     | [x]    |          |        |        |
| Magnitude | [x]    | [x]      | [x]    | [x]    |
| Euclidean |        | [x]      | [x]    |        |
| Infinity  |        | [x]      | [x]    | [x]    |
| P-Norm    |        | [x]      | [x]    | [x]    |
| Lpq-Norm  |        |          |        | [x]    |
| Frobenius |        |          |        | [x]    |
| Trace     |        |          |        | [x]    |
| Index     |        |          | [x]    | [x]    |
| Component |        | [x]      |        |        |

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

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{White}&space;||A||_F&space;=&space;\sqrt{\sum_i&space;{\sum_j{||a_{ij}||^2}}}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{White}&space;||A||_F&space;=&space;\sqrt{\sum_i&space;{\sum_j{||a_{ij}||^2}}}}" title="{\color{White} ||A||_F = \sqrt{\sum_i {\sum_j{||a_{ij}||^2}}}}" /></a>

### Infinity

This is infinity (i.e. max norm) which is available for `Array 3D`, `Vector` and `Matrix` type.
For an `Array 3D` variable following equation is used.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{White}&space;||\underline{V}||_\infty&space;=&space;\max\left&space;\lbrace&space;{|V_X|,&space;|V_Y|,&space;|V_Z|}\right&space;\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{White}&space;||\underline{V}||_\infty&space;=&space;\max\left&space;\lbrace&space;{|V_X|,&space;|V_Y|,&space;|V_Z|}\right&space;\rbrace}" title="{\color{White} ||\underline{V}||_\infty = \max\left \lbrace {|V_X|, |V_Y|, |V_Z|}\right \rbrace}" /></a>

For a `Vector` variable following equation is used.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{White}&space;||\underline{V}||_\infty&space;=&space;\max\left&space;\lbrace&space;{|V_0|,&space;|V_1|,&space;|V_2|,...,{V_{n-1}}}\right&space;\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{White}&space;||\underline{V}||_\infty&space;=&space;\max\left&space;\lbrace&space;{|V_0|,&space;|V_1|,&space;|V_2|,...,{V_{n-1}}}\right&space;\rbrace}" title="{\color{White} ||\underline{V}||_\infty = \max\left \lbrace {|V_0|, |V_1|, |V_2|,...,{V_{n-1}}}\right \rbrace}" /></a>

For a `Matrix` variable following equation is used.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{A}||_\infty&space;=&space;\max_{0&space;\leq&space;i&space;<&space;n}\left(&space;\sum_{j=0}^{n-1}&space;|a_{ij}|\right&space;)}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{A}||_\infty&space;=&space;\max_{0&space;\leq&space;i&space;<&space;n}\left(&space;\sum_{j=0}^{n-1}&space;|a_{ij}|\right&space;)}" title="{\color{Black} ||\underline{A}||_\infty = \max_{0 \leq i < n}\left( \sum_{j=0}^{n-1} |a_{ij}|\right )}" /></a>

### Trace

Trace is only for `Matrix` type variables. If the given matrix is not squre, an error is thrown.

<a href="https://www.codecogs.com/eqnedit.php?latex={\color{Black}&space;||\underline{A}||_{trace}&space;=&space;\sum_{i=0}^{n-1}a_{ij}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?{\color{Black}&space;||\underline{A}||_{trace}&space;=&space;\sum_{i=0}^{n-1}a_{ij}}" title="{\color{Black} ||\underline{A}||_{trace} = \sum_{i=0}^{n-1}a_{ij}}" /></a>

### P norm

### Lpq norm

### Index based

### Component based


## Spatial methods

Add supported spatial methods, with input arguments, output arguments

### Spatial method containers

Add spatial method containers with examples from c++ and python

### Spatial statistics process

## Temporal methods

Add supported temporal methods, with arguments, params

### Temporal method containers

Add temporal method containers and exampels from c++ and python

### Temporal statistics process

Add a description about the process