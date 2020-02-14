# Statistics Application

Statistics application consist of widely used methods to calculate statistics in various containers of KratosMultiphysics. There are mainly two groups of statistical methods namely **Spatial** and **Temporal**. **Spatial** methods calculate statistics on spatial containers and output the values whenever they are called. **Temporal** methods calculate statistics on the fly in a transient simulation. All the temporal methods gurantee that, the resultant statistics will be same as the ones if one calculates accumulating all the data upto that time instance and calculate the same statistics. All of these methods in each group is `OpenMP` and `MPI` compatible, and tested.

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

Distribution methods calculates distribution of a given variable with respect to given norm in a given container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will be a tuple with followings in the same order:

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

## Norm definitions

### Value

### Magnitude

### Euclidean

### Infinity

### Trace

### P norm

### Lpq norm

### Index based

## Variable types and supported norms

ADD the table here...

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