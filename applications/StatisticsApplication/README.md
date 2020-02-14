# Statistics Application

Statistics application consist of widely used methods to calculate statistics in various containers of KratosMultiphysics. There are mainly two groups of statistical methods namely **Spatial** and **Temporal**. **Spatial** methods calculate statistics on spatial containers and output the values whenever they are called. **Temporal** methods calculate statistics on the fly in a transient simulation. All the temporal methods gurantee that, the resultant statistics will be same as the ones if one calculates accumulating all the data upto that time instance and calculate the same statistics. All of these methods in each group is `OpenMP` and `MPI` compatible, and tested.

## Method definitions

There are two types of methods under each **Spatial** and **Temporal** method groups. They are namely **Value** and **Norm** methods.

### Value methods

#### Sum

In the case of spatial domain, it adds up all the variable values for a given container and returns summed up value as shown in following equation. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{sum&space;=&space;\sum_1^n{\underline{\phi}_i}}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{sum&space;=&space;\sum_1^n{\underline{\phi}_i}}" title="\color{Black}{sum = \sum_1^n{\underline{\phi}_i}}" /></a>

Following is an example of summation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Sum(model_part, Kratos.VELOCITY)
```

In the case of temporal domain, **Sum** methods is the time integrated quantity for a specific variable. It will be stored each element under user specified variable and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{sum&space;=&space;\sum_{i=1}^{n}{\underline{x}_i\Delta&space;t_i}&space;\quad&space;where&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{sum&space;=&space;\sum_{i=1}^{n}{\underline{x}_i\Delta&space;t_i}&space;\quad&space;where&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{sum = \sum_{i=1}^{n}{\underline{x}_i\Delta t_i} \quad where \quad \Delta t_i = T_{i} - T_{i-1} \quad \forall T_i \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

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

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{mean&space;=&space;\frac{1}{N}\sum_{i=1}^{N}\underline{x}_i}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{mean&space;=&space;\frac{1}{N}\sum_{i=1}^{N}\underline{x}_i}" title="\color{Black}{mean = \frac{1}{N}\sum_{i=1}^{N}\underline{x}_i}" /></a>

Following is an example of mean calculation of non historical `VELOCITY` over the whole model part's nodes

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
mean = KratosStats.SpatialMethods.NonHistorical.Nodes.ValueMethods.Mean(model_part, Kratos.VELOCITY)
```

In the case of temporal domain, **Mean** methods is the time integrated quantity's mean for a specific variable. It will be stored each element under user specified variable and a user specified container. x<sub>i</sub> is the i<sup>th</sup> element's variable value of the corresponding container. Result will have the same type as the type of the variable specified by the user preserving the dimensionality as in the spatial case.

<a href="https://www.codecogs.com/eqnedit.php?latex=\color{Black}{mean&space;=&space;\frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}_i\Delta&space;t_i}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" target="_blank"><img src="https://latex.codecogs.com/svg.latex?\color{Black}{mean&space;=&space;\frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}_i\Delta&space;t_i}&space;\quad&space;where&space;\quad&space;T_{total}&space;=&space;T_{end}&space;-&space;T_{initial}&space;\quad&space;and&space;\quad&space;\Delta&space;t_i&space;=&space;T_{i}&space;-&space;T_{i-1}&space;\quad&space;\forall&space;T_i&space;\in&space;\left\lbrace&space;T_{initial},&space;...,&space;T_{end}&space;\right\rbrace}" title="\color{Black}{mean = \frac{1}{T_{total}}\sum_{i=1}^{n}{\underline{x}_i\Delta t_i} \quad where \quad T_{total} = T_{end} - T_{initial} \quad and \quad \Delta t_i = T_{i} - T_{i-1} \quad \forall T_i \in \left\lbrace T_{initial}, ..., T_{end} \right\rbrace}" /></a>

Following is an example of mean calculation of non historical velocity. Input variable is node's non-historical container's `VELOCITY` and output variable is same containers `DISPLACEMENT` where mean will be stored for each node. The `0` represents echo level for this method object. Blank "" indicates that value method is used.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
sum_method = KratosStats.TemporalMethods.NonHistorical.Nodes.ValueMethods.Mean.Array(model_part, "", Kratos.VELOCITY, 0, Kratos.DISPLACEMENT)
sum_method.InitializeStatisticsMethod()
previous_time = 2
for t in range(3, 6):
    delta_time = t - previous_time
    sum_method.CalculateStatistics(delta_time)
    previous_time = t
```

#### Root mean square

#### Variance

#### Min

#### Max

#### Median

#### Distribution

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