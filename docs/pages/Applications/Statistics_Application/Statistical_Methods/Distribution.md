---
title: Distribution
keywords: 
tags: [Distribution.md]
sidebar: statistics_application
summary: 
---
## Definition

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

## Example

Following is an example of median method of non historical `VELOCITY`'s magnitude over the whole model part's nodes.

```python
import KratosMultiphysics as Kratos
import KratosMultiphysics.StatisticsApplication as KratosStats
model = Kratos.Model()
model_part = model.CreateModelPart("test_model_part")
min_value, max_value, group_upper_values, group_histogram, group_percentage_distribution, group_means, group_variances = KratosStats.SpatialMethods.NonHistorical.Nodes.NormMethods.Distribution(model_part, Kratos.VELOCITY, "magnitude")
```