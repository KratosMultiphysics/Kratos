---
title: Elements and Conditions
keywords: 
tags: [Python Script Tutorial ModelPart Elements Conditions]
sidebar: kratos_for_users
summary: 
---

The elements and conditions are the main extension points of the *Kratos*. In input script accessing them improves the flexibility in dealing with complex problems.

## Starting
First of all we need to create a python file with following code to import the *Kratos*, create a `ModelPart` and read it from input as described in the here :

```python
from KratosMultiphysics import *
import KratosMultiphysics.FluidDynamicsApplication

this_model = Model()
fluid_model_part = this_model.CreateModelPart("FluidPart")

fluid_model_part.AddNodalSolutionStepVariable(VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)
fluid_model_part.AddNodalSolutionStepVariable(TEMPERATURE)

fluid_model_part_io = ModelPartIO("path/to/file/example")
fluid_model_part_io.ReadModelPart(fluid_model_part)

fluid_model_part.SetBufferSize(3)
```
{: data-lang="Python"}

## Accessing `Elements`
The elements stored in the ModelPart can be accessed using the Elements parameter:

```python
model_part_elements = fluid_model_part.Elements
```
{: data-lang="Python"}
 
Iteration over all elements in a model part is very similar to the nodes. For example writing the elements in a model part can be done as follow:

```python
for element in fluid_model_part.Elements:
    print(element)
```
{: data-lang="Python"}

and printing the ID for all of the elements:

```python
for element in fluid_model_part.Elements:
    print(element.Id)
```
{: data-lang="Python"}

## Accessing `Conditions`
Conditions parameter of model part provides access to the conditions it stores:

```python
model_part_conditions = fluid_model_part.Conditions
```
{: data-lang="Python"}

Iteration over conditions is very similar to the elements. In the same way printing conditions is as follow:

```python
for condition in fluid_model_part.Conditions:
    print(condition)
```
{: data-lang="Python"}

and printing the ID for all of the conditions:

```python
for condition in fluid_model_part.Conditions:
    print(condition.Id)
```
{: data-lang="Python"}

**Prev** [Nodes and Nodal Data](Nodes_And_Data)