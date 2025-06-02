---
title: Nodes and Nodal Data
keywords: 
tags: [Python Script Tutorial Nodes Nodal Data]
sidebar: kratos_for_users
summary: 
---

In this tutorial the access to the nodes stored in a `ModelPart` and their nodal data will be described. More information about nodes and nodal data can be found here.

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

## Accessing Nodes
The nodes stored in the ModelPart can be accessed using the Nodes parameter:

```python
model_part_nodes = fluid_model_part.Nodes
```
{: data-lang="Python"}

Having access to the nodes make iteration over all nodes very easy. For example to print all nodes in the model part:

```python
for node in fluid_model_part.Nodes:
    print(node)
```
{: data-lang="Python"}

Here is a loop over all of the nodes in a model part, which prints the ID for all of the nodes:

```python
for node in fluid_model_part.Nodes:
    print(node.Id)
```
{: data-lang="Python"}

## Node Coordinates
The coordinates can be accessed by X,Y,Z parameters of the node:

```python
node_x = node.X
node_y = node.Y
node_z = node.Z
```
{: data-lang="Python"}

Or we can extend the previous example writing also the coordinates of all the nodes in the ModelPart:

```python
for node in fluid_model_part.Nodes:
    print(node.Id, node.X, node.Y, node.Z)
```
{: data-lang="Python"}

This access is very useful in order to classify the nodes due to their position. For example we can extend the previous loop to write node information exclusively on the nodes with positive X

```python
for node in fluid_model_part.Nodes:
    if(node.X > 0.0): # Printing the ID of all of the nodes with positive X
        print(node.Id, node.X, node.Y)
```
{: data-lang="Python"}

## Nodal Data
The Python interface provides full access to the nodal database. The access to the historical variables is given by `GetSolutionStepValue` and `SetSolutionStepValue` passing the variable you want:

```python
node_velocity = node.GetSolutionStepValue(VELOCITY) # node's velocity at the current time step
```
{: data-lang="Python"}

We can write the velocities of all the nodes:

```python
for node in fluid_model_part.Nodes:
    node_velocity = node.GetSolutionStepValue(VELOCITY) # node's velocity at the current time step
    print(node_velocity)
```
{: data-lang="Python"}

you can also get a value for n time step ago, where n is the buffer size:

```python
node_previous_velocity = node.GetSolutionStepValue(VELOCITY, 1) # node's velocity at 1 time step ago 
node_earlier_velocity = node.GetSolutionStepValue(VELOCITY, 2) # node's velocity at 2 time step ago
```
{: data-lang="Python"}
 
For getting the previous time step velocities of all the nodes:  

```python
for node in fluid_model_part.Nodes:
   print(node.GetSolutionStepValue(VELOCITY, 1)) # node's velocity at 1 time step ago
```
{: data-lang="Python"}

To set the historical value for a variable in a node we can use the `SetSolutionStepValue`. To make an example 
let's assume that we want to set the variable `TEMPERATURE` to the value of 100.0 on the nodes in our `ModelPart`. This is obtained immediately by typing

```python
for node in fluid_model_part.Nodes:
    node.SetSolutionStepValue(TEMPERATURE,0,100.0)
```
{: data-lang="Python"}

The command above should be interpreted as: for the node pointed by iterator "it" assign to the variable `TEMPERATURE` at the current step (the current step is identified by 0) the value of 100.0.

**Next** [ModelPart Elements and Conditions](Elems_And_Conds)<br>
**Prev** [Writing Output File](../Input_Output_and_Visualization/Writing_output)