---
title: Nodes and Nodal Data
keywords: 
tags: [Nodes NodalData Data]
sidebar: kratos_for_developers
summary: 
---

# Nodes and Nodal Data

In this tutorial the access to the nodes stored in a `ModelPart` and their nodal data will be described. More information about nodes and nodal data can be found here.

## Starting
First of all we need to create a python file with following code to import the *Kratos*, create a `ModelPart` and read it from input as described in the here :

```python
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication

this_model = KratosMultiphysics.Model()
structural_model_part = this_model.CreateModelPart("StructuralPart", 3)

structural_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
structural_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

structural_model_part_io = KratosMultiphysics.ModelPartIO("KratosWorkshop2019_high_rise_building_CSM")
structural_model_part_io.ReadModelPart(structural_model_part)
```

## Accessing Nodes
The nodes stored in the ModelPart can be accessed using the Nodes parameter:

```python
model_part_nodes = structural_model_part.Nodes
```

Having access to the nodes make iteration over all nodes very easy. For example to print all nodes in the model part:

```python
for node in model_part_nodes:
    print(node)
```

Here is a loop over all of the nodes in a model part, which prints the ID for all of the nodes:

```python
for node in model_part_nodes:
    print(node.Id)
```

## Node Coordinates
The coordinates can be accessed by X,Y,Z parameters of the node:

```
node_x = node.X
node_y = node.Y
node_z = node.Z
```

Or we can extend the previous example writing also the coordinates of all the nodes in the ModelPart:

```python
for node in model_part_nodes:
    print(node.Id, node.X, node.Y, node.Z)
```

This access is very useful in order to classify the nodes due to their position. For example we can extend the previous loop to write node information exclusively on the nodes with positive X

```python
for node in model_part_nodes:
    if(node.X > 0.0): # Printing the ID of all of the nodes with positive X
        print(node.Id, node.X, node.Y)
```

## Nodal Data
The Python interface provides full access to the nodal database. The access to the historical variables is given by `GetSolutionStepValue` and `SetSolutionStepValue` passing the variable you want:

```python
node_displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) # node's displacement at the current time step
```

We can write the displacements of all the nodes:

```python
for node in model_part_nodes:
    node_displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT) # node's displacement at the current time step
    print(node_displacement)
```
you can also get a value for n time step ago, where n is the buffer size:

```python
node_previous_displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1) # node's displacement at 1 time step ago 
node_earlier_displacement = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 2) # node's displacement at 2 time step ago
```
 
For getting the previous time step displacements of all the nodes: 

```python
for node in model_part_nodes:
   print(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 1)) # node's displacement at 1 time step ago
```

To set the historical value for a variable in a node we can use the `SetSolutionStepValue`. To make an example  let's assume that we want to set the variable `DISPLACEMENT_X` to the value of 1.0e-6 on the nodes in our `ModelPart`. This is obtained immediately by typing

```python
for node in model_part_nodes:
    node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0,1.0e-6)
```

To set the non-historical value for a variable in a node we can use the `SetValue`. Later this can be accessed with `GetValue`. For example:

```python
for node in model_part_nodes:
    node.SetValue(KratosMultiphysics.TEMPERATURE,100.0)
    print(node.GetValue(KratosMultiphysics.TEMPERATURE))
```