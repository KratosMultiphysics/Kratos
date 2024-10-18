---
title: Entities creation from scratch
keywords: 
tags: [Entities Creation From Scratch]
sidebar: kratos_for_users
summary: 
---

## Entities creation from scratch

Here we present the interface in order to create new entities in an empty model from *Python*.

### Creating nodes

We can create a node by doing. If we Try to create a new node with the same Id and different coordinates we would get an error.

```python
main_model_part.CreateNewNode(1, 1.0,0.0,0.0)
#main_model_part.CreateNewNode(1, 0.0,0.0,0.0)  # Here an error is thrown
```
{: data-lang="Python"}

However if we try to create a node with the same coordinates twice nothing is actually done (and no error is thrown)

```python
main_model_part.CreateNewNode(1, 1.00,0.00,0.00) 
print(main_model_part.NumberOfNodes()) # This still returns 1!!
```
{: data-lang="Python"}

We can now access the node as needed, for example:

```python
print(main_model_part.GetNode(1).Id) # Gives 1
print(main_model_part.GetNode(1,0).X) # Gives 1.0
```
{: data-lang="Python"}

`Nodes` can be created in every order

```python
main_model_part.CreateNewNode(2000, 2.00,0.00,0.00)
main_model_part.CreateNewNode(2, 2.00,0.00,0.00)
```
{: data-lang="Python"}

We could then loop over all the nodes:

```python
for node in main_model_part.Nodes:
    print(node.Id, node.X, node.Y, node.Z)
```
{: data-lang="Python"}

Or eventually remove nodes one by one by doing:

```python
main_model_part.RemoveNode(2000)
```
{: data-lang="Python"}

Let's now see what happens if we add a node to a submodelpart. Here the node will be both in root `ModelPart` and "BC", but for example not in derived submodelparts from "BC" or root `ModelPart`.

```python
bc_model_part = main_model_part.GetSubModelPart("BC")
bc_model_part.CreateNewNode(3, 3.00,0.00,0.00)
```
{: data-lang="Python"}

Multiple nodes can be removed at once (and from all levels) by flagging them:

```python
for node in main_model_part.Nodes:
    if node.Id < 3:
         node.Set(KratosMultiphysics.TO_ERASE,True)   
main_model_part.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
```
{: data-lang="Python"}

One could call simply the function `RemoveNodes` and remove them from the current level down.

### Creating Elements and Conditions
`Elements` and `Conditions` can be created from the python interface by providing their connectivity as well as the `Properties` to be employed in the creation. The string to be provided is the name by which the element is registered in *Kratos*. An error is thrown if i try to create an element with the same Id

```python
main_model_part.AddProperties(KratosMultiphysics.Properties(1))
main_model_part.CreateNewElement("Element2D3N", 1, [1,2,3], main_model_part.GetProperties()[1])
```
{: data-lang="Python"}

An identical interface is provided for Conditions, as well as functions equivalent to the nodes for removing from one level or from all the levels.

### Creating MasterSlaveConstraints

The case of the constraints, like it requires the existance of `DoF` we will create a new model part in order to avoid problems. This example will advance some methods that will be introduced in future sections, as `AddNodalSolutionStepVariable`. The current interface for python works directly with nodes, it is not possible to manipulate `DoF` manually.

```python
mp = this_model.CreateModelPart("constraint_example")

mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

mp.CreateNewNode(1, 0.00000, 0.00000, 0.00000)
mp.CreateNewNode(2, 0.00000, 1.00000, 0.00000)

KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)

mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, mp.Nodes[1], KratosMultiphysics.DISPLACEMENT_X, mp.Nodes[2], KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)
```
{: data-lang="Python"}