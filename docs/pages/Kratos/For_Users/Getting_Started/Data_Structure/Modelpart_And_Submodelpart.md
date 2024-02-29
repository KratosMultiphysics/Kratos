---
title: ModelPart and SubModelPart
keywords: 
tags: [Python Script Tutorial ModelPart SubModelPart]
sidebar: kratos_for_users
summary: 
---

In the previous part of the tutorial, we already saw how the `ModelPart` is the object containing `Element`, `Conditions`, `Nodes` and `Properties`.

A fundamental feature is that it can also hierarchically contain **"SubModelParts"** intended as other `ModelParts` which belong to the same parent. This relation can be repeated recursively, so that each "root" `ModelPart` can actually own a tree of SubModelParts.

A quite extensive testing can be found [here](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/tests/test_model_part.py)

However let's try to make an example to explain this better.

```python        
# Create a ModelPart root
current_model = Model()
model_part = current_model.CreateModelPart("Main")

# Now create a SubModelPart
model_part.CreateSubModelPart("Inlets")

# Let's output what is there:
print(model_part)
```
{: data-lang="Python"}

the output is:

```console      
-Main- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 1
    Current solution step index : 0

    Mesh 0 :
        Number of Nodes      : 0
        Number of Properties : 0
        Number of Elements   : 0
        Number of Conditions : 0

    -Inlets- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes      : 0
            Number of Properties : 0
            Number of Elements   : 0
            Number of Conditions : 0
```
{: data-lang="Python Output"}

We could now verify if a given submodelpart exists, or how many SubModelParts exist as

```python        
model_part.HasSubModelPart("Inlets") #returns True
model_part.NumberOfSubModelParts() #returns 1
model_part.GetSubModelPart("Inlets").Name #returns the name --> Inlets
```
{: data-lang="Python"}

Let's now create some other SubModelParts

```python
# On the first level
model_part.CreateSubModelPart("Temp")
model_part.CreateSubModelPart("Outlet")

# On the second level --> "sub-sub modelparts"
sub_model_part_1 = model_part.GetSubModelPart("Inlets")
sub_model_part_1.CreateSubModelPart("Inlet1")
sub_model_part_1.CreateSubModelPart("Inlet2")

# Output
print(model_part)
```
{: data-lang="Python"}

to give

```console
-Main- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 3
    Current solution step index : 0

    Mesh 0 :
        Number of Nodes      : 0
        Number of Properties : 0
        Number of Elements   : 0
        Number of Conditions : 0

    -Outlet- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes      : 0
            Number of Properties : 0
            Number of Elements   : 0
            Number of Conditions : 0
    -Temp- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes      : 0
            Number of Properties : 0
            Number of Elements   : 0
            Number of Conditions : 0
    -Inlets- model part
        Number of tables : 0
        Number of sub model parts : 2

        Mesh 0 :
            Number of Nodes      : 0
            Number of Properties : 0
            Number of Elements   : 0
            Number of Conditions : 0
        -Inlet2- model part
            Number of tables : 0
            Number of sub model parts : 0

            Mesh 0 :
                Number of Nodes      : 0
                Number of Properties : 0
                Number of Elements   : 0
                Number of Conditions : 0
        -Inlet1- model part
            Number of tables : 0
            Number of sub model parts : 0

            Mesh 0 :
                Number of Nodes      : 0
                Number of Properties : 0
                Number of Elements   : 0
                Number of Conditions : 0
```
{: data-lang="Python Output"}

Each `ModelPart` is only directly aware of its first level siblings. that is

```python
model_part.HasSubModelPart("Inlet1") #--> returns False
```
{: data-lang="Python"}

However

```python
model_part.GetSubModelPart("Inlets").HasSubModelPart("Inlet1") #--> returns True
```
{: data-lang="Python"}

Eventually we can loop on all the SubModelParts of a given submodelpart by doing

```python
for part in model_part.SubModelParts:
    print(part.Name)
```
{: data-lang="Python"}

Which would output

```Console
Outlet
Inlets
Temp
```
{: data-lang="Python Output"}

## Data Ownership
The parent-son relation is such that **anything that belongs to a given SubModelPart also belongs to the parent ModelPart**.

This implies that the ultimate "owner" of any `Node`, `Element`, etc, will be the "root" `ModelPart`. The consistency of the tree is ensured by the `ModelPart` **API**, which provides the tools needed for creating or removing anything any of the contained objects.

As usually let's try to explain this by examples. We can create a node by doing

```python
model_part.CreateNewNode(1, 1.00,0.00,0.00)
```
{: data-lang="Python"}

If we Try to create a new node with the same Id and different coordinates we would get an error

```python
model_part.CreateNewNode(1, 0.00,0.00,0.00)  # Here an error is thrown
```
{: data-lang="Python"}

However if we try to create a node with the same coordinates twice nothing is actually done (and no error is thrown)

```python
model_part.CreateNewNode(1, 1.00,0.00,0.00) 
print(model_part.NumberOfNodes()) # This still returns 1!!
```
{: data-lang="Python"}

# We can now access the node as needed, for example

```python
model_part.GetNode(1).Id #gives 1
model_part.GetNode(1,0).X #gives 1.0
```
{: data-lang="Python"}

`Nodes` can be created in every order

```python
model_part.CreateNewNode(2000, 2.00,0.00,0.00)
model_part.CreateNewNode(2, 2.00,0.00,0.00)
```
{: data-lang="Python"}

We could then loop over all the nodes

```python
for node in model_part.Nodes:
    print(node.Id, node.X, node.Y, node.Z)
```
{: data-lang="Python"}

Or eventually remove nodes one by one by doing

```python
model_part.RemoveNode(2000)
```
{: data-lang="Python"}

Let's now see what happens if we add a node to a submodelpart.
here the node will be both in root `ModelPart` and in Inlets, but for example not in "Temp" or "Outlet"

```python
inlets_model_part = model_part.GetSubModelPart("Inlets")
inlets_model_part.CreateNewNode(3, 3.00,0.00,0.00)
```
{: data-lang="Python"}

If we add to a sub-sub modelpart, it will belong to the root and the parent, but for example not to Inlet1

```python
inlet2_model_part = inlets_model_part.GetSubModelPart("Inlet2")
inlet2_model_part.CreateNewNode(4, 4.00,0.00,0.00)
```
{: data-lang="Python"}

Multiple nodes can be removed at once (and from all levels) by flagging them

```python
for node in model_part.Nodes:
    if(node.Id < 3):
        node.Set(TO_ERASE,True)

model_part.RemoveNodesFromAllLevels(TO_ERASE)
```
{: data-lang="Python"}

One could call simply the function `RemoveNodes` and remove them from the current level down.

# Creating Elements and Conditions
`Elements` and `Conditions` can be created from the python interface by providing their connectivity as well as the `Properties` to be employed in the creation. The string to be provided is the name by which the element is registered in **Kratos*.

```python
model_part.AddProperties(Properties(1))
model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])
```
{: data-lang="Python"}

an error is thrown if i try to create an element with the same Id

```python
model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1]) # Here an error is thrown
```
{: data-lang="Python"}

An identical interface is provided for Conditions, as well as functions equivalent to the nodes for removing from one level or from all the levels.
