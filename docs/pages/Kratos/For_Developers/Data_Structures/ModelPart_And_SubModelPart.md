---
title: ModelPart and SubModelPart
keywords: 
tags: [ModelPart SubModelPart]
sidebar: kratos_for_developers
summary: 
---

# Create and edit ModelParts

Here we introduce the basics to create and access to model parts and submodelparts.

## Starting

First of all we need to create a python file with following code to import the `KratosMultiphysics`.

```python
import KratosMultiphysics
```

## Create Model

Then we need to create the `Model`, which will be the resposible to manage the different `ModelPart`s that we will create.

```python
this_model = KratosMultiphysics.Model()
```

## Create ModelPart

Now we can create a `ModelPart`. The `ModelPart` is the object containing `Element`, `Conditions`, `Nodes` and `Properties`. For now we create the *Main* model part, which will store the successive submodelparts.

```python
main_model_part = this_model.CreateModelPart("Main")
```

We can create a new model part with a certain *buffer size* using the following (we need to delete the model part to avoid errors):

```python
this_model.DeleteModelPart("Main")
main_model_part = this_model.CreateModelPart("Main", 2)
```

We can now execute different operations with the `Model`:

```python
print(this_model.HasModelPart("Main")) # It will return True
print(this_model.GetModelPartNames()) # It will return ['Main']
main_model_part_again = this_model.GetModelPart("Main") # Getting again
```

Let's output what is there:

```python
print(main_model_part)
```

```console
-Main- model part
    Buffer Size : 2
    Number of tables : 0
    Number of sub model parts : 0
    Current solution step index : 0

    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0
```

Some other operations we can do are:

```python
print(main_model_part.NumberOfNodes()) # It will return 0
print(main_model_part.NumberOfElements()) # It will return 0
print(main_model_part.NumberOfConditions()) # It will return 0
print(main_model_part.NumberOfMasterSlaveConstraints()) # It will return 0
print(main_model_part.NumberOfProperties()) # It will return 0
print(main_model_part.NumberOfMeshes()) # It will return 1
print(main_model_part.GetBufferSize()) # It will return 2
main_model_part.SetBufferSize(3) # Set the buffer size to 3 instead of 2
```

## Create sub ModelPart

A fundamental feature is that it can also hierarchically contain *SubModelParts* intended as other `ModelParts` which belong to the same parent. This relation can be repeated recursively, so that each "root" `ModelPart` can actually own a tree of SubModelParts. We can create a submodelpart with the following:

```python
bc_model_part = main_model_part.CreateSubModelPart("BC")
```

Let's output what is there:

```python
print(main_model_part)
```

```console
-Main- model part
    Buffer Size : 3
    Number of tables : 0
    Number of sub model parts : 1
    Current solution step index : 0

    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0

    -BC- model part
        Number of tables : 0
        Number of sub model parts : 0

        Mesh 0 :
            Number of Nodes       : 0
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 0
            Number of Constraints : 0

```

Now we can do several operations with this:

```python
print(main_model_part.HasSubModelPart("BC")) #returns True
print(main_model_part.NumberOfSubModelParts()) #returns 1
print(main_model_part.GetSubModelPart("BC").Name) #returns the name --> BC
```

## Data Ownership

The parent-son relation is such that **anything that belongs to a given SubModelPart also belongs to the parent ModelPart**. 

This implies that the ultimate *owner* of any `Node`, `Element`, etc, will be the *root* `ModelPart`. The consistency of the tree is ensured by the `ModelPart` **API**, which provides the tools needed for creating or removing anything any of the contained objects.