---
title: ModelPart Entities
keywords: 
tags: [ModelPart Entities]
sidebar: kratos_for_developers
summary: 
---

# Starting
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

# Accessing `Elements`
The elements stored in the ModelPart can be accessed using the Elements parameter:

```python
model_part_elements = structural_model_part.Elements
```
 
Iteration over all elements in a model part is very similar to the nodes. For example writing the ID elements in a model part can be done as follow:

```python
for element in model_part_elements:
    print(element.Id)
```

Additionally we can access for example the geometry of the element and ask the area of each element:

```python
for element in structural_model_part.Elements:
    print("ID", element.Id, " AREA: ", element.GetGeometry().Area())
```

# Accessing `Conditions`
Conditions parameter of model part provides access to the conditions it stores:

```python
model_part_conditions = structural_model_part.Conditions
```

Iteration over conditions is very similar to the elements. In the same way printing the ID conditions is as follow:

```python
for condition in model_part_conditions:
    print(condition.Id)
```