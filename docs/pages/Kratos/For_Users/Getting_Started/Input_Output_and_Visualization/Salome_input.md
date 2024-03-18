---
title: Reading ModelPart From Salome
keywords: 
tags: [Python Script Tutorial Reading ModelPart Input File]
sidebar: kratos_for_users
summary: 
---

The `ModelPart` represents an arbitrary part of the `Model` to be simulated and stores the mesh and additional data for it. Most of the *Kratos* routines take a `ModelPart` as their argument. So always is necessary to create and fill a `ModelPart`. 

In order to load a file from Salome, we will make use of the Med Application to convert it from its native format to Kratos' mdpa. 

In this example we will show how to read a Salome model and use it for a structural analysis:

# Setup
First of all we need to create a python file with following code to import the *Kratos*:

```python
# Kratos and MedApplication
import KratosMultiphysics as KMP
import KratosMultiphysics.MedApplication import MED

# StructuralMechanicAnalysis from StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
```
{: data-lang="Python"}

# Creating a ModelPart

First, we will read our parameters:

```python
# Load Parameters
with open("ProjectParameters.json",'r') as parameter_file:
    parameters = KMP.Parameters(parameter_file.read())
```
{: data-lang="Python"}

To create a `ModelPart`, one has to create a new `Model` first, and then to call its constructor passing the `ModelPart``s name as its argument:

```python
# Create the model and modelpart
model = KMP.Model()
model_part = model.CreateModelPart("Structure")
```
{: data-lang="Python"}

After that, we must create our simulation. Since this is a structural analysis, we will use the StructuralMechanicsAnalysis

```python
# Create the Analysis Stage
simulation = StructuralMechanicsAnalysis(model,parameters)
```
{: data-lang="Python"}

Then you must read the `.med` file:

```python
# Read the med file
KratosMED.MedModelPartIO("plate_with_hole.med", KratosMED.MedModelPartIO.READ).ReadModelPart(model_part)
```
{: data-lang="Python"}

Once the file is processed, the geometry from the `.med` file will be loaded into the kratos modelpart. In order to give the correct elements and condition we will use the `CreateEntitiesFromGeometriesModeler`:

```python
# apply the elements and conditions
params = KMP.Parameters("""{
    "elements_list" : [
        {
            "model_part_name" : "Structure.solid",
            "element_name" : "SmallDisplacementElement2D3N"
        }
    ],
    "conditions_list" : [
        {
            "model_part_name" : "Structure.load",
            "condition_name" : "LineLoadCondition2D2N"
        {
    ]
}""")

modeler = KMP.CreateEntitiesFromGeometriesModeler(model, params)
modeler.SetupModelPart()
```
{: data-lang="Python"}

Currently the modeler does not assign any property to our conditions and elements, so we will assign one:

```python
# Assign a fresh properties container to the model
properties = model_part.CreateNewProperties(1)
for cond in model_part.Conditions:
    cond.Properties = properties

for elem in model_part.Elements:
    elem.Properties = properties
```
{: data-lang="Python"}

We can print the modelpart to check that everything worked correctly:

```python
>>> print(model_part)
-Structure- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 3
    Current solution step index : 0

    Number of Geometries  : 3758
    Mesh 0 :
        Number of Nodes       : 1879
        Number of Properties  : 1
        Number of Elements    : 3447
        Number of Conditions  : 160
        Number of Constraints : 0

    -fixture- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 104
        Mesh 0 :
            Number of Nodes       : 105
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 0
            Number of Constraints : 0
    -load- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 160
        Mesh 0 :
            Number of Nodes       : 161
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 160
            Number of Constraints : 0
    -solid- model part
        Number of tables : 0
        Number of sub model parts : 0

        Number of Geometries  : 3447
        Mesh 0 :
            Number of Nodes       : 1879
            Number of Properties  : 0
            Number of Elements    : 3447
            Number of Conditions  : 0
            Number of Constraints : 0
```
{: data-lang="Python Output"}
