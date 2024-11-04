---
title: Kratos Data Management
keywords: 
tags: [Kratos-data-data-management]
sidebar: kratos_for_users
summary: 
---

# 1. Introduction

In the past section we have seen how to manipulate our geometry and perform several actions such as creating nodes, elements and writing or reading them.
If we go back to section [1.2 Creating the database]() there was one action that needed to be done before performing this operation and that was adding our variables.

In this section we will cover in detail this process and show how Kratos handles the databases and how to query the info from several structures.

# 2. Variables

When speaking about data in kratos we typically reffer to variables and we normally distigush three different types of data that one may use in a simulation.

## 2.1 Regular Variables

We understand by regular variables the ones offered by the programing language. In this case anything provided by python or C++ such an integer, a double or a string. 

These variables are used to represent generic properties that do not represent a physical magnitude or are common enough that are usefull to have in a generic way. 

For example, the coordinates of a node are a triplet of doubles that represents the physical position of that node in the space:

We can print the information about one of the nodes that we created in the previous section here:

```python
print(f"Node Info | ID: {node_3.Id}, X: {node_3.X}")
```

## 2.2. Kratos Variables

On the contraty, we understand by a Kratos variables a space that the Kratos frameworks provides to you in which you can store a value and are tied to a node, element, condition or similar. While you can create your own, there is a set of predefined variables that you can find in X and Y.

The main advantage of this variables is that can be mixed and shared among all components and applications in Kratos and its one of the pillars that allow the multiphyscis capabilites of Kratos.

In other words, this means that if your Strcutrual Solver modifies the value of the pressure in one node, everyone in Kratos will be able to account for that change, without explictly needing to communicate it.

Kratos variables are divided into two main categories:

### 2.2.1. Non Historical Variables: 

The non historical variables hold information for information that typically does not change over time, or that if it does, does it in a way that renders the old values obsolete.

These variables can be used at any point in time without needing to tell the Model which ones are needed before the simulation starts. typicaly in kratos `Nodes`, `Elements`. `Conditions`, `Properties` and `ProcessInfo` can have non-historical variables.

Non Historical Variables can be interacted with their own providded `GetValue` and `SetValue` functions. For example, setting and getting a value of PRESSURE of a node:

```python
import KratosMultiphysics

# ...

node_3.SetValue(KMP.PRESSURE, 1.5)

print(node.GetValue(KratosMultiphysics.PRESSURE))
``` 

will produce:

```
> 1.5
```
        

### 2.2.2. Historical Variables: 
Historical variables on the other hand hold not only information about the their current value but also about the value they had in previous time steps of the simulation. The number of historical values that different historical variables have is the same and can be set with the `SetBufferSize` for every modelpart.

For performance reasons the list of historical variables needs to be known before the simulation starts (and for that matter, before reading a mdpa) and cannot change. The historical variables are currently only allowed in the `nodes`.

As with the non-historical variables, there are a set of functions provided to interact with them: 
- **`SetBufferSize`**: To specify the size of the buffer (by default 1). Note that this needs to be set inside the modelpart, and no for individual entities. 

    ```Python
    modelpart = Model.CreateModelPart("Test")
    modelpart.SetBufferSize(3) # Indicates that I want to keep the value for the last three time steps
    ```

- **`AddNodalSolutionStepVariable`**: To inform the `ModelPart` that a variable will need to have an historical database. Also note that this is done at the `ModelPart` level.

    ```Python
    modelpart = Model.CreateModelPart("Test")
    modelpart.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE) # Indicates that I want to use PRESSURE as historical variable
    ```

- **`SetSolutionStepValue`**: Sets a value to the entity containing the variable. If not specified it will do it by adding it to the current time step, but that can be changed specifying how many times steps back you want to set the value:

    ```Python
    node_3.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 2.0)    # Sets PRESSURE to 2 in the current time step
    node_3.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 2.0, 0) # Sets PRESSURE to 3 in the time step 0 (current time step)
    node_3.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 4.0, 1) # Sets PRESSURE to 4 in the time step 1 (previous time step)
    ```

- **`GetSolutionStepValue`**: Gets the variable in the current entity. If not specified it will do it by adding it to the current time step, but that can be changed specifying how many times steps back you want to set the value:

    ```Python
    print(node_3.GetSolutionStepValue(KratosMultiphysics.PRESSURE))    # Gets PRESSURE in the current time step
    print(node_3.GetSolutionStepValue(KratosMultiphysics.PRESSURE, 0)) # Gets PRESSURE in the time step 0 (current time step)
    print(node_3.GetSolutionStepValue(KratosMultiphysics.PRESSURE, 1)) # Gets PRESSURE in the time step 1 (previous time step)
    ```

    will produce:

    ```Python
    > 2.0
    > 2.0
    > 4.0
    ```

# 3 Model and ModelPart structures

We are already very familiar with the modelpart and we have seen how to read and add entities to it, but we still have to cover a very important part which is how to access its data effectively

## 3.1 Model

As we have explained the model, modelpart and submodelparts are structure hierarchicly and the root of the is the Model. There three main actions that one may want to do with a model:

- **Create a ModelPart**: This is a very basic operation that we have already seen several times during the course. An interesting point to add is that you may directly create a tree structure of modelparts and submodelparts by using a name spearated by `.`:

    ```Python
    import KratosMultiphysics

    model = KratosMultiphysics.Model()
    root  = model.CreateModelPart("Root")
    child = model.CreateModelPart("Root.Child")
    leaf  = model.CreateModelPart("Root.Child.Leaf")

    print(root)
    ```

    will produce:

    ```
    -Root- model part
    Buffer Size : 1
    Number of tables : 0
    Number of sub model parts : 1
    Current solution step index : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0

    -Child- model part
        Number of tables : 0
        Number of sub model parts : 1

        Number of Geometries  : 0
        Mesh 0 :
            Number of Nodes       : 0
            Number of Properties  : 0
            Number of Elements    : 0
            Number of Conditions  : 0
            Number of Constraints : 0
        -Leaf- model part
            Number of tables : 0
            Number of sub model parts : 0

            Number of Geometries  : 0
            Mesh 0 :
                Number of Nodes       : 0
                Number of Properties  : 0
                Number of Elements    : 0
                Number of Conditions  : 0
                Number of Constraints : 0
    ```


- **Query a ModelPart**: One may want to check the status of a modelpart or get its reference. We can do this with the `HasModelPart` and `GetModelPart` functions. We can even obtain the list of all abailable modelparts using `GetModelPartNames`:

    ```Python
    import KratosMultiphysics

    model = KratosMultiphysics.Model()
    model.CreateModelPart("Root")
    model.CreateModelPart("Root.Child")
    model.CreateModelPart("Root.Child.Leaf")

    # Get the list of existing modelparts
    print(model.GetModelPartNames())

    # Query for the existance of various ModelParts
    print(model.HasModelPart("Child"))
    print(model.HasModelPart("Root.Child"))

    # Retrieve a given ModelPart
    print(model.GetModelPart("Root.Child.Leaf"))
    ```

    will produce:

    ```
    ['Root', 'Root.Child', 'Root.Child.Leaf'] # List of ModelParts

    False   # "Child"
    True    # "Root.Child"

    -Leaf- model part
    Number of tables : 0
    Number of sub model parts : 0

    Number of Geometries  : 0
    Mesh 0 :
        Number of Nodes       : 0
        Number of Properties  : 0
        Number of Elements    : 0
        Number of Conditions  : 0
        Number of Constraints : 0
    ```

- **Delete a ModelPart**: Finally after a modelpart becomes useless, we can delete it to free space with `DeleteModelPart`:

    ```Python
    import KratosMultiphysics

    model = KratosMultiphysics.Model()
    model.CreateModelPart("Root")
    model.CreateModelPart("Root.Child")
    model.CreateModelPart("Root.Child.Leaf")

    # Delete ModelPart
    if (model.HasModelPart("Root"))
        model.DeleteModelPart("Root")

    # Query for the existance of various ModelParts
    print(model.HasModelPart("Root"))
    ```

    will produce:

    ```
    False
    ```

## 3.2 ModelPart

Modepart becomes a more interesting data structre. We have seen in the previous sections some of the very basic operations needed to have a insight of the structure, like adding nodes and elements, or preparing the database as we commented in the data section, but much more can be done. In this section we will go throug some of the typical operations you will see in Kratos scripts:

## 3.2.1 Creating Entities

We have already seen this, but for the sake of completness we will present it again.

Given a ModelPart:

```Python
import KratosMultiphysics

model = KratosMultiphysics.Model()
model_part = model.CreateModelPart("ModelPart")

model_part.AddNodalSolutionStepVariable(KratosMultiphyscis.PRESSURE)
model_part.AddNodalSolutionStepVariable(KratosMultiphyscis.DISPLACEMENT)
model_part.AddNodalSolutionStepVariable(KratosMultiphyscis.REACTION)
model_part.SetBufferSize(2)
```

One can create nodes:

```python
node_1 = model_part.CreateNewNode(1, 0.00, 0.00, 0.00)
node_2 = model_part.CreateNewNode(2, 1.00, 0.00, 0.00)
node_3 = model_part.CreateNewNode(3, 1.00, 1.00, 0.00)
node_4 = model_part.CreateNewNode(4, 1.00, 0.00, 0.00)
```

Elements

```python
element_1 = model_part.CreateNewElement("Element2D3N", 1, [1,3,2], model_part.GetProperties(1))
element_2 = model_part.CreateNewElement("Element2D3N", 2, [1,4,3], model_part.GetProperties(1))
```

Conditions:

```python
condition_1 = model_part.CreateNewElement("LineCondition2D2N", 1, [1,2], model_part.GetProperties(1))
```

And master slave constrains:

```python
msc_1 = model.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, 
    node_1, KratosMultiphysics.DISPLACEMENT_X, 
    node_2, KratosMultiphysics.DISPLACEMENT_X, 
    1.0, 0
)
```

## 3.2.2 Accessing Entities

Once a modelpart contains entities, we can be interested in accessing them for any given reason. The prefered way to access entities is trough iterators. For example:

```Python
for element in model_part.Elements:
    print(element.Id)
```

will produce:

```
```

As you can se, we have accessed all elements, and for performance resons we cannot assume that the iterator order will be the same as the natural order in which they were stored in the modelpart.

Similarly, nodes can be queried the same way:

```Python
for node in model_part.Nodes:
    print(node.Id)
```

will produce:

```
```

One may be interested in accessing entities belonging to a particular area of the ModelPart or which share a common given property. If your remember submodelparts, they are particulary suited for that need. Lets create a submodelpart, add a set of nodes and iterate over them:

```Python
# You can chose any of the following lines to create the submodelpart
sub_model_part = model.CreateModelPart("ModelPart.SubModelPart")
# sub_model_part = model_part.CreateSubModelPart("SubModelPart")

# Add the nodes
sub_model_part.AddNodes([1,3])

# Print the id's of the SubModelPart Nodes
for node in sub_model_part.Nodes:
    print(node.Id)
```

will produce:

```
1
3
```

⚠️  Of course, we may have interest in accessing a particular individual entity. We can access it by id using the `[]` operator, but bear in mind that this is a very constly operation and should be avoided if possible.

```Python
my_node = model_part.Nodes[2]

print(my_node)
```

will produce:

```
Node #2 :  (1, 0, 0)
```