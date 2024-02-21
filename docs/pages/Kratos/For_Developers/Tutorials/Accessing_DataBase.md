---
title: How to Access DataBase
keywords: 
tags: [How-to-Access-DataBase.md]
sidebar: kratos_for_developers
summary: 
---

# Content
* [Overview][overview]
* [Nodal database][nodal]
   * [Historial Database][historical]
   * [Non-Historial Database][nonhistorical]
* [Nodal database][nodal]
   * [Non-Historial Database][nonhistoricalelem]

[overview]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Access-DataBase#overview
[nodal]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Access-DataBase#nodal-database
[historical]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Access-DataBase#historical-database
[nonhistorical]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Access-DataBase#non-historical-database
[nonhistoricalelem]: https://github.com/KratosMultiphysics/Kratos/wiki/How-to-Access-DataBase#elemental-database

# Overview

Accessing to the database is one of the most frequent operations in a Finite Element context.

In a Finite Element context we often need to store information on the `Nodes` and on the `Elements`. This data may take different forms (scalars, vectors, matrices ... ) and needs to be widely accessible inside a finite element program. In the _Kratos_, this data is stored directly in the object from which it will be accessed, meaning that the nodes will contain the nodal data, and the elements the elemental data.

The _Kratos_ provides both a nodal and elemental `Variable Database`, which allows to store in a flexible way `VARIABLES` of different type. A detailed description of the design concepts together with a study of the performance issues can be found in [Pooyan's Thesis](https://futur.upc.edu/3239147) . The aim of current "How To" is just to provide a brief introduction on the usage of the `C++` interface. 

# Nodal database

The `Nodes` constitute one of the building blocks of the _Kratos_ interface and provide the richest database structure. Conceptually each nodes allocates a portion of memory sufficient to contain a database with all of its data. In the practice the memory needed is defined by the user which at the beginning of the program has to select and add to the model part the list of variables that will be needed during the analysis.

As a concept we should also note that for some variables a "_history_" (the value of this variables in the past steps) will be needed, while for others only the current value will be necessary. In order to make this possible the _Kratos_ provides two different implementations for the "_Historical Database_" and for the "_Non-Historical one_".

The two databases are independent from each other and are **not synchronized**. 

## Historical database

The "_Historical Database_" stores the present value of a variable and the value it took in the previews steps. The Number of steps stored depends on the size of `Buffer`. This can be modified by the function `SetBufferSize()` while the present value can be obtained by `GetBufferSize()`. To make an example :

```cpp
model_part.SetBufferSize(3)
unsigned int model_part.GetBufferSize(); ##now outputs 3
```

Implies storing the current step plus the 2 last steps.

The list of the variables to be stored and the buffer size needs to be known in order to use the "_Historical Database_". This is obtained by providing at the beginning of the analysis the list of all the variables that will be involved in the calculation. to make an example if a given solver will need the variables `TEMPERATURE` and `VELOCITY`, the user should provide **before creating the list of nodes** the two commands.

```cpp
model_part.AddNodalSolutionStepVariable(TEMPERATURE)
model_part.AddNodalSolutionStepVariable(VELOCITY)
```

In the practice each solver provides the list of its variables through the function `AddVariables(...)` (examples of this can be found in any of the `tests`). Note that trying to access for example to `PRESSURE` will lead to an error as the database did not allocate any memory to hold the variable `PRESSURE`

Once this preliminary operations are performed, the memory for the _Kratos_ is allocated and can be used efficiently. Two Functions exist to access to the solution step data: 

* `GetSolutionStepData `
*  `FastGetSolutionStepData` 

Their functionality is identical but they differ in the checks that are performed internally in the _Release_ version (They are identical in _Debug_). The `Fast...` operator, does not perform any check, meaning that if, in the example before, we try to access to `PRESSURE` (not in the variable list) we will get a **segmentation fault** without any further information.

The syntax is as follows: 

```cpp
  Node::iterator inode = model_part.NodesBegin(); //to make example let's take the first node
  
  //here  we get REFERENCES to the database!!
  array_1d<double,3>& vel = inode->FastGetSolutionStepValue(VELOCITY); //vel here has the value of velocity at the current step
  const array_1d<double,3>& oldvel = inode->FastGetSolutionStepValue(VELOCITY,1); //vel 1 step in the past
  const array_1d<double,3>& veryoldvel = inode->FastGetSolutionStepValue(VELOCITY,2); //vel 2 steps in the past
  //
  //predicting the vel as     vel = 2*oldvel - veryoldvel;
  noalias(vel) = 2.0 * oldvel;
  noalias(vel) -= veryoldvel;
  //NOTE: "vel" is a REFERENCE to the database => changing "vel" we change the database
```

##  Non-Historical Database 

On the side of the "_historical_" database it is also possible to store variables which do not change as the time step evolves. As an important difference with respect to the historical database, the list of variables included nor the "_size_" of the objects does not need to be prescribed at the beginning of the program and **can be changed dynamically during the run**.

**ANY** element type can be included in this sort of database. As an example, a part for the standard types double `array_1d<double,3>` `Matrix` we can use user defined types. Two useful examples of this are:

* `NEIGHBOUR_NODES`
* `NEIGHBOUR_ELEMENTS` 

The usage of such entities is described elsewhere (..) the important point is here that an array of pointers can be stored without any problems together with the standard datatypes. The implementation details and the discussion of the implementation can be found in [Pooyan's Thesis](https://scholar.google.es/citations?view_op=view_citation&hl=es&user=Z3lsvS8AAAAJ&citation_for_view=Z3lsvS8AAAAJ:hqOjcs7Dif8C).

The access is through the functions `GetValue()` and `SetValue()`. To make an example: 

```cpp
  double temperature = inode->GetValue(TEMPERATURE)
```
or
```cpp
  double aaa = 10.0
  inode->SetValue(TEMPERATURE,aaa);
```

# Elemental Database

The Elemental database is identical to the "_non-historical_" database (it shares exactly the same implementation). Non historical database is implemented for elements and conditions.

The user should nevertheless note that the user can implement a user-defined storage through the elemental interface (see for example the function `Calculate` or `CalculateOnIntegrationPoints`) 