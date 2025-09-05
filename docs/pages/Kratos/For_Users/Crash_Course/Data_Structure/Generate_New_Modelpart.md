---
title: Generating New Modelparts
keywords: 
tags: [Python Script Tutorial Generating New Modelpart]
sidebar: kratos_for_users
summary: 
---

`ModelParts` are the essential data structure to hold **FEM** objects in Kratos.

Since the "physics" of a problem is provided by the `Element` and `Condition` which implement it, in order to describe a new physical problem one should provide a new `ModelPart` describing the connectivity but also made of the relevant element technology.

In *Kratos*, essentially all of the FEM objects (`Nodes`, `Elements`, `Conditions`, `Properties`, `ProcessInfo`) are managed by shared pointers, and can hence have multiple owners. This mean in the practice that a given, say, `Node` may belong at the same time to multiple ModelParts. 

A typical use case is that in which multiple physical problems should be solved on a single discretization of the problem. This is the case for example of the Fluid-Thermal problem that we will describe at the end of the tutorial. *Kratos* has a special `Modeler` named `ConnectivityPreserveModeler` that fills a modelpart by preserving the same connectivity as in the source model part while changing the element technology.

The resultant modelpart is a "free standing" root modelpart, completely independent of the original one. Nevertheless it shares:
- Pointers to the same nodes
- Same `ProcessInfo`
- Same `Properties`
- Same `Tables`
as in the source modelpart.

The difference thus lays in the Element type being employed, which substitute in the new modelpart the element of the old modelpart with a new element implementing the desired physics.

The usage of the modeler is as follows:

```python
modeler = KratosMultiphysics.ConnectivityPreserveModeler()
modeler.GenerateModelPart(self.main_model_part, self.thermal_model_part, "Element2D3N", "Condition2D2N")
```
{: data-lang="Python"}





