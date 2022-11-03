---
title: Design variables
keywords: 
tags: [Design_variables.md]
sidebar: shape_optimization_application
summary: 
---

This section defines the freedom of the optimization problem to achieve the optimum design surface optimiying the objectives and sstosfying the given constraints.

## Design variables dialog

<p align="center">
    <img src="images/design_variables.png" alt="Design variables dialog box"/>
</p>
<p align="center">Figure 1: Design variables dialog box</p>

|Option|Description|
|------|-----------|
|Selection method| Selection method to identify allowed surfaces to change to obtain optimum design surface. Three types are supported. **Manual**: Selects the design surfaces manually. **Element face group**: Selects surfaces of a selected element group. **Auto: All_Faces**: Selects all the outer faces of the model|
|Filter shape| Filter function type (refer [vertex morphing filters](../Technologies/Vertex_morphing.html#vertex-morphing-options)) |
|Filter radius| Filter radius (i.e. $$r$$) in vertex morphing methodology (refer [vertex morphing filter radius](../Technologies/Vertex_morphing.html#effect-of-filter-radii)) |
|Use damping| Enable damped regions |
|Damping shape| The damping function to be used. (refer [damping functions](../Technologies/Damping.html#damping-functions)) |
|Damping radius| Radius of damping effect around a given single node in the design surface |
|Damping directions| Selection of damping directions. Either "X", "Y" or "Z" direction of the sensitivities can be damped.|
