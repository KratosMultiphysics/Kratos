---
title: Import results
keywords: 
tags: [Import_results.md]
sidebar: shape_optimization_application
summary: 
---

This section of the ribbon is used to import the results produced by KratosMultiphysics optmization procedure for the specified problem.

This imports two types of results as shown in figure 1.

<p align="center">
    <img src="images/imported_results.png" alt="Imported results"/>
</p>
<p align="center">Figure 1: Imported results</p>

First import is always the optimization result import where all the updated designs and their respective nodal values are imported. Following table briefly describe the functionality of each variable.

|Variable name|Description|
|-------------|-----------|
|SHAPE_CHANGE| The change of the design surface computed considering objectives and constraints for all the iterations. This represents all the shape updates occured up until selected iteration.|
|MESH_CHANGE| The change of the mesh nodal position due to the use of mesh motion solver when changing the design surface nodes to obtain the optimized designs|
|DF1DX| Objective gradient w.r.t. nodal coordinates (refer $$\frac{df}{d\underline{s}}$$ in [objective values](../Technologies/Objectives.html))|
|DC1DX| Constraint gradient w.r.t. nodal coordinates (refer $$\frac{dg}{d\underline{s}}$$ in [constraint values](../Technologies/Constraints.html))|
|DF1DX_MAPPED| Vertex morphed objective gradients w.r.t. nodal coordinates (refer $$\left(\frac{df}{d\underline{s}}\right)_{i, morphed}$$ in [vertex morphing](../Technologies/Vertex_morphing.html#algorithm))|
|DC1DX_MAPPED| Vertex morphed constraint gradients w.r.t. nodal coordinates (refer $$\left(\frac{df}{d\underline{s}}\right)_{i, morphed}$$ in [vertex morphing](../Technologies/Vertex_morphing.html#algorithm))|
|SHAPE_UPDATE| Current shape update computed by selected algorithm (refer $$U$$ in [steepest descent](../Technologies/Algorithms/steepest_descent.html) or [penalized projection](../Technologies/Algorithms/penalized_projection.html))|

Next few imports corresponds to primal result imports computed in each design iteration. The always have the "primal" keyword in the import names followed by the design iteration number.

|Variable name|Description|
|-------------|-----------|
|DISPLACEMENT| The solved displacement for the primal problem with the selected design iterations design surface|
|POINT_LOAD| Loads applied on the primal problem in the selected design iteration|
|VON_MISES_STRESS| Computed Von Mises stress for each element|
|REACTION| Residuals of each nodal equaton. Corresponds to reaction forces at fixed nodes|
