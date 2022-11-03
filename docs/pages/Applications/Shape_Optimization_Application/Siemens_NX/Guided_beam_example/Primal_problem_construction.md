---
title: Fixed beam primal problem construction
keywords: 
tags: [guided_example.md]
sidebar: shape_optimization_application
summary: 
---

## Introduction

This example guides the user through a creation of an primal problem in **Siemens NX** plugin to be used in an optimization problem.

## Primal problem description

Figure 1 illustrates the fixed beam configuration.

<p align="center">
    <img src="images/fixed_beam.png" alt="Fixed beam configuration"/>
</p>
<p align="center">Figure 1: Fixed beam configuration</p>

Following table illustrates the values used in this problem.

|Symbol| Value|
|------|------|
|L| $$100\,mm$$|
|T| $$10\,mm$$|
|P| $$ 100\,N$$|

## Primal problem configuration steps

### Creating a new project

First of all the geometry needs to be created. For that, you can open a new model by clicking "File" and then "New" as indicated in figure 2.

<p align="center">
    <img src="images/step_1_new.png" alt="Choosing new"/>
</p>
<p align="center">Figure 2: Choosing new</p>


Then you will be given the "New" dialog box as shown in the figure 3. The 3rd step is to select "Model" in the "Templates" section. Then a proper "Name" and "Folder" should be given in 4th and 5th steps.
<p align="center">
    <img src="images/step_2_model.png" alt="Creating a new model"/>
</p>
<p align="center">Figure 3: Creating a new model</p>

### Creating geometry

Then you can create the geometry to given dimensions. Illustration of final geometry is given in figure 4.
<p align="center">
    <img src="images/step_3_geometry.png" alt="Creating geometry"/>
</p>
<p align="center">Figure 4 Creating geometry</p>

### Creating FEM Simulation

Then "PrePost" can be enabled by following the steps illustrated in figure 5 by first clicking on "Application" and then clicking on "PrePost" on the top ribbon.
<p align="center">
    <img src="images/step_4_prepost.png" alt="Activating pre-post"/>
</p>
<p align="center">Figure 5: Activating pre-post</p>

Thereafter, a FEM simulation can be created by clicking on "New FEM and Simulaion" as shown in the figure 6.
<p align="center">
    <img src="images/step_5_femsim.png" alt="Selecting FEM Simulation"/>
</p>
<p align="center">Figure 6: Selecting FEM Simulation</p>

When you select "New FEM Simulation", following dialog box will appear as shown in figure 7. In there, "Solver" needs to be set to "Simcenter Nastran" and "Analysis Type" needs to be set to "Structural" as shown in step 10. Thereafter, select "Solution type" as "SOL 101: Linear Statics" to obtain a linear statics solution as in step 12. Then click "Ok" on the next dialog box.
<p align="center">
    <img src="images/step_6_structfem.png" alt="Creating FEM simulation"/>
</p>
<p align="center">Figure 6: Creating FEM simulation</p>

### Creating mesh

To create the mesh, first, the FEM model needs to be selected from "Simulation File View" by double clicking on the "model1_fem1" as shown in figure 7 step 13. Then "3D Tetrahedral" can be selected from the ribbon as in step 14. Then the the dialog box to the right in figure 7 appears. In there, "Object to Mesh" needs to selected with the whole body of your constructed geometry. Then tetrahedral element with 4 nodes needs to be selected as the element type in "Element Properties -> Type". This can be done by selecting "CTETRA(4)" in the options box as in step 16. Thereafter, "Element size" is put to "1mm" as shown in step 17. Onceyou press "Ok", then the mesh will be generated.
<p align="center">
    <img src="images/step_7_mesh.png" alt="Creating mesh"/>
</p>
<p align="center">Figure 7: Creating mesh</p>

The generated mesh will look like in figure 8.
<p align="center">
    <img src="images/step_8_generated_mesh.png" alt="Generated mesh"/>
</p>
<p align="center">Figure 8: Generated mesh</p>

### Applying material properties

Figure 9 illstrates steps required to apply material "Auminium 2014" to the whole solid part of the fem model.
<p align="center">
    <img src="images/step_9_materials.png" alt="Applying materials"/>
</p>
<p align="center">Figure 9: Applying materials</p>

### Applying boundary conditions

Figure 10 illustrates the steps required to apply fixed boundary condition on the two end surfaces of the beam. Right figure of figure 10 shows the illustration after applying the fixed boundary conditions to two end square surfaces.

<p align="center">
    <img src="images/step_10_bc_fixed.png" alt="Applying fixed boundary condition"/>
</p>
<p align="center">Figure 10: Applying fixed boundary condition</p>

Figure 11 illustrates the steps required to apply force boundary condition on the top surface (normal to the "Z" direction). Right figure of figure 11 shows the illustration after applying the fixed boundary condition and force boundary condition.

<p align="center">
    <img src="images/step_11_bc_load.png" alt="Applying force boundary condition"/>
</p>
<p align="center">Figure 11: Applying force boundary condition</p>

### Solving primal problem

Figure 12 illustrates steps required to solve the primal problem.
<p align="center">
    <img src="images/step_12_run.png" alt="Running simulation"/>
</p>
<p align="center">Figure 12: Running simulation</p>

### Post processing

Figure 13 illustrates steps required to post process the primal problem.
<p align="center">
    <img src="images/step_13_postprocess.png" alt="Post processing"/>
</p>
<p align="center">Figure 13: Post processing</p>