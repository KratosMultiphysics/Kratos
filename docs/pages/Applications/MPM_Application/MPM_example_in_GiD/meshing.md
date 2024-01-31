---
title: Meshing
keywords: 
tags: [meshing.md]
sidebar: mpm_application
summary: 
---

## 3. Meshing
Before we start the meshing of the problem, we save our project by clicking on **`Files`&rightarrow;`Save`** in the top toolbar. 
 
For the Meshing we click on **`Mesh`** in the upper toolbar. Then a dropdown menu opens with different meshing options. In MPM the background grid and the body have different meshes, so two meshes are to be defined. For the body mesh we choose a structured mesh of quadrilaterals with a mesh size of **0.05**. To assign it to the body, one clicks on **`Mesh`&rightarrow;`Structured`&rightarrow;`Surfaces`&rightarrow;`Assign sizes`**. Subsequently, select the area of the cantilever by clicking on the purple rectangle within it. The inner rectangle of the cantilever changes its color from purple to red once it is selected. This can be seen in the picture below.

![assign_body_mesh_1](https://user-images.githubusercontent.com/51473791/168785222-28198d44-991a-4d50-9c42-892504e7b806.png)


Then a pop-up window appears where we have to enter the mesh size of lines. Here we enter **0.05**, too (the mesh size of the assigned lines creates the mesh size of the area). 

![bodymesh_enter_sizes](https://user-images.githubusercontent.com/51473791/190999754-934c1ee9-e118-4ddc-b9b7-828a5fe6ae1a.jpg)

To assign a structured mesh to the area of the cantilever, we have to specify the edges of the area. Therefore, we click firstly on one of the horizontal (blue) lines of the area of the cantilever.

![bodymesh_enter_surrounding_lines_1](https://user-images.githubusercontent.com/51473791/191000193-4f8d1395-8570-4e1f-8e28-db843173b354.jpg)

To select the vertical lines of the cantilever, draw a rectangle (click into the drawing plane, hold the left button of the mouse and move the pointer), as indicated in the picture below, over the area of the cantilever. 

![bodymesh_enter_surrounding_lines_2](https://user-images.githubusercontent.com/51473791/191000262-fe2cf0a1-5b6e-4d56-a8cb-376a54a54359.jpg)

As the vertical edges of the larger and the smaller rectangle (boundaries of the background- and body-domain) lie above each other, both of them are selected by the created rectangle. However, since only the edges of the smaller one belong to the inner surface, an error-message appears:

![bodymesh_enter_surrounding_lines_3](https://user-images.githubusercontent.com/51473791/191000297-1a883ceb-d3c1-4782-a969-edf62fad323c.jpg)

After confirming the error-message, the following lines should be selected (in case the horizontal lines aren't red anymore, just select them again):

![final_lines](https://user-images.githubusercontent.com/51473791/191000366-28d219af-b8ae-470c-8dde-bc3379175dba.jpg)

Subsequently, press *'ESC'* and assign the mesh size **0.05** to the surrounding lines. 

![size_on_lines](https://user-images.githubusercontent.com/51473791/191000597-291537e2-2b47-46c7-b10e-9c0450828c9d.jpg)

The background grid is created in an analogous manner. Instead of selecting the area of the cantilever, select now the area of the background domain (large purple rectangle). Assign a mesh size of **0.2** to the area and the surrounding edges of the background domain.

**Note:** It is generally recommended to choose the size of the body mesh smaller than half the size of the background mesh. This leads to more meaningful results and enhances the stability of the calculation process. 

After adjusting the meshing parameters, create the mesh by clicking on **`Mesh`&rightarrow;`Generate Mesh`**. Another pop-up window will appear that asks for a mesh size. 

![generate_mesh_1](https://user-images.githubusercontent.com/51473791/168788863-66d50217-3464-4a06-9467-4221f8817bfb.png)

The meshing parameters, which were defined in the previous step, override the size that is displayed in the pop-up window. So simply click **`Ok`**. A further pop-up window shows the number of nodes, elements, etc. that were generated. Finally, the meshed geometry is displayed.

![refined_mesh](https://user-images.githubusercontent.com/51473791/191000975-fb7b4628-4276-4f12-8e06-6d0c055f22f0.jpg)

