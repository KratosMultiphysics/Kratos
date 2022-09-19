---
title: Meshing
keywords: 
tags: [meshing.md]
sidebar: particle_mechanics_application
summary: 
---

## 3. Meshing
Before we start the meshing of the problem, we save our project by clicking on **`Files`&rightarrow;`Save`** in the top toolbar. 
 
For the Meshing we click on **`Mesh`** in the upper toolbar. Then a dropdown menu opens with different meshing options. In MPM the background grid and the body have different meshes, so two meshes are to be defined. For the body mesh we choose a structured mesh of quadrilaterals with a mesh size of **0.05**. To assign it to the body, one clicks on **`Mesh`&rightarrow;`Structured`&rightarrow;`Surfaces`&rightarrow;`Assign sizes`**.

![assign_body_mesh_1](https://user-images.githubusercontent.com/51473791/168785222-28198d44-991a-4d50-9c42-892504e7b806.png)


In a first step, we assign a structured mesh to the cantilever (body mesh). Therefore select the area of the cantilever by clicking on the purple rectangle within it that changes its color from purple to red as a consequence. One can see this in the picture above. Thereafter, a pop-up window appears where we have to enter the mesh size of lines. Here we enter in our example **0.05**. 

![bodymesh_enter_sizes](https://user-images.githubusercontent.com/51473791/190999754-934c1ee9-e118-4ddc-b9b7-828a5fe6ae1a.jpg)

As we assign a strucutured mesh to the area of the cantilever we also have to specify the edges of the area to create the structure of the mesh. The mesh size of the assigned lines creates the mesh size of the area. 

In the example we choose for the cantilever its the upper and lower edge. Usually one should provide four edges for a good mesh generation, but in this case the outer edges of cantilever area are part of the background rectangle. So we can't select them.

Therefore, we click firstly on one of the horizontal (blue) lines of the area of the cantilever.

![bodymesh_enter_surrounding_lines_1](https://user-images.githubusercontent.com/51473791/191000193-4f8d1395-8570-4e1f-8e28-db843173b354.jpg)

To select the vertical lines of the cantilever, draw a rectangle, as indicated in the picture below, over the area of the cantilever. 

![bodymesh_enter_surrounding_lines_2](https://user-images.githubusercontent.com/51473791/191000262-fe2cf0a1-5b6e-4d56-a8cb-376a54a54359.jpg)

As the vertical edges of the larger and the smaller rectangle (boundaries of the background- and body-domain) lie above each other, both of them are selected by the created rectangle. But since only the edges of the smaller one belong to the inner surface, an error-message appears:

![bodymesh_enter_surrounding_lines_3](https://user-images.githubusercontent.com/51473791/191000297-1a883ceb-d3c1-4782-a969-edf62fad323c.jpg)

After confirming the error-message, the following lines should be selected (in case the horizontal lines aren't red anymore, just click on one of the once again):

![final_lines](https://user-images.githubusercontent.com/51473791/191000366-28d219af-b8ae-470c-8dde-bc3379175dba.jpg)

Subsequently, press *'ESC'* and enter the meshsize of the surrounding lines. Here we enter **0.05** once more.

![size_on_lines](https://user-images.githubusercontent.com/51473791/191000597-291537e2-2b47-46c7-b10e-9c0450828c9d.jpg)

The background grid is created in an analogous manner. Instead of the area of the cantilever we choose the area of the background domain and here we can assign a size to the four edges of the rectangle for the mesh generation.

**Note:** Altough we chose in this example the same size for body- and background mesh, it is generally recommended to choose the size of the body mesh half as large as the size of the background mesh. This leads to more meaningful results and enhances the stability of the calculation process. 

When we are done with adjusting the meshing parameters, we create the mesh by clicking on **`Mesh`&rightarrow;`Generate Mesh`**. Another pop-up window will appear that asks for a mesh size. 

![generate_mesh_1](https://user-images.githubusercontent.com/51473791/168788863-66d50217-3464-4a06-9467-4221f8817bfb.png)

The meshing parameters we defined in the previous step override the size that is displayed in the pop-up window. So we simply click **`Ok`**. A further pop-up window informs us about the number of elements, etc. that were generated. Thereafter, the meshed geometry is displayed.

![refined_mesh](https://user-images.githubusercontent.com/51473791/191000975-fb7b4628-4276-4f12-8e06-6d0c055f22f0.jpg)

