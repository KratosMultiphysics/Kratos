---
title: Meshing
keywords: 
tags: [meshing.md]
sidebar: particle_mechanics_application
summary: 
---

## 3. Meshing
Before we start the meshing of the problem, we save our project by clicking on **`Files`&rightarrow;`Save`** in the top toolbar. 
 
For the Meshing we click on **`Mesh`** in the upper toolbar. Then a dropdown menu opens with different meshing options. In MPM the background grid and the body have different meshes, so two meshes are to be defined. For the body mesh we choose a structured mesh of quadrilaterals with a mesh size of **0.2**. To assign it to the body, one clicks on **`Mesh`&rightarrow;`Structured`&rightarrow;`Surfaces`&rightarrow;`Assign sizes`**.

![assign_body_mesh_1](https://user-images.githubusercontent.com/51473791/168785222-28198d44-991a-4d50-9c42-892504e7b806.png)


In a first step, we assign a structured mesh to the cantilever (body mesh). Therefore select the area of the cantilever by clicking on the purple rectangle within it that changes its color from purple to red as a consequence. One can see this in the picture above. Thereafter, a pop-up window appears where we have to enter the mesh size of lines. Here we enter in our example **0.2**. 

![assign_body_mesh_11](https://user-images.githubusercontent.com/51473791/168787212-de5a83b8-7c5a-4da3-b870-349ddf3d6294.png)

As we assign a strucutured mesh to the area of the cantilever we also have specify edges of the area to create the structure of the mesh. The mesh size of the assigned lines creates the mesh size of the area. In the example we choose for the cantilever its the upper and lower edge. Usually one should provide four edges for a good mesh generation, but in this case the outer edges of cantilever area are part of the background rectangle. So we can't select them.

![assign_body_mesh_2](https://user-images.githubusercontent.com/51473791/168787234-3f996da7-f688-449b-a7e3-5703debf3888.png)

The background grid is created in an analogous manner. Instead of the area of the cantilever we choose the area of the background domain and here we can assign a size to the four edges of the rectangle for the mesh generation.

**Note:** Altough we chose in this example the same size for body- and background mesh, it is generally recommended to choose the size of the body mesh half as large as the size of the background mesh. This leads to more meaningful results and enhances the stability of the calculation process. 

When we are done with adjusting the meshing parameters, we create the mesh by clicking on **`Mesh`&rightarrow;`Generate Mesh`**. Another pop-up window will appear that asks for a mesh size. 

![generate_mesh_1](https://user-images.githubusercontent.com/51473791/168788863-66d50217-3464-4a06-9467-4221f8817bfb.png)

The meshing parameters we defined in the previous step override the size that is displayed in the pop-up window. So we simply click **`Ok`**. A further pop-up window informs us about the number of elements, etc. that were generated. Thereafter, the meshed geometry is displayed.

![meshed_geometry](https://user-images.githubusercontent.com/51473791/168790908-22e3b59b-2ef8-4346-8035-031df06e0639.png)
