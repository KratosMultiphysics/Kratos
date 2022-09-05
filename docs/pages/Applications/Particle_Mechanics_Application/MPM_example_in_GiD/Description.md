# Setting up a MPM example in GiD

## 1.Introduction
This tutorial is intended to give an overview on how to pre- and post-process a problem with the material point method (MPM) in [GiD](https://www.gidhome.com/). The general procedure how to set up a problem geometry, specify the boundary conditions, generate the mesh, calculate the solution and display the results, is explained subsequently by going through an example. For this purpose a cantilever beam under dead load is regarded. The corresponding geometric and cross-sectional properties are displayed in the picture below. For the calculation of the problem [KRATOS Multiphysics](https://github.com/KratosMultiphysics/Kratos) is called internally from GiD.

![Structural system](https://user-images.githubusercontent.com/51473791/168762544-750d2f29-6ed7-409d-8205-a6257e6a72ac.png) 



## 2.Set up a Material Point Method problem in GiD
### 2.1 Load the MPM GUI
The very first step is to open GiD. Then we need to load the Kratos GiD GUI. Therefore click on **`Data` &rightarrow; `Problem type` &rightarrow; `Kratos`** in the top toolbar. In the opening window, the Kratos Application Market, select then the **MPM** application. 

![Kratos_Application_Market](https://user-images.githubusercontent.com/51473791/167376588-366ee16a-1ecc-4c00-8b0c-77e903fead76.png)

<p align="center">
<img src="https://user-images.githubusercontent.com/51473791167375311-61660521-3a93-456c-8bf6-7c194efbcd02.png" alt="Kratos application Market" title="Kratos application market"/>
</p>

Following this, another window is going to appear where one has to choose the dimensions of the MPM problem that is to be modelled. For this example we choose **2D** by clicking on the left field.

![GUI_GiD_wa](https://user-images.githubusercontent.com/51473791/168768013-80e01bcd-c7c1-44a9-afd0-337247c7f060.png)

In the sidebar with the title 'Particle mechanics' on the left, we can adjust different parameters for the solution of the problem. The most important are:

- **Solver Type:** Since the regarded problem is static, we select the corresponding solver type.
- **Parts:** Within 'parts' the domains for the body mesh (particles) and the background grid are chosen and important parameters are defined. **Solid** contains the data for the body mesh; **Grid** the data for the background mesh.
- **Boundary conditions:** Here we can assign Dirichlet boundary conditions to the problem domain.
- **Loads and other conditions:** In this section different types of loads and slip conditions can be defined.

### 2.2 Set up the example geometry
In our example we will regard a cantilever beam under gravity forces. As already mentioned above, we choose therefore the solver type **static**.

To set up the geometry of the problem, two domains have to be defined: one for the background grid and one for the body domain. For that purpose we define at first the background domain. By clicking on the **create object**-icon and choosing a rectangle, we define the background domain. In the next step we define another, smaller rectangle that represents the cantilever. 

![create_object_classic_wa](https://user-images.githubusercontent.com/51473791/168771191-6eb60514-6000-40bf-8399-0d094f1a8a4d.png)

After clicking on the indicated button a menu opens. To create a rectangle click on the marked rectangle:

![create_rectangle_classic_wa](https://user-images.githubusercontent.com/51473791/168771227-7d563d15-f525-48ef-aafe-0a11675a52f6.png)

Now follow the instructions given in the command line: 

![command_window_classic_wa](https://user-images.githubusercontent.com/51473791/168772080-2eee839d-3416-49ae-87e8-69133f5573b9.png)

We create two rectangles by entering the coordinates of the left lower and the right upper corner point in the command line. We enter the coordinates in the format shown in the image above and confirm them by pressing **`Enter`**. The first rectangle that is going to represent the background domain has its corner points at **(0,0)** and **(5,2)**. The second one that represents the cantilever has its corner points at **(0, 0.5)** and **(5,1.5)**. Now the geometry should be 
displayed in GiD's main window as in the picture below:

![rectangles_created](https://user-images.githubusercontent.com/51473791/168774252-b19afebe-cd08-4105-80ad-6165e43714aa.png)

Subsequently, we assign the geometries to the different domains, choose element formulations and specify the material parameters of the structures. We start with the body structure. Click on **`Solid`**, then a menu opens above the command line. 

![solid_1](https://user-images.githubusercontent.com/51473791/168551778-89cefe5a-cb50-4bde-ae09-fc29cdd57579.png)

![menu_above_command_line](https://user-images.githubusercontent.com/51473791/168551845-64139237-cb73-4ef6-8c8d-c4651b7b6e52.png)

Here we enter the material parameters. The cantilever in our problem is made of steel, so we enter the corresponding material parameters.
Then we have to assign the data that we just entered to a domain of our problem. This step is very important and error prone. To assign our entered data to the rectangle that represents the cantilever beam click on **`Select`**.

![solid_1_new_wa](https://user-images.githubusercontent.com/51473791/168775682-e3298c63-9d94-4897-bfa3-b1059c2bd094.png)

With the mouse that has now a different cursor, a little quadrat with a smaller quadrat inside, click on the purple rectangle that lies within the area of the 'cantilever-rectangle'. After clicking on it, this inner purple line should change its color to red, as indicated in the image below. 

![solid_3_new_wa](https://user-images.githubusercontent.com/51473791/168777223-ce0e9b92-0889-480e-b822-56c64f730ef2.png)

Then leave the selection of areas by clicking **`ESC`** and confirm your choice by clicking **`OK`**. The same procedure we repeat with the 'background rectangle'. We start by clicking on **`Grid`**. In the opening menu we click on **`Select`** and choose now the inner purple rectangle of the 'background rectangle'. Then we hit **`ESC`** and confirm our choice by clicking **`OK`**. Now we assigned element formulations, material parameters and the domains to our problem geometry. In the next step we will deal with the boundary conditions.

### 2.2 Definition of the boundary conditions
The cantilever has a fixed support on the left. So we set there all displacements permanently to 0. On the right side, the cantilever can deflect freely, so we assign there a slip condition. To impose the Dirichlet boundary condition on the left side, we click on **`Boundary Conditions` &rightarrow; `Displacements`**. 

![boundary_conditions_1_new_wa](https://user-images.githubusercontent.com/51473791/168780032-7d0650fe-5d81-4c9a-a31c-b04a36c93a30.png)


Afterwards, a new menu opens above the command line. Here we choose, as already indicated in the picture, line conditions. We click again on **`Select`** and choose the left edge of the large rectangle. This edge changes its color from blue to red, then we end the selection by hitting **`ESC`** and clicking **`OK`**. 

![slip_cond_1_wa](https://user-images.githubusercontent.com/51473791/168783010-dbfb2fe0-d351-4648-877d-1ba725e1005d.png)

To assign a slip condition on the right side of the background rectangle, click on **`Loads and other conditions`&rightarrow;`Slip`**. Then select the right edge of the background rectangle in the same way as we selected the left edge. Now we imposed all the necessary boundary conditions on the example problem. Note, that the boundary conditions were imposed in this case on the background mesh. In MPM it is also possible to impose boundary conditions on the body mesh. 

### 2.3 Application of Loads

The next step is to add loads to the model. As we regard the cantilever beam under dead load in our example, we don't have to apply additional loads on the structure. GiD resp. KRATOS takes the dead load automatically into account.

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

## 4. Calculation
By now, the setup of the model of the cantilever is complete. To calculate the displacement of the cantilever, click in the top toolbar on **`Calculate`&rightarrow;`Calculate`**, as shown in the image below. 

![calculate_cantilever_wa](https://user-images.githubusercontent.com/51473791/170933305-179a5ab5-c8ff-4b6d-b4f9-6f6c4f8965aa.png)

Then the calculation of the problem is carried out. As already mentioned in the first paragraph, the calculation is done within the [KRATOS Multiphysics](https://github.com/KratosMultiphysics/Kratos) framework. More information on the material point method and its implementation within KRATOS is available under [ParticleMechanicsApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/ParticleMechanicsApplication).

After the calculation, one can find in the project folder of our project the following files, among others:

- MainKratos.py
- ProjectParamters.json
- ParticleMaterials.json
- project_name_Body.mdpa
- project_name_Grid.mdpa

- project_name_Body.post.res
- project_name_Grid.post.bin 

The files 'MainKratos.py'-'project_name_Grid.mdpa' contain the necessary data for the calculation of the problem with KRATOS. The last two files 'project_name_Body.post.res' and 'project_name_Grid.post.bin' contain essential results of the calculation and are used in the post-processing.

## 5. Post-Processing
When the calculation finished, we examine the results in the post-processing mode of GiD. For this purpose, navigate in windows explorer to the folder that contains your project files. Then select the file 'project_name_Body.post.res' and had it via 'Drag and Drop'over to GiD (drop it over the main window of GiD, e.g. over the meshed problem).

![drag_and_drop_wa](https://user-images.githubusercontent.com/51473791/170936451-6f8caa41-a339-4d8c-b251-561a813525fa.png)

Then GiD automatically switches to its post-processing surface and loads the results file.

![results_particles](https://user-images.githubusercontent.com/51473791/170936811-77ba39c9-a142-4a19-8921-2aaa48345643.png)

Now, we can see the particles that we assigned to the domain of the cantilever in a previous step. Often it is useful to compare the behaviour of the particles with the geometry of the background mesh. To add the background mesh to the post-processing, click on **`Files`&rightarrow;`Merge`** and select the file 'project_name_Grid.post.bin' to open.
 
![merge](https://user-images.githubusercontent.com/51473791/170937564-ecc48ed7-73ee-419c-83d9-fbf9e31c105d.png)

![merge2](https://user-images.githubusercontent.com/51473791/170937565-e41df791-0ad3-4935-b6ea-f32418f62ac7.png)

Then we can see the particles as well as the background mesh.

![results_2_wa](https://user-images.githubusercontent.com/51473791/170946323-d54778cd-3a37-4a92-b249-2be14f58be64.png)

To adapt the depiction of the background mesh, we click in the right toolbar in the line of the backgroundmesh on the yellow quadrat that contains a black cross. Subsequently, another dropdown menue opens and we click on the first entry, as depicted in the image below. 

![right_toolbar_wa](https://user-images.githubusercontent.com/51473791/170946238-f77bdc44-1705-4bdb-bdc9-e2fbc28b3a98.png)

Now we can see the particles within the outer edges of the background mesh:

![results_3](https://user-images.githubusercontent.com/51473791/170947386-106bd57f-40d2-491b-a2bd-f5b7695354e3.png)

To see the deformation of the particles due to the dead load of the structure, one has to select the depiction of the displacement of the 
material points (particles). This is done by clicking on the following symbol and selecting **`MP Displacement`**.

![results_4_wa](https://user-images.githubusercontent.com/51473791/170948919-06811936-7dd6-4c78-844b-534ca7dc3bfc.png)

Subsequently, another pop-up window opens where we have to specify the deformation (scaling) factor. In this example, the real deformation is so small that it cannot be seen with a scaling of 1. Therefore we choose a deformation factor of 100. To do this, we click on `1x` and enter in the next pop-up window a scaling factor of 100.

![deformation_factor_wa](https://user-images.githubusercontent.com/51473791/170950810-f1dd3e8c-c739-46a3-b5b6-5e88b0e1b4d9.png)
![new_deformation_factor](https://user-images.githubusercontent.com/51473791/170950838-7cda9d63-4b23-4429-a28b-aa414fb4a8e3.png)

Then we can finally see the deformation of the cantilever beam:
 
![results_5](https://user-images.githubusercontent.com/51473791/170951287-bf2bcb66-ae6f-4afa-8b49-5b65277af171.png)

In some cases it may be useful to see the change of the deformation over the time (e.g. for a dynamic problem). To see the change of deformation over time, click on the play button in the upper toolbar, as shown in the following picture:

![results_6_wa](https://user-images.githubusercontent.com/51473791/170952039-65d49339-d936-4aab-8f3b-2331f5d38f91.png)

Note: Since we calculated here a static example, only the initial state of the system and the final displacement are depicted within the deformation over time depiction.






