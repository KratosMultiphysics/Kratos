---
title: Set up of example
keywords: 
tags: [set_up.md]
sidebar: mpm_application
summary: 
---
## 2.Set up a Material Point Method problem in GiD
### 2.1 Load the MPM GUI
The very first step is to open GiD. Then we need to load the Kratos GiD GUI. Therefore click on **`Data` &rightarrow; `Problem type` &rightarrow; `Kratos`** in the top toolbar. In the opening window, the Kratos Application Market, select then the **MPM** application.

![Kratos_Application_Market](https://user-images.githubusercontent.com/51473791/167376588-366ee16a-1ecc-4c00-8b0c-77e903fead76.png)

<p align="center">
<img src="https://user-images.githubusercontent.com/51473791167375311-61660521-3a93-456c-8bf6-7c194efbcd02.png" alt="Kratos application Market" title="Kratos application market"/>
</p>

Following this, another window is going to appear where one has to choose the dimensions of the MPM problem that is to be modelled. For this example we choose **2D** by clicking on the left field.

![GUI_GiD_wa](https://user-images.githubusercontent.com/51473791/168768013-80e01bcd-c7c1-44a9-afd0-337247c7f060.png)

In the sidebar with the title 'MPM' on the left, we can adjust different parameters for the solution of the problem. The most important are:

- **Solver Type:** Since the regarded problem is static, we select the corresponding solver type.
- **Parts:** Within 'parts' the domains for the body mesh (material points) and the background grid are chosen and important parameters are defined. **Solid** contains the data for the body mesh; **Grid** the data for the background mesh.
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
