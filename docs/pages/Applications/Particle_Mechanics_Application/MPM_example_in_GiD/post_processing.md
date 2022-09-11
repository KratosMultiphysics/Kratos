---
title: Post processing
keywords: 
tags: [post_processing.md]
sidebar: particle_mechanics_application
summary: 
---

## 5. Post-Processing
After the calculation finished, we examine the results in the post-processing mode of GiD. For this purpose, navigate in windows explorer to the folder that contains your project files. Then select the file *project_name_Body.post.res* and hand it via *Drag and Drop* over to GiD (drop it over the main window of GiD, e.g. over the meshed problem).

![drag_and_drop_wa](https://user-images.githubusercontent.com/51473791/170936451-6f8caa41-a339-4d8c-b251-561a813525fa.png)

Then GiD automatically switches to its post-processing surface and loads the results file.

![results_particles](https://user-images.githubusercontent.com/51473791/170936811-77ba39c9-a142-4a19-8921-2aaa48345643.png)

Now, we can see the particles that we assigned to the domain of the cantilever in a previous step. Often it is useful to compare the behaviour of the particles with the geometry of the background mesh. To add the background mesh to the post-processing, click on **`Files`&rightarrow;`Merge`** and select the file *project_name_Grid.post.bin* to open.
 
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

