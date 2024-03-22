---
title: Post processing
keywords: 
tags: [post_processing.md]
sidebar: mpm_application
summary: 
---

## 5. Post-Processing
When the calculation finished, examine the results in the post-processing mode of GiD. For this purpose, navigate in windows explorer to the folder that contains your project files. Then select the file *project_name_Body.post.res* and hand it via *Drag and Drop* over to GiD (drop it over the main window of GiD, e.g. over the meshed problem).

![drag_and_drop](https://user-images.githubusercontent.com/51473791/190993773-bd58206b-a5ff-4c8f-9591-7f411916ecff.jpg)


Then GiD automatically switches to its post-processing surface and loads the results file. Now, we can see the material points that we assigned to the domain of the cantilever in a previous step:

![updated_particles_result](https://user-images.githubusercontent.com/51473791/190993829-0a32d698-038c-4deb-a980-ea12e896ab28.jpg)

Often it is useful to compare the behaviour of the material points with the geometry of the background mesh. To add the background mesh to the post-processing, click on **`Files`&rightarrow;`Merge`** and select the file *project_name_Grid.post.bin* to open.

![merge](https://user-images.githubusercontent.com/51473791/170937564-ecc48ed7-73ee-419c-83d9-fbf9e31c105d.png)

![merge2](https://user-images.githubusercontent.com/51473791/170937565-e41df791-0ad3-4935-b6ea-f32418f62ac7.png)

Then we can see the material points as well as the background domain.

![merged_result_files](https://user-images.githubusercontent.com/51473791/190993887-72c795f9-c4d3-463b-8122-9b70c3c60883.jpg)

To adapt the depiction of the background mesh, we click in the right toolbar in the line of the backgroundmesh on the yellow quadrat that contains a black cross. Subsequently, another dropdown menu opens and we click on the first entry, as depicted in the image below.

![right_toolbar_wa](https://user-images.githubusercontent.com/51473791/170946238-f77bdc44-1705-4bdb-bdc9-e2fbc28b3a98.png)

Now we can see the material points within the outer edges of the background mesh:

![merged_image_2](https://user-images.githubusercontent.com/51473791/191246995-0dce65d8-3f90-48aa-aedb-2139f02612fd.jpg)



To see the deformation of the structure, one has to select the depiction of the displacement of the
material points. This is done by clicking on the following symbol and selecting **`MP Displacement`**.

![choose_mp_displacement](https://user-images.githubusercontent.com/51473791/191246922-99cd2f61-0993-48d2-9520-c8d116e1c741.jpg)

Subsequently, another pop-up window opens where we have to specify the deformation (scaling) factor. In this example, the real deformation is so small that it cannot be seen with a scaling of 1. Therefore we choose a deformation factor of 100. To do this, we click on `1x` and enter in the next pop-up window a scaling factor of 100.

![deformation_factor_wa](https://user-images.githubusercontent.com/51473791/170950810-f1dd3e8c-c739-46a3-b5b6-5e88b0e1b4d9.png)
![new_deformation_factor](https://user-images.githubusercontent.com/51473791/170950838-7cda9d63-4b23-4429-a28b-aa414fb4a8e3.png)

Then we can finally see the deformation of the cantilever beam:

![deformed_structure](https://user-images.githubusercontent.com/51473791/190995668-4b21da90-b2f4-498d-94ef-c454efebab72.jpg)

In some cases it may be useful to see the change of the deformation over the time (e.g. for a dynamic problem). To see the change of deformation over time, click on the play button in the upper toolbar, as shown in the following picture:

![show_video_button](https://user-images.githubusercontent.com/51473791/190995834-88570b1a-3a17-420a-ab1f-e4ae021e079b.jpg)


Note: Since we calculated here a static example, only the initial state of the system and the final displacement are depicted within the deformation over time depiction.

