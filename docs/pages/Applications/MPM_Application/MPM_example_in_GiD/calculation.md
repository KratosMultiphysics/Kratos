---
title: Calculation
keywords: 
tags: [Description.md]
sidebar: mpm_application
summary: 
---

## 4. Calculation
By now, the setup of the model of the cantilever is complete. To calculate the displacement of the cantilever, click in the top toolbar on **`Calculate`&rightarrow;`Calculate`**, as shown in the image below.

![calculation_refined_mesh](https://user-images.githubusercontent.com/51473791/191244080-393633b9-19c8-4123-969d-c9c4ef4431aa.jpg)


Then the calculation of the problem is carried out. As already mentioned in the first paragraph, the calculation is done within the [KRATOS Multiphysics](https://github.com/KratosMultiphysics/Kratos) framework. More information on the material point method and its implementation within KRATOS is available under [MPMApplication](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/MPMApplication).

After the calculation, one can find in the project folder the following files (among others):

- MainKratos.py
- ProjectParamters.json
- ParticleMaterials.json
- project_name_Body.mdpa
- project_name_Grid.mdpa

- project_name_Body.post.res
- project_name_Grid.post.bin

The files *MainKratos.py*, *ProjectParamters.json*, *ParticleMaterials.json*, *project_name_Body.mdpa* and *project_name_Grid.mdpa* contain the data  that is required for the calculation of the problem with KRATOS.
<br/><br/>
The last two files *project_name_Body.post.res* and *project_name_Grid.post.bin* contain essential results of the calculation and are used in the post-processing.
