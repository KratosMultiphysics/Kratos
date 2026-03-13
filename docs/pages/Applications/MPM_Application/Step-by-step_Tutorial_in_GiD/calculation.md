---
title: 4. Calculation
keywords: 
tags: [Description.md]
sidebar: mpm_application
summary: 
---

The setup of the model of the cantilever is complete. As already mentioned
before, the calculation of the displacement of the cantilever is carried out
using the [Kratos Multiphysics](https://github.com/KratosMultiphysics/Kratos)
framework.

The computation can be performed directly in GiD, [if Kratos is linked to
it](https://github.com/KratosMultiphysics/GiDInterface/blob/master/README.md),
or it can be started from the command line. We proceed by showing both
the approaches.

## Running Kratos within GiD

To perform the simulation by running Kratos within GiD, click in the top
toolbar on **`Calculate` &#8594; `Calculate`**, as shown in the image below.

<img src="https://user-images.githubusercontent.com/51473791/191244080-393633b9-19c8-4123-969d-c9c4ef4431aa.jpg" style="display: block; margin: auto; margin-bottom: 20px">

After the calculation is completed, one can find in the project folder the
following files (among others):

- `MainKratos.py`
- `ProjectParamters.json`
- `ParticleMaterials.json`
- `project_name_Body.mdpa`
- `project_name_Grid.mdpa`

- `project_name_Body.post.res`
- `project_name_Grid.post.bin`

The files `MainKratos.py`, `ProjectParamters.json`, `ParticleMaterials.json`,
`project_name_Body.mdpa` and `project_name_Grid.mdpa` contain the [data that is
required for the calculation of the problem with Kratos](../Input_Files/overview).

The last two files `project_name_Body.post.res` and
`project_name_Grid.post.bin` contain essential results of the calculation and
are used in the post-processing.

## Running Kratos from the command line

In the top toolbar, click on **`Kratos` &#8594; `Write Calculation Files - No Run`**.

The following files are generated in the project folder:

- `MainKratos.py`
- `ProjectParamters.json`
- `ParticleMaterials.json`
- `project_name_Body.mdpa`
- `project_name_Grid.mdpa`

These files contain the [data and settings required for the calculation of the
problem with Kratos](../Input_Files/overview).

Then, open the terminal, navigate to the project folder and execute the command

```bash
python MainKratos.py
```

During the simulation process, the files `project_name_Body.post.res` and
`project_name_Grid.post.bin` are generated. These files contain essential
calculation results and are used for post-processing.

