---
title: RomApplication Tutorial with RomManager
keywords:
tags: []
sidebar: rom_application
summary:
---

The objective of this tutorial is to show Kratos users how to build ROMs using the Kratos RomApplication.

In Katos, the [RomManager](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/RomApplication/python_scripts/rom_manager.py) is the class that seamlessly orschestrates the simulations involved in both stages. ![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/RomApp_AddFigures/RomApp_Tutorial/Figures/RomTutorial_1_1.png)

In the subsequent sections, the parameters passed to the RomManager will be introduced.


## Setting up a Structural Mechanics Parametric Simulation

The first example is a structural mechanics simulation with three pressure loads applied as shown in the figure

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/RomApp_Tutorial/Figures/RomTutorial_2_1.png)

The material follows a Neo-Hookean Hyperelastic constitutive law, and the bottom part is fixed. The step-by-step procedure for generating this geometry in GiD is explained in [this video](https://youtu.be/3gJIHf5gQ88?si=5gPumMJTlYwBL0e3). Moreover, the geometry files can be obtained [here](https://github.com/KratosMultiphysics/Documentation/tree/master/RomApp_Tutorial/RomAppTutorial_Part2).



## Setting up a Structural Mechanics ROM

A video version of this section is available [here](https://youtu.be/KtO-XxbgLwU?si=MEmMpKlW1LOJaCTW), while the files used are available [here](https://github.com/KratosMultiphysics/Documentation/tree/master/RomApp_Tutorial/RomAppTutorial_Part3).


The [RomManager](https://github.com/KratosMultiphysics/Kratos/blob/master/applications/RomApplication/python_scripts/rom_manager.py) accepts a KratosParameters object as the key word argument * general_rom_manager_parameters*. In this tutorial section we consider the following subset of parameters:

```python
def rom_manager_parameters():

    default_settings = KratosMultiphysics.Parameters("""{
        "rom_stages_to_train" : ["ROM", "HROM"],
        "rom_stages_to_test" : [],
        "save_gid_output": true,
        "save_vtk_output": true
        "ROM":{
            "use_non_converged_sols": true,
            "model_part_name": "Structure",
            "nodal_unknowns": ["DISPLACEMENT_X","DISPLACEMENT_Y"]
        }
    }""")

    return default_settings
```

We setup the simplest workflow by creating a rom_manager object and passing to it the KratosParameters returned by the function created above, that is:

```python
if __name__ == "__main__":
    rom_manager = RomManager(general_rom_manager_parameters=rom_manager_parameters())
    rom_manager.Fit()
    rom_manager.PrintErrors()
```

Indeed, by running the above code, a unique case is launched, obtaining a single column vector, which is used for creating a ROM. Then, a ROM is also launched and the results of FOM and ROM are compared with the PrintErrors() method. This workflow looks as follows:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/RomApp_Tutorial/Figures/RomTutorial_3_1.png)


The parameters vector $\mu$ contains the pressure loads speficied in the ProjectParameters.json, that is [10, 10, 10].

There are 3 ways of imposing different parameters vectors $\mu$ using the RomManager, they are:

- Specify a function that modifies the KratosParameters object from the ProjectParameters.json. This is useful for impossing boundary conditions, initial conditions, loads, or other parameter variation imposed via a Kratos Process.
- Specify a function that modifies the json file containg the material. This re-writes, with the desired material property, the json file that is loaded by each simulation.
- Speficy a custom analysis stage. This allows to either impose a parameter variation or acces any of the methods of the derived Analysis stage class.


We now focus on the first option. The function modifying the KratosParameters object to impose different values of pressure is the following:

```python
def UpdateProjectParameters(parameters, mu):
    parameters["processes"]["loads_process_list"][0]["Parameters"]["value"].SetDouble(mu[0])
    parameters["processes"]["loads_process_list"][1]["Parameters"]["value"].SetDouble(mu[1])
    parameters["processes"]["loads_process_list"][2]["Parameters"]["value"].SetDouble(mu[2])

    return parameters

```


This function should be passed as the key word argument *UpdateProjectParameters* to the rom_manager object as follows:

```python
rom_manager = RomManager(general_rom_manager_parameters=rom_manager_parameters(), UpdateProjectParameters=UpdateProjectParameters)

```

Finally, the $\mu$ vector is passed to the Fit() method as a list of lists.

```python
mu_train = [[20,45,60]]
rom_manager.Fit(mu_train)
rom_manager.PrintErrors()
```

The workflow for the above script is the one presented in Sec. Overview of ROMs in Kratos.


The Test() method of the RomManager also accepts a list of lists specifying the testing parameters. This method runs both the FOM and ROM, but the snapshots from the FOM are not used for improving the ROM. An example snippet of the Test() method of the RomManager is the following:

```python
mu_test = [[200,450,600]]
rom_manager.Test(mu_test)
rom_manager.PrintErrors()
```

The workflow, when calling the Test() method looks like this:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/RomApp_Tutorial/Figures/RomTutorial_3_2.png)

The RunFOM() method of the RomManager allows to launch the FOM simulations without introducing overheads. No extra data (besides the Results files if they are chosen to be kept with the flags *"save_gid_output"* and *"save_vtk_output"*) is generated.  An example snippet of the RunROM() method of the RomManager is the following:

```python
mu_run = [[2,4,6]]
rom_manager.RunFOM(mu_run)
```

The workflow, when calling the RunFOM() method looks like this:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/RomApp_Tutorial/Figures/RomTutorial_3_3.png)


Finally, the RunROM() method of the RomManager launches the ROM simulations for the parameters passed, without storing extra data. The workflow, when calling the RunROM() method looks like this:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/RomApp_Tutorial/Figures/RomTutorial_3_4.png)



## Setting up a Fluid Dynamics Parametric Simulation
to be added
