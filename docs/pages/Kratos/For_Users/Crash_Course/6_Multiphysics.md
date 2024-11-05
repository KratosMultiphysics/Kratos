---
title: 6 - Multiphysics
keywords: 
tags: [Kratos Crash Course Multiphyscis]
sidebar: kratos_for_users
summary: 
---

# Multiphysics example

# 1. Introduction

This will be the final chapter of the Crash course and you will put in practice everything you have learnt about Kratos.

The goal of this section is to become familiar with key aspects of multiphysics simulations using a prototypical example and also take tha chance to take a close loop to the ProjectParameters.json file. The chosen problem is a case of fluid-thermal interaction. One should become familiar with the various components involved in setting up and running such simulations, as well as the necessary steps to ensure the quality and physical relevance of results.

This tutorial will should one of the most simple multiphyscis applications, but we have tutorials for more complex examples with solver-level coupling in the [wiki]() 

Due to the complexity of the task and the time constraints, we will provide you with all the material. Please download the case from [here]() so you can follow the rest of the sections.

# 2 Description of files

If you look at the files that you downloaded you will see some familiar faces:

- **MainKratos.py**: Our main script
- **ProjectParameters.json**: The configuration file
- **multiphysics_example.mdpa**: The geometry of our problem
- **FluidMaterials.json** & **ConvectionDiffusionMaterials.json**: The json files for the properties of the materials during the fluid and thermal stages of the problem.

and some new files that we will explain:

- **custom_transfer_operation.py**: A custom operation which will perform the information transfer between stages 

## 2.1. MainKratos.py

To begin, let's take a look at the MainKratos.py. You should already be familiar with this file:

```python
import sys
import importlib
import KratosMultiphysics
from KratosMultiphysics.project import Project

if __name__ == "__main__":

    # Check if a custom project parameters filename is provided
    if len(sys.argv) == 1:
        project_parameters_filename = "ProjectParameters.json"
    else:
        project_parameters_filename = str(sys.argv[1])

    # Parse simulation settings and run simulation
    with open(project_parameters_filename, 'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        project = Project(project_parameters)

        orchestrator_reg_entry = KratosMultiphysics.Registry[project.GetSettings()["orchestrator"]["name"].GetString()]
        orchestrator_module = importlib.import_module(orchestrator_reg_entry["ModuleName"])
        orchestrator_class = getattr(orchestrator_module, orchestrator_reg_entry["ClassName"])
        orchestrator_instance = orchestrator_class(project)
        orchestrator_instance.Run()% 
```

As you see, there are some key differences from the one we saw at the begining of the coure. The key change is that instead of creating an `AnalysisStage` we are creating something called `Project` and `Orchestrator`. These classes add one level of abstraction to the script and are used when multiple stages are needed. 

The role of the `Orchestrator` is to decide which `AnalysisStage` and when will be used, as well as executing the code that needs to transfer data between them.

All this information is now contained in the `ProjectParameters.json`

## 2.2. ProjectParameters.json

As with a normal simulation, this file contains the configuration of all components of Kratos using the JSON sintax. Looking at it closely will reveral that it has some differences from the files we have been using as examples until now:

First take a look at the first key:

```json
"orchestrator" : {
    "name" : "Orchestrators.KratosMultiphysics.SequentialOrchestrator",
    "settings" : {
        "echo_level" : 0,
        "execution_list" : ["fluid_stage", "thermal_stage"],
        "load_from_checkpoint" : null,
        "stage_checkpoints" : false
    }
}
```
 
As you can see, there is a new `orchestrator`. As we have seen this will be in charge of scheduling and launching our different `AnalsisStages`. For this particular cas, we are using the `SequentialOrchestrator` because one will be executed after the other.

The code in the main file:

```python
orchestrator_reg_entry = KratosMultiphysics.Registry[project.GetSettings()["orchestrator"]["name"].GetString()]
orchestrator_module = importlib.import_module(orchestrator_reg_entry["ModuleName"])
orchestrator_class = getattr(orchestrator_module, orchestrator_reg_entry["ClassName"])
orchestrator_instance = orchestrator_class(project)
```

Takes this `string` and creates an instance of the selected orchestrato, so you can change it directly in the `ProjectParameters.json` file maintaining a generic `MainKratos.py` script.

Our study case consists on a fluid and a thermal part, and as you can see in the `execution_list` key, those are present:

```json
"stages" : {
    "fluid_stage" : {
        ...
    },
    "thermal_stage" : {
        ...
    }
}
```

Following the orchestrator, we have the `stages` key. As you may have guesses each stages is an `AnalysisStage`. Note that the names of the stages are the same as the ones we added in the `orchestrator`. Inspecting any of the stages, will show us what are their settings.

Let's take a look at the Fluid stage:

```json
"stage_preprocess" : {
    "modelers" : [{
        "name" : "Modelers.KratosMultiphysics.ImportMDPAModeler",
        "parameters" : {
            "input_filename"  : "multiphysics_example",
            "model_part_name" : "FluidModelPart"
        }
    },{
        "name" : "Modelers.KratosMultiphysics.CreateEntitiesFromGeometriesModeler",
        "parameters" : {
            "elements_list"   : [{
                "model_part_name" : "FluidModelPart.Fluid",
                "element_name"    : "Element2D3N"
            }],
            "conditions_list" : [{
                "model_part_name" : "FluidModelPart.Walls",
                "condition_name"  : "WallCondition2D2N"
            }]
        }
    }],
    "operations" : []
}
```

The first section is called `stage_preprocess` and holds the properties of the ojbects that will be executed before the stage is executed. As you can see this is very similar to how `AnalysisStages` and `Processes` are executed.

In this case we use this section to execute two modelers, one for reading the `ModelPart` and one for setting up the type of `Elements` and `Conditions`. We do this because this particular mdpa file has no information about the types that we will use, as will be defferents for the different `AnalysisStages` 

Now let's take a look at the Thermal stage:

```json
"stage_preprocess" : {
    "modelers" : [{
        "name" : "Modelers.KratosMultiphysics.ImportMDPAModeler",
        "parameters" : {
            "input_filename"  : "multiphysics_example",
            "model_part_name" : "ThermalModelPart"
        }
    },{
        "name" : "Modelers.KratosMultiphysics.CreateEntitiesFromGeometriesModeler",
        "parameters" : {
            "elements_list"   : [{
                "model_part_name" : "ThermalModelPart.Fluid",
                "element_name"    : "Element2D3N"
            }],
            "conditions_list" : [{
                "model_part_name" : "ThermalModelPart.Walls",
                "condition_name"  : "WallCondition2D2N"
            }]
        }
    }],
    "operations" : [{
        "name" : "custom_transfer_operation.CustomTransferOperation",
        "parameters" : {
            "input_model_part_name" : "FluidModelPart.Fluid",
            "output_model_part_name" : "ThermalModelPart.Fluid"
        }
    }]
}
```

Most of the preprocess is the same, this time with a different set of `elements`, but this time the `operations` section is not emtpy. This is because we use the `stage_preprocess` to transfer the results from the Fluid stage, which has already finished to the Thermal stage, which is about to begin. If you look at the name, it will be familiar to you: `custom_transfer_operation.CustomTransferOperation`, which is the name of one of the extra files that were present when you downloaded the example: `custom_transfer_operation.py`

Both the `orchestrator` and the `AnalysisStage` are able to inspect your execution folder in order to execute custom `operators`, `modelers` and `processes` that you define.

Let's see what our custom operation does:

## 2.3. Custom Operation

```python
class CustomTransferOperation(KratosMultiphysics.Operation):

    def __init__(self, Model, Settings):
        KratosMultiphysics.Operation.__init__(self)

        self.model = Model
        self.input_model_part_name = Settings["input_model_part_name"].GetString()
        self.output_model_part_name = Settings["output_model_part_name"].GetString()

    def Execute(self):

        input_model_part = self.model.GetModelPart(self.input_model_part_name)
        output_model_part = self.model.GetModelPart(self.output_model_part_name)

        for node_in, node_out in zip(input_model_part.Nodes, output_model_part.Nodes):
            velocity = node_in.GetSolutionStepValue(KratosMultiphysics.VELOCITY, 0)
            node_out.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, velocity)
```

As you can see, the code is very simple: We get two `ModelParts`, the input (Fluid) and the output (Thermal), and we iterate over the nodes of both. For every node, we read the `VELOCITY` variable from the historical database of the `node` in the input `ModelPart` and set its value to the historical database of the corresponding node of the output `ModelPart` this case there is 1-1 relation because we have used the same geometry in both `AnalysisStages`, but you could use a more complex mapping function.

# 3 Executing and checking the Results

You have all the necessary script and data to run the simulation, so let's execute it. As with the first example, the only thing you need to do is to call the `MainKratos.py` script:

```bash
python MainKratos.py
```

As this is a multiphyscis example, it may take a little more time to complete. Once it finished you should see the results.
