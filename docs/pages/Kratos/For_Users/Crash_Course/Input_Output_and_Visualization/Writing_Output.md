---
title: Writing Output File
keywords: 
tags: [Python Script Tutorial Writing Output File]
sidebar: kratos_for_users
summary: 
---

In this tutorial the procedure for writing the mesh and data from a `ModelPart` to output files will be described briefly. More information can be found here

## Starting
First of all we need to create a python file with following code to import the *Kratos*, create a `ModelPart` and read it from input as described in the previous tutorial :

```python
from KratosMultiphysics import *
import KratosMultiphysics.FluidDynamicsApplication

this_model = Model()
fluid_model_part = this_model.CreateModelPart("FluidPart")

fluid_model_part.AddNodalSolutionStepVariable(VELOCITY)
fluid_model_part.AddNodalSolutionStepVariable(PRESSURE)

fluid_model_part_io = ModelPartIO("path/to/file/example")
fluid_model_part_io.ReadModelPart(fluid_model_part)
```
{: data-lang="Python"}

## Creating an instance of the *GiD* output utility

The `GiDOutputProcess` helper class can be used to write GiD output from Python. It works mostly as a standard *Kratos* `Process` but it has an extra methods to write the to the output file.

Before we begin, we will need to import the module that defines the `GiDOutputProcess` class:

```python
from gid_output_process import GiDOutputProcess
```
{: data-lang="Python"}

The GiD output class can be instantiated by passing the `ModelPart` you intend to print, the name of the file, and the configuration `Parameters` to it. The example `*.json` file for this tutorial already contains the required block, so we just need to retrieve it:

```Pyhon
gid_output = GiDOutputProcess(
    fluid_model_part,
    'FluidModelPart',
    ProjectParameters["output_configuration"]
)
```
{: data-lang="Python"}

We can check the settings passed as `Parameters` to see a few common configuration options:

```Pyhon
print(ProjectParameters["output_configuration"].PrettyPrintJsonString())
```
{: data-lang="Python"}

```json
"result_file_configuration" : {
    "gidpost_flags"       : {
        "GiDPostMode"           : "GiD_PostBinary",
        "WriteDeformedMeshFlag" : "WriteDeformed",
        "WriteConditionsFlag"   : "WriteConditions",
        "MultiFileFlag"         : "SingleFile"
    },
    "file_label"          : "time",
    "output_control_type" : "step",
    "output_frequency"    : 1,
    "body_output"         : true,
    "node_output"         : false,
    "skin_output"         : false,
    "plane_output"        : [],
    "nodal_results"       : ["VELOCITY","PRESSURE"],
    "gauss_point_results" : []
},
"point_data_configuration"  : []
```
{: data-lang="Json"}

### Output configuration

- The `GiDPostMode` option controls whether the output will be text-based (`"GiD_PostAscii"`) or binary (`"GiD_PostBinary"`). The **ASCII** format is only recommended for debugging, where the user needs to manually access and read the output data. In general, the binary format should be preferred as it yields much smaller output files and can be read and written much faster.

- `WriteDeformedMeshFlag` specifies whether the meshes are written in deformed or undeformed configuration. Setting it to `"WriteDeformedMesh"` will cause the meshes to be written in deformed state, while `"WriteUndeformedMesh"` causes the meshes to be written in their original configuration.

- `WriteConditionsFlag` is used to switch between the option to write only the elements of each mesh (option `"WriteElementsOnly"`) or to write additionally the Condition objects in the current model part (option `"WriteConditions"`).

- `MultiFileFlag` determines whether all are written in a single file or to individual files for each time step. If the flag is set to `"SingleFile"`, all meshes are written into one file. Note that this is only possible in Binary mode. If instead the flag is set to `"MultipleFiles"`, a new file (or, in **ASCII** mode, set of mesh and result files) is generated for each step. If multiple files are used, the names for the files are chosen automatically according to the `file_label` option (both `"time"` and `"step"` based labels are supported).

### Choosing output variables

The configuration also gives us the option to select which results will be printed on the nodes (and integration points) of the mesh. Additional results can be selected by adding the corresponding variable name (as a string argument) to the `nodal_results` or `gauss_point_results` list arguments.

 > Note that it is currently not possible to print the same variable both on nodes and on integration points.

## Using the *GiD* output utility

The usage of the *GiD* output class follows the structure of `Kratos` processes, with one exception. The following methods should be called before starting the solution loop:

```python
gid_output.ExecuteInitialize()
gid_output.ExecuteBeforeSolutionLoop()
```
{: data-lang="Python"}

During the solution loop, one should call, at each time step:

```python
gid_output.ExecuteInitializeSolutionStep()

## The actual solution goes here

if gid_output.IsOutputStep():

    # Call ExecuteBeforeOutputStep on any auxiliar processes you have here

    gid_output.PrintOutput()

    # Call AfterBeforeOutputStep on any auxiliar processes you have here

gid_output.ExecuteFinalizeSolutionStep()
```
{: data-lang="Python"}

Finally, once the solution is done, call
```python
gid_output.ExcecuteFinalize()
```
{: data-lang="Python"}

to close any remaining open files.

**Next** [Nodes and Nodal Data](../Data_Structure/Nodes_And_Data)<br>
**Prev** [Reading ModelPart From Input File](Reading_Input)