 
## Adjoint Sensitivity Analysis

### General remarks:

This feature provides the framework to compute sensitivities of structural responses (e.g. displacements, strain energy or stresses) with respect to different types of design variables (e.g. nodal coordinates, material or cross-sectional properties or load intensity) with the adjoint approach. Therefore for each response function an adjoint problem has to be solved. The sensitivies are than computed in a post-processing step. The implemented sensitivity analysis uses a so called semi-analytic approach which means that the derivatives at element level are computed by finite differences.

*Please note:* 
- This feature currently only works for linear problems
- This feature makes use of the HDF5Application

### Features:  

- Response utilities (response functions):
    * Base class of structural response functions
    * Strain energy
    * Displacement or rotation of a node 
    * Stress resultant of a single element
  
- Schemes:
	* Scheme to solve the adjoint problem

- Processes:
    * replacement process (replaces all elements and conditions of a model with its adjoint equivalents and vice versa)

- Adjoint *Neumann* conditions:
    * Point load (derived from PointLoadCondition)
    * Surface load (derived from SurfaceLoadCondition3D)
   
- Structural adjoint elements:
    * Uni-dimensional elements:
       	* Linear 3D beam element (derived from CrBeamElementLinear3D2N)
    * Two-dimensional elements:
        * Thin triangular shell (derived from ShellThinElement3D3N)

*Please note:* 
The adjoint elements and conditions are derived from elements/conditions of the Structural Mechanics Application and can not be seen independently from them. Rather they have to be traced as additions of their parents in the context of adjoint sensitivity analysis. So basic tasks like the computation of the stiffness matrix are not overwritten. The main task of the adjoint elements/conditions is to derive different quantities with respect to the design variable or state (e.g. the right hand side or post-processing results like stresses).

### Usage: 
In order to perform a sensitivity analysis for one response function, the solutions of two linear static problems are necessary: The primal and the adjoint problem. 

*Please note:* 
For the solution of the two problems different kind of variables are used in order to store the results. For the primal problem the usual variables ```DISPLACEMENT``` and ```ROTATION``` and for the adjoint problem ```ADJOINT_DISPLACEMENT``` and ```ADJOINT_ROTATION``` are used.

#### Definition of the Primal Problem
The primal problem can be defind by the regular input files which are needed for an usual linear static analysis. As only difference the output process of the HDF5Application has to be added to the ```list_other_processes``` in the project parameters:

```python
    "list_other_processes" :[{
        "kratos_module" : "KratosMultiphysics.HDF5Application",
        "python_module" : "single_mesh_primal_output_process",
        "help"          : "",
        "process_name"  : "",
        "Parameters" : {
            "model_part_name" : "Structure",
            "file_settings" : {
                "file_access_mode" : "truncate"
            },
            "model_part_output_settings" : {
                "prefix" : "/ModelData"
            },
            "nodal_results_settings" : {
                "list_of_variables": ["DISPLACEMENT", "ROTATION"]
            }
        }
    } 
```

#### Definition of the Adjoint Problem
In order to define the adjoint problem an additional *.json-file for the adjoint project parameters is necessary. This input file is in principle very similar to the respective file of the primal analysis. In comparsion to a regular file for a linear static analysis three points have to be modified:
- ```solver_settings``` by using the ```adjoint_structural_solver``` as ```solver_type``` and by the definition of the ```response_function_settings```
- The input process of the HDF5Application has to be added to the ```list_other_processes``` in order to read the primal solution
- When defining *Dirichlet* conditions in the ```constraints_process_list``` the ```variable_name``` has to be modified to ```ADJOINT_DISPLACEMENT``` respective ```ADJOINT_ROTATION```

For example the ```solver_settings``` can be look like this (Hints for the ```response_function_settings``` are given below):

```python
    "solver_settings"                  : {
        "solver_type"                  : "adjoint_structural_solver",
        "scheme_settings" : {
            "scheme_type"              : "structural"
            },
        "response_function_settings" : {
                "response_type"     : "adjoint_nodal_displacement",
                "gradient_mode"     : "semi_analytic",
                "sensitivity_model_part_name" : "Parts_Beam",
                "nodal_sensitivity_variables"  : ["SHAPE_SENSITIVITY"],
                "element_sensitivity_variables"  : ["I22"],
                "condition_sensitivity_variables"  : [],
                "step_size"         : 1e-6,
                "traced_node"       : 6,
                "traced_dof"        : "DISPLACEMENT_Z"

            },
        "echo_level"                   : 0,
        "problem_domain_sub_model_part_list" : ["Parts_Beam"],
        "processes_sub_model_part_list"      : ["DISPLACEMENT_support","ROTATION_support"],
        "computing_model_part_name" : "computing_domain",
        "rotation_dofs"                      : true,
        "linear_solver_settings"       : {
            "solver_type"         : "Super_LU"
        },
        "model_import_settings"        : {
            "input_type"     : "mdpa",
            "input_filename" : "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/Beam_structure"
        },
        "material_import_settings" :{
            "materials_filename": "adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/materials_beam.json"
        }
    }
```

and the ```list_other_processes``` like:

```python
 "loads_process_list"       : [],
    "list_other_processes" :[{
        "kratos_module" : "KratosMultiphysics.HDF5Application",
        "python_module" : "single_mesh_primal_input_process",
        "help"          : "",
        "process_name"  : "",
        "Parameters" : {
	        "model_part_name" : "Structure",
            "file_settings" : {
                "file_access_mode" : "read_only"
            },
            "nodal_results_settings" : {
                "list_of_variables": ["DISPLACEMENT", "ROTATION"]
            }
        }
     }
```
#### How to run the analysis:

Are all necessary input-files defined the analysis can be performed with a simple python script by calling the primal and afterwards the adjoint analysis. The sensitivties are computed in a post-processing of the adjoint problem.

A possible python code can look like this:
```python 
    # Solve the primal problem     
    with open("concrete_building_parameters.json",'r') as parameter_file:
        ProjectParametersPrimal = Parameters( parameter_file.read())
    primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal)
    primal_analysis.Run()
    # Solve adjoint problem and compute sensitivities
    with open("concrete_building_adjoint_parameters.json",'r') as parameter_file:
        ProjectParametersAdjoint = Parameters( parameter_file.read())
    adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
    adjoint_analysis.Run()
```

#### Possible ```response_function_settings```:

Independet from the chosen response function the following definitions are always necessary:

- ```gradient_mode```: Currently there is only ```semi_analytic``` availible.
- ```step_size``` is the perturbation measure for finite difference computations within the semi-analytic apporach. 
- ```sensitivity_model_part_name```: Add here the name of the model part for which components sensitivities has to be computed (e.g. if the chosen design parameter is ```THICKNESS``` than for each element in this model part the sensitivity w.r.t. this variable is calculated).
- ```nodal_sensitivity_variables```: Currently only ```SHAPE_SENSITIVITY``` is availible. Doing this the sensitivities w.r.t. to the x-, y- and z-coordinate of all nodes in the ```sensitivity_model_part_name``` are computed.
- ```element_sensitivity_variables```: Here are sensitivities with respect to the properties of the elements computed. For that the respective name of the Kratos-Variable has to be given (e.g. ```THICKNESS```, ```I22``` or ```YOUNG_MODULUS```)
- ```condition_sensitivity_variables```: Here are sensitivities with respect to the properties of the conditions computed. Currenty there is only ```POINT_LOAD``` availible.

**Important, please note:**
In order to use an element or condition design variable one has to ensure that a corresponding Kratos-Variable for the sensitivity is defined (e.g. for the design variable ```THICKNESS``` a corresponding variable called ```THICKNESS_SENSITIVITY``` is necessary). This additional variable is necessary to store the results of the sensitivity analysis.

There are currently three different types of response functions availible which can be chosen as ```response_type```. For each of them specific settings are neccessary:
- ```adjoint_nodal_displacement```: The response is the displacement or rotation of a single node. Necessary additional settings are:
    * ```traced_node```: ID of the traced node
    * ```traced_dof```: Define the traced DOF (e.g. ```DISPLACEMENT_Z``` or ```ROTATION_X```)

- ```adjoint_strain_energy```: The response is the linear strain energy. No additional settings are necessary.

- ```adjoint_local_stress```: The response is the stress or stress-resultant of a single element. Necessary additional settings are:
    * ```traced_element```: ID of the element which should be traced
    * ```stress_type```: Stress type which should be traced (e.g. FX, MY, FXX or MXX)
    * ```stress_treatment```: There are three possibilities: ```node``` (Takes the response value at the position of a defined node of the element. Only availible for beam element.), ```GP``` (Takes the response value at a defined Gauss-Point of the element.) and ```mean``` (The response is the mean value of all Gauss-Point results of the traced stress type.)
    * ```stress_location```: Only necessary if ```node``` or ```GP``` is chosen as ```stress_treatment```. Define here the local ID of the position where the stress has to be traced (e.g. if the stress resultant of a beam element should be traced at one of his two nodes ```stress_location``` has to be 1 or 2)

Examples:
```python
    "response_function_settings" : {
            "response_type"     : "adjoint_nodal_displacement",
            "gradient_mode"     : "semi_analytic",
            "sensitivity_model_part_name" : "Parts_Beam",
            "nodal_sensitivity_variables"  : ["SHAPE_SENSITIVITY"],
            "element_sensitivity_variables"  : ["I22"],
            "condition_sensitivity_variables"  : ["POINT_LOAD"],
            "step_size"         : 1e-6,
            "traced_node"       : 6,
            "traced_dof"        : "DISPLACEMENT_Z"
        }
```   

```python
    "response_function_settings" : {
            "response_type"     : "adjoint_strain_energy",
            "gradient_mode"     : "semi_analytic",
            "sensitivity_model_part_name" : "Parts_Beam",
            "nodal_sensitivity_variables"  : ["SHAPE_SENSITIVITY"],
            "element_sensitivity_variables"  : ["I22"],
            "condition_sensitivity_variables"  : ["POINT_LOAD"],
            "step_size"         : 1e-6
        }
```   

```python
    "response_function_settings" : {
            "response_type"     : "adjoint_local_stress",
            "gradient_mode"     : "semi_analytic",
            "sensitivity_model_part_name" : "Parts_Beam",
            "nodal_sensitivity_variables"  : ["SHAPE_SENSITIVITY"],
            "element_sensitivity_variables"  : ["I22"],
            "condition_sensitivity_variables"  : ["POINT_LOAD"],
            "step_size"         : 1e-6,
            "traced_element"    : 6,
            "stress_type"       : "MY",
            "stress_treatment"  : "node",
            "stress_location"   : 1
        }
    }
```

### Post-Processing

The results of the sensitivity analysis are accessible in the post-processing as ```nodal_results``` (currently only ```SHAPE_SENSITIVITY``` and ```POINT_LOAD_SENSITIVITY```) and ```gauss_point_results``` (sensitivities for elemental design variables like ```THICKNESS_SENSITIVITY``` or ```I22_SENSITIVITY```).
```python
    "nodal_results"       : ["DISPLACEMENT","ADJOINT_DISPLACEMENT", "SHAPE_SENSITIVITY", "POINT_LOAD_SENSITIVITY"],
    "gauss_point_results" : ["THICKNESS_SENSITIVITY"]
```    







        
                
         



