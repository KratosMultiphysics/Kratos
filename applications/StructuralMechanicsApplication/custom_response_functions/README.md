
## Adjoint Sensitivity Analysis

### General remarks:

This feature provides the framework to compute sensitivities of structural responses (e.g. displacements, strain energy or stresses) with respect to different types of design variables (e.g. nodal coordinates, material or cross-sectional properties or load intensity) with the adjoint approach. An adjoint system has to be solved for each responsefunciton. The sensitivities are then computed in a post-processing step.

*Please note:*
- This feature currently only works for linear problems (exception: geometric nonlinear problems exclusively modeled with TrussElement3D2N)
- This feature makes use of the HDF5Application

### Features:

Provided is a **scheme** to solve the adjoint problem and a **replacement process** which replaces all elements and conditions of a model with its adjoint equivalents and vice versa.



&nbsp;

| Response Function | Semi-analytic | Analytic | Static | Transient |
| ----------------- |  ---------------- | -------- | ------- | ----------|
|  Linear strain energy|    x   |          |      x   |         |
|  Displacement or rotation of a node      |       |    x    |      x    |         |
|  Stress resultant of a single element |    x   |          |      x   |         |
*In the table above, the first two columns (Semi-analytic, Analytic) refer to the way, the partial  derivative of the response is calculated.*


&nbsp;

| Structural Condition | Adjoint Condition | Design Variables | Semi-analytic | Analytic | Static | Transient |
| -------------------- | ----------------- | ------------------- |  ---------------- | -------- | ------- | ----------|
| PointLoadCondition | PointLoadAdjointCondition¹ | POINT_LOAD |    x   |          |      x   |         |
|                      |                             | SHAPE     |    x   |          |      x    |     |



&nbsp;

| Structural Element | Adjoint Element | Design Variables | Semi-analytic | Analytic | Static | Transient |
| -------------------- | ----------------- | ------------------- |  ---------------- | -------- | ------- | ----------|
| ShellThinElement3D3N | AdjointFiniteDifferencingShellElement¹ | THICKNESS² |    x   |          |      x   |         |
|                      |                             | SHAPE     |    x   |          |      x    |
| CrBeamElementLinear3D2N | AdjointFiniteDifferenceCrBeamElement¹ | I22²|    x   |          |      x   |         |
|                      |                             | SHAPE     |    x   |          |      x   |          |
| TrussElement3D2N | AdjointFiniteDifferenceTrussElement¹ | CROSS_AREA²|    x   |   x³       |      x   |         |
|                      |                             | SHAPE     |    x   |          |      x   |          |
| TrussLinearElement3D2N | AdjointFiniteDifferenceTrussLinearElement¹ | CROSS_AREA²|    x   |         |      x   |         |
|                      |                             | SHAPE     |    x   |          |      x   |          |

¹ The adjoint elements and conditions wrap elements/conditions of the Structural Mechanics Application and can call its public functions.  The main task of the adjoint elements/conditions is to derive different quantities (e.g. the right hand side or post-processing results like stresses) with respect to the design variable or state.

² In principle the adjoint elements are able to compute sensitivities w.r.t. to all properties which are available at a specific element. ```THICKNESS``` and ```I22``` are possible examples. One only has to ensure that a corresponding Kratos-Variable for the sensitivity result is defined (e.g. for the design variable ```THICKNESS``` a corresponding variable called ```THICKNESS_SENSITIVITY``` is necessary). This additional variable is necessary to store the results of the sensitivity analysis. See also the paragraph Post-Processing.

³ the nonlinear adjoint truss element is able to compute its stress-displacement-derivative analytically. All other derivatives are computed by finite-difference approximation on element level.


### Usage:
In order to perform a sensitivity analysis for one response function, the solution of two problems is necessary: The primal and the adjoint problem.

*Please note:*
For the solution of the two problems different kind of variables are used in order to store the results. For the primal problem the usual variables ```DISPLACEMENT``` and ```ROTATION``` and for the adjoint problem ```ADJOINT_DISPLACEMENT``` and ```ADJOINT_ROTATION``` are used.

#### Definition of the Primal Problem
The primal problem can be defined by the regular input files which are needed for standard analysis. As only difference the output process of the HDF5Application has to be added in the project parameters. An example for a possible file is:
*StructuralMechanicsApplication/tests/adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_parameters.json*


#### Definition of the Adjoint Problem
In order to define the adjoint problem an additional \*.json-file for the adjoint project parameters is necessary. This input file is in principle very similar to the respective file of the primal analysis. In comparison to a regular file for a linear static analysis three points have to be modified:
- ```solver_settings``` by using the ```adjoint``` as ```solver_type``` and by the definition of the ```response_function_settings```
- The input process of the HDF5Application has to be added in order to read the primal solution
- When defining *Dirichlet* conditions ```ADJOINT_DISPLACEMENT``` respective ```ADJOINT_ROTATION``` instead of ```DISPLACEMENT``` respective ```ROTATION``` have to be used.

An example for a possible input file is:
*StructuralMechanicsApplication/tests/adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_local_stress_adjoint_parameters.json*

#### How to run the analysis:

If all necessary input files are defined, the analysis can be performed with a simple python script by calling the primal and afterwards the adjoint analysis. The sensitivities are computed in a post-processing of the adjoint problem.

A possible python workflow can be found in:
*StructuralMechanicsApplication/tests/test_adjoint_sensitity_analysis_beam_3d2n_structure.py*

#### Possible ```response_function_settings```:

Independent from the chosen response function the following definitions are always necessary:

- ```gradient_mode```: Currently there is only ```semi_analytic``` available.
- ```step_size``` is the perturbation measure for finite difference computations within the semi-analytic approach.
- ```sensitivity_model_part_name```: The name of the model part for which components sensitivities have to be computed (e.g. if the chosen design parameter is ```THICKNESS``` then for each element in this model part the sensitivity w.r.t. this variable is calculated).
- ```nodal_sensitivity_variables```: Currently only ```SHAPE``` is available. Doing this the sensitivities w.r.t. to the x-, y- and z-coordinate of all nodes in the ```sensitivity_model_part_name``` are computed.
- ```element_sensitivity_variables```: Here sensitivities with respect to the properties of the elements are computed. For that the respective name of the Kratos-Variable has to be given (e.g. ```THICKNESS```, ```I22``` or ```YOUNG_MODULUS```)
- ```condition_sensitivity_variables```: Here sensitivities with respect to the properties of the conditions are computed. Currently there is only ```POINT_LOAD``` available.

**Important, please note:**
In order to use an element or condition design variable one has to ensure that a corresponding Kratos-Variable for the sensitivity is defined (e.g. for the design variable ```THICKNESS``` a corresponding variable called ```THICKNESS_SENSITIVITY``` is necessary). This additional variable is necessary to store the results of the sensitivity analysis.

There are currently three different types of response functions available which can be chosen as ```response_type```. For each of them specific settings are necessary:
- ```adjoint_nodal_displacement```: The response is the displacement or rotation of a single node. Necessary additional settings are:
    * ```traced_node_id```: ID of the traced node
    * ```traced_dof```: Define the traced DOF (e.g. ```DISPLACEMENT_Z``` or ```ROTATION_X```)

- ```adjoint_linear_strain_energy```: The response is the linear strain energy. No additional settings are necessary.

- ```adjoint_local_stress```: The response is the stress or stress-resultant of a single element. Necessary additional settings are:
    * ```traced_element_id```: ID of the element which should be traced
    * ```stress_type```: Stress type which should be traced (e.g. FX, MY, FXX or MXX)
    * ```stress_treatment```: There are three possibilities:
        * ```node``` (Takes the response value at the position of a defined node of the element. Only available for beam element.),
        * ```GP``` (Takes the response value at a defined Gauss-Point of the element.)
        * ```mean``` (The response is the mean value of all Gauss-Point results of the traced stress type.)
    * ```stress_location```: Only necessary if ```node``` or ```GP``` is chosen as ```stress_treatment```. Local index of the gauss point or node where the stress should be traced respectively. E.g. if the stress resultant of a beam element should be traced at the first of the two nodes ```stress_location``` has to be 1)

Exemplary input files:
- StructuralMechanicsApplication/tests/adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_local_stress_adjoint_parameters.json
- StructuralMechanicsApplication/tests/adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_nodal_disp_adjoint_parameters.json
- StructuralMechanicsApplication/tests/adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_strain_energy_adjoint_parameters.json



### Post-Processing

The results of the sensitivity analysis are accessible in the post-processing as ```nodal_results``` (currently only ```SHAPE_SENSITIVITY``` and ```POINT_LOAD_SENSITIVITY```) and ```gauss_point_results``` (sensitivities for elemental design variables like ```THICKNESS_SENSITIVITY``` or ```I22_SENSITIVITY```).
```python
    "nodal_results"       : ["ADJOINT_DISPLACEMENT", "SHAPE_SENSITIVITY", "POINT_LOAD_SENSITIVITY"],
    "gauss_point_results" : ["THICKNESS_SENSITIVITY"]
```













