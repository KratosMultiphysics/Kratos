# Custom processes

This folder contains the custom processes that are used in the GeoMechanicsApplication. Processes are commands that can be executed during several stages of a calculation. In the GeoMechanicsApplication, their main use is to apply boundary conditions and/or constraints , but there are more specific examples (e.g. the K0 procedure, nodal extrapolation, or (de)activation of certain parts of the model). Since this effort is a work in progress, not all processes are documented yet. If you have any questions, please contact the maintainers (@KratosMultiphysics/geomechanics on GitHub).

Documented processes:
- [$c-\phi$ reduction process](#c-phi-reduction-process)
- [GeoExtrapolateIntegrationPointValuesToNodesProcess](#extrapolation-of-integration-values-to-nodes)
- [ApplyFinalStressesOfPreviousStageToInitialState](#apply-final-stresses-of-previous-stage-to-initial-state-process)
- [$K_0$ procedure process](#K_0-procedure-process)
- [ApplyInitialUniformStress](#apply-initial-uniform-stress)
- [FindNeighboursOfInterfaces](#find-neighbours-of-interfaces)

## $c-\phi$ reduction process
For the assessment of a safety factor to characterize slope stability, a Mohr-Coulomb material based $c-\phi$ reduction 
scheme is implemented. The apex of the Mohr-Coulomb cone shaped failure surface is kept in the same position, 
therefore both $c$ and $\tan \phi$ will diminish at the same rate.

### Incrementation scheme
The $c-\phi$ reduction process requires the existence of a stress state in your model and the use of a Mohr-Coulomb material 
(in a UDSM or UMAT formulation). Preferably this stress state is an equilibrium state, such that no stresses in integration 
points violate the given Mohr-Coulomb failure surface. During the stage with the active $c-\phi$ reduction process, 
$c$ and $\tan \phi$ will be incrementally reduced in steps with an initial size of 10%. For each reduction step stresses are 
mapped back onto the reduced Mohr-Coulomb yield surface and equilibrium is found if possible. When equilibrium is no longer 
found in the given number of iterations for the Newton-Raphson scheme, the step is retried with a halved reduction increment.
This is repeated until the allowed number of cycles for a step is reached or until the reduction factor drops below 0.01.
As the stepping is "mis"using the time-stepping of the GeoMechanics application, reaching the end_time also terminates the
reduction process.

### Safety factor
The safety factor is computed as the inverse of the reduction factor [[1]](#1):

$$SF = \frac{1}{\alpha}$$

where the reduction factor $\alpha$ is the ratio of the found critical values for cohesion and friction angle $c_c$ or $\phi_c$ and original material parameters.

$$\alpha = \frac{c_c}{c} = \frac{\tan \phi_c}{\tan \phi}$$ 

## Extrapolation of integration values to nodes 
The `GeoExtrapolateIntegrationPointValuesToNodesProcess` can be used as a post-processing step to acquire nodal data for variables that are stored at the integration points. This is useful for visualization services which expect nodal data.

Conceptually the process consists of the following steps:
1. Determine a count for each node, to keep track of how many elements will contribute to the nodal value.  Note that only the elements of the given model part(s) are considered!
2. Calculate the extrapolation matrix, to distribute the integration values to the nodes.
3. Calculate the integration point values of the variables of interest, by using the `CalculateOnIntegrationPoints` function of the `Element` class.
4. For each element, distribute the integration point values to their respective nodes by multiplying the extrapolation matrix with the integration point values.
5. Divide the nodal values by the count to get the average value.

### Limitations
The process supports floating-point scalar, vector, and matrix variables.  The supported element shapes include lines, triangles, quadrilaterals, tetrahedra, and hexahedra as well as line and plane interfaces.  Furthermore, the order of the element's shape functions must be either linear or quadratic.  The extrapolation is always done linearly.  For the quadratic elements, this means the corner nodes are extrapolated as usual, but the mid-side nodes are extrapolated using linear combinations of the extrapolation contributions for the corner nodes.

### Usage
The process is defined as follows in json (also found in some of the [integration tests](../tests/test_integration_node_extrapolation)):
```json
{
  "python_module": "geo_extrapolate_integration_point_values_to_nodes_process",
  "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
  "process_name":  "GeoExtrapolateIntegrationPointValuesToNodesProcess",
  "Parameters":    {
    "model_part_name":   "ModelPartName",
    "list_of_variables": ["Variable1", "Variable2", "Variable3"]
  }
}
```
The process receives either a single model part name (when using `model_part_name`) or a list of model part names (when using `model_part_name_list`).  Note that any elements that are not part of the given model part(s) are **not** considered by the extrapolation process.  Inactive elements are automatically discarded by the process.  In general, the variables that are to be extrapolated (supplied through `list_of_variables`) should be valid for all elements of the supplied model part(s), or else errors may occur.  For instance, extrapolation of bending moments only makes sense for structural elements that have curvatures.  Similarly, extrapolation of traction vectors only makes sense in the context of interface elements.  Therefore, it is recommended to have one extrapolation process per group of variables that can be calculated for all of the elements of the given model part(s).  As mentioned in the [Section Limitations](#limitations), these variables could be of any type, as long as the `Element` class has an implementation of the `CalculateOnIntegrationPoints` function for them.


When this process is added to the `ProjectParameters.json`, the variables specified in `list_of_variables` can be exported as nodal output (e.g. as `nodal_results` in the `GiDOutputProcess`). 

## Apply Final Stresses Of Previous Stage To Initial State Process
The `ApplyFinalStressesOfPreviousStageToInitialState` process can be used to change the reference point of the displacements to the displacement at the start of that stage. This process only needs to be applied to structural and interface elements, to convert the displacements from a total displacement to a staged displacement.

### Requirements
For this process to work, the following requirements have to be met:
1. The elements in the model part that the process is applied to should have an implementation for `CalculateOnIntegrationPoints` that calculates the PK2_STRESS_VECTOR as well as an overload of `CalculateOnIntegrationPoints` that returns a list of ConstitutiveLaw::Pointer objects for each integration point.
2. The ConstitutiveLaw used in the elements this process is applied to should use the `InitialState` to apply the initial stresses to the calculated stresses.


## $K_0$ procedure process
For the initialization of an in-situ stress field, the $K_0$ procedure derives the horizontal effective stresses from a field of vertical effective stresses.
Pre-requisite is a computed stress field with the desired normal effective stresses in the direction indicated with "K0_MAIN_DIRECTION".
The normal effective stress in "K0_MAIN_DIRECTION" remains as is.
Effective normal stresses in the other two directions are affected by the $K_0$ value, all shear stresses are erased.
A specialized method for computation of the stress field with normal effective stresses is steered with "use_standard_procedure".
When set to true, the material used for the modelpart used for the $K_0$ procedure process is changed to an incremental elastic material with constitutive tensor that has only zero off-diagonal terms and zero shear terms.
When the stress computation is completed, the original constitutive law is restored.

Depending on the given input parameters, the following scheme is adapted for computation of the $K_0$ value.
$K_0^{nc}$ is gotten from either "K0_NC" the material input file or by computation from input of "INDEX_OF_UMAT_PHI_PARAMETER" and "UMAT_PARAMETERS":

$$K_0^{nc} = 1.0 - \sin \phi$$

When the overconsolidation ratio ("OCR") and optionally the unloading-reloading Poisson's ratio $\nu_{ur}$ ("POISSON_UNLOADING_RELOADING") are supplied, the normal consolidation value $K_0^{nc}$ is modified:

$$K_0 = OCR.K_0^{nc} +  \frac{\nu_{ur}}{1 - \nu_{ur}} ( OCR - 1 )$$

```math
\sigma^{'}_{initial} = \begin{bmatrix} {K_0 \sigma^{'}_{zz}} & 0 & 0 \\
                                        0 & {K_0 \sigma^{'}_{zz}} & 0 \\
                                        0 & 0 & {\sigma^{'}_{zz}} \end{bmatrix}
```

Alternatively, when the pre-overburden pressure "POP" is specified, the initial stress tensor becomes:

```math
\sigma^{'}_{initial} = \begin{bmatrix} {K_0^{nc} (\sigma^{'}_{zz} + POP )} - \frac{\nu_{ur}}{1 - \nu_{ur}} POP & 0 & 0 \\
                                        0 & {K_0^{nc} (\sigma^{'}_{zz} + POP)} - \frac{\nu_{ur}}{1 - \nu_{ur}} POP & 0 \\
                                        0 & 0 & {\sigma^{'}_{zz}} \end{bmatrix}
```

When the optional unloading-reloading Poisson's ratio is omitted, a default value $\nu_{ur} = 0$ is used such that the correction term drops to 0.

### Note:
After the stress adaptation by the $K_0$ procedure, the stress state may not be in equilibrium with the present external forces anymore. Equilibrium may be reached by performing a step without applying additional load. Reaching equilibrium may then be accomplished by movement.

### Usage
The process is defined as follows in "ProjectParameters.json" (also found in some of the [integration tests](../tests/test_k0_procedure_process)). Without the addition of this process, no adaptation of the horizontal stresses takes place.
```json
{
  "auxiliary_process_list": [
    {
      "python_module": "apply_k0_procedure_process",
      "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
      "process_name": "ApplyK0ProcedureProcess",
      "Parameters": {
        "model_part_name": "PorousDomain.porous_computational_model_part",
        "use_standard_procedure": true
      }
    }
  ]
}
```
Next to specifying a single model part, it is also possible to provide a list:
```json
{
  "auxiliary_process_list": [
    {
      "python_module": "apply_k0_procedure_process",
      "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
      "process_name": "ApplyK0ProcedureProcess",
      "Parameters": {
        "model_part_name_list": ["PorousDomain.Clay", "PorousDomain.Sand"],
        "use_standard_procedure": true
      }
    }
  ]
}
```

The "apply_k0_procedure_process" needs the following material parameter input to be added in the "MaterialParameters.json".
```json
{
  "Variables": {
    "K0_MAIN_DIRECTION":           1,
    "K0_NC":                       0.6,
    "UDSM_NAME"                :  "MohrCoulomb64.dll",
    "IS_FORTRAN_UDSM"          :  true,
    "INDEX_OF_UMAT_PHI_PARAMETER": 4,
    "UMAT_PARAMETERS"          :  [30000000,
                                   0.2,
                                   1000.0,
                                   30,
                                   0.0,
                                   1000],
    "OCR":                         1.4,
    "POISSON_UNLOADING_RELOADING": 0.35,
    "POP":                         800.0
  }
}
```

## Apply Initial Uniform Stress

This process applies an initial uniform stress field to all elements in a model part.
The elements in the model part need to be able to calculate and set the CAUCHY_STRESS_VECTOR variable.
The Parameters object should contain a "value" field, which is a vector representing the stress components.
The vector should have a length equal to the strain size (e.g. 4 for plane strain and axisymmetric cases, 6 for 3D).
Note that this means that if you want to apply a uniform stress field to 
elements with different strain sizes, you will need to apply the process multiple times with separate model parts.

Example usage for a case with 3D elements in a ProjectParameters.json file:
```json
{
    "loads_process_list": [
      {
        "python_module": "apply_initial_uniform_stress_field",
        "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
        "process_name":  "ApplyInitialUniformStressField",
        "Parameters":    {
        "model_part_name": "PorousDomain.Soil",
        "value": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        }
      }
    ]
}
```


## Find Neighbours Of Interfaces

This process finds the neighbouring elements of interface elements in a model part. These neighbours are then used to calculate and apply a prestress to the interfaces based on the stress state of the neighbouring elements. Typically, this process is used in a multi-stage analysis, where in a specific stage the interfaces are installed (along with a structural element that models, for instance, a sheet pile wall). To avoid deformations due to differences in stress between already existing soil elements and the newly added interface elements, equilibrium is ensured by prestressing the interfaces using the stresses of the surrounding soil.

The process of applying prestress to the interfaces consists of the following steps:
1. The neighbouring elements of the interface elements are found using this process.
2. The stresses at the integration points of the neighbouring elements are extrapolated to their respective nodes.
3. The nodal stresses are interpolated to the integration points of the interface elements.

Note that steps 2 and 3 are not part of `FindNeighboursOfInterfacesProcess`, but they are taken care of by the interface element itself when neighbours are known. The process only finds neighbouring elements with a higher local dimension than the interface elements, to avoid prestressing the element with stresses of non-continuum elements (e.g. structural elements or other interface elements). 

Example usage for a case in a ProjectParameters.json file:

```json
{
  "python_module": "find_neighbours_of_interfaces_process",
  "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
  "process_name": "FindNeighboursOfInterfacesProcess",
  "Parameters": {
    "model_part_name": "PorousDomain.Interface",
    "model_part_name_for_neighbouring_elements": "PorousDomain.porous_computational_model_part"
  }
}
```

The `model_part_name_for_neighbouring_elements` is used to specify the model part that contains the elements which can be neighbours to the interface elements. Typically, this is the entire computational domain.

Next to specifying a single model part, it is also possible to provide a list of model parts containing interface elements:
```json
{
  "auxiliary_process_list": [
    {
      "python_module": "find_neighbours_of_interfaces_process",
      "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
      "process_name": "FindNeighboursOfInterfacesProcess",
      "Parameters": {
        "model_part_name_list": ["Interfaces_Left", "Interfaces_Right"],
        "model_part_name_for_neighbouring_elements": "PorousDomain.porous_computational_model_part"
      }
    }
  ]
}
```
## References
<a id="1">[1]</a> Brinkgreve, R.B.J., Bakker, H.L., 1991. Non-linear finite element analysis of safety factors, Computer Methods and Advances in Geomechanics, Beer, Booker & Carterr (eds), Balkema, Rotterdam.
