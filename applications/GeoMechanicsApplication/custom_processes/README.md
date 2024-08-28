# Custom processes

This folder contains the custom processes that are used in the GeoMechanicsApplication. Processes are commands that can be executed during several stages of a calculation. In the GeoMechanicsApplication, their main use is to apply boundary conditions and/or constraints , but there are more specific examples (e.g. the K0 procedure, nodal extrapolation, or (de)activation of certain parts of the model). Since this effort is a work in progress, not all processes are documented yet. If you have any questions, please contact the maintainers (@KratosMultiphysics/geomechanics on GitHub).

Documented processes:
- $c-\phi$ reduction process
- [GeoExtrapolateIntegrationPointValuesToNodesProcess](#extrapolation-of-integration-values-to-nodes)

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
found in the given number of iterations for the Newton-Raphson scheme, the step is retried with a halved reduction increment. This is repeated until the allowed number of cycles for a step is reached.   

### Safety factor
The safety factor is computed as the inverse of the reduction factor [[1]](#1):

$$SF = \frac{1}{\alpha}$$

where the reduction factor $\alpha$ is the ratio of the found critical values for cohesion and friction angle $c_c$ or $\phi_c$ and original material parameters.

$$\alpha = \frac{c_c}{c} = \frac{\tan \phi_c}{\tan \phi}$$ 

## Extrapolation of integration values to nodes 
The `GeoExtrapolateIntegrationPointValuesToNodesProcess` can be used as a post-processing step to acquire nodal data for variables that are stored at the integration points. This is useful for visualization services which expect nodal data.

Conceptually the process consists of the following steps:
1. Determine a count for each node, to keep track of how many elements will contribute to the nodal value.
2. Calculate the extrapolation matrix, to distribute the integration values to the nodes.
3. Calculate the integration point values of the variables of interest, by using the `CalculateOnIntegrationPoints` function of the `Element` class.
4. For each element, distribute the integration point values to their respective nodes by multiplying the extrapolation matrix with the integration point values.
5. Divide the nodal values by the count to get the average value.

### Restrictions
Currently, this process is only implemented for 3-noded or 6-noded `Triangle` and 4-noded or 8-noded `Quadrilateral` elements in 2D. The extrapolation is always done linearly. For the higher order 6-noded and 8-noded elements, this means the corner nodes are extrapolated as usual, but the mid-side nodes are extrapolated using linear combinations of the extrapolation contributions for the corner nodes.

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
Where the `model_part_name` should contain the name of the model part where the extrapolation is to be performed for the variables in `list_of_variables`. These variables could be of any type, as long as the `Element` class has an implementation of the `CalculateOnIntegrationPoints` function for them.


When this process is added to the `ProjectParameters.json`, the variables specified in `list_of_variables` can be exported as nodal output (e.g. as `nodal_results` in the `GiDOutputProcess`). 

## $K0$ procedure process
For the initialization of an in-situ stress field, the $K0$ procedure derives the horizontal effective stresses from a field of vertical effective stresses.
Pre-requisite is a computed stress field with the desired normal effective stresses in the direction indicated with "K0_MAIN_DIRECTION". The normal effective stress in "K0_MAIN_DIRECTION" remains as is. Effective normal stresses in the other two directions are affected by the $K_0$ value, all shear stresses are erased.


Depending on the given input parameters, the following scheme is adapted for computation of the $K0$ value.
$K0_{NC}$ is gotten from either "K0_NC" the material input file or by computation from input of "INDEX_OF_UMAT_PHI_PARAMETER" and "UMAT_PARAMETERS"

$$K0_{NC} = 1.0 - \sin \phi$$

When "OCR" and optionally "POISSON_UNLOADING_RELOADING" are supplied, the normal consolidation value $K0_NC$ is modified:

$$K0 = OCR \cdot K0_{NC} +  \frac{\nu_{ur}}{1 - \nu_{ur}} ( OCR - 1 )$$

$$\sigma^{'}_{initial} = \begin{bmatrix} K0 \cdot \sigma^{'}_{zz} & 0 & 0 \\ 0 & K0 \cdot \sigma^{'}_{zz} & 0 \\ 0 & 0 & \sigma^{'}_{zz} \end{bmatrix}$$

Alternaternively, when the pre-overburden pressure "POP" is specified, the initial stress tensor becomes:

$$\sigma^{'}_{initial} = \begin{bmatrix} K0_{NC} \cdot (\sigma^{'}_{zz} + POP ) & 0 & 0 \\ 0 & K0_{NC} \cdot (\sigma^{'}_{zz} + POP) & 0 \\ 0 & 0 & \sigma^{'}_{zz} \end{bmatrix}$$

###Note:
After the stress adaptation by the $K0$ procedure, the stress state may not be in equilibrium with the present external forces anymore. Equilibrium may be reached by performing a step without applying additional load. Reaching equilibrium may then be accomplished by movement.  




### Usage
The process is defined as follows in json (also found in some of the [integration tests](../tests/test_k0_procedure_process)):
```json
{
  "auxilliary_process_list": [
    {
      "python_module": "apply_k0_procedure_process",
      "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
      "process_name": "ApplyK0ProcedureProcess",
      "Parameters": {
        "model_part_name": "PorousDomain.porous_computational_model_part",
        "variable_name": "CAUCHY_STRESS_TENSOR"
      }
    }
  ]
}
```
The apply_k0_procedure_process need the following material parameter input.
```json
{
  "Variables": {
    "K0_MAIN_DIRECTION":           1,
    "K0_NC":                       0.6,
    "UDSM_NAME"                :  "MohrCoulomb64.dll",
    "IS_FORTRAN_UDSM"          :  true,
    "NUMBER_OF_UMAT_PARAMETERS":  6,
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
  },
}
```

## References
<a id="1">[1]</a> Brinkgreve, R.B.J., Bakker, H.L., 1991. Non-linear finite element analysis of safety factors, Computer Methods and Advances in Geomechanics, Beer, Booker & Carterr (eds), Balkema, Rotterdam.
