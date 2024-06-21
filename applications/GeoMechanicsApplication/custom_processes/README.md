# Processes
This folder contains the custom processes that are used in the GeoMechanicsApplication. Processes are commands that can be executed during several stages of a calculation. In the GeoMechanicsApplication, their main use is to apply boundary conditions and/or constraints , but there are more specific examples (e.g. the K0 procedure, nodal extrapolation, or (de)activation of certain parts of the model).

In this document, we will explain the different processes. Since this effort is a work in progress, not all processes are documented yet. If you have any questions, please contact the mainainers (@KratosMultiphysics/geomechanics on GitHub).

Documented processes:
- [GeoIntegrationValuesExtrapolationToNodesProcess](#extrapolation-of-integration-values-to-nodes)

## Extrapolation of integration values to nodes
The `GeoIntegrationValuesExtrapolationToNodesProcess` can be used as a post-processing step to acquire nodal data for variables that are stored at the integration points. This is useful for visualization services which expect nodal data.

Conceptually the process consists of the following steps:
1. Determine a count for each node, to keep track of how many elements will contribute to the nodal value.
2. Calculate the extrapolation matrix, to distribute the integration values to the nodes.
3. Calculate the integration point values of the variables of interest, by using the `CalculateOnIntegrationPoints` function of the `Element` class.
4. For each element, distribute the integration values to their respective nodes by multiplying the extrapolation matrix with the integration values.
5. Divide the nodal values by the count to get the average value.

### Restrictions
Currently, this process is only implemented for 3-noded or 6-noded `Triangle` and 4-noded or 8-noded `Quadrilateral` elements in 2D. The extrapolation is always done linearly. For the higher order 6-noded and 8-noded elements, this means the corner nodes are extrpolated as usual, but the mid-side nodes are extrapolated using linear combinations of the extrapolation contributions for the corner nodes.

### Usage
The process is defined as follows in json (also found in some of the [integration tests](../tests/test_integration_node_extrapolation)):
```json
{
  "process_name": "GeoIntegrationValuesExtrapolationToNodesProcess",
  "Parameters":   {
    "model_part_name":   "ModelPartName",
    "list_of_variables": ["Variable1", "Variable2", "Variable3"],
    "average_variable":  "NODAL_AREA"
  }
}
```
Where the model_part_name should contain the name of the model part where the extrapolation is to be performed for the variables in `list_of_variables`. These variables could be of any type, as long as the `Element` class has an implementation of the `CalculateOnIntegrationPoints` function for them.

The `average_variable` is used to divide the nodal values by the count, and it is recommended to use `NODAL_AREA` here.

When this process is added to the `ProjectParameters.json`, the variables specified in `list_of_variables` can be exported as nodal output (e.g. as `nodal_results` in the `GiDOutputProcess`). 

