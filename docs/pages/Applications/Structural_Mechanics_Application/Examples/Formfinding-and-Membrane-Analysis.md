---
title: Formfinding and Membrane Analysis
keywords: 
tags: [Formfinding-and-Membrane-Analysis.md]
sidebar: structural_mechanics_application
summary: 
---

This page descibes the simulation of prestressed membranes (isotropic and anisotropic) for analysis and formfinding which is implemented in the StructuralMechanicsApplication. The implentation is based on the following sources:
* [Numerical Methods for the Design and Analysis of Tensile Structures, Falko Hartmut Dieringer](https://mediatum.ub.tum.de/doc/1197480/880619.pdf)
*  [Mechanik und Numerik der Formfindung und Fluid-Struktur-Interaktion von Membrantragwerken, Roland Wüchner](https://mediatum.ub.tum.de/doc/601102/601102.pdf)

# Prestressed Membrane Element
The prestressed membrane element (`prestressed_membrane_element.h`) is designed in such a way that it is suitable for both formfinding and membrane analysis. The prestress can either be read from a .mdpa file (which could for example result from a formfinding analysis, example see following code)
```
Begin ElementalData MEMBRANE_PRESTRESS
1	[3,1]((1.0),(1.0),(1.0))
2	[3,1]((2.0),(1.0),(-1.0))
3	[3,1]((3.0),(1.0),(6.0))
4	[3,1]((5.0),(2.0),(1.0))
5	[3,1]((6.0),(3.0),(-5.0))
End ElementalData 
```
or from a projection (planar or rotational, see picture). The prestress, the projection type and the projection direction(s) are defined in the material. In the rotational case the variable PRESTRESS_AXIS_1_GLOBAL defines the axis of rotation. The entries in the PRESTRESS_VECTOR define the components for ["radial prestress", "tangential prestress", 0]:
```
{
    "properties" : [{
        "Material"        : {
            "Variables"        : {
                "PRESTRESS_VECTOR": [6000,5000,0]
                "PROJECTION_TYPE_COMBO": "planar"
                "PRESTRESS_AXIS_1_GLOBAL: [1 0 0],
                "PRESTRESS_AXIS_2_GLOBAL: [0 1 0]
}}}]}
```
```
{
    "properties" : [{
        "Material"        : {
            "Variables"        : {
                "PRESTRESS_VECTOR"      : [100,0,0],
                "PROJECTION_TYPE_COMBO" : "radial",
                "PRESTRESS_AXIS_1_GLOBAL" : [0,0,1],
}}}]}
```

![planar projection of prestress](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Formfinding_Membrane_Analysis/planar_rotational.png)

The projection of the prestress is executed by the method `ProjectPrestress`. 
The distinction between formfinding and membrane analysis is set by the value `IS_FORMFINDING`. If `IS_FORMFINDING = True`, the update reference strategy is applied where the base vectors are updated in each nonlinear iteration (Method: `ComputeBaseVectors`). In the case of prestress adaption, the prestress is updated in the method `UpdatePrestress`.

# Formfinding Analysis
The formfinding analysis is executed in `formfinding_updated_reference_strategy.hpp<´. In order to use the formfinding strategy the following parameters have to be included in the parameter input file:
```
"solver_settings"          : {
        "solver_type"                        : "Static",
        "echo_level"                         : 2,
        "analysis_type"                      : "formfinding",
        "line_search"                        : false,
        "print_formfinding_iterations"       : true
    }
```
In order to use the result of a formfinding analysis for further simulations, the necessary information (whole modelpart or prestress data) can be printed in an output file using the `formfinding_io_process`. The same process can be used in subsequent simulations to read the prestress data.
```
"list_other_processes": [
        {
            "python_module": "formfinding_IO_process",
            "kratos_module": "KratosMultiphysics.StructuralMechanicsApplication",
            "help": "This process is for input and output of prestress data",
            "process_name": "FormfindingIOProcess",
            "Parameters": {
                "model_part_name": "Structure",
                "print_mdpa": true,
                "print_prestress": false,
                "read_prestress": false
            }}]
```
# Formfinding example
The following pictures show an example where the formfinding is executed for a 4-point-sail. Isotropic prestress is applied by a planar projection. The first picture shows the original geometry and the second picture the formfinding result.
![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Formfinding_Membrane_Analysis/formfinding_original.png)![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/Application_cases/Formfinding_Membrane_Analysis/formfinding_result.png)