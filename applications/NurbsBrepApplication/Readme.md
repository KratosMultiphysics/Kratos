<img src="https://github.com/KratosMultiphysics/Kratos/tree/Brep_Application/applications/NurbsBrepMechanicsApplication/readme_application_description/nurbs_symbol.png">

# Nurbs Brep Application

The application is designed to handle CAD designed geometries inside KRATOS. It is enhanced with several CAD-functionalities, which are needed for simulation.

## Getting Started

The application is part of Kratos Multiphysics. Instructions how to install Kratos for development and testing purposes are available for both [Linux](http://kratos-wiki.cimne.upc.edu/index.php/LinuxInstall) and [Windows](http://kratos-wiki.cimne.upc.edu/index.php/Windows_7_Download_and_Installation) systems. The application is available for both platforms.

Please add following instructions to the configure file for the correct compilation:
``` cmake
-DNURBS_BREP_APPLICATION=ON  ^
```
If it is wished to do structural mechanics simulation with Isogeometric Analysis the *IGAStructuralMechanicsApplication* provides necessary element and condition formulations for thin walled structures.

## Generation of Integration Domain
The application can create a proper Isogeometric B-Rep Analysis integration domain for numerical simulation.

After reading in the geometries, the generation of the integration domain can be called inside python with the following command:
```
modeler.CreateIntegrationDomain(2, model_part_integration_domain)
```
It has to be passed a model part which is used as data holder for the integration points. Additionally the wished number of derivatives of shape functions can be given. If it is chosen as -1, no shape functions are enhanced.
<img src="https://github.com/KratosMultiphysics/Kratos/tree/Brep_Application/applications/NurbsBrepApplication/readme_application_description/integration_domain.png">

## Meshless formulations
Meshless formulations require a precomputation of numerical data. This means in general, the definition of shape functions, shape function derivatives and the integration weight. Some conditions might need extra information.
The application can provide those entities. Thus, afterwards no NURBS-functionalities are needed anymore.

The outcome of the application is structured as follows. For each different topology item different data is necessary:
```
**SURFACE_POINT**
INTEGRATION_WEIGHT
SHAPE_FUNTION_VALUES
SHAPE_FUNCTION_LOCAL_DERIVATIVES
```
Edge points for application of weak conditions contain additionally the directions of tangents:
```
**EDGE_POINT**
INTEGRATION_WEIGHT
SHAPE_FUNTION_VALUES
SHAPE_FUNCTION_LOCAL_DERIVATIVES
TANGENTS
```
For the application of continuities the shape functions of both patches are necessary:
```
**COUPLING_EDGE_POINT**
INTEGRATION_WEIGHT
SHAPE_FUNTION_VALUES
SHAPE_FUNCTION_LOCAL_DERIVATIVES
TANGENTS
SHAPE_FUNTION_VALUES_SLAVE
SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE
TANGENTS_SLAVE
```