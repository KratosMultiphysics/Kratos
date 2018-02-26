<img src="https://github.com/KratosMultiphysics/Kratos/tree/Brep_Application/applications/NurbsBrepApplication/readme_application_description/nurbs_symbol.png">

# Nurbs Brep Application

The application is designed to handle CAD designed geometries inside KRATOS. It is enhanced with several CAD-functionalities, which are needed for simulation.

## Getting Started

The application is part of Kratos Multiphysics. Instructions how to install Kratos for development and testing purposes are available for both [Linux](http://kratos-wiki.cimne.upc.edu/index.php/LinuxInstall) and [Windows](http://kratos-wiki.cimne.upc.edu/index.php/Windows_7_Download_and_Installation) systems. The application is available for both platforms.

Please add following instructions to the configure file for the correct compilation:
``` cmake
-DNURBS_BREP_APPLICATION=ON  ^
```
If it is wished to do structural mechanics simulation with Isogeometric Analysis the *IGAStructuralMechanicsApplication* provides necessary element and condition formulations for thin walled structures.

## Preparation of Geometries
Generally the geometries are stored in the geometry.json file. This file is a standardized B-Rep based NURBS geometry file. It is recommended to use this standard. Anyways read in of geometries with different input formats is possible, too.
```
#import define_output
cad_geometry_file = open("geometry.json",'r')
cad_geometry = Parameters( cad_geometry_file.read())
```
The json file can be read by the *BrepModelGeometryReader*. It generates a structure to handle full CAD simulation possibilities.
```
geometry_reader = BrepModelGeometryReader(cad_geometry)
```
The geometry structure is handled by the *NurbsBrepModeler*. Here, the creation of the Integration Domain, projections towards the surface and additional processes are possible.
```
modeler = NurbsBrepModeler(geometry_reader, model_part)
```
## Generation of Integration Domain
The application can create a proper Isogeometric B-Rep Analysis integration domain for numerical simulation.

After reading in the geometries, the generation of the integration domain can be called inside python with the following command:
```
modeler.CreateIntegrationDomain(order, model_part_integration_domain)
```
It has to be passed a model part which is used as data holder for the integration points. Additionally the wished number of derivatives of shape functions can be given. If it is chosen as -1, no shape functions are enhanced. With the order being 0 only the shape functions are given.

[Creation of integration domain]<img src="https://github.com/KratosMultiphysics/Kratos/tree/Brep_Application/applications/NurbsBrepApplication/readme_application_description/integration_domain.png">

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

## Project parameter set up

The application holds the possibility to specify special parameters for the geometrical interpretation. 
The possible parameters are mentioned in the section:
``` 
"nurbs_geometry_configuration": {
  ...
}
``` 

### Description of integration domain
```
"integration_domain_parameter": {
  "integration_scheme": "triangulation", //"AGIP"
  "accuracy": 10e-8,
  "number_of_projection_iterations": 100,
  "integration_domains": {
    "faces": true,
    "faces_embedded": false,
    "faces_reversed": false,
    "edges": true,
    "coupling_edges": true
  }
},
```
### Geometry refinement
If it is wished to manually refine the geometries that were obtained from CAD those refinement operations can be specified in the following section. Here, either all patches or specific patches can be selected.
The order of those operations is line-wise.
```
"refinement": [
  {
    "selection": "ALL", //"SELECTION"
    // only if "SELECTION": "patch_id": [int, ...],
    "parameters": {
      "knot_insertions_u": [ int ... ],
      "knot_insertions_v": [ 1, 2 ],
      "multiply_knots_u": 4,
      "multiply_knots_v": 4,
      "max_element_size_u": 0.3,
      "max_element_size_v": 0.3,
      "order_elevation_p": 1,
      "order_elevation_q": 1,
      "min_order_p": 1,
      "min_order_q": 1
    }
  }
]
```