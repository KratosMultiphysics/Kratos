<img src="https://github.com/KratosMultiphysics/Kratos/tree/Brep_Application/applications/IGAStructuralMechanicsApplication/readme_application_description/IGA_symbol.png">

# IGA Structural Mechanics Application

The Application handles a couple of necessary element and condition formulations for structural mechanics simulations. The application is completely based on meshless formulations, thus additionally it could be used with standard FEM shape functions.

## Meshless formulations
Meshless formulations require a precomputation of numerical data. This means in general the definition of shape functions, shape function derivatives and the integration weight. Some conditions might extra information. This information can be obtained by the use of the NurbsBrepApplication.

## Element formulations
The application is optimized for the use of thin-walled structures, but can be enhanced with other formulations. Hence the conditions are formulated towards a Kirchhoff-Love shell.