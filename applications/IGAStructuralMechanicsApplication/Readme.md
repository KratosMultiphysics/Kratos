<img src="https://github.com/KratosMultiphysics/Kratos/tree/Brep_Application/applications/IGAStructuralMechanicsApplication/readme_application_description/IGA_symbol.png">

# IGA Structural Mechanics Application

The Application handles a couple of necessary element and condition formulations for structural mechanics simulations optimized for the use of NURBS-basis functions. The application is completely based on meshless formulations, thus additionally it could be used with standard FEM shape functions.

## Getting Started

The application is part of Kratos Multiphysics. Instructions how to install Kratos for development and testing purposes are available for both [Linux](http://kratos-wiki.cimne.upc.edu/index.php/LinuxInstall) and [Windows](http://kratos-wiki.cimne.upc.edu/index.php/Windows_7_Download_and_Installation) systems. The application is available for both platforms.

Please add following instructions to the configure file for the correct compilation:
``` cmake
-DIGA_STRUCTURAL_MECHANICS_APPLICATION=ON  ^
```
It is recommended to use the *NurbsBrepApplication* for preparation of geometry and generation of the numerical integration domain.

## Meshless formulations
Meshless formulations require a precomputation of numerical data. This means in general, the definition of shape functions, shape function derivatives and the integration weight. Some conditions might need extra information. This information can be obtained by the use of the *NurbsBrepApplication*.

## Element formulations
The application is optimized for the use of thin-walled structures, but can be enhanced with other formulations. Hence the conditions are formulated for Kirchhoff-Love shell (see [#Kiendl2011]). The structure of the application is 

[#Kiendl2011]: http://nbn-resolving.de/urn/resolver.pl?urn:nbn:de:bvb:91-diss-20110321-1002634-1-5   "Isogeometric Analysis and Shape Optimal Design of Shell Structures."
