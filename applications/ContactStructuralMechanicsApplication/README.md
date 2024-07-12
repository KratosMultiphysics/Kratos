# Contact Structural Mechanics Application 
 
 |             **Application**             |                                                                                    **Description**                                                                                    |                              **Status**                              | **Authors** |
|:---------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------:|:-----------:|
| `ContactStructuralMechanicsApplication` | The *Contact Structural Mechanics Application* contains the contact mechanics implementations that can be used by the *Structural Mechanics Application* and *Constitutive Laws Application* within *Kratos Multiphysics* | <img src="https://img.shields.io/badge/Status-%F0%9F%94%A7Maintained-blue"  width="300px"> | [*Vicente Mataix Ferrándiz*](mailto:vmataix@altair.com)  <br /> [*Alejandro Cornejo Velázquez*](mailto:acornejo@cimne.upc.edu)  |
 
<p align="center">
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/validation/double_arch/data/result.gif" alt="Solution" style="width: 300px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/validation/double_arch/data/result_frictional.gif" alt="Solution" style="width: 300px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/use_cases/in_ring/data/animation.gif" alt="Solution" style="width: 300px;"/>
 <img src="https://github.com/KratosMultiphysics/Examples/raw/master/contact_structural_mechanics/use_cases/hyperelastic_tubes/data/half_cylinders.gif" alt="Solution" style="width: 300px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/mmg_remeshing_examples/use_cases/contacting_cylinders/data/nodal_h.gif" alt="Solution" style="width: 300px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/use_cases/self_contact/data/animation.gif" alt="Solution" style="width: 300px;"/>
</p>

The application includes tests to check the proper functioning of the application.
 
## 😎 Features: 
 
- *Mesh tying conditions based in mortar formulation*
 
- *Augmented Lagrangian contact conditions based in mortar formulation*
 
    * Frictionless formulation

    * Frictional formulation

- *Penalty contact conditions based in mortar formulation*

     * Frictionless formulation

     * Frictional

- *Simplified **MPC** conditions based in mortar formulation. With the mortar formulation the weight are computed, allowing to compute a Simplified *NTN* and a simplified *NTS**

     * Frictionless formulation

     * Frictional formulation

     * Mesh tying formulation, with tension checking
 
- *Self-contact compatible*

- *Strategies, processes, solvers and convergence criterias used by the contact formulation*

- *Several strategies for adaptive remeshing*
 
- *The application includes search utilities in order to create the contact conditions*

- *Frictional laws (**WIP**) in order to consider different types of frictional behaviour*

- *+115 Python unittest, including Validation tests, and +85 cpp tests*

## ⚙️ Examples:

Examples can be found [here](https://github.com/KratosMultiphysics/Examples/tree/master/contact_structural_mechanics), and [here](https://github.com/KratosMultiphysics/Examples/tree/master/mmg_remeshing_examples/) for several contact adaptive remeshing examples.

## 🗎 Documentation:

Further information regarding the formulation can be accessed in Chapter 4 of the *PhD thesis* authored by [Vicente Mataix Ferrándiz](mailto:vmataix@altair.com), available on [UPC Commons](https://upcommons.upc.edu/bitstream/2117/328952/1/TVMF1de1.pdf).
