## Contact Structural Mechanics Application 
 
The Contact Structural Mechanics Application contains the contact mechanics implementations that can be used by the Structural Mechanics Application within Kratos Multiphysics. 
 
<p align="center">
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/validation/double_arch/data/result.gif" alt="Solution" style="width: 600px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/validation/double_arch/data/result_frictional.gif" alt="Solution" style="width: 600px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/use_cases/in_ring/data/animation.gif" alt="Solution" style="width: 600px;"/>
  <img src="https://github.com/KratosMultiphysics/Examples/raw/master/contact_structural_mechanics/use_cases/hyperelastic_tubes/data/half_cylinders.gif" alt="Solution" style="width: 600px;"/>
 <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/mmg_remeshing_examples/use_cases/contacting_cylinders/data/nodal_h.gif" alt="Solution" style="width: 600px;"/>
  <img src="https://raw.githubusercontent.com/KratosMultiphysics/Examples/master/contact_structural_mechanics/use_cases/self_contact/data/animation.gif" alt="Solution" style="width: 600px;"/>
</p>


The application includes tests to check the proper functioning of the application
 
### Features: 
 
- Mesh tying conditions based in mortar formulation
 
- Augmented Lagrangian contact conditions based in mortar formulation
 
    * Frictionless formulation

    * Frictional formulation

- Penalty contact conditions based in mortar formulation

     * Frictionless formulation

     * Frictional

- Simplified **MPC** conditions based in mortar formulation. With the mortar formulation the weight are computed, allowing to compute a Simplified *NTN* and a simplified *NTS*

     * Frictionless formulation

     * Frictional formulation

     * Mesh tying formulation, with tension checking
 
- Self-contact compatible

- Strategies, processes, solvers and convergence criterias used by the contact formulation

- Several strategies for adaptive remeshing
 
- The application includes search utilities in order to create the contact conditions

- Frictional laws (WIP) in order to consider different types of frictional behaviour 

- +115 *Python* unittest, including Validation tests, and +85 cpp tests

### Examples:

Examples can be found [here](https://github.com/KratosMultiphysics/Examples/tree/master/contact_structural_mechanics), and [here](https://github.com/KratosMultiphysics/Examples/tree/master/mmg_remeshing_examples/) for several contact adaptive remeshing examples
