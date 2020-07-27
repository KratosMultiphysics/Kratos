## Particle Mechanics Application Frequently Asked Questions (FAQs)

The aim of this page is to pull together some common tips and tricks for managing issues that arrise when using the Particle Mechanics application (mainly MPM). 

This page isn't an exhaustive list of fixes. If you encounter a problem that you don't find on here please:
1. Check the [issues page](https://github.com/KratosMultiphysics/Kratos/issues)
2. Raise an [issue](https://github.com/KratosMultiphysics/Kratos/issues)
3. Contact the team members at the bottom of the [Particle Mechanics homepage](https://github.com/KratosMultiphysics/Kratos/tree/MPM/linear_implicit/applications/ParticleMechanicsApplication)


## Quick fixes
1. `RuntimeError: Error: DEPRECATION: The ModelPart "XXX" is retrieved from the Model by using the flat-map! This was removed end of November 2019 Please prepend the Parent-ModelPart-Names like this: "Initial_MPM_Material.XXX"`~

Put `Initial_MPM_Material.` in front of the model_part_name in `ParticleMaterials.json`. You may also have to put `Background_Grid` in front of the model_part_name if you are apply some boundary condition process.



2. Error using `LinearElasticPlaneStress2DLaw`.

Use `LinearElasticIsotropicPlaneStress2DLaw` instead. This may occur with other material laws too.



3. Error using `mpm_gid_output_process`.

Use `particle_gid_output_process` instead as per issue [#5641](https://github.com/KratosMultiphysics/Kratos/issues/5641)

## MPM modelling in GiD
MPM modelling in GiD is relatively straightforward if you keep the following pointers in mind:
- Use the layers and groups to control the visibility of the background grid and the object of analysis. These always lay on top of each other and selecting items can be painful without hiding either the grid or object.
- Select background grid nodes for restraint groups after you have meshed. If you mesh again you might have to re-assign nodes to your restraint group.
