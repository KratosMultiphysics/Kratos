---
title: KratosAPI MultiPoint Constraint (MPC)
keywords: 
tags: [KratosAPI MultiPoint Constraint MPC]
sidebar: kratos_for_developers
summary: 
---

## Overview
MultiFreedom Constraints (**MFC**) are functional equations that connect two or more degrees of freedom components:

```
F(nodal degree of freedom)  = prescribed value
```

An **MFC** of this form is called multipoint or multinode (**MPC**) if it involves *DoF* components at different nodes. The constraint is called *linear* if all displacement components appear linearly on the left-hand-side, and *nonlinear* otherwise.

There are basically three methods for treating **MFC**. 

1. **Master-Slave Elimination**. The degrees of freedom involved in each MFC are separated into master and slave freedoms. The slave freedoms are then explicitly eliminated. The modified equations do not contain the slave freedoms.
2. **Penalty Augmentation**. Also called the penalty function method. Each MFC is viewed as
the presence of a fictitious elastic structural element called penalty element that enforces it approximately. This element is parametrized by a numerical weight. The exact constraint is recovered if the weight goes to infinity. The **MFC**s are imposed by augmenting the finite element model with the penalty elements.
3. **Lagrange Multiplier Adjunction**. For each MFC an additional unknown is adjoined to the master stiffness equations. Physically this set of unknowns represent constraint forces that would enforce the constraints exactly should they be applied to the unconstrained system.

These methods can be compared, the following table resumes the differences between them.

![Comparison](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/mfc_comparison.png)

The **MFC** implemented in *Kratos* is the *Master-Slave Elimination*, this will be the one presented in the [*Formulation*](#formulation) section.

## Formulation

### Problem introduction

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/bar_sctruture.png)

We will use a very simple structure example in order to explain how the formulation works. This structure consists of six bar elements connected by seven nodes that can only displace in the x direction. Before imposing various multifreedom constraints discussed below, the master stiffness equations for this problem are assumed to be

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/kuf.png)

Now let us specify a multifreedom constraint that states that nodes 2 and 6 are to move by the same
amount:

```
u2 = u6
```

Passing all node displacements to the right hand side gives the canonical form:

```
u2-u6=0
```

Constraintconditionsofthistypearesometimescalled rigid links becausetheycanbemechanically
interpreted as forcing node points 2 and 6 to move together as if they were tied by a rigid member.

### The Master-Slave method

To apply this method by hand, the *MFC*s are taken one at a time. For each constraint a slave degree of freedom is chosen. The freedoms remaining in that constraint are labeled master. A new set of degrees of freedom **u** is established by removing all slave freedoms from *u*. This new vector contains master freedoms as well as those that do not appear in the **MFC**s. A matrix transformation equation that relates *u* to **u** is generated. This equation is used to apply a congruential transformation to the master stiffness equations. This procedure yields a set of modified stiffness equations that are expressed in terms of the new freedom set **u**. Because the modified system does not contain the slave freedoms, these have been effectively eliminated.

The mechanics of the process is best seen by going through an example.

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/eq1.png)

This is the required transformation relation. In compact form:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/eq2.png)

Replacing and premultiplying by **T** yields the modified system:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/eq3.png)

Carrying out the indicated matrix multiplications yields:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/eq4.png)

## Kratos interface (how to use it)

The **MPC** are implemented in *Kratos* via the [*MasterSlaveConstraint* class](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/includes/master_slave_constraint.h), and particularly the linear kinematic relationship can be found in in the class [*LinearMasterSlaveConstraint* class](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/includes/linear_master_slave_constraint.h).

The *LinearMasterSlaveConstraint* is defined by the following:

* **Index**: The Id of the new created constraint
* **Vector of master DoF**: The vector containing the *DoF* of the master side
* **Vector of slave DoF**: The vector containing the *DoF* of the slave side
* **A relation matrix** : The relation matrix between the master/slave *DoF*. This is **T** in the theoretical formulation.
* **A constant vector**: The vector containing the additional kinematic relationship

Additionally to this, in order to consider the **MPC** into ths system assembly the *ResidualBasedBlockBuilderAndSolverWithConstraints* must be used as builder and solver.

The *MasterSlaveConstraint* are part of the *Mesh* class, so they are stored in the [*ModelPart*](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/includes/model_part.h) as the elements and conditions do. Then the *ModelPart* has an interface to access and create them. Particularly the methods `CreateNewMasterSlaveConstraint` are of interest, there are three options:

```cpp
MasterSlaveConstraint::Pointer CreateNewMasterSlaveConstraint(
	const std::string& ConstraintName,
	IndexType Id, 
	DofsVectorType& rMasterDofsVector,
	DofsVectorType& rSlaveDofsVector,
	const MatrixType& RelationMatrix,
	const VectorType& ConstantVector,
	IndexType ThisIndex = 0);

MasterSlaveConstraint::Pointer CreateNewMasterSlaveConstraint(
	const std::string& ConstraintName, 
	IndexType Id, 
	NodeType& rMasterNode,
	const DoubleVariableType& rMasterVariable,
	NodeType& rSlaveNode,
	const DoubleVariableType& rSlaveVariable,
	const double Weight,
	const double Constant,
	IndexType ThisIndex = 0);

MasterSlaveConstraint::Pointer CreateNewMasterSlaveConstraint(
	const std::string& ConstraintName, 
	IndexType Id, 
	NodeType& rMasterNode,
	const VariableComponentType& rMasterVariable,
	NodeType& rSlaveNode,
	const VariableComponentType& rSlaveVariable,
	double Weight,
	double Constant,
	IndexType ThisIndex = 0);
```

This can be seen in a simple example:

```python
import KratosMultiphysics

# Define a Model
current_model = KratosMultiphysics.Model()
mp = current_model.CreateModelPart("Main")

mp.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
mp.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

mp.CreateNewNode(1, 0.00000, 0.00000, 0.00000)
mp.CreateNewNode(2, 0.00000, 1.00000, 0.00000)

KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, mp)
KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, mp)
KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, mp)

mp.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, mp.Nodes[1], KratosMultiphysics.DISPLACEMENT_X, mp.Nodes[2], KratosMultiphysics.DISPLACEMENT_X, 1.0, 0)

# We impose the BC
mp.Nodes[1].Fix(KratosMultiphysics.DISPLACEMENT_X)
mp.Nodes[1].SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0.1)

#define a minimal newton raphson solver
linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()
builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
convergence_criterion = KratosMultiphysics.ResidualCriteria(1e-10, 1e-12)
convergence_criterion.SetEchoLevel(0)

max_iters = 100
compute_reactions = False
reform_step_dofs = True
move_mesh_flag = False
strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
    mp, scheme, linear_solver, convergence_criterion,
    builder_and_solver, max_iters, compute_reactions,
    reform_step_dofs, move_mesh_flag)
strategy.SetEchoLevel(0)
strategy.Initialize()

strategy.Check()
strategy.Solve()

print(mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X))

assert mp.Nodes[2].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X) == mp.Nodes[1].GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X)
```

**IMPORTANT: The slaves nodes must have the flag SLAVE set to true**

## How to create en example
In this section a brief tutorial of "How To" create an example is explained.

### Set up of the interface
In order to generate the required input data files we use the geometry pre and post processor GiD (https://www.gidhome.com/download/official-versions). Additionally we need to download the KratosMultiphysics GiD interface (https://github.com/KratosMultiphysics/GiDInterface - > download .zip).
Once we open the GiD software we must load the Kratos problem type by:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto1.png)

And then select the StructuralMechanicsApplication as:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto2.png)

And then we select 3D for this particular case:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto3.png)

Then a new menu will appear on the left with all the available options. In order to create MPC contact problems you should switch the Kratos Problemtype to _Developer Mode_ in `Kratos->Kratos Preferences`.

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto4.png)

### Creating a problem
Now we have generated an arbitrary geometry like the one below.

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto5.png)

It is **Very important** to check the initial orientation of the normal of the surfaces (lines in 2D). This can be done by `Utilities->Swap Normals->Surfaces->Select`. One has to ensure that **all the normals point outwards** the bodies like in the folowing picture.

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto6.PNG)

Next we assign the material properties to each of the bodies (we recommend to assign separately the mat props to each body in order to see them separately in the post-process) as:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto7.PNG)

Then we must apply the Boundary conditions (in thiscase we fix the tips of the beam) as:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto8.PNG)

And the **MASTER-SLAVE CONDITIONS**. The assignation of this Master-Slave conditions is mandatory to detect the contact between the potential contact surfaces of the bodies. We should assign one surface of the body 1 as MASTER and the corresponding contacting body (or body 2) the slave flag.
In our case, we have assigned as master the contact surfaces o the grey body and as slave the ones of the blue body.

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto9.PNG)

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/tuto10.PNG)

Afterwards we could apply several loading conditions with the _Loads_ folder. Once we have finished assigning loads, we generate the mesh (**we strongly recommend to generate Structured meshes**) and create the input data files for kratos by `Kratos->Write Calculation Files`. Several files will be generated in your problem folder:

* ProjectParameters.json (**to be slightly modified for the MPC**)
* problem_name.mdpa
* StructuralMaterials.json
* MainKratos.py

### Modifications to the ProjectParameters.json

You should change:

```json
"convergence_criterion": "contact_residual_criterion"
```

by 

```json
"convergence_criterion": "residual_criterion"
```

You should change:

```json
"contact_settings": {
    "mortar_type": "ALMContactFrictionlessComponents"
}
```

by 

```json
"mpc_contact_settings": {
    "contact_type": "MeshTying"
}
```

and replace 

```json
{
    "python_module" : "alm_contact_process",
    "kratos_module" : "KratosMultiphysics.ContactStructuralMechanicsApplication",
    "process_name"  : "ALMContactProcess",
    "Parameters"    : {
        "model_part_name"     : "Structure",
        "contact_model_part"  : {
            "0" : ["CONTACT_Contact_slave_Auto1","CONTACT_Contact_master_Auto1"]
        },
        "assume_master_slave" : {
            "0" : ["CONTACT_Contact_slave_Auto1"]
        },
        "contact_type"        : "FrictionlessComponents"
    }
}
```

by 

```json
{
    "python_module" : "mpc_contact_process",
    "kratos_module" : "ContactStructuralMechanicsApplication",
    "process_name"  : "MPCContactProcess",
    "Parameters"    : {
        "model_part_name"     : "Structure",
        "contact_model_part"  : {
            "0" : ["CONTACT_Contact_slave_Auto1","CONTACT_Contact_master_Auto1"]
        },
        "assume_master_slave" : {
            "0" : ["CONTACT_Contact_slave_Auto1"]
        },
        "contact_type"        : "MeshTying",
        "reaction_check_stiffness_factor" : 2.0e-3
    }
}
```

### Final result
In our case we applied a body force over the whole volume, obtaining the following deformed shape:

![](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Wiki_files/MPC_wiki/post.png)

## References
[1] "MultiFreedom Constraint I " by C. Felippa. [Link](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Resources_files/Colorado%20IFEM%20course/IFEM.Ch08.pdf)

[2] "MultiFreedom Constraint II " by C. Felippa. [Link](https://raw.githubusercontent.com/KratosMultiphysics/Documentation/master/Resources_files/Colorado%20IFEM%20course/IFEM.Ch09.pdf)

[3] "Multi-Point Constraints". [Link](https://mashayekhi.iut.ac.ir/sites/mashayekhi.iut.ac.ir/files//files_course/lesson_16.pdf)