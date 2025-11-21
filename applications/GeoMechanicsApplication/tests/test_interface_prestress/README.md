# Interface Prestress Test

This test validates the functionality of applying prestress to interface elements based on neighbouring continuum soil elements.

For details on how neighbouring elements of interfaces are found using the FindNeighboursOfInterfacesProcess, see the [process documentation](../../custom_processes/README.md#find-neighbours-of-interfaces).

## Setup
The test consists of two stages. In the first stage, two U-Pw differential order continuum elements are created, connected by master-slave constraints at their coincident nodes. A uniform load of 1 kPa is applied on the top surface of the model. In the second stage, a line interface element is activated between the two continuum elements (and the master-slave constraints are deactivated). Due to the prestress applied to the interface element from the first stage, no additional displacements should occur in the second stage. A schematic representation of the setup is shown below:

```text
=== Stage 1 ===

+---+---+
|   |   |   (Uniform load load)
v   v   v
o---o---o
|       |
o       o   (Second U-Pw diff order element)
|       |
o---o---o   (These three nodes are coincident with...
o---o---o      ...these three nodes, and connected using master-slave constraints)
|       |
o       o   (First U-Pw diff order element)
|       |
o---o---o


=== Stage 2 ===

+---+---+
|   |   |   (Uniform load load)
v   v   v
o---o---o
|       |
o       o   (Second U-Pw diff order element)
|       |
o---o---o
:   :   :   (Line interface element becomes active; it replaces the master-slave constraints)
o---o---o
|       |
o       o   (First U-Pw diff order element)
|       |
o---o---o
```

Next to the aforementioned elements and top load, the following constraints are applied:
- The bottom nodes of the continuum elements are fixed in all directions.
- The side nodes are fixed in the horizontal direction.
- The water pressure is kept at zero in the entire domain.

A linear elastic material is used for both the continuum and the interface elements.

## Assertions

The following assertions are done to validate the prestress functionality:
- In the second stage, all stage displacements of the nodes are checked to be zero (within a small tolerance). This would not be the case if the interfaces were not prestressed based on the first stage stresses in the continuum elements.
- The stresses at the integration points of the interface elements are checked to be equal to the expected stresses at the integration points of the neighbouring continuum element. These are expected to be 1 kPa in the normal direction (following from the vertical equilibrium) and zero in the shear direction.
- The relative (normal and shear) displacements of the interface elements are checked to be zero (within a small tolerance).