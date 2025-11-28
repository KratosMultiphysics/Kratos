# Submerged construction of an excavation

This document describes excavation modeling test.

## Geometry
The picture shows the modelling set-up.
![Set up](setup.svg)

The geometry is 2D, with a vertical symmetry line at the center of the excavated area, allowing only half of the domain to be modeled. The excavated area is shown as a shaded region.

- Excavation dimensions: total width $`30~m`$, final depth $`20~m`$.
- Boundary conditions: vertical diaphragm walls extend $`30~m`$ into the ground, with struts near their tops. Walls and struts are depicted in red and blue, respectively.
- Constraints: the domain bottom is fully fixed when domain sides are fixed only in the horizontal direction.
- External load: distributed load of $`5~Kn/m^2/m`$ near the wall (orange in the figure).

The ground consists of two homogeneous layers:
1. Top 20 m: clay
2. Bottom: sand 

Both layers have homogenous properties, which are listed in [Material properties section](#material-properties).

The analysis is performed in seven stages, described in the [Staged analysis section](#staged-analysis).

## Mesh

The mesh is generated in GiD and shown below:

![Mesh](mesh.png)

It is refined near the walls and the interface between the two soil layers.

### Mesh preparation workflow

1. Assign unique element IDs for interface elements
- Prepend left-side interface element IDs with 1 and right-side IDs with 2.
- Update sub-model part element IDs accordingly.

2. Convert line elements to interface elements

   The line elements that represent the soil-structure interfaces are converted to actual line interface elements with 
```shell
# Run the following commands in the directory containing the `.mdpa` file
python insert_interfaces_in_excavation_model.py  
```
This produces a new `.mdpa` file (i.e., it doesn't overwrite the existing one): `submerged_excavation_gid_project_with_interfaces.mdpa`.  **Note**, that node IDs in sub-model parts that were **not** given to the interface inserter, won't be updated.  This typically occurs for "overarching" sub-model parts that put other sub-model parts together.

3. Create master-slave constraints using the FEA-Tools script:

```shell
python <path-to-FEA-Tools-repo>\MiscellaneousTools\InterfaceInserter\insert_master_slave_constraints.py Clay_Left,Sand_Left Clay_Upper_Right,Clay_Middle_Right,Clay_Lower_Right,Sand_Right submerged_excavation_gid_project_with_interfaces.mdpa 1_Initial_stage.json --search_radius 0.01
```
This script produces two new files: 
- `out_submerged_excavation_gid_project_with_interfaces.mdpa` and 
- `out_1_Initial_stage.json`.  

After inspecting the differences between the original files and the new files, the new files can be renamed to the original file names.

4. Adjust truss elements

GiD generates 3-noded trusses; convert them to 2-noded by removing the final node and setting the type to LinearTrussElement2D2N.

5. Update line load conditions

For the line load condition (that is being used for the surface load),
- change its type from `LineLoadCondition2D3N` to `LineLoadDiffOrderCondition2D3N`.
- Swap the end nodes of these conditions  to ensure proper orientation. Prior to swapping, the process that finds neighboring elements would not work, because of the reversed orientations.

6. Adjust water pressure conditions

- Replace the corner node of the clay bottom next to the diaphragm wall with the duplicated node.


## Material properties

Material properties for soil, interfaces, diaphragm walls, and struts are summarized below.

For the initial stage (which uses a $`K_0`$ procedure), linear elastic material properties are defined by a Young's modulus $`E`$ and a Poisson's ratio $`\nu`$. 

The following table lists the material properties of the soil layers that have been adopted by the Kratos model.

| Property                                                  | Kratos input parameter | Clay                             | Sand                             | Unit                           |
|-----------------------------------------------------------|------------------------|----------------------------------|----------------------------------|--------------------------------|
| Grain density $`\rho_{\mathrm{g}}`$                       | `DENSITY_SOLID`        | 2038.74                          | 2475.61                          | $`\mathrm{kg} / \mathrm{m}^3`$ |
| Water density $`\rho_{\mathrm{w}}`$                       | `DENSITY_WATER`        | 1019.37                          | 1019.37                          | $`\mathrm{kg} / \mathrm{m}^3`$ |
| Porosity $`n`$                                            | `POROSITY`             | 0.20                             | 0.30                             | $`[-]`$                        |
| Retention law type                                        | `RETENTION_LAW`        | `SaturatedBelowPhreaticLevelLaw` | `SaturatedBelowPhreaticLevelLaw` | N/A                            |
| Residual saturation $`S_{\mathrm{res}}`$                  | `RESIDUAL_SATURATION`  | $`1 \cdot 10^{-10}`$             | $`1 \cdot 10^{-10}`$             | $`[-]`$                        |
| Saturated saturation $`S_{\mathrm{sat}}`$                 | `SATURATED_SATURATION` | 1.0                              | 1.0                              | $`[-]`$                        |
| Young's modulus $`E`$                                     | `YOUNG_MODULUS`        | $`12 \cdot 10^3`$                | $`120 \cdot 10^3`$               | $`\mathrm{kN} / \mathrm{m}^2`$ |
| Poisson's ratio $`\nu`$                                   | `POISSON_RATIO`        | 0.15                             | 0.20                             | $`[-]`$                        |
| Cohesion $`c`$                                            | `GEO_COHESION`         | 1000.0                           | 0.0                              | $`\mathrm{N} / \mathrm{m}^2`$  |
| Friction angle $`\phi`$                                   | `GEO_FRICTION_ANGLE`   | 25.0                             | 32.0                             | $`^{\circ}`$                   |
| Dilatancy angle $`\psi`$                                  | `GEO_DILATANCY_ANGLE`  | 0.0                              | 2.0                              | $`^{\circ}`$                   |
| K0-value for normal consolidation $`K_{\mathrm{0}}^{nc}`$ | `K0_NC`                | 0.5774                           | 0.4701                           | $`[-]`$                        |


The following table lists the material properties of the interfaces that have been adopted by the Kratos model.

| Property                                 | Kratos input parameter | Clay      | Sand    | Unit                           |
|------------------------------------------|------------------------|-----------|---------|--------------------------------|
| Normal stiffness $`k_{\mathrm{n}}`$      | `NORMAL_STIFFNESS`     | 48000     | 480000  | $`\mathrm{kN} / \mathrm{m}^3`$ |
| Shear stiffness $`k_{\mathrm{s}}`$       | `SHEAR_STIFFNESS`      | 20869.565 | 200000  | $`\mathrm{kN} / \mathrm{m}^3`$ |
| Cohesion $`c`$                           | `GEO_COHESION`         | 1000.0    | 0.0     | $`\mathrm{N} / \mathrm{m}^2`$  |
| Friction angle $`\phi`$                  | `GEO_FRICTION_ANGLE`   | 25.0      | 32.0    | $`^{\circ}`$                   |
| Dilatancy angle $`\psi`$                 | `GEO_DILATANCY_ANGLE`  | 0.0       | 2.0     | $`^{\circ}`$                   |
| Tensile strength $`\sigma_{\mathrm{t}}`$ | `GEO_TENSILE_STRENGTH` | 2.1445    | 0.0     | $`\mathrm{N} / \mathrm{m}^2`$  |

The following table lists the material properties of the diaphragm wall that have been adopted by the Kratos model.

| Property                                       | Kratos input parameter  | Diaphragm wall      | Unit                           |
|------------------------------------------------|-------------------------|---------------------|--------------------------------|
| Young's modulus $`E`$                          | `YOUNG_MODULUS`         | $`5.93 \cdot 10^6`$ | $`\mathrm{kN} / \mathrm{m}^2`$ |
| Poisson's ratio $`\nu`$                        | `POISSON_RATIO`         | 0.0                 | $`[-]`$                        |
| Thickness $`t`$                                | `THICKNESS`             | 1.265               | $`[m]`$                        |
| Effective shear thickness y $`t_{\mathrm{y}}`$ | `THICKNESS_EFFECTIVE_Y` | 0.025               | $`[m]`$                        |
| Density $`\rho`$                               | `DENSITY`               | 1019.368            | $`\mathrm{kg} / \mathrm{m}^3`$ |

The following table lists the material properties of the strut that have been adopted by the Kratos model.

| Property                  | Kratos input parameter | Strut               | Unit                           |
|---------------------------|------------------------|---------------------|--------------------------------|
| Young's modulus $`E`$     | `YOUNG_MODULUS`        | $`2.1 \cdot 10^8`$  | $`\mathrm{kN} / \mathrm{m}^2`$ |
| Density $`\rho`$          | `DENSITY`              |  7850.0             | $`\mathrm{kg} / \mathrm{m}^3`$ |
| Cross-sectional area $`A` | `CROSS_SECTIONAL_AREA` | $`1.9 \cdot 10^-3`$ | $`[m^2]`$                      |
| Truss pre-stress PK2      | `TRUSS_PRESTRESS_PK2`  | 0.0                 | $`\mathrm{kN} / \mathrm{m}^3`$ |


## Staged analysis

The excavation is modeled in seven stages:

1. **Initial stage:**

- Entire soil domain is active, the diaphragm wall as well the interfaces at both sides of it are inactive.
- Master-slave constraints are applied where the soil will later be separated by the diaphragm wall. This ensures continuity of the displacement field in this early stage of analysis.  
- The only load that is being applied is self-weight.  
- At the end of the stage, a $`K_0`$ procedure is performed to initialize the horizontal stress field.  **Note**, the $`K_0`$ procedure requires the use of linear elastic materials for all soil parts.

2. **Null step:**

- The material models change from linear elastic to Mohr-Coulomb model. This triggers a stiffness redistribution, and hence a stress redistribution.

3. **Diaphragm wall installation and applying the external load:**

- Activation of the diaphragm wall and the interfaces that are attached to its left and right sides.
- Deactivate the master-slave constraints, since there is no need in the continuous displacement field across the entire domain.
- Apply a surface load to a part of the top of the soil on the left hand side.

4. **First excavation stage:**

- Excavate the top part of the clay to the right of the diaphragm wall (deactivate corresponding model part).
- Deactivate the interface elements that connect the now excavated clay layer to the diaphragm wall.

5. **Installation of a strut:**

- Activate the strut.

6. **Second excavation stage:**

- Excavate the next top portion of clay to the right of the wall.
- Deactivate corresponding model parts and interface elements.
- Apply a normal contact stress to the exposed part of the diaphragm wall as well as the bottom of the excavation pit.

7. **Third excavation stage:**

- Excavate the final top portion of clay to the right of the wall.
- Deactivate corresponding model parts and interface elements.
- Apply a normal contact stress to the exposed part of the diaphragm wall as well as the bottom of the excavation pit.