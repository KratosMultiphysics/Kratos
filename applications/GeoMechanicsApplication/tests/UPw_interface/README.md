# UPw_interface

This set contains compact 2D column benchmarks for validating line-interface behavior in the GeoMechanics `U_Pw` formulation.

The tests compare:
- a continuous soil column (no explicit interface), 
- two stacked soil blocks connected by an explicit horizontal line interface, and
- a continuous soil column that has a fixture to the left boundary via an interface.

Both same-order and diff-order element families are included. The figure below shows the column and constraints. The dashed line shows the position of the horizontal interface case.

![Setup](setup.svg)

## Cases

### Same-order element cases

- `column`
	- Soil elements: `UPwSmallStrainElement2D4N` (2 quadrilaterals)
	- No interface element

- `column_horizontal_interface`
	- Soil elements: `UPwSmallStrainElement2D4N` (2 quadrilaterals with duplicated nodes at the internal boundary)
	- Interface element: `Geo_ULineInterfacePlaneStrainElement2Plus2N` (a single line element based on the duplicated nodes)
	
- `column_vertical_interface`
	- Soil elements: `UPwSmallStrainElement2D4N` (2 quadrilaterals)
	- Interface elements: `Geo_ULineInterfacePlaneStrainElement2Plus2N` (2 line elements that connect left side nodes of the soil elements with the boundary nodes)

### Different-order element cases

- `column_diff_order_elements`
	- Soil elements: `SmallStrainUPwDiffOrderElement2D6N` (4 triangles with midside nodes)
	- No interface element

- `column_horizontal_interface_diff_order_elements`
	- Soil elements: `SmallStrainUPwDiffOrderElement2D6N` (4 triangles with duplicated nodes at the internal boundary)
	- Interface element: `Geo_UPwLineInterfacePlaneStrainDiffOrderElement3Plus3N` (a single line element)

## Geometry and boundary conditions

- Column size: $1 \times 2\ \mathrm{[m]}$
- Bottom: fixed in both directions.
- Left and right sides: fixed in the horizontal direction.
- Stage 1: gravity loading + K0 initialization
- Stage 2:
	- prescribed top vertical displacement (table to $-0.01\ \mathrm{[m]}$)
	- prescribed top water pressure table (which varies linearly in time from $0$ (at $`t = 0`$) to $1000\ \mathrm{[Pa]}$ (at $`t = 1`$))

## Material models

### Soil (`PorousDomain.Soil`)

- Constitutive law: `GeoLinearElasticPlaneStrain2DLaw`
- $E = 1.0e9\ \mathrm{[Pa]}$, 
- $\nu = 0.2$
- $\rho_s = 2000\ \mathrm{[kg/m^3]}$, $\rho_w = 1000\ \mathrm{[kg/m^3]}$
- $n = 0.3$, $K_s = 1.0e12\ \mathrm{[Pa]}$, $K_w = 2.2e9\ \mathrm{[Pa]}$
- $k_{xx}=k_{yy}=4.5e-10\ \mathrm{[m^2]}$, 
- $\mu=1.0e-3\ \mathrm{[Pa\cdot s]}$

### Interface (`PorousDomain.Interface`)

- Constitutive law: `GeoIncrementalLinearElasticInterfaceLaw`
- $k_{\mathrm{n}} = 6.4e10\ \mathrm{[N/m^3]}$
- $k_{\mathrm{t}} = 2.6000869565e10\ \mathrm{[N/m^3]}$
- $k_{\perp} = 5.0e-4\ \mathrm{[m^2]}$
- $\mu=1.0e-3\ \mathrm{[Pa\cdot s]}$

## Assertions

Displacement and water pressure values obtained with and without the interface are compared for all nodes. The following pairs are used:

- `column_horizontal_interface` and `column`,
- `column_vertical_interface` and `column`,
- `column_horizontal_interface_diff_order_elements` and `column_diff_order_elements`.
