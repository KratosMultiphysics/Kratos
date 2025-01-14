# Thermal DEM Application

This application is an extension of the [DEM Application](https://github.com/KratosMultiphysics/Kratos/tree/master/applications/DEMApplication) to include **thermal effects** such as:

- Heat transfer between particle-particle, particle-wall, and particle-surrounding fluid.
- Heat transfer mechanisms by conduction, convection, and radiation.
- Heat generation by energy dissipation and internal sources.
- Temperature dependent material properties.

Theoretical information on thermal DEM analysis can be found [here](https://doi.org/10.5281/zenodo.13846126).

A [Matlab version](https://github.com/rlrangel/DEMLab) of this application is also available.

## Table of Contents
- [Authorship](#authorship)
- [License](#license)
- [Getting Started](#getting-started)
- [Instructions](#instructions)
    - [MainKratos](#mainkratos-python-file)
    - [Project Parameters](#project-parameters-json-file)
    - [Materials](#materials-json-file)
    - [DEM Model Parts](#dem-model-parts-mdpa-file)
    - [DEM-FEM Boundary Model Parts](#dem-fem-boundary-model-parts-mdpa-file)
- [Input Explanations](#input-explanations)
- [Testing](#testing)

## Authorship

- **Rafael Rangel** - (<rrangel@cimne.upc.edu>)

International Center for Numerical Methods in Engineering ([CIMNE](https://www.cimne.com/))

## License

The Thermal DEM application is OPEN SOURCE.
The main code and program structure is available and aimed to grow with the need of any users willing to expand it.
The BSD (Berkeley Software Distribution) licence allows to use and distribute the existing code without any restriction,
but with the possibility to develop new parts of the code on an open or close basis depending on the developers.

## Getting Started

This application is part of the ***Kratos Multiphysics*** framework.
Instructions on how to get a copy of the project and run on your local machine for development and testing purposes are available for both
[Linux](http://kratos-wiki.cimne.upc.edu/index.php/LinuxInstall) and [Windows](http://kratos-wiki.cimne.upc.edu/index.php/Windows_7_Download_and_Installation) systems.

Before building *Kratos Multiphysics*, make sure to add the following applications to your configure file: 

- ThermalDEMApplication
- DelaunayMeshingApplication

## Instructions

To create a model for the Thermal DEM Application, the following adaptations must be done to the input files of the DEM Application.

### MainKratos (python file)

Replace the import of the *AnalysisStage* of the DEM Application with the *AnalysisStage* of the Thermal DEM Application, which is imported as follows:

	from KratosMultiphysics.ThermalDEMApplication.thermal_dem_analysis import ThermalDEMAnalysis

You can check the template of the *MainKratos* file [here](https://github.com/KratosMultiphysics/Kratos/blob/4c8a07592cf4056557be0e5dff1c19839a5e3b98/applications/ThermalDEMApplication/python_scripts/MainKratosThermalDEMAnalysis.py).

### Project Parameters (json file)

Change **solver strategy**.\
Available options: *thermal_sphere_strategy*.

	"solver_settings": {
		"strategy": "thermal_sphere_strategy"
	}

Change **element type**.\
Available options: *ThermalSphericPartDEMElement3D*.

	"ElementType": "ThermalSphericPartDEMElement3D"

Add **thermal settings** with desired options:

	"thermal_settings" : {
		"thermal_integration_scheme"     : "forward_euler",
		"numerical_integration_method"   : "adaptive_simpson",
		"thermal_solve_frequency"        : 1,
		"voronoi_tesselation_frequency"  : 1000,
		"porosity_update_frequency"      : 1000,
		"automatic_solve_frequency"      : true or false,
		"compute_motion"                 : true or false,
		"compute_direct_conduction"      : true or false,
		"compute_indirect_conduction"    : true or false,
		"compute_convection"             : true or false,
		"compute_radiation"              : true or false,
		"compute_heat_generation"        : true or false,
		"compute_adjusted_contact"       : true or false,
		"direct_conduction_model"        : "batchelor_obrien_simple" or "batchelor_obrien_complete" or "batchelor_obrien_modified" or "thermal_pipe" or "collisional",
		"indirect_conduction_model"      : "surrounding_layer" or "voronoi_a" or "voronoi_b" or "vargas_mccarthy",
		"nusselt_correlation"            : "sphere_hanz_marshall" or "sphere_whitaker" or "sphere_gunn" or "sphere_li_mason",
		"radiation_model"                : "continuum_zhou" or "continuum_krause",
		"heat_generation_model"          : ["sliding_friction","rolling_friction","contact_damping"],
		"adjusted_contact_model"         : "zhou" or "lu" or "morris",
		"voronoi_method"                 : "tesselation" or "porosity",
		"porosity_method"                : "global" or "average_convex_hull" or "average_alpha_shape",
		"min_conduction_distance"        : 0.0000000275,
		"max_conduction_distance"        : 1.0,
		"conduction_radius"              : 1.0,
		"fluid_layer_thickness"          : 0.4,
		"isothermal_core_radius"         : 0.5,
		"max_radiation_distance"         : 2.0,
		"heat_generation_ratio"          : 1.0,
		"global_porosity"                : 0.0,
		"alpha_shape_parameter"          : 1.2,
		"integral_tolerance"             : 0.000001,
		"heat_map_corners"               : [[0,0,0],[1,1,1]],
		"heat_map_subdivisions"          : [10,10,10],
		"global_fluid_properties"        : {
			"fluid_density"              : 1.0,
			"fluid_viscosity"            : 1.0,
			"fluid_thermal_conductivity" : 1.0,
			"fluid_heat_capacity"        : 1.0,
			"fluid_temperature"          : 0.0,
			"fluid_velocity_X"           : 0.0,
			"fluid_velocity_Y"           : 0.0,
			"fluid_velocity_Z"           : 0.0
		}
	}

Add **post options** with desired options:

	"PostTemperature"                      : true or false,
	"PostHeatFlux"                         : true or false,
	"PostGraphParticleTempMin"             : true or false,
	"PostGraphParticleTempMax"             : true or false,
	"PostGraphParticleTempAvg"             : true or false,
	"PostGraphParticleTempAvgVol"          : true or false,
	"PostGraphParticleTempDev"             : true or false,
	"PostGraphMechanicalEnergy"            : true or false,
	"PostGraphDissipatedEnergy"            : true or false,
	"PostGraphThermalEnergy"               : true or false,
	"PostGraphHeatFluxContributions"       : true or false,
	"PostGraphHeatGenerationValues"        : true or false,
	"PostGraphHeatGenerationContributions" : true or false,
	"PostHeatMapGeneration"                : true or false

### Materials (json file)

Add **material variables and tables**:

	"materials": [{
		"Variables": {
			"THERMAL_CONDUCTIVITY"          : 0.0,
			"SPECIFIC_HEAT"                 : 0.0,
			"EMISSIVITY"                    : 0.0,
			"THERMAL_EXPANSION_COEFFICIENT" : 0.0
		},
		"Tables": {
			"TABLE_NAME": {
				"input_variable"  : "TEMPERATURE",
				"output_variable" : "PARTICLE_DENSITY" or "YOUNG_MODULUS" or "POISSON_RATIO" or "SPECIFIC_HEAT" or "THERMAL_CONDUCTIVITY" or "EMISSIVITY" or "THERMAL_EXPANSION_COEFFICIENT",
				"data"            : [[0.0,0.0],[0.0,0.0]]
			}
		}
	}]

Add **material relations variables**:

	"material_relations": [{
		"Variables": {
			"DYNAMIC_FRICTION" : 0.0
		}
	}]

### DEM Model Parts (mdpa file)

Change **elements type**.\
Available options: *ThermalSphericParticle*.

	Begin Elements ThermalSphericParticle

Add **SubModelPartData** to sub model parts with desired options:

	Begin SubModelPart PART_NAME
		Begin SubModelPartData
			TEMPERATURE               0.0
			FIXED_TEMPERATURE         0
			HEATFLUX                  0.0
			HEATSOURCE                0.0
			ADIABATIC                 0
			REAL_YOUNG_MODULUS_RATIO  0.0
		End SubModelPartData
	End SubModelPart

### DEM-FEM Boundary Model Parts (mdpa file)

Add **NodalData TEMPERATURE** to all nodes:

	Begin NodalData TEMPERATURE
		1  0  0.0
		2  0  0.0
		3  0  0.0
		...
	End NodalData

Add **SubModelPartData** to sub model parts with desired options:

	Begin SubModelPart PART_NAME
		Begin SubModelPartData
			ADIABATIC   0
		End SubModelPartData
	End SubModelPart

## Input Explanations

**Thermal settings**
- *"thermal_integration_scheme"*:\
  Selected scheme for time integration of thermal problem.\
  Default: "forward_euler"

- *"numerical_integration_method"*:\
  Selected metohd for solving integration expressions numerically.\
  Default: "adaptive_simpson"

- *"thermal_solve_frequency"*:\
  Number of steps in which thermal problem is solved.\
  Default: 1

- *"voronoi_tesselation_frequency"*:\
  Number of steps in which Voronoi diagram is built, in case it is required.\
  Default: 1000

- *"porosity_update_frequency"*:\
  Number of steps in which porosity is computed, in case it is required.\
  Default: 1000

- *"automatic_solve_frequency"*:\
  Boolean for automatically setting the thermal solve frequency based on the maximum allowed time step (it overrides the value set for thermal_solve_frequency).\
  Default: false

- *"compute_motion"*:\
  Boolean for solving mechanical problem.\
  Default: true

- *"compute_direct_conduction"*:\
  Boolean for computing heat transfer between elements by direct conduction.\
  Default: true

- *"compute_indirect_conduction"*:\
  Boolean for computing heat transfer between elements by indirect conduction.\
  Default: false

- *"compute_convection"*:\
  Boolean for computing heat transfer from fluid to particles by convection.\
  Default: false

- *"compute_radiation"*:\
  Boolean for computing heat transfer between elements by radiation.\
  Default: false

- *"compute_heat_generation"*:\
  Boolean for computing heat generation by energy dissipation between elements.\
  Default: false

- *"compute_adjusted_contact"*:\
  Boolean for adjusting contact geometry between elements according to the real value of Young modulus.\
  Default: false

- *"direct_conduction_model"*:\
  Selected model for simulating heat transfer by direct conduction.\
  Default: "batchelor_obrien_simple"

- *"indirect_conduction_model"*:\
  Selected model for simulating heat transfer by indirect conduction.\
  Default: "surrounding_layer"

- *"nusselt_correlation"*:\
  Selected correlation for computing Nusselt number and thus simulating heat convection between solids and fluid.\
  Default: "sphere_hanz_marshall"

- *"radiation_model"*:\
  Selected model for simulating heat transfer by radiation.\
  Default: "continuum_zhou"

- *"heat_generation_model"*:\
  List of selected models for simulating heat generation by energy dissipation.\
  Default: ["sliding_friction"]

- *"adjusted_contact_model"*:\
  Selected model for adjusting contact geometry.\
  Default: "zhou"

- *"voronoi_method"*:\
  Selected method for obtaining data from Voronoi diagram, in case it is required.\
  Default: "tesselation"

- *"porosity_method"*:\
  Selected method for computing porosity, in case it is required.\
  Default: "average_alpha_shape"

- *"min_conduction_distance"*:\
  Minimum heat conduction distance required for indirect conduction model "surrounding_layer".\
  Default: 0.0000000275

- *"max_conduction_distance"*:\
  Maximum distance for heat conduction (ratio of particles radii) required for conduction model "batchelor_obrien_complete", "batchelor_obrien_modified", "voronoi_a" and "voronoi_b".\
  Default: 1.0

- *"conduction_radius"*:\
  Radius of cylindrical conductive region (ratio of particles radii) required for conduction model "batchelor_obrien_complete" and "batchelor_obrien_modified".\
  Default: 1.0

- *"fluid_layer_thickness"*:\
  Thickness of particle fluid layer (ratio of particles radii) required for conduction model "surrounding_layer".\
  Default: 0.4

- *"isothermal_core_radius"*:\
  Radius of particle isothermal core (ratio of particles radii) required for conduction model "voronoi_b".\
  Default: 0.5

- *"max_radiation_distance"*:\
  Maximum distance for heat radiation (ratio of particles radii) required for all radiation models.\
  Default: 2.0

- *"heat_generation_ratio"*:\
  Ratio of dissipated energy that is converted into heat.\
  Default: 1.0

- *"global_porosity"*:\
  Prescribed value for the porosity, required for porosity computation method "global".\
  Default: 0.0

- *"alpha_shape_parameter"*:\
  Alpha radius for the alpha-shape method, required for porosity computation method "average_alpha_shape".\
  Default: 1.2

- *"integral_tolerance"*:\
  Numerical tolerance for numerical integration, in case it is required.\
  Default: 0.000001

- *"heat_map_corners"*:\
  Coordinates of two opposite corners of a cuboid that define the region where heat map is evaluated.\
  Default: [[0,0,0],[1,1,1]]

- *"heat_map_subdivisions"*:\
  Number of subdivisions in X,Y,Z directions to defined the grid of the heat map.\
  Default: [10,10,10]

- *"global_fluid_properties"*:\
  Prescribed values for the properties of the interstitial fluid, assumed as constant throughout all the analysis (fluid behavior does not change as it is not simulated).

- *"global_fluid_properties.fluid_density"*:\
  Fixed value of fluid density.\
  Default: 1.0

- *"global_fluid_properties.fluid_viscosity"*:\
  Fixed value of fluid viscosity.\
  Default: 1.0

- *"global_fluid_properties.fluid_thermal_conductivity"*:\
  Fixed value of fluid thermal conductivity.\
  Default: 1.0

- *"global_fluid_properties.fluid_heat_capacity"*:\
  Fixed value of fluid heat capacity.\
  Default: 1.0

- *"global_fluid_properties.fluid_temperature"*:\
  Fixed value of fluid temperature.\
  Default: 0.0

- *"global_fluid_properties.fluid_velocity_X"*:\
  Fixed value of fluid velocity in global X direction.\
  Default: 0.0

- *"global_fluid_properties.fluid_velocity_Y"*:\
  Fixed value of fluid velocity in global Y direction.\
  Default: 0.0

- *"global_fluid_properties.fluid_velocity_Z"*:\
  Fixed value of fluid velocity in global Z direction.\
  Default: 0.0

**Post options**
- *"PostTemperature"*:\
  Boolean for showing elements temperature in post processing.\
  Default: false

- *"PostHeatFlux"*:\
  Boolean for showing elements heat flux in post processing.\
  Default: false

- *"PostGraphParticleTempMin"*:\
  Boolean for writing a graph with the minimum particle temperature.\
  Default: false

- *"PostGraphParticleTempMax"*:\
  Boolean for writing a graph with the maximum particle temperature.\
  Default: false

- *"PostGraphParticleTempAvg"*:\
  Boolean for writing a graph with the average temperature of all particles.\
  Default: false

- *"PostGraphParticleTempAvgVol"*:\
  Boolean for writing a graph with the volume-weighted average temperature of all particles.\
  Default: false

- *"PostGraphParticleTempDev"*:\
  Boolean for writing a graph with the standard deviation of the temperature of all particles.\
  Default: false

- *"PostGraphMechanicalEnergy"*:\
  Boolean for writing a graph with the mechanical energy components of all partilces.\
  Default: false

- *"PostGraphDissipatedEnergy"*:\
  Boolean for writing a graph with the accumulated energy dissipation of all partilces.\
  Default: false

- *"PostGraphThermalEnergy"*:\
  Boolean for writing a graph with the accumulated thermal energy (from heat generation by energy dissipation) of all partilces.\
  Default: false

- *"PostGraphHeatFluxContributions"*:\
  Boolean for writing a graph with the contribution of each heat transfer mechanism to the total heat transfer.\
  Default: false

- *"PostGraphHeatGenerationValues"*:\
  Boolean for writing a graph with the total values of  each heat generation mechanism.\
  Default: false

- *"PostGraphHeatGenerationContributions"*:\
  Boolean for writing a graph with the contribution of each heat generation mechanism to the total heat generation.\
  Default: false

- *"PostHeatMapGeneration"*:\
  Boolean for assemblying and writing the heat map of heat generation.\
  Default: false

**Material properties**
- *"materials.Variables.THERMAL_CONDUCTIVITY"*:\
  Thermal conductivity of material (always required).

- *"materials.Variables.SPECIFIC_HEAT"*:\
  Specific heat of material (always required).

- *"materials.Variables.EMISSIVITY"*:\
  Emissivity of material surface, from 0.0 to 1.0, in which 1.0 is a perfect black body (required by radiative heat transfer).

- *"materials.Variables.THERMAL_EXPANSION_COEFFICIENT"*:\
  Thermal expansion coefficient of material (default: 0.0).

- *"materials.Tables.input_variable"*:\
  Input (X-axis) variable of the curve given by the table.

- *"materials.Tables.output_variable"*:\
  Output (Y-axis) variable of the curve given by the table.

- *"materials.Tables.data"*:\
  Points of the curve given by the table.

- *"material_relations.Variables.DYNAMIC_FRICTION"*:\
  Dynamic friction coefficient between materials (required by heat generation by sliding friction).

**SubModelPartData**
- *TEMPERATURE*:\
  Initial temperature of the elements.\
  Default: 0.0

- *FIXED_TEMPERATURE*:\
  Boolean for fixing the temperature of the elements.\
  Default: 0

- *HEATFLUX*:\
  Prescribed value of heat flux applied over elements surface (thermal power by area).\
  Default: 0.0

- *HEATSOURCE*:\
  Prescribed value of heat source applied over elements volume (thermal power by volume).\
  Default: 0.0

- *ADIABATIC*:\
  Boolean for setting elements as adiabatic (no heat transfer with other elements).\
  Default: 0

- *REAL_YOUNG_MODULUS_RATIO*:\
  Ratio between the real value of the Young modulus and the value used in the simulation.\
  Default: 1.0

## Testing

To test if the application is working correctly, run the *test_ThermalDEMApplication.py* file, located in the *tests* folder.
