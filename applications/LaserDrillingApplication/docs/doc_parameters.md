In this file we describe the parameters that appear in the JSON files that determine the simulation parameters.

## ProjectParameters.json

- problem_name:
- parallel_type: Which kind of parallelisation form to use
- time_step (ms): Ellapsed time between each time step
- start_time (ms): The time at the start of the simulation 
- end_time (ms): The time at the end of the simulation
- average_laser_power (W): The average laser power over a period
- pulse_frequency (Hz): The amount of pulses per unit time. Also called the pulse repetition rate
- focus_Z_offset (mm): 
- Rayleigh_length (mm): The Rayleigh length of the gaussian mode of the laser
- beam_waist_diameter (mm): The diameter of the gaussian beam at the waist
- consider_material_refraction (bool): Toggle refraction at the interface between air and the sample
- adjust_T_field_after_ablation (bool):
- reference_T_after_laser (K):
- mesh_size (string): The size of the elements of the mesh of the included .mdpa geometry files. 
- mesh_type (string): The type (structured or unstructured) of the mesh of the included .mdpa geometry files.
- print_hole_geometry_files (bool): Toggle to print the hole geometry files.
- compute_vaporisation (bool): (unused) Toggle the computation of vaporisation due to the temperature of an element exceeding the vaporisation temprature of the material.
- vaporisation_temperature (K): The vaporisation temperature of the material.
- echo_level (int):
- print_debug_info (bool):

## LaserDrillingMaterials.json
- compute_optical_penetration_depth_using_refractive_index (bool): Toggle to choose between setting a value for the penetration depth as a parameter or computing it from the refractive index.
- compute_energy_per_unit_volume_threshold_using_enthalpy_and_ionization (bool):Toggle to choose between setting a value for the energy per unit volume threshold as a parameter or computing it from the enthalpy and the ionization energy.
- DENSITY:
- CONDUCTIVITY:
- SPECIFIC_HEAT:
- IONIZATION_ALPHA:
- PENETRATION_DEPTH:
- ABLATION_THRESHOLD:
- THERMAL_DEPTH:
- ENTHALPY:
- OPTICAL_PENETRATION_DEPTH:
- ENERGY_PER_VOLUME_THRESHOLD:
- REFRACTIVE_INDEX:
