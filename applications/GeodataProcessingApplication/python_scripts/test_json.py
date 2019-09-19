"""
# IT'S ONLY A TEST. SCRIPT NOT COMPLETE 
"""


import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication

from geo_importer import GeoImporter
from geo_mesher import GeoMesher
from geo_preprocessor import GeoPreprocessor
from geo_building import GeoBuilding


# we read the JSON file for the CFD analysis
input_string = open("data/parameters/test/ProjectParameters.json",'r').read()
settings = KratosMultiphysics.Parameters(input_string)
print("01")

# we set "problem_name"
problem_name = "test_00"
settings["problem_data"].AddEmptyValue("problem_name")
settings["problem_data"]["problem_name"].SetString(problem_name)

print("02")

# we set "output_name"
problem_name = "ABC"
gid_output = settings["output_processes"]["gid_output"]		# gid_output is an array with only one element (where there all the entries are)
gid_output[0]["Parameters"].AddEmptyValue("output_name")
gid_output[0]["Parameters"]["output_name"].SetString(problem_name)

print("03")

# we set "model_import_settings"
input_filename = "CUSTOM_NAME"
settings["solver_settings"]["model_import_settings"].AddEmptyValue("input_filename")
settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(input_filename)

print("04")

# we set "volume_model_part_name"
volume_name = "Parts_Fluid"
settings["solver_settings"].AddEmptyValue("volume_model_part_name")
settings["solver_settings"]["volume_model_part_name"].SetString(volume_name)

print("05")

# we set "skin_parts"
skin_name_custom = ["Inlet", "Outlet", "BottomModelPart", "TopModelPart", "SKIN_ISOSURFACE"]
settings["solver_settings"].RemoveValue("skin_parts")	# to avoid duplicates
settings["solver_settings"].AddEmptyArray("skin_parts")

for i in range(len(skin_name_custom)):
	settings["solver_settings"]["skin_parts"].Append(skin_name_custom[i])

print("06")


# we set "boundary_conditions_process_list" array
inlet_name = "INLET"
outlet_name = "OUTLET"
slip_name = ["Slip", "Slip_2"]		# it can be a string or a list of strings
noslip_name = "NoSlip"				# it can be a string or a list of strings

settings["processes"].RemoveValue("boundary_conditions_process_list")	# to avoid duplicates
settings["processes"].AddEmptyArray("boundary_conditions_process_list")
# inlet
inlet_process_1 = KratosMultiphysics.Parameters("""
				{
					"python_module" : "apply_inlet_process",
					"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
					"Parameters"    : {
						"model_part_name" : "FluidModelPart.""" + str(inlet_name) + """",
						"variable_name"   : "VELOCITY",
						"modulus"         : 1.0,
						"direction"       : "automatic_inwards_normal",
						"interval"        : [0,1]
					}
				} """)
settings["processes"]["boundary_conditions_process_list"].Append(inlet_process_1)

inlet_process_2 = KratosMultiphysics.Parameters("""
				{
					"python_module" : "apply_inlet_process",
					"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
					"Parameters"    : {
						"model_part_name" : "FluidModelPart.""" + str(inlet_name) + """",
						"variable_name"   : "VELOCITY",
						"modulus"         : 0.5,
						"direction"       : "automatic_inwards_normal",
						"interval"        : [1,"End"]
					}
				} """)
settings["processes"]["boundary_conditions_process_list"].Append(inlet_process_2)

# outlet
outlet_process = KratosMultiphysics.Parameters("""
				{
					"python_module" : "apply_outlet_process",
					"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
					"Parameters"    : {
						"model_part_name" : "FluidModelPart.""" + str(outlet_name) + """",
						"variable_name"   : "PRESSURE",
						"constrained"        : true,
						"value"              : 0.0,
						"hydrostatic_outlet" : false,
						"h_top"              : 0.0
					}
				} """)
settings["processes"]["boundary_conditions_process_list"].Append(outlet_process)

# slip (it can be a string or a list of strings)
if (isinstance(slip_name, list)):
	for i_slip in slip_name:
		slip_process = KratosMultiphysics.Parameters("""
						{
							"python_module" : "apply_slip_process",
							"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
							"Parameters"    : {
								"model_part_name" : "FluidModelPart.""" + str(i_slip) + """"
							}
						} """)
		settings["processes"]["boundary_conditions_process_list"].Append(slip_process)
elif (isinstance(slip_name, str)):
	slip_process = KratosMultiphysics.Parameters("""
					{
						"python_module" : "apply_slip_process",
						"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
						"Parameters"    : {
							"model_part_name" : "FluidModelPart.""" + str(slip_name) + """"
						}
					} """)
	settings["processes"]["boundary_conditions_process_list"].Append(slip_process)

# noslip (it can be a string or a list of strings)
if (isinstance(noslip_name, list)):
	for i_noslip in noslip_name:
		noslip_process = KratosMultiphysics.Parameters("""
						{
							"python_module" : "apply_noslip_process",
							"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
							"Parameters"    : {
								"model_part_name" : "FluidModelPart.""" + str(i_noslip) + """"
							}
						} """)
		settings["processes"]["boundary_conditions_process_list"].Append(noslip_process)
elif (isinstance(noslip_name, str)):
	noslip_process = KratosMultiphysics.Parameters("""
					{
						"python_module" : "apply_noslip_process",
						"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
						"Parameters"    : {
							"model_part_name" : "FluidModelPart.""" + str(noslip_name) + """"
						}
					} """)
	settings["processes"]["boundary_conditions_process_list"].Append(noslip_process)

# # we set processes/boundary_conditions_process_list/Parameters/model_part_name (Inlet)
# inlet_name = "AAAAAAAAA"
# bc_list = settings["processes"]["boundary_conditions_process_list"]		# bc_list is an array with only one element (where there all the entries are)
# bc_list[0]["Parameters"].AddEmptyValue("model_part_name")
# bc_list[0]["Parameters"]["model_part_name"].SetString("FluidModelPart.{}".format(inlet_name))

# bc_list[1]["Parameters"].AddEmptyValue("model_part_name")
# bc_list[1]["Parameters"]["model_part_name"].SetString("FluidModelPart.{}".format(inlet_name))

# outlet_name = "BBBBBBBB"
# bc_list[2]["Parameters"].AddEmptyValue("model_part_name")
# bc_list[2]["Parameters"]["model_part_name"].SetString("FluidModelPart.{}".format(outlet_name))

# print(settings.PrettyPrintJsonString())
f = open("data/parameters/test/ProjectParameters_MOD.json",'w')
f.write(settings.PrettyPrintJsonString())
f.close()