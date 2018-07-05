from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

def CreateSolver(model, custom_settings):

#    if (type(model) != KratosMultiphysics.Model):
#        raise Exception("input is expected to be provided as a Kratos Model object")

#    if (type(custom_settings) != KratosMultiphysics.Parameters):
#        raise Exception("input is expected to be provided as a Kratos Parameters object")

    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@") 
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@") 
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@") 
    print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@") 
    print(custom_settings) 
    
     
    parallelism = "OpenMP"
    solver_type = "Transient"
    default_settings = KratosMultiphysics.Parameters("""
            {
                "solver_type" : "ThermallyCoupled",
                "fluid_solver_settings": {
                        "solver_type": "navier_stokes_solver_vmsmonolithic",
                        "model_import_settings": {
                                "input_type": "mdpa",
                                "input_filename": "unknown_name"
                        }
                },
                "thermal_solver_settings": {
                        "solver_type": "Transient",
                        "analysis_type": "linear",
                        "model_import_settings": {
                                "input_type": "use_input_model_part",
                                "input_filename": "unknown_name"
                        },
                        "computing_model_part_name": "Thermal",
                        "material_import_settings": {
                                "materials_filename": "ThermicMaterials.json"
                        },
                        "convection_diffusion_variables": {
                                "density_variable": "DENSITY",
                                "diffusion_variable": "CONDUCTIVITY",
                                "unknown_variable": "TEMPERATURE",
                                "volume_source_variable": "HEAT_FLUX",
                                "surface_source_variable": "FACE_HEAT_FLUX",
                                "projection_variable": "PROJECTED_SCALAR1",
                                "convection_variable": "CONVECTION_VELOCITY",
                                "mesh_velocity_variable": "MESH_VELOCITY",
                                "transfer_coefficient_variable": "",
                                "velocity_variable": "VELOCITY",
                                "specific_heat_variable": "SPECIFIC_HEAT",
                                "reaction_variable": "REACTION_FLUX"
                        }
                }
            }
            """)
    print(parallelism) 
    print(solver_type) 
    
    # Solvers for OpenMP parallelism
    if (parallelism == "OpenMP"):
        if (solver_type == "transient" or solver_type == "Transient"):
            solver_module_name = "convection_diffusion_transient_solver"
        elif (solver_type == "stationary" or solver_type == "Stationary"):
            solver_module_name = "convection_diffusion_stationary_solver"
        else:
            err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            err_msg += "Available options are: \"transient\", \"stationary\""
            raise Exception(err_msg)

    # Solvers for MPI parallelism
    elif (parallelism == "MPI"):
        err_msg =  "The requested parallel type MPI is not yet available!\n"
        raise Exception(err_msg)
        #if (solver_type == "transient" or solver_type == "Transient"):
            #solver_module_name = "trilinos_convection_diffusion_transient_solver"

        #elif (solver_type == "stationary" or solver_type == "Stationary"):
            #solver_module_name = "trilinos_convection_diffusion_stationary_solver"

        #else:
            #err_msg =  "The requested solver type \"" + solver_type + "\" is not in the python solvers wrapper\n"
            #err_msg += "Available options are: \"transient\", \"stationary\""
            #raise Exception(err_msg)
    else:
        err_msg =  "The requested parallel type \"" + parallelism + "\" is not available!\n"
        err_msg += "Available options are: \"OpenMP\", \"MPI\""
        raise Exception(err_msg)

    # Remove settings that are not needed any more
    #custom_settings["solver_settings"].RemoveValue("solver_type")
    #custom_settings["solver_settings"].RemoveValue("time_integration_method") # does not throw even if the value is not existing
    
    solver_module = __import__(solver_module_name)
    solver = solver_module.CreateSolver(model, default_settings["thermal_solver_settings"])

    return solver
