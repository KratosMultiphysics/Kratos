from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMesh
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
from KratosMultiphysics.gid_output_process import GiDOutputProcess

def Factory(settings, Model): 
    if(type(settings) != KratosMultiphysics.Parameters): 
        raise Exception("expected input shall be a Parameters object, encapsulating a json string") 
    return ApplyErrorCalculatorProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process" 
class ApplyErrorCalculatorProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(""" 
            {
                "model_part_name"                       : "MainModelPart",
                "minimal_size"                        : 0.005,
                "maximal_size"                        : 2.0,
                "nodal_averaging"                     : true,
                "reference_variable_name"             : "ERROR_RATIO",
                "historical_results"                  : true,
                "echo_level"                          : 0
            }
            """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.historical_results = settings["historical_results"].GetBool()
        print(self.model_part.NumberOfNodes())
        print(self.model_part.NumberOfElements())


        #Process Parameters
        error_calc_params = KratosMultiphysics.Parameters("""{}""")
        error_calc_params.AddValue("minimal_size", settings["minimal_size"])
        error_calc_params.AddValue("maximal_size", settings["maximal_size"])
        error_calc_params.AddValue("nodal_averaging", settings["nodal_averaging"])
        error_calc_params.AddValue("reference_variable_name", settings["reference_variable_name"])
        error_calc_params.AddValue("echo_level", settings["echo_level"])
        error_calc_params.AddValue("historical_results", settings["historical_results"])

        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.error_calc = KratosConvDiff.MetricsTemperatureGradientProcess2D(self.model_part, error_calc_params)
        else:
            self.error_calc = KratosConvDiff.MetricsTemperatureGradientProcess3D(self.model_part, error_calc_params)
        
        remesh_params = KratosMultiphysics.Parameters("""{}""")
        if (self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.remesh_process = KratosMesh.MmgProcess2D(self.model_part, remesh_params)
        else:
            self.remesh_process = KratosMesh.MmgProcess3D(self.model_part, remesh_params)
    
    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        if self.historical_results:
            # We execute our process
            self.error_calc.Execute()
        # if self.model_part.ProcessInfo[KratosMultiphysics.TIME] == 10.0:
        #     self.remesh_process.Execute()

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        if not self.historical_results:
            # We execute our process
            self.error_calc.Execute()
        self.remesh_process.Execute()
        # Fix Model before Printing
        # self.model_part.RemoveSubModelPart("thermal_computing_domain")

        model_parts = self.model_part.SubModelParts

        for r_sub_model_part in model_parts:
            if r_sub_model_part.NumberOfElements() == 0 and r_sub_model_part.NumberOfConditions() == 0:
                r_sub_model_part.Nodes.clear()
                for node in self.model_part.Nodes:
                    r_sub_model_part.AddNode(node, 0)
        
        print(self.model_part)

        output_params = KratosMultiphysics.Parameters("""{
            "output_configuration"                  : {
                "result_file_configuration" : {
                    "gidpost_flags" : {
                        "GiDPostMode"           : "GiD_PostBinary",
                        "WriteDeformedMeshFlag" : "WriteDeformed",
                        "WriteConditionsFlag"   : "WriteConditions",
                        "MultiFileFlag"         : "SingleFile"
                    },
                    "file_label"          : "time",
                    "output_control_type" : "time",
                    "output_frequency"    : 0.1,
                    "body_output"         : true,
                    "node_output"         : false,
                    "skin_output"         : false,
                    "nodal_results"       : []
                },
                "point_data_configuration"  : []
            }
        }""")

        KratosMultiphysics.ModelPartIO("Revised_fluid_buoyancy", KratosMultiphysics.IO.WRITE | KratosMultiphysics.IO.MESH_ONLY).WriteModelPart(self.model_part)
        gid_output = GiDOutputProcess(self.model_part, "Revised_fluid_buoyancy", output_params["output_configuration"])
        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()

        vtk_settings = KratosMultiphysics.Parameters("""{
            "condition_data_value_variables": [],
            "condition_flags": [],
            "custom_name_postfix": "",
            "custom_name_prefix": "",
            "element_data_value_variables": [],
            "element_flags": [],
            "file_format": "ascii",
            "folder_name": "VTK_Output",
            "gauss_point_variables_extrapolated_to_nodes": [],
            "gauss_point_variables_in_elements": [],
            "model_part_name": "FluidPart",
            "nodal_data_value_variables": [],
            "nodal_flags": [],
            "nodal_solution_step_data_variables": ["VELOCITY","PRESSURE","TEMPERATURE","NODAL_TEMP_GRADIENT"],
            "output_control_type": "step",
            "output_frequency": 1.0,
            "output_precision": 7,
            "output_sub_model_parts": false,
            "save_output_files_in_folder": false,
            "write_deformed_configuration": false,
            "write_ids": false
        }""")
        #vtk_settings["model_part_name"] = str(model_input)
        vtk_io = KratosMultiphysics.VtkOutput(self.model_part, vtk_settings)
        vtk_io.PrintOutput()

    def Clear(self):
        pass



