from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

#import numpy
#from stl import mesh
import KratosMultiphysics
#import KratosMultiphysics.MeshingApplication as MeshingApplication
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import math

class SkinVisualization():
    
    def __init__(self, model1, model2, settings=KratosMultiphysics.Parameters("""{}""")):
        model_A = KratosMultiphysics.Model()
        self.model_fluid = model_A.CreateModelPart("FluidModelPart")
        self.model_fluid.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        self.model_fluid.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.model_fluid.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        KratosMultiphysics.ModelPartIO(model1).ReadModelPart(self.model_fluid)

        model_B = KratosMultiphysics.Model()
        self.model_solid = model_B.CreateModelPart("SolidModelPart")
        self.model_solid.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        #self.model_solid.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        #self.model_solid.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        KratosMultiphysics.ModelPartIO(model2).ReadModelPart(self.model_solid)

        model_C = KratosMultiphysics.Model()
        self.model_visualization = model_C.CreateModelPart("VisualizationModel")
        self.model_visualization.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        default_params = KratosMultiphysics.Parameters( """
        {
            "shape_functions"                       : "standard",
            "reform_model_part_at_each_time_step"   : false,
            "visualization_variables"               : ["DISTANCE","DISTANCE_GRADIENT"],
            "help"                                  : "This process detects the skin from a given submodelpart and it generates the correspong conditions",
            "model_part_name"                       : "SolidModelPart",
            "computing_model_part_name"             : "ComputingDomain",
            "recursive_detection"                   : true,
            "name_auxiliar_model_part"              : "SolidModelSkinPart",
            "name_auxiliar_condition"               : "Condition",
            "list_model_parts_to_assign_conditions" : [],
            "echo_level"                            : 0,
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
        } """ )

        settings.ValidateAndAssignDefaults(default_params)

        # Add the visualization model part variables to the visualization model part.
        # Add them to the nodal_results GiD output process list as well.
        for i_var in range(0, settings["visualization_variables"].size()):
            variable_name = settings["visualization_variables"][i_var].GetString()
            settings["output_configuration"]["result_file_configuration"]["nodal_results"].Append(variable_name)
            self.model_visualization.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable_name))

        # Set an auxilar Kratos Parameters object to build the skin visualization process
        aux_params = KratosMultiphysics.Parameters("""{}""")
        aux_params.AddValue("shape_functions",settings["shape_functions"])
        aux_params.AddValue("visualization_variables",settings["visualization_variables"])
        aux_params.AddValue("reform_model_part_at_each_time_step",settings["reform_model_part_at_each_time_step"])

        self.EmbeddedSkinVisualizationProcess = KratosFluid.EmbeddedSkinVisualizationProcess(
            self.model_fluid,
            self.model_visualization,
            aux_params)

        #self.computing_model_part = Model[settings["computing_model_part_name"].GetString()]
        self.recursive_detection = settings["recursive_detection"].GetBool()
        self.name_auxiliar_model_part = settings["name_auxiliar_model_part"].GetString()
        
        # Process parameters
        detect_skin_parameters = KratosMultiphysics.Parameters("""{}""")
        detect_skin_parameters.AddValue("name_auxiliar_model_part", settings["name_auxiliar_model_part"])
        detect_skin_parameters.AddValue("name_auxiliar_condition", settings["name_auxiliar_condition"])
        detect_skin_parameters.AddValue("list_model_parts_to_assign_conditions", settings["list_model_parts_to_assign_conditions"])
        detect_skin_parameters.AddValue("echo_level", settings["echo_level"])
        if (self.model_solid.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            self.detect_skin = KratosMultiphysics.SkinDetectionProcess2D(self.model_solid, detect_skin_parameters)
        else:
            self.detect_skin = KratosMultiphysics.SkinDetectionProcess3D(self.model_solid, detect_skin_parameters)
    
        self.gid_output = GiDOutputProcess(self.model_visualization, "VisualizationModel", settings["output_configuration"])
    
    def ExecuteInitialize(self):
         # We execute the process
        self.detect_skin.Execute()
        # We copy the conditions to the contact model part
        skin_model_part = self.model_solid.GetSubModelPart(self.name_auxiliar_model_part)
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(self.model_fluid, skin_model_part).Execute()
        #transfer_process = KratosMultiphysics.FastTransferBetweenModelPartsProcess(self.computing_model_part, skin_model_part, KratosMultiphysics.FastTransferBetweenModelPartsProcess.EntityTransfered.CONDITIONS)
        #transfer_process.Execute()
        
        self.EmbeddedSkinVisualizationProcess.ExecuteInitialize()
        self.gid_output.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        # Set time in case the GiD process control output is time
        self.model_visualization.ProcessInfo[KratosMultiphysics.TIME] = self.model_fluid.ProcessInfo[KratosMultiphysics.TIME]
        # If recurive we detect each time step
        if (self.recursive_detection is True):
            self.detect_skin.Execute()
            # We copy the conditions to the contact model part
            skin_model_part = self.model_solid.GetSubModelPart(self.name_auxiliar_model_part)
            KratosMultiphysics.CalculateDistanceToSkinProcess2D(self.model_fluid, skin_model_part).Execute()

        self.EmbeddedSkinVisualizationProcess.ExecuteInitializeSolutionStep()
        self.gid_output.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteInitializeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeOutputStep()
        if (self.gid_output.IsOutputStep()):
            self.gid_output.PrintOutput()

    def ExecuteFinalize(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalize()
        self.gid_output.ExecuteFinalize()


if __name__ == "__main__":
    vis_tool = SkinVisualization("cavity0_6_fluid_water", "solid_circle")
    vis_tool.ExecuteInitialize()
    vis_tool.ExecuteBeforeSolutionLoop()
    vis_tool.ExecuteInitializeSolutionStep()
    vis_tool.ExecuteFinalizeSolutionStep()
    vis_tool.ExecuteBeforeOutputStep()
    vis_tool.ExecuteFinalize()
