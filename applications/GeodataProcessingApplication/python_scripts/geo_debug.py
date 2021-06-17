import KratosMultiphysics as Kratos
import KratosMultiphysics.MeshingApplication as KratosMesh
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo

from geo_processor import GeoProcessor
import os

class GeoDebug( GeoProcessor ):

    def __init__(self):
        pass


    def ExportSubModelPart(self, dir_out, gid_post_mode="GiD_PostBinary"):
        for smp in self.ModelPart.SubModelParts:
            print("[DEBUG] SubModelPart name: ", smp.Name)

            output_name = dir_out + smp.Name
            self.CreateGidOutput(output_name, smp, gid_post_mode)


    def CreateGidOutput( self, output_name, model_part, gid_post_mode="GiD_PostBinary"):
        from KratosMultiphysics.gid_output_process import GiDOutputProcess
        gid_parameters = Kratos.Parameters("""
            {
                "result_file_configuration" : {
                    "gidpost_flags": {
                        "WriteDeformedMeshFlag": "WriteDeformed",
                        "WriteConditionsFlag": "WriteConditions",
                        "MultiFileFlag": "SingleFile"
                    },
                    "node_output"         : true,
                    "nodal_results"       : ["DISTANCE","DISTANCE_GRADIENT"],
                    "nodal_nonhistorical_results": [],
                    "nodal_flags_results": []
                }
            }
            """)
        gid_parameters["result_file_configuration"]["gidpost_flags"].AddEmptyValue("GiDPostMode")
        gid_parameters["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].SetString(gid_post_mode)      # the options are: GiD_PostAscii / GiD_PostBinary
        
        gid_output = GiDOutputProcess(model_part, output_name, gid_parameters)
        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()


    def PrintInfo(self):
        pass
