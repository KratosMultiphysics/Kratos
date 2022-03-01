import KratosMultiphysics as Kratos
import KratosMultiphysics.MeshingApplication as KratosMesh
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo

from KratosMultiphysics.GeodataProcessingApplication.geo_processor import GeoProcessor

import os
import time

class GeoDebug( GeoProcessor ):

    def __init__(self):
        super(GeoDebug, self).__init__()
        
        self.start_time = time.time()
        self.curr_time = time.time()
        self.dict_time = dict()         # dictionary with intervalls. key=name_intervall; value=time


    def ExportSubModelPart(self, dir_out, gid_post_mode="GiD_PostBinary"):
        if not self.HasModelPart:
            Kratos.Logger.PrintWarning("GeoDebug", "Model part has to be set, first.")
            return
        
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


    def PrintTime(self, msg="debug"):
        """
            time between steps
        """
        print("[DEBUG] {} OK in {} seconds -> {} minutes".format(msg, time.time()-self.curr_time, (time.time()-self.curr_time)/60))
        self.curr_time = time.time()


    def PrintTimeIntervals(self, msg="debug", time_end=None, time_start=None):
        """
            time between 2 intervals
        """
        # if (time_start == None): time_start = self.start_time
        time_start = self.start_time if (time_start == None) else time_start
        
        # if (time_end == None): time_end = time.time()
        time_end = time.time() if (time_end == None) else time_end
        
        print("[DEBUG] {} OK in {} seconds -> {} minutes".format(msg, time_end-time_start, (time_end-time_start)/60))
        self.curr_time = time.time()


    def SaveIntervall(self, name="interval"):
        """
            function to save intervals
        """
        self.dict_time[name] = time.time()
