import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo

# Superclass

class GeoProcessor:

    def __init__( self ):
        self.HasModelPart = False

    # set geo model part
    def SetGeoModelPart( self, modelPartIn ):
        if self.HasModelPart:
            Kratos.Logger.PrintWarning("GeoProcessor", "Model part was over-written")
        self.ModelPart = modelPartIn
        self.HasModelPart = True


    # get geo model part
    def GetGeoModelPart( self ):
        if self.HasModelPart:
            return self.ModelPart
        else:
            Kratos.Logger.PrintWarning("GeoProcessor", "No model part can be returned")


    # delete geo model part
    def DelGeoModelPart(self):
        if self.HasModelPart:
            model = self.ModelPart.GetModel()
            model.DeleteModelPart(self.ModelPart.Name)
            # self.ModelPart.DeleteModelPart("ModelPart")
            self.HasModelPart = False
        else:
            Kratos.Logger.PrintWarning("GeoProcessor", "No model part is present")


    #### operation on the model part (cleaning, information etc.)

    def ShowModelPartQuality( self ):
	    # only to be applied if the more elegant C++ routine should have problems
        required_node_list = []

        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()
            for node in nodes:
                required_node_list.append( node.Id )

        # avoiding redundancies
        required_node_list = list( set( required_node_list ) )
        difference = len( self.ModelPart.Nodes ) - len( required_node_list )

        Kratos.Logger.PrintWarning("ModelPartQuality", str( len( self.ModelPart.Nodes ) ) + " node(s)")
        Kratos.Logger.PrintWarning("ModelPartQuality", str( len( self.ModelPart.Conditions ) ) + " condition(s)")
        Kratos.Logger.PrintWarning("ModelPartQuality", str( len( self.ModelPart.Elements ) ) + " element(s)")

        if difference > 0:
            Kratos.Logger.PrintWarning("ModelPartQuality", str( difference) + " nodes are present that do not belong to an element! Critical. \n")
        elif difference == 0:
            Kratos.Logger.PrintWarning("ModelPartQuality", "No isolated nodes. Model part ready for further steps. \n")

    def CleanIsolatedNodes( self ):
        # execute the utility for the sake of performance
        utility = KratosGeo.CleaningUtilities( self.ModelPart )
        utility.CleanIsolatedNodes()


    def ReBuildModelPart( self ):
        # execute the utility for the sake of performance
        utility = KratosGeo.CleaningUtilities( self.ModelPart )
        self.ModelPart = utility.ReBuildModelPart()


    def CreateGidControlOutput( self, output_name, gid_post_mode="GiD_PostBinary" ):
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
        
        gid_output = GiDOutputProcess(self.ModelPart, output_name, gid_parameters)
        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()


    def WriteMdpaOutput( self, mdpa_file_name="model_part"):

        if self.HasModelPart:
            model_part = self.ModelPart
            model_part_io = Kratos.ModelPartIO(mdpa_file_name, Kratos.IO.WRITE)
            model_part_io.WriteModelPart( model_part )
        else:
            Kratos.Logger.PrintWarning("WriteMdpaOutput", "No model part is present. \n")

    def ReadMdpaInput( self, mdpa_file_name):

        if self.HasModelPart:
            Kratos.Logger.PrintWarning("ReadMdpaInput", "Reading an MDPA input file replaced the ModelPart  \n")
        else:
            Kratos.Logger.PrintWarning("ReadMdpaInput", "No model part is present. \n")
    

    # model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_write.out"), KratosMultiphysics.IO.WRITE)
    # model_part_io.WriteModelPart(model_part)