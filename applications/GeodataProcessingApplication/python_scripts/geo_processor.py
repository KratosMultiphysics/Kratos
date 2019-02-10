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


    def CreateGidControlOutput( self, output_name ):
        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(
                self.ModelPart,
	            output_name,
	            Kratos.Parameters("""
	                {
	                    "result_file_configuration" : {
	                        "gidpost_flags": {
	                            "GiDPostMode": "GiD_PostBinary",
	                            "MultiFileFlag": "SingleFile"
	                        },
	                        "nodal_results"       : ["DISTANCE","DISTANCE_GRADIENT"],
	                        "nodal_nonhistorical_results": [],
	                        "nodal_flags_results": []
	                    }
	                }
	                """)
	            )
        gid_output.ExecuteInitialize()
        gid_output.ExecuteBeforeSolutionLoop()
        gid_output.ExecuteInitializeSolutionStep()
        gid_output.PrintOutput()
        gid_output.ExecuteFinalizeSolutionStep()
        gid_output.ExecuteFinalize()