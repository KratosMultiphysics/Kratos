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
        self.HasModelPart

    # get geo model part
    def GetGeoModelPart( self ):
        if self.HasModelPart:
            return self.ModelPart
        else:
            Kratos.Logger.PrintWarning("GeoProcessor", "No model part can be returned")


    # operation on the model part (cleaning, information etc.)

    def ShowModelPartQuality( self ):
	    # only to be applied if the more elegant C++ routine should have problems
        required_node_list = []

        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()
            required_node_list.extend( [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id] )

        # avoiding redundancies
        required_node_list = list( set( required_node_list ) )
        difference = len( self.ModelPart.Nodes ) - len( required_node_list )

        Kratos.Logger.PrintWarning("ModelPartQuality", str( len( self.ModelPart.Nodes ) ) + " nodes \n")
        Kratos.Logger.PrintWarning("ModelPartQuality", str( len( self.ModelPart.Conditions ) ) + " conditions \n")
        Kratos.Logger.PrintWarning("ModelPartQuality", str( len( self.ModelPart.Elements ) ) + " elements \n")

        if difference > 0:
            Kratos.Logger.PrintWarning("ModelPartQuality", str( difference) + " nodes are present that do not belong to an element! Critical. \n")


    def CleanIsolatedNodes( self ):
        # execute the utility for the sake of performance
        utility = KratosGeo.CleaningUtilities( self.ModelPart )
        utility.CleanIsolatedNodes()