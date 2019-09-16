import KratosMultiphysics
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.MeshingApplication as KratosMesh

from geo_processor import GeoProcessor

class GeoModel( GeoProcessor ):

    def __init__( self ):
        super(GeoMesher, self).__init__()

        self.HasModelPart = False
        self.HasExtrusionHeight = False

    def DefineElementType( self, type_of_element ):

        if ( not self.HasModelPart ):
            KratosMultiphysics.Logger.PrintWarning("GeoModel", "No model part has been set for this function")
            return



        # for elem in self.ModelPart: