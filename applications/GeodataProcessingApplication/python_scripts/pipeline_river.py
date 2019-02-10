import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo

### --- An exemplaric pipeline for the meshing a river bed --- ###

from geo_importer import GeoImporter
from geo_mesher import GeoMesher

tool1 = GeoImporter()
tool1.xyz_import( "data/toy_terrain.xyz" )
tool1.ShowModelPartQuality()
model_part = tool1.GetGeoModelPart()

tool2 = GeoMesher()
tool2.SetGeoModelPart( model_part )
tool2.ComputeExtrusionHeight( 5.0, 4.0, 1.0, 10 )
tool2.MeshConcaveHullWithTerrainPoints()
tool2.ShowModelPartQuality()
tool2.CreateGidControlOutput( "out_of_pipeline_1" )