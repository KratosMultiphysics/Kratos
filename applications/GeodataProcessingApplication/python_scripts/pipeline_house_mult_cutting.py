import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo

### --- An exemplaric pipeline for the meshing a town--- ###

# REMARK: The building are process one after another and cut
# out immediately. This reduces memory consumption.

from geo_importer import GeoImporter
from geo_mesher import GeoMesher
from geo_building import GeoBuilding

### IMPORTER ### -------------------------------------------------

tool1 = GeoImporter()
tool1.XyzImport( "data/toy_terrain.xyz" )

tool1.ShowModelPartQuality()

model_part = tool1.GetGeoModelPart()

### TERRAIN MESHING ### ------------------------------------------

tool2 = GeoMesher()
tool2.SetGeoModelPart( model_part )

tool2.ComputeExtrusionHeight( 5.0, 4.0, 1.0, 10 )

tool2.MeshConcaveHullWithTerrainPoints()

# 1st iteration
tool2.ComputeDistanceFieldFromGround()
tool2.RefineMeshNearGround( 7.0 )
# 2nd iteration
tool2.ComputeDistanceFieldFromGround()
tool2.RefineMeshNearGround( 5.0 )
# Ready
tool2.CreateGidControlOutput( "out_of_pipeline_1" )
tool2.CleanIsolatedNodes()
tool2.ShowModelPartQuality()
model_part = tool2.GetGeoModelPart()

### BUILDING MESHING ### ------------------------------------------

tool3 = GeoBuilding()
tool3.SetGeoModelPart( model_part )

tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
tool3.ShiftBuildingHull(1, 0.0, 0 )

tool3.ComputeDistanceFieldFromHull(False, 0.3)
tool3.SubtractBuilding( 0.05, 0.5, 0.05)
tool3.CreateGidControlOutput( "out_of_pipeline_2a" )

tool3.ComputeDistanceFieldFromHull()
tool3.SubtractBuilding( 0.1, 0.5, 0.05)

tool3.CreateGidControlOutput( "out_of_pipeline_2b" )
model_part = tool3.GetGeoModelPart()

### Adding another building --------------------------------------

tool3 = GeoBuilding()
tool3.SetGeoModelPart( model_part )

tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
tool3.ShiftBuildingHull(3.6, 0.0, 0 )

tool3.ComputeDistanceFieldFromHull(False, 0.3)
tool3.SubtractBuilding( 0.05, 0.5, 0.05)

tool3.ComputeDistanceFieldFromHull()
tool3.SubtractBuilding( 0.1, 0.5, 0.1)

tool3.CreateGidControlOutput( "out_of_pipeline_3" )
model_part = tool3.GetGeoModelPart()

### Adding another building --------------------------------------

tool3 = GeoBuilding()
tool3.SetGeoModelPart( model_part )

tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
tool3.ShiftBuildingHull(9, 0.3, 0.7 )

tool3.ComputeDistanceFieldFromHull(False, 0.3)
tool3.SubtractBuilding( 0.05, 0.5, 0.05)

tool3.ComputeDistanceFieldFromHull()
tool3.SubtractBuilding( 0.1, 0.5, 0.1)

tool3.CreateGidControlOutput( "out_of_pipeline_4" )
model_part = tool3.GetGeoModelPart()

print( model_part )

print( "CP - final" )