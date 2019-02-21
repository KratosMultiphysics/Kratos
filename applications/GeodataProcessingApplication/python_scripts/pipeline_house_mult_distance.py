import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo

### --- An exemplaric pipeline for the meshing a town--- ###

# REMARK: The buildings are summarized in one common distance field.
# Only in a final step, the building volumes are substracted all together

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
tool3.ShiftBuildingHull(1, 0, 0 )
tool3.ComputeDistanceFieldFromHull()

### Adding another building --------------------------------------

tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
tool3.ShiftBuildingHull(5, 0, 0 )
tool3.AddDistanceFieldFromHull()

### Adding another building --------------------------------------

tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
tool3.ShiftBuildingHull(9, 0.0, 0 )
tool3.AddDistanceFieldFromHull()

### Refinement for all -------------------------------------------

tool3.CreateGidControlOutput( "AllDist1" )
tool3.RefineMeshNearBuilding( 0.2 )

### Resetting complete distance ----------------------------------

tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
tool3.ShiftBuildingHull(1, 0, 0 )
tool3.ComputeDistanceFieldFromHull()

### Adding another building --------------------------------------

tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
tool3.ShiftBuildingHull(5, 0, 0 )
tool3.AddDistanceFieldFromHull()

### Adding another building --------------------------------------

tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
tool3.ShiftBuildingHull(9, 0.0, 0 )
tool3.AddDistanceFieldFromHull()

tool3.CreateGidControlOutput( "AllDist2" )

tool3.SubtractBuilding( 0.1, 0.5, 0.05)

tool3.CreateGidControlOutput( "AllDist3" )
model_part = tool3.GetGeoModelPart()

### Adding another building --------------------------------------

print( model_part )

print( "CP - final" )