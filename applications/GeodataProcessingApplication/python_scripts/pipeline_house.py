import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo

### --- An exemplary pipeline a single toy house --- ###

from geo_importer import GeoImporter
from geo_mesher import GeoMesher
from geo_building import GeoBuilding

### GEO IMPORTER ### -------------------------------------------------

tool1 = GeoImporter()
tool1.XyzImport( "data/toy_terrain.xyz" )
tool1.ShowModelPartQuality()
model_part = tool1.GetGeoModelPart()


### GEO MESHER ### ------------------------------------------

tool2 = GeoMesher()
tool2.SetGeoModelPart( model_part )
tool2.ComputeExtrusionHeight( 5.0, 4.0, 1.0, 10 )
tool2.MeshConcaveHullWithTerrainPoints( 1.0 )
tool2.CreateGidControlOutput( "terrain_1" )

# 1st ground refinement iteration
tool2.ComputeDistanceFieldFromGround()
tool2.RefineMeshNearGround( 7.0 )
tool2.CreateGidControlOutput( "terrain_2" )
# 2nd ground refinement iteration
tool2.ComputeDistanceFieldFromGround()
tool2.RefineMeshNearGround( 4.0 )
tool2.CreateGidControlOutput( "terrain_3" )
# terrain volume mesh ready
model_part = tool2.GetGeoModelPart()


### GEO BUILDING ### ------------------------------------------

strategy_options = ["RefineBeforeCutting", "PerformMultipleCuts"]
strategy = strategy_options[1]

tool3 = GeoBuilding()
tool3.SetGeoModelPart( model_part )
tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )

if ( strategy == "RefineBeforeCutting" ):
    # 1st zero-level refinement iteration
    tool3.ComputeDistanceFieldFromHull()
    tool3.CreateGidControlOutput( "refine_1a" )
    tool3.RefineMeshNearBuilding( 0.05 )
    tool3.CreateGidControlOutput( "refine_1b" )
    # 2nd zero-level refinement iteration
    tool3.ComputeDistanceFieldFromHull()
    tool3.CreateGidControlOutput( "refine_2a" )
    tool3.RefineMeshNearBuilding( 0.02 )
    tool3.CreateGidControlOutput( "refine_2b" )

if ( strategy == "PerformMultipleCuts" ):
    # 1st cutting iteration (shrinkage of hull by 0.2 m)
    tool3.ComputeDistanceFieldFromHull( False, 0.2 )
    tool3.CreateGidControlOutput( "cut_1a" )
    tool3.SubtractBuilding( 0.5, 1.0, 0.05 )
    tool3.CreateGidControlOutput( "cut_1b" )
    # 2nd cutting iteration (shrinkage of hull by 0.05 m)
    tool3.ComputeDistanceFieldFromHull( False, 0.05 )
    tool3.CreateGidControlOutput( "cut_2a" )
    tool3.SubtractBuilding( 0.05, 1.0, 0.01 )
    tool3.CreateGidControlOutput( "cut_2b" )

tool3.ShowModelPartQuality()
tool3.ComputeDistanceFieldFromHull()
tool3.CreateGidControlOutput( "result_distance" )
tool3.SubtractBuilding( 0.1, 0.5, 0.1)

if ( strategy == "RefineBeforeCutting" ):
    tool3.CreateGidControlOutput( "result_refine" )
if ( strategy == "PerformMultipleCuts" ):
    tool3.CreateGidControlOutput( "result_cut" )

model_part = tool3.GetGeoModelPart()
print( model_part )