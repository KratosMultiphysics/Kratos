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
tool1.XyzImport( "dataTown/winterberg_clipped_2-5m.xyz" )
# tool1.ShowModelPartQuality()
model_part = tool1.GetGeoModelPart()

### TERRAIN MESHING ### ------------------------------------------

tool2 = GeoMesher()
tool2.SetGeoModelPart( model_part )
tool2.SetUniformExtrusionHeight( 710.0 )
tool2.MeshConcaveHullWithTerrainPoints( 3.0 )

# 1st iteration
tool2.ComputeDistanceFieldFromGround()
# tool2.RefineMeshNearGround( 40.0 )
# 2nd iteration
# tool2.ComputeDistanceFieldFromGround()
# tool2.RefineMeshNearGround( 10.0 )
# Ready
tool2.CreateGidControlOutput( "town_1" )
tool2.CleanIsolatedNodes()
tool2.ShowModelPartQuality()
model_part = tool2.GetGeoModelPart()

### BUILDING MESHING ### ------------------------------------------

# tool3 = GeoBuilding()
# tool3.SetGeoModelPart( model_part )

# tool3.ImportBuildingHullMDPA( "/dataTown/MdpaFromObj_10" )
# tool3.ShiftBuildingHull(1, 0.0, 0 )

# tool3.ComputeDistanceFieldFromHull(False, 1)
# tool3.CreateGidControlOutput( "town_2" )

# tool3.SubtractBuilding( 1.0, 7.0, 0.5)

# tool3.ComputeDistanceFieldFromHull()
# tool3.CreateGidControlOutput( "town_2a" )
# tool3.SubtractBuilding( 2.0, 5.0, 0.3 )

# tool3.CreateGidControlOutput( "town_2b" )
# model_part = tool3.GetGeoModelPart()

### Adding another building --------------------------------------

tool3 = GeoBuilding()
tool3.SetGeoModelPart( model_part )

tool3.ImportBuildingHullMDPA( "/dataTown/MdpaFromObj_4" )
# tool3.ShiftBuildingHull(1, 0.0, 0 )

tool3.ComputeDistanceFieldFromHull(False, 1)
tool3.CreateGidControlOutput( "town_3a" )

tool3.SubtractBuilding( 2.0, 12.0, 0.5)
tool3.CreateGidControlOutput( "town_3b" )

tool3.ComputeDistanceFieldFromHull()
tool3.SubtractBuilding( 2.0, 5.0, 0.3 )

tool3.CreateGidControlOutput( "town_3c" )
model_part = tool3.GetGeoModelPart()

### Adding another building --------------------------------------

tool3 = GeoBuilding()
tool3.SetGeoModelPart( model_part )

tool3.ImportBuildingHullMDPA( "/dataTown/MdpaFromObj_10" )
# tool3.ShiftBuildingHull(1, 0.0, 0 )

tool3.ComputeDistanceFieldFromHull(False, 1)
tool3.SubtractBuilding( 2.0, 12.0, 0.5)
tool3.CreateGidControlOutput( "town_4a" )

tool3.ComputeDistanceFieldFromHull()
tool3.SubtractBuilding( 2.0, 5.0, 0.3 )

tool3.CreateGidControlOutput( "town_4b" )
model_part = tool3.GetGeoModelPart()

### Adding another building --------------------------------------

tool3 = GeoBuilding()
tool3.SetGeoModelPart( model_part )

tool3.ImportBuildingHullMDPA( "/dataTown/MdpaFromObj_2" )
# tool3.ShiftBuildingHull(1, 0.0, 0 )

tool3.ComputeDistanceFieldFromHull(False, 1)
tool3.SubtractBuilding( 2.0, 12.0, 0.5)
tool3.CreateGidControlOutput( "town_5a" )

tool3.ComputeDistanceFieldFromHull()
tool3.SubtractBuilding( 2.0, 5.0, 0.3 )

tool3.CreateGidControlOutput( "town_5b" )
model_part = tool3.GetGeoModelPart()

### Adding another building --------------------------------------

tool3 = GeoBuilding()
tool3.SetGeoModelPart( model_part )

tool3.ImportBuildingHullMDPA( "/dataTown/MdpaFromObj_5" )
# tool3.ShiftBuildingHull(1, 0.0, 0 )

tool3.ComputeDistanceFieldFromHull(False, 1)
tool3.SubtractBuilding( 2.0, 12.0, 0.5)
tool3.CreateGidControlOutput( "town_6a" )

tool3.ComputeDistanceFieldFromHull()
tool3.SubtractBuilding( 2.0, 5.0, 0.3 )

tool3.CreateGidControlOutput( "town_6b" )
model_part = tool3.GetGeoModelPart()


print( model_part )

print( "CP - final" )