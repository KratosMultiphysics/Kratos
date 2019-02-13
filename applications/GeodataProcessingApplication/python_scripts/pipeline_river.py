import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication
import os

### --- An exemplaric pipeline for the meshing a river bed --- ###

from geo_importer import GeoImporter
from geo_mesher import GeoMesher
from geo_building import GeoBuilding

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

### IMPORTER ### -------------------------------------------------

tool1 = GeoImporter()
# tool1.xyz_import( "data/torrent1_qgis.xyz" )
tool1.xyz_import( "data/torrent2_qgis.xyz" )

tool1.ShowModelPartQuality()

model_part = tool1.GetGeoModelPart()

### TERRAIN MESHING ### ------------------------------------------

tool2 = GeoMesher()
tool2.SetGeoModelPart( model_part )

### for section 1 (torrent1_qgis.xyz)
# tool2.ComputeExtrusionHeight( 13, 8, 1, 10 )

### for section 1 (torrent1_qgis.xyz)
tool2.ComputeExtrusionHeight( 7, 10, 1, 10 )

tool2.MeshConcaveHullWithTerrainPoints()

tool2.CreateGidControlOutput( "out_of_pipeline_0" )

# 1st iteration
tool2.ComputeDistanceFieldFromGround()
tool2.RefineMeshNearGround( 15.0 )
print( "End of Refinement 1")
# 2nd iteration
# tool2.ComputeDistanceFieldFromGround()
# tool2.RefineMeshNearGround( 8.0 )
# Ready
tool2.CreateGidControlOutput( "out_of_pipeline_1" )

tool2.CleanIsolatedNodes()
tool2.ShowModelPartQuality()

tool2.CreateGidControlOutput( "out_of_pipeline_2" )

main_model_part = tool2.GetGeoModelPart()
print( "Generation of main_model_part")
print( main_model_part )
print( "Generation of main_model_part")

###################################################################################################

current_model = Kratos.Model()

if current_model.HasModelPart( "NewModelPart" ):
    # clear the existing model part
    new_model_part = current_model.GetModelPart( "NewModelPart" )
    new_model_part.Elements.clear()
    new_model_part.Conditions.clear()
    new_model_part.Nodes.clear()
else:
    new_model_part = current_model.CreateModelPart( "NewModelPart" )

new_model_part.AddProperties(Kratos.Properties(1))
new_model_part.AddProperties(Kratos.Properties(2))

fluidSubModelPart = new_model_part.CreateSubModelPart("Parts_Fluid")
wallSubModelPart = new_model_part.CreateSubModelPart("Wall")
inletSubModelPart = new_model_part.CreateSubModelPart("Inlet")
outletSubModelPart = new_model_part.CreateSubModelPart("Outlet")
print( "Generation of Sub Model Parts finsished")

# finding only nodes that belong to an element
required_node_list = []
for elem in main_model_part.Elements:
    nodes = elem.GetNodes()
    required_node_list.extend( [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id] )

print( "Required node list finished")

# avoiding redundancies
required_node_list = list( set( required_node_list ) )
for n_id in required_node_list:
    node = main_model_part.GetNode( n_id )
    n = new_model_part.CreateNewNode( node.Id, node.X, node.Y, node.Z )
    fluidSubModelPart.AddNode(n, 0)
    if node.Y < 0.1:
        inletSubModelPart.AddNode(n, 0)
    if node.Y > 99.9:
        outletSubModelPart.AddNode(n, 0 )

print( "Node distribution finished")

for elem in main_model_part.Elements:
    nodes = elem.GetNodes()
    e = new_model_part.CreateNewElement("Element3D4N", elem.Id,  [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id],new_model_part.GetProperties()[1])
    fluidSubModelPart.AddElement( e, 0 )

print( "ELement distribution finished")

for cond in main_model_part.Conditions:
    nodes = cond.GetNodes()
    c = new_model_part.CreateNewCondition("WallCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id],new_model_part.GetProperties()[0])
    print( cond.Id )

    if ( nodes[0].Y < 0.1 and nodes[1].Y < 0.1 and nodes[2].Y < 0.1 ):
        inletSubModelPart.AddCondition( c, 0 )

    elif ( nodes[0].Y > 99.9 and nodes[1].Y > 99.9 and nodes[2].Y > 99.9 ):
        outletSubModelPart.AddCondition( c, 0 )

    else:
        wallSubModelPart.AddCondition( c, 0 )
        for n in c.GetNodes():
            wallSubModelPart.AddNode( n, 0 )

print( "Condition distribution finished")

model_part_io = Kratos.ModelPartIO("out", Kratos.IO.WRITE)
model_part_io.WriteModelPart( new_model_part )


### BUILDING MESHING ### ------------------------------------------

# tool3 = GeoBuilding()
# tool3.SetGeoModelPart( model_part )

# tool3.ImportBuildingHullMDPA( "/data/toy_house_skin" )
# # 1st iteration
# tool3.ComputeDistanceFieldFromHull()
# tool3.RefineMeshNearBuilding( 0.2 )
# # 2nd iteration
# # tool3.ComputeDistanceFieldFromHull()
# # tool3.RefineMeshNearBuilding( 0.05 )
# # 3rd iteration
# # tool3.ComputeDistanceFieldFromHull()
# # tool3.RefineMeshNearBuilding( 0.03 )
# # Ready
# tool3.ShowModelPartQuality()
# tool3.ComputeDistanceFieldFromHull()
# tool3.SubtractBuilding( 0.1, 0.5, 0.1)

# tool3.CreateGidControlOutput( "out_of_pipeline_2" )
# model_part = tool3.GetGeoModelPart()

# model = Kratos.Model()
# new_model_part = model.CreateModelPart("new_model_part")

# file_path = os.path.dirname(os.path.realpath(__file__))
# Kratos.ModelPartIO( file_path + "/out" ).ReadModelPart( new_model_part )

# print( new_model_part )

print( "CP - final" )