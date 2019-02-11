import KratosMultiphysics
from geo_processor import GeoProcessor
import os

class GeoBuilding( GeoProcessor ):

    def __init__( self ):
        super(GeoBuilding, self).__init__()

        self.HasBuildingHull = False


    def ImportBuildingHullSTL( self, file_name ):

        print("To be done...")


    def ImportBuildingHullMDPA( self, file_name ):

        self._generate_building_model_part()

        # function to load form MDPA to MODEL_PART
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO( file_path + file_name ).ReadModelPart(self.building_hull_model_part)

        self.HasBuildingHull = True


    def ComputeDistanceFieldFromHull( self, invert_distance_field = False ):

	    aux_model_part_name = "AuxModelPart"
	    current_model = complete_domain_model_part.GetModel()

	    if current_model.HasModelPart(aux_model_part_name):
		    # clear the existing model part (to be sure)
		    aux_model_part = current_model.GetModelPart( aux_model_part_name )
		    aux_model_part.Elements.clear()
		    aux_model_part.Conditions.clear()
		    aux_model_part.Nodes.clear()

	    else:
		    # create the model part from scratch
		    aux_model_part = current_model.CreateModelPart(aux_model_part_name)
		    aux_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
		    aux_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
		    aux_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)

	    # populating the empty model part (copy, do not just assign!)
	    prop = self.ModelPart.Properties[0]
	    for node in self.ModelPart.Nodes:
		    n = aux_model_part.CreateNewNode( node.Id, node.X, node.Y, node.Z )
	    for elem in self.ModelPart.Elements:
		    nodes = elem.GetNodes()
		    e = aux_model_part.CreateNewElement("Element3D4N", elem.Id,  [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], prop)

	    KratosMultiphysics.CalculateDistanceToSkinProcess3D(aux_model_part, self.building_hull_model_part ).Execute()
	    _test_output_dist( aux_model_part, "Terrain_with_ZeroLevel_1" )

	    # inversion of the distance field
	    if invert_distance_field:
		    for node in aux_model_part.Nodes:
			    node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, -node.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) )

	    # getting rid of the +/- inf values around the zero level
	    for node in aux_model_part.Nodes:
		    if ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) > 1000.0 ):
			    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.0 )
		    elif ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) < -1000.0 ):
			    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -1.0 )

	    variational_distance_process = _set_variational_distance_process_serial( aux_model_part, "DistanceFromSkin" )
	    variational_distance_process.Execute()

	    for node in aux_model_part.Nodes:
		    source = node
		    destination = complete_domain_model_part.GetNodes()[source.Id]
		    destination.SetSolutionStepValue( KratosMultiphysics.DISTANCE, source.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) + size_reduction )


### --- auxiliary functions --- ### -------------------------------------

    def _generate_building_model_part( self ):

        current_model = self.ModelPart.GetModel()
        self.building_hull_model_part = current_model.CreateModelPart( "BuildingModelPart" )
        self.building_hull_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.building_hull_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)

    def _test_output_dist( self, model_part, name ):

        from gid_output_process import GiDOutputProcess
        gid_output = GiDOutputProcess(model_part,
                                    name,
                                    KratosMultiphysics.Parameters("""
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


    def _set_variational_distance_process_serial(self, complete_model, aux_name):
        # Construct the variational distance calculation process
        serial_settings = KratosMultiphysics.Parameters("""
            {
                "linear_solver_settings"   : {
                    "solver_type" : "AMGCL"
                }
            }
        """)
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(serial_settings["linear_solver_settings"])
        maximum_iterations = 5
        if complete_model.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                complete_model,
                linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.CALCULATE_EXACT_DISTANCES_TO_PLANE,
        		aux_name )
        else:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                complete_model,
                linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.CALCULATE_EXACT_DISTANCES_TO_PLANE,
        		aux_name )
        return variational_distance_process