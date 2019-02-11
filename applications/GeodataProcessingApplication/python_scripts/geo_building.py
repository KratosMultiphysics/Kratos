import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMesh

from geo_processor import GeoProcessor
import os

class GeoBuilding( GeoProcessor ):

    def __init__( self ):
        super(GeoBuilding, self).__init__()

        self.HasBuildingHull = False



    def ImportBuildingHullSTL( self, file_name ):

        print("To be done...")



    def ImportBuildingHullMDPA( self, file_name ):

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Model part has to be set, first.")
            return

        self._generate_building_model_part()

        # function to load form MDPA to MODEL_PART
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO( file_path + file_name ).ReadModelPart(self.building_hull_model_part)

        self.HasBuildingHull = True



    def ComputeDistanceFieldFromHull( self, invert_distance_field = False, size_reduction = 0.0 ):

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Model part has to be set, first.")
            return

        if not self.HasBuildingHull:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Function ComputeDistanceFieldFromHull requires to import a building hull, first.")
            return

        aux_model_part_name = "AuxModelPart"
        current_model = self.ModelPart.GetModel()

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

        variational_distance_process = self._set_variational_distance_process_serial( aux_model_part, "DistanceFromSkin" )
        variational_distance_process.Execute()

        for node in aux_model_part.Nodes:
            source = node
            destination = self.ModelPart.GetNodes()[source.Id]
            destination.SetSolutionStepValue( KratosMultiphysics.DISTANCE, source.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) + size_reduction )


    def RefineMeshNearBuilding( self, single_parameter ):

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess( self.ModelPart )
        find_nodal_h.Execute()

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.ModelPart.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.ModelPart, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero (or unit) the metric
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 1.0; ZeroVector[1] = 1.0; ZeroVector[2] = 1.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in self.ModelPart.Nodes:
        	node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        min_size = single_parameter
        max_dist = 1.25 * single_parameter
        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
        	{
        		"minimal_size"                         : """ + str(min_size) + """,
        		"enforce_current"                      : false,
        		"anisotropy_remeshing"                 : true,
        		"anisotropy_parameters": {
        			"hmin_over_hmax_anisotropic_ratio"      : 0.9,
        			"boundary_layer_max_distance"           : """ + str(max_dist) + """,
        			"interpolation"                         : "Linear" }
        	}
        	""")
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(self.ModelPart, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        # We create the remeshing process
        remesh_param = KratosMultiphysics.Parameters("""{ }""")
        MmgProcess = KratosMesh.MmgProcess3D(self.ModelPart, remesh_param)
        MmgProcess.Execute()

        print( "next step")



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