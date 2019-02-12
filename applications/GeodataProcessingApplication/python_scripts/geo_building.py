import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as KratosMesh
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo

from geo_processor import GeoProcessor
import os

class GeoBuilding( GeoProcessor ):

    def __init__( self ):
        super(GeoBuilding, self).__init__()

        self.HasBuildingHull = False
        self.HasDistanceField = False



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

        self.HasDistanceField = True


    def RefineMeshNearBuilding( self, single_parameter ):

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Model part has to be set, first.")
            return

        if not self.HasDistanceField:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Refinement around the building not possible. Please compute the distance field from the hull.")
            return

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



    def SubtractBuilding( self, min_size, max_size, hausdorff_value ):

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Model part has to be set, first.")
            return

        if not self.HasDistanceField:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Building subtraction not possible. Please compute the distance field from the hull.")
            return

        ### moving procedure to another model part
        current_model = self.ModelPart.GetModel()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        properties_0 = main_model_part.GetProperties()[0]
        properties_1 = main_model_part.GetProperties()[1]
        main_model_part.CreateSubModelPart("AuxSubModelPart")

        # Copying the content of the model
        for node in self.ModelPart.Nodes:
            n = main_model_part.CreateNewNode( node.Id, node.X, node.Y, node.Z )
            main_model_part.GetSubModelPart("AuxSubModelPart").AddNode( n, 0 )
        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()
            e = main_model_part.CreateNewElement("Element3D4N", elem.Id,  [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], properties_1)
            main_model_part.GetSubModelPart("AuxSubModelPart").AddElement( e, 0 )
        # coping (only!) the distance field
        for node in self.ModelPart.Nodes:
            receiver = main_model_part.GetNode( node.Id )
            receiver.SetSolutionStepValue( KratosMultiphysics.DISTANCE, node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) )

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
        find_nodal_h.Execute()

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in main_model_part.Nodes:
        	node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
        	{
        		"minimal_size"                         : 0.05,
        		"enforce_current"                      : false,
        		"anisotropy_remeshing"                 : true,
        		"anisotropy_parameters": {
        			"hmin_over_hmax_anisotropic_ratio"      : 0.1,
        			"boundary_layer_max_distance"           : 0.2,
        			"interpolation"                         : "Linear" }
        	}
        	""")
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(main_model_part, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        # The Hausdorff parameter is an important parameter to decrease the mesh size
        remesh_param = KratosMultiphysics.Parameters("""{
            "advanced_parameters": {
                "deactivate_detect_angle": false,
                "force_gradation_value": true,
                "force_hausdorff_value": true,
                "gradation_value": 1.3,
                "hausdorff_value": """ + str( hausdorff_value ) + """,
                "no_insert_mesh": false,
                "no_move_mesh": false,
                "no_surf_mesh": false,
                "no_swap_mesh": false,
                "normal_regularization_mesh": false
            },
            "buffer_size": 0,
            "debug_result_mesh": false,
            "discretization_type": "IsoSurface",
            "echo_level": 3,
            "extrapolate_contour_values": true,
            "filename": "out",
            "force_sizes": {
                "force_max": true,
                "force_min": true,
                "maximal_size": """ + str(max_size) + """,
                "minimal_size": """ + str(min_size) + """
            },
            "framework": "Eulerian",
            "initialize_entities": true,
            "internal_variables_parameters": {
                "allocation_size": 1000,
                "bucket_size": 4,
                "internal_variable_interpolation_list": [],
                "interpolation_type": "LST",
                "search_factor": 2
            },
            "interpolate_non_historical": true,
            "isosurface_parameters": {
                "isosurface_variable": "DISTANCE",
                "nonhistorical_variable": false,
                "remove_regions": true
            },
            "max_number_of_searchs": 1000,
            "remesh_at_non_linear_iteration": false,
            "save_external_files": false,
            "save_mdpa_file": false,
            "search_parameters": {
                "allocation_size": 1000,
                "bucket_size": 4,
                "search_factor": 2.0
            },
            "step_data_size": 0,
            "surface_elements": false
        }""")
        MmgProcess = KratosMesh.MmgProcess3D(main_model_part, remesh_param)
        MmgProcess.Execute()

        for node in main_model_part.Nodes:
            node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, 1.0 )
        for cond in main_model_part.Conditions:
            nodes = cond.GetNodes()
            for node in nodes:
                if ( node in main_model_part.Nodes ):
                    node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, -1.0e-7 )

        self.ModelPart.Nodes.clear()
        self.ModelPart.Conditions.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.CreateSubModelPart("Parts_Fluid")
        self.ModelPart.CreateSubModelPart("Complete_Boundary")

        # Copying the content of the model
        for node in main_model_part.Nodes:
            n = self.ModelPart.CreateNewNode( node.Id, node.X, node.Y, node.Z )
            self.ModelPart.GetSubModelPart("Parts_Fluid").AddNode( n, 0 )
        for elem in main_model_part.Elements:
            nodes = elem.GetNodes()
            e = self.ModelPart.CreateNewElement("Element3D4N", elem.Id,  [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], properties_1)
            self.ModelPart.GetSubModelPart("Parts_Fluid").AddElement( e, 0 )
        for cond in main_model_part.Conditions:
            nodes = cond.GetNodes()
            c = self.ModelPart.CreateNewCondition("Condition3D3N", cond.Id,  [nodes[0].Id, nodes[1].Id, nodes[2].Id], properties_0)
            self.ModelPart.GetSubModelPart("Complete_Boundary").AddCondition( c, 0 )
            for node in c.GetNodes():
                self.ModelPart.GetSubModelPart("Complete_Boundary").AddNode( node, 0 )



### --- auxiliary functions --- ### -------------------------------------

    def _generate_building_model_part( self ):

        current_model = self.ModelPart.GetModel()
        self.building_hull_model_part = current_model.CreateModelPart( "BuildingModelPart" )
        self.building_hull_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.building_hull_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)


    def _set_variational_distance_process_serial(self, complete_model, aux_name):
        # Construct the variational distance calculation process
        serial_settings = KratosMultiphysics.Parameters("""
            {
                "linear_solver_settings"   : {
                    "solver_type" : "amgcl"
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