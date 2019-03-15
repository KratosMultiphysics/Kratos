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

        self._generate_building_model_part("ModelPart"+file_name)

        # function to load form MDPA to MODEL_PART
        file_path = os.path.dirname(os.path.realpath(__file__))
        KratosMultiphysics.ModelPartIO( file_path + file_name ).ReadModelPart(self.building_hull_model_part)

        self.HasBuildingHull = True


    def ShiftBuildingHull( self, dx, dy, dz ):

        if not self.HasModelPart:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Model part has to be set, first.")
            return

        if not self.HasBuildingHull:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Function SiftBuildingHull requires to import a building hull, first.")
            return

        for node in self.building_hull_model_part.Nodes:
            node.X += dx
            node.Y += dy
            node.Z += dz


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

        # check the distance field
        pos = 0; neg = 0
        for node in aux_model_part.Nodes:
            if node.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) > 0.0:
                pos = 1
            if node.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) < 0.0:
                neg = 1

        if ( pos == 0 or neg == 0):
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "The distance field from the building hull does not have its zero-level inside the domain.")
            return False

        # inversion of the distance field
        if invert_distance_field:
            for node in aux_model_part.Nodes:
        	    node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, -node.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) )

        # getting rid of the +/- inf values around the zero level
        for node in aux_model_part.Nodes:
            if ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) > 1000.0 ):
        	    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 20.0 )
            elif ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) < -1000.0 ):
        	    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -20.0 )

        if (size_reduction != 0.0):
            variational_distance_process = self._set_variational_distance_process_serial( aux_model_part, "DistanceFromSkin" )
            variational_distance_process.Execute()

        for node in aux_model_part.Nodes:
            source = node
            destination = self.ModelPart.GetNodes()[source.Id]
            destination.SetSolutionStepValue( KratosMultiphysics.DISTANCE, source.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) + size_reduction )

        self.HasDistanceField = True
        return True


    def AddDistanceFieldFromHull( self, invert_distance_field = False, size_reduction = 0.0 ):

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

        # check the distance field
        pos = 0; neg = 0
        for node in aux_model_part.Nodes:
            if node.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) > 0.0:
                pos = 1
            if node.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) < 0.0:
                neg = 1

        if ( pos == 0 or neg == 0):
            return

        # inversion of the distance field
        if invert_distance_field:
            for node in aux_model_part.Nodes:
        	    node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, -node.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) )

        # getting rid of the +/- inf values around the zero level
        for node in aux_model_part.Nodes:
            if ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) > 1000.0 ):
        	    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 100.0 )
            elif ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) < -1000.0 ):
        	    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -100.0 )

        # full distance field only needed if shifts are necessary
        if (size_reduction != 0):
            variational_distance_process = self._set_variational_distance_process_serial( aux_model_part, "DistanceFromSkin" )
            variational_distance_process.Execute()

        for node in aux_model_part.Nodes:
            source = node
            destination = self.ModelPart.GetNodes()[source.Id]
            # finding the minimum
            distance = min( source.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) + size_reduction, destination.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) )
            destination.SetSolutionStepValue( KratosMultiphysics.DISTANCE, distance )

        self.HasDistanceField = True


    def SetInitialDistanceField( self, distance_value = 1.0 ):

        for node in self.ModelPart.Nodes:
            node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, distance_value )

        self.HasDistanceField = True


    def ShiftGlobalBuildingDistanceField( self, distance_shift_value = 0.0 ):

        if not self.HasDistanceField:
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "Refinement around the building not possible. Please compute the distance field from the hull.")
            return

        for node in self.ModelPart.Nodes:
            distance_value = node.GetSolutionStepValue( KratosMultiphysics.DISTANCE ) + distance_shift_value
            node.SetSolutionStepValue( KratosMultiphysics.DISTANCE, distance_value )


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
        		"enforce_current"                      : true,
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
        if current_model.HasModelPart( "MainModelPart" ):
            current_model.DeleteModelPart( "MainModelPart" )

        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        properties_0 = main_model_part.GetProperties()[0]
        properties_1 = main_model_part.GetProperties()[1]
        main_model_part.CreateSubModelPart("AuxSubModelPart")

        # # Copying the content of the model --- VERY SLOW --- NOW IN C++!
        # for node in self.ModelPart.Nodes:
        #     n = main_model_part.CreateNewNode( node.Id, node.X, node.Y, node.Z )
        #     main_model_part.GetSubModelPart("AuxSubModelPart").AddNode( n, 0 )
        # for elem in self.ModelPart.Elements:
        #     nodes = elem.GetNodes()
        #     e = main_model_part.CreateNewElement("Element3D4N", elem.Id,  [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], properties_1)
        #     main_model_part.GetSubModelPart("AuxSubModelPart").AddElement( e, 0 )
        # # coping (only!) the distance field
        # for node in self.ModelPart.Nodes:
        #     receiver = main_model_part.GetNode( node.Id )
        #     receiver.SetSolutionStepValue( KratosMultiphysics.DISTANCE, node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) )

        main_model_part = KratosGeo.CleaningUtilities( self.ModelPart ).HardCopyBeforeSurfaceDiscretization( self.ModelPart, main_model_part )

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
        		"enforce_current"                      : true,
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
                "gradation_value": 1.2,
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
            "interpolate_non_historical": false,
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
        self.ModelPart.RemoveSubModelPart("Parts_Fluid")
        self.ModelPart.CreateSubModelPart("Parts_Fluid")
        self.ModelPart.RemoveSubModelPart("Complete_Boundary")
        self.ModelPart.CreateSubModelPart("Complete_Boundary")

        self.ModelPart = KratosGeo.CleaningUtilities( self.ModelPart ).HardCopyAfterSurfaceDiscretization( main_model_part, self.ModelPart )

        # # Copying the content of the model --- VERY SLOW!
        # for node in main_model_part.Nodes:
        #     n = self.ModelPart.CreateNewNode( node.Id, node.X, node.Y, node.Z )
        #     self.ModelPart.GetSubModelPart("Parts_Fluid").AddNode( n, 0 )
        #     print( "Node " + str(node.Id) )
        # for elem in main_model_part.Elements:
        #     nodes = elem.GetNodes()
        #     e = self.ModelPart.CreateNewElement("Element3D4N", elem.Id,  [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], properties_1)
        #     self.ModelPart.GetSubModelPart("Parts_Fluid").AddElement( e, 0 )
        #     print( "Element " + str(elem.Id) )
        # for cond in main_model_part.Conditions:
        #     nodes = cond.GetNodes()
        #     c = self.ModelPart.CreateNewCondition("Condition3D3N", cond.Id,  [nodes[0].Id, nodes[1].Id, nodes[2].Id], properties_0)
        #     self.ModelPart.GetSubModelPart("Complete_Boundary").AddCondition( c, 0 )
        #     for node in c.GetNodes():
        #         self.ModelPart.GetSubModelPart("Complete_Boundary").AddNode( node, 0 )
        #     print( "Condition " + str(cond.Id) )


    def ShiftBuildingOnTerrain(self, buildings_model_part, terrain_model_part):
        # function to shift Buildings on terrain
        import numpy as np      # TODO: REMOVE NUMPY

        N = KratosMultiphysics.Vector(3)
        coords =  KratosMultiphysics.Array3()

        locate_on_background = KratosMultiphysics.BinBasedFastPointLocator2D(terrain_model_part)
        locate_on_background.UpdateSearchDatabase()

        already_moved = []            # list with all nodes already moved
        height = []

        self.ModelPart = buildings_model_part
        # KratosMultiphysics.ModelPartIO("data/building_model_part_ORG", KratosMultiphysics.IO.WRITE).WriteModelPart(self.ModelPart)
        ID_vertex = self.ModelPart.NumberOfNodes() + 1

        for num_building in range(self.ModelPart.NumberOfSubModelParts()):
            current_sub_model = self.ModelPart.GetSubModelPart("Building_{}".format(num_building+1))
            
            displacement = []        # the vector where we save the Z coordinate to evaluate the minimum
            for node_building in current_sub_model.Nodes:
                # take only nodes with z = 0.0 which are the nodes at the base of the Building
                if node_building.Z != 0.0:
                    continue

                # fill coords array
                coords[0] = node_building.X
                coords[1] = node_building.Y
                coords[2] = node_building.Z
                
                # "found" is a boolean variable (True if "node" is inside the mesh element)
                # "pelem" is a pointer to element that contain the node of the Building
                [found, N, pelem] = locate_on_background.FindPointOnMesh(coords)        # here we have the terrain element (but it is on xy plane)

                if not isinstance(pelem, KratosMultiphysics.Element):            # if "pelem" is not a "Kratos.Element", we go to the next node
                    continue
                
                # calculation of the intersection point between Building and terrain
                Norm = self._normal_triangle(pelem)

                # define plane
                planeNormal = np.array(Norm)
                planePoint = np.array(pelem.GetNode(0))

                # define ray
                rayDirection = np.array([0, 0, 1])
                rayPoint = np.array(coords)

                ndotu = planeNormal.dot(rayDirection)

                w = rayPoint - planePoint
                si = -planeNormal.dot(w) / ndotu
                Psi = w + si *rayDirection + planePoint
                
                displacement.append(Psi[2])                # append just coordinate Z
            
            # check if "displacement" is empty
            if displacement:
                height.append(min(displacement))                # quantity to move the n-th Building
            else:
                height.append(0)
            
            # check on the nodes shared by multiple elements
            node_mod = {}    # dictionary where the key is the old node Id and the "value" is the new node Id
            for node in current_sub_model.Nodes:
                if node.Id in already_moved:                    # if it is one of the already moved nodes
                    self.ModelPart.CreateNewNode(ID_vertex, node.X, node.Y, (node.Z + height[num_building]))    # a new node with different Id and Z coordinate moved
                    node_mod[node.Id] = ID_vertex                # fill the dictionary with the node to be removed (the key) and the node to be added (the value) in this current_sub_model
                    node.Id = ID_vertex                            # update the node of the element
                    ID_vertex += 1
                else:
                    already_moved.append(node.Id)                # fill "already_moved" with current Id node
                    node.Z = node.Z + height[num_building]            # move current node of the quantity in list "height"
            
            for old_node, new_node in node_mod.items():
                current_sub_model.AddNodes([new_node])            # add the new nodes in the current sub model part
                current_sub_model.RemoveNode(old_node)            # remove the old nodes that are now replaced

        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,True)

        for sub_model in self.ModelPart.SubModelParts:
            for node in sub_model.Nodes:
                node.Set(KratosMultiphysics.TO_ERASE,False)

        # to erase unused nodes
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

        # KratosMultiphysics.ModelPartIO("data/building_model_part_MOD", KratosMultiphysics.IO.WRITE).WriteModelPart(self.ModelPart)


    # FUNCTION UNDER CONSTRUCTION
    def DeleteBuildingsUnderValue(self, z_value):
        # this function delete the buildings that are under a certain Z value

        del_buildings = []
        for sub_model in self.ModelPart.SubModelParts:
            for elem in sub_model.Elements:
                for node in elem.GetNodes():
                    if (node.Z <= z_value):
                        del_buildings.append(sub_model.Name)
                        # print("save the information to will delete the element")
        
        del_buildings = list(set(del_buildings))
        
        for sub_model_to_del in del_buildings:
            self.ModelPart.RemoveSubModelPart(sub_model_to_del)


### --- auxiliary functions --- ### -------------------------------------

    def _generate_building_model_part( self, name ):

        current_model = self.ModelPart.GetModel()

        if current_model.HasModelPart( name ):
            current_model.DeleteModelPart( name )

        self.building_hull_model_part = current_model.CreateModelPart( name )
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
        maximum_iterations = 3
        if complete_model.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess2D(
                complete_model,
                linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess2D.NOT_CALCULATE_EXACT_DISTANCES_TO_PLANE,
        		aux_name )
        else:
            variational_distance_process = KratosMultiphysics.VariationalDistanceCalculationProcess3D(
                complete_model,
                linear_solver,
                maximum_iterations,
                KratosMultiphysics.VariationalDistanceCalculationProcess3D.NOT_CALCULATE_EXACT_DISTANCES_TO_PLANE,
        		aux_name )

        return variational_distance_process


    def _normal_triangle(self, elem):
        P1 = elem.GetNode(0)
        P2 = elem.GetNode(1)
        P3 = elem.GetNode(2)

        Nx = (P2.Y-P1.Y)*(P3.Z-P1.Z) - (P3.Y-P1.Y)*(P2.Z-P1.Z)
        Ny = (P2.Z-P1.Z)*(P3.X-P1.X) - (P2.X-P1.X)*(P3.Z-P1.Z)
        Nz = (P2.X-P1.X)*(P3.Y-P1.Y) - (P3.X-P1.X)*(P2.Y-P1.Y)

        return [Nx, Ny, Nz]
