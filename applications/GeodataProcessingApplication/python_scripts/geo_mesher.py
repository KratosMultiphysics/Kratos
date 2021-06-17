import KratosMultiphysics
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.MeshingApplication as KratosMesh
if KratosMultiphysics.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI
    print("[DEBUG][geo_mesher] import KratosMultiphysics.mpi")

from geo_processor import GeoProcessor
import triangle
import numpy as np
import math
from meshpy.tet import MeshInfo, build

class GeoMesher( GeoProcessor ):

    def __init__( self ):
        super(GeoMesher, self).__init__()

        # self.HasModelPart = False     # [NG] reduntant!!! we define it in the GeoProcessor class
        self.HasExtrusionHeight = False

        # useful variables for volumetric mesh and for division into sectors
        self.x_center = float()     # X coordinate of the center of the domain
        self.y_center = float()     # Y coordinate of the center of the domain
        self.r_boundary = float()   # radius of cylindrical domain
        self.r_ground = float()     # radius of ground (smoothing between r_ground and r_boundary)
        self.r_buildings = float()  # radius of the portion on which the buildings will be placed

    """ TEST FUNCTIONS """
    def Mesh_2D(self):
        all_points = []     # list with all coordinates of nodes
        coord_2D = []
        for node in self.ModelPart.Nodes:
            all_points.append((node.X, node.Y, node.Z))
            coord_2D.append([node.X, node.Y])

        num_id_bot = len(all_points)    # number of id up to this moment

        ### TRIANGLE ###
        A = dict(vertices = coord_2D)
        triangle_dict = triangle.triangulate(A)

        all_facets = []		# list with all node ids of faces
        all_markers = []	# list of integers [1 = bottom, 2 = topper and 3 = lateral]
        list_id_bottom = []

        # we fill list for "bottom" faces
        for f in triangle_dict["triangles"]:
            all_facets.append([f[0], f[1], f[2]])
            all_markers.append(1)

        list_id_top = []

        # we fill list for "top" faces
        for f in triangle_dict["triangles"]:
            all_facets.append([f[0], f[1], f[2]])
            all_markers.append(2)

        # we fill all_facets list with "lateral" faces
        for i in range(bottom_elem):
            if (i == bottom_elem-1):
                # here we have the last value
                all_facets.append([list_id_circle_bottom[i], list_id_circle_bottom[0], list_id_circle_top[0], list_id_circle_top[i]])
                all_markers.append(3)
                break
            all_facets.append([list_id_circle_bottom[i], list_id_circle_bottom[i+1], list_id_circle_top[i+1], list_id_circle_top[i]])
            all_markers.append(3)

        pass


    ### --- functions to determine an extrusion height --- ###

    def ComputeExtrusionHeight( self, radius, height, free_board, iterations ):

        tool = KratosGeo.ExtrusionHeightUtilities( self.ModelPart )
        tool.SetExtrusionHeight( radius, height, free_board )
        tool.SmoothExtrusionHeight( radius/5.0, iterations, 0.8*free_board )
        self.HasExtrusionHeight = True

    def SetUniformExtrusionHeight( self, height ):

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable( KratosGeo.EXTRUSION_HEIGHT, height, self.ModelPart.Nodes )
        self.HasExtrusionHeight = True

    ### --- functions to create a 3D domain mesh --- ###

    def MeshConcaveHullWithTerrainPoints( self, average_point_dist ):

        if ( not self.HasExtrusionHeight ):
            KratosMultiphysics.Logger.PrintWarning("GeoMesher", "Extrusion height must be determined first.")
            return

        ### reading points from model part
        X = []; Y = []; Z = []      # lists are generated to use predefined operations
        for node in self.ModelPart.Nodes:
            X.append(node.X)
            Y.append(node.Y)
            Z.append(node.Z)

        ### preparing nodes to be handed over to triangle
        pts = np.vstack((X, Y)).T
        A = dict(vertices = pts)
        triangle_dict = triangle.triangulate(A)		# triangle_dict {vertices: [...], triangles: [...], vertex_markers: [...]}

        all_bottom_points = []			# vector with ALL coordinates of nodes on the terrain
        all_bottom_facets = []			# vector with ALL node ids of faces
        all_top_points = []				# vector with ALL coordinates of nodes on the terrain

        for node in self.ModelPart.Nodes:
			# nodes should now stay were the are (!!!)
            all_bottom_points.append( (node.X, node.Y, node.Z) )
            all_top_points.append( (node.X, node.Y, node.GetValue( KratosGeo.EXTRUSION_HEIGHT )) )

        for face in triangle_dict["triangles"]:
			# facets are further checked
            all_bottom_facets.append([face[0], face[1], face[2]])

        maximal_initial_id_terrain_point = len( all_bottom_points )		# id up to this moment for the actual imported terrain
        maximal_initial_id_terrain_facet = len( all_bottom_facets )		# id up to this moment for the actual imported terrain

        ### reading back into the model part from the triangle library
        properties = self.ModelPart.Properties[0]
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()
        self.ModelPart.Nodes.clear()

        for i in range( 0, maximal_initial_id_terrain_point ):
            # node id is shifted up by 1
            node = self.ModelPart.CreateNewNode( (i+1), (all_bottom_points[i])[0], (all_bottom_points[i])[1] , 0.0)

        for i in range( 0, maximal_initial_id_terrain_facet ):
            # node id is shifted up by 1
            nodes_for_element = [ (all_bottom_facets[i])[0]+1, (all_bottom_facets[i])[1]+1, (all_bottom_facets[i])[2]+1 ]
            # facet id is shifted up by 1
            elem = self.ModelPart.CreateNewElement("Element2D3N", (i+1), nodes_for_element, properties )

        ### generation of the concave hull
        Mesher = KratosMesh.TriGenPFEMModeler()
        for node in self.ModelPart.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.NODAL_H, 0, 0.1)
            node.SetSolutionStepValue(KratosMultiphysics.IS_FLUID, 0, 1)

        node_erase_process = KratosMultiphysics.NodeEraseProcess(self.ModelPart)
        neigh_finder = KratosMultiphysics.FindNodalNeighboursProcess(self.ModelPart, 9, 18)
        neigh_finder.Execute()
        Mesher.ReGenerateMesh("Element2D3N", "Condition2D", self.ModelPart, node_erase_process, False, False, 1.4*average_point_dist, 0.5)
        
	    ### read back from "element + condition" to replace the old "facets"  ---
        all_bottom_facets = []
        all_bottom_markers = []				# BOTTOM = marker 1
        all_top_facets = []
        all_top_markers = []				# TOP = marker 2
        all_lateral_facets = []
        all_lateral_markers = []				# LATERAL = marker 3

        ### preparing a closed volume
        for elem in self.ModelPart.Elements:
            # the shift must be performed in reverse
            node1 = (elem.GetNodes()[0]).Id - 1
            node2 = (elem.GetNodes()[1]).Id - 1
            node3 = (elem.GetNodes()[2]).Id - 1
            # facets at the - bottom - (Marker 1)
            all_bottom_facets.append( [ node1, node2, node3 ] )
            all_bottom_markers.append(1)
            # facets at the - top - (Marker 2)
            all_top_facets.append( [ 	node1 + maximal_initial_id_terrain_point,
            						    node2 + maximal_initial_id_terrain_point,
            						    node3 + maximal_initial_id_terrain_point ] )
            all_top_markers.append(2)

        # creating the - lateral - limitting facets from the conditions in the model_part
        for cond in self.ModelPart.Conditions:
		# the shift must be performed in reverse
            node1 = (cond.GetNodes()[0]).Id - 1
            node2 = (cond.GetNodes()[1]).Id - 1
            node3 = node1 + maximal_initial_id_terrain_point
            node4 = node2 + maximal_initial_id_terrain_point
            all_lateral_facets.append([ node1, node2, node4, node3 ])
            all_lateral_markers.append(3)

        ### putting all components together
        # putting the node lists together
        all_points = []
        all_points.extend( all_bottom_points )
        all_points.extend( all_top_points )
        # putting facet lists togehter (same order as markers)
        all_facets = []
        all_facets.extend( all_bottom_facets )
        all_facets.extend( all_top_facets )
        all_facets.extend( all_lateral_facets )
        print( len(all_facets) )
        # putting the markers together (same order as facets)
        all_markers = []
        all_markers.extend( all_bottom_markers )
        all_markers.extend( all_top_markers )
        all_markers.extend( all_lateral_markers )
        print( len(all_markers) )
        self.nodes = all_points
        self.facets = all_facets

        ### using Meshpy for the creation of an initial mesh https://mathema.tician.de/software/meshpy/
        points = np.array( all_points )
        facets = np.array( all_facets )

        mesh_info = MeshInfo()
        mesh_info.set_points(points)
        mesh_info.set_facets(facets, markers=all_markers )

        # build the mesh
        mesh = build(mesh_info)
        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()

        lateral_model_part = self.ModelPart.CreateSubModelPart("LateralModelPart")
        bottom_model_part = self.ModelPart.CreateSubModelPart("BottomModelPart")
        top_model_part = self.ModelPart.CreateSubModelPart("TopModelPart")
        properties = self.ModelPart.Properties[0]

        bottom_cond = []; bottom_points = []
        top_cond = []; top_points = []
        lateral_cond = []; lateral_points = []

        # AWARE: Shift by 1 in index!!!
        ### Nodes
        for i in range( 0, len( mesh.points ) ):
        	# node id is shifted up by 1
        	coords = mesh.points[i]
        	node = self.ModelPart.CreateNewNode( (i+1), coords[0], coords[1], coords[2] )

        ### Conditions and Nodes in SubModelParts
        for j in range( 0, len( mesh.faces ) ):
            points = mesh.faces[j]
            marker = mesh.face_markers[j]
            if ( marker == 1 ):
                # cond = self.ModelPart.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                cond = self.ModelPart.CreateNewCondition("WallCondition3D3N", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                bottom_cond.append( j+1 )
            elif ( marker == 2 ):
                # cond = self.ModelPart.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                cond = self.ModelPart.CreateNewCondition("WallCondition3D3N", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                top_cond.append( j+1 )
            elif ( marker == 3 ):
                # cond = self.ModelPart.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                cond = self.ModelPart.CreateNewCondition("WallCondition3D3N", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                lateral_cond.append( j+1 )

        bottom_model_part.AddConditions( bottom_cond )
        top_model_part.AddConditions( top_cond )
        lateral_model_part.AddConditions( lateral_cond )

        for cond in bottom_model_part.Conditions:
        	for node in cond.GetNodes():
        		bottom_points.append( node.Id )
        for cond in top_model_part.Conditions:
        	for node in cond.GetNodes():
        		top_points.append( node.Id )
        for cond in lateral_model_part.Conditions:
        	for node in cond.GetNodes():
        		lateral_points.append( node.Id )

        bottom_model_part.AddNodes( bottom_points )
        top_model_part.AddNodes( top_points )
        lateral_model_part.AddNodes( lateral_points )

        ### Elements
        for k in range( 0, len( mesh.elements ) ):
        	points = mesh.elements[k]
        	elem = self.ModelPart.CreateNewElement("Element3D4N", (k+1), [ points[0]+1, points[1]+1, points[2]+1, points[3]+1 ], properties )


    def ComputeDistanceFieldFromGround( self ):

        if( self.ModelPart.HasSubModelPart("BottomModelPart") ):

            for node in self.ModelPart.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1.0)

            for node in self.ModelPart.GetSubModelPart("BottomModelPart").Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -1.0e-7)

        else:

            KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.ModelPart.Conditions, 3)
            # assigning the initial distances
            for node in self.ModelPart.Nodes:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, 1.0)

            for cond in self.ModelPart.Conditions:
                n = cond.GetValue(KratosMultiphysics.NORMAL)
                if( n[2] / math.sqrt(n[0]**2 + n[1]**2 + n[2]**2) < -0.0001 ):
                    for node in cond.GetNodes():
                        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -1.0e-7)


        # computing the distance field
        variational_distance_process = self._set_variational_distance_process_serial( self.ModelPart, "DistanceFromGround1" )
        variational_distance_process.Execute()


    def RefineMeshNearGround( self, single_parameter ):

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess( self.ModelPart )
        find_nodal_h.Execute()

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.ModelPart.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.ModelPart, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero the metric
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 1.0; ZeroVector[1] = 1.0; ZeroVector[2] = 1.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in self.ModelPart.Nodes:
        	node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        min_size = single_parameter
        max_dist = 1.5 * single_parameter
        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
        	{
        		"minimal_size"                         : """ + str(min_size) + """,
        		"enforce_current"                      : true,
        		"anisotropy_remeshing"                 : true,
        		"anisotropy_parameters": {
        			"hmin_over_hmax_anisotropic_ratio"      : 1.0,
        			"boundary_layer_max_distance"           : """ + str(max_dist) + """,
        			"interpolation"                         : "linear" }
        	}
        	""")
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(self.ModelPart, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        # We create the remeshing process
        remesh_param = KratosMultiphysics.Parameters("""{ }""")
        MmgProcess = KratosMesh.MmgProcess3D(self.ModelPart, remesh_param)
        MmgProcess.Execute()


    def RefineMesh_test(self, min_size, max_size, max_dist=1.0):
        # test function to refine the mesh

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess( self.ModelPart )
        find_nodal_h.Execute()

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.ModelPart.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.ModelPart, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero the metric
        # the metric tensor is a symmetric matrix; therefore we can provide just the upper triangular part
        # the convention MMG API is: (m11, m12, m13, m22, m23, m33)
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in self.ModelPart.Nodes:
            node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
            {
                "minimal_size"                         : """ + str(min_size) + """,
                "maximal_size"                         : """ + str(max_size) + """,
                "sizing_parameters":
                {
                    "reference_variable_name"          : "DISTANCE",
                    "boundary_layer_max_distance"      : """ + str(max_dist) + """,
                    "interpolation"                    : "linear"
                },
                "enforce_current"                      : true,
                "anisotropy_remeshing"                 : true,
                "anisotropy_parameters":
                {
                    "reference_variable_name"              : "DISTANCE",
                    "hmin_over_hmax_anisotropic_ratio"      : 1.0,
                    "boundary_layer_max_distance"           : """ + str(max_dist) + """,
                    "interpolation"                         : "linear"
        }
            }
            """)
        # level_set_param = KratosMultiphysics.Parameters("""
        #     {
        #         "minimal_size"                         : """ + str(min_size) + """,
        #         "maximal_size"                         : """ + str(max_size) + """,
        #         "sizing_parameters":
        #         {
        #             "reference_variable_name"          : "DISTANCE",
        #             "boundary_layer_max_distance"      : """ + str(max_dist) + """,
        #             "interpolation"                    : "linear"
        #         }
        # 	}
        # 	""")
        # "boundary_layer_max_distance"      : 1.0,
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(self.ModelPart, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        # We create the remeshing process
        remesh_param = KratosMultiphysics.Parameters("""{ }""")

        if KratosMultiphysics.IsDistributedRun():
            # MPI -> ParMMG
            # KratosMultiphysics.mpi.ParallelFillCommunicator(self.ModelPart).Execute()
            KratosMPI.ParallelFillCommunicator(self.ModelPart).Execute()
            pmmg_process = KratosMesh.ParMmgProcess3D(self.ModelPart, remesh_param)
            pmmg_process.Execute()
            print("[DEBUG][geo_mesher][RefineMesh_test] ParMMG DONE!")
        else:
            # serial -> MMG
            mmg_process = KratosMesh.MmgProcess3D(self.ModelPart, remesh_param)
            mmg_process.Execute()
            print("[DEBUG][geo_mesher][RefineMesh_test] MMG DONE!")

#######################################################################################################################
    def RefineMesh_mpi(self, min_size, max_size, max_dist=1.0):
        # test function to refine the mesh
        
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        size = comm.Get_size()

        find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess( self.ModelPart )
        find_nodal_h.Execute()

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, self.ModelPart.Nodes)
        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.ModelPart, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        # We set to zero the metric
        # the metric tensor is a symmetric matrix; therefore we can provide just the upper triangular part
        # the convention MMG API is: (m11, m12, m13, m22, m23, m33)
        ZeroVector = KratosMultiphysics.Vector(6)
        ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
        ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

        for node in self.ModelPart.Nodes:
            node.SetValue(KratosMesh.METRIC_TENSOR_3D, ZeroVector)

        # We define a metric using the ComputeLevelSetSolMetricProcess
        level_set_param = KratosMultiphysics.Parameters("""
            {
                "minimal_size"                         : """ + str(min_size) + """,
                "maximal_size"                         : """ + str(max_size) + """,
                "sizing_parameters":
                {
                    "reference_variable_name"          : "DISTANCE",
                    "boundary_layer_max_distance"      : """ + str(max_dist) + """,
                    "interpolation"                    : "linear"
                },
                "enforce_current"                      : true,
                "anisotropy_remeshing"                 : true,
                "anisotropy_parameters":
                {
                    "reference_variable_name"              : "DISTANCE",
                    "hmin_over_hmax_anisotropic_ratio"      : 1.0,
                    "boundary_layer_max_distance"           : """ + str(max_dist) + """,
                    "interpolation"                         : "linear"
                }
            }
            """)
        # level_set_param = KratosMultiphysics.Parameters("""
        #     {
        #         "minimal_size"                         : """ + str(min_size) + """,
        #         "maximal_size"                         : """ + str(max_size) + """,
        #         "sizing_parameters":
        #         {
        #             "reference_variable_name"          : "DISTANCE",
        #             "boundary_layer_max_distance"      : """ + str(max_dist) + """,
        #             "interpolation"                    : "linear"
        #         }
        # 	}
        # 	""")
        # "boundary_layer_max_distance"      : 1.0,
        metric_process = KratosMesh.ComputeLevelSetSolMetricProcess3D(self.ModelPart, KratosMultiphysics.DISTANCE_GRADIENT, level_set_param)
        metric_process.Execute()

        # We create the remeshing process
        remesh_param = KratosMultiphysics.Parameters("""{ }""")


        if KratosMultiphysics.IsDistributedRun():
            # MPI -> ParMMG
            # KratosMPI.ParallelFillCommunicator(self.ModelPart).Execute()
            pmmg_process = KratosMesh.ParMmgProcess3D(self.ModelPart, remesh_param)
            pmmg_process.Execute()
            print("[DEBUG][geo_mesher][RefineMesh_test] ParMMG DONE!")
        else:
            # serial -> MMG
            mmg_process = KratosMesh.MmgProcess3D(self.ModelPart, remesh_param)
            mmg_process.Execute()
            print("[DEBUG][geo_mesher][RefineMesh_test] MMG DONE!")
# END MPI
#############################################################################################################################################

    def MeshCircleWithTerrainPoints( self, h_value=0.0, circ_division=60, extract_center=False ):

        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()

        # sub model part to store the nodes of the bottom
        node_inner_bottom = self.ModelPart.CreateSubModelPart("NodeInnerBottom")

        x_min = float("inf");   x_max = float("-inf")   # we initialize the variables
        y_min = float("inf");   y_max = float("-inf")   # we initialize the variables
        for node in self.ModelPart.Nodes:
            # we add nodes in the sub model part
            node_inner_bottom.AddNode(node, 0)
            # we compute min and max X coordinate
            if (node.X < x_min): x_min = node.X
            elif (node.X > x_max): x_max = node.X
            # we compute min and max Y coordinate
            if (node.Y < y_min): y_min = node.Y
            elif (node.Y > y_max): y_max = node.Y

        # the centre of the circle
        x_center = (x_min + x_max)/2
        y_center = (y_min + y_max)/2

        # radius calculation considering the smaller side of the rectangular domain
        r_boundary = min((x_max-x_min), (y_max-y_min))/2    # radius of the domain
        r_ground = r_boundary * 2.0 / 3.0                   # 2/3 of the radius of the domain
        r_buildings = r_boundary / 3.0                      # 1/3 of the radius of the domain

        # buffering zone
        delta = r_boundary/1000         # evaluate this value!!!

        """ DELETING NODES PROCEDURE (nodes outside the r_boundary) """
        del_id = []             # list where there are the ids to be deleted because are outside of r_boundary
        z_min = float("inf");   z_max = float("-inf")   # we initialize the variables
        # for index in idx_to_smooth:
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE, False)    # we set all nodes as TO_ERASE=False
            
            x_curr = node.X             # current X coordinate
            y_curr = node.Y             # current Y coordinate
            dist = math.sqrt(((x_center-x_curr)**2)+((y_center-y_curr)**2))     # calculate the distance between current node and center of circle
            # we check that the nodes are inside the cylindrical domain
            if (dist > r_boundary-delta):
                del_id.append(node.Id)      # we will delete the nodes outside the domain
                continue
            else:
                # we compute min and max Z coordinate (we considering only nodes inside the domain)
                if (node.Z < z_min): z_min = node.Z
                elif (node.Z > z_max): z_max = node.Z

        for id in del_id:
            self.ModelPart.GetNode(id).Set(KratosMultiphysics.TO_ERASE, True)
        
        # we delete all nodes positioned outside the cylindrical domain
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

        # renumerate process [CHECK IT]
        node_Id = 1
        for node in self.ModelPart.Nodes:
            node.Id = node_Id
            node_Id += 1
        print("\n*** Renumbered Nodes DONE! ***\n")

        """ SMOOTHING PROCEDURE (bottom) """
        for node in self.ModelPart.Nodes:
            x_curr = node.X
            y_curr = node.Y
            dist = math.sqrt(((x_center-x_curr)**2)+((y_center-y_curr)**2))
            if dist > r_ground:
                Z_beta = (-(dist-r_ground) / (r_boundary-r_ground))+1       # we calculate beta
                self.ModelPart.GetNode(node.Id).Z = (node.Z - z_min) * Z_beta + z_min

        # a list with 360 degree (from 0 to 360)
        # theta = self._custom_range(0.0, 2*math.pi, math.pi/90)
        theta = self._custom_range(0.0, 2*math.pi, 2*math.pi/circ_division)  # circ_division is the number of division of the 2*pi

        # lists with nodes on circle boundary
        x_circle = []; y_circle = []
        node_Id = self.ModelPart.NumberOfNodes() + 1
        node_circle_bottom = self.ModelPart.CreateSubModelPart("NodeCircleBottom")  # sub model part in which we will store all the nodes along the circle boundary
        for th in theta:
            # TODO: WE CAN USE A VARIABLE INSTEAD OF A LIST
            x_circle.append(r_boundary * math.cos(th) + x_center)   # we move the node along X to consider that the circle have not the centre in (0.0, 0.0)
            y_circle.append(r_boundary * math.sin(th) + y_center)   # we move the node along Y to consider that the circle have not the centre in (0.0, 0.0)
            # we add the circle nodes in the ModelPart...
            node = self.ModelPart.CreateNewNode(node_Id, x_circle[-1], y_circle[-1], z_min)
            # ...and we add it also in the sub model part
            node_circle_bottom.AddNode(node, 0)
            node_Id += 1

        """ TRIANGLE (bottom) """
        # preparing nodes to be handed over to triangle
        coord_2D = []       # vertically list with X and Y coordinates (TO AVOID NUMPY VSTACK)
        all_points = []     # we need this list for TetGen process
        for node in self.ModelPart.Nodes:
            coord_2D.append([node.X, node.Y])   # we add coordinates X and Y of terrain nodes (bottom)
            all_points.append((node.X, node.Y, node.Z))

        vertices = dict(vertices = coord_2D)
        triangle_dict_bottom = triangle.triangulate(vertices)     # triangle_dict_bottom {vertices: [...], triangles: [...], vertex_markers: [...]}

        # TODO: CHECK IF WE CAN USE THE ELEMENTS IN MODEL PART INSTEAD OF THIS LIST
        all_facets = []		# list with all node ids of faces
        all_markers = []	# list of integers [1 = bottom, 2 = topper and 3 = lateral]

        # we fill the lists for "bottom"
        for face in triangle_dict_bottom["triangles"]:
            all_facets.append([face[0], face[1], face[2]])
            all_markers.append(1)
        
        # [VALUATE IF NECESSARY] only for a test
        elem_Id = 1
        for face in triangle_dict_bottom["triangles"]:
            self.ModelPart.CreateNewElement("Element2D3N", elem_Id, [face[0]+1, face[1]+1, face[2]+1], self.ModelPart.Properties[0])
            elem_Id += 1

        # domain height
        volume_height = z_max + h_value             # TODO: this value will be set by users

        # we copy the nodes along the circle boundary (with different Z coordinate)
        node_circle_top = self.ModelPart.CreateSubModelPart("NodeCircleTop")
        for node in self.ModelPart.GetSubModelPart("NodeCircleBottom").Nodes:
            node_top = self.ModelPart.CreateNewNode(node_Id, node.X, node.Y, volume_height)
            node_circle_top.AddNode(node_top, 0)
            node_Id += 1
            all_points.append((node.X, node.Y, volume_height))
        
        # [VALUATE IT]
        # we add a node in (x_center, y_center, volume_height) for a best triangular mesh
        node_top = self.ModelPart.CreateNewNode(node_Id, x_center, y_center, volume_height)
        node_circle_top.AddNode(node_top, 0)
        node_Id += 1
        all_points.append((x_center, y_center, volume_height))
        
        """ TRIANGLE (top) """
        # preparing nodes to be handed over to triangle
        coord_2D = []       # vertically list with X and Y coordinates (TO AVOID NUMPY VSTACK)
        for node in self.ModelPart.GetSubModelPart("NodeCircleTop").Nodes:
            coord_2D.append([node.X, node.Y])   # we add coordinates X and Y of top nodes

        vertices = dict(vertices = coord_2D)
        triangle_dict_top = triangle.triangulate(vertices)     # triangle_dict_top {vertices: [...], triangles: [...], vertex_markers: [...]}

        nun_node_bottom = self.ModelPart.GetSubModelPart("NodeCircleBottom").NumberOfNodes()
        nun_node_bottom += self.ModelPart.GetSubModelPart("NodeInnerBottom").NumberOfNodes()

        # we fill the lists for "top"
        for face in triangle_dict_top["triangles"]:
            all_facets.append([face[0]+nun_node_bottom+1, face[1]+nun_node_bottom+1, face[2]+nun_node_bottom+1])
            all_markers.append(2)
        
        # [VALUATE IF NECESSARY] only for a test
        for face in triangle_dict_top["triangles"]:
            self.ModelPart.CreateNewElement("Element2D3N", elem_Id, [face[0]+nun_node_bottom+1, face[1]+nun_node_bottom+1, face[2]+nun_node_bottom+1], self.ModelPart.Properties[0])
            elem_Id += 1

        # we fill the lists for "lateral"
        # TODO: USE DIRECTLY THE SUB MODEL PART WITHOUT THE LISTS
        node_bottom = []
        for node in self.ModelPart.GetSubModelPart("NodeCircleBottom").Nodes:
            node_bottom.append(node.Id)
        node_top = []
        for node in self.ModelPart.GetSubModelPart("NodeCircleTop").Nodes:
            node_top.append(node.Id)
        
        for i in range(len(node_bottom)):
            if (i == len(node_bottom)-1):
                # here we have the last value of the list
                all_facets.append([node_bottom[i], node_bottom[0], node_top[0], node_top[i]])
                all_markers.append(3)

                # [VALUATE IF NECESSARY] only for a test
                self.ModelPart.CreateNewElement("Element2D4N", elem_Id, [node_bottom[i], node_bottom[0], node_top[0], node_top[i]], self.ModelPart.Properties[0])
                elem_Id += 1
            else:
                all_facets.append([node_bottom[i], node_bottom[i+1], node_top[i+1], node_top[i]])
                all_markers.append(3)

                # [VALUATE IF NECESSARY] only for a test
                self.ModelPart.CreateNewElement("Element2D4N", elem_Id, [node_bottom[i], node_bottom[i+1], node_top[i+1], node_top[i]], self.ModelPart.Properties[0])
                elem_Id += 1


        """ TETGEN (VIA MESHPY) """
        ### using Meshpy for the creation of an initial mesh https://mathema.tician.de/software/meshpy/
        mesh_info = MeshInfo()
        mesh_info.set_points(all_points)
        mesh_info.set_facets(all_facets, markers=all_markers)

        # build the mesh
        mesh = build(mesh_info)

        # we clear the ModelPart
        self.ModelPart.RemoveSubModelPart("NodeCircleBottom")
        self.ModelPart.RemoveSubModelPart("NodeInnerBottom")
        self.ModelPart.RemoveSubModelPart("NodeCircleTop")

        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()

        lateral_model_part = self.ModelPart.CreateSubModelPart("LateralModelPart")
        bottom_model_part = self.ModelPart.CreateSubModelPart("BottomModelPart")
        top_model_part = self.ModelPart.CreateSubModelPart("TopModelPart")
        properties = self.ModelPart.Properties[0]

        bottom_cond = []; bottom_points = []
        top_cond = []; top_points = []
        lateral_cond = []; lateral_points = []

        # AWARE: Shift by 1 in index!!!
        ### Nodes
        for i in range(len(mesh.points)):
            # node id is shifted up by 1
            coords = mesh.points[i]
            self.ModelPart.CreateNewNode((i+1), coords[0], coords[1], coords[2])

        ### Conditions and Nodes in SubModelParts
        for j in range(len(mesh.faces)):
            points = mesh.faces[j]
            marker = mesh.face_markers[j]
            if (marker == 1):
                self.ModelPart.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties)
                bottom_cond.append(j+1)
            elif (marker == 2):
                self.ModelPart.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties)
                top_cond.append(j+1)
            elif ( marker == 3 ):
                self.ModelPart.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties)
                lateral_cond.append(j+1)

        bottom_model_part.AddConditions(bottom_cond)
        top_model_part.AddConditions(top_cond)
        lateral_model_part.AddConditions(lateral_cond)

        for cond in bottom_model_part.Conditions:
            for node in cond.GetNodes():
                bottom_points.append(node.Id)
        for cond in top_model_part.Conditions:
            for node in cond.GetNodes():
                top_points.append(node.Id)
        for cond in lateral_model_part.Conditions:
            for node in cond.GetNodes():
                lateral_points.append(node.Id)

        bottom_model_part.AddNodes(bottom_points)
        top_model_part.AddNodes(top_points)
        lateral_model_part.AddNodes(lateral_points)

        ### Elements
        for k in range(len(mesh.elements)):
            points = mesh.elements[k]
            self.ModelPart.CreateNewElement("Element3D4N", (k+1), [points[0]+1, points[1]+1, points[2]+1, points[3]+1], properties)

        print("\nModelPart\n", self.ModelPart)

        """ OPTIONAL VALUES """
        if extract_center:
            # we extract the center of the geometry and we add them in a sub model part
            center_elem = self.ModelPart.CreateSubModelPart("CenterElement")

            # elem_inside = {}        # dictionary with the information about the element that are in r_buildings
            # elem_id = 1

            elem_list = []  #########################################################################################################################

            for node in triangle_dict["triangles"]:
                for i in range(3):
                    x_curr = self.ModelPart.GetNode(node[i]+1).X
                    y_curr = self.ModelPart.GetNode(node[i]+1).Y
                    dist = math.sqrt(((x_center-x_curr)**2) + ((y_center-y_curr)**2))

                    if dist > r_buildings:
                        break   # if at least one node is external to r_buildings, we reject the entire element

                else:
                    # we get here if all nodes of the element are inside the r_buildings
                    node1 = self.ModelPart.GetNode(node[0]+1)
                    node2 = self.ModelPart.GetNode(node[1]+1)
                    node3 = self.ModelPart.GetNode(node[2]+1)

                    # elem_inside[elem_id] = [[node1.Id, node1.X, node1.Y, node1.Z],
                    #                         [node2.Id, node2.X, node2.Y, node2.Z],
                    #                         [node3.Id, node3.X, node3.Y, node3.Z]]

                    elem_list.append([node1.X, node1.Y, node1.Z]) ###################################################################################
                    elem_list.append([node2.X, node2.Y, node2.Z]) ###################################################################################
                    elem_list.append([node3.X, node3.Y, node3.Z]) ###################################################################################

                    # elem_id += 1

            return elem_list    # JUST FOR THE CHECK

        """
        # [NG] new code under construction
        # we create a sub model part (if extract_center = True) and fill it with conditions that are inside r_buildings
        if extract_center:
            # we extract the center of the geometry and we add them in a sub model part
            center_cond = self.ModelPart.CreateSubModelPart("CenterCondition")

            for cond in self.ModelPart.GetSubModelPart("BottomModelPart").Conditions:
                list_nodes = []     # node ids of the condition
                for node in cond.GetNodes():
                    list_nodes.append(node.Id)
                    dist = math.sqrt(((x_center-node.X)**2) + ((y_center-node.Y)**2))
                    
                    if (dist > r_buildings):
                        break   # if at least one node is external to r_buildings, we reject the entire condition
                        
                else:
                    # we get here if all nodes of the condition are inside the r_buildings
                    center_cond.AddNodes(list_nodes)
                    center_cond.AddCondition(cond, 0)
        """



#########################################################################################################################################################################
    def MeshCircleWithTerrainPoints_old(self, h_value=0.0, circ_division=60, num_sector=12, extract_center=False):
        "A function that creates a cylindrical volumetric mesh"
        # TODO: check the differences with "MeshCircleWithTerrainPoints"

        ### reading points from model part
        X = []; Y = []; Z = []      # lists are generated to use predefined operations
        ids_all = []                # id node with Z!=0
        coord_2D = []               # list with XY coordinates of vertices [X, Y]
        for node in self.ModelPart.Nodes:
            X.append(node.X)
            Y.append(node.Y)
            Z.append(node.Z)
            ids_all.append(node.Id)
            coord_2D.append([node.X, node.Y])

        # calculating minimum and maximum value of X, Y and Z coordinates of the entire initial domain
        x_min = min(X);    x_max = max(X)
        y_min = min(Y);    y_max = max(Y)
        # z_min = min(Z);    z_max = max(Z)     # EDIT 02 SEPTEMBER 2019

        # the centre of the circle
        self.x_center = (x_min + x_max)/2
        self.y_center = (y_min + y_max)/2

        # TODO: DEVONO ESSERE IMPOSTATI TRAMITE IL FILE JSON! NON QUI
        # radius calculation, considering the smaller side of the rectangle
        self.r_boundary = min((x_max-x_min), (y_max-y_min))/2       # radius of the domain
        # self.r_ground = self.r_boundary * 2.0 / 3.0                 # 2/3 of the radius of the domain (this value can be different)
        # self.r_buildings = self.r_boundary / 3.0                    # 1/3 of the radius of the domain (this value can be different)
        self.r_ground = self.r_boundary * 8.0 / 10.0                 # test values
        self.r_buildings = self.r_boundary * 6.0 / 10.0              # test values

        # buffering zone
        delta = self.r_boundary/20         # evaluate this value!!!
        print("[DEBUG][MeshCircleWithTerrainPoints_old] delta: ", delta)

        """ we delete nodes outside the circumference """
        del_id = []             # list where there are the ids to be deleted because are outside of r_boundary
        Z = []
        # for index in idx_to_smooth:
        for node in self.ModelPart.Nodes:
            x_curr = node.X             # current X coordinate
            y_curr = node.Y             # current Y coordinate

            dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))     # calculate the distance between current node and center of circle

            if (dist > self.r_boundary-delta):
                del_id.append(node.Id)
                continue
            else:
                # we fill the Z list to compute the minimum with the only nodes inside the domain
                Z.append(node.Z)        # this list is useful to compute the min(Z) after

        # the nodes from SubModelPart("bottom_terrain") that are outside the r_boundary-delta are deleted
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,False)

        for id in del_id:
            self.ModelPart.GetNode(id).Set(KratosMultiphysics.TO_ERASE,True)

        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

        # JUST FOR THE TESTS (CHECK IT)
        if Z:
            z_min = min(Z)
            z_max = max(Z)      # EDIT 02 SEPTEMBER 2019
        else:
            z_min = 0

        """ SMOOTHING PROCEDURE """
        for node in self.ModelPart.Nodes:
            x_curr = node.X
            y_curr = node.Y
            dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))
            if dist > self.r_ground:
                Z_beta = (-(dist-self.r_ground) / (self.r_boundary-self.r_ground))+1       # we calculate beta
                self.ModelPart.GetNode(node.Id).Z = (node.Z - z_min) * Z_beta + z_min

        # this step is useful to the n sectors
        circ_division = num_sector * round(circ_division / num_sector)
        # a list with the division points of the circumference
        theta = self._custom_range(0.0, 2*math.pi, 2*math.pi/circ_division)  # circ_division is the number of division of the 2*pi

        # lists with nodes on circle boundary
        x_circle = []; y_circle = []
        for th in theta:
            x_circle.append(self.r_boundary * math.cos(th) + self.x_center)   # we move the node along X to consider that the circle have not the centre in (0.0, 0.0)
            y_circle.append(self.r_boundary * math.sin(th) + self.y_center)   # we move the node along Y to consider that the circle have not the centre in (0.0, 0.0)

        """ TRIANGLE """
        coord_2D = []       # vertically list with X and Y coordinates (TO AVOID NUMPY VSTACK)
        all_points = []     # list with all coordinates of nodes
        # Bottom
        for node in self.ModelPart.Nodes:
            coord_2D.append([node.X, node.Y])     # we add coordinates X and Y of terrain nodes
            all_points.append((node.X, node.Y, node.Z))

        # Bottom: we add circumference nodes into coord_2D list
        for i in range(len(x_circle)):
            coord_2D.append([x_circle[i], y_circle[i]])

        ### preparing nodes to be handed over to triangle
        A = dict(vertices = coord_2D)
        triangle_dict = triangle.triangulate(A)     # triangle_dict {vertices: [...], triangles: [...], vertex_markers: [...]}

        inner_id_terrain = len(all_points)		# id up to this moment
        list_id_circle_bottom = []				# list with only ids of the bottom circle

        # Bottom
        for i in range(len(x_circle)):
            all_points.append((x_circle[i], y_circle[i], z_min))
            list_id_circle_bottom.append(i + inner_id_terrain)

        all_facets = []		# list with all node ids of faces
        all_markers = []	# list of integers [1 = bottom, 2 = topper and 3 = lateral]

        # Bottom: fill in the "facets" and "markers" lists
        for face in triangle_dict["triangles"]:
            all_facets.append([face[0], face[1], face[2]])
            all_markers.append(1)

        number_node_terrain = len(all_points)		# number of points up to now
        volume_height = z_max + h_value             # CHECK THIS VALUE

        # we fill all_points list for "top"
        for node in self.ModelPart.Nodes:
            all_points.append((node.X, node.Y, volume_height))

        inner_id_top = len(all_points)              # id up to this moment
        list_id_circle_top = []                     # list with only ids of the top circle

        for i in range(len(x_circle)):
            all_points.append((x_circle[i], y_circle[i], volume_height))
            list_id_circle_top.append(i + inner_id_top)

        # we fill all_facets list with "top" faces
        for face in triangle_dict["triangles"]:
            all_facets.append([face[0]+number_node_terrain, face[1]+number_node_terrain, face[2]+number_node_terrain])		# list of lists. "...+number_node_terrain" is UGLY! improve it
            all_markers.append(2)
        
        # we fill all_facets list with "lateral" faces
        for i in range(len(list_id_circle_bottom)):
            if (i == len(list_id_circle_bottom)-1):
                # here we have the last value
                all_facets.append([list_id_circle_bottom[i], list_id_circle_bottom[0], list_id_circle_top[0], list_id_circle_top[i]])
                all_markers.append(3)
                break
            all_facets.append([list_id_circle_bottom[i], list_id_circle_bottom[i+1], list_id_circle_top[i+1], list_id_circle_top[i]])
            all_markers.append(3)

        """ TETGEN (VIA MESHPY) """
        ### using Meshpy for the creation of an initial mesh https://documen.tician.de/meshpy/
        mesh_info = MeshInfo()
        mesh_info.set_points(all_points)
        mesh_info.set_facets(all_facets, markers=all_markers)

        # build the mesh
        mesh = build(mesh_info)
        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()

        lateral_model_part = self.ModelPart.CreateSubModelPart("LateralModelPart")
        bottom_model_part = self.ModelPart.CreateSubModelPart("BottomModelPart")
        top_model_part = self.ModelPart.CreateSubModelPart("TopModelPart")
        properties = self.ModelPart.Properties[0]

        bottom_cond = [];   bottom_points = []
        top_cond = [];      top_points = []
        lateral_cond = [];  lateral_points = []

        # AWARE: Shift by 1 in index!!!
        ### Nodes
        for i in range(len(mesh.points)):
            # node id is shifted up by 1
            coords = mesh.points[i]
            self.ModelPart.CreateNewNode((i+1), coords[0], coords[1], coords[2])

        ### Conditions and Nodes in SubModelParts
        for j in range(len(mesh.faces)):
            points = mesh.faces[j]
            marker = mesh.face_markers[j]
            self.ModelPart.CreateNewCondition("WallCondition3D3N", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties)
            if (marker == 1):
                bottom_cond.append(j+1)
            elif (marker == 2):
                top_cond.append(j+1)
            elif (marker == 3):
                lateral_cond.append(j+1)
            else:
                KratosMultiphysics.Logger.PrintWarning("GeoMesher", "Marker {} not valid. 1 = bottom, 2 = topper and 3 = lateral".format(marker))

        bottom_model_part.AddConditions(bottom_cond)
        top_model_part.AddConditions(top_cond)
        lateral_model_part.AddConditions(lateral_cond)

        for cond in bottom_model_part.Conditions:
            for node in cond.GetNodes():
                bottom_points.append(node.Id)
        for cond in top_model_part.Conditions:
            for node in cond.GetNodes():
                top_points.append(node.Id)
        for cond in lateral_model_part.Conditions:
            for node in cond.GetNodes():
                lateral_points.append(node.Id)

        bottom_model_part.AddNodes(bottom_points)
        top_model_part.AddNodes(top_points)
        lateral_model_part.AddNodes(lateral_points)

        ### Elements
        for k in range(len(mesh.elements)):
            points = mesh.elements[k]
            self.ModelPart.CreateNewElement("Element3D4N", (k+1), [points[0]+1, points[1]+1, points[2]+1, points[3]+1], properties)


        # """ OPTIONAL VALUES """
        # # [NG] code under construction
        # # we create a sub model part (if extract_center = True) and fill it with conditions that are inside r_buildings
        # if extract_center:
        #     # we extract the center of the geometry and we add them in a sub model part
        #     center_cond = self.ModelPart.CreateSubModelPart("CenterCondition")

        #     for cond in self.ModelPart.GetSubModelPart("BottomModelPart").Conditions:
        #         list_nodes = []     # node ids of the condition\
        #         for node in cond.GetNodes():
        #             list_nodes.append(node.Id)
        #             dist = math.sqrt(((self.x_center-node.X)**2) + ((self.y_center-node.Y)**2))
                    
        #             # if (dist > r_ground):       # Only for a test. The value correcr is r_buildings
        #             if (dist > self.r_buildings):
        #                 break   # if at least one node is external to r_buildings, we reject the entire condition
                        
        #         else:
        #             # we get here if all nodes of the condition are inside the r_buildings
        #             center_cond.AddNodes(list_nodes)
        #             center_cond.AddCondition(cond, 0)





    def MeshCircleWithTerrainPoints_MOD(self, h_value=0.0, circ_division=60, num_sector=12):
        "A function that creates a cylindrical volumetric mesh"
        # TODO: check the differences with "MeshCircleWithTerrainPoints"

        ### reading points from model part
        X = []; Y = []; Z = []      # lists are generated to use predefined operations
        ids_all = []                # id node with Z!=0
        coord_2D = []               # list with XY coordinates of vertices [X, Y]
        for node in self.ModelPart.Nodes:
            X.append(node.X)
            Y.append(node.Y)
            Z.append(node.Z)
            ids_all.append(node.Id)
            coord_2D.append([node.X, node.Y])

        # calculating minimum and maximum value of X, Y and Z coordinates of the entire initial domain
        x_min = min(X);    x_max = max(X)
        y_min = min(Y);    y_max = max(Y)
        # z_min = min(Z);    z_max = max(Z)     # EDIT 02 SEPTEMBER 2019

        # the centre of the circle
        self.x_center = (x_min + x_max)/2
        self.y_center = (y_min + y_max)/2

        # TODO: DEVONO ESSERE IMPOSTATI TRAMITE IL FILE JSON! NON QUI
        # radius calculation, considering the smaller side of the rectangle
        # self.r_boundary = min((x_max-x_min), (y_max-y_min))/2       # radius of the domain
        # self.r_ground = self.r_boundary * 2.0 / 3.0                 # 2/3 of the radius of the domain (this value can be different)
        # self.r_buildings = self.r_boundary / 3.0                    # 1/3 of the radius of the domain (this value can be different)
        # self.r_ground = self.r_boundary * 8.0 / 10.0                 # test values
        # self.r_buildings = self.r_boundary * 6.0 / 10.0              # test values
        print("r_boundary: ", self.r_boundary)
        print("r_ground: ", self.r_ground)
        print("r_buildings: ", self.r_buildings)


        # buffering zone
        delta = self.r_boundary/20         # evaluate this value!!!
        print("[DEBUG][MeshCircleWithTerrainPoints_old] delta: ", delta)

        """ we delete nodes outside the circumference """
        del_id = []             # list where there are the ids to be deleted because are outside of r_boundary
        Z = []
        # for index in idx_to_smooth:
        for node in self.ModelPart.Nodes:
            x_curr = node.X             # current X coordinate
            y_curr = node.Y             # current Y coordinate

            dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))     # calculate the distance between current node and center of circle

            if (dist > self.r_boundary-delta):
                del_id.append(node.Id)
                continue
            else:
                # we fill the Z list to compute the minimum with the only nodes inside the domain
                Z.append(node.Z)        # this list is useful to compute the min(Z) after

        # the nodes from SubModelPart("bottom_terrain") that are outside the r_boundary-delta are deleted
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,False)

        for id in del_id:
            self.ModelPart.GetNode(id).Set(KratosMultiphysics.TO_ERASE,True)

        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

        # JUST FOR THE TESTS (CHECK IT)
        if Z:
            z_min = min(Z)
            z_max = max(Z)      # EDIT 02 SEPTEMBER 2019
        else:
            z_min = 0
            z_max = 0

        """ SMOOTHING PROCEDURE """
        for node in self.ModelPart.Nodes:
            x_curr = node.X
            y_curr = node.Y
            dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))
            if dist > self.r_ground:
                Z_beta = (-(dist-self.r_ground) / (self.r_boundary-self.r_ground))+1       # we calculate beta
                self.ModelPart.GetNode(node.Id).Z = (node.Z - z_min) * Z_beta + z_min

        # this step is useful to the n sectors
        circ_division = num_sector * round(circ_division / num_sector)
        # a list with the division points of the circumference
        theta = self._custom_range(0.0, 2*math.pi, 2*math.pi/circ_division)  # circ_division is the number of division of the 2*pi

        # lists with nodes on circle boundary
        x_circle = []; y_circle = []
        for th in theta:
            x_circle.append(self.r_boundary * math.cos(th) + self.x_center)   # we move the node along X to consider that the circle have not the centre in (0.0, 0.0)
            y_circle.append(self.r_boundary * math.sin(th) + self.y_center)   # we move the node along Y to consider that the circle have not the centre in (0.0, 0.0)

        """ TRIANGLE """
        coord_2D_bottom = []    # list with X and Y coordinates (only bottom)
        coord_2D_top = []       # list with X and Y coordinates (only top)
        all_points_3D = []      # list with all 3D coordinates of nodes

        all_facets = []         # list with all node ids of faces
        all_markers = []        # list of integers [1 = bottom, 2 = topper and 3 = lateral]

        # circumference nodes
        for i in range(len(x_circle)):
            coord_2D_bottom.append([x_circle[i], y_circle[i]])          # Bottom
            coord_2D_top.append([x_circle[i], y_circle[i]])             # Top
            all_points_3D.append((x_circle[i], y_circle[i], z_min))     # Bottom

        # Bottom: inner nodes
        for node in self.ModelPart.Nodes:
            coord_2D_bottom.append([node.X, node.Y])     # we add coordinates X and Y of bottom
            all_points_3D.append((node.X, node.Y, node.Z))

        # Bottom: triangle
        A = dict(vertices = coord_2D_bottom)
        triangle_dict_bottom = triangle.triangulate(A)     # triangle_dict_bottom {vertices: [...], triangles: [...], vertex_markers: [...]}

        # Bottom: number of bottom nodes
        n_bottom_nodes = len(triangle_dict_bottom["vertices"])

        # Bottom: fill in the "facets" and "markers" lists
        for face in triangle_dict_bottom["triangles"]:
            all_facets.append([face[0], face[1], face[2]])
            all_markers.append(1)

        # Top: triangle
        volume_height = z_max + h_value
        B = dict(vertices = coord_2D_top)
        triangle_dict_top = triangle.triangulate(B)     # triangle_dict_top {vertices: [...], triangles: [...], vertex_markers: [...]}
        
        # Top: updated all_points_3D with top nodes
        for node in triangle_dict_top["vertices"]:
            all_points_3D.append((node[0], node[1], volume_height))
        
        # Top: fill in the "facets" and "markers" lists
        for face in triangle_dict_top["triangles"]:
            all_facets.append([ face[0]+n_bottom_nodes,
                                face[1]+n_bottom_nodes,
                                face[2]+n_bottom_nodes])        # list of lists. "...+n_bottom_nodes"
            all_markers.append(2)

        for n, coord in enumerate(coord_2D_top[:-1]):
            all_facets.append([n, n+1, n+n_bottom_nodes+1, n+n_bottom_nodes])
            all_markers.append(3)
        # the last value
        all_facets.append([n+1, 0, n_bottom_nodes, n+n_bottom_nodes+1])
        all_markers.append(3)


        """ TETGEN (VIA MESHPY) """
        ### using Meshpy for the creation of an initial mesh https://documen.tician.de/meshpy/
        mesh_info = MeshInfo()
        mesh_info.set_points(all_points_3D)
        mesh_info.set_facets(all_facets, markers=all_markers)

        # build the mesh
        mesh = build(mesh_info)
        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()

        lateral_model_part = self.ModelPart.CreateSubModelPart("LateralModelPart")
        bottom_model_part = self.ModelPart.CreateSubModelPart("BottomModelPart")
        top_model_part = self.ModelPart.CreateSubModelPart("TopModelPart")
        properties = self.ModelPart.Properties[0]

        bottom_cond = [];   bottom_points = []
        top_cond = [];      top_points = []
        lateral_cond = [];  lateral_points = []

        # AWARE: Shift by 1 in index!!!
        ### Nodes
        for i in range(len(mesh.points)):
            # node id is shifted up by 1
            coords = mesh.points[i]
            self.ModelPart.CreateNewNode((i+1), coords[0], coords[1], coords[2])

        ### Conditions and Nodes in SubModelParts
        for j in range(len(mesh.faces)):
            points = mesh.faces[j]
            marker = mesh.face_markers[j]
            self.ModelPart.CreateNewCondition("WallCondition3D3N", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties)
            if (marker == 1):
                bottom_cond.append(j+1)
            elif (marker == 2):
                top_cond.append(j+1)
            elif (marker == 3):
                lateral_cond.append(j+1)
            else:
                KratosMultiphysics.Logger.PrintWarning("GeoMesher", "Marker {} not valid. 1 = bottom, 2 = topper and 3 = lateral".format(marker))

        bottom_model_part.AddConditions(bottom_cond)
        top_model_part.AddConditions(top_cond)
        lateral_model_part.AddConditions(lateral_cond)

        for cond in bottom_model_part.Conditions:
            for node in cond.GetNodes():
                bottom_points.append(node.Id)
        for cond in top_model_part.Conditions:
            for node in cond.GetNodes():
                top_points.append(node.Id)
        for cond in lateral_model_part.Conditions:
            for node in cond.GetNodes():
                lateral_points.append(node.Id)

        bottom_model_part.AddNodes(bottom_points)
        top_model_part.AddNodes(top_points)
        lateral_model_part.AddNodes(lateral_points)

        ### Elements
        for k in range(len(mesh.elements)):
            points = mesh.elements[k]
            self.ModelPart.CreateNewElement("Element3D4N", (k+1), [points[0]+1, points[1]+1, points[2]+1, points[3]+1], properties)











    def MeshSectors(self, n_sectors=12):
        """ function to split the LateralModelPart into n SubModelParts

        Args:
            n_sectors: number of sectors into which the ModelPart must be divided

        Returns:
            ModelPart updated with n SubModelPart for each sector
        """

        if not (self.ModelPart.HasSubModelPart("LateralModelPart")):
            KratosMultiphysics.Logger.PrintWarning("GeoMesher", "LateralModelPart not found!")
            return
        
        # creation of the n sub model parts
        for n in range(1, n_sectors+1):
            if not (self.ModelPart.HasSubModelPart("LateralSector_{}".format(n))):
                self.ModelPart.CreateSubModelPart("LateralSector_{}".format(n))
        
        center_domain = (self.x_center, self.y_center)

        # sector with all intervals
        list_sector = self._frange(0.0, 360.0, 360.0/n_sectors)
        sects = list_sector + [360]     # this is necessary for the last sector. for example in 12 sectors, the last one will be 330.0, 360.0
        
        # we populate a dictionary with all sectors. key: sector id; value: range of angles (first angle:included; second angle: not included)
        sect_id = 1
        dict_sects = {}
        for i in range(len(list_sector)):
            dict_sects[sect_id] = (sects[i], sects[i+1])
            sect_id += 1

        for cond in self.ModelPart.GetSubModelPart("LateralModelPart").Conditions:
            coords = cond.GetGeometry().Center()    # center of the n-th condition
            center_coords = (coords.X, coords.Y)
            Comp_x, Comp_y = self._vector_components(center_coords, center_domain)
            angle = self._angle_vectors(Comp_x, Comp_y)

            for n_sect, angles in dict_sects.items():
                if angles[0] <= angle < angles[1]:
                    # we add the current Condition in the specific sub model part
                    sub_model_part_sect = self.ModelPart.GetSubModelPart("LateralSector_{}".format(n_sect))
                    # sub_model_part_sect.AddCondition(cond.Id)
                    sub_model_part_sect.AddCondition(cond)

                    # we add nodes in the specific sub model part
                    points = []
                    for node in cond.GetNodes():
                        points.append(node.Id)
                    sub_model_part_sect.AddNodes(points)

                    break
        
        # LateralModelPart has been split and therefore can be removed
        self.ModelPart.RemoveSubModelPart("LateralModelPart")
 
#########################################################################################################################################################################


    def extract_center(self, terrain_model_part):
        """ function to create a ModelPart with "elements" of the central portion

        Args:
            terrain_model_part: ModelPart of the terrain

        Returns:
            ModelPart with the elements belonging to the central portion
        """
        
        model = KratosMultiphysics.Model()
        self.center_model_part = model.CreateModelPart("center_model_part")
        prop = self.center_model_part.Properties[0]
        smp_bottom = terrain_model_part.GetSubModelPart("BottomModelPart")

        for cond in smp_bottom.Conditions:
            nodes = cond.GetNodes()
            for node in nodes:
                dist = math.sqrt(((self.x_center-node.X)**2) + ((self.y_center-node.Y)**2))
                if (dist > self.r_buildings*1.2):
                    break
            else:
                for node in nodes:
                    self.center_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                self.center_model_part.CreateNewElement("Element2D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)




############################ --- Auxiliary functions --- ###########################################

    def _set_variational_distance_process_serial(self, complete_model, aux_name):
        "Construct the variational distance calculation process"

        serial_settings = KratosMultiphysics.Parameters("""
            {
                "linear_solver_settings"   : {
                    "solver_type" : "amgcl"
                }
            }
        """)
        from KratosMultiphysics import python_linear_solver_factory
        linear_solver = python_linear_solver_factory.ConstructSolver(serial_settings["linear_solver_settings"])

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


    def _custom_range(self, start, stop, step=1.0):
        "function to replace the numpy range function"
        list_floats = [float(start)]
        while (list_floats[-1]+step) < stop:
            list_floats.append(list_floats[-1]+step)
        return list_floats


    def _vector_components(self, P1, P2):
        "function to compute the vector components. return a tuble"
        if len(P1) == 2:
            return (P1[0]-P2[0], P1[1]-P2[1])                  # 2D case
        else:
            return (P1[0]-P2[0], P1[1]-P2[1], P1[2]-P2[2])     # 3D case


    def _angle_vectors(self, Px, Py):
        "angle between two vectors. return an angle in degree"
        angle = math.degrees(math.atan2(Py, Px))
        return angle + 360 if angle < 0 else angle


    def _frange(self, start, stop=None, step=None):
        "floating range. Like 'range' function but with float. return a list"

        start = round(start, 5)      # we make sure we have only 5 decimal digits
        if stop == None:
            stop = start
            start = 0.0
        else:
            stop = round(stop, 5)    # we make sure we have only 5 decimal digits
        
        if step == None:
            step = 1.0
        else:
            step = round(step, 5)    # we make sure we have only 5 decimal digits

        list_range = [start]
        last = start
        
        while True:
            last = round(last+step, 5)   # we make sure we have only 5 decimal digits
            if last >= stop:
                break
            list_range.append(last)
        
        return list_range