import enum
from os import curdir
from time import perf_counter

from triangle.plot import holes, vertices
import KratosMultiphysics
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.MeshingApplication as KratosMesh

if KratosMultiphysics.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI
    print("[DEBUG][geo_mesher] import KratosMultiphysics.mpi")

from KratosMultiphysics.GeodataProcessingApplication.geo_processor import GeoProcessor
from KratosMultiphysics.GeodataProcessingApplication.geo_building import GeoBuilding

import math
from meshpy.tet import MeshInfo, Options, build
import numpy as np
import triangle

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
        self.height = float()       # height value of the domain


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
                # we fill the Z list to compute the min/max with the only nodes inside the domain
                Z.append(node.Z)        # this list is useful to compute the min(Z)/max(Z) after

        # the nodes from SubModelPart("bottom_terrain") that are outside the r_boundary-delta are deleted
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,False)

        for id in del_id:
            self.ModelPart.GetNode(id).Set(KratosMultiphysics.TO_ERASE,True)

        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

        # JUST FOR THE TESTS (CHECK IT)
        if Z:
            z_min = min(Z)
            z_max = max(Z)
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
        # reference: https://rufat.be/triangle/examples.html
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


    # TODO: ELIMINARE QUESTA FUNZIONE
    # TEST!
    def MeshCircle(self, h_value=0.0, circ_division=60, num_sector=12):
        "mesh circle from terrain surface (with buildings)"
        # TODO: it is just a test at the moment (04 aug 2021)

        ### reading points from model part
        X = []; Y = []; Z = []      # lists are generated to use predefined operations
        ids_all = []                # id node with Z!=0
        coord_2D = []               # list with XY coordinates of vertices [X, Y]
        for node in self.ModelPart.Nodes:
            X.append(node.X)
            Y.append(node.Y)
            Z.append(node.Z)
            ids_all.append(node.Id)
            coord_2D.append([node.X, node.Y])   # TODO: CONTROLLARE SE I NODI DEGLI EDIFICI DEVONO ESSERE CONSIDERATI O MENO

        # calculating minimum and maximum value of X, Y and Z coordinates of the entire initial domain
        x_min = min(X);    x_max = max(X)
        y_min = min(Y);    y_max = max(Y)

        # the centre of the circle
        self.x_center = (x_min + x_max)/2
        self.y_center = (y_min + y_max)/2


        print("r_boundary: ", self.r_boundary)
        print("r_ground: ", self.r_ground)
        print("r_buildings: ", self.r_buildings)


        # buffering zone
        delta = self.r_boundary/20         # evaluate this value!!!
        print("[DEBUG][MeshCircle] delta: ", delta)

        Z = []

        # we set all nodes and elements as NOT to erase
        for elem in self.ModelPart.Elements:
            for node in elem.GetNodes():
                node.Set(KratosMultiphysics.TO_ERASE,False)
            elem.Set(KratosMultiphysics.TO_ERASE,False)
        
        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()     # list of nodes
            for node in nodes:
                x_curr = node.X         # current X coordinate
                y_curr = node.Y         # current Y coordinate

                dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))     # calculate the distance between current node and center of circle

                if (dist > self.r_boundary-delta):
                    node.Set(KratosMultiphysics.TO_ERASE,True)
                    elem.Set(KratosMultiphysics.TO_ERASE,True)
                else:
                    # we fill the Z list to compute the min/max with the only nodes inside the domain
                    Z.append(node.Z)        # this list is useful to compute the min(Z)/max(Z) after
        
        # we delete nodes and element marked as TO_ERASE
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.ModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)

        # we update z_min and z_max
        if Z:
            z_min = min(Z)
            z_max = max(Z)
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
        for th in theta[:-1]:   # the first and last values are the same
            x_circle.append(self.r_boundary * math.cos(th) + self.x_center)   # we move the node along X to consider that the circle have not the centre in (0.0, 0.0)
            y_circle.append(self.r_boundary * math.sin(th) + self.y_center)   # we move the node along Y to consider that the circle have not the centre in (0.0, 0.0)

        list_edges = list()     # list with all edges (an edge is a couple odf nodes)
        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()

            # we save edges with ordered nodes
            list_edges.append((nodes[0].Id, nodes[1].Id)) if (nodes[0].Id < nodes[1].Id) else list_edges.append((nodes[1].Id, nodes[0].Id))
            list_edges.append((nodes[1].Id, nodes[2].Id)) if (nodes[1].Id < nodes[2].Id) else list_edges.append((nodes[2].Id, nodes[1].Id))
            list_edges.append((nodes[2].Id, nodes[0].Id)) if (nodes[2].Id < nodes[0].Id) else list_edges.append((nodes[0].Id, nodes[2].Id))

        ext_edges = list()      # list with external edges
        ext_nodes = set()       # list with external nodes
        for edge in list_edges:
            if (list_edges.count(edge) == 1) and (not edge in ext_edges):
                # if the edge is unique and not yet in ext_edges
                ext_edges.append(edge)      # we update ext_edges with unique edge
                ext_nodes.add(edge[0])      # we update ext_nodes with nodes
                ext_nodes.add(edge[1])

        print("[DEBUG] external edges: ", ext_edges)    # DEBUG. TODO: DELETE IT
        print("[DEBUG] external nodes: ", ext_nodes)    # DEBUG. TODO: DELETE IT


        """ TRIANGLE """
        # reference: https://rufat.be/triangle/examples.html
        coord_2D_bottom = list()    # list with X and Y coordinates (only the external portion on bottom)
        coord_2D_top = list()       # list with X and Y coordinates (only top)
        all_points_3D = list()      # list with all 3D coordinates of nodes (useful for meshpy-TetGen)

        all_facets = list()         # list with all node ids of faces
        all_markers = list()        # list of integers [1 = bottom, 2 = topper and 3 = lateral]

        segments = list()           # list of segments delimiting triangulations (Bottom)

        len_circ = len(x_circle)    # x_circle and y_circle have the same number of nodes

        # circumference nodes and segments on circumference
        for id, (x_value, y_value) in enumerate(zip(x_circle, y_circle)):
            coord_2D_bottom.append([x_value, y_value])          # Bottom
            coord_2D_top.append([x_value, y_value])             # Top
            all_points_3D.append([x_value, y_value, z_min])
            
            if (id == len_circ-1):     # the last one
                segments.append([id, 0])    # the last value with the first to close the circle
                continue
            segments.append([id, id+1])     # the i-th segment
        

        temp_dict = dict()          # SALVO LA CORRISPONDENTE POSIZIONE IN coord_2D_bottom E L'ID NEL MODELPART. key: id nel model part; value: id in coord_2D_bottom
        # inv_dict = dict()           # DIZIONARIO INVERSO RISPETTO A temp_dict. key: id in coord_2D_bottom; value: id nel model part
        points_3D_to_mdpa = dict()  # key: posizione in all_points_3D; value: id in model part

        # AGGIUNGO TUTTI I NODI DEL MODEL PART IN coord_2D_bottom
        # SE IL NODO I-ESIMO È IN ext_nodes ALLORA MI SALVO LA POSIZIONE CHE OCCUPA NELLA LISTA coord_2D_bottom (MI SERVIRÀ PER I SEGMENTI)
        x_hole = 0; y_hole = 0      # coordinates of hole (useful in Triangle)
        for node in self.ModelPart.Nodes:
            coord_2D_bottom.append([node.X, node.Y])     # we add coordinates X and Y of bottom
            all_points_3D.append([node.X, node.Y, node.Z])
            points_3D_to_mdpa[len(all_points_3D)-1] = node.Id
            
            if (node.Id in ext_nodes):
                temp_dict[node.Id] = len(coord_2D_bottom)-1
                # inv_dict[len(coord_2D_bottom)-1] = node.Id
                
                x_hole += node.X
                y_hole += node.Y

        # we calculate the centroid of the "hole"
        len_ext_nodes = len(ext_nodes)
        x_hole /= len_ext_nodes
        y_hole /= len_ext_nodes

        # # in ext_edges there are the node pairs belonging to the segments
        # # nodes in ext_nodes have been added to the coord_2D_bottom list
        # # the position (in ext_nodes) corresponding to the pair in ext_edges is read
        # ext_nodes = list(ext_nodes)     # we transform the set into list because we need the "index" function
        # for id_1, id_2 in ext_edges:
        #     node_1 = ext_nodes.index(id_1)
        #     node_2 = ext_nodes.index(id_2)
        #     segments.append([node_1+len_circ, node_2+len_circ])

        # in ext_edges there are the node pairs belonging to the segments
        # LEGGO L'ID CORRISPONDENTE A coord_2D_bottom E CREO I SEGMENTI
        for id_1, id_2 in ext_edges:
            node_1 = temp_dict[id_1]
            node_2 = temp_dict[id_2]
            segments.append([node_1, node_2])


        # Bottom: triangle
        A = dict(vertices=coord_2D_bottom, segments=segments, holes=[[x_hole,y_hole]])
        triangle_dict_bottom = triangle.triangulate(A, "p")     # triangle_dict_bottom {vertices: [...], triangles: [...], vertex_markers: [...]}


        # TODO: DELETE IT
        import matplotlib.pyplot as plt
        triangle.compare(plt, A, triangle_dict_bottom)
        plt.show()
        input("STOP 1")

        print(triangle_dict_bottom)
        input("STOP 2")


        max_node_id = self._find_max_node_id() + 1
        max_elem_id = self._find_max_elem_id() + 1
        properties = self.ModelPart.Properties[0]

        # AGGIUNGO GLI ELEMENTI MANCANTI NEL MODEL PART (GLI ELEMENTI APPENA CREATI CON triangle)
        # INOLTRE, AGGIUNGO I NODI SULLA CIRCONFERENZA NEL MODEL PART
        for nodes in triangle_dict_bottom["triangles"]:
            nodes_for_element = list()
            for node in nodes:
                if node in points_3D_to_mdpa:           # SE IL NODO È IN points_3D_to_mdpa VUOL DIRE CHE È UN NODO INTERNO (NON SULLA CIRCONFERENZA)
                    node_id = points_3D_to_mdpa[node]   # LEGGO IL CORRISPONDENTE ID PER IL MODEL PART
                    nodes_for_element.append(node_id)   # nodes useful to create a new element
                    continue
                else:   # SE NEL MODEL PART IL NODO NON C'È, VUOL DIRE CHE È UN NODO SULLA CIRCONFERENZA. LO AGGIUNGO NEL MODEL PART
                    coords = triangle_dict_bottom["vertices"][node]   # LEGGO IL NODO NELLA LISTA IN vertices
                    self.ModelPart.CreateNewNode(max_node_id, coords[0], coords[1], z_min)  # SE IL NODO NON È NEL MODELPART VUOL DIRE CHE È SULLA CIRCONFERENZA E QUINDI HA LA Z PARI A z_min
                    # all_points_3D.append([coords[0], coords[1], z_min])
                    # points_3D_to_mdpa[len(all_points_3D)-1] = max_node_id
                    nodes_for_element.append(max_node_id)   # nodes useful to create a new element
                    max_node_id += 1    # we update the Id
            
            self.ModelPart.CreateNewElement("Element2D3N", max_elem_id, nodes_for_element, properties)
            max_elem_id += 1
        
        # reverse dictionary of points_3D_to_mdpa 
        mdpa_to_points_3D = {points_3D_to_mdpa[k]: k for k in points_3D_to_mdpa}   # key: node id in model part; value: node position in all_points_3D list
        return

        # Bottom: fill in the "facets" and "markers" lists
        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()
            faces = list()
            for node in nodes:
                faces.append(mdpa_to_points_3D[node.Id])    # LEGGO IL CORRISPONDETE node.Id DEL MODEL PART IN all_points_3D
            
            all_facets.append(faces)
            all_markers.append(1)       # marker: 1 bottom; 2 top; 3 lateral
        
        # Bottom: number of bottom nodes
        n_bottom_nodes = len(all_points_3D)     # so far there are only bottom nodes

        # Top: triangle
        volume_height = z_max + h_value
        B = dict(vertices=coord_2D_top)
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

        ##############

        # DEVO SISTEMARE LA CREAZIONE DELLE FACCE LATERALI
        # NON CONOSCO PIÙ LA POSIZIONE DEI NODI SULLA CIRCONFERENZA IN BOTTOM PERCHÉ GLI HO AGGIUNTI INSIEME ALLE FACCE CREATE CON triangle

        for id, val in enumerate(all_points_3D):
            print(id, val)
        input("STOP 3")

        for n, coord in enumerate(coord_2D_top[:-1]):
            print(all_points_3D.index([coord[0], coord[1], z_min]))
        
        input("STOP 4")


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


    # TODO: tradurre commenti e aggiungerne altri per spiegare bene quello che viene fatto
    def MeshCircle_2(self, h_value=0.0, circ_division=60, num_sector=12):
        " VERSIONE MIGLIORATA E FUNZIONANTE DI MeshCircle "
        ### reading points from model part
        X = list(); Y = list()              # lists are generated to use predefined operations
        for node in self.ModelPart.Nodes:
            node.X = self.truncate(node.X)       # we truncate the coordinate to the fifth decimal digit
            node.Y = self.truncate(node.Y)       # we truncate the coordinate to the fifth decimal digit
            X.append(node.X)
            Y.append(node.Y)
        
        # calculating minimum and maximum value of X, Y and Z coordinates of the entire initial domain
        x_min = min(X);    x_max = max(X)
        y_min = min(Y);    y_max = max(Y)

        # the centre of the circle
        self.x_center = (x_min + x_max)/2
        self.y_center = (y_min + y_max)/2

        # buffering zone
        delta = self.r_boundary/20         # evaluate this value!!!
        print("[DEBUG][MeshCircle] delta: ", delta)

        # we set all nodes and elements as NOT to erase
        for elem in self.ModelPart.Elements:
            for node in elem.GetNodes():
                node.Set(KratosMultiphysics.TO_ERASE,False)
            elem.Set(KratosMultiphysics.TO_ERASE,False)

        Z = list()

        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()     # list of nodes
            for node in nodes:
                x_curr = node.X         # current X coordinate
                y_curr = node.Y         # current Y coordinate

                dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))     # calculate the distance between current node and center of circle

                if (dist > self.r_boundary-delta):
                    node.Set(KratosMultiphysics.TO_ERASE,True)
                    elem.Set(KratosMultiphysics.TO_ERASE,True)
                else:
                    # we fill the Z list to compute the min/max with the only nodes inside the domain
                    Z.append(node.Z)        # this list is useful to compute the min(Z)/max(Z) after
        
        # we delete nodes and element marked as TO_ERASE
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.ModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)

        # we update z_min and z_max
        if Z:
            z_min = min(Z)
            z_max = max(Z)
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
        for th in theta[:-1]:   # the first and last values are the same
            x = self.r_boundary * math.cos(th) + self.x_center      # we move the node along X to consider that the circle have not the centre in (0.0, 0.0)
            y = self.r_boundary * math.sin(th) + self.y_center      # we move the node along Y to consider that the circle have not the centre in (0.0, 0.0)
            x_circle.append(self.truncate(x))    # we add truncated coordinates
            y_circle.append(self.truncate(y))    # we add truncated coordinates

        list_edges = list()     # list with all edges (an edge is a couple odf nodes)
        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()
            # we save edges with ordered nodes
            list_edges.append((nodes[0].Id, nodes[1].Id)) if (nodes[0].Id < nodes[1].Id) else list_edges.append((nodes[1].Id, nodes[0].Id))
            list_edges.append((nodes[1].Id, nodes[2].Id)) if (nodes[1].Id < nodes[2].Id) else list_edges.append((nodes[2].Id, nodes[1].Id))
            list_edges.append((nodes[2].Id, nodes[0].Id)) if (nodes[2].Id < nodes[0].Id) else list_edges.append((nodes[0].Id, nodes[2].Id))

        ext_edges = list()      # list with external edges
        ext_nodes = set()       # list with external nodes
        for edge in list_edges:
            if (list_edges.count(edge) == 1) and (not edge in ext_edges):
                # if the edge is unique and not yet in ext_edges
                ext_edges.append(edge)      # we update ext_edges with unique edge
                ext_nodes.add(edge[0])      # we update ext_nodes with nodes
                ext_nodes.add(edge[1])

        # print("[DEBUG] external edges: ", ext_edges)    # DEBUG. TODO: DELETE IT
        # print("[DEBUG] external nodes: ", ext_nodes)    # DEBUG. TODO: DELETE IT


        """ TRIANGLE """
        # reference: https://rufat.be/triangle/examples.html

        # *** NB: coord_2D_bottom E all_points_3D HANNO GLI STESSI INDICI. L'UNICA COSA DIVERSA È CHE NEL 3D HO LA COORDINATA Z ANCHE
        coord_2D_bottom = list()    # list with X and Y coordinates (bottom)
        coord_2D_top = list()       # list with X and Y coordinates (top)
        all_points_3D = list()      # list with all 3D coordinates of nodes (useful for meshpy-TetGen)
        segments = list()           # list of segments delimiting triangulations (bottom)

        len_circ = len(x_circle)    # x_circle and y_circle have the same number of nodes

        # circumference nodes and segments on it
        for id, (x_value, y_value) in enumerate(zip(x_circle, y_circle)):
            coord_2D_bottom.append([x_value, y_value])          # Bottom
            coord_2D_top.append([x_value, y_value])             # Top
            all_points_3D.append([x_value, y_value, z_min])
            
            if (id == len_circ-1):     # the last one
                segments.append([id, 0])    # the last value with the first to close the circle
                continue
            segments.append([id, id+1])     # the i-th segment

        
        # AGGIUNGO TUTTI I NODI DEL MODEL PART IN coord_2D_bottom
        # SE IL NODO I-ESIMO È IN ext_nodes ALLORA MI SALVO LA POSIZIONE CHE OCCUPA NELLA LISTA coord_2D_bottom (MI SERVIRÀ PER I SEGMENTI)
        mdpa_to_coord_2D = dict()       # SALVO LA CORRISPONDENTE POSIZIONE IN coord_2D_bottom E L'ID NEL MODELPART. key: id nel model part; value: id in coord_2D_bottom
        coord_3D_to_mdpa = dict()       # key: POSIZIONE DEL NODO I-ESIMO NELLA LISTA all_points_3D; value: NODO ID NEL MODEL PART
        for node in self.ModelPart.Nodes:
            coord_2D_bottom.append([node.X, node.Y])     # we add coordinates X and Y of bottom
            all_points_3D.append([node.X, node.Y, node.Z])
            
            mdpa_to_coord_2D[node.Id] = len(coord_2D_bottom)-1
            coord_3D_to_mdpa[len(all_points_3D)-1] = node.Id

        # in ext_edges there are the node pairs belonging to the segments
        # LEGGO L'ID CORRISPONDENTE A coord_2D_bottom E CREO I SEGMENTI
        for id_1, id_2 in ext_edges:
            node_1 = mdpa_to_coord_2D[id_1]
            node_2 = mdpa_to_coord_2D[id_2]
            segments.append([node_1, node_2])
        
        # Bottom: triangle
        A = dict(vertices=coord_2D_bottom, segments=segments, holes=[[self.x_center,self.y_center]])
        triangle_dict_bottom = triangle.triangulate(A, "p")     # triangle_dict_bottom {vertices: [...], triangles: [...], vertex_markers: [...]}

        # # TODO: DELETE IT
        # import matplotlib.pyplot as plt
        # triangle.compare(plt, A, triangle_dict_bottom)
        # plt.show()
        # input("STOP 1")

        # print(triangle_dict_bottom)
        # input("STOP 2")

        # # PER VISUALIZZARE I NODI PRESENTI IN triangle_dict_bottom["vertices"]
        # p_x = []
        # p_y = []
        # for x,y in triangle_dict_bottom["vertices"]:
        #     p_x.append(x)
        #     p_y.append(y)
        #     plt.text(x, y+0.5, "({}, {})".format(x, y))
        # plt.plot(p_x, p_y, "r*")
        # plt.show()
        # input("PAUSE")

        max_node_id = self._find_max_node_id() + 1
        max_elem_id = self._find_max_elem_id() + 1
        properties = self.ModelPart.Properties[0]

        # AGGIUNGO GLI ELEMENTI MANCANTI NEL MODEL PART (GLI ELEMENTI APPENA CREATI CON triangle)
        # INOLTRE, AGGIUNGO I NODI SULLA CIRCONFERENZA NEL MODEL PART
        for nodes in triangle_dict_bottom["triangles"]:
            nodes_for_element = list()
            for node in nodes:
                if node in coord_3D_to_mdpa:           # SE IL NODO È IN coord_3D_to_mdpa VUOL DIRE CHE È UN NODO INTERNO (NON SULLA CIRCONFERENZA)
                    node_id = coord_3D_to_mdpa[node]   # LEGGO IL CORRISPONDENTE ID PER IL MODEL PART
                    nodes_for_element.append(node_id)   # nodes useful to create a new element
                else:   # SE NEL MODEL PART IL NODO NON C'È, VUOL DIRE CHE È UN NODO SULLA CIRCONFERENZA. LO AGGIUNGO NEL MODEL PART
                    coords = triangle_dict_bottom["vertices"][node]   # LEGGO LE COORDINATE DEL NODO NELLA LISTA IN vertices
                    self.ModelPart.CreateNewNode(max_node_id, coords[0], coords[1], z_min)  # SE IL NODO NON È NEL MODELPART VUOL DIRE CHE È SULLA CIRCONFERENZA E QUINDI HA LA Z PARI A z_min
                    nodes_for_element.append(max_node_id)   # nodes useful to create a new element
                    coord_3D_to_mdpa[node] = max_node_id
                    mdpa_to_coord_2D[max_node_id] = node
                    max_node_id += 1    # we update the Id
            
            self.ModelPart.CreateNewElement("Element2D3N", max_elem_id, nodes_for_element, properties)
            max_elem_id += 1

        all_facets = list()         # list with all node ids of faces
        all_markers = list()        # list of integers [1 = bottom, 2 = topper and 3 = lateral]

        # Bottom: fill in the "facets" and "markers" lists
        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()
            faces = list()
            for node in nodes:
                faces.append(mdpa_to_coord_2D[node.Id])    # LEGGO IL CORRISPONDETE node.Id DEL MODEL PART IN all_points_3D
            
            all_facets.append(faces)
            all_markers.append(1)       # marker: 1 bottom; 2 top; 3 lateral
        
        # Bottom: number of bottom nodes
        n_bottom_nodes = len(all_points_3D)     # so far there are only bottom nodes

        # Top: triangle
        volume_height = z_max + h_value
        B = dict(vertices=coord_2D_top)
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

        # Lateral: quadrilateral elements
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


    def MeshHole(self, building_model_part):
        # bisogna settare il terreno in mesher e passare il model part con gli edifici
        # terrain_model_part = self.ModelPart
        # building_model_part: ogni edificio è in un sub model part
        # funzione basata su ShiftBuildingOnTerrain in geo_building
        
        import numpy as np

        coords = KratosMultiphysics.Array3()

        # we set terrain_model_part in BuildingUtilities
        BU = KratosGeo.BuildingUtilities(self.ModelPart)

        for current_building in building_model_part.SubModelParts:
            print("[DEBUG][geo_building][ShiftBuildingOnTerrain] Building name: ", current_building.Name)
            # each Building is in a different sub_model_part
            
            displacement = []        # the vector where we save the Z coordinate to evaluate the minimum for each building
            shift = 0.0
            for node_building in current_building.Nodes:
                # take only nodes with z = 0.0 which are the nodes at the base of the Building
                if node_building.Z != 0.0:
                    continue

                # fill coords array
                coords[0] = node_building.X
                coords[1] = node_building.Y
                coords[2] = node_building.Z

                # return the ID of the element containing the located point (from coords)
                # if (id_elem = 0) the point is outside from the elements
                id_elem = BU.CheckIfInternal(coords)

                # we check if all the nodes are inside terrain_model_part. If at least one node is outside, the whole building must be outside (and deleted)
                if (id_elem == 0):
                    break

                # we extract the element with the id_elem
                pelem = self.ModelPart.GetElement(id_elem)
                
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

                print(Psi)
                
                displacement.append(Psi[2])                # append just coordinate Z

            else:
                # we check if "displacement" is empty
                if displacement:
                    # shift = min(displacement)
                    shift = max(displacement)

            for node in current_building.Nodes:
                if (node.Z > 0.0):  # extrusion
                    node.Z += shift

        """ CHECK THIS """
        self.HasBuildingHull = True


    def ForoTerreno(self, building_model_part, h_value=0.0, circ_division=60, num_sector=12):
        # terrain_model_part = self.ModelPart
        # building_model_part: ogni edificio è in un sub model part e sono già posizionati correttamente sul terreno

        # ### reading points from model part
        # X = []; Y = []      # lists are generated to use predefined operations
        # for node in self.ModelPart.Nodes:
        #     node.X = self.truncate(node.X)       # we truncate the coordinate to the fifth decimal digit
        #     node.Y = self.truncate(node.Y)       # we truncate the coordinate to the fifth decimal digit
        #     X.append(node.X)
        #     Y.append(node.Y)
        
        # # calculating minimum and maximum value of X, Y and Z coordinates of the entire initial domain
        # x_min = min(X);    x_max = max(X)
        # y_min = min(Y);    y_max = max(Y)

        # # the centre of the circle
        # self.x_center = (x_min + x_max)/2.0
        # self.y_center = (y_min + y_max)/2.0

        # buffering zone
        delta = self.r_boundary/20         # evaluate this value!!!
        print("[DEBUG][ForoTerreno] delta: ", delta)

        # we set all nodes and elements as NOT to erase
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,False)
        for elem in self.ModelPart.Elements:
            elem.Set(KratosMultiphysics.TO_ERASE,False)

        Z = []

        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()     # list of nodes
            for node in nodes:
                x_curr = node.X         # current X coordinate
                y_curr = node.Y         # current Y coordinate

                dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))     # calculate the distance between current node and center of circle

                if (dist > self.r_boundary-delta):
                    node.Set(KratosMultiphysics.TO_ERASE,True)
                    elem.Set(KratosMultiphysics.TO_ERASE,True)
                else:
                    # we fill the Z list to compute the min/max with the only nodes inside the domain
                    Z.append(node.Z)        # this list is useful to compute the min(Z)/max(Z) after
        
        # we delete nodes and element marked as TO_ERASE
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.ModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)

        # we update z_min and z_max
        z_min = min(Z) if Z else 0.0    # z_min = 0.0 if Z is empty
        z_max = max(Z) if Z else 0.0    # z_max = 0.0 if Z is empty


        """ SMOOTHING PROCEDURE """
        for node in self.ModelPart.Nodes:
            x_curr = node.X
            y_curr = node.Y
            dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))
            if (dist > self.r_ground):
                # node_Z = (node_Z - node_Z_min) * Z_beta + node_Z_min
                Z_beta = (-(dist-self.r_ground) / (self.r_boundary-self.r_ground))+1       # we calculate beta
                self.ModelPart.GetNode(node.Id).Z = (node.Z - z_min) * Z_beta + z_min

        # this step is useful to the n sectors
        circ_division = num_sector * round(circ_division / num_sector)
        # a list with the division points of the circumference
        theta = self._custom_range(0.0, 2*math.pi, 2*math.pi/circ_division)  # circ_division is the number of division of the 2*pi

        # lists with nodes on circle boundary
        x_circle = []; y_circle = []
        for th in theta[:-1]:   # the first and last values are the same
            x = self.r_boundary * math.cos(th) + self.x_center      # we move the node along X to consider that the circle have not the centre in (0.0, 0.0)
            y = self.r_boundary * math.sin(th) + self.y_center      # we move the node along Y to consider that the circle have not the centre in (0.0, 0.0)
            x_circle.append(self.truncate(x))    # we add truncated coordinates
            y_circle.append(self.truncate(y))    # we add truncated coordinates


        """ TRIANGLE """
        # reference: https://rufat.be/triangle/examples.html
        coord_2D = []
        coord_2D_top = []   # coordinate 2D dei nodi sulla circonferenza (utile per creare la parte superiore del dominio)
        coord_3D = []       # lista con le coordinate 3D di tutti i nodi
        segments = []       # list of segments delimiting triangulations (bottom)

        all_facets = []     # list with all node ids of faces
        all_markers = []    # list of integers [1 = bottom, 2 = topper and 3 = lateral]

        len_circ = len(x_circle)    # x_circle and y_circle have the same number of nodes

        # circumference nodes and segments on it
        for id, (x_value, y_value) in enumerate(zip(x_circle, y_circle)):
            coord_2D.append([x_value, y_value])         # Bottom
            coord_2D_top.append([x_value, y_value])     # Top
            coord_3D.append([x_value, y_value, z_min])
            
            if (id == len_circ-1):     # the last one
                segments.append([id, 0])    # the last value with the first to close the circle
                continue
            segments.append([id, id+1])     # the i-th segment

        # aggiungo tutti i nodi del Model Part in coord_2D e coord_3D
        for node in self.ModelPart.Nodes:
            coord_2D.append([node.X, node.Y])
            coord_3D.append([node.X, node.Y, node.Z])

        # we delete buildings out of r_buildings
        building = GeoBuilding()
        building.SetGeoModelPart(building_model_part)
        building.DeleteBuildingsOutsideBoundary(x_center=self.x_center,
                                                y_center=self.y_center,
                                                radius=self.r_buildings)

        # we estract the center_model_part
        self.extract_center(self.ModelPart)
        # lists with the coordinates of base and top of the buildings (in 3D)
        coord_3D_base_building, coord_3D_top_building = building.FindDistanceFromTerrain(self.center_model_part)

        len_only_terrain = len(coord_2D)    # terrain nodes only
        n_nodes_per_building = []           # number of nodes for each building (base only)

        # reference: https://stackoverflow.com/questions/48462044/shapely-parallel-offset-sometimes-does-not-generate-closed-ring
        from shapely.geometry.polygon import LinearRing

        holes = []      # ogni foro rappresenta un edificio sul terreno
        # print("NumberOfSubModelParts: ", building_model_part.NumberOfSubModelParts()); input()
        print(10*"1")   # 1111111111
        for current_sub_model in building_model_part.SubModelParts:
            temp_segment = []   # TODO: delete it
            
            n_nodes_per_building.append(int(current_sub_model.NumberOfNodes() / 2))     # we only need the base nodes
            first_id_building = len(coord_2D)   # numero dei nodi fino ad ora (utile per l'id dei nodi per il calcolo dei segmenti)
            coord_2D_building = []              # coordinate 2D dell'edificio corrente

            # CREO UNA NUOVA LISTA CON I SOLI NODI A QUOTA Z=0.0
            # PERCHÉ, ALTRIMENTI, NON SAPREI QUANDO STO PROCESSANDO L'ULTIMO NODO
            # (UTILE PER POTER CHIUDERE IL SEGMENTO IN Triangle)
            nodes_list = []     # lista con i nodi dell'edificio a quota Z=0.0
            for node in current_sub_model.Nodes:
                if (node.Z >= 1e-5):    # se la coordinata Z è diversa da 0.0 vado avanti
                    continue
                nodes_list.append((node.X, node.Y, node.Z))
            
            for id, node in enumerate(nodes_list):
                coord_2D_building.append([node[0], node[1]])

                if (id == len(nodes_list)-1):
                    segments.append([first_id_building+id, first_id_building])
                    temp_segment.append([id, 0])    # TODO: delete it
                    continue
                segments.append([first_id_building+id, first_id_building+id+1])
                temp_segment.append([id, id+1])

            # """
            # NON POSSO SKIPPARE IN QUESTO MODO SE coord_2D_building È VUOTO
            # """
            # if not coord_2D_building:
            #     n_nodes_per_building = n_nodes_per_building[:-1]
            #     print("*********** coord_2D_building è vuoto.\nPREMERE UN TASTO PER CONTINUARE"); input()
            #     continue    # if coord_2D_building is empty

            coord_2D.extend(coord_2D_building)      # aggiungiamo coord_2D_building in coord_2D

            poly_line = LinearRing(coord_2D_building)

            # reference: https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
            if (self._towards(coord_2D_building) == "clockwise"):
                side = "right"
            else:
                side = "left"
            
            poly_line_offset = poly_line.parallel_offset(0.05, side=side, resolution=16, 
                                                         join_style=2, mitre_limit=1)
            # print()
            # print("poly_line: ", poly_line)
            # print("type poly_line_offset: ", poly_line_offset.type)
            # print("type poly_line: ", poly_line.type)
            # print("poly_line_offset: ", poly_line_offset)
            # print("coord_2D_building: ", coord_2D_building)
            # # print("poly_line_offset.xy: ", poly_line_offset.xy)
            # print()
            # # input(poly_line_offset.type)
            if (poly_line_offset.type == "MultiLineString"):
                poly_line_offset = poly_line_offset[0]
                print("poly_line_offset new: ", poly_line_offset)
                input("continua...")
            x_hole = poly_line_offset.xy[0][0]
            y_hole = poly_line_offset.xy[1][0]

            holes.append([x_hole, y_hole])

            # visualizzo la geometria dell'edificio da 420m in poi (con r_building pari a 450m -> 450*0.93=418.5m)
            for node in current_sub_model.Nodes:
                dist = math.sqrt(((self.x_center-node.X)**2)+((self.y_center-node.Y)**2))
                # if (dist >= self.r_buildings*0.93):
                #     print("dist: ", dist)
                #     temp = dict(vertices=coord_2D_building, segments=temp_segment)
                #     triangle_temp = triangle.triangulate(temp, "p")
                    
                #     import matplotlib.pyplot as plt
                #     triangle.compare(plt, temp, triangle_temp)
                #     # manager = plt.get_current_fig_manager()
                #     # manager.full_screen_toggle()
                #     plt.show()
                break
        print(10*"2")   # 2222222222
        
        # print("n_nodes_per_building: ", n_nodes_per_building); input()
        # print("holes: ", holes); input()

        # print("\ncoord_2D: ", coord_2D)
        # print("\nsegments: ", segments)
        # print("\nholes:", holes)
        
        # Bottom: triangle
        C = dict(vertices=coord_2D, segments=segments, holes=holes)
        triangle_dict_bottom = triangle.triangulate(C, "p")     # triangle_dict_bottom {vertices: [...], triangles: [...], vertex_markers: [...]}
        print(10*"3")   # 3333333333

        # import matplotlib.pyplot as plt
        # triangle.compare(plt, C, triangle_dict_bottom)
        # plt.show()

        # aggiungo prima tutti i nodi della base degli edifici e poi tutti i nodi della parte superiore degli edifici
        coord_3D.extend(coord_3D_base_building)
        coord_3D.extend(coord_3D_top_building)
        print(10*"4")   # 4444444444

        # Bottom: number of bottom nodes
        n_bottom_nodes = len(coord_3D)      # so far there are only bottom nodes
        print(10*"5")   # 5555555555

        # Bottom: fill in the "facets" and "markers" lists
        for face in triangle_dict_bottom["triangles"]:
            all_facets.append([face[0], face[1], face[2]])
            all_markers.append(1)
        print(10*"6")   # 6666666666


        """
        in n_nodes_per_building ci sono i numeri dei nodi per ogni edificio
        
        ESEMPIO DI OUTPUT
        n_nodes_per_building:  [4, 4, 12, 4, 4, 12, 4, 8, 8, 4, 6, 8, 4, 11, 36,
                                4, 4, 13, 4, 7, 4, 4, 6, 4, 4, 4, 7, 6, 4, 4, 4,
                                4, 4, 4, 4, 4, 4, 4, 4, 6, 4, 4, 6, 4, 4, 4, 4, 4]
        """

        # Bottom: fill "facets" of buildings ("markers" are still 1)
        # conosco il numero di nodi per ogni edificio;
        # quindi sposto il "puntatore" dei nodi già processati.
        # di volta in volta aggiungo il numero di nodi dell'edificio appena processato
        target_down = len_only_terrain  # aggiungerò il numero dei nodi degli edifici già processati
        target_up = len(coord_2D)       # target dei nodi superiori dell'edificio. all'inizio è pari alla lunghezza di coord_2D (cioè tutti i nodi del terreno e della base degli edifici)
        for n_nodes in n_nodes_per_building:
            coord_2D_top_building = []  # coordinate dei nodi superiori dell'edificio che stiamo processando
            segments_top_building = []  # segmenti sulla parte superiore dell'edificio corrente
            for n in range(n_nodes):
                # 
                # n3         n6 ------- n5
                # |`\         `\  face  | 
                # |  `\         `\  2   | 
                # |    `\         `\    | 
                # | face `\         `\  | 
                # |   1    `\         `\| 
                # n1 ------- n2         n4

                # node ids face lateral 1
                n1 = target_down + n
                n2 = target_down if (n == n_nodes-1) else target_down + n + 1 
                n3 = target_up + n   # target_up è il numero dei nodi del terreno e base edifici. dopo target_up iniziano i nodi nella parte superiore dell'edificio che stiamo processando
                all_facets.append([n1, n2, n3])
                all_markers.append(1)

                # node ids face lateral 2
                n4 = n2
                n5 = target_up if (n == n_nodes-1) else target_up + n + 1
                n6 = n3
                all_facets.append([n4, n5, n6])
                all_markers.append(1)

                # in n3 c'è l'id del nodo superiore. alla fine del ciclo avrò tutti i nodi superiori dell'edificio corrente
                coord_2D_top_building.append([coord_3D[n3][0], coord_3D[n3][1]])
                if (n == n_nodes-1):
                    segments_top_building.append([n, 0])
                    continue
                segments_top_building.append([n, n+1])
            
            dict_top_building = dict(vertices=coord_2D_top_building,
                                     segments=segments_top_building)
            triangle_dict_top_building = triangle.triangulate(dict_top_building, "p")

            # triangle.compare(plt, dict_top_building, triangle_dict_top_building)
            # plt.show()

            for node in triangle_dict_top_building["triangles"]:
                all_facets.append([ node[0]+target_up,
                                    node[1]+target_up,
                                    node[2]+target_up])
                all_markers.append(1)

            # aggiorno i nodi già processati
            target_down += n_nodes
            target_up += n_nodes

        print(10*"7")   # 7777777777


        # Top: triangle
        volume_height = z_max + h_value
        D = dict(vertices=coord_2D_top)
        triangle_dict_top = triangle.triangulate(D)     # triangle_dict_top {vertices: [...], triangles: [...], vertex_markers: [...]}

        # Top: updated coord_3D with top nodes
        for node in triangle_dict_top["vertices"]:
            coord_3D.append((node[0], node[1], volume_height))

        # Top: fill in the "facets" and "markers" lists
        for face in triangle_dict_top["triangles"]:
            all_facets.append([ face[0]+n_bottom_nodes,
                                face[1]+n_bottom_nodes,
                                face[2]+n_bottom_nodes])        # list of lists. "...+n_bottom_nodes"
            all_markers.append(2)

        # Lateral: quadrilateral elements
        for n, coord in enumerate(coord_2D_top[:-1]):
            all_facets.append([n, n+1, n+n_bottom_nodes+1, n+n_bottom_nodes])
            all_markers.append(3)
        # the last value
        all_facets.append([n+1, 0, n_bottom_nodes, n+n_bottom_nodes+1])
        all_markers.append(3)

        print(10*"8")   # 8888888888

        """ TETGEN (VIA MESHPY) """
        ### using Meshpy for the creation of an initial mesh https://documen.tician.de/meshpy/
        mesh_info = MeshInfo()
        mesh_info.set_points(coord_3D)
        mesh_info.set_facets(all_facets, markers=all_markers)
        print(10*"9")   # 9999999999

        # print("coord_3D: ", coord_3D)
        print("len coord_3D: ", len(coord_3D))
        # print("all_facets: ", all_facets)
        print("len all_facets: ", len(all_facets))
        # print("all_markers: ", all_markers)
        print("len all_markers: ", len(all_markers))
        print("mesh_info: ", mesh_info)

        # build the mesh
        mesh = build(mesh_info)
        print(10*"10")  #10101010101010101010
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


    def MeshSoloForo(self, building_model_part, dict_ext_edges=None, circ_division=60, num_sector=12):
        """
            [ok] viene estratta una porzione circolare;
            [ok] viene eseguita la procedura di smoothing;
            [ok] viene creato un sub model part "bottom_model_part" in cui c'è il terreno meno l'impronta degli edifici;
            [ok] viene creato un sub model part "top_model_part" in cui c'è il tappo superiore del dominio;
            [ok] viene creato un sub model part "lateral_model_part" in cui ci sono le superfici laterali del dominio;
            [ok] viene creato un dizionario in cui si sono le informazioni sui nodi degli edifici (utile per posizionare gli edifici successivamente)
        """
        # terrain_model_part = self.ModelPart
        # building_model_part: ogni edificio è in un sub model part e sono già posizionati correttamente sul terreno

        # buffering zone
        # delta = self.r_boundary/20        # evaluate this value!!!
        delta = self.r_boundary/100         # evaluate this value!!!
        print("[DEBUG][ForoTerreno] delta: ", delta)

        # we set all nodes and elements as NOT to erase
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,False)
        for elem in self.ModelPart.Elements:
            elem.Set(KratosMultiphysics.TO_ERASE,False)

        Z = []

        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()     # list of nodes
            for node in nodes:
                x_curr = node.X         # current X coordinate
                y_curr = node.Y         # current Y coordinate

                dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))     # calculate the distance between current node and center of circle

                if (dist > self.r_boundary-delta):
                    node.Set(KratosMultiphysics.TO_ERASE,True)
                    elem.Set(KratosMultiphysics.TO_ERASE,True)
                else:
                    # we fill the Z list to compute the min/max with the only nodes inside the domain
                    Z.append(node.Z)        # this list is useful to compute the min(Z)/max(Z) after
        
        # we delete nodes and element marked as TO_ERASE
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.ModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)

        # we update z_min and z_max
        z_min = min(Z) if Z else 0.0    # z_min = 0.0 if Z is empty
        z_max = max(Z) if Z else 0.0    # z_max = 0.0 if Z is empty


        """ SMOOTHING PROCEDURE """
        for node in self.ModelPart.Nodes:
            x_curr = node.X
            y_curr = node.Y
            dist = math.sqrt(((self.x_center-x_curr)**2)+((self.y_center-y_curr)**2))
            if (dist > self.r_ground):
                # node_Z = (node_Z - node_Z_min) * Z_beta + node_Z_min
                Z_beta = (-(dist-self.r_ground) / (self.r_boundary-self.r_ground))+1       # we calculate beta
                self.ModelPart.GetNode(node.Id).Z = (node.Z - z_min) * Z_beta + z_min

        # this step is useful to the n sectors
        circ_division = num_sector * round(circ_division / num_sector)
        # a list with the division points of the circumference
        theta = self._custom_range(0.0, 2*math.pi, 2*math.pi/circ_division)  # circ_division is the number of division of the 2*pi

        # lists with nodes on circle boundary
        x_circle = []; y_circle = []
        for th in theta[:-1]:   # the first and last values are the same
            x = self.r_boundary * math.cos(th) + self.x_center      # we move the node along X to consider that the circle have not the centre in (0.0, 0.0)
            y = self.r_boundary * math.sin(th) + self.y_center      # we move the node along Y to consider that the circle have not the centre in (0.0, 0.0)
            x_circle.append(self.truncate(x))    # we add truncated coordinates
            y_circle.append(self.truncate(y))    # we add truncated coordinates


        """ TRIANGLE """
        # reference: https://rufat.be/triangle/examples.html
        coord_2D = []
        coord_2D_top = []   # coordinate 2D dei nodi sulla circonferenza (utile per creare la parte superiore del dominio)
        coord_3D = []       # lista con le coordinate 3D di tutti i nodi
        segments = []       # list of segments delimiting triangulations (bottom)

        # TODO: posso creare solo coord_3D e poi estrarre coord_2D per usarlo in "triangle"
        #       coord_2D = [[x,y] for x,y,z in coord_3D]

        len_circ = len(x_circle)    # x_circle and y_circle have the same number of nodes

        # circumference nodes and segments on it
        for id, (x_value, y_value) in enumerate(zip(x_circle, y_circle)):
            coord_2D.append([x_value, y_value])         # Bottom
            coord_2D_top.append([x_value, y_value])     # Top
            coord_3D.append([x_value, y_value, z_min])
            
            if (id == len_circ-1):     # the last one
                segments.append([id, 0])    # the last value with the first to close the circle
                continue
            segments.append([id, id+1])     # the i-th segment

        # aggiungo tutti i nodi del Model Part in coord_2D e coord_3D
        for node in self.ModelPart.Nodes:
            coord_2D.append([node.X, node.Y])
            coord_3D.append([node.X, node.Y, node.Z])


        # reference: https://stackoverflow.com/questions/48462044/shapely-parallel-offset-sometimes-does-not-generate-closed-ring
        from shapely.geometry.polygon import LinearRing

        holes = []      # ogni foro rappresenta un edificio sul terreno
        coord_base_building = []
        n_nodes_per_building = []           # number of nodes for each building (base only)

        self.dict_buildingMP_to_terrainMP = {}  # key: node.Id in building Model Part; value: node.Id in terrain Model Part


        print(10*"1")   # 1111111111
        for current_sub_model in building_model_part.SubModelParts:
            
            n_nodes_per_building.append(int(current_sub_model.NumberOfNodes() / 2))     # we only need the base nodes
            first_id_building = len(coord_2D)   # numero dei nodi fino ad ora (utile per l'id dei nodi per il calcolo dei segmenti)
            coord_2D_building = []              # coordinate 2D dell'edificio corrente
            # print("prima first_id_building: ", first_id_building)

            if (dict_ext_edges is None):
                # aggiungo solo i nodi alla base in coord_2D_building e creo i segmenti (perimetro dell'edificio)
                for pos, node in enumerate(current_sub_model.Nodes):    # posizione e nodo nell'i-esimo edificio
                    if (pos < n_nodes_per_building[-1]):
                        # nodi alla base dell'edificio
                        self.dict_buildingMP_to_terrainMP[node.Id] = first_id_building+pos+1    # mappo l'id che il nodo ha nel MP degli edifici con il nodo che avrà nel MP del terreno
                        # print("node.Id: ", node.Id, " ---> first_id_building+pos: ", first_id_building+pos)
                        coord_base_building.append([node.X, node.Y, node.Z])
                        coord_2D_building.append([node.X, node.Y])
                        if (pos == n_nodes_per_building[-1] - 1):   # l'ultimo nodo alla base dell'i-esimo edificio
                            segments.append([first_id_building+pos, first_id_building])
                            continue
                        segments.append([first_id_building+pos, first_id_building+pos+1])
            else:
                # we get the boundary edges and we sort them
                edges = dict_ext_edges[current_sub_model.Name]
                edges = self._sort_edges(edges)
                
                for pos, nodes in enumerate(edges):     # edges = [(1,2), (2,3), ...]
                    self.dict_buildingMP_to_terrainMP[node.Id] = first_id_building+pos    # TODO: ho tolto il +1.     mappo l'id che il nodo ha nel MP degli edifici con il nodo che ha nel MP del terreno (il +1 finale è perché la lista parte da 0 mentre il MP parte da 1)
                    node = current_sub_model.GetNode(nodes[0])
                    coord_base_building.append([node.X, node.Y, node.Z])
                    coord_2D_building.append([node.X, node.Y])
                    if (pos == len(edges)-1):
                        segments.append([first_id_building+pos, first_id_building])
                        continue
                    segments.append([first_id_building+pos, first_id_building+pos+1])

            
            coord_2D.extend(coord_2D_building)      # aggiungiamo coord_2D_building in coord_2D
            # print("len(coord_2D_building): ", len(coord_2D_building))
            # print("dopo len(coord_2D): ", len(coord_2D))
            # input("PAUSE")

            poly_line = LinearRing(coord_2D_building)

            # reference: https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order
            if (self._towards(coord_2D_building) == "clockwise"):
                side = "right"
            else:
                side = "left"
            
            poly_line_offset = poly_line.parallel_offset(0.05, side=side, resolution=16, 
                                                         join_style=2, mitre_limit=1)
            
            # if poly_line_offset is "MultiLineString" we get the first value
            if (poly_line_offset.type == "MultiLineString"):
                poly_line_offset = poly_line_offset[0]

            # hole coordinates
            x_hole = poly_line_offset.xy[0][0]
            y_hole = poly_line_offset.xy[1][0]

            holes.append([x_hole, y_hole])

        print(10*"2")   # 2222222222
        
        
        # Bottom: triangle
        C = dict(vertices=coord_2D, segments=segments, holes=holes)
        triangle_dict_bottom = triangle.triangulate(C, "p")     # triangle_dict_bottom {vertices: [...], triangles: [...], vertex_markers: [...]}
        print(10*"3")   # 3333333333


        # aggiungo tutti i nodi della base degli edifici
        coord_3D.extend(coord_base_building)
        

        #######################################################
        # ripulisco il model part
        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()

        bottom_model_part = self.ModelPart.CreateSubModelPart("BottomModelPart")
        top_model_part = self.ModelPart.CreateSubModelPart("TopModelPart")
        lateral_model_part = self.ModelPart.CreateSubModelPart("LateralModelPart")
        properties = self.ModelPart.Properties[0]

        bottom_cond = [];   bottom_nodes = []
        top_cond = [];      top_nodes = []
        lateral_cond = [];  lateral_nodes = []

        ### Bottom Nodes
        n_start_nodes = self._find_max_node_id() + 1    # at the beginning _find_max_node_id() return 0 because the Model Part is empty
        id_node = n_start_nodes
        for node in coord_3D:
            # node id is shifted up by 1
            self.ModelPart.CreateNewNode(id_node, node[0], node[1], node[2])
            bottom_nodes.append(id_node)
            id_node += 1
        
        ### Bottom Conditions
        id_cond = self._find_max_cond_id() + 1
        for nodes in triangle_dict_bottom["triangles"]:
            self.ModelPart.CreateNewCondition(  "WallCondition3D3N",
                                                id_cond,
                                                [nodes[0]+n_start_nodes,
                                                nodes[1]+n_start_nodes,
                                                nodes[2]+n_start_nodes],
                                                properties)
            bottom_cond.append(id_cond)
            id_cond += 1
        
        # Bottom: add Nodes and Conditions to bottom_model_part
        bottom_model_part.AddNodes(bottom_nodes)
        bottom_model_part.AddConditions(bottom_cond)


        # Top: triangle
        volume_height = z_max + self.height
        D = dict(vertices=coord_2D_top)
        triangle_dict_top = triangle.triangulate(D)     # triangle_dict_top {vertices: [...], triangles: [...], vertex_markers: [...]}

        ### Top Nodes
        n_bottom_nodes = self._find_max_node_id() + 1   # utile per la creazione delle condizioni
        id_node = n_bottom_nodes
        for node in triangle_dict_top["vertices"]:
            self.ModelPart.CreateNewNode(id_node, node[0], node[1], volume_height)
            top_nodes.append(id_node)
            id_node += 1

        ### Top Conditions
        # id_cond = self._find_max_cond_id() + 1
        for nodes in triangle_dict_top["triangles"]:
            self.ModelPart.CreateNewCondition(  "WallCondition3D3N",
                                                id_cond,
                                                [nodes[0]+n_bottom_nodes,
                                                nodes[1]+n_bottom_nodes,
                                                nodes[2]+n_bottom_nodes],
                                                properties)
            top_cond.append(id_cond)
            id_cond += 1

        # Top: add Nodes and Conditions to top_model_part
        top_model_part.AddNodes(top_nodes)
        top_model_part.AddConditions(top_cond)


        # Lateral Conditions
        len_coord_2D_top = len(coord_2D_top)
        for i in range(len_coord_2D_top):
            # 
            # n3         n6 ------- n5
            # |`\         `\  face  | 
            # |  `\         `\  2   | 
            # |    `\         `\    | 
            # | face `\         `\  | 
            # |   1    `\         `\| 
            # n1 ------- n2         n4
            
            # node ids face 1
            n1 = n_start_nodes + i
            n2 = n_start_nodes if (i == len_coord_2D_top-1) else n_start_nodes + i + 1
            n3 = n_bottom_nodes + i
            lateral_nodes.extend([n1, n2, n3])

            self.ModelPart.CreateNewCondition(  "WallCondition3D3N",
                                                id_cond,
                                                [n1, n2, n3],
                                                properties)
            lateral_cond.append(id_cond)
            id_cond += 1

            # node ids face 2
            n4 = n2
            n5 = n_bottom_nodes if (i == len_coord_2D_top-1) else n_bottom_nodes + i + 1
            n6 = n3
            lateral_nodes.extend([n4, n5, n6])

            self.ModelPart.CreateNewCondition(  "WallCondition3D3N",
                                                id_cond,
                                                [n4, n5, n6],
                                                properties)
            lateral_cond.append(id_cond)
            id_cond += 1

        # Lateral: add Nodes and Conditions to lateral_model_part
        lateral_model_part.AddNodes(lateral_nodes)
        lateral_model_part.AddConditions(lateral_cond)

        # set all nodes as to erase
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,True)
        
        # set as NOT to erase Nodes belonging to Conditions and Elements
        for cond in self.ModelPart.Conditions:
            nodes = cond.GetNodes()
            for node in nodes:
                node.Set(KratosMultiphysics.TO_ERASE,False)
        for elem in self.ModelPart.Elements:
            nodes = elem.GetNodes()
            for node in nodes:
                node.Set(KratosMultiphysics.TO_ERASE,False)

        # delete all free Nodes
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)


    """" FUNZIONE CREATA PRIMA DI RAFFITTIRE LA MESH DEGLI EDIFICI """
    # def AddBuildingsOnTerrain(self, building_model_part):
    #     """ aggiungo gli edifici al model part del terreno """

    #     if not self.ModelPart.HasSubModelPart("BottomModelPart"):
    #         KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "\"BottomModelPart\" not found!")
    #         return
        
    #     if not self.dict_buildings_on_terrain:
    #         KratosMultiphysics.Logger.PrintWarning("GeoMesher", "\"dict_buildings_on_terrain\" must be set first!")
    #         return
        
    #     skin_isosurface = self.ModelPart.CreateSubModelPart("SKIN_ISOSURFACE")  # TODO: VALUTARE QUESTO NOME
    #     properties = self.ModelPart.Properties[0]

    #     buildings_cond = [];    building_nodes = []

    #     # DEVO CONTARE IL NUMERO DELLE CONDIZIONI E IL NUMERO DI NODI IN self.ModelPart
    #     # CREO LE NUOVE CONDIZIONI CON LE INFO DAGLI EDIFICI
    #     # PER I NODI, DEVO SOSTITUIRE GLI ID ALLA BASE CON QUELLI NEL DIZIONARIO; MENTRE PER I
    #     #   NODI DI COPERTURA, DEVO CREARNE DI NUOVI (SOMMARE IL NUMERO DEGLI ID GIÀ PRESENTI)

    #     # current Id (for new Nodes). we add the number of Nodes in building_model_part to avoid creating Nodes with already existing Ids
    #     curr_node_id = self._find_max_node_id() + building_model_part.NumberOfNodes() + 1
    #     # current Id (for new Conditions). we add the number of Conditions in building_model_part to avoid creating Conditions with already existing Ids
    #     curr_cond_id = self._find_max_cond_id() + building_model_part.NumberOfConditions() + 1

    #     # gli Id dei nodi alla base sono sostituiti con gli Id nel dizionario (perché già esistenti nel model part del terreno);
    #     # gli Id della copertura degli edifici vengono aggiunti al sub model part SKIN_ISOSURFACE e, contemporaneamente,
    #     #   viene aggiornato il sub model part degli edifici per poter estrarre le facce dell'edificio con gli Id dei nodi già aggiornati
    #     # aggiungo gli Id dei nodi al sub model part SKIN_ISOSURFACE e aggiungo anche le facce (che mano a mano sono create processando il model part degli edifici)
    #     for current_sub_model in building_model_part.SubModelParts:
    #         ### Nodes
    #         for pos, node in enumerate(current_sub_model.Nodes):    # posizione e nodo nell'i-esimo edificio
    #             if (pos < current_sub_model.NumberOfNodes()/2):     # processo i nodi alla base
    #                 # nodi alla base dell'edificio
    #                 node.Id = self.dict_buildings_on_terrain[current_sub_model.Name][pos]   # aggiorno l'Id del nodo leggendolo dal dizionario (perché già presente nel model part del terreno)
    #                 building_nodes.append(node.Id)  # aggiungo il nodo con l'Id appena aggiornato
    #             else:
    #                 # nodi della copertura dell'edificio
    #                 self.ModelPart.CreateNewNode(curr_node_id, node.X, node.Y, node.Z)   # creo un nuovo nodo nel model part del terreno
    #                 node.Id = curr_node_id          # aggiorno l'Id per poter estrarre le facce dell'edificio con gli Id già sistemati
    #                 building_nodes.append(node.Id)  # aggiungo il nodo con l'Id appena aggiornato
    #                 curr_node_id += 1               # aggiorno l'Id corrente

    #         ### Conditions
    #         # in building_model_part sono elementi ma qui devono essere delle condizioni
    #         for elem in current_sub_model.Elements:
    #             nodes = elem.GetNodes()
    #             self.ModelPart.CreateNewCondition(  "WallCondition3D3N",
    #                                                 curr_cond_id,
    #                                                 [nodes[0].Id, nodes[1].Id, nodes[2].Id],
    #                                                 properties)
    #             buildings_cond.append(curr_cond_id)
    #             curr_cond_id += 1

    #         # add Nodes and Conditions to bottom_model_part
    #         skin_isosurface.AddNodes(building_nodes)
    #         skin_isosurface.AddConditions(buildings_cond)


    # def AddBuildingsOnTerrain(self, building_model_part):
    #     """ aggiungo gli edifici al model part del terreno """

    #     if not self.ModelPart.HasSubModelPart("BottomModelPart"):
    #         KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "\"BottomModelPart\" not found!")
    #         return
        
    #     if not self.dict_buildings_on_terrain:
    #         KratosMultiphysics.Logger.PrintWarning("GeoMesher", "\"dict_buildings_on_terrain\" must be set first!")
    #         return
        
    #     skin_isosurface = self.ModelPart.CreateSubModelPart("SKIN_ISOSURFACE")  # TODO: VALUTARE QUESTO NOME
    #     properties = self.ModelPart.Properties[0]

    #     buildings_cond = [];    building_nodes = []

    #     # DEVO CONTARE IL NUMERO DELLE CONDIZIONI E IL NUMERO DI NODI IN self.ModelPart
    #     # CREO LE NUOVE CONDIZIONI CON LE INFO DAGLI EDIFICI
    #     # PER I NODI, DEVO SOSTITUIRE GLI ID ALLA BASE CON QUELLI NEL DIZIONARIO; MENTRE PER I
    #     #   NODI DI COPERTURA, DEVO CREARNE DI NUOVI (SOMMARE IL NUMERO DEGLI ID GIÀ PRESENTI)

    #     # current Id (for new Nodes). we add the number of Nodes in building_model_part to avoid creating Nodes with already existing Ids
    #     curr_node_id = self._find_max_node_id() + building_model_part.NumberOfNodes() + 1
    #     # current Id (for new Conditions). we add the number of Conditions in building_model_part to avoid creating Conditions with already existing Ids
    #     curr_cond_id = self._find_max_cond_id() + building_model_part.NumberOfConditions() + 1

    #     # gli Id dei nodi alla base sono sostituiti con gli Id nel dizionario (perché già esistenti nel model part del terreno);
    #     # gli Id della copertura degli edifici vengono aggiunti al sub model part SKIN_ISOSURFACE e, contemporaneamente,
    #     #   viene aggiornato il sub model part degli edifici per poter estrarre le facce dell'edificio con gli Id dei nodi già aggiornati
    #     # aggiungo gli Id dei nodi al sub model part SKIN_ISOSURFACE e aggiungo anche le facce (che mano a mano sono create processando il model part degli edifici)
    #     for current_sub_model in building_model_part.SubModelParts:
    #         ### Nodes
    #         for node in current_sub_model.Nodes:
    #             if (node.Id in self.dict_buildingMP_to_terrainMP):
    #                 # è un nodo alla base dell'edificio
    #                 node.Id = self.dict_buildingMP_to_terrainMP[node.Id]
    #                 building_nodes.append(node.Id)  # aggiungo il nodo con l'Id appena aggiornato
    #             else:
    #                 # tutti gli altri nodi dell'edificio
    #                 self.ModelPart.CreateNewNode(curr_node_id, node.X, node.Y, node.Z)   # creo un nuovo nodo nel model part del terreno
    #                 node.Id = curr_node_id          # aggiorno l'Id per poter estrarre le facce dell'edificio con gli Id già sistemati
    #                 building_nodes.append(node.Id)  # aggiungo il nodo con l'Id appena aggiornato
    #                 curr_node_id += 1               # aggiorno l'Id corrente

    #         ### Conditions
    #         # in building_model_part sono elementi ma qui devono essere delle condizioni
    #         nodes_visited = {}      # Nodes already visited. key: building node.Id, value: node id in terrain Model Part
    #         for elem in current_sub_model.Elements:
    #             nodes = elem.GetNodes()
    #             # list_nodes = []     # Nodes useful for the new Condition
    #             # for node in nodes:
    #             #     if (node.Id in self.dict_buildingMP_to_terrainMP):
    #             #         n_id = self.dict_buildingMP_to_terrainMP[node.Id]
    #             #         list_nodes.append(n_id)
    #             #         building_nodes.append(n_id)
    #             #     else:
    #             #         if (node.Id in nodes_visited):
    #             #             n_id = nodes_visited[node.Id]
    #             #             list_nodes.append(n_id)
    #             #             building_nodes.append(n_id)
    #             #         else:
    #             #             self.ModelPart.CreateNewNode(curr_node_id, node.X, node.Y, node.Z)   # creo un nuovo nodo nel model part del terreno
    #             #             list_nodes.append(curr_node_id)
    #             #             nodes_visited[node.Id] = curr_node_id
    #             #             building_nodes.append(curr_node_id)
    #             #             curr_node_id += 1


    #             # [list_nodes[0], list_nodes[1], list_nodes[2]],
    #             self.ModelPart.CreateNewCondition(  "WallCondition3D3N",
    #                                                 curr_cond_id,
    #                                                 [nodes[0].Id, nodes[1].Id, nodes[2].Id],
    #                                                 properties)
    #             buildings_cond.append(curr_cond_id)
    #             curr_cond_id += 1

    #         # add Nodes and Conditions to bottom_model_part
    #         skin_isosurface.AddNodes(building_nodes)
    #         skin_isosurface.AddConditions(buildings_cond)

    ### PROVO AD AGGIORNARE GLI ID DEI NODI CREANDO UN DIZIONARIO SENZA MODIFICARE IL MODEL PART DEGLI EDIFICI ###
    def AddBuildingsOnTerrain(self, building_model_part):
        """ aggiungo gli edifici al model part del terreno """

        if not self.ModelPart.HasSubModelPart("BottomModelPart"):
            KratosMultiphysics.Logger.PrintWarning("GeoBuilding", "\"BottomModelPart\" not found!")
            return
        
        skin_isosurface = self.ModelPart.CreateSubModelPart("SKIN_ISOSURFACE")  # TODO: VALUTARE QUESTO NOME
        properties = self.ModelPart.Properties[0]

        buildings_cond = [];    building_nodes = []

        # DEVO CONTARE IL NUMERO DELLE CONDIZIONI E IL NUMERO DI NODI IN self.ModelPart
        # CREO LE NUOVE CONDIZIONI CON LE INFO DAGLI EDIFICI
        # PER I NODI, DEVO SOSTITUIRE GLI ID ALLA BASE CON QUELLI NEL DIZIONARIO; MENTRE PER I
        #   NODI DI COPERTURA, DEVO CREARNE DI NUOVI (SOMMARE IL NUMERO DEGLI ID GIÀ PRESENTI)

        # current Id (for new Nodes). we add the number of Nodes in building_model_part to avoid creating Nodes with already existing Ids
        curr_node_id = self._find_max_node_id() + building_model_part.NumberOfNodes() + 1
        # current Id (for new Conditions). we add the number of Conditions in building_model_part to avoid creating Conditions with already existing Ids
        curr_cond_id = self._find_max_cond_id() + building_model_part.NumberOfConditions() + 1

        # gli Id dei nodi alla base sono sostituiti con gli Id nel dizionario (perché già esistenti nel model part del terreno);
        # gli Id della copertura degli edifici vengono aggiunti al sub model part SKIN_ISOSURFACE e, contemporaneamente,
        #   viene aggiornato il sub model part degli edifici per poter estrarre le facce dell'edificio con gli Id dei nodi già aggiornati
        # aggiungo gli Id dei nodi al sub model part SKIN_ISOSURFACE e aggiungo anche le facce (che mano a mano sono create processando il model part degli edifici)
        nodes_new_orders = {}   # key: old id (id in building Model Part); value: updated id (node id already in terrain Model Part or the new one)
        for current_sub_model in building_model_part.SubModelParts:
            ### Nodes
            for node in current_sub_model.Nodes:
                if (node.Id in self.dict_buildingMP_to_terrainMP):
                    # è un nodo alla base dell'edificio
                    n_id = self.dict_buildingMP_to_terrainMP[node.Id]   # id del nodo nel Model Part del terreno
                    nodes_new_orders[node.Id] = n_id
                    building_nodes.append(n_id)  # aggiungo il nodo con l'Id appena aggiornato
                else:
                    # tutti gli altri nodi dell'edificio
                    self.ModelPart.CreateNewNode(curr_node_id, node.X, node.Y, node.Z)   # creo un nuovo nodo nel model part del terreno
                    nodes_new_orders[node.Id] = curr_node_id
                    building_nodes.append(curr_node_id)  # aggiungo il nodo con l'Id appena aggiornato
                    curr_node_id += 1               # aggiorno l'Id corrente

            ### Conditions
            # in building_model_part sono elementi ma qui devono essere delle condizioni
            for elem in current_sub_model.Elements:
                nodes = elem.GetNodes()
                n1 = nodes_new_orders[nodes[0].Id]
                n2 = nodes_new_orders[nodes[1].Id]
                n3 = nodes_new_orders[nodes[2].Id]
                # print(nodes[0].Id, " ---> ", n1)
                # print(nodes[1].Id, " ---> ", n2)
                # print(nodes[2].Id, " ---> ", n3)
                self.ModelPart.CreateNewCondition(  "WallCondition3D3N",
                                                    curr_cond_id,
                                                    [n1, n2, n3],
                                                    properties)
                buildings_cond.append(curr_cond_id)
                curr_cond_id += 1

            # add Nodes and Conditions to bottom_model_part
            skin_isosurface.AddNodes(building_nodes)
            skin_isosurface.AddConditions(buildings_cond)


    def BuildMesh3D(self, metric_file):
        """ costruiamo mesh volumetrica """

        if not (self.ModelPart.HasSubModelPart("BottomModelPart")):
            KratosMultiphysics.Logger.PrintWarning("GeoMesher", "BottomModelPart not found!")
            return

        if not (self.ModelPart.HasSubModelPart("TopModelPart")):
            KratosMultiphysics.Logger.PrintWarning("GeoMesher", "TopModelPart not found!")
            return

        if not (self.ModelPart.HasSubModelPart("LateralModelPart")):
            KratosMultiphysics.Logger.PrintWarning("GeoMesher", "LateralModelPart not found!")
            return

        if not (self.ModelPart.HasSubModelPart("SKIN_ISOSURFACE")):
            KratosMultiphysics.Logger.PrintWarning("GeoMesher", "SKIN_ISOSURFACE not found!")
            return

        all_nodes = []      # list with all node coordinates
        all_facets = []     # list with all node ids of faces
        all_markers = []    # list of integers [1=bottom, 2=top, 3=lateral, 4=buildings]

        ###########################################################################################
        # INIZIO NUOVO BLOCCO
        # # all Nodes as NOT visited
        # for node in self.ModelPart.Nodes:
        #     node.Set(KratosMultiphysics.VISITED, False)
        
        # # we set all Nodes of buildings as visited
        # for cond in self.ModelPart.GetSubModelPart("SKIN_ISOSURFACE").Conditions:
        #     for node in cond.GetNodes():
        #         node.Set(KratosMultiphysics.VISITED, True)
        
        # metric file
        str_mtr = "{} 1\n".format(self.ModelPart.NumberOfNodes())
        mtr_buildings = 3.0     # metric of the buildings
        mtr_lateral = 20.0      # metric of the lateral surfaces
        mtr_top = 100.0         # metric of the topper
        # min and max values
        mtr_terrain_min = 5.0
        mtr_terrain_max = 20.0
        mtr_middle = 0
        # FINE NUOVO BLOCCO
        ###########################################################################################

       
        # Conditions
        mp_to_list_nodes = {}   # key: node id in Model Part; value: node position in all_nodes list
        for smp in self.ModelPart.SubModelParts:
            markers = 0     # inizializzo il marker con un id non valido
            if (smp.Name == "BottomModelPart"):
                markers = 1
                for node in smp.Nodes:
                    if not (node.Id in mp_to_list_nodes):
                        # we check if the current Node is already in the dictionary
                        all_nodes.append([node.X, node.Y, node.Z])
                        mp_to_list_nodes[node.Id] = len(all_nodes) - 1  # salviamo la posizione del nodo i-esimo nel dizionario
                        dist = math.sqrt(((self.x_center-node.X)**2)+((self.y_center-node.Y)**2))   # calculate the distance between current node and center of circle
                        if (dist <= (self.r_buildings*1.2)):
                            mtr_middle = mtr_terrain_min
                        else:
                            coeff = (mtr_terrain_max-mtr_terrain_min) / (self.r_boundary-(self.r_buildings*1.2))
                            mtr_middle = coeff * (dist - (self.r_buildings*1.2)) + mtr_terrain_min
                        str_mtr += "{}\n".format(mtr_middle)
            elif (smp.Name == "TopModelPart"):
                markers = 2
                for node in smp.Nodes:
                    if not (node.Id in mp_to_list_nodes):
                        # we check if the current Node is already in the dictionary
                        all_nodes.append([node.X, node.Y, node.Z])
                        mp_to_list_nodes[node.Id] = len(all_nodes) - 1  # salviamo la posizione del nodo i-esimo nel dizionario
                        str_mtr += "{}\n".format(mtr_top)
            elif (smp.Name == "LateralModelPart"):
                markers = 3
                for node in smp.Nodes:
                    if not (node.Id in mp_to_list_nodes):
                        # we check if the current Node is already in the dictionary
                        all_nodes.append([node.X, node.Y, node.Z])
                        mp_to_list_nodes[node.Id] = len(all_nodes) - 1  # salviamo la posizione del nodo i-esimo nel dizionario
                        if (node.Z > (self.height*0.9)):    # if Z is on the top
                            str_mtr += "{}\n".format(mtr_top)
                        else:
                            str_mtr += "{}\n".format(mtr_lateral)
            elif (smp.Name == "SKIN_ISOSURFACE"):
                markers = 4
                for node in smp.Nodes:
                    if not (node.Id in mp_to_list_nodes):
                        # we check if the current Node is already in the dictionary
                        all_nodes.append([node.X, node.Y, node.Z])
                        mp_to_list_nodes[node.Id] = len(all_nodes) - 1  # salviamo la posizione del nodo i-esimo nel dizionario
                        str_mtr += "{}\n".format(mtr_buildings)
            else:
                KratosMultiphysics.Logger.PrintWarning("GeoMesher", "Unknown SubModelPart: {}".format(smp.Name))
            
            for cond in smp.Conditions:
                nodes = cond.GetNodes()
                n1 = nodes[0].Id
                n2 = nodes[1].Id
                n3 = nodes[2].Id
                all_facets.append([ mp_to_list_nodes[n1],
                                    mp_to_list_nodes[n2],
                                    mp_to_list_nodes[n3]])
                all_markers.append(markers)


        """ TETGEN (VIA MESHPY) """
        ### using Meshpy for the creation of an initial mesh https://documen.tician.de/meshpy/
        mesh_info = MeshInfo()
        mesh_info.set_points(all_nodes)
        mesh_info.set_facets(all_facets, markers=all_markers)
        print(10*"9")   # 9999999999


        ###########################################################################################
        # INIZIO NUOVO BLOCCO
        with open(metric_file, "w") as fo:
            fo.write(str_mtr)
        
        metric_file = metric_file[:-4]
        mesh_info.load_mtr(metric_file)
        print(10*"A")   # AAAAAAAAAA
        # FINE NUOVO BLOCCO
        ###########################################################################################


        # build the mesh
        # mesh = build(mesh_info)
        # mesh = build(mesh_info, options=Options("qpm"))
        mesh = build(mesh_info, options=Options("qm"))
        print(10*"10")  #10101010101010101010
        self.ModelPart.Nodes.clear()
        self.ModelPart.Elements.clear()
        self.ModelPart.Conditions.clear()

        for smp in self.ModelPart.SubModelParts:
            smp.Nodes.clear()
            smp.Elements.clear()
            smp.Conditions.clear()

        bottom_model_part = self.ModelPart.GetSubModelPart("BottomModelPart")
        top_model_part = self.ModelPart.GetSubModelPart("TopModelPart")
        lateral_model_part = self.ModelPart.GetSubModelPart("LateralModelPart")
        building_model_part = self.ModelPart.GetSubModelPart("SKIN_ISOSURFACE")
        fluid_model_part = self.ModelPart.CreateSubModelPart("Parts_Fluid")
        properties = self.ModelPart.Properties[0]

        bottom_cond   = [];   bottom_nodes = []
        top_cond      = [];      top_nodes = []
        lateral_cond  = [];  lateral_nodes = []
        building_cond = []; building_nodes = []
        fluid_elem    = [];    fluid_nodes = []

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
            elif (marker == 4):
                building_cond.append(j+1)
            else:
                KratosMultiphysics.Logger.PrintWarning("GeoMesher", "Marker {} not valid. 1=bottom, 2=top, 3=lateral, 4=buildings".format(marker))

        bottom_model_part.AddConditions(bottom_cond)
        top_model_part.AddConditions(top_cond)
        lateral_model_part.AddConditions(lateral_cond)
        building_model_part.AddConditions(building_cond)

        for cond in bottom_model_part.Conditions:
            for node in cond.GetNodes():
                bottom_nodes.append(node.Id)
        for cond in top_model_part.Conditions:
            for node in cond.GetNodes():
                top_nodes.append(node.Id)
        for cond in lateral_model_part.Conditions:
            for node in cond.GetNodes():
                lateral_nodes.append(node.Id)
        for cond in building_model_part.Conditions:
            for node in cond.GetNodes():
                building_nodes.append(node.Id)

        bottom_model_part.AddNodes(bottom_nodes)
        top_model_part.AddNodes(top_nodes)
        lateral_model_part.AddNodes(lateral_nodes)
        building_model_part.AddNodes(building_nodes)

        ### Elements
        for k in range(len(mesh.elements)):
            nodes = mesh.elements[k]
            self.ModelPart.CreateNewElement("Element3D4N",
                                            (k+1),
                                            [nodes[0]+1, nodes[1]+1, nodes[2]+1, nodes[3]+1],
                                            properties)
            fluid_nodes.extend([nodes[0]+1, nodes[1]+1, nodes[2]+1, nodes[3]+1])
            fluid_elem.append(k+1)

        fluid_model_part.AddElements(fluid_elem)
        fluid_model_part.AddNodes(fluid_nodes)


    #################################################################
    def Mesh3D_buildings(self):
        """ raffinamento delle superfici degli edifici e restituisce i segmenti di bordo inferiori
            utili per inserire gli edifici sul terreno [ogni segmento è formato da 2 nodi]
            (bisogna considerare anche i nuovi nuovi che vengono creati con il raffinamento)
        """

        import sys

        n_buildings = self.ModelPart.NumberOfSubModelParts()
        n_th = 1            # n-th building
        bar_length = 30     # length of progress bar

        # TODO: potrei salvarlo come self.dict_ext_edges
        dict_ext_edges = {}     # key: name of building; value: list of tuple of external edges [(n1,n2), (n2,n3), ...]

        for current_sub_model in self.ModelPart.SubModelParts:
            # progress bar. Reference: https://gist.github.com/sibosutd/c1d9ef01d38630750a1d1fe05c367eb8
            percent = 100.0*n_th/n_buildings
            sys.stdout.write('\r')
            sys.stdout.write("Mesh3D_buildings: [{:{}}] {:>3}%".format('='*int(percent/(100.0/bar_length)),
                                                                                 bar_length,
                                                                                 int(percent)))
            sys.stdout.flush()

            segments = []
            coord_2D = []   # useful for trianle
            coord_3D = []   # useful for tetgen

            all_facets = []
            all_markers = []    # 1: top and lateral faces (building); 2: bottom faces (building)

            list_to_mp_nodes = {}   # key: position in coord_2D; value: node.Id in model part
            mp_to_list_nodes = {}   # key: node id in Model Part; value: node position in coord_3D list
            
            n_nodes_base = current_sub_model.NumberOfNodes() / 2
            first_id = len(segments)
            for pos, node in enumerate(current_sub_model.Nodes):    # posizione e nodo nell'i-esimo edificio
                coord_3D.append([node.X, node.Y, node.Z])
                mp_to_list_nodes[node.Id] = len(coord_3D) - 1   # salviamo la posizione del nodo i-esimo nel dizionario
                if (pos >= n_nodes_base):
                    # we exclude all Nodes on the top of the buildings
                    continue
                
                list_to_mp_nodes[first_id+pos] = node.Id
                coord_2D.append([node.X, node.Y])
                if (pos == n_nodes_base - 1):
                    segments.append([first_id+pos, first_id])
                    continue
                segments.append([first_id+pos, first_id+pos+1])

            # Bottom buildings: triangle
            A = dict(vertices=coord_2D, segments=segments)
            triangle_dict_bottom_building = triangle.triangulate(A, "p")

            # import matplotlib.pyplot as plt
            # triangle.compare(plt, A, triangle_dict_bottom_building)
            # plt.show()
            
            properties = self.ModelPart.Properties[0]


            # top and lateral faces (building)
            for elem in current_sub_model.Elements:
                nodes = elem.GetNodes()
                
                n1 = mp_to_list_nodes[nodes[0].Id]
                n2 = mp_to_list_nodes[nodes[1].Id]
                n3 = mp_to_list_nodes[nodes[2].Id]
                all_facets.append([n1, n2, n3])
                all_markers.append(1)

            # bottom faces (building)
            id_elem = self._find_max_elem_id() + 1
            for nodes in triangle_dict_bottom_building["triangles"]:
                n1 = list_to_mp_nodes[nodes[0]]
                n2 = list_to_mp_nodes[nodes[1]]
                n3 = list_to_mp_nodes[nodes[2]]
                current_sub_model.CreateNewElement("Element2D3N",
                                                   id_elem,
                                                   [n1, n2, n3],
                                                   properties)
                id_elem += 1
                
                n1 = mp_to_list_nodes[n1]
                n2 = mp_to_list_nodes[n2]
                n3 = mp_to_list_nodes[n3]
                all_facets.append([n1, n2, n3])
                all_markers.append(2)


            """ TETGEN (VIA MESHPY) """
            ### using Meshpy for the creation of an initial mesh https://documen.tician.de/meshpy/
            mesh_info = MeshInfo()
            mesh_info.set_points(coord_3D)
            mesh_info.set_facets(all_facets, markers=all_markers)
            
            # build the mesh
            mesh = build(mesh_info, options=Options("pa50.0"))

            # we clean the current Sub Model Part
            # set all Nodes as to erase
            for node in current_sub_model.Nodes:
                node.Set(KratosMultiphysics.TO_ERASE,True)
            
            # set all Elements as to erase
            for elem in current_sub_model.Elements:
                elem.Set(KratosMultiphysics.TO_ERASE,True)
            
            # set all Conditions as to erase
            for cond in current_sub_model.Conditions:
                cond.Set(KratosMultiphysics.TO_ERASE,True)
            
            # delete all Nodes, Elements and Conditions
            self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
            self.ModelPart.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)
            self.ModelPart.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)


            ### Nodes
            node_id = self._find_max_node_id()
            for i in range(len(mesh.points)):
                # node id is shifted up by 1
                coords = mesh.points[i]
                current_sub_model.CreateNewNode((node_id+i+1), coords[0], coords[1], coords[2])
            
            ### Elements and Nodes in SubModelParts
            elem_id = self._find_max_elem_id()
            list_edges = []     # list with all edges at the bottom of the building
            for j in range(len(mesh.faces)):
                points = mesh.faces[j]
                marker = mesh.face_markers[j]
                if (marker == 2):
                    # we save all edges at the bottom (building)
                    # we saves the ids in ascending order to count the number of equal edges
                    list_edges.append((points[0]+node_id+1, points[1]+node_id+1)) if (points[0]+node_id+1 < points[1]+node_id+1) else list_edges.append((points[1]+node_id+1, points[0]+node_id+1))
                    list_edges.append((points[1]+node_id+1, points[2]+node_id+1)) if (points[1]+node_id+1 < points[2]+node_id+1) else list_edges.append((points[2]+node_id+1, points[1]+node_id+1))
                    list_edges.append((points[2]+node_id+1, points[0]+node_id+1)) if (points[2]+node_id+1 < points[0]+node_id+1) else list_edges.append((points[0]+node_id+1, points[2]+node_id+1))
                    continue
                current_sub_model.CreateNewElement("Element3D3N",
                                                    (elem_id+j+1),
                                                    [points[0]+node_id+1, points[1]+node_id+1, points[2]+node_id+1],
                                                    properties)
            
            ext_edges = []      # list with external edges
            ext_nodes = set()   # set with external nodes (to avoid duplicates)
            for edge in list_edges:
                if (list_edges.count(edge) == 1) and (not edge in ext_edges):
                    # if the edge is unique and not yet in ext_edges is an external edges
                    ext_edges.append(edge)      # we update ext_edges with unique edge
                    ext_nodes.add(edge[0])      # we update ext_nodes with nodes
                    ext_nodes.add(edge[1])

            dict_ext_edges[current_sub_model.Name] = ext_edges



            # ### Elements
            # elem_id = self._find_max_elem_id()
            # for k in range(len(mesh.elements)):
            #     nodes = mesh.elements[k]
            #     current_sub_model.CreateNewElement("Element3D4N",
            #                                        (elem_id+k+1),
            #                                        [nodes[0]+node_id+1, nodes[1]+node_id+1, nodes[2]+node_id+1, nodes[3]+node_id+1],
            #                                        properties)


            n_th += 1
        sys.stdout.write("\n")

        # set all Nodes as to erase
        for node in self.ModelPart.Nodes:
            node.Set(KratosMultiphysics.TO_ERASE,True)
        for elem in self.ModelPart.Elements:
            for node in elem.GetNodes():
                node.Set(KratosMultiphysics.TO_ERASE,False)
        for cond in self.ModelPart.Conditions:
            for node in cond.GetNodes():
                node.Set(KratosMultiphysics.TO_ERASE,False)
        
        # delete all free Nodes
        self.ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

        return dict_ext_edges
    #################################################################





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
                    # print("LateralSector_{}".format(n_sect))
                    # nnn = cond.GetNodes()
                    # print(cond.Id)
                    # for n in nnn:
                    #     print(f"\t{n.X} {n.Y} {n.Z}")
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
 

    def extract_center(self, terrain_model_part=None):
        """ function to create a ModelPart with "elements" of the central portion

        Args:
            terrain_model_part: ModelPart of the terrain. If "None" we use self.ModelPart

        Returns:
            ModelPart with the elements belonging to the central portion
        """

        if (terrain_model_part == None):
            terrain_model_part = self.ModelPart
        
        model = KratosMultiphysics.Model()
        self.center_model_part = model.CreateModelPart("center_model_part")
        prop = self.center_model_part.Properties[0]

        if (terrain_model_part.HasSubModelPart("BottomModelPart")):
            smp_bottom = terrain_model_part.GetSubModelPart("BottomModelPart")
            items = smp_bottom.Conditions
        else:
            items = self.ModelPart.Elements

        # for cond in smp_bottom.Conditions:
        for item in items:
            # nodes = cond.GetNodes()
            nodes = item.GetNodes()
            for node in nodes:
                dist = math.sqrt(((self.x_center-node.X)**2) + ((self.y_center-node.Y)**2))
                if (dist > self.r_buildings*1.2):
                    break
            else:
                for node in nodes:
                    self.center_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
                # self.center_model_part.CreateNewElement("Element2D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
                self.center_model_part.CreateNewElement("Element2D3N", item.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)




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


    def _find_max_node_id(self):
        " returns the maximum Node id "
        if (self.ModelPart.NumberOfNodes()):
            return max((node.Id for node in self.ModelPart.Nodes))
        else:
            return 0


    def _find_max_elem_id(self):
        " returns the maximum Element id "
        if (self.ModelPart.NumberOfElements()):
            return max((elem.Id for elem in self.ModelPart.Elements))
        else:
            return 0


    def _find_max_cond_id(self):
        " returns the maximum Condition id "
        if (self.ModelPart.NumberOfConditions()):
            return max((cond.Id for cond in self.ModelPart.Conditions))
        else:
            return 0


    def _find_top_building(self, current_building):
        " returns the Z coordinate of the top of the current building "
        return max((node.Z for node in current_building.Nodes))


    def truncate(self, num, ndigits=5):
        " returns a truncated number of ndigits"
        return int(num*(10**ndigits)) / (10**ndigits)


    def _towards(self, coords):
        " check if coords are clockwise or counterclockwise. return a string "
        # coords must be a list like [(x1,y1), (x2,y2), ...]

        # if _sum is positive, the coordinates are clockwise; otherwise counterclockwise
        _sum = 0
        _len_coords = len(coords)
        for i in range(_len_coords):
            if (i == _len_coords-1):
                # the last one: sum of (x1 − xn)(y1 + yn)
                _sum += ((coords[0][0]-coords[i][0])*(coords[0][1]+coords[i][1]))
                continue
            
            # sum of (x2 − x1)(y2 + y1)
            _sum += ((coords[i+1][0]-coords[i][0])*(coords[i+1][1]+coords[i][1]))
        
        if (_sum >= 0):
            return "clockwise"
        else:
            return "counterclockwise"


    def _sort_edges(self, edges):
        """ return a sorted list
            \"edges\" must be a list of list [(n1,n2), (n2,n3), ...]
        """
        
        len_target = len(edges)
        edges_sort = []     # new list with sorted edges
        edges_sort.append(edges[0])
        del edges[0]    # we delete edge alreay sorted

        for i in range(1, len_target):
            rif = edges_sort[-1][1]     # second node in the last couple of nodes
            for j in range(len(edges)):
                node1 = edges[j][0]
                node2 = edges[j][1]
                if (rif == node1):
                    edges_sort.append((node1, node2))
                    break
                elif (rif == node2):
                    edges_sort.append((node2, node1))
                    break
            del edges[j]
        
        return edges_sort
