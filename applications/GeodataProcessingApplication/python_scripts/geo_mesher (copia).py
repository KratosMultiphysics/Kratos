"""
    # edit:     02 September 2019 -> z_min and z_max are calculated only in the nodes of the central region (in MeshCircleWithTerrainPoints function)
    
"""

import KratosMultiphysics
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.MeshingApplication as KratosMesh

from geo_processor import GeoProcessor
import triangle
import numpy as np
import math
from meshpy.tet import MeshInfo, build

class GeoMesher( GeoProcessor ):

    def __init__( self ):
        super(GeoMesher, self).__init__()

        self.HasModelPart = False
        self.HasExtrusionHeight = False

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

    def ExtrudeBox(self, height):

        pass
    """ END TEST FUNCTIONS """


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
                cond = self.ModelPart.CreateNewCondition("SurfaceCondition3D3N", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                bottom_cond.append( j+1 )
            elif ( marker == 2 ):
                # cond = self.ModelPart.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                cond = self.ModelPart.CreateNewCondition("SurfaceCondition3D3N", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                top_cond.append( j+1 )
            elif ( marker == 3 ):
                # cond = self.ModelPart.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
                cond = self.ModelPart.CreateNewCondition("SurfaceCondition3D3N", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
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
        MmgProcess = KratosMesh.MmgProcess3D(self.ModelPart, remesh_param)
        MmgProcess.Execute()


    def MeshCircleWithTerrainPoints( self, h_value=0.0, circ_division=60, extract_center=False ):

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
        x_center = (x_min + x_max)/2
        y_center = (y_min + y_max)/2

        # radius calculation, considering the smaller side of the rectangle
        r_boundary = min((x_max-x_min), (y_max-y_min))/2    # radius of the domain
        r_ground = r_boundary * 2.0 / 3.0                   # 2/3 of the radius of the domain
        r_buildings = r_boundary / 3.0                      # 1/3 of the radius of the domain

        # buffering zone
        delta = r_boundary/1000         # evaluate this value!!!

        """ SMOOTHING PROCEDURE """
        del_id = []             # list where there are the ids to be deleted because are outside of r_boundary
        Z = []
        # for index in idx_to_smooth:
        for node in self.ModelPart.Nodes:
            x_curr = node.X             # current X coordinate
            y_curr = node.Y             # current Y coordinate

            dist = math.sqrt(((x_center-x_curr)**2)+((y_center-y_curr)**2))     # calculate the distance between current node and center of circle

            if (dist > r_boundary-delta):
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
        for th in theta:
            x_circle.append(r_boundary * math.cos(th) + x_center)   # we move the node along X to consider that the circle have not the centre in (0.0, 0.0)
            y_circle.append(r_boundary * math.sin(th) + y_center)   # we move the node along Y to consider that the circle have not the centre in (0.0, 0.0)

        """ TRIANGLE """
        coord_2D = []       # vertically list with X and Y coordinates (TO AVOID NUMPY VSTACK)
        all_points = []     # list with all coordinates of nodes
        for node in self.ModelPart.Nodes:
            coord_2D.append([node.X, node.Y])     # we add coordinates X and Y of terrain nodes
            all_points.append((node.X, node.Y, node.Z))

        # we add circumference nodes into coord_2D list
        for i in range(len(x_circle)):
            coord_2D.append([x_circle[i], y_circle[i]])

        ### preparing nodes to be handed over to triangle
        A = dict(vertices = coord_2D)
        triangle_dict = triangle.triangulate(A)     # triangle_dict {vertices: [...], triangles: [...], vertex_markers: [...]}

        inner_id_terrain = len(all_points)		# id up to this moment
        list_id_circle_bottom = []				# list with only ids of the bottom circle

        for i in range(len(x_circle)):
            all_points.append((x_circle[i], y_circle[i], z_min))
            list_id_circle_bottom.append(i + inner_id_terrain)

        all_facets = []		# list with all node ids of faces
        all_markers = []	# list of integers [1 = bottom, 2 = topper and 3 = lateral]

        # we fill the list for "bottom"
        for f in triangle_dict["triangles"]:
            all_facets.append([f[0], f[1], f[2]])
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
        for f in triangle_dict["triangles"]:
            all_facets.append([f[0]+number_node_terrain, f[1]+number_node_terrain, f[2]+number_node_terrain])		# list of lists. "...+number_node_terrain" is UGLY! improve it
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
        ### using Meshpy for the creation of an initial mesh https://mathema.tician.de/software/meshpy/
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


        """ OPTIONAL VALUES """
        # if extract_center:
        #     # we extract the center of the geometry and we add them in a sub model part
        #     center_elem = self.ModelPart.CreateSubModelPart("CenterElement")

        #     # elem_inside = {}        # dictionary with the information about the element that are in r_buildings
        #     # elem_id = 1

        #     elem_list = []  #########################################################################################################################

        #     for node in triangle_dict["triangles"]:
        #         for i in range(3):
        #             x_curr = self.ModelPart.GetNode(node[i]+1).X
        #             y_curr = self.ModelPart.GetNode(node[i]+1).Y
        #             dist = math.sqrt(((x_center-x_curr)**2) + ((y_center-y_curr)**2))

        #             if dist > r_buildings:
        #                 break   # if at least one node is external to r_buildings, we reject the entire element

        #         else:
        #             # we get here if all nodes of the element are inside the r_buildings
        #             node1 = self.ModelPart.GetNode(node[0]+1)
        #             node2 = self.ModelPart.GetNode(node[1]+1)
        #             node3 = self.ModelPart.GetNode(node[2]+1)

        #             # elem_inside[elem_id] = [[node1.Id, node1.X, node1.Y, node1.Z],
        #             #                         [node2.Id, node2.X, node2.Y, node2.Z],
        #             #                         [node3.Id, node3.X, node3.Y, node3.Z]]

        #             elem_list.append([node1.X, node1.Y, node1.Z]) ###################################################################################
        #             elem_list.append([node2.X, node2.Y, node2.Z]) ###################################################################################
        #             elem_list.append([node3.X, node3.Y, node3.Z]) ###################################################################################

        #             # elem_id += 1

        #     return elem_list    # JUST FOR THE CHECK

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


############################ --- Auxiliary functions --- ###########################################

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


    def _custom_range(self, start, stop, step=1.0):
        # function to replace the numpy range function
        list_floats = [float(start)]
        while (list_floats[-1]+step) < stop:
            list_floats.append(list_floats[-1]+step)
        return list_floats
