import KratosMultiphysics
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.MeshingApplication as KratosMesh

from geo_processor import GeoProcessor
import triangle
import numpy as np
from meshpy.tet import MeshInfo, build

class GeoMesher( GeoProcessor ):

    def __init__( self ):
        super(GeoMesher, self).__init__()

        self.HasExtrusionHeight = False

    ### --- functions to determine an extrusion height --- ###

    def ComputeExtrusionHeight( self, radius, height, free_board, iterations ):

        tool = KratosGeo.ExtrusionHeightUtilities( self.ModelPart )
        tool.SetExtrusionHeight( radius, height, free_board )
        tool.SmoothExtrusionHeight( radius/5.0, iterations, 0.8*free_board )
        self.HasExtrusionHeight = True

    def SetUniformExtrusionHeight( self, radius, height ):

        KratosMultiphysics.VariableUtils().SetNonHistoricalVariable( KratosGeo.EXTRUSION_HEIGHT, height, self.ModelPart.Nodes )
        self.HasExtrusionHeight = True


    ### --- functions to create a 3D domain mesh --- ###

    def MeshConcaveHullWithTerrainPoints( self ):

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
        Mesher.ReGenerateMesh("Element2D3N", "Condition2D", self.ModelPart, node_erase_process, False, False, 1.4, 0.5)

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

        # output (requries some vtk module, do you know more, Nicola?)
        # vtk_file_out = "mesh.vtk"
        # mesh.write_vtk(vtk_file_out)

		## create a new model part
        current_model = KratosMultiphysics.Model()
        main_model_part = current_model.CreateModelPart("MainModelPart")
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, 0.0)
        main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME, 1.0)

        # We add the variables needed (needed for the mmg process)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FLAG_VARIABLE)
        main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_BOUNDARY)

        lateral_model_part = main_model_part.CreateSubModelPart("LateralModelPart")
        bottom_model_part = main_model_part.CreateSubModelPart("BottomModelPart")
        top_model_part = main_model_part.CreateSubModelPart("TopModelPart")

        properties = main_model_part.Properties[0]

        bottom_cond = []; bottom_points = []
        top_cond = []; top_points = []
        lateral_cond = []; lateral_points = []

        # AWARE: Shift by 1 in index!!!

        ### Nodes
        for i in range( 0, len( mesh.points ) ):
        	# node id is shifted up by 1
        	coords = mesh.points[i]
        	node = main_model_part.CreateNewNode( (i+1), coords[0], coords[1], coords[2] )

        ### Conditions and Nodes in SubModelParts
        for j in range( 0, len( mesh.faces ) ):
        	points = mesh.faces[j]
        	marker = mesh.face_markers[j]
        	if ( marker == 1 ):
        		cond = main_model_part.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
        		bottom_cond.append( j+1 )
        	elif ( marker == 2 ):
        		cond = main_model_part.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
        		top_cond.append( j+1 )
        	elif ( marker == 3 ):
        		cond = main_model_part.CreateNewCondition("Condition3D", (j+1), [points[0]+1, points[1]+1, points[2]+1], properties )
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
        	elem = main_model_part.CreateNewElement("Element3D4N", (k+1), [ points[0]+1, points[1]+1, points[2]+1, points[3]+1 ], properties )

        print( bottom_model_part )
        print( top_model_part )
        print( lateral_model_part )

        self.ModelPart = main_model_part







    def MeshCircleWithTerrainPoints( self ):

        print("Not interfereing here, Nicola - you are the expert!")
