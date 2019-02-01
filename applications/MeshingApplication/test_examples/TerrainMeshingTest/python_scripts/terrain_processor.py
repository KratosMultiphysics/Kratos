from KratosMultiphysics import *

class terrain_processor:

	def __init__ ( self ):

		# definition of the data structure

		### dictionaries
		self.vertex = {}				# dictionary of all vertices - key = (x, y, z); value = id_vertex
		self.coord_3D = {}				# dictionary of all vertices - key = id_vertex; value = (x, y, z)

		### lists
		self.id_vertex = []				# list with IDs of vertices with Z != 0
		self.coord_2D = []				# DEPRECATED # list with 2D (x,y) coordinates

		self.nodes = []
		self.facets = []


	def FileCropXYZ(self, infile_name, outfile_name, bounding_box = (-1E20, -1E20, 1E20, 1E20) ):
		from terrain_preprocessor import xyz_cut
		xyz_cut( infile_name, outfile_name, bounding_box )
		print("Cropping process finished")

	def FileShiftCoordsXYZ(self, infile_name, outfile_name, x_shift, y_shift):
		from terrain_preprocessor import xyz_shift
		xyz_shift( infile_name, outfile_name, x_shift, y_shift )
		print("Coordinate shift process finished")

	def FileExtractHighRegionsXYZ(self, infile_name, outfile_name, min_height):
		from terrain_preprocessor import xyz_mountain
		xyz_mountain( infile_name, outfile_name, min_height )
		print("High regions extaction process finished")

	def FileExtractLowRegionsXYZ(self, infile_name, outfile_name, max_height):
		from terrain_preprocessor import xyz_valley
		xyz_valley( infile_name, outfile_name, max_height )
		print("Low regions extaction process finished")

	def FileExtractRiver(self, infile_name, outfile_name, radius, start_x_y):
		from terrain_preprocessor import xyz_river
		xyz_river( infile_name, outfile_name, radius, start_x_y )
		print("River bed extaction process finished")



	def Import(self, file_name_in ):

		import os
		file_name, extension_in = os.path.splitext( file_name_in )

		if extension_in.capitalize() == ".stl":

			from terrain_importer import stl_import
			self.id_vertex, self.vertex, self.coord_2D, self.coord_3D = stl_import( file_name_in )

		elif extension_in.capitalize() == ".xyz":

			from terrain_importer import xyz_import
			self.id_vertex, self.vertex, self.coord_2D, self.coord_3D = xyz_import( file_name_in )

		else:

			print( "No importer module found for file format " + extension_in.capitalize() )
			print( "Import of the following file types is possible: \n - *.stl, \n - *.xyz")


	def MeshConcaveHull( self ):

		###	  P R E P A R A T I O N   ###
		# evaluation of imported data
		X = []; Y = []; Z = []      # lists are generated to use predefined operations

		for _, coords in self.coord_3D.items():
			X.append(coords[0])
			Y.append(coords[1])
			Z.append(coords[2])

		# calculate minimum and maximum value of X, Y and Z coordinate
		x_min = min(X); x_max = max(X)
		y_min = min(Y); y_max = max(Y)
		z_min = min(Z); z_max = max(Z)

		vertex_id_list = self.id_vertex		# filling the id list

		import triangle
		import numpy as np

		# fill pts list with the couple x, y for 2D triangulation
		pts = np.vstack((X, Y)).T
		A = dict(vertices = pts)
		triangle_dict = triangle.triangulate(A)		# triangle_dict {vertices: [...], triangles: [...], vertex_markers: [...]}

		print( triangle_dict )

		height = 5.0

		all_bottom_points = []			# vector with ALL coordinates of nodes on the terrain
		all_bottom_facets = []			# vector with ALL node ids of faces
		all_top_points = []				# vector with ALL coordinates of nodes on the terrain


    	# fill all_xxx_points list
		for _, coords in self.coord_3D.items():
			## nodes stay as the were red in
			all_bottom_points.append((coords[0], coords[1], coords[2]))
			all_top_points.append((coords[0], coords[1], coords[2]+height))

		for f in triangle_dict["triangles"]:
			## facets are further checked
			all_bottom_facets.append([f[0], f[1], f[2]])

		maximal_initial_id_terrain_point = len( all_bottom_points )		# id up to this moment for the actual imported terrain
		maximal_initial_id_terrain_facet = len( all_bottom_facets )		# id up to this moment for the actual imported terrain

		import KratosMultiphysics

		## create a new model part
		current_model = KratosMultiphysics.Model()
		model_part = current_model.CreateModelPart("ModelPart")
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FLUID)
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FREE_SURFACE)
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_BOUNDARY)

		properties = model_part.Properties[0]
		### AWARE: Shift by 1 in index!!!
		for i in range( 0, maximal_initial_id_terrain_point ):
			# node id is shifted up by 1
			node = model_part.CreateNewNode( (i+1), (all_bottom_points[i])[0], (all_bottom_points[i])[1] , 0.0)

		for i in range( 0, maximal_initial_id_terrain_facet ):
			# node id is shifted up by 1
			nodes_for_element = [ (all_bottom_facets[i])[0]+1, (all_bottom_facets[i])[1]+1, (all_bottom_facets[i])[2]+1 ]
			# facet id is shifted up by 1
			elem = model_part.CreateNewElement("Element2D3N", (i+1), nodes_for_element, properties )

		## the mesh has a convex hull and many distorted elements exist
		self.test_output( model_part, "out1" )

		# defining the mesher
		import KratosMultiphysics.MeshingApplication
		Mesher = KratosMultiphysics.MeshingApplication.TriGenPFEMModeler()

		for node in model_part.Nodes:
			node.SetSolutionStepValue(KratosMultiphysics.NODAL_H, 0, 0.1)
			node.SetSolutionStepValue(KratosMultiphysics.IS_FLUID, 0, 1)

		node_erase_process = NodeEraseProcess(model_part)
		neigh_finder = FindNodalNeighboursProcess(model_part, 9, 18)
		neigh_finder.Execute()

		Mesher.ReGenerateMesh("Element2D3N", "Condition2D", model_part, node_erase_process, False, False, 1.4, 0.5)

		## the mesh has a good quality
		## conditions are introduced and allow an identification of the boundary
		self.test_output( model_part, "out2" )

		## read back from "element + condition" to repalce the old "facets"  ---
		all_bottom_facets = []
		all_bottom_markers = []				# BOTTOM = marker 1
		all_top_facets = []
		all_top_markers = []				# TOP = marker 2

		for elem in model_part.Elements:
			# the shift must be performed in reverse
			node1 = (elem.GetNodes()[0]).Id - 1
			node2 = (elem.GetNodes()[1]).Id - 1
			node3 = (elem.GetNodes()[2]).Id - 1
			# facets at the bottom
			all_bottom_facets.append( [ node1, node2, node3 ] )
			all_bottom_markers.append(1)
			# facets at the top
			all_top_facets.append( [ 	node1 + maximal_initial_id_terrain_point,
										node2 + maximal_initial_id_terrain_point,
										node3 + maximal_initial_id_terrain_point ] )
			all_top_markers.append(2)

		all_lateral_facets = []
		all_lateral_markers = []				# LATERAL = marker 3

		## creating the lateral limitting facets from the conditions in the model_part
		for cond in model_part.Conditions:
			# the shift must be performed in reverse
			node1 = (cond.GetNodes()[0]).Id - 1
			node2 = (cond.GetNodes()[1]).Id - 1
			node3 = node1 + maximal_initial_id_terrain_point
			node4 = node2 + maximal_initial_id_terrain_point
			all_lateral_facets.append([ node1, node2, node4, node3 ])
			all_lateral_markers.append(3)

		## putting the node lists together
		all_points = []
		all_points.extend( all_bottom_points )
		all_points.extend( all_top_points )

		## putting facet lists togehter (same order as markers)
		all_facets = []
		all_facets.extend( all_bottom_facets )
		all_facets.extend( all_top_facets )
		all_facets.extend( all_lateral_facets )
		print( len(all_facets) )

		## putting the markers together (same order as facets)
		all_markers = []
		all_markers.extend( all_bottom_markers )
		all_markers.extend( all_top_markers )
		all_markers.extend( all_lateral_markers )
		print( len(all_markers) )

		self.nodes = all_points
		self.facets = all_facets

		## create a new model part
		current_model = KratosMultiphysics.Model()
		model_part = current_model.CreateModelPart("ModelPart")
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FLUID)
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_FREE_SURFACE)
		model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_BOUNDARY)

		properties = model_part.Properties[0]
		### AWARE: Shift by 1 in index!!!
		for i in range( 0, len( all_points ) ):
			# node id is shifted up by 1
			node = model_part.CreateNewNode( (i+1), (all_points[i])[0], (all_points[i])[1] , (all_points[i])[2])

		for i in range( 0, len( all_facets ) ):
			# node id is shifted up by 1
			nodes_for_element = [ (all_facets[i])[0]+1, (all_facets[i])[1]+1, (all_facets[i])[2]+1 ]
			# facet id is shifted up by 1
			elem = model_part.CreateNewElement("Element2D3N", (i+1), nodes_for_element, properties )

		## the mesh has a convex hull and many distorted elements exist
		self.test_output( model_part, "out3" )


#######################################################################################################################
###														meshpy.tet													###
#######################################################################################################################
# https://mathema.tician.de/software/meshpy/

		from meshpy.tet import MeshInfo, build

		points = np.array( self.nodes )
		facets = np.array( self.facets )

		mesh_info = MeshInfo()
		mesh_info.set_points(points)
		mesh_info.set_facets(facets, markers=all_markers )

		# build the mesh
		mesh = build(mesh_info)

		# output
		vtk_file_out = "mesh.vtk"
		mesh.write_vtk(vtk_file_out)

### --- re-import into Kratos --------------------------------------------------------

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


		## DISTANCE on terrain ( No 1 ) ------------------------------------------------------------------------------------
		from KratosMultiphysics import MeshingApplication

		# assigning the initial distances
		for node in main_model_part.Nodes:
			node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,1.0)
		for node in bottom_model_part.Nodes:
			node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,-1e-7)

		# computing the distance field
		variational_distance_process = self._set_variational_distance_process_serial( main_model_part, "auxName1" )
		variational_distance_process.Execute()

		find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
		find_nodal_h.Execute()

		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)
		local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
		local_gradient.Execute()

		# We set to zero the metric
		ZeroVector = KratosMultiphysics.Vector(6)
		ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
		ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

		for node in main_model_part.Nodes:
			node.SetValue(MeshingApplication.METRIC_TENSOR_3D, ZeroVector)

		# We define a metric using the ComputeLevelSetSolMetricProcess
		level_set_param = KratosMultiphysics.Parameters("""
			{
				"minimal_size"                         : 2.0,
				"enforce_current"                      : false,
				"anisotropy_remeshing"                 : true,
				"anisotropy_parameters": {
					"hmin_over_hmax_anisotropic_ratio"      : 0.9,
					"boundary_layer_max_distance"           : 1.5,
					"interpolation"                         : "Linear" }
			}
			""")
		metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part,KratosMultiphysics.DISTANCE_GRADIENT,level_set_param)
		metric_process.Execute()

		# We create the remeshing process
		remesh_param = KratosMultiphysics.Parameters("""{ }""")
		MmgProcess = MeshingApplication.MmgProcess3D(main_model_part, remesh_param)
		MmgProcess.Execute()

		self.test_output_dist( main_model_part, "TerrainMesh1" )


		## DISTANCE on terrain ( No 2 ) ------------------------------------------------------------------------------------

		# assigning the initial distances
		for node in main_model_part.Nodes:
			node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,1.0)
		for node in bottom_model_part.Nodes:
			node.SetSolutionStepValue(KratosMultiphysics.DISTANCE,0,-1e-7)
		# computing the distance field
		variational_distance_process = self._set_variational_distance_process_serial( main_model_part, "auxName2" )
		variational_distance_process.Execute()

		find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
		find_nodal_h.Execute()

		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)
		local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
		local_gradient.Execute()

		# We set to zero the metric
		ZeroVector = KratosMultiphysics.Vector(6)
		ZeroVector[0] = 0.0; ZeroVector[1] = 0.0; ZeroVector[2] = 0.0
		ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

		for node in main_model_part.Nodes:
			node.SetValue(MeshingApplication.METRIC_TENSOR_3D, ZeroVector)

		# We define a metric using the ComputeLevelSetSolMetricProcess
		# level_set_param = KratosMultiphysics.Parameters("""
		# 	{
		# 		"minimal_size"                         : 0.2,
		# 		"enforce_current"                      : false,
		# 		"anisotropy_remeshing"                 : true,
		# 		"anisotropy_parameters": {
		# 			"hmin_over_hmax_anisotropic_ratio"      : 0.05,
		# 			"boundary_layer_max_distance"           : 0.1,
		# 			"interpolation"                         : "Linear" }
		# 	}
		# 	""")
		metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part,KratosMultiphysics.DISTANCE_GRADIENT,level_set_param)
		metric_process.Execute()

		# We create the remeshing process
		remesh_param = KratosMultiphysics.Parameters("""{ }""")
		MmgProcess = MeshingApplication.MmgProcess3D(main_model_part, remesh_param)
		MmgProcess.Execute()

		self.test_output_dist( main_model_part, "TerrainMesh2" )

############  --- under construction --- ##########################

		# function to load form MDPA to MODEL_PART
		file_path = os.path.dirname(os.path.realpath(__file__))
		building_model_part = current_model.CreateModelPart("BuildingModelPart")
		building_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
		building_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
		KratosMultiphysics.ModelPartIO(file_path + "/toy_house_skin").ReadModelPart(building_model_part)

		# self.test_output_dist( building_model_part, "outTest" )
		self.test_output_dist( main_model_part, "Building_raw" )

		# for node in main_model_part.Nodes:
		# 	node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.X )

		## SURFACE MESH
		mmg_parameters = KratosMultiphysics.Parameters("""
		{
		    "filename"                         : "toy_house_skin",
		    "save_external_files"              : true,
		    "echo_level"                       : 0
		}
		""")

		# We calculate the gradient of the distance variable
		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, building_model_part.Nodes)
		find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(building_model_part)
		find_nodal_h.Execute()

		# We set to zero the metric
		metric_vector = KratosMultiphysics.Vector(6)
		metric_vector[0] = 1.0; metric_vector[1] = 1.0; metric_vector[2] = 1.0
		metric_vector[3] = 0.0; metric_vector[4] = 0.0; metric_vector[5] = 0.0

		for node in building_model_part.Nodes:
		    node.SetValue(MeshingApplication.METRIC_TENSOR_3D, metric_vector)

		for node in building_model_part.Nodes:
		    node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.X )

		mmg_parameters["filename"].SetString(file_path + "/" + mmg_parameters["filename"].GetString())
		mmg_process = MeshingApplication.MmgProcess3DSurfaces(building_model_part, mmg_parameters)
		mmg_process.Execute()


		self.test_output_dist( main_model_part, "Building_remeshed" )


		# NOT a complete distance filed but ruther a specification of the 0-level only
		KratosMultiphysics.CalculateDistanceToSkinProcess3D(main_model_part, building_model_part).Execute()


		self.test_output_dist( main_model_part, "Terrain_with_ZeroLevel_1" )

		for node in main_model_part.Nodes:
			if ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) > 1000.0 ):
				node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.0 )
			elif ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) < -1000.0 ):
				node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -1.0 )


		# NOW creating the full distance field
		variational_distance_process = self._set_variational_distance_process_serial( main_model_part, "auxName3" )
		variational_distance_process.Execute()



		self.test_output_dist( main_model_part, "Terrain_with_FullDist_1" )



		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, main_model_part.Nodes)
		find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
		find_nodal_h.Execute()

		print("CP1")

		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)
		local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
		local_gradient.Execute()
		self.test_output_dist( main_model_part, "DistGrad" )

		print("CP2")

		# We set to zero the metric
		ZeroVector = KratosMultiphysics.Vector(6)
		ZeroVector[0] = 1.0; ZeroVector[1] = 1.0; ZeroVector[2] = 1.0
		ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

		for node in main_model_part.Nodes:
			node.SetValue(MeshingApplication.METRIC_TENSOR_3D, ZeroVector)

		print("CP3")

		# We define a metric using the ComputeLevelSetSolMetricProcess
		level_set_param = KratosMultiphysics.Parameters("""
			{
				"minimal_size"                         : 0.2,
				"enforce_current"                      : false,
				"anisotropy_remeshing"                 : true,
				"anisotropy_parameters": {
					"hmin_over_hmax_anisotropic_ratio"      : 0.3,
					"boundary_layer_max_distance"           : 0.4,
					"interpolation"                         : "Linear" }
			}
			""")
		## works without!!!
		metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part,KratosMultiphysics.DISTANCE_GRADIENT,level_set_param)
		metric_process.Execute()

		print("CP4")

		# We create the remeshing process
		remesh_param = KratosMultiphysics.Parameters("""{ }""")
		MmgProcess = MeshingApplication.MmgProcess3D(main_model_part, remesh_param)
		MmgProcess.Execute()


		self.test_output_dist( main_model_part, "Terrain_with Building_1" )


		### ------------- WORKING SO FAR -------------------------------------------------------------

		# NOT a complete distance filed but ruther a specification of the 0-level only
		KratosMultiphysics.CalculateDistanceToSkinProcess3D(main_model_part, building_model_part).Execute()


		self.test_output_dist( main_model_part, "Terrain_with_ZeroLevel_2" )

		for node in main_model_part.Nodes:
			if ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) > 1000.0 ):
				node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 1.0 )
			elif ( node.GetSolutionStepValue( KratosMultiphysics.DISTANCE) < -1000.0 ):
				node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, -1.0 )


		# NOW creating the full distance field
		variational_distance_process = self._set_variational_distance_process_serial( main_model_part, "auxName4" )
		variational_distance_process.Execute()

		self.test_output_dist( main_model_part, "Terrain_with_FullDist_2" )

		variable_utils = VariableUtils()
		variable_utils.SaveScalarVar( DISTANCE, DISTANCE, main_model_part.Nodes)

		############################################################################### TEST ###########################################

		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, main_model_part.Nodes)
		find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
		find_nodal_h.Execute()

		print("CP1")

		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)
		local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
		local_gradient.Execute()
		self.test_output_dist( main_model_part, "DistGrad" )

		print("CP2")

		# We set to zero the metric
		ZeroVector = KratosMultiphysics.Vector(6)
		ZeroVector[0] = 1.0; ZeroVector[1] = 1.0; ZeroVector[2] = 1.0
		ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

		for node in main_model_part.Nodes:
			node.SetValue(MeshingApplication.METRIC_TENSOR_3D, ZeroVector)

		print("CP3")

		# We define a metric using the ComputeLevelSetSolMetricProcess
		level_set_param = KratosMultiphysics.Parameters("""
			{
				"minimal_size"                         : 0.1,
				"enforce_current"                      : false,
				"anisotropy_remeshing"                 : true,
				"anisotropy_parameters": {
					"hmin_over_hmax_anisotropic_ratio"      : 0.1,
					"boundary_layer_max_distance"           : 0.05,
					"interpolation"                         : "Linear" }
			}
			""")
		## works without!!!
		metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part,KratosMultiphysics.DISTANCE_GRADIENT,level_set_param)
		metric_process.Execute()

		print("CP4")

		# We create the remeshing process
		remesh_param = KratosMultiphysics.Parameters("""{
				"discretization_type": "Isosurface"
		}""")
		MmgProcess = MeshingApplication.MmgProcess3D(main_model_part, remesh_param)
		MmgProcess.Execute()

		self.test_output_dist( main_model_part, "Terrain_ISOmesh" )


		ERRRRRRR   ### End due to error

		############################################################################### TEST ###########################################

		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, main_model_part.Nodes)
		find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
		find_nodal_h.Execute()

		print("CP1")

		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_AREA, 0.0, main_model_part.Nodes)
		local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(main_model_part, KratosMultiphysics.DISTANCE, KratosMultiphysics.DISTANCE_GRADIENT, KratosMultiphysics.NODAL_AREA)
		local_gradient.Execute()
		self.test_output_dist( main_model_part, "DistGrad" )

		print("CP2")

		# We set to zero the metric
		ZeroVector = KratosMultiphysics.Vector(6)
		ZeroVector[0] = 1.0; ZeroVector[1] = 1.0; ZeroVector[2] = 1.0
		ZeroVector[3] = 0.0; ZeroVector[4] = 0.0; ZeroVector[5] = 0.0

		for node in main_model_part.Nodes:
			node.SetValue(MeshingApplication.METRIC_TENSOR_3D, ZeroVector)

		print("CP3")

		# We define a metric using the ComputeLevelSetSolMetricProcess
		level_set_param = KratosMultiphysics.Parameters("""
			{
				"minimal_size"                         : 0.1,
				"enforce_current"                      : false,
				"anisotropy_remeshing"                 : true,
				"anisotropy_parameters": {
					"hmin_over_hmax_anisotropic_ratio"      : 0.1,
					"boundary_layer_max_distance"           : 0.05,
					"interpolation"                         : "Linear" }
			}
			""")
		## works without!!!
		metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part,KratosMultiphysics.DISTANCE_GRADIENT,level_set_param)
		metric_process.Execute()

		print("CP4")

		# We create the remeshing process
		remesh_param = KratosMultiphysics.Parameters("""{ }""")
		MmgProcess = MeshingApplication.MmgProcess3D(main_model_part, remesh_param)
		MmgProcess.Execute()


		self.test_output_dist( main_model_part, "Terrain_with Building_2" )



		KratosMultiphysics.CalculateDistanceToSkinProcess3D(main_model_part, building_model_part).Execute()


		self.test_output_dist( main_model_part, "Terrain_with_ZeroLevel_3" )

		# REMESH No 2
		for node in main_model_part.Nodes:
			node.SetValue(MeshingApplication.METRIC_TENSOR_3D, ZeroVector)

		# We define a metric using the ComputeLevelSetSolMetricProcess
		level_set_param = KratosMultiphysics.Parameters("""
								{
									"minimal_size"                         : 3.0,
									"enforce_current"                      : false,
									"anisotropy_remeshing"                 : true,
									"anisotropy_parameters":
									{
										"hmin_over_hmax_anisotropic_ratio"      : 0.9,
										"boundary_layer_max_distance"           : 1.5,
										"interpolation"                         : "Linear"
									}
								}
								""")

		KratosMultiphysics.VariableUtils().SetNonHistoricalVariable(KratosMultiphysics.NODAL_H, 0.0, main_model_part.Nodes)
		find_nodal_h = KratosMultiphysics.FindNodalHNonHistoricalProcess(main_model_part)
		find_nodal_h.Execute()

		metric_process = MeshingApplication.ComputeLevelSetSolMetricProcess3D(main_model_part,KratosMultiphysics.DISTANCE_GRADIENT,level_set_param)
		metric_process.Execute()

		# KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BLOCKED, True, building_sub_model_part.Conditions)

		remesh_param = KratosMultiphysics.Parameters("""{ }""")

		MmgProcess = MeshingApplication.MmgProcess3D(main_model_part, remesh_param)
		MmgProcess.Execute()

		## DISTANCE TO SKIN
		KratosMultiphysics.CalculateDistanceToSkinProcess3D(main_model_part, building_model_part).Execute()
		print( "CP6" )

		self.test_output_dist( main_model_part, "out11" )



	def test_output( self, model_part, name ):

		import KratosMultiphysics
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
		                                        "nodal_results"       : ["IS_BOUNDARY"],
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



	def test_output_dist( self, model_part, name ):

		import KratosMultiphysics
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
		import KratosMultiphysics
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


	def _generate_building_model_part( self, current_model ):

		a = -5.0
		building_model_part = current_model.CreateModelPart("BuildingModelPart")
		properties = building_model_part.Properties[0]
		building_model_part.CreateNewNode( 1  , 30.0000000000 , 50.0000000000 , 660.0000000000 + a)
		building_model_part.CreateNewNode( 2  , 24.8715445833 , 44.7978514583 , 665.0758117317 + a)
		building_model_part.CreateNewNode( 3  , 30.0000000000 , 40.0000000000 , 660.0000000000 + a)
		building_model_part.CreateNewNode( 4  , 20.0000000000 , 50.0000000000 , 660.0000000000 + a)
		building_model_part.CreateNewNode( 5  , 30.0000000000 , 50.0000000000 , 670.0000000000 + a)
		building_model_part.CreateNewNode( 6  , 20.0000000000 , 40.0000000000 , 660.0000000000 + a)
		building_model_part.CreateNewNode( 7  , 30.0000000000 , 40.0000000000 , 670.0000000000 + a)
		building_model_part.CreateNewNode( 8  , 20.0000000000 , 50.0000000000 , 670.0000000000 + a)
		building_model_part.CreateNewNode( 9  , 20.0000000000 , 40.0000000000 , 670.0000000000 + a)
		building_model_part.CreateNewNode( 10 ,  22.8333333333,  44.1181423743,  677.1666666667+ a )
		building_model_part.CreateNewNode( 11 ,  30.0000000000,  50.0000000000,  680.0000000000+ a )
		building_model_part.CreateNewNode( 12 ,  30.0000000000,  40.0000000000,  680.0000000000+ a )
		building_model_part.CreateNewNode( 13 ,  20.0000000000,  50.0000000000,  680.0000000000+ a )
		building_model_part.CreateNewNode( 14 ,  20.0000000000,  40.0000000000,  680.0000000000+ a )
		building_model_part.CreateNewElement("Element3D3N", 1,  [4, 6, 1], properties )
		building_model_part.CreateNewElement("Element3D3N", 2,  [1, 6, 3], properties )
		building_model_part.CreateNewElement("Element3D3N", 3,  [4, 1, 5], properties )
		building_model_part.CreateNewElement("Element3D3N", 4,  [5, 11, 8], properties )
		building_model_part.CreateNewElement("Element3D3N", 5,  [8, 11, 13], properties )
		building_model_part.CreateNewElement("Element3D3N", 6,  [5, 8, 4], properties )
		building_model_part.CreateNewElement("Element3D3N", 7,  [1, 3, 7], properties )
		building_model_part.CreateNewElement("Element3D3N", 8,  [7, 12, 11], properties )
		building_model_part.CreateNewElement("Element3D3N", 9,  [11, 5, 7], properties )
		building_model_part.CreateNewElement("Element3D3N", 10, [7, 5, 1], properties )
		building_model_part.CreateNewElement("Element3D3N", 11, [3, 6, 9], properties )
		building_model_part.CreateNewElement("Element3D3N", 12, [9, 14, 12], properties )
		building_model_part.CreateNewElement("Element3D3N", 13, [12, 7, 9], properties )
		building_model_part.CreateNewElement("Element3D3N", 14, [9, 7, 3], properties )
		building_model_part.CreateNewElement("Element3D3N", 15, [6, 4, 8], properties )
		building_model_part.CreateNewElement("Element3D3N", 16, [8, 13, 9], properties )
		building_model_part.CreateNewElement("Element3D3N", 17, [9, 13, 14], properties )
		building_model_part.CreateNewElement("Element3D3N", 18, [8, 9, 6], properties )
		building_model_part.CreateNewElement("Element3D3N", 19, [13, 14, 11], properties )
		building_model_part.CreateNewElement("Element3D3N", 20, [11, 14, 12], properties )
		print( building_model_part )

		return building_model_part