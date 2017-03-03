# importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.IGAStructuralMechanicsApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from beam_sections_python_utility import SetProperties
import json
import math
#import KratosMultiphysics.model_part

def Factory(JsonFileName, settings, projectparameters):
	return ModelPartIOIGA(JsonFileName, settings, projectparameters)


class ModelPartIOIGA(KratosMultiphysics.IO):
	def __init__(self, JsonFileName, settings, projectparameters):
		KratosMultiphysics.IO.__init__(self)
		self.JsonFileName = JsonFileName
		self.settings = settings
		self.projectparameters = projectparameters

	def ReadModelPart(self, model_part):
		dimension = 2

		with open(self.JsonFileName) as JsonFile:
			nodes_list = json.load(JsonFile)['nodes']
			
		with open(self.JsonFileName) as JsonFile:
			dim2_elements_list = json.load(JsonFile)['2d_elements']

		with open(self.JsonFileName) as JsonFile:
			brep_elements_list = json.load(JsonFile)['brep_elements']

		self.ReadSubModelPart(model_part)

		#self.ReadProperties(model_part)

		for property in model_part.Properties:
			properties = property

		weights = {}
		for node in nodes_list:
			node1 = model_part.CreateNewNode(int(node[0]),node[1][0],node[1][1], node[1][2])
			weights.update({node[0] : node[1][3]})
			#print("Node: ", node[0],node[1][0],node[1][1], node[1][2])


		IGAModeler = {}
		ElementNodeIds = {}

		for patch in dim2_elements_list:
			patchID = patch[0]
			
			for dim2_element in patch[1]:
				nodes_for_element_ids = []
				weights_for_element = []
				for i in range (0, len(dim2_element[3])):
					nodes_for_element_ids.append(dim2_element[3][i])
					weights_for_element.append(weights[dim2_element[3][i]])


				ElementNodeIds.update({dim2_element[0]:nodes_for_element_ids})

				PolynomialDegreeXi = dim2_element[1][0]
				PolynomialDegreeEta = dim2_element[1][1]
				KnotsXi = dim2_element[2][0]
				KnotsEta = dim2_element[2][1]
		
				IGAModeler.update({dim2_element[0] : NurbsShapeFunctionModeler()})

				IGAModeler[dim2_element[0]].SetUp(model_part, dimension, nodes_for_element_ids, weights_for_element, PolynomialDegreeXi, PolynomialDegreeEta, KnotsXi, KnotsEta)
				for quatrature_point in dim2_element[5]:
					QuadraturePointID = quatrature_point[0]
					QuadraturePointWeightingFactor = quatrature_point[1]

					LocalCoordinatesOfQuadraturePoint = KratosMultiphysics.Vector(3)
					LocalCoordinatesOfQuadraturePoint[0] = quatrature_point[2][0]
					LocalCoordinatesOfQuadraturePoint[1] = quatrature_point[2][1]
					LocalCoordinatesOfQuadraturePoint[2] = 0.0000000000

					#local variable definition
					ShapeFunctionValues = KratosMultiphysics.Vector()
					ShapeFunctionLocalDerivatives = KratosMultiphysics.Matrix()
					ElementType = 1
					if (ElementType == 0):
						IGAModeler[dim2_element[0]].EvaluateShapeFunction(LocalCoordinatesOfQuadraturePoint, ShapeFunctionValues, ShapeFunctionLocalDerivatives)
						elem = model_part.CreateNewElement("MeshlessMembraneElement", QuadraturePointID, nodes_for_element_ids, properties)
		
						elem.SetValue(SHAPE_FUNCTION_VALUES, ShapeFunctionValues)
						elem.SetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES, ShapeFunctionLocalDerivatives)
						elem.SetValue(INTEGRATION_WEIGHT, QuadraturePointWeightingFactor)

					elif (ElementType == 1):
						#IGAModeler[-1].EvaluateShapeFunction(LocalCoordinatesOfQuadraturePoint, ShapeFunctionValues, ShapeFunctionLocalDerivatives)
						IGAModeler[dim2_element[0]].EvaluateShapeFunctionSecondOrder(LocalCoordinatesOfQuadraturePoint, ShapeFunctionValues, ShapeFunctionLocalDerivatives)

						elem = model_part.CreateNewElement("MeshlessShellElement", QuadraturePointID, nodes_for_element_ids, properties)

						elem.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_VALUES, ShapeFunctionValues)
						elem.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES, ShapeFunctionLocalDerivatives)
						elem.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT, QuadraturePointWeightingFactor)


		#--- CONDITIONS ---#
		with open("ProjectParameters.json") as JsonFile:
			load_list = json.load(JsonFile)["loads_process_list"]

		for load in load_list:
			loadID = load["parameters"]["brep_id"]
			for element in dim2_elements_list:
				if (loadID == element[0]):
					self.ReadConditions(model_part, element, load, IGAModeler, ElementNodeIds)
			for element in brep_elements_list:
				if (loadID == element[0]):
					self.ReadConditions(model_part, element, load, IGAModeler, ElementNodeIds)


		for brep_element in brep_elements_list:
			brepID = brep_element[0]

			parameter_file = open("ProjectParameters.json",'r')
			ProjectParameters = KratosMultiphysics.Parameters( parameter_file.read())
			
			Condition = 0
		
			with open("ProjectParameters.json") as JsonFile:
				support_list = json.load(JsonFile)["constraints_process_list"]
			with open("ProjectParameters.json") as JsonFile:
				coupling_list = json.load(JsonFile)["coupling_condition_list"]

			for support in support_list:
				#if ProjectParameters["loads_process_list"]["Parameters"]["mesh_id"] in property:
				if support["parameters"]["brep_id"] == brepID:
					Condition = "support"
					brepProperty = support

			for coupling in coupling_list:
				#if ProjectParameters["loads_process_list"]["Parameters"]["mesh_id"] in property:
				if coupling["parameters"]["brep_id"] == brepID:
					Condition = "coupling"
					brepProperty = coupling


			for brep in brep_element[1]:
				
				for quadrature_point in brep[1]:
					cond_id = quadrature_point[1][0]
					integration_weight = quadrature_point[1][1]
					
					if Condition == "coupling":
						LocalCoordinatesOfQuadraturePointMaster = KratosMultiphysics.Vector(3)
						LocalCoordinatesOfQuadraturePointMaster[0] = quadrature_point[1][2][0]
						LocalCoordinatesOfQuadraturePointMaster[1] = quadrature_point[1][2][1]
						LocalCoordinatesOfQuadraturePointMaster[2] = 0.0000000000
	
						LocalCoordinatesOfQuadraturePointSlave = KratosMultiphysics.Vector(3)
						LocalCoordinatesOfQuadraturePointSlave[0] = quadrature_point[1][4][0]
						LocalCoordinatesOfQuadraturePointSlave[1] = quadrature_point[1][4][1]
						LocalCoordinatesOfQuadraturePointSlave[2] = 0.0000000000

						#local variable definition
						ShapeFunctionValuesMaster = KratosMultiphysics.Vector()
						ShapeFunctionLocalDerivativesMaster = KratosMultiphysics.Matrix()

						#local variable definition
						ShapeFunctionValuesSlave = KratosMultiphysics.Vector()
						ShapeFunctionLocalDerivativesSlave = KratosMultiphysics.Matrix()

						IGAModeler[quadrature_point[0][0]].EvaluateShapeFunction(LocalCoordinatesOfQuadraturePointMaster, ShapeFunctionValuesMaster, ShapeFunctionLocalDerivativesMaster)
						IGAModeler[quadrature_point[0][1]].EvaluateShapeFunction(LocalCoordinatesOfQuadraturePointSlave, ShapeFunctionValuesSlave, ShapeFunctionLocalDerivativesSlave)


						nodeIdsMaster = []
						nodeIdsSlave = []
		
						nodeIdsMaster = ElementNodeIds[quadrature_point[0][0]]
						nodeIdsSlave = ElementNodeIds[quadrature_point[0][1]]

						n = len(ShapeFunctionValuesMaster)
						m = len(ShapeFunctionValuesSlave)

						#if (brepProperty["parameters"]["condition_name"] == "ContinuityConditionPenalty"):
						#	penalty = math.sqrt(brepProperty["parameters"]["penalty_factor"]) # math.pow(10,7)
						#	nodeIds = []
						#	penaltyVector = KratosMultiphysics.Vector(n+m) 
						#	for index in range (0, n):
						#		nodeIds.append(nodeIdsMaster[index])
						#		penaltyVector[index] = ShapeFunctionValuesMaster[index]*penalty
						#	for index in range (0, m):
						#		nodeIds.append(nodeIdsSlave[index])
						#		penaltyVector[n+index] = - ShapeFunctionValuesSlave[index]*penalty


						#	sub_model_part = model_part.GetSubModelPart(brepProperty["parameters"]["model_part_name"])
						#	for nodeID in nodeIds:
						#		sub_model_part.AddNode(model_part.GetNode(nodeID),0)
						#	Master = sub_model_part.CreateNewCondition("ContinuityConditionPenalty", cond_id, nodeIds, properties)
						#	Master.SetValue(KratosMultiphysics.PENALTY, penaltyVector)
						#	Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER, ShapeFunctionLocalDerivativesMaster)
						#	Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT, integration_weight)

						#	#Penalty
						#	penalty = brepProperty["parameters"]["penalty_factor"] # math.pow(10,7)
						#	Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.PENALTY_FACTOR, penalty)
							
						#	LocalTangents = KratosMultiphysics.Vector(4)
						#	LocalTangents[0] = quadrature_point[1][3][0]
						#	LocalTangents[1] = quadrature_point[1][3][1]
						#	LocalTangents[2] = quadrature_point[1][5][0]
						#	LocalTangents[3] = quadrature_point[1][5][1]
						#	Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.TANGENTS, LocalTangents)


						#elif (brepProperty["parameters"]["condition_name"] == "ContinuityConditionLagrange"):
						#	nodeIds = []
						#	LagrangeVector = KratosMultiphysics.Vector(2*(n+m)) 
						#	for index in range (0, n):
						#		nodeIds.append(nodeIdsMaster[index])
						#		LagrangeVector[index] = -ShapeFunctionValuesMaster[index]
						#	for index in range (0, m):
						#		nodeIds.append(nodeIdsSlave[index])
						#		LagrangeVector[n+index] = ShapeFunctionValuesSlave[index]
						#	for index in range (0, n):
						#		LagrangeVector[index+n+m] = ShapeFunctionValuesMaster[index]
						#	for index in range (n, m+n):
						#		LagrangeVector[index+n+m] = 0#- ShapeFunctionValuesSlave[index-n]#0#ShapeFunctionValuesMaster[index-n]

						#	Master = model_part.CreateNewCondition("ContinuityConditionLagrange", cond_id, nodeIds, properties)
						#	Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT, integration_weight)
						#	Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_VALUES, LagrangeVector)

						if (brepProperty["parameters"]["condition_name"] == "MeshlessPenaltyCouplingRotationCondition"):
							nodeIds = []
							ShapeFunctionVector = KratosMultiphysics.Vector(n+m) 
							for index in range (0, n):
								nodeIds.append(nodeIdsMaster[index])
								ShapeFunctionVector[index] = ShapeFunctionValuesMaster[index]
							for index in range (0, m):
								nodeIds.append(nodeIdsSlave[index])
								ShapeFunctionVector[n+index] = - ShapeFunctionValuesSlave[index]

							Master = model_part.CreateNewCondition("MeshlessPenaltyCouplingRotationCondition", cond_id, nodeIds, properties)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT, integration_weight)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_VALUES, ShapeFunctionVector)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER, ShapeFunctionLocalDerivativesMaster)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE, ShapeFunctionLocalDerivativesSlave)

							#PENALTY
							penalty = brepProperty["parameters"]["penalty_factor"] # math.pow(10,7)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.PENALTY_FACTOR, penalty)

							#TANGENTS
							LocalTangents = KratosMultiphysics.Vector(4)
							LocalTangents[0] = quadrature_point[1][3][0]
							LocalTangents[1] = quadrature_point[1][3][1]
							LocalTangents[2] = quadrature_point[1][5][0]
							LocalTangents[3] = quadrature_point[1][5][1]
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.TANGENTS, LocalTangents)

							DisplacementRotationFix = 0 #defined by rot, dispx, dispy, dispz
							if (brepProperty["parameters"]["is_fixed_rot"]):
								DisplacementRotationFix += 1000
							if (brepProperty["parameters"]["is_fixed_x"]):
								DisplacementRotationFix += 100
							if (brepProperty["parameters"]["is_fixed_y"]):
								DisplacementRotationFix += 10
							if (brepProperty["parameters"]["is_fixed_z"]):
								DisplacementRotationFix += 1

							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.DISPLACEMENT_ROTATION_FIX, DisplacementRotationFix)

						elif (brepProperty["parameters"]["condition_name"] == "MeshlessLagrangeCouplingCondition"):
							nodeIds = []
							LagrangeVector = KratosMultiphysics.Vector(2*(n+m)) 
							for index in range (0, n):
								nodeIds.append(nodeIdsMaster[index])
								LagrangeVector[index] = ShapeFunctionValuesMaster[index]
							for index in range (0, m):
								nodeIds.append(nodeIdsSlave[index])
								LagrangeVector[n+index] = -ShapeFunctionValuesSlave[index]
							for index in range (0, n):
								LagrangeVector[index+n+m] = ShapeFunctionValuesMaster[index]
							for index in range (n, m+n):
								LagrangeVector[index+n+m] = 0

							Master = model_part.CreateNewCondition("MeshlessLagrangeCouplingCondition", cond_id, nodeIds, properties)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT, integration_weight)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_VALUES, LagrangeVector)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES_MASTER, ShapeFunctionLocalDerivativesMaster)
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES_SLAVE, ShapeFunctionLocalDerivativesSlave)


							LocalTangents = KratosMultiphysics.Vector(4)
							LocalTangents[0] = quadrature_point[1][3][0]
							LocalTangents[1] = quadrature_point[1][3][1]
							LocalTangents[2] = quadrature_point[1][5][0]
							LocalTangents[3] = quadrature_point[1][5][1]
							Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.TANGENTS, LocalTangents)


					if Condition == "support":
						LocalCoordinatesOfQuadraturePointMaster = KratosMultiphysics.Vector(3)
						LocalCoordinatesOfQuadraturePointMaster[0] = quadrature_point[1][2][0]
						LocalCoordinatesOfQuadraturePointMaster[1] = quadrature_point[1][2][1]
						LocalCoordinatesOfQuadraturePointMaster[2] = 0.0000000000
						#local variable definition
						ShapeFunctionValuesMaster = KratosMultiphysics.Vector()
						ShapeFunctionLocalDerivativesMaster = KratosMultiphysics.Matrix()

						IGAModeler[quadrature_point[0][0]].EvaluateShapeFunction(LocalCoordinatesOfQuadraturePointMaster, ShapeFunctionValuesMaster, ShapeFunctionLocalDerivativesMaster)


						nodeIdsMaster = []
						nodeIdsMaster = ElementNodeIds[quadrature_point[0][0]]
						n = len(ShapeFunctionValuesMaster)

						nodeIds = []
						Penalty = math.sqrt(brepProperty["parameters"]["penalty_factor"])
						for index in range (0, n):
							nodeIds.append(nodeIdsMaster[index])
						sub_model_part = model_part.GetSubModelPart(brepProperty["parameters"]["model_part_name"])

						for nodeID in nodeIds:
							sub_model_part.AddNode(model_part.GetNode(nodeID),0)

						QuadratureWeighting = quadrature_point[1][1]				
						Master = sub_model_part.CreateNewCondition("MeshlessSupportRotationCondition", cond_id, nodeIds, properties)
						Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_VALUES, ShapeFunctionValuesMaster)
						Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES, ShapeFunctionLocalDerivativesMaster)
						Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT, QuadratureWeighting)
						
						LocalTangents = KratosMultiphysics.Vector(2)
						LocalTangents[0] = quadrature_point[1][3][0]
						LocalTangents[1] = quadrature_point[1][3][1]
						Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.TANGENTS, LocalTangents)

						#PENALTY
						penalty = brepProperty["parameters"]["penalty_factor"] # math.pow(10,7)
						Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.PENALTY_FACTOR, penalty)
		print("Read Model Part finished")

	def ReadProperties(self, model_part):
		with open(self.PatchFileName) as JsonFile:
			properties_list = json.load(JsonFile)['materials']
		for property in properties_list:
			self.ReadPropertiesBlock(model_part, property)

	def ReadPropertiesBlock(self, model_part, property):
		property_id = property['material_id']
		properties = model_part.Properties[property_id]

		if 'youngs_modulus' in property:
			properties.SetValue( KratosMultiphysics.YOUNG_MODULUS, property['youngs_modulus'] )
		if 'nue' in property:
			properties.SetValue( KratosMultiphysics.POISSON_RATIO, property['nue'] )
		if 'thickness' in property:
			properties.SetValue( KratosMultiphysics.THICKNESS, property['thickness'] )
		
		properties.SetValue( KratosMultiphysics.THICKNESS, 1 )
		properties.SetValue( KratosMultiphysics.BODY_FORCE, [0.0, 0.0, 0.0] )
		properties.SetValue( KratosMultiphysics.DENSITY, 1 )

		mat = LinearElasticPlaneStress2DLaw();
		#mat.Initialize();
		properties.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone());

	def ReadSubModelPart(self, model_part):
		for i in range(self.settings["processes_sub_model_part_list"].size()):
			part_name = self.settings["processes_sub_model_part_list"][i].GetString()
			model_part.CreateSubModelPart(part_name)

	def ReadConditions(self, model_part, condition, property, IGAModeler, ElementNodeIds):
		#if (property["description"] == "Edge"):

		for element in condition[1]:
			if (property["description"] == "Surface"):
				qpID = 5
				elementID = int(element[0])
			if (property["description"] == "Edge"):
				qpID = 1

			for quadrature_point_system in element[qpID]:
				if (property["description"] == "Surface"):
					quadrature_point = quadrature_point_system
					cond_id = int(quadrature_point[0])
				
				if (property["description"] == "Edge"):
					quadrature_point = quadrature_point_system[1]
					elementID = quadrature_point_system[0][0]
					cond_id = int(quadrature_point_system[1][0])

				LocalCoordinatesOfQuadraturePointMaster = KratosMultiphysics.Vector(3)
				LocalCoordinatesOfQuadraturePointMaster[0] = quadrature_point[2][0]
				LocalCoordinatesOfQuadraturePointMaster[1] = quadrature_point[2][1]
				LocalCoordinatesOfQuadraturePointMaster[2] = 0.0000000000

				#local variable definition
				ShapeFunctionValuesMaster = KratosMultiphysics.Vector()
				ShapeFunctionLocalDerivativesMaster = KratosMultiphysics.Matrix()


				IGAModeler[elementID].EvaluateShapeFunction(LocalCoordinatesOfQuadraturePointMaster, ShapeFunctionValuesMaster, ShapeFunctionLocalDerivativesMaster)

				nodeIdsMaster = []
				nodeIdsMaster = ElementNodeIds[elementID]

				n = len(ShapeFunctionValuesMaster)

				nodeIds = []
				ShapeFunctionVector = KratosMultiphysics.Vector(n) 
				for index in range (0, n):
					nodeIds.append(nodeIdsMaster[index])
					ShapeFunctionVector[index] = ShapeFunctionValuesMaster[index]
				sub_model_part = model_part.GetSubModelPart(property["parameters"]["model_part_name"])

				for nodeID in nodeIds:
					sub_model_part.AddNode(model_part.GetNode(nodeID),0)

				QuadratureWeighting = quadrature_point[1]
				Master = sub_model_part.CreateNewCondition(property["condition_name"], cond_id, nodeIds, model_part.Properties[0])
				Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_VALUES, ShapeFunctionVector)
				Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES, ShapeFunctionLocalDerivativesMaster)
				Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT, QuadratureWeighting)
				
				if (property["description"] == "Edge"):
					LocalTangents = KratosMultiphysics.Vector(2)
					LocalTangents[0] = quadrature_point[3][0]
					LocalTangents[1] = quadrature_point[3][1]
					Master.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.TANGENTS, LocalTangents)