# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.IGAStructuralMechanicsApplication
import KratosMultiphysics.NURBSBRepApplication
from KratosMultiphysics.SolidMechanicsApplication import *
from beam_sections_python_utility import SetProperties
import json
import math
#import KratosMultiphysics.model_part

def Factory(model_part_elements, settings):
	return ModelPartIOIGA(model_part_elements, settings)


class ModelPartIOIGA(KratosMultiphysics.IO):
	def __init__(self, model_part_elements, settings):
		KratosMultiphysics.IO.__init__(self)
		self.model_part_elements = model_part_elements
		self.settings = settings

	def ReadModelPart(self, model_part):
		dimension = 2

		nodes_list = self.model_part_elements.Nodes
		
		self.ReadSubModelPart(model_part)

		for property in model_part.Properties:
			properties = property

		cplist = model_part.Nodes
		for node in cplist:
			print(node.Id)
		
		for node in nodes_list:
			cps = node.GetValue(KratosMultiphysics.NURBSBRepApplication.CONTROL_POINT_IDS)
			intcps = []
			for cp in cps:
				intcps.append(int(cp))
			print(node.Id)
			element = model_part.CreateNewElement("MeshlessShellElement", node.Id, intcps, properties)
			element.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_VALUES, node.GetValue(KratosMultiphysics.NURBSBRepApplication.SHAPE_FUNCTION_VALUES))
			element.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES, node.GetValue(KratosMultiphysics.NURBSBRepApplication.SHAPE_FUNCTION_DERIVATIVES))
			element.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES, node.GetValue(KratosMultiphysics.NURBSBRepApplication.SHAPE_FUNCTION_SECOND_DERIVATIVES))
			element.SetValue(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT, node.GetValue(KratosMultiphysics.NURBSBRepApplication.INTEGRATION_WEIGHT))
		

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