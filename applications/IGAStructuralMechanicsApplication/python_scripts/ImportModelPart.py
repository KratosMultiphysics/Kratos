# importing the Kratos Library
from KratosMultiphysics import *
import KratosMultiphysics.IGAStructuralMechanicsApplication
import KratosMultiphysics.NurbsBrepApplication
from beam_sections_python_utility import SetProperties
import json
import math
#import KratosMultiphysics.model_part

def Factory(model_part_elements, project_parameters):
	return ModelPartIOIGA(model_part_elements, project_parameters)


class ModelPartIOIGA():
	def __init__(self, model_part_elements, project_parameters):
		self.model_part_elements = model_part_elements
		self.project_parameters = project_parameters

	def ReadModelPart(self, model_part):
		dimension = 2

		for element_id in range(0,self.project_parameters["element_list"].size()): #project_parameters["element_list"]:
			element_type = self.project_parameters["element_list"][element_id]
			element_name = element_type["element_name"].GetString()
			if(model_part.HasSubModelPart(element_type["iga_model_part"].GetString())):
				iga_model_part = model_part.GetSubModelPart(element_type["iga_model_part"].GetString())
			else:
				iga_model_part = model_part.CreateSubModelPart(element_type["iga_model_part"].GetString())
			#iga_model_part = model_part.CreateSubModelPart(element_type["iga_model_part"].GetString())
			#cad_model_part_names = element_type["cad_model_part"]
			#properties = model_part.GetProperties(element_type["parameters"]["properties_id"].GetInt())

			for property in model_part.Properties:
				if (property.Id == element_type["parameters"]["properties_id"].GetInt()):
					properties = property
					break
			for cad_model_part_name_id in range (0,element_type["cad_model_part"].size()):
				cad_model_part_name = element_type["cad_model_part"][cad_model_part_name_id]
				#cad_model_part  = self.model_part_elements.GetSubModelPart(cad_model_part_name.GetString())
				#model = KratosMultiphysics.Model()
				#model.AddModelPart(self.model_part_elements)
				cad_model_part  = model[cad_model_part_name.GetString()]

			var = KratosMultiphysics.KratosGlobals.GetVariable(element_type["parameters"]["control_points"].GetString())
			list_var_nurbs = []
			list_var_iga = []
			for variable in range(0,element_type["parameters"]["variables"].size()):
				list_var_nurbs.append(KratosMultiphysics.KratosGlobals.GetVariable(element_type["parameters"]["variables"][variable]["cad_variable"].GetString()))
				list_var_iga.append(KratosMultiphysics.KratosGlobals.GetVariable(element_type["parameters"]["variables"][variable]["iga_variable"].GetString()))
			#for node in cad_model_part.Nodes:
			#	print(node.Id)
			#	print(node.X)
			#	print(node.Y)
			#	print(node.Z)
			#	print(node.GetValue(KratosMultiphysics.INTEGRATION_WEIGHT))
			#	print(node.GetValue(KratosMultiphysics.NurbsBrepApplication.LOCAL_PARAMETERS))
			#	print(node.GetValue(KratosMultiphysics.NurbsBrepApplication.NURBS_SHAPE_FUNCTIONS))
			for node in cad_model_part.Nodes:
				control_points = node.GetValue(var)
				int_control_points = []
				for cp in control_points:
					int_control_points.append(int(cp))
					iga_model_part.AddNode(model_part.GetNode(int(cp)),0)

				element = iga_model_part.CreateNewElement(element_name, node.Id, int_control_points, properties)
				i = 1
				for iga_var, brep_var in zip(list_var_iga, list_var_nurbs):# in range(0,element_type["parameters"]["variables"].size()):#element_type["parameters"]["variables"]:
					element.SetValue(iga_var, node.GetValue(brep_var))

		for condition_id in range(0,self.project_parameters["condition_list"].size()): #project_parameters["condition_list"]:
			condition_type = self.project_parameters["condition_list"][condition_id]
			condition_name = condition_type["condition_name"].GetString()
			if(model_part.HasSubModelPart(condition_type["iga_model_part"].GetString())):
				iga_model_part = model_part.GetSubModelPart(condition_type["iga_model_part"].GetString())
			else:
				iga_model_part = model_part.CreateSubModelPart(condition_type["iga_model_part"].GetString())
			#cad_model_part_names = condition_type["cad_model_part"]
			#properties = model_part.GetProperties(condition_type["parameters"]["properties_id"].GetInt())
			for property in model_part.Properties:
				if (property.Id == condition_type["parameters"]["properties_id"].GetInt()):
					properties = property
					break

			for cad_model_part_name_id in range (0,condition_type["cad_model_part"].size()):
				cad_model_part_name = condition_type["cad_model_part"][cad_model_part_name_id]
				model = KratosMultiphysics.Model()
				model.AddModelPart(self.model_part_elements)
				cad_model_part  = model[cad_model_part_name.GetString()]

				list_var_nurbs = []
				list_var_iga = []
				for variable in range(0,condition_type["parameters"]["variables"].size()):
					list_var_nurbs.append(KratosMultiphysics.KratosGlobals.GetVariable(condition_type["parameters"]["variables"][variable]["cad_variable"].GetString()))
					list_var_iga.append(KratosMultiphysics.KratosGlobals.GetVariable(condition_type["parameters"]["variables"][variable]["iga_variable"].GetString()))

				for node in cad_model_part.Nodes:
					int_control_points = []
					for i in range(0,condition_type["parameters"]["control_points"].size()):
						control_points = node.GetValue(KratosMultiphysics.KratosGlobals.GetVariable(condition_type["parameters"]["control_points"][i].GetString()))

						for cp in control_points:
							int_control_points.append(int(cp))
							iga_model_part.AddNode(model_part.GetNode(int(cp)),0)

					condition = iga_model_part.CreateNewCondition(condition_name, node.Id, int_control_points, properties)
					for iga_var, brep_var in zip(list_var_iga, list_var_nurbs):# in range(0,element_type["parameters"]["variables"].size()):#element_type["parameters"]["variables"]:
						condition.SetValue(iga_var, node.GetValue(brep_var))
					#for variable in range(0,condition_type["parameters"]["variables"].size()):#condition_type["parameters"]["variables"]:
					#	condition.SetValue(eval(condition_type["parameters"]["variables"][variable]["iga_variable"].GetString()), node.GetValue(eval(condition_type["parameters"]["variables"][variable]["cad_variable"].GetString())))