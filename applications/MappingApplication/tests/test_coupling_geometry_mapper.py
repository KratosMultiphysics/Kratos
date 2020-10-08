import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
from KratosMultiphysics import KratosUnittest
from KratosMultiphysics.MappingApplication import Mapper
data_comm = KM.DataCommunicator.GetDefault()
import mapper_test_case
import os

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), "mdpa_files", file_name)

class TestCouplingGeometryMapper(KratosUnittest.TestCase):


    @classmethod
    def setUpClass(self):
        mapper_params = KM.Parameters("""{
            "mapper_type": "coupling_geometry",
            "echo_level" : 0,
            "is_precompute_mapping_matrix" : false,
			"dual_mortar": false,
			"consistency_scaling" : true,
			"modeler_name" : "MappingGeometriesModeler",
            "modeler_parameters":{
						"origin_model_part_name" : "origin",
						"destination_model_part_name" : "destination",
						"is_interface_sub_model_parts_specified" : true,
						"origin_interface_sub_model_part_name" : "origin.line_tri",
						"destination_interface_sub_model_part_name" : "destination.line_quad"
					}
        }""")

        SetupModelParts(self)
        CreateMapper(self,mapper_params)

    def test_map_displacements(self):
        reference_displacement = 1.0
        SetConstantVariable(self.interface_model_part_origin,KM.DISPLACEMENT,reference_displacement)
        self.mapper.Map(KM.DISPLACEMENT, KM.DISPLACEMENT)
        mapped_results = GetNodalVariable(self.interface_model_part_destination,KM.DISPLACEMENT)
        for nodal_vector_result in mapped_results:
            for nodal_component in nodal_vector_result:
                self.assertAlmostEqual(nodal_component,reference_displacement)




def SetupModelParts(self):
    self.model = KM.Model()
    self.model_part_origin = self.model.CreateModelPart("origin")
    self.model_part_destination = self.model.CreateModelPart("destination")

    self.model_part_origin.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    self.model_part_origin.AddNodalSolutionStepVariable(KM.FORCE)

    self.model_part_destination.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    self.model_part_destination.AddNodalSolutionStepVariable(KM.REACTION)

    origin_mdpa_file_name = "cube_tri"
    destination_mdpa_file_name = "cube_quad"

    ReadModelPart(self.model_part_origin, origin_mdpa_file_name)
    ReadModelPart(self.model_part_destination, destination_mdpa_file_name)


def ReadModelPart(model_part, mdpa_file_name):
    import_flags = KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER
    KM.ModelPartIO(GetFilePath(mdpa_file_name), import_flags).ReadModelPart(model_part)


def CreateMapper(self,mapper_parameters):
    self.mapper_type = mapper_parameters["mapper_type"].GetString()
    self.mapper_parameters = mapper_parameters.Clone()

    origin_interface_string = mapper_parameters["modeler_parameters"]["origin_interface_sub_model_part_name"].GetString()
    self.interface_model_part_origin =self.model.GetModelPart(origin_interface_string)

    dest_interface_string = mapper_parameters["modeler_parameters"]["destination_interface_sub_model_part_name"].GetString()
    self.interface_model_part_destination =self.model.GetModelPart(dest_interface_string)

    if data_comm.IsDistributed():
        self.mapper = KratosMapping.MapperFactory.CreateMPIMapper(
            self.model_part_origin, self.model_part_destination, mapper_parameters)
    else:
        self.mapper = KratosMapping.MapperFactory.CreateMapper(
            self.model_part_origin, self.model_part_destination, mapper_parameters)

def SetConstantVariable(model_part,variable,reference_value):
    for node in model_part.Nodes:
        var_x = reference_value
        var_y = reference_value
        var_z = reference_value
        node.SetSolutionStepValue(variable, KM.Vector([var_x, var_y, var_z]))

def SetLinearDisplacements(model_part):
    for node in model_part.Nodes:
        disp_x = 1.0*node.X
        disp_y = 1.0*node.Y
        disp_z = 1.0*node.Z
        node.SetSolutionStepValue(KM.DISPLACEMENT, KM.Vector([disp_x, disp_y, disp_z]))

def GetNodalVariable(model_part,variable):
    var_vector = []
    for node in model_part.Nodes:
        var_vector.append(node.GetSolutionStepValue(variable))
    return var_vector


if __name__ == '__main__':
    KratosUnittest.main()
