import KratosMultiphysics as Kratos
from KratosMultiphysics import Node, Parameters

class ParallelepipedRegularMesher:
    def __init__(self,
                 model_part_to_be_filled,
                 lower_corner_coordinates,
                 higher_corner_coordinates,
                 number_of_divisions_per_dimension,
                 element_name='Element3D4N',
                 condition_name='WallCondition3D'):
        self.lc = lower_corner_coordinates
        self.hc = higher_corner_coordinates
        self.n_divisions = number_of_divisions_per_dimension
        self.mp = model_part_to_be_filled

        parameters = Parameters("{}")
        parameters.AddEmptyValue("element_name").SetString(element_name)
        parameters.AddEmptyValue("condition_name").SetString(condition_name)
        parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        parameters.AddEmptyValue("number_of_divisions").SetVector(self.n_divisions)

        self.domain_geometry = Kratos.Hexahedra3D8(Node(1, self.hc[0], self.hc[1], self.lc[2]),
                                          Node(2, self.lc[0], self.hc[1], self.lc[2]),
                                          Node(3, self.lc[0], self.lc[1], self.lc[2]),
                                          Node(4, self.hc[0], self.lc[1], self.lc[2]),
                                          Node(5, self.hc[0], self.hc[1], self.hc[2]),
                                          Node(6, self.lc[0], self.hc[1], self.hc[2]),
                                          Node(7, self.lc[0], self.lc[1], self.hc[2]),
                                          Node(8, self.hc[0], self.lc[1], self.hc[2]))

        self.mesh_generator_process = Kratos.StructuredMeshGeneratorProcess(self.domain_geometry,
                                                                        self.mp,
                                                                        parameters)
    def FillModelPartWithNewMesh(self):
        n_elements = len(self.mp.Elements)
        if n_elements > 0:
            raise Exception('The model_part already contains elements ( '
                            + str(n_elements)
                            + ' in total). Regular mesh cannot be generated.)')
        else:
            if not self.mp.HasSubModelPart("Skin"):
                self.mp.CreateSubModelPart("Skin")

            self.mesh_generator_process.Execute()

class RectangularRegularMesher:
    def __init__(self,
                 model_part_to_be_filled,
                 lower_corner_coordinates,
                 higher_corner_coordinates,
                 number_of_divisions_per_dimension,
                 element_name='Element2D3N',
                 condition_name='WallCondition2D'):
        self.lc = lower_corner_coordinates
        self.hc = higher_corner_coordinates
        self.n_divisions = number_of_divisions_per_dimension
        self.mp = model_part_to_be_filled

        parameters = Parameters("{}")
        parameters.AddEmptyValue("element_name").SetString(element_name)
        parameters.AddEmptyValue("condition_name").SetString(condition_name)
        parameters.AddEmptyValue("create_skin_sub_model_part").SetBool(False)
        parameters.AddEmptyValue("number_of_divisions").SetVector(self.n_divisions)

        self.domain_geometry = Kratos.Quadrilateral2D4(Node(1, self.lc[0], self.hc[1], 0),
                                          Node(2, self.hc[0], self.hc[1], 0),
                                          Node(3, self.hc[0], self.lc[1], 0),
                                          Node(4, self.lc[0], self.lc[1], 0))

        self.mesh_generator_process = Kratos.StructuredMeshGeneratorProcess(self.domain_geometry,
                                                                        self.mp,
                                                                        parameters)
    def FillModelPartWithNewMesh(self):
        n_elements = len(self.mp.Elements)
        if n_elements > 0:
            raise Exception('The model_part already contains elements ( '
                            + str(n_elements)
                            + ' in total). Regular mesh cannot be generated.)')
        else:
            if not self.mp.HasSubModelPart("Skin"):
                self.mp.CreateSubModelPart("Skin")

            self.mesh_generator_process.Execute()
