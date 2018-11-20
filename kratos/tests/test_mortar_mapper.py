from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import os
import math

class TestMortarMapperCore(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def __base_test_mapping(self, input_filename, num_nodes, master_num_nodes, pure_implicit, inverted):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        self.StructureModel = KratosMultiphysics.Model()

        self.main_model_part = self.StructureModel.CreateModelPart("Structure",2)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        self.main_model_part.CloneTimeStep(1.01)

        KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TEMPERATURE, self.main_model_part)

        if (inverted is True):
            self.model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")
            self.model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        else:
            self.model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
            self.model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")

        for node in self.model_part_master.Nodes:
            x = node.X
            y = node.Y
            z = node.Z
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, z)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, x)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, y)
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, z)
        del(node)

        map_parameters = KratosMultiphysics.Parameters("""
        {
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2
        }
        """)

        if (pure_implicit == True):
            #linear_solver = ExternalSolversApplication.SuperLUSolver()
            linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()

            if (num_nodes == 3):
                if (master_num_nodes == 3):
                    self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D3NDouble(self.model_part_master, self.model_part_slave, KratosMultiphysics.TEMPERATURE, map_parameters, linear_solver)
                    self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D3NVector(self.model_part_master, self.model_part_slave, KratosMultiphysics.DISPLACEMENT, map_parameters, linear_solver)
                else:
                    self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D3N4NDouble(self.model_part_master, self.model_part_slave, KratosMultiphysics.TEMPERATURE, map_parameters, linear_solver)
                    self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D3N4NVector(self.model_part_master, self.model_part_slave, KratosMultiphysics.DISPLACEMENT, map_parameters, linear_solver)
            else:
                if (master_num_nodes == 4):
                    self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D4NDouble(self.model_part_master, self.model_part_slave, KratosMultiphysics.TEMPERATURE, map_parameters, linear_solver)
                    self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D4NVector(self.model_part_master, self.model_part_slave, KratosMultiphysics.DISPLACEMENT, map_parameters, linear_solver)
                else:
                    self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D4N3NDouble(self.model_part_master, self.model_part_slave, KratosMultiphysics.TEMPERATURE, map_parameters, linear_solver)
                    self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D4N3NVector(self.model_part_master, self.model_part_slave, KratosMultiphysics.DISPLACEMENT, map_parameters, linear_solver)
        else:
            if (num_nodes == 3):
                if (master_num_nodes == 3):
                    self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D3NDouble(self.model_part_master, self.model_part_slave, KratosMultiphysics.TEMPERATURE, map_parameters)
                    self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D3NVector(self.model_part_master, self.model_part_slave, KratosMultiphysics.DISPLACEMENT, map_parameters)
                else:
                    self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D3N4NDouble(self.model_part_master, self.model_part_slave, KratosMultiphysics.TEMPERATURE, map_parameters)
                    self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D3N4NVector(self.model_part_master, self.model_part_slave, KratosMultiphysics.DISPLACEMENT, map_parameters)
            else:
                if (master_num_nodes == 4):
                    self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D4NDouble(self.model_part_master, self.model_part_slave, KratosMultiphysics.TEMPERATURE, map_parameters)
                    self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D4NVector(self.model_part_master, self.model_part_slave, KratosMultiphysics.DISPLACEMENT, map_parameters)
                else:
                    self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D4N3NDouble(self.model_part_master, self.model_part_slave, KratosMultiphysics.TEMPERATURE, map_parameters)
                    self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D4N3NVector(self.model_part_master, self.model_part_slave, KratosMultiphysics.DISPLACEMENT, map_parameters)

    def _mapper_tests(self, input_filename, num_nodes, master_num_nodes, pure_implicit = False, inverted = False):

        self.__base_test_mapping(input_filename, num_nodes, master_num_nodes, pure_implicit, inverted)

        self.mortar_mapping_double.Execute()
        self.mortar_mapping_vector.Execute()

        # Debug postprocess file
        #self.__post_process()

        import from_json_check_result_process

        check_parameters = KratosMultiphysics.Parameters("""
        {
            "check_variables"      : ["TEMPERATURE","DISPLACEMENT"],
            "input_file_name"      : "",
            "model_part_name"      : "Structure",
            "sub_model_part_name"  : "Parts_Parts_Auto1"
        }
        """)

        if (inverted is True):
            check_parameters["input_file_name"].SetString(input_filename+"_inverted.json")
        else:
            check_parameters["input_file_name"].SetString(input_filename+".json")

        check = from_json_check_result_process.FromJsonCheckResultProcess(self.StructureModel, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

        ## The following is used to create the solution database
        #import json_output_process

        #out_parameters = KratosMultiphysics.Parameters("""
        #{
            #"output_variables"     : ["TEMPERATURE","DISPLACEMENT"],
            #"output_file_name"     : "",
            #"model_part_name"      : "Structure",
            #"sub_model_part_name"  : "Parts_Parts_Auto1"
        #}
        #""")

        #if (inverted is True):
            #out_parameters["output_file_name"].SetString(input_filename+"_inverted.json")
        #else:
            #out_parameters["output_file_name"].SetString(input_filename+".json")

        #out = json_output_process.JsonOutputProcess(self.StructureModel, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()

    def test_less_basic_mortar_mapping_triangle_pure_implicit(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/mortar_mapper_python_tests/test_integration_several_triangles"
        self._mapper_tests(input_filename, 3, 3, True)

    def test_less_basic_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files/mortar_mapper_python_tests/test_integration_several_triangles"
        self._mapper_tests(input_filename, 3, 3)

    def test_simple_curvature_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files/mortar_mapper_python_tests/test_simple_curvature"
        self._mapper_tests(input_filename, 3, 3)

    def test_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files/mortar_mapper_python_tests/test_double_curvature_integration_triangle"
        self._mapper_tests(input_filename, 3, 3)

    def test_mortar_mapping_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files/mortar_mapper_python_tests/test_double_curvature_integration_quadrilateral"
        self._mapper_tests(input_filename, 4, 4)

    def test_mortar_mapping_quad_tri(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files/mortar_mapper_python_tests/test_double_curvature_integration_triangle_quadrilateral"
        self._mapper_tests(input_filename, 4, 3)

    def test_mortar_mapping_tri_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files/mortar_mapper_python_tests/test_double_curvature_integration_triangle_quadrilateral"
        self._mapper_tests(input_filename, 3, 4, False, True)

    def __post_process(self):
        from gid_output_process import GiDOutputProcess
        self.gid_output = GiDOutputProcess(self.main_model_part,
                                    "gid_output",
                                    KratosMultiphysics.Parameters("""
                                        {
                                            "result_file_configuration" : {
                                                "gidpost_flags": {
                                                    "GiDPostMode": "GiD_PostBinary",
                                                    "WriteDeformedMeshFlag": "WriteUndeformed",
                                                    "WriteConditionsFlag": "WriteConditionsOnly",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["DISPLACEMENT","NORMAL","TEMPERATURE"],
                                                "nodal_nonhistorical_results": ["NODAL_AREA","NODAL_MAUX","NODAL_VAUX"]
                                            }
                                        }
                                        """)
                                    )

        self.gid_output.ExecuteInitialize()
        self.gid_output.ExecuteBeforeSolutionLoop()
        self.gid_output.ExecuteInitializeSolutionStep()
        self.gid_output.PrintOutput()
        self.gid_output.ExecuteFinalizeSolutionStep()
        self.gid_output.ExecuteFinalize()

    def __sci_str(self, x):
        from decimal import Decimal
        s = 10*Decimal(str(x))
        s = ('{:.' + str(len(s.normalize().as_tuple().digits) - 1) + 'E}').format(s)
        s = s.replace('E+','D0')
        s = s.replace('E-','D0-')
        s = s.replace('.','')
        if s.startswith('-'):
            return '-.' + s[1:]
        else:
            return '.' + s

if __name__ == '__main__':
    KratosUnittest.main()
