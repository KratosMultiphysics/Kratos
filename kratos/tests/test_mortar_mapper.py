from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.basic_mapping_process import BasicMappingProcess

import os
import math

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMortarMapperCore(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def __base_test_mapping(self, input_filename, num_nodes, master_num_nodes, pure_implicit, inverted, discontinuous, origin_are_conditions, destination_are_conditions):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
        self.model = KratosMultiphysics.Model()

        self.main_model_part = self.model.CreateModelPart("Main",2)

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        self.main_model_part.CloneTimeStep(1.01)

        KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

        ## DEBUG Generate discontinous case
        #self.__generate_discontinous_case(inverted)

        if inverted:
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

        # Double parameters
        double_map_parameters = KratosMultiphysics.Parameters("""
        {
            "origin_model_part_name"           : "please_specify_model_part_name",
            "destination_model_part_name"      : "please_specify_model_part_name",
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2,
            "origin_variable"                  : "TEMPERATURE",
            "discontinuous_interface"          : false,
            "origin_are_conditions"            : true,
            "destination_are_conditions"       : true
        }
        """)
        double_map_parameters["discontinuous_interface"].SetBool(discontinuous)
        double_map_parameters["origin_are_conditions"].SetBool(origin_are_conditions)
        double_map_parameters["destination_are_conditions"].SetBool(destination_are_conditions)
        double_map_parameters["origin_model_part_name"].SetString(self.model_part_master.Name)
        double_map_parameters["destination_model_part_name"].SetString(self.model_part_slave.Name)

        if pure_implicit:
            linear_solver_settings = KratosMultiphysics.Parameters("""
            {
                "solver_type" : "skyline_lu_factorization"
            }
            """)
            double_map_parameters.AddValue("linear_solver_settings", linear_solver_settings)

        self.mortar_mapping_double = BasicMappingProcess(self.model, double_map_parameters)

        # Vector parameters
        vector_map_parameters = KratosMultiphysics.Parameters("""
        {
            "origin_model_part_name"           : "please_specify_model_part_name",
            "destination_model_part_name"      : "please_specify_model_part_name",
            "echo_level"                       : 0,
            "absolute_convergence_tolerance"   : 1.0e-9,
            "relative_convergence_tolerance"   : 1.0e-4,
            "max_number_iterations"            : 10,
            "integration_order"                : 2,
            "origin_variable"                  : "DISPLACEMENT",
            "discontinuous_interface"          : false,
            "origin_are_conditions"            : true,
            "destination_are_conditions"       : true
        }
        """)
        vector_map_parameters["discontinuous_interface"].SetBool(discontinuous)
        vector_map_parameters["origin_are_conditions"].SetBool(origin_are_conditions)
        vector_map_parameters["destination_are_conditions"].SetBool(destination_are_conditions)
        vector_map_parameters["origin_model_part_name"].SetString(self.model_part_master.Name)
        vector_map_parameters["destination_model_part_name"].SetString(self.model_part_slave.Name)

        if pure_implicit:
            linear_solver_settings = KratosMultiphysics.Parameters("""
            {
                "solver_type" : "skyline_lu_factorization"
            }
            """)
            vector_map_parameters.AddValue("linear_solver_settings", linear_solver_settings)

        self.mortar_mapping_vector = BasicMappingProcess(self.model, vector_map_parameters)

    def _mapper_tests(self, input_filename, num_nodes, master_num_nodes, pure_implicit = False, inverted = False, discontinuous = False, origin_are_conditions = True, destination_are_conditions = True):

        self.__base_test_mapping(input_filename, num_nodes, master_num_nodes, pure_implicit, inverted, discontinuous, origin_are_conditions, destination_are_conditions)

        self.mortar_mapping_double.ExecuteInitializeSolutionStep()
        self.mortar_mapping_vector.ExecuteInitializeSolutionStep()

        # Debug postprocess file
        #self.__post_process()

        import from_json_check_result_process

        check_parameters = KratosMultiphysics.Parameters("""
        {
            "check_variables"      : ["TEMPERATURE","DISPLACEMENT"],
            "input_file_name"      : "",
            "model_part_name"      : "Main",
            "sub_model_part_name"  : "Parts_Parts_Auto1"
        }
        """)

        if inverted:
            check_parameters["input_file_name"].SetString(input_filename+"_inverted.json")
        else:
            check_parameters["input_file_name"].SetString(input_filename+".json")

        check = from_json_check_result_process.FromJsonCheckResultProcess(self.model, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

        ## The following is used to create the solution database
        #import json_output_process

        #out_parameters = KratosMultiphysics.Parameters("""
        #{
            #"output_variables"     : ["TEMPERATURE","DISPLACEMENT"],
            #"output_file_name"     : "",
            #"model_part_name"      : "Main",
            #"sub_model_part_name"  : "Parts_Parts_Auto1"
        #}
        #""")

        #if inverted:
            #out_parameters["output_file_name"].SetString(input_filename+"_inverted.json")
        #else:
            #out_parameters["output_file_name"].SetString(input_filename+".json")

        #out = json_output_process.JsonOutputProcess(self.model, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()

    def test_less_basic_mortar_mapping_triangle_pure_implicit(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unnitest/mortar_mapper_python_tests/test_integration_several_triangles"
        self._mapper_tests(input_filename, 3, 3, True)

    def test_less_basic_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unnitest/mortar_mapper_python_tests/test_integration_several_triangles"
        self._mapper_tests(input_filename, 3, 3)

    def test_simple_curvature_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unnitest/mortar_mapper_python_tests/test_simple_curvature"
        self._mapper_tests(input_filename, 3, 3)

    def test_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unnitest/mortar_mapper_python_tests/test_double_curvature_integration_triangle"
        self._mapper_tests(input_filename, 3, 3)

    def test_mortar_mapping_triangle_discontinous_interface(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unnitest/mortar_mapper_python_tests/test_double_curvature_integration_triangle_discontinous_interface"
        self._mapper_tests(input_filename, 3, 3, False, False, True)

    def test_mortar_mapping_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unnitest/mortar_mapper_python_tests/test_double_curvature_integration_quadrilateral"
        self._mapper_tests(input_filename, 4, 4)

    def test_mortar_mapping_quad_tri(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unnitest/mortar_mapper_python_tests/test_double_curvature_integration_triangle_quadrilateral"
        self._mapper_tests(input_filename, 4, 3, False, False, False, False, True)

    def test_mortar_mapping_tri_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/auxiliar_files_for_python_unnitest/mortar_mapper_python_tests/test_double_curvature_integration_triangle_quadrilateral"
        self._mapper_tests(input_filename, 3, 4, False, True, False, True, False)

    def __post_process(self, debug = "GiD"):

        if debug == "GiD":
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
        elif debug == "VTK":
            from vtk_output_process import VtkOutputProcess
            self.vtk_output_process = VtkOutputProcess(self.model,
                                        KratosMultiphysics.Parameters("""{
                                                "model_part_name"                    : "Main",
                                                "nodal_solution_step_data_variables" : ["DISPLACEMENT","NORMAL","TEMPERATURE"],
                                                "nodal_data_value_variables": ["NODAL_AREA","NODAL_MAUX","NODAL_VAUX"]
                                            }
                                            """)
                                        )

            self.vtk_output_process.ExecuteInitialize()
            self.vtk_output_process.ExecuteBeforeSolutionLoop()
            self.vtk_output_process.ExecuteInitializeSolutionStep()
            self.vtk_output_process.PrintOutput()
            self.vtk_output_process.ExecuteFinalizeSolutionStep()
            self.vtk_output_process.ExecuteFinalize()

    def __generate_discontinous_case(self, inverted):
        counter_nodes = 0
        for node in self.main_model_part.Nodes:
            counter_nodes += 1

        counter_conditions = 0
        for cond in self.main_model_part.Conditions:
            counter_conditions += 1

        if inverted:
            self.model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")
            self.model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        else:
            self.model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
            self.model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")

        for cond in self.model_part_slave.Conditions:

            counter_conditions += 1
            length = cond.GetGeometry().Area() * 0.01
            list_nodes = []
            for node in cond.GetNodes():
                counter_nodes += 1
                list_nodes.append(counter_nodes)
                self.model_part_slave.CreateNewNode(counter_nodes, node.X + length, node.Y - length, node.Z)
                node.Set(KratosMultiphysics.TO_ERASE)
            cond.Set(KratosMultiphysics.TO_ERASE)

            self.model_part_slave.CreateNewCondition("SurfaceCondition3D3N", counter_conditions, list_nodes, self.main_model_part.GetProperties()[1])

        self.main_model_part.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)
        self.main_model_part.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

        # Debug
        #self.__post_process()

        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_double_curvature_integration_triangle_discontinous_interface"), KratosMultiphysics.IO.WRITE)
        model_part_io.WriteModelPart(self.main_model_part)

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
