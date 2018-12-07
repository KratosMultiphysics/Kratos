from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication as ExternalSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os
import math

class TestMortarMapping(KratosUnittest.TestCase):
    def setUp(self):
        pass

    def __base_test_mapping(self, input_filename, num_nodes, pure_implicit):
        KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

        self.model = KratosMultiphysics.Model()
        self.main_model_part = self.model.CreateModelPart("Structure", 2)

        ## Creation of the Kratos model (build sub_model_parts or submeshes)
        self.StructureModel = {"Structure": self.main_model_part}

        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_GAP)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)

        self.main_model_part.CloneTimeStep(1.01)

        KratosMultiphysics.ModelPartIO(input_filename).ReadModelPart(self.main_model_part)

        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TEMPERATURE, self.main_model_part)

        if (self.main_model_part.HasSubModelPart("Contact")):
            interface_model_part = self.main_model_part.GetSubModelPart("Contact")
        else:
            interface_model_part = self.main_model_part.CreateSubModelPart("Contact")

        self.mapping_model_part = self.main_model_part.GetSubModelPart("DISPLACEMENT_Displacement_Auto2")

        self.model_part_slave = self.main_model_part.GetSubModelPart("Parts_Parts_Auto1")
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLAVE, True, self.model_part_slave.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.MASTER, False, self.model_part_slave.Nodes)
        self.model_part_master = self.main_model_part.GetSubModelPart("Parts_Parts_Auto2")
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.SLAVE, False, self.model_part_master.Nodes)
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.MASTER, True, self.model_part_master.Nodes)

        for prop in self.main_model_part.GetProperties():
            prop[ContactStructuralMechanicsApplication.INTEGRATION_ORDER_CONTACT] = 3

        self.main_model_part.ProcessInfo[ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR] = 3.0e-1

        for node in self.mapping_model_part.Nodes:
            node.Set(KratosMultiphysics.INTERFACE, True)

        Preprocess = ContactStructuralMechanicsApplication.InterfacePreprocessCondition(self.main_model_part)

        interface_parameters = KratosMultiphysics.Parameters("""{"simplify_geometry": false}""")
        Preprocess.GenerateInterfacePart3D(self.mapping_model_part, interface_parameters)

        # We compute NODAL_H that can be used in the search and some values computation
        find_nodal_h = KratosMultiphysics.FindNodalHProcess(self.mapping_model_part)
        find_nodal_h.Execute()

        # We copy the conditions to the ContactSubModelPart
        for cond in self.mapping_model_part.Conditions:
            interface_model_part.AddCondition(cond)
        del(cond)
        for node in self.mapping_model_part.Nodes:
            interface_model_part.AddNode(node, 0)
        del(node)

        # We initialize the conditions
        alm_init_var = ContactStructuralMechanicsApplication.ALMFastInit(self.mapping_model_part)
        alm_init_var.Execute()

        search_parameters = KratosMultiphysics.Parameters("""
        {
            "search_factor"               : 3.5,
            "allocation_size"             : 1000,
            "check_gap"                   : "NoCheck",
            "type_search"                 : "InRadius"
        }
        """)
        if (num_nodes == 3):
            contact_search = ContactStructuralMechanicsApplication.TreeContactSearch3D3N(self.main_model_part, search_parameters)
        else:
            contact_search = ContactStructuralMechanicsApplication.TreeContactSearch3D4N(self.main_model_part, search_parameters)

        # We initialize the search utility
        contact_search.CreatePointListMortar()
        contact_search.InitializeMortarConditions()
        contact_search.UpdateMortarConditions()

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

        # Copy respective conditions
        for cond in self.main_model_part.Conditions:
            if (cond.Is(KratosMultiphysics.SLAVE)):
                self.model_part_slave.AddCondition(cond)
            if (cond.Is(KratosMultiphysics.MASTER)):
                self.model_part_master.AddCondition(cond)

        if (pure_implicit == True):
            #linear_solver = ExternalSolversApplication.SuperLUSolver()
            linear_solver = KratosMultiphysics.SkylineLUFactorizationSolver()

            if (num_nodes == 3):
                self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D3NDouble(self.model_part_master, self.model_part_slave,  KratosMultiphysics.TEMPERATURE, map_parameters, linear_solver)
                self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D3NVector(self.model_part_master, self.model_part_slave,  KratosMultiphysics.DISPLACEMENT, map_parameters, linear_solver)
            else:
                self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D4NDouble(self.model_part_master, self.model_part_slave,  KratosMultiphysics.TEMPERATURE, map_parameters, linear_solver)
                self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D4NVector(self.model_part_master, self.model_part_slave,  KratosMultiphysics.DISPLACEMENT, map_parameters, linear_solver)
        else:
            if (num_nodes == 3):
                self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D3NDouble(self.model_part_master, self.model_part_slave,  KratosMultiphysics.TEMPERATURE, map_parameters)
                self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D3NVector(self.model_part_master, self.model_part_slave,  KratosMultiphysics.DISPLACEMENT, map_parameters)
            else:
                self.mortar_mapping_double = KratosMultiphysics.SimpleMortarMapperProcess3D4NDouble(self.model_part_master, self.model_part_slave,  KratosMultiphysics.TEMPERATURE, map_parameters)
                self.mortar_mapping_vector = KratosMultiphysics.SimpleMortarMapperProcess3D4NVector(self.model_part_master, self.model_part_slave,  KratosMultiphysics.DISPLACEMENT, map_parameters)

    def _mapper_tests(self, input_filename, num_nodes, pure_implicit = False):

        self.__base_test_mapping(input_filename, num_nodes, pure_implicit)

        self.mortar_mapping_double.Execute()
        self.mortar_mapping_vector.Execute()

        ### DEBUG
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

        check_parameters["input_file_name"].SetString(input_filename+".json")

        check = from_json_check_result_process.FromJsonCheckResultProcess(self.StructureModel, check_parameters)
        check.ExecuteInitialize()
        check.ExecuteBeforeSolutionLoop()
        check.ExecuteFinalizeSolutionStep()

        #import json_output_process

        #out_parameters = KratosMultiphysics.Parameters("""
        #{
            #"output_variables"     : ["TEMPERATURE","DISPLACEMENT"],
            #"output_file_name"     : "",
            #"model_part_name"      : "Structure",
            #"sub_model_part_name"  : "Parts_Parts_Auto1"
        #}
        #""")

        #out_parameters["output_file_name"].SetString(input_filename+".json")

        #out = json_output_process.JsonOutputProcess(self.StructureModel, out_parameters)
        #out.ExecuteInitialize()
        #out.ExecuteBeforeSolutionLoop()
        #out.ExecuteFinalizeSolutionStep()

    def test_less_basic_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/integration_tests/test_integration_several_triangles"
        self._mapper_tests(input_filename, 3, False)

    def test_simple_curvature_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/integration_tests/test_simple_curvature"
        self._mapper_tests(input_filename, 3, False)

    def test_mortar_mapping_triangle(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/integration_tests/test_double_curvature_integration_triangle"
        self._mapper_tests(input_filename, 3, False)

    def test_mortar_mapping_quad(self):
        input_filename = os.path.dirname(os.path.realpath(__file__)) + "/integration_tests/test_double_curvature_integration_quadrilateral"
        self._mapper_tests(input_filename, 4, False)

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
                                                    "WriteConditionsFlag": "WriteElementsOnly",
                                                    "MultiFileFlag": "SingleFile"
                                                },
                                                "nodal_results"       : ["DISPLACEMENT","NORMAL","TEMPERATURE"],
                                                "nodal_nonhistorical_results": ["NODAL_AREA","NODAL_MAUX","NODAL_VAUX"],
                                                "nodal_flags_results": ["MASTER","SLAVE"]
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
