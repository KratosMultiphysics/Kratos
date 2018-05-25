from KratosMultiphysics import *
from KratosMultiphysics.FSIApplication import *

import KratosMultiphysics.KratosUnittest as UnitTest

from os import remove

try:
    from KratosMultiphysics.StructuralMechanicsApplication import *
    from KratosMultiphysics.FluidDynamicsApplication import *
    missing_external_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)

import NonConformant_OneSideMap

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipIf(missing_external_dependencies,"Missing required application: "+ missing_application)
class NonConformantOneSideMapTest(UnitTest.TestCase):

    def setUp(self):
        # Test case settings
        self.domain_size = None
        self.fluid_input_file = None
        self.solid_input_file = None
        self.work_folder = "NonConformantOneSideMapTest"

        # Mapper settings
        self.search_radius_factor = 2.0
        self.mapper_max_iterations = 200
        self.mapper_tolerance = 1e-12

        # Testing settings
        self.check_tolerance = None
        self.print_output = False

    def tearDown(self):
        self.deleteOutFile(self.fluid_input_file+'.time')
        self.deleteOutFile(self.solid_input_file+'.time')

    def deleteOutFile(self,filename):
        with WorkFolderScope(self.work_folder):
            try:
                remove(filename)
            except FileNotFoundError as e:
                pass

    def test2D_1(self):
        self.domain_size = 2

        self.fluid_input_file = "2D_fluid_mesh_test1"
        self.solid_input_file = "2D_solid_mesh_test1"

        self.fluid_interface_name = "Fluid_interface"
        self.solid_interface_name = "Solid_interface"

        self.check_tolerance = 0.1

        def SetFluidData():
            fluid_interface = self.GetFluidInterfaceModelPart()

            # Constant pressure distribution
            for node in fluid_interface.Nodes:
                node.SetSolutionStepValue(PRESSURE, 0, 2.0)

            # Linear x-reaction distribution
            fluid_interface.Nodes[1].SetSolutionStepValue(REACTION_X, 0, 2.0)
            fluid_interface.Nodes[3].SetSolutionStepValue(REACTION_X, 0, 4.7)
            fluid_interface.Nodes[5].SetSolutionStepValue(REACTION_X, 0, 5.6)
            fluid_interface.Nodes[7].SetSolutionStepValue(REACTION_X, 0, 6.6)
            fluid_interface.Nodes[9].SetSolutionStepValue(REACTION_X, 0, 3.6)

            # Linear y-reaction distribution
            fluid_interface.Nodes[1].SetSolutionStepValue(REACTION_Y, 0, 3.6)
            fluid_interface.Nodes[3].SetSolutionStepValue(REACTION_Y, 0, 6.6)
            fluid_interface.Nodes[5].SetSolutionStepValue(REACTION_Y, 0, 5.6)
            fluid_interface.Nodes[7].SetSolutionStepValue(REACTION_Y, 0, 4.7)
            fluid_interface.Nodes[9].SetSolutionStepValue(REACTION_Y, 0, 2.0)

        def SetSolidData():
            solid_interface = self.GetSolidInterfaceModelPart()

            # Parabolic velocity distribution
            for node in solid_interface.Nodes:
                node.SetSolutionStepValue(VELOCITY_X, 0, 2*node.Y - node.Y*node.Y/2.0)
                node.SetSolutionStepValue(VELOCITY_Y, 0, 2*node.Y - node.Y*node.Y/2.0)
                node.SetSolutionStepValue(VELOCITY_Z, 0, 0.0)

        def CheckFluidResults():
            fluid_interface = self.GetFluidInterfaceModelPart()

            # Mapped VELOCITY check
            for node in fluid_interface.Nodes:
                obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
                expected_velocity_value = [2*node.Y-node.Y*node.Y/2.0, 2*node.Y-node.Y*node.Y/2.0, 0.0]
                for i in range(0,3):
                    self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=self.check_tolerance)

        def CheckSolidResults():
            solid_interface = self.GetSolidInterfaceModelPart()

            # Mapped PRESSURE check
            for node in solid_interface.Nodes:
                obtained_pressure_value = node.GetSolutionStepValue(PRESSURE,0)
                expected_pressure_value = 2.0
                self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value)

            # Force equilibrium check (mapped POINT_LOAD)
            expected_sum_fx = 22.5
            expected_sum_fy = 22.5
            expected_sum_fz = 0.0

            obtained_sum_fx = 0.0
            obtained_sum_fy = 0.0
            obtained_sum_fz = 0.0

            for node in solid_interface.Nodes:
                f_solid = node.GetSolutionStepValue(POINT_LOAD)
                obtained_sum_fx += f_solid[0]
                obtained_sum_fy += f_solid[1]
                obtained_sum_fz += f_solid[2]

            self.assertAlmostEqual(obtained_sum_fx, expected_sum_fx, delta=self.check_tolerance)
            self.assertAlmostEqual(obtained_sum_fy, expected_sum_fy, delta=self.check_tolerance)
            self.assertAlmostEqual(obtained_sum_fz, expected_sum_fz, delta=self.check_tolerance)

        # Run test
        self.RunTestCase(SetFluidData, SetSolidData, CheckFluidResults, CheckSolidResults)

    def test2D_2(self):
        self.domain_size = 2

        self.fluid_input_file = "2D_fluid_mesh_test2"
        self.solid_input_file = "2D_solid_mesh_test2"

        self.fluid_interface_name = "Fluid_interface"
        self.solid_interface_name = "Solid_interface"

        self.check_tolerance = 0.1

        def SetFluidData():
            fluid_interface = self.GetFluidInterfaceModelPart()

            for node in fluid_interface.Nodes:
                node.SetSolutionStepValue(PRESSURE, 0, node.Y*20)               # Linear PRESSURE variation
                node.SetSolutionStepValue(REACTION_X, 0, node.Y*30)             # Linear REACTION_X variation
                node.SetSolutionStepValue(REACTION_Y, 0, node.Y*node.Y*200)     # Parabolic REACTION_Y variation

        def SetSolidData():
            solid_interface = self.GetSolidInterfaceModelPart()

            for node in solid_interface.Nodes:
                node.SetSolutionStepValue(VELOCITY_X, 0, node.Y*10)             # Linear VELOCITY_X distribution
                node.SetSolutionStepValue(VELOCITY_Y, 0, node.Y*node.Y*10000)   # Parabolic VELOCITY_Y distribution
                node.SetSolutionStepValue(VELOCITY_Z, 0, 0.0)

        def CheckFluidResults():
            fluid_interface = self.GetFluidInterfaceModelPart()

            # Mapped VELOCITY check
            for node in fluid_interface.Nodes:
                obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
                expected_velocity_value = [node.Y*10, node.Y*node.Y*10000, 0.0]
                for i in range(0,3):
                    self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=0.15)

        def CheckSolidResults():
            solid_interface = self.GetSolidInterfaceModelPart()

            # Mapped PRESSURE check
            for node in solid_interface.Nodes:
                obtained_pressure_value = node.GetSolutionStepValue(PRESSURE,0)
                expected_pressure_value = node.Y*20
                self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value, delta=self.check_tolerance)

            # Force equilibrium check (mapped POINT_LOAD)
            expected_sum_fx = 765.0
            expected_sum_fy = 858.375
            expected_sum_fz = 0.0

            obtained_sum_fx = 0.0
            obtained_sum_fy = 0.0
            obtained_sum_fz = 0.0

            for node in solid_interface.Nodes:
                f_solid = node.GetSolutionStepValue(POINT_LOAD)
                obtained_sum_fx += f_solid[0]
                obtained_sum_fy += f_solid[1]
                obtained_sum_fz += f_solid[2]

            self.assertAlmostEqual(obtained_sum_fx, expected_sum_fx, delta=self.check_tolerance)
            self.assertAlmostEqual(obtained_sum_fy, expected_sum_fy, delta=self.check_tolerance)
            self.assertAlmostEqual(obtained_sum_fz, expected_sum_fz, delta=self.check_tolerance)

        # Run test
        self.RunTestCase(SetFluidData, SetSolidData, CheckFluidResults, CheckSolidResults)

    def test3D_1(self):
        self.domain_size = 3

        self.fluid_input_file = "3D_fluid_mesh_test1"
        self.solid_input_file = "3D_solid_mesh_test1"

        self.fluid_interface_name = "Fluid_interface"
        self.solid_interface_name = "Solid_interface"

        self.check_tolerance = 5e-3

        def SetFluidData():
            fluid_interface = self.GetFluidInterfaceModelPart()

            # Linear pressure distribution
            for node in fluid_interface.Nodes:
                node.SetSolutionStepValue(PRESSURE, 0, node.X*node.Y*node.Z)

            # Linear reaction distribution
            for node in fluid_interface.Nodes:
                node.SetSolutionStepValue(REACTION_X, 0,  node.Y)
                node.SetSolutionStepValue(REACTION_Y, 0, 2*node.Z)
                node.SetSolutionStepValue(REACTION_Z, 0, 3*node.X)

        def SetSolidData():
            solid_interface = self.GetSolidInterfaceModelPart()

            # Parabolic velocity distribution
            for node in solid_interface.Nodes:
                node.SetSolutionStepValue(VELOCITY_X, 0, -(2*node.Z - node.Y*node.Z/2.0))
                node.SetSolutionStepValue(VELOCITY_Y, 0, 2*node.Z - node.Y*node.Z/2.0)
                node.SetSolutionStepValue(VELOCITY_Z, 0, node.X*node.Y)

        def CheckFluidResults():
            fluid_interface = self.GetFluidInterfaceModelPart()

            # Mapped VELOCITY check
            for node in fluid_interface.Nodes:
                obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
                expected_velocity_value = [-(2*node.Z-node.Y*node.Z/2.0) , 2*node.Z-node.Y*node.Z/2.0, node.X*node.Y]
                for i in range(0,3):
                    self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=self.check_tolerance)

        def CheckSolidResults():
            solid_interface = self.GetSolidInterfaceModelPart()

            # Mapped PRESSURE check
            for node in solid_interface.Nodes:
                obtained_pressure_value = node.GetSolutionStepValue(PRESSURE,0)
                expected_pressure_value = node.X*node.Y*node.Z
                self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value, delta=self.check_tolerance)

            # Force equilibrium check (mapped POINT_LOAD)
            expected_sum_fx = 69.54492675560007
            expected_sum_fy = 187.00000000000003
            expected_sum_fz = 352.3652168853

            obtained_sum_fx = 0.0
            obtained_sum_fy = 0.0
            obtained_sum_fz = 0.0

            for node in solid_interface.Nodes:
                f_solid = node.GetSolutionStepValue(POINT_LOAD)
                obtained_sum_fx += f_solid[0]
                obtained_sum_fy += f_solid[1]
                obtained_sum_fz += f_solid[2]

            self.assertAlmostEqual(obtained_sum_fx, expected_sum_fx, delta=2.5)
            self.assertAlmostEqual(obtained_sum_fy, expected_sum_fy, delta=5.0)
            self.assertAlmostEqual(obtained_sum_fz, expected_sum_fz, delta=10.0)

        # Run test
        self.RunTestCase(SetFluidData, SetSolidData, CheckFluidResults, CheckSolidResults)

    def test3D_two_faces(self):
        self.domain_size = 3

        self.fluid_input_file = "3D_fluid_mesh_two_faces_test"
        self.solid_input_file = "3D_solid_mesh_two_faces_test"

        self.fluid_positive_interface_name = "Fluid_interface_pos"
        self.fluid_negative_interface_name = "Fluid_interface_neg"
        self.solid_interface_name = "Solid_interface"

        self.check_tolerance = 5e-3

        def SetFluidData():
            positive_interface_model_part = self.fluid_main_model_part.GetSubModelPart(self.fluid_positive_interface_name)
            negative_interface_model_part = self.fluid_main_model_part.GetSubModelPart(self.fluid_negative_interface_name)

            # Assign the PRESSURE values
            for node_pos, node_neg in zip(positive_interface_model_part.Nodes, negative_interface_model_part.Nodes):
                if (node_neg.Id == node_pos.Id):
                    node_pos.SetSolutionStepValue(PRESSURE, 0, 0.0) # If the node is shared, set a unique value
                else:
                    node_pos.SetSolutionStepValue(PRESSURE, 0, -(-36.07*node_pos.Y*node_pos.Y+36.04*node_pos.Y-8)*(8*node_pos.Z-16*node_pos.Z*node_pos.Z)) # Positive face nodal values
                    node_neg.SetSolutionStepValue(PRESSURE, 0, (-36.07*node_neg.Y*node_neg.Y+36.04*node_neg.Y-8)*(8*node_neg.Z-16*node_neg.Z*node_neg.Z))  # Negative face nodal values

        def SetSolidData():
            solid_interface = self.GetSolidInterfaceModelPart()

            # Parabolic velocity distribution
            for node in solid_interface.Nodes:
                node.SetSolutionStepValue(VELOCITY_X, 0, -(2*node.Z - node.Y*node.Z/2.0))
                node.SetSolutionStepValue(VELOCITY_Y, 0, 2*node.Z - node.Y*node.Z/2.0)
                node.SetSolutionStepValue(VELOCITY_Z, 0, node.X*node.Y)

        def CheckFluidResults():
            positive_interface_model_part = self.fluid_main_model_part.GetSubModelPart(self.fluid_positive_interface_name)
            negative_interface_model_part = self.fluid_main_model_part.GetSubModelPart(self.fluid_negative_interface_name)

            # Mapped VELOCITY check in positive interface
            for node in positive_interface_model_part.Nodes:
                obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
                expected_velocity_value = [-(2*node.Z-node.Y*node.Z/2.0) , 2*node.Z-node.Y*node.Z/2.0, node.X*node.Y]
                for i in range(0,3):
                    self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=self.check_tolerance)

            # Mapped VELOCITY check in negative interface
            for node in negative_interface_model_part.Nodes:
                obtained_velocity_value = node.GetSolutionStepValue(VELOCITY,0)
                expected_velocity_value = [-(2*node.Z-node.Y*node.Z/2.0) , 2*node.Z-node.Y*node.Z/2.0, node.X*node.Y]
                for i in range(0,3):
                    self.assertAlmostEqual(obtained_velocity_value[i], expected_velocity_value[i], delta=self.check_tolerance)

        def CheckSolidResults():
            solid_interface = self.GetSolidInterfaceModelPart()

            # Mapped PRESSURE check
            for node in solid_interface.Nodes:
                obtained_pressure_value = node.GetSolutionStepValue(POSITIVE_FACE_PRESSURE,0)
                expected_pressure_value = -(-36.07*node.Y*node.Y+36.04*node.Y-8)*(8*node.Z-16*node.Z*node.Z)
                self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value, delta=2.5e-2)

            for node in solid_interface.Nodes:
                obtained_pressure_value = node.GetSolutionStepValue(NEGATIVE_FACE_PRESSURE,0)
                expected_pressure_value = (-36.07*node.Y*node.Y+36.04*node.Y-8)*(8*node.Z-16*node.Z*node.Z)
                self.assertAlmostEqual(obtained_pressure_value, expected_pressure_value, delta=2.5e-2)

        # Run test
        self.RunTwoFacesTestCase(SetFluidData, SetSolidData, CheckFluidResults, CheckSolidResults)

    def RunTestCase(self, set_fluid_data, set_solid_data, check_fluid_results, check_solid_results):
        with WorkFolderScope(self.work_folder):
            # Problem set up
            self.SetUpProblem()

            # Set the interface flag in the fluid and solid interfaces
            self.GenerateInterface()

            # Construct the mapper object
            self.mapper = NonConformant_OneSideMap.NonConformant_OneSideMap(self.fluid_main_model_part,
                                                                            self.solid_main_model_part,
                                                                            self.search_radius_factor,
                                                                            self.mapper_max_iterations,
                                                                            self.mapper_tolerance)

            # Set the data to be mapped
            set_fluid_data()
            set_solid_data()

            # Map information between fluid and solid domains
            self.PerformMapping()

            # If required, print output
            if self.print_output:
                self.InitializeOutput()
                self.PrintOutput()
                self.FinalizeOutput()

            # Check the mapped results
            check_fluid_results()
            check_solid_results()

    def RunTwoFacesTestCase(self, set_fluid_data, set_solid_data, check_fluid_results, check_solid_results):
        with WorkFolderScope(self.work_folder):
            # Problem set up
            self.SetUpProblem()

            # Set the interface flag in the fluid and solid interfaces
            self.GenerateTwoFacesInterface()

            # Construct the mapper object
            self.mapper = NonConformant_OneSideMap.NonConformantTwoFaces_OneSideMap(self.fluid_main_model_part.GetSubModelPart(self.fluid_positive_interface_name),
                                                                                    self.fluid_main_model_part.GetSubModelPart(self.fluid_negative_interface_name),
                                                                                    self.solid_main_model_part.GetSubModelPart(self.solid_interface_name),
                                                                                    self.search_radius_factor,
                                                                                    self.mapper_max_iterations,
                                                                                    self.mapper_tolerance)

            # Set the data to be mapped
            set_fluid_data()
            set_solid_data()

            # Map information between fluid and solid domains
            self.PerformTwoFacesMapping()

            # If required, print output
            if self.print_output:
                self.InitializeOutput()
                self.PrintOutput()
                self.FinalizeOutput()

            # Check the mapped results
            check_fluid_results()
            check_solid_results()

    def SetUpProblem(self):
        # Defining a model part for the fluid and one for the structure
        self.fluid_main_model_part = ModelPart("fluid_part")
        self.solid_main_model_part = ModelPart("solid_part")

        # Set the domain size (2D or 3D test)
        self.solid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.domain_size)
        self.fluid_main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.domain_size)

        # Fluid model part variables addition
        self.fluid_main_model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(PRESSURE)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(REACTION)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(MAPPER_SCALAR_PROJECTION_RHS)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(MAPPER_VECTOR_PROJECTION_RHS)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(VAUX_EQ_TRACTION)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(NORMAL)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(NODAL_MAUX)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED)
        self.fluid_main_model_part.AddNodalSolutionStepVariable(VECTOR_PROJECTED)

        # Structure model part variables addition
        self.solid_main_model_part.AddNodalSolutionStepVariable(VELOCITY)
        self.solid_main_model_part.AddNodalSolutionStepVariable(PRESSURE)
        self.solid_main_model_part.AddNodalSolutionStepVariable(POINT_LOAD)
        self.solid_main_model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
        self.solid_main_model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
        self.solid_main_model_part.AddNodalSolutionStepVariable(MAPPER_SCALAR_PROJECTION_RHS)
        self.solid_main_model_part.AddNodalSolutionStepVariable(MAPPER_VECTOR_PROJECTION_RHS)
        self.solid_main_model_part.AddNodalSolutionStepVariable(VAUX_EQ_TRACTION)
        self.solid_main_model_part.AddNodalSolutionStepVariable(NORMAL)
        self.solid_main_model_part.AddNodalSolutionStepVariable(NODAL_MAUX)
        self.solid_main_model_part.AddNodalSolutionStepVariable(SCALAR_PROJECTED)
        self.solid_main_model_part.AddNodalSolutionStepVariable(VECTOR_PROJECTED)

        # Model parts reading
        ModelPartIO(self.fluid_input_file).ReadModelPart(self.fluid_main_model_part)
        ModelPartIO(self.solid_input_file).ReadModelPart(self.solid_main_model_part)

        # Buffer size set
        self.fluid_main_model_part.SetBufferSize(1)
        self.solid_main_model_part.SetBufferSize(1)

    def GenerateInterface(self):
        # Set fluid and structure interfaces
        for node in self.fluid_main_model_part.GetSubModelPart(self.fluid_interface_name).Nodes:
            node.Set(INTERFACE, True)
        for node in self.solid_main_model_part.GetSubModelPart(self.solid_interface_name).Nodes:
            node.Set(INTERFACE, True)

    def GenerateTwoFacesInterface(self):
        # Set fluid and structure interfaces
        for node in self.fluid_main_model_part.GetSubModelPart(self.fluid_positive_interface_name).Nodes:
            node.Set(INTERFACE, True)
        for node in self.fluid_main_model_part.GetSubModelPart(self.fluid_negative_interface_name).Nodes:
            node.Set(INTERFACE, True)
        for node in self.solid_main_model_part.GetSubModelPart(self.solid_interface_name).Nodes:
            node.Set(INTERFACE, True)

    def PerformMapping(self):
        self.mapper.FluidToStructure_ScalarMap(PRESSURE, PRESSURE,   True)
        self.mapper.FluidToStructure_VectorMap(REACTION, POINT_LOAD, True, True)
        self.mapper.StructureToFluid_VectorMap(VELOCITY, VELOCITY,   True, False)

    def PerformTwoFacesMapping(self):
        self.mapper.PositiveFluidToStructure_ScalarMap(PRESSURE, POSITIVE_FACE_PRESSURE,   True)
        self.mapper.NegativeFluidToStructure_ScalarMap(PRESSURE, NEGATIVE_FACE_PRESSURE,   True)
        self.mapper.StructureToPositiveFluid_VectorMap(VELOCITY, VELOCITY,   True, False)
        self.mapper.StructureToNegativeFluid_VectorMap(VELOCITY, VELOCITY,   True, False)

    def GetFluidInterfaceModelPart(self):
        return self.fluid_main_model_part.GetSubModelPart(self.fluid_interface_name)

    def GetSolidInterfaceModelPart(self):
        return self.solid_main_model_part.GetSubModelPart(self.solid_interface_name)

    def InitializeOutput(self):
        gid_mode = GiDPostMode.GiD_PostBinary
        multifile = MultiFileFlag.SingleFile
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteDeformed
        write_conditions = WriteConditionsFlag.WriteConditions
        mesh_name = 0.0

        self.fluid_gid_io = GidIO(self.fluid_input_file, gid_mode, multifile, deformed_mesh_flag, write_conditions)
        self.fluid_gid_io.InitializeMesh( mesh_name)
        self.fluid_gid_io.WriteMesh( self.fluid_main_model_part.GetMesh() )
        self.fluid_gid_io.FinalizeMesh()
        self.fluid_gid_io.InitializeResults(mesh_name,(self.fluid_main_model_part).GetMesh())

        self.solid_gid_io = GidIO(self.solid_input_file, gid_mode, multifile, deformed_mesh_flag, write_conditions)
        self.solid_gid_io.InitializeMesh( mesh_name)
        self.solid_gid_io.WriteMesh( self.solid_main_model_part.GetMesh() )
        self.solid_gid_io.FinalizeMesh()
        self.solid_gid_io.InitializeResults(mesh_name,(self.solid_main_model_part).GetMesh())

    def FinalizeOutput(self):
        self.fluid_gid_io.FinalizeResults()
        self.solid_gid_io.FinalizeResults()

    def PrintOutput(self):
        fluid_label = self.fluid_main_model_part.ProcessInfo[TIME]
        solid_label = self.solid_main_model_part.ProcessInfo[TIME]

        self.fluid_gid_io.WriteNodalResults(VELOCITY, self.fluid_main_model_part.Nodes, fluid_label, 0)
        self.fluid_gid_io.WriteNodalResults(PRESSURE, self.fluid_main_model_part.Nodes, fluid_label, 0)
        self.fluid_gid_io.WriteNodalResults(REACTION, self.fluid_main_model_part.Nodes, fluid_label, 0)
        self.fluid_gid_io.WriteNodalResults(NODAL_MAUX, self.fluid_main_model_part.Nodes, fluid_label, 0)
        self.fluid_gid_io.WriteNodalResults(VAUX_EQ_TRACTION, self.fluid_main_model_part.Nodes, fluid_label, 0)

        self.solid_gid_io.WriteNodalResults(VELOCITY, self.solid_main_model_part.Nodes, solid_label, 0)
        self.solid_gid_io.WriteNodalResults(PRESSURE, self.solid_main_model_part.Nodes, solid_label, 0)
        self.solid_gid_io.WriteNodalResults(POINT_LOAD, self.solid_main_model_part.Nodes, solid_label, 0)
        self.solid_gid_io.WriteNodalResults(NODAL_MAUX, self.solid_main_model_part.Nodes, solid_label, 0)
        self.solid_gid_io.WriteNodalResults(VAUX_EQ_TRACTION, self.solid_main_model_part.Nodes, solid_label, 0)

if __name__ == '__main__':
    test = NonConformantOneSideMapTest()
    test.setUp()
    test.print_output = True
    # test.test2D_1()
    # test.test2D_2()
    # test.test3D_1()
    test.test3D_two_faces()
    test.tearDown()
