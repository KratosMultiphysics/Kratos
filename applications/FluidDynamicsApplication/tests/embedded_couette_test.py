import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError as e:
    have_external_solvers = False


import KratosMultiphysics.KratosUnittest as UnitTest

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class EmbeddedCouetteTest(UnitTest.TestCase):

    # Embedded element tests
    def testEmbeddedCouette2D(self):
        self.distance = 0.25
        self.slip_flag = False
        self.work_folder = "EmbeddedCouette2DTest"
        self.reference_file = "reference_couette_embedded_2D"
        self.settings = "EmbeddedCouette2DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedCouette3D(self):
        self.distance = 0.25
        self.slip_flag = False
        self.work_folder = "EmbeddedCouette3DTest"
        self.reference_file = "reference_couette_embedded_3D"
        self.settings = "EmbeddedCouette3DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedSlipCouette2D(self):
        self.distance = 0.25
        self.slip_flag = True
        self.work_folder = "EmbeddedCouette2DTest"
        self.reference_file = "reference_couette_embedded_slip_2D"
        self.settings = "EmbeddedCouette2DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedSlipCouette3D(self):
        self.distance = 0.25
        self.slip_flag = True
        self.work_folder = "EmbeddedCouette3DTest"
        self.reference_file = "reference_couette_embedded_slip_3D"
        self.settings = "EmbeddedCouette3DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    # Embedded Ausas element tests
    def testEmbeddedAusasCouette2D(self):
        self.distance = 0.25
        self.slip_flag = True
        self.work_folder = "EmbeddedCouette2DTest"
        self.reference_file = "reference_couette_ausas_2D"
        self.settings = "EmbeddedAusasCouette2DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedAusasCouette3D(self):
        self.distance = 0.25
        self.slip_flag = True
        self.work_folder = "EmbeddedCouette3DTest"
        self.reference_file = "reference_couette_ausas_3D"
        self.settings = "EmbeddedAusasCouette3DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    # Embedded development element tests
    def testEmbeddedDevelopmentCouette2D(self):
        self.distance = 0.25
        self.slip_flag = False
        self.work_folder = "EmbeddedCouette2DTest"
        self.reference_file = "reference_couette_development_2D"
        self.settings = "EmbeddedDevelopmentCouette2DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedSlipDevelopmentCouette2D(self):
        self.distance = 0.25
        self.slip_flag = True
        self.work_folder = "EmbeddedCouette2DTest"
        self.reference_file = "reference_couette_development_slip_2D"
        self.settings = "EmbeddedDevelopmentCouette2DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedDevelopmentCouette3D(self):
        self.distance = 0.25
        self.slip_flag = False
        self.work_folder = "EmbeddedCouette3DTest"
        self.reference_file = "reference_couette_development_3D"
        self.settings = "EmbeddedDevelopmentCouette3DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedSlipDevelopmentCouette3D(self):
        self.distance = 0.25
        self.slip_flag = True
        self.work_folder = "EmbeddedCouette3DTest"
        self.reference_file = "reference_couette_development_slip_3D"
        self.settings = "EmbeddedDevelopmentCouette3DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    # Embedded Ausas development element tests
    def testEmbeddedAusasDevelopmentCouette2D(self):
        self.distance = 0.25
        self.slip_flag = True
        self.work_folder = "EmbeddedCouette2DTest"
        self.reference_file = "reference_couette_ausas_development_2D"
        self.settings = "EmbeddedAusasDevelopmentCouette2DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedAusasDevelopmentCouette3D(self):
        self.distance = 0.25
        self.slip_flag = True
        self.work_folder = "EmbeddedCouette3DTest"
        self.reference_file = "reference_couette_ausas_development_3D"
        self.settings = "EmbeddedAusasDevelopmentCouette3DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def ExecuteEmbeddedCouetteTest(self):
        with WorkFolderScope(self.work_folder):
            self.setUp()
            self.setUpProblem()
            self.setUpDistanceField()
            if (self.slip_flag):
                self.setUpSlipInitialCondition()
                self.setUpSlipBoundaryConditions()
            else:
                self.setUpNoSlipInitialCondition()
                self.setUpNoSlipBoundaryConditions()
            self.runTest()
            self.tearDown()
            self.checkResults()

    def setUp(self):
        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

    def tearDown(self):
        with WorkFolderScope(self.work_folder):
            try:
                os.remove(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')
            except FileNotFoundError as e:
                pass

    def setUpProblem(self):
        with WorkFolderScope(self.work_folder):
            with open(self.settings, 'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.model = KratosMultiphysics.Model()

            ## Solver construction
            import python_solvers_wrapper_fluid
            self.solver = python_solvers_wrapper_fluid.CreateSolver(self.model, self.ProjectParameters)

            ## Set the "is_slip" field in the json settings (to avoid duplication it is set to false in all tests)
            if self.slip_flag:
                self.solver.settings["formulation"]["is_slip"].SetBool(True)

            self.solver.AddVariables()

            ## Read the model - note that SetBufferSize is done here
            self.solver.ImportModelPart()
            self.solver.PrepareModelPart()

            ## Add AddDofs
            self.solver.AddDofs()

            ## Solver initialization
            self.solver.Initialize()

            ## Processes construction
            import process_factory
            self.list_of_processes  = process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["gravity"] )
            self.list_of_processes += process_factory.KratosProcessFactory(self.model).ConstructListOfProcesses( self.ProjectParameters["boundary_conditions_process_list"] )

            self.main_model_part = self.model.GetModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

            ## Processes initialization
            for process in self.list_of_processes:
                process.ExecuteInitialize()

    def setUpDistanceField(self):
        # Set the distance function
        if (self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            for node in self.main_model_part.Nodes:
                distance = node.Y-self.distance
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance)
        elif (self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
            for node in self.main_model_part.Nodes:
                distance = node.Z-self.distance
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance)

        # Set the ELEMENTAL_DISTANCES value
        n_nodes = len(self.main_model_part.Elements[1].GetNodes())
        for element in self.main_model_part.Elements:
            elem_dist = KratosMultiphysics.Vector(n_nodes)
            elem_nodes = element.GetNodes()
            for i_node in range(0,n_nodes):
                elem_dist[i_node] = elem_nodes[i_node].GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            element.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, elem_dist)

    def setUpSlipBoundaryConditions(self):
        # Set the inlet function
        for node in self.main_model_part.GetSubModelPart("Inlet").Nodes:
            if (node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0):
                aux_vel = KratosMultiphysics.Vector([1.0,0.0,0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
            else:
                aux_vel = KratosMultiphysics.Vector([0.0,0.0,0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)

        # Set the SLIP elemental flag (only used in Winter's formulation)
        if (self.slip_flag):
            for element in self.main_model_part.Elements:
                element.Set(KratosMultiphysics.SLIP, True)

    def setUpSlipInitialCondition(self):
        for node in self.main_model_part.Nodes:
            if (node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0.0):
                init_v = KratosMultiphysics.Vector([1.0,0.0,0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, init_v)

    def setUpNoSlipBoundaryConditions(self):
        # Set the inlet function
        for node in self.main_model_part.GetSubModelPart("Inlet").Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if (dist > 0.0):
                aux_vel = KratosMultiphysics.Vector([dist,0.0,0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
            else:
                aux_vel = KratosMultiphysics.Vector([0.0,0.0,0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)

        # Set and fix the top boundary velocity
        for node in self.main_model_part.GetSubModelPart("Top").Nodes:
            aux_vel = KratosMultiphysics.Vector([node.GetSolutionStepValue(KratosMultiphysics.DISTANCE),0.0,0.0])
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)

    def setUpNoSlipInitialCondition(self):
        for node in self.main_model_part.Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if (dist > 0.0):
                init_v = KratosMultiphysics.Vector([dist,0.0,0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, init_v)

    def runTest(self):
        with WorkFolderScope(self.work_folder):
            if (self.print_output):
                gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
                multifile = KratosMultiphysics.MultiFileFlag.SingleFile
                deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
                write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
                gid_io = KratosMultiphysics.GidIO(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString(),gid_mode,multifile,deformed_mesh_flag, write_conditions)

                mesh_name = 0.0
                gid_io.InitializeMesh( mesh_name)
                gid_io.WriteMesh( self.main_model_part.GetMesh() )
                gid_io.FinalizeMesh()
                gid_io.InitializeResults(mesh_name,(self.main_model_part).GetMesh())

            end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

            time = 0.0

            for process in self.list_of_processes:
                process.ExecuteBeforeSolutionLoop()

            while(time <= end_time):

                time = self.solver.AdvanceInTime(time)

                for process in self.list_of_processes:
                    process.ExecuteInitializeSolutionStep()

                self.solver.InitializeSolutionStep()
                self.solver.Predict()
                self.solver.SolveSolutionStep()
                self.solver.FinalizeSolutionStep()

                for process in self.list_of_processes:
                    process.ExecuteFinalizeSolutionStep()

                for process in self.list_of_processes:
                    process.ExecuteBeforeOutputStep()

                if (self.print_output):
                    gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY,self.main_model_part.Nodes,time,0)
                    gid_io.WriteNodalResults(KratosMultiphysics.PRESSURE,self.main_model_part.Nodes,time,0)
                    gid_io.WriteNodalResults(KratosMultiphysics.DISTANCE,self.main_model_part.Nodes,time,0)

                for process in self.list_of_processes:
                    process.ExecuteAfterOutputStep()

            for process in self.list_of_processes:
                process.ExecuteFinalize()

            if (self.print_output):
                gid_io.FinalizeResults()

    def checkResults(self):
        with WorkFolderScope(self.work_folder):
            ## 2D results check
            if (self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                if self.print_reference_values:
                    with open(self.reference_file+'.csv','w') as ref_file:
                        ref_file.write("#ID, VELOCITY_X, VELOCITY_Y\n")
                        for node in self.main_model_part.Nodes:
                            vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                            ref_file.write("{0}, {1}, {2}\n".format(node.Id, vel[0], vel[1]))
                else:
                    with open(self.reference_file+'.csv','r') as reference_file:
                        reference_file.readline() # skip header
                        line = reference_file.readline()

                        for node in self.main_model_part.Nodes:
                            values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                            node_id = values[0]
                            reference_vel_x = values[1]
                            reference_vel_y = values[2]

                            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                            self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                            self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)

                            line = reference_file.readline()
                        if line != '': # If we did not reach the end of the reference file
                            self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")
            ## 3D results check
            elif (self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                if self.print_reference_values:
                    with open(self.reference_file+'.csv','w') as ref_file:
                        ref_file.write("#ID, VELOCITY_X, VELOCITY_Y, VELOCITY_Z\n")
                        for node in self.main_model_part.Nodes:
                            vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                            ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], vel[2]))
                else:
                    with open(self.reference_file+'.csv','r') as reference_file:
                        reference_file.readline() # skip header
                        line = reference_file.readline()

                        for node in self.main_model_part.Nodes:
                            values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                            node_id = values[0]
                            reference_vel_x = values[1]
                            reference_vel_y = values[2]
                            reference_vel_z = values[3]

                            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                            self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                            self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)
                            self.assertAlmostEqual(reference_vel_z, velocity[2], delta = self.check_tolerance)

                            line = reference_file.readline()
                        if line != '': # If we did not reach the end of the reference file
                            self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

if __name__ == '__main__':
    test = EmbeddedCouetteTest()
    test.setUp()
    test.distance = 0.25
    test.slip_flag = False
    test.print_output = True
    test.print_reference_values = True
    test.work_folder = "EmbeddedCouette3DTest"
    test.reference_file = "reference_couette_development_3D"
    test.settings = "EmbeddedDevelopmentCouette3DTestParameters.json"
    test.setUpProblem()
    test.setUpDistanceField()
    if (test.slip_flag):
        test.setUpSlipInitialCondition()
        test.setUpSlipBoundaryConditions()
    else:
        test.setUpNoSlipInitialCondition()
        test.setUpNoSlipBoundaryConditions()
    test.runTest()
    test.tearDown()
    test.checkResults()