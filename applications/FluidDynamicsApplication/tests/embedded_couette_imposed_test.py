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
class EmbeddedCouetteImposedTest(UnitTest.TestCase):

    # Embedded element tests
    def testEmbeddedCouetteImposed2D(self):
        self.slip_flag = False
        self.distance = 1.65
        self.embedded_velocity = 2.56
        self.work_folder = "EmbeddedCouetteImposed2DTest"
        self.reference_file = "reference_couette_imposed_2D"
        self.settings = "EmbeddedCouetteImposed2DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def testEmbeddedCouetteImposed3D(self):
        self.slip_flag = False
        self.distance = 1.65
        self.embedded_velocity = 2.56
        self.work_folder = "EmbeddedCouetteImposed3DTest"
        self.reference_file = "reference_couette_imposed_3D"
        self.settings = "EmbeddedCouetteImposed3DTestParameters.json"
        self.ExecuteEmbeddedCouetteTest()

    def ExecuteEmbeddedCouetteTest(self):
        with WorkFolderScope(self.work_folder):
            self.setUp()
            self.setUpProblem()
            self.setUpDistanceField()
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
            if self.slip_flag and self.solver.settings.Has("is_slip"):
                self.ProjectParameters["solver_settings"]["is_slip"].SetBool(True)

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

            ## Processes initialization
            for process in self.list_of_processes:
                process.ExecuteInitialize()

            self.main_model_part = self.model.GetModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())

    def setUpDistanceField(self):
        # Set the distance function
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if (domain_size == 2):
            for node in self.main_model_part.Nodes:
                distance = self.distance - node.Y
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance)
        elif (domain_size == 3):
            for node in self.main_model_part.Nodes:
                distance = self.distance - node.Z
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, distance)

        # Set the ELEMENTAL_VELOCITY in the intersected elements
        for elem in self.main_model_part.Elements:
            n_pos = 0
            n_neg = 0
            for node in elem.GetNodes():
                if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                    n_neg += 1
                else:
                    n_pos += 1

            if ((n_pos != 0) and (n_neg != 0)):
                embedded_velocity = KratosMultiphysics.Vector([self.embedded_velocity,0.0,0.0])
                elem.SetValue(KratosMultiphysics.EMBEDDED_VELOCITY,embedded_velocity)

    def setUpNoSlipBoundaryConditions(self):
        # Set the inlet function
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if (domain_size == 2):
            for node in self.main_model_part.GetSubModelPart("Inlet").Nodes:
                dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if (dist > 0.0):
                    vel_x = (self.embedded_velocity/self.distance)*node.Y
                    aux_vel = KratosMultiphysics.Vector([vel_x,0.0,0.0])
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
                else:
                    aux_vel = KratosMultiphysics.Vector([0.0,0.0,0.0])
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)
        else:
            for node in self.main_model_part.GetSubModelPart("Inlet").Nodes:
                dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if (dist > 0.0):
                    vel_x = (self.embedded_velocity/self.distance)*node.Z
                    aux_vel = KratosMultiphysics.Vector([vel_x,0.0,0.0])
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
                else:
                    aux_vel = KratosMultiphysics.Vector([0.0,0.0,0.0])
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, aux_vel)
                node.Fix(KratosMultiphysics.VELOCITY_X)
                node.Fix(KratosMultiphysics.VELOCITY_Y)
                node.Fix(KratosMultiphysics.VELOCITY_Z)

    def setUpNoSlipInitialCondition(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if (domain_size == 2):
            for node in self.main_model_part.Nodes:
                dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if (dist > 0.0):
                    vel_x = (self.embedded_velocity/self.distance)*node.Y
                    init_v = KratosMultiphysics.Vector([vel_x,0.0,0.0])
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, init_v)
        else:
            for node in self.main_model_part.Nodes:
                dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                if (dist > 0.0):
                    vel_x = (self.embedded_velocity/self.distance)*node.Z
                    init_v = KratosMultiphysics.Vector([vel_x,0.0,0.0])
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
            step = 0

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
    test = EmbeddedCouetteImposedTest()
    test.setUp()
    test.distance = 1.65
    test.embedded_velocity = 2.56
    test.print_output = False
    test.print_reference_values = False
    test.work_folder = "EmbeddedCouetteImposed2DTest"
    test.reference_file = "reference_couette_imposed_2D"
    test.settings = "EmbeddedCouetteImposed2DTestParameters.json"
    test.setUpProblem()
    test.setUpDistanceField()
    test.setUpNoSlipInitialCondition()
    test.setUpNoSlipBoundaryConditions()
    test.runTest()
    test.tearDown()
    test.checkResults()
