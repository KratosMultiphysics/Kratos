import math
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

class CustomFluidDynamicsAnalysis(FluidDynamicsAnalysis):
    def __init__(self, model, project_parameters, A = 1.5, w = 2.0 * math.pi, print_output = False):
        self.A = A # Piston amplitude
        self.w = w # Piston angular frequency
        self.print_output = print_output # Print output flag
        super(CustomFluidDynamicsAnalysis, self).__init__(model,project_parameters)

    def ModifyInitialGeometry(self):
        # Call the parent ModifyInitialGeometry()
        super(CustomFluidDynamicsAnalysis,self).ModifyInitialGeometry()

        # Create and read the structure model part
        structure_model_part = self.model.CreateModelPart("Structure")
        structure_model_part.SetBufferSize(2)
        structure_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)
        structure_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        structure_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        structure_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        KratosMultiphysics.ModelPartIO("embedded_piston_2d_test_structure").ReadModelPart(structure_model_part)

        # Initialize the level-set function and VELOCITY field
        A = 1.5
        w = 2.0 * math.pi
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            dis = node.X0 - 2.5
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, dis)
            if dis > 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, [w * A, 0.0, 0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 1, [w * A, 0.0, 0.0])
                node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 2, [w * A, 0.0, 0.0])

        # Set initial condition in the negative nodes of the intersected elements as well
        for element in self._GetSolver().GetComputingModelPart().Elements:
            n_pos = 0
            n_neg = 0
            for node in element.GetNodes():
                if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                    n_neg += 1
                else:
                    n_pos += 1
            if n_pos != 0 and n_neg != 0:
                for node in element.GetNodes():
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, [w * A, 0.0, 0.0])
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 1, [w * A, 0.0, 0.0])
                    node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 2, [w * A, 0.0, 0.0])

    def ApplyBoundaryConditions(self):
        # Do the structure advance in time
        structure_model_part = self.model.GetModelPart("Structure")
        dt = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.DELTA_TIME]
        current_time = self._GetSolver().GetComputingModelPart().ProcessInfo[KratosMultiphysics.TIME]
        new_time = current_time + dt
        structure_model_part.CloneTimeStep(new_time)
        structure_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        # Move the structure
        dis = self.A * math.sin(self.w * current_time)
        vel = self.w * self.A * math.cos(self.w * current_time)
        acc = - self.w * self.w * self.A * math.sin(self.w * current_time)
        for node in self.model.GetModelPart("Structure").Nodes:
            node.X = node.X0 + dis
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, [dis,0.0,0.0])
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, [vel,0.0,0.0])
            node.SetSolutionStepValue(KratosMultiphysics.ACCELERATION, [acc,0.0,0.0])

        # Recompute the level-set function
        x_0 = 2.5
        for node in self._GetSolver().GetComputingModelPart().Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, node.X0 - x_0 - dis)

        # Apply base BCs
        # Note it is done after the level-set update to correct it with the distance modification process
        super(CustomFluidDynamicsAnalysis,self).ApplyBoundaryConditions()

    def OutputSolutionStep(self):
        if self.print_output:
            super(CustomFluidDynamicsAnalysis,self).OutputSolutionStep()

@UnitTest.skipIfApplicationsNotAvailable("MeshMovingApplication")
class EmbeddedPistonTest(UnitTest.TestCase):

    # Embedded element tests
    def testEmbeddedPiston2D(self):
        self.A = 1.5
        self.w = 2.0 * math.pi
        self.work_folder = "EmbeddedPiston2DTest"
        self.reference_file = "reference_embedded_piston_2D"
        self.settings = "EmbeddedPiston2DTestParameters.json"
        self.ExecuteEmbeddedPistonTest()

    def ExecuteEmbeddedPistonTest(self):
        self.setUp()
        self.setUpProblem()
        self.runTest()
        self.tearDown()
        self.checkResults()

    def setUp(self):
        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False
        self.work_folder = "EmbeddedPiston2DTest"
        self.reference_file = "reference_embedded_piston_2D"
        self.settings = "EmbeddedPiston2DTestParameters.json"

    def tearDown(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            filename_fluid = self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()
            filename_structure = filename_fluid + '_structure'
            KratosUtilities.DeleteFileIfExisting(filename_fluid + '.time')
            KratosUtilities.DeleteFileIfExisting(filename_structure + '.time')
            if not self.print_output:
                KratosUtilities.DeleteFileIfExisting(filename_fluid + '.post.bin')
                KratosUtilities.DeleteFileIfExisting(self.work_folder + '.post.lst')

    def setUpProblem(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            with open(self.settings, 'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())
            self.model = KratosMultiphysics.Model()
            self.simulation = CustomFluidDynamicsAnalysis(self.model, self.ProjectParameters, self.A, self.w, self.print_output)

    def runTest(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            self.simulation.Run()

    def checkResults(self):
        with UnitTest.WorkFolderScope(self.work_folder, __file__):
            results_model_part = self.simulation._GetSolver().GetComputingModelPart()
            ## 2D results check
            if (results_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                if self.print_reference_values:
                    with open(self.reference_file+'.csv','w') as ref_file:
                        ref_file.write("#ID, VELOCITY_X, VELOCITY_Y\n")
                        for node in results_model_part.Nodes:
                            vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                            ref_file.write("{0}, {1}, {2}\n".format(node.Id, vel[0], vel[1]))
                else:
                    with open(self.reference_file+'.csv','r') as reference_file:
                        reference_file.readline() # skip header
                        line = reference_file.readline()

                        for node in results_model_part.Nodes:
                            values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                            reference_vel_x = values[1]
                            reference_vel_y = values[2]

                            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                            self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                            self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)

                            line = reference_file.readline()
                        if line != '': # If we did not reach the end of the reference file
                            self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")
            ## 3D results check
            elif (results_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                if self.print_reference_values:
                    with open(self.reference_file+'.csv','w') as ref_file:
                        ref_file.write("#ID, VELOCITY_X, VELOCITY_Y, VELOCITY_Z\n")
                        for node in results_model_part.Nodes:
                            vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                            ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], vel[2]))
                else:
                    with open(self.reference_file+'.csv','r') as reference_file:
                        reference_file.readline() # skip header
                        line = reference_file.readline()

                        for node in results_model_part.Nodes:
                            values = [ float(i) for i in line.rstrip('\n ').split(',') ]
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
    UnitTest.main()
