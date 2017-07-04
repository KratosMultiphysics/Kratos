from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *

import KratosMultiphysics.KratosUnittest as UnitTest

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class SourceTermTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.input_file = "source_test"
        self.reference_file = "source_reference"

        self.convection_diffusion_solver = "eulerian"
        self.dt = 1e10 # This is steady state test
        self.nsteps = 3

        self.check_tolerance = 1e-6
        self.print_output = True
        self.print_reference_values = True

    def tearDown(self):
        import os
        with WorkFolderScope("SourceTermTest"):
            try:
                os.remove(self.input_file+'.time')
            except FileNotFoundError as e:
                pass

    def testEulerian(self):
        #self.reference_file = "reference10_eulerian"
        self.testSourceTerm()

    def testSourceTerm(self):

        with WorkFolderScope("SourceTermTest"):
            self.setUpModel()
            self.setUpSolvers()
            self.setUpProblem()

            self.runTest()

            self.checkResults()
            if self.print_output:
                self.printOutput()

    def setUpModel(self):

        self.model_part = ModelPart("TestModelPart")

        thermal_settings = ConvectionDiffusionSettings()
        thermal_settings.SetUnknownVariable(TEMPERATURE)
        thermal_settings.SetDensityVariable(DENSITY)
        thermal_settings.SetSpecificHeatVariable(SPECIFIC_HEAT)
        thermal_settings.SetDiffusionVariable(CONDUCTIVITY)
        thermal_settings.SetVolumeSourceVariable(HEAT_FLUX)
        #thermal_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
        thermal_settings.SetVelocityVariable(VELOCITY)
        thermal_settings.SetMeshVelocityVariable(MESH_VELOCITY)
        #thermal_settings.SetProjectionVariable(PROJECTED_SCALAR1)

        self.model_part.ProcessInfo.SetValue(CONVECTION_DIFFUSION_SETTINGS,thermal_settings)

    def setUpSolvers(self):

        import convection_diffusion_solver as thermal_solver
        thermal_solver.AddVariables(self.model_part)

        thermal_solver.AddVariables(self.model_part)

        model_part_io = ModelPartIO(self.input_file)
        model_part_io.ReadModelPart(self.model_part)

        self.model_part.SetBufferSize(2)
        thermal_solver.AddDofs(self.model_part)

        # thermal solver
        self.thermal_solver = thermal_solver.ConvectionDiffusionSolver(self.model_part,self.domain_size)
        self.thermal_solver.Initialize()


    def setUpProblem(self):
        xmin = 0.0
        xmax = 10.0
        ymin = 0.0
        ymax = 1.0

        rho = 1.0
        c = 1.0
        k = 1.0
        ux = 1.0
        velocity = Array3()
        velocity[0] = ux
        velocity[1] = 0.0
        velocity[2] = 0.0

        source = ux*rho / (xmax-xmin)

        ## Set initial and boundary conditions
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(DENSITY,rho)
            node.SetSolutionStepValue(CONDUCTIVITY,k)
            node.SetSolutionStepValue(SPECIFIC_HEAT,c)

            if node.X == xmin:
                node.Fix(TEMPERATURE)
                node.SetSolutionStepValue(TEMPERATURE,T_xmin)
            elif node.X == xmax:
                node.Fix(TEMPERATURE)
                node.SetSolutionStepValue(TEMPERATURE,T_xmax)

    def runTest(self):
        time = 0.0

        for step in range(self.nsteps):
            time = time+self.dt
            self.model_part.CloneTimeStep(time)
            self.buoyancy_process.ExecuteInitializeSolutionStep()
            self.fluid_solver.Solve()
            self.thermal_solver.Solve()

    def checkResults(self):

        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, VELOCITY_X, VELOCITY_Y, TEMPERATURE\n")
                for node in self.model_part.Nodes:
                    vel = node.GetSolutionStepValue(VELOCITY,0)
                    temp = node.GetSolutionStepValue(TEMPERATURE,0)
                    ref_file.write("{0}, {1}, {2}, {3}\n".format(node.Id, vel[0], vel[1], temp))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()
                node_iter = self.model_part.Nodes

                for node in self.model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    node_id = values[0]
                    reference_vel_x = values[1]
                    reference_vel_y = values[2]
                    reference_temp = values[3]

                    velocity = node.GetSolutionStepValue(VELOCITY)
                    self.assertAlmostEqual(reference_vel_x, velocity[0], delta=self.check_tolerance)
                    self.assertAlmostEqual(reference_vel_y, velocity[1], delta=self.check_tolerance)
                    temperature = node.GetSolutionStepValue(TEMPERATURE)
                    self.assertAlmostEqual(reference_temp, temperature, delta=self.check_tolerance)

                    line = reference_file.readline()
                if line != '': # If we did not reach the end of the reference file
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

    def printOutput(self):
        gid_mode = GiDPostMode.GiD_PostBinary
        multifile = MultiFileFlag.SingleFile
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteElementsOnly
        gid_io = GidIO(self.input_file,gid_mode,multifile,deformed_mesh_flag, write_conditions)

        mesh_name = 0.0
        gid_io.InitializeMesh( mesh_name)
        gid_io.WriteMesh( self.model_part.GetMesh() )
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name,(self.model_part).GetMesh())

        label = self.model_part.ProcessInfo[TIME]
        gid_io.WriteNodalResults(VELOCITY,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(PRESSURE,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(TEMPERATURE,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(DENSITY,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(VISCOSITY,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(CONDUCTIVITY,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(SPECIFIC_HEAT,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(BODY_FORCE,self.model_part.Nodes,label,0)

        gid_io.FinalizeResults()

if __name__ == '__main__':
    a = SourceTermTest()
    a.setUp()
    a.testEulerian()
    a.tearDown()