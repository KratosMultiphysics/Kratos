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

class TestCaseConfiguration(object):
    def __init__(self):
        self.xmin = 0.0
        self.xmax = 10.0
        self.ymin = 0.0
        self.ymax = 1.0

        self.rho = 1.0
        self.c = 1.0
        self.k = 1.0

        self.ux = 0.0
        self.source = 0.0
        

class SourceTermTest(UnitTest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.input_file = "source_test"
        self.reference_file = "source_reference"

        self.dt = 1e10 # This is steady state test
        self.nsteps = 1
        self.theta = 1.0 # Since it is steady state, use backward Euler
        # Note: Crank-Nicolson (theta=0.5) won't converge in a single iteration (or at all, for huge dt)

        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False

        self.config = TestCaseConfiguration()

    def tearDown(self):
        import os
        with WorkFolderScope("SourceTermTest"):
            try:
                os.remove(self.input_file+'.time')
            except FileNotFoundError as e:
                pass

    def testPureDiffusion(self):
        self.reference_file = "pure_diffusion"
        self.config.T_xmin = 0.0
        self.config.T_xmax = 0.0
        self.config.source = 1.0
        self.config.ux = 0.0

        # What follows is the analytical solution
        #def F(x):
        #    L = self.config.xmax - self.config.xmin
        #    return self.config.source * x * (L-x) / 2.
        #print( [ F(float(i)) for i in range(0, 10) ] )

        self.testSourceTerm()

    def testDiffusionDominated(self):
        self.reference_file = "diffusion_dominated"
        self.config.T_xmin = 0.0
        self.config.T_xmax = 0.0
        self.config.ux = 0.1
        self.config.source = self.config.ux * self.config.rho / (self.config.xmax - self.config.xmin)

        # What follows is the analytical solution
        #def F(x):
        #    import math
        #    L = self.config.xmax - self.config.xmin
        #    a = self.config.rho*self.config.ux / self.config.k
        #    return x / L - (1.0 - math.exp( a*x ) ) / ( 1.0 - math.exp( a*L ) ) 
        #print( [ F(float(i)) for i in range(0, 10) ] )            
        
        self.testSourceTerm()


    def testConvectionDominated(self):
        self.reference_file = "convection_dominated"
        self.config.T_xmin = 0.0
        self.config.T_xmax = 0.0
        self.config.ux = 2.0
        self.config.source = self.config.ux * self.config.rho / (self.config.xmax - self.config.xmin)

        # What follows is the analytical solution
        #def F(x):
        #    import math
        #    L = self.config.xmax - self.config.xmin
        #    a = self.config.rho*self.config.ux / self.config.k
        #    return x / L - (1.0 - math.exp( a*x ) ) / ( 1.0 - math.exp( a*L ) ) 
        #print( [ F(float(i)) for i in range(0, 10) ] )  
        
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
        self.thermal_solver.theta = self.theta
        self.thermal_solver.Initialize()


    def setUpProblem(self):

        velocity = Array3()
        velocity[0] = self.config.ux
        velocity[1] = 0.0
        velocity[2] = 0.0

        ## Set initial and boundary conditions
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(DENSITY,self.config.rho)
            node.SetSolutionStepValue(CONDUCTIVITY,self.config.k)
            node.SetSolutionStepValue(SPECIFIC_HEAT,self.config.c)
            node.SetSolutionStepValue(HEAT_FLUX,self.config.source)
            node.SetSolutionStepValue(VELOCITY,velocity)

            if node.X == self.config.xmin:
                node.Fix(TEMPERATURE)
                node.SetSolutionStepValue(TEMPERATURE,self.config.T_xmin)
            elif node.X == self.config.xmax:
                node.Fix(TEMPERATURE)
                node.SetSolutionStepValue(TEMPERATURE,self.config.T_xmax)

    def runTest(self):
        time = 0.0

        for step in range(self.nsteps):
            time = time+self.dt
            self.model_part.CloneTimeStep(time)
            self.thermal_solver.Solve()

    def checkResults(self):

        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, TEMPERATURE\n")
                for node in self.model_part.Nodes:
                    temp = node.GetSolutionStepValue(TEMPERATURE,0)
                    ref_file.write("{0}, {1}\n".format(node.Id, temp))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()
                node_iter = self.model_part.Nodes

                for node in self.model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    node_id = values[0]
                    reference_temp = values[1]

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
        gid_io.WriteNodalResults(TEMPERATURE,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(DENSITY,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(CONDUCTIVITY,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(SPECIFIC_HEAT,self.model_part.Nodes,label,0)
        gid_io.WriteNodalResults(HEAT_FLUX,self.model_part.Nodes,label,0)

        gid_io.FinalizeResults()

if __name__ == '__main__':
    a = SourceTermTest()
    a.setUp()
    a.testPureDiffusion()
    #a.testConvectionDominated()
    #a.testDiffusionDominated()
    a.tearDown()