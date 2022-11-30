import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_analysis

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities

class TestCaseConfiguration(object):
    """Auxiliary class to configure the test

    This auxiliary class customizes the material propertiesand
    boundary conditions for each test variant.

    Public member variables:
    xmin -- Minimum x-coordinate
    xmax -- Maximum x-coordinate
    ymin -- Minimum y-coordinate
    ymax -- Maximum y-coordinate
    T_xmin -- Temperature at minimum x-coordinate
    T_xmax -- Temperature at maximum x-coordinate
    rho -- Density
    c -- Specific heat
    k -- Thermal conductivity
    ux -- Convective velocity
    source -- Source term value
    """

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

class SourceTermTestConvectionDiffusionAnalysis(convection_diffusion_analysis.ConvectionDiffusionAnalysis):
    """Derived convection-diffusion analysis stage to set the test material properties and boundary conditions."""

    def __init__(self, model, project_parameters, test_config):
        super().__init__(model, project_parameters)
        self.config = test_config

    def ModifyInitialProperties(self):
        super().ModifyInitialProperties()

        ## Set the material properties according to the test case configuration
        for node in self.model.GetModelPart("ThermalModelPart").Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, self.config.rho)
            node.SetSolutionStepValue(KratosMultiphysics.CONDUCTIVITY, self.config.k)
            node.SetSolutionStepValue(KratosMultiphysics.SPECIFIC_HEAT, self.config.c)

    def ApplyBoundaryConditions(self):
        super().ApplyBoundaryConditions()

        velocity = KratosMultiphysics.Vector(3)
        velocity[0] = self.config.ux
        velocity[1] = 0.0
        velocity[2] = 0.0

        ## Set initial and boundary conditions according to the test configuration
        for node in self.model.GetModelPart("ThermalModelPart").Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.HEAT_FLUX, self.config.source)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, velocity)
            if node.X == self.config.xmin:
                node.Fix(KratosMultiphysics.TEMPERATURE)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.config.T_xmin)
            elif node.X == self.config.xmax:
                node.Fix(KratosMultiphysics.TEMPERATURE)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, self.config.T_xmax)

class SourceTermTest(KratosUnittest.TestCase):

    def setUp(self):
        self.domain_size = 2
        self.input_file = "source_test"
        self.theta = 1.0 # Since it is steady state, use backward Euler
        # Note: Crank-Nicolson (theta=0.5) won't converge in a single iteration (or at all, for huge dt)

        self.checked_variable = KratosMultiphysics.TEMPERATURE
        self.check_tolerance = 1e-6
        self.print_output = False
        self.print_reference_values = False
        self.calculate_reactions = False

        self.config = TestCaseConfiguration()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope("SourceTermTest", __file__):
            KratosUtilities.DeleteFileIfExisting(self.input_file+'.time')

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

    def testReaction(self):
        self.reference_file = "reaction_test"
        self.config.T_xmin = 0.0
        self.config.T_xmax = 0.0
        self.config.source = 1.0
        self.config.ux = 0.0
        self.checked_variable = KratosMultiphysics.REACTION_FLUX
        self.calculate_reactions = True

        self.testSourceTerm()

    def testSourceTerm(self):
        with KratosUnittest.WorkFolderScope("SourceTermTest", __file__):
            ## Set up the custom source term input
            self.model = KratosMultiphysics.Model()
            with open("SourceTermTestProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
                parameters["solver_settings"]["compute_reactions"].SetBool(self.calculate_reactions)
                parameters["solver_settings"]["transient_parameters"]["theta"].SetDouble(self.theta)

            ## Run test
            self.source_term_analysis = SourceTermTestConvectionDiffusionAnalysis(self.model, parameters, self.config)
            self.source_term_analysis.Run()

            ## Check results
            self.checkResults()

            ## If required, print output
            if self.print_output:
                self.printOutput()

    def checkResults(self):
        model_part = self.model.GetModelPart("ThermalModelPart")
        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                ref_file.write("#ID, {0}\n".format(self.checked_variable.Name()))
                for node in model_part.Nodes:
                    value = node.GetSolutionStepValue(self.checked_variable,0)
                    ref_file.write("{0}, {1}\n".format(node.Id, value))
        else:
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                line = reference_file.readline()

                for node in model_part.Nodes:
                    values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                    node_id = values[0]
                    reference_value = values[1]

                    value = node.GetSolutionStepValue(self.checked_variable)
                    self.assertAlmostEqual(reference_value, value, delta=self.check_tolerance)

                    line = reference_file.readline()
                if line != '': # If we did not reach the end of the reference file
                    self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

    def printOutput(self):
        gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
        multifile = KratosMultiphysics.MultiFileFlag.SingleFile
        deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions =KratosMultiphysics. WriteConditionsFlag.WriteElementsOnly
        gid_io = KratosMultiphysics.GidIO(self.input_file, gid_mode, multifile, deformed_mesh_flag, write_conditions)

        mesh_name = 0.0
        model_part = self.model.GetModelPart("ThermalModelPart")
        gid_io.InitializeMesh(mesh_name)
        gid_io.WriteMesh(model_part.GetMesh())
        gid_io.FinalizeMesh()
        gid_io.InitializeResults(mesh_name, model_part.GetMesh())

        label = model_part.ProcessInfo[KratosMultiphysics.TIME]
        gid_io.WriteNodalResults(KratosMultiphysics.VELOCITY, model_part.Nodes,label,0)
        gid_io.WriteNodalResults(KratosMultiphysics.TEMPERATURE, model_part.Nodes,label,0)
        gid_io.WriteNodalResults(KratosMultiphysics.DENSITY, model_part.Nodes,label,0)
        gid_io.WriteNodalResults(KratosMultiphysics.CONDUCTIVITY, model_part.Nodes,label,0)
        gid_io.WriteNodalResults(KratosMultiphysics.SPECIFIC_HEAT, model_part.Nodes,label,0)
        gid_io.WriteNodalResults(KratosMultiphysics.HEAT_FLUX, model_part.Nodes,label,0)
        gid_io.WriteNodalResults(KratosMultiphysics.REACTION_FLUX, model_part.Nodes,label,0)

        gid_io.FinalizeResults()

if __name__ == '__main__':
    a = SourceTermTest()
    a.setUp()
    a.print_reference_values = False
    a.print_output = True
    #a.testPureDiffusion()
    #a.testConvectionDominated()
    #a.testDiffusionDominated()
    a.testReaction()
    a.tearDown()
