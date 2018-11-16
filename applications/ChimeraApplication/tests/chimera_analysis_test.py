import KratosMultiphysics as km
import KratosMultiphysics.FluidDynamicsApplication as kfd
import KratosMultiphysics.ChimeraApplication as kchim
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError:
    have_external_solvers = False

from fluid_chimera_analysis import FluidChimeraAnalysis

import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.kratos_utilities as kratos_utilities

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, exc_type, exc_value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class ChimeraAnalysisTest(UnitTest.TestCase):

    def setUp(self):
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        self.print_reference_values = False
        self.reference_file = "reference_chimera_monolithic_simple_test"

    def testMonolithic(self):
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        self.print_reference_values = False
        self.reference_file = "reference_chimera_monolithic_simple_test"
        work_folder = "chimera_monolithic_simple_test"
        settings_file_name = "test_chimera_monolithic_simple_ProjectParameters.json"
        with WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
            kratos_utilities.DeleteFileIfExisting("test_chimera_monolithic_simple_input.time")
            kratos_utilities.DeleteFileIfExisting("ModifiedPatch_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("GENERIC_domainboundary_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("ModifiedPatchBoundary_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("chimera_monolithic_simple_test.post.lst")
            kratos_utilities.DeleteFileIfExisting("test_chimera_monolithic_simple_input_0.post.msh")
            kratos_utilities.DeleteFileIfExisting("test_chimera_monolithic_simple_input_0.post.res")
            #self._checkResults()
    
    def testFractionalStep(self):
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        self.print_reference_values = False
        self.reference_file = "reference_chimera_fractionalstep_simple_test"
        work_folder = "chimera_fractionalstep_simple_test"
        settings_file_name = "test_chimera_fractionalstep_simple_ProjectParameters.json"
        with WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
            kratos_utilities.DeleteFileIfExisting("test_chimera_fractionalstep_simple_input.time")
            kratos_utilities.DeleteFileIfExisting("ModifiedPatch_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("GENERIC_domainboundary_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("ModifiedPatchBoundary_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("chimera_fractionalstep_simple_test.post.lst")
            kratos_utilities.DeleteFileIfExisting("test_chimera_fractionalstep_simple_input_0.post.msh")
            kratos_utilities.DeleteFileIfExisting("test_chimera_fractionalstep_simple_input_0.post.res")
            
    def testMultipleOverlappingPatch(self):
        self.check_tolerance = 1e-6
        # Set to true to get post-process files for the test
        self.print_reference_values = False
        self.reference_file = "reference_chimera_multiple_overlapping_patch_simple_test"
        work_folder = "chimera_multiple_overlapping_patch"
        settings_file_name = "test_chimera_multiple_overlapping_simple_ProjectParameters.json"
        with WorkFolderScope(work_folder):
            self._run_test(settings_file_name)
            kratos_utilities.DeleteFileIfExisting("test_chimera_multiple_overlapping_simple_input.time")
            kratos_utilities.DeleteFileIfExisting("ModifiedPatch_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("GENERIC_domainboundary_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("ModifiedPatchBoundary_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("NoSlip2D_structure_1_1_0_0.vtk")
            kratos_utilities.DeleteFileIfExisting("chimera_multiple_overlapping_patch.post.lst")
            kratos_utilities.DeleteFileIfExisting("test_chimera_multiple_overlapping_simple_0.post.msh")
            kratos_utilities.DeleteFileIfExisting("test_chimera_multiple_overlapping_simple_0.post.res")
    
    def _run_test(self,settings_file_name):
        self.model = km.Model()
        with open(settings_file_name,'r') as settings_file:
            settings = km.Parameters(settings_file.read())

        analysis = FluidChimeraAnalysis(self.model,settings)
        analysis.Run()

        print("Finished!!!!!!!")
        print("started to check results")

        print("started to check results")
        if self.print_reference_values:
            with open(self.reference_file+'.csv','w') as ref_file:
                print("started to check results  2" )
                ref_file.write("#ID, VELOCITY_X, VELOCITY_Y\n")
                print("started to check results  3" )
                for node in self.model["FluidModelPart"].Nodes:
                    vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY,0)
                    ref_file.write("{0}, {1}, {2}\n".format(node.Id, vel[0], vel[1]))
        else:
            print("started to check results")
            with open(self.reference_file+'.csv','r') as reference_file:
                reference_file.readline() # skip header
                print("started to check results1")
                line = reference_file.readline()
                print("started to check results12")
                for node in self.model["FluidModelPart"].Nodes:
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

        #with open(self.reference_file+'.csv','r') as reference_file:
            #reference_file.readline() # skip header
            #print("started to check results1")
            #line = reference_file.readline()
            #print("started to check results12")
            #for node in self.model["MainModelPart"].Nodes:
                #print("inside")
                #values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                #node_id = values[0]
                #reference_vel_x = values[1]
                #reference_vel_y = values[2]

                #velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
                #self.assertAlmostEqual(reference_vel_x, velocity[0], delta = self.check_tolerance)
                #self.assertAlmostEqual(reference_vel_y, velocity[1], delta = self.check_tolerance)
                #print("test")

                #line = reference_file.readline()
            #if line != '': # If we did not reach the end of the reference file
                #self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")


    #def _checkResults(self):
        



if __name__ == '__main__':
    test_case = ChimeraAnalysisTest()
    #test_case.setUp()
    test_case.testMonolithic()
    test_case.testFractionalStep()
    #test_case.testMultipleOverlappingPatch()
    print("completed all tests")


