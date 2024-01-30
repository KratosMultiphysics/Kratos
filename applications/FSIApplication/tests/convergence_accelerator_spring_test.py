import KratosMultiphysics
import KratosMultiphysics.FSIApplication as KratosFSI

import KratosMultiphysics.KratosUnittest as KratosUnittest

import os, glob
from KratosMultiphysics.FSIApplication import convergence_accelerator_factory

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/AcceleratorSpringTests/" + fileName

class ConvergenceAcceleratorSpringTest(KratosUnittest.TestCase):

    def constant_force(self,model_part,guess,k,reference_z):
        f = KratosMultiphysics.Vector(3*len(model_part.Nodes))
        for i,node in enumerate(model_part.Nodes):
            j = 3*i
            position = guess[j+2] + node.Z
            if node.Z > reference_z:
                delta = position - reference_z
            else:
                delta = reference_z - position
            f[j]   = 0.0
            f[j+1] = 0.0
            f[j+2] = k*delta
        return f

    def variable_stiffness(self,model_part,guess,k,reference_z):
        f = KratosMultiphysics.Vector(3*len(model_part.Nodes))
        for i,node in enumerate(model_part.Nodes):
            j = 3*i
            position = guess[j+2] + node.Z
            if node.Z > reference_z:
                delta = position - reference_z
            else:
                delta = reference_z - position
            f[j]   = 0.0
            f[j+1] = 0.0
            f[j+2] = k * ( 1 + node.X*node.Y ) * delta
        return f

    def ComputeResidual(self,model_part,guess,force1,force2):

        f = force1(model_part,guess)
        g = force2(model_part,guess)

        if self.print_gid_output:
            for i,node in enumerate(model_part.Nodes):
                j = 3*i
                node.SetSolutionStepValue(KratosMultiphysics.FSI_INTERFACE_RESIDUAL_X,0, f[j]  -g[j]   )
                node.SetSolutionStepValue(KratosMultiphysics.FSI_INTERFACE_RESIDUAL_Y,0, f[j+1]-g[j+1] )
                node.SetSolutionStepValue(KratosMultiphysics.FSI_INTERFACE_RESIDUAL_Z,0, f[j+2]-g[j+2] )

        return f - g

    def InitializeGuess(self,model_part):
        guess = KratosMultiphysics.Vector(3*len(model_part.Nodes))
        for i,node in enumerate(model_part.Nodes):
            j = 3*i
            guess[j]   = 0.0
            guess[j+1] = 0.0
            guess[j+2] = 0.0

        return guess

    def PrintGuess(self,model_part,guess):
        if self.print_gid_output:
            for i,node in enumerate(model_part.Nodes):
                j = 3*i
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X,0, guess[j]   )
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y,0, guess[j+1] )
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0, guess[j+2] )

    def ComputeResidualNorm(self,residual):
        norm = 0.0
        for r in residual:
            norm += r*r

        return norm**0.5

    def ReadModelPart(self,filename):
        self.model = KratosMultiphysics.Model()
        model_part = self.model.CreateModelPart("Test ModelPart")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_RESIDUAL)

        model_part_io = KratosMultiphysics.ModelPartIO( filename )
        model_part_io.ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        return model_part

    def setUp(self):
        self.print_gid_output = False
        self.accelerator_tolerance = 1e-10
        self.accelerator_iterations = 50
        self.assert_delta = 1e-7

        self.model_part = self.ReadModelPart(GetFilePath("box_fluid"))

    def tearDown(self):
        if self.print_gid_output:
            local_nodes = self.model_part.GetCommunicator().LocalMesh().Nodes

            gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
            multifile = KratosMultiphysics.MultiFileFlag.SingleFile
            deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
            write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

            gid_io = KratosMultiphysics.GidIO(GetFilePath("box_fluid"),
                                              gid_mode, multifile,
                                              deformed_mesh_flag, write_conditions)

            gid_io.InitializeMesh(0)
            gid_io.WriteMesh(self.model_part.GetMesh())
            gid_io.FinalizeMesh()

            gid_io.InitializeResults(0.0,self.model_part.GetMesh())
            gid_io.WriteNodalResults(KratosMultiphysics.PARTITION_INDEX,self.model_part.Nodes,0.0,0)
            gid_io.WriteNodalResults(KratosMultiphysics.DISPLACEMENT,self.model_part.Nodes,0.0,0)
            gid_io.WriteNodalResults(KratosMultiphysics.FSI_INTERFACE_RESIDUAL,local_nodes,0.0,0)
            gid_io.FinalizeResults()

        # clean temporary files
        for f in glob.glob(GetFilePath('*.time')):
            os.remove(f)

    def _test_accelerator(self,force1,force2,solution,accelerator_settings):

        print("")
        print("Testing accelerator: ",accelerator_settings["solver_type"].GetString())

        # Construct the accelerator strategy
        coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(accelerator_settings)

        top_part = self.model_part.GetSubModelPart("Top")

        coupling_utility.Initialize()
        coupling_utility.InitializeSolutionStep()

        nl_it = 0
        x_guess = self.InitializeGuess(top_part)
        residual = self.ComputeResidual(top_part,x_guess,force1,force2)
        res_norm = self.ComputeResidualNorm(residual)

        while (nl_it <= self.accelerator_iterations):

            if res_norm > self.accelerator_tolerance:
                coupling_utility.InitializeNonLinearIteration()
                coupling_utility.UpdateSolution(residual, x_guess)
                coupling_utility.FinalizeNonLinearIteration()
            else:
                coupling_utility.FinalizeSolutionStep()
                break

            nl_it += 1
            residual = self.ComputeResidual(top_part,x_guess,force1,force2)
            res_norm = self.ComputeResidualNorm(residual)

        # Check the obtained solution
        expected_x = solution(top_part)

        if self.print_gid_output:
            self.PrintGuess(top_part,x_guess)

        for i in range(len(expected_x)):
            self.assertAlmostEqual(expected_x[i],x_guess[i],delta=self.assert_delta)

    # MVQN recursive accelerator test
    def _test_mvqn_recursive_accelerator(self,force1,force2,solution):

        mvqn_recursive_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "MVQN_recursive",
            "w_0": 0.5,
            "buffer_size": 5
        }""")

        self._test_accelerator(force1,force2,solution,mvqn_recursive_settings)

    def test_mvqn_recursive_accelerator_constant_forces(self):

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def force1(model_part,guess):
            return self.constant_force(model_part,guess,k1,z_equilibrium_1)

        def force2(model_part,guess):
            return self.constant_force(model_part,guess,k2,z_equilibrium_2)

        def analytical_solution(model_part,k1,k2,z_equilibrium_1,z_equilibrium_2):
            s = KratosMultiphysics.Vector(3*len(model_part.Nodes))
            dtot = z_equilibrium_2 - z_equilibrium_1
            z_solution = dtot * k2 / (k1+k2) + z_equilibrium_1
            for i,node in enumerate(model_part.Nodes):
                j = 3*i
                s[j]   = 0.0
                s[j+1] = 0.0
                s[j+2] = z_solution - node.Z

            return s

        def solution(model_part):
            return analytical_solution(model_part,k1,k2,z_equilibrium_1,z_equilibrium_2)

        self._test_mvqn_recursive_accelerator(force1,force2,solution)

    def test_mvqn_recursive_accelerator_variable_stiffness(self):

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def forceA(model_part,guess):
            return self.constant_force(model_part,guess,k1,z_equilibrium_1)

        def forceB(model_part,guess):
            return self.variable_stiffness(model_part,guess,k2,z_equilibrium_2)

        def analytical_solution(model_part,k1,k2,z_equilibrium_1,z_equilibrium_2):
            s = KratosMultiphysics.Vector(3*len(model_part.Nodes))
            dtot = z_equilibrium_2 - z_equilibrium_1
            for i,node in enumerate(model_part.Nodes):
                k2_variable = k2 * ( 1 + node.X*node.Y )
                z_solution = dtot * k2_variable / (k1+k2_variable) + z_equilibrium_1
                j = 3*i
                s[j]   = 0.0
                s[j+1] = 0.0
                s[j+2] = z_solution - node.Z

            return s

        def solution(model_part):
            return analytical_solution(model_part,k1,k2,z_equilibrium_1,z_equilibrium_2)

        self._test_mvqn_recursive_accelerator(forceA,forceB,solution)

    # Aitken accelerator test
    def _test_aitken_accelerator(self,force1,force2,solution):

        aitken_settings = KratosMultiphysics.Parameters("""{
                                                            "solver_type"        : "Relaxation",
                                                            "acceleration_type"  : "Aitken",
                                                            "w_0"                : 0.825
                                                           }""")

        self._test_accelerator(force1,force2,solution,aitken_settings)

    def test_aitken_accelerator_constant_forces(self):

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def force1(model_part,guess):
            return self.constant_force(model_part,guess,k1,z_equilibrium_1)

        def force2(model_part,guess):
            return self.constant_force(model_part,guess,k2,z_equilibrium_2)

        def analytical_solution(model_part,k1,k2,z_equilibrium_1,z_equilibrium_2):
            s = KratosMultiphysics.Vector(3*len(model_part.Nodes))
            dtot = z_equilibrium_2 - z_equilibrium_1
            z_solution = dtot * k2 / (k1+k2) + z_equilibrium_1
            for i,node in enumerate(model_part.Nodes):
                j = 3*i
                s[j]   = 0.0
                s[j+1] = 0.0
                s[j+2] = z_solution - node.Z

            return s

        def solution(model_part):
            return analytical_solution(model_part,k1,k2,z_equilibrium_1,z_equilibrium_2)

        self._test_aitken_accelerator(force1,force2,solution)

    def test_aitken_accelerator_variable_stiffness(self):

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def forceA(model_part,guess):
            return self.constant_force(model_part,guess,k1,z_equilibrium_1)

        def forceB(model_part,guess):
            return self.variable_stiffness(model_part,guess,k2,z_equilibrium_2)

        def analytical_solution(model_part,k1,k2,z_equilibrium_1,z_equilibrium_2):
            s = KratosMultiphysics.Vector(3*len(model_part.Nodes))
            dtot = z_equilibrium_2 - z_equilibrium_1
            for i,node in enumerate(model_part.Nodes):
                k2_variable = k2 * ( 1 + node.X*node.Y )
                z_solution = dtot * k2_variable / (k1+k2_variable) + z_equilibrium_1
                j = 3*i
                s[j]   = 0.0
                s[j+1] = 0.0
                s[j+2] = z_solution - node.Z

            return s

        def solution(model_part):
            return analytical_solution(model_part,k1,k2,z_equilibrium_1,z_equilibrium_2)

        self._test_aitken_accelerator(forceA,forceB,solution)

if __name__ == '__main__':
    test = ConvergenceAcceleratorSpringTest()
    test.setUp()
    test.test_mvqn_recursive_accelerator_variable_stiffness()
    test.tearDown()