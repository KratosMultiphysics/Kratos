from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.FSIApplication as KratosFSI

try:
    from KratosMultiphysics.mpi import *
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos
    import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
except ImportError:
    pass

import KratosMultiphysics.KratosUnittest as KratosUnittest

import math
import convergence_accelerator_factory

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/AcceleratorSpringTests/" + fileName

def GetPartitionedFilePath(fileName):
    return GetFilePath( "{0}_{1}".format(fileName,mpi.rank ) )

class ConvergenceAcceleratorSpringMPITest(KratosUnittest.TestCase):

    def constant_force(self,model_part,variable,k,reference_z):
        for i,node in enumerate(model_part.Nodes):
            j = 3*i
            position = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0) + node.Z
            if node.Z > reference_z:
                delta = position - reference_z
            else:
                delta = reference_z - position
            forces = KratosMultiphysics.Array3()
            forces[0] = 0.0
            forces[1] = 0.0
            forces[2] = k*delta
            node.SetSolutionStepValue(variable,0,forces)

    def variable_stiffness(self,model_part,variable,k,reference_z):
        for i,node in enumerate(model_part.Nodes):
            j = 3*i
            position = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0) + node.Z
            if node.Z > reference_z:
                delta = position - reference_z
            else:
                delta = reference_z - position
            forces = KratosMultiphysics.Array3()
            forces[0] = 0.0
            forces[1] = 0.0
            forces[2] = k * ( 1 + node.X*node.Y ) * delta
            node.SetSolutionStepValue(variable,0,forces)


    def ComputeResidual(self,model_part,guess,force1,force2):

        force1(model_part,KratosMultiphysics.FORCE)
        force2(model_part,KratosMultiphysics.REACTION)

        residual = self.partitioned_utilities.SetUpInterfaceVector(model_part)
        self.partitioned_utilities.ComputeInterfaceResidualVector(model_part,KratosMultiphysics.FORCE,
                                                                  KratosMultiphysics.REACTION,
                                                                  residual)
        return residual

    def ComputeResidualNorm(self,residual):
        return self.space.TwoNorm(residual.GetReference())

    def ReadModelPart(self,filename):

        model_part = KratosMultiphysics.ModelPart("Test ModelPart")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FSI_INTERFACE_RESIDUAL)
        mpi.world.barrier()

        KratosMetis.SetMPICommunicatorProcess(model_part).Execute()

        model_part_io = KratosMultiphysics.ModelPartIO( filename )
        model_part_io.ReadModelPart(model_part)
        model_part.SetBufferSize(2)

        KratosTrilinos.ParallelFillCommunicator(model_part).Execute()

        return model_part

    def InitializeNodalEpetraVector(self,model_part,rows_per_node):
        v = self.space.CreateEmptyVectorPointer(self.epetra_comm)
        local_size = rows_per_node * len(model_part.GetCommunicator().LocalMesh().Nodes)
        global_size = model_part.GetCommunicator().SumAll(local_size)
        self.space.ResizeVector(v,int(global_size))
        return v

    def setUp(self):

        # So far, the MPI convergence accelerator tests must be run with 2 processes
        if (mpi.size != 2):
            raise Exception("The MPI convergence accelerator tests must be run with 2 processes.")

        self.print_gid_output = False
        self.accelerator_tolerance = 1e-10
        self.accelerator_iterations = 50
        self.assert_delta = 1e-7

        self.model_part = self.ReadModelPart(GetPartitionedFilePath("box_fluid"))
        self.model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)

        self.space = KratosTrilinos.TrilinosSparseSpace()
        self.epetra_comm = KratosTrilinos.CreateCommunicator()

        self.partitioned_utilities = KratosTrilinos.TrilinosPartitionedFSIUtilities3D(self.epetra_comm)

    def tearDown(self):
        if self.print_gid_output:
            local_nodes = self.model_part.GetCommunicator().LocalMesh().Nodes

            gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostBinary
            multifile = KratosMultiphysics.MultiFileFlag.SingleFile
            deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
            write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteConditions

            gid_io = KratosMultiphysics.GidIO(GetPartitionedFilePath("box_fluid"),
                                              gid_mode, multifile,
                                              deformed_mesh_flag, write_conditions)

            gid_io.InitializeMesh(0)
            gid_io.WriteMesh(self.model_part.GetMesh())
            gid_io.FinalizeMesh()

            gid_io.InitializeResults(0.0,self.model_part.GetMesh())
            gid_io.WriteNodalResults(KratosMultiphysics.PARTITION_INDEX,self.model_part.Nodes,0.0,0)
            gid_io.WriteNodalResults(KratosMultiphysics.DISPLACEMENT,self.model_part.Nodes,0.0,0)
            gid_io.WriteNodalResults(KratosMultiphysics.FSI_INTERFACE_RESIDUAL,local_nodes,0.0,0)
            gid_io.WriteNodalResults(KratosMultiphysics.FORCE,local_nodes,0.0,0)
            gid_io.WriteNodalResults(KratosMultiphysics.REACTION,local_nodes,0.0,0)
            gid_io.FinalizeResults()

        # clean temporary files
        import os,glob
        if mpi.rank == 0:
            for f in glob.glob(GetFilePath('*.time')):
                os.remove(f)

    def test_accelerator(self,accelerator_settings,force1,force2,solution):

        print("")
        print("Testing accelerator: ",accelerator_settings["solver_type"].GetString())

        top_part = self.model_part.GetSubModelPart("Top")

        # Construct the accelerator strategy
        coupling_utility = convergence_accelerator_factory.CreateTrilinosConvergenceAccelerator(top_part, self.epetra_comm, accelerator_settings)

        coupling_utility.Initialize()

        nl_it = 0
        convergence = False

        coupling_utility.InitializeSolutionStep()

        x_guess = self.partitioned_utilities.SetUpInterfaceVector(top_part)
        residual = self.ComputeResidual(top_part,x_guess,force1,force2)
        res_norm = self.ComputeResidualNorm(residual)

        while (nl_it <= self.accelerator_iterations):

            print(mpi.rank,": Iteration: ", nl_it," residual norm: ", res_norm, file=sys.stderr)

            if res_norm > self.accelerator_tolerance:
                coupling_utility.InitializeNonLinearIteration()
                coupling_utility.UpdateSolution(residual, x_guess)
                self.partitioned_utilities.UpdateInterfaceValues(top_part, KratosMultiphysics.DISPLACEMENT, x_guess)
                coupling_utility.FinalizeNonLinearIteration()
            else:
                coupling_utility.FinalizeSolutionStep()
                convergence = True
                break

            nl_it += 1
            residual = self.ComputeResidual(top_part,x_guess,force1,force2)
            res_norm = self.ComputeResidualNorm(residual)

        # Check the obtained solution
        expected_x = solution(top_part)

        for i,node in enumerate(top_part.Nodes):
            expected = expected_x[3*i+2]
            obtained = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0)
            self.assertAlmostEqual(expected,obtained,delta=self.assert_delta)

    # Aitken accelerator tests
    def test_aitken_accelerator(self,force1,force2,solution):

        aitken_settings = KratosMultiphysics.Parameters("""{
                                                            "solver_type"        : "Relaxation",
                                                            "acceleration_type"  : "Aitken",
                                                            "w_0"                : 0.825
                                                           }""")

        self.test_accelerator(aitken_settings, force1, force2, solution)

    def test_aitken_accelerator_constant_forces(self):
        self.print_gid_output = False

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def force1(model_part,variable):
            return self.constant_force(model_part,variable,k1,z_equilibrium_1)

        def force2(model_part,variable):
            return self.constant_force(model_part,variable,k2,z_equilibrium_2)

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

        self.test_aitken_accelerator(force1,force2,solution)

    def test_aitken_accelerator_variable_stiffness(self):

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def forceA(model_part,variable):
            return self.constant_force(model_part,variable,k1,z_equilibrium_1)

        def forceB(model_part,variable):
            return self.variable_stiffness(model_part,variable,k2,z_equilibrium_2)

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

        self.test_aitken_accelerator(forceA,forceB,solution)

    def test_aitken_accelerator_ghost_nodes(self):
        self.print_gid_output = False

        # relax tolerance requirements to force differences between processors
        self.accelerator_tolerance = 1e-2
        self.accelerator_iterations = 10
        self.assert_delta = 1e-3

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def forceA(model_part,variable):
            return self.constant_force(model_part,variable,k1,z_equilibrium_1)

        def forceB(model_part,variable):
            return self.variable_stiffness(model_part,variable,k2,z_equilibrium_2)

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

        self.test_aitken_accelerator(forceA,forceB,solution)

        ghost_displacements = [ node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0) for node in self.model_part.GetCommunicator().GhostMesh(0).Nodes ]
        self.model_part.GetCommunicator().SynchronizeNodalSolutionStepsData()
        owner_displacements = [ node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0) for node in self.model_part.GetCommunicator().GhostMesh(0).Nodes ]

        for vold,vnew in zip(ghost_displacements,owner_displacements):
            self.assertAlmostEqual(vold,vnew,delta=1e-7)

    # MVQN recursive accelerator tests
    def test_mvqn_recursive_accelerator(self,force1,force2,solution):

        mvqn_recursive_settings = KratosMultiphysics.Parameters("""{
                                                                    "solver_type"        : "MVQN_recursive",
                                                                    "w_0"                : 0.825,
                                                                    "buffer_size"        : 5
                                                                }""")

        self.test_accelerator(mvqn_recursive_settings, force1, force2, solution)

    def test_mvqn_recursive_accelerator_constant_forces(self):
        self.print_gid_output = False

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def force1(model_part,variable):
            return self.constant_force(model_part,variable,k1,z_equilibrium_1)

        def force2(model_part,variable):
            return self.constant_force(model_part,variable,k2,z_equilibrium_2)

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

        self.test_mvqn_recursive_accelerator(force1,force2,solution)

    def test_mvqn_recursive_accelerator_variable_stiffness(self):

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def forceA(model_part,variable):
            return self.constant_force(model_part,variable,k1,z_equilibrium_1)

        def forceB(model_part,variable):
            return self.variable_stiffness(model_part,variable,k2,z_equilibrium_2)

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

        self.test_mvqn_recursive_accelerator(forceA,forceB,solution)

    def test_mvqn_recursive_accelerator_ghost_nodes(self):
        self.print_gid_output = False

        # relax tolerance requirements to force differences between processors
        self.accelerator_tolerance = 1e-2
        self.accelerator_iterations = 15
        self.assert_delta = 1e-3

        k1 = 100
        k2 = 500
        z_equilibrium_1 = 0.5
        z_equilibrium_2 = 1.5

        def forceA(model_part,variable):
            return self.constant_force(model_part,variable,k1,z_equilibrium_1)

        def forceB(model_part,variable):
            return self.variable_stiffness(model_part,variable,k2,z_equilibrium_2)

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

        self.test_mvqn_recursive_accelerator(forceA,forceB,solution)

        ghost_displacements = [ node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0) for node in self.model_part.GetCommunicator().GhostMesh(0).Nodes ]
        self.model_part.GetCommunicator().SynchronizeNodalSolutionStepsData()
        owner_displacements = [ node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z,0) for node in self.model_part.GetCommunicator().GhostMesh(0).Nodes ]

        for vold,vnew in zip(ghost_displacements,owner_displacements):
            self.assertAlmostEqual(vold,vnew,delta=1e-7)
