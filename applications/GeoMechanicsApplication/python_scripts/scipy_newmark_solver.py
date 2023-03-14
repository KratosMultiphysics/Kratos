import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication
import KratosMultiphysics.scipy_conversion_tools
from KratosMultiphysics.gid_output_process import GiDOutputProcess
from KratosMultiphysics.GeoMechanicsApplication.geomechanics_U_Pw_solver import UPwSolver

import numpy as np
import scipy
import scipy.sparse as sparse
import scipy.sparse.linalg as la


def CreateSolver(model, custom_settings):
    return ScipyNewmarkSolver(model, custom_settings)


class ScipyNewmarkSolver(UPwSolver):

    def __init__(self, model, custom_settings):
        super(ScipyNewmarkSolver, self).__init__(model, custom_settings)
        self.beta = 0.25
        self.gamma = 0.5

        self.dt = 0

        self.max_iter = 100

    def Assemble(self):
        KratosMultiphysics.Timer.Start("Kratos Assemble")

        # construct graph of system of equations
        self.dof_list = self.builder_and_solver.GetDofSet()
        self.Agraph = KratosMultiphysics.SparseContiguousRowGraph(len(self.dof_list))

        for elem in self.main_model_part.Elements:
            elem.Initialize(self.GetComputingModelPart().ProcessInfo)
            equation_ids = elem.EquationIdVector(self.GetComputingModelPart().ProcessInfo)
            self.Agraph.AddEntries(equation_ids)

        for cond in self.main_model_part.Conditions:
            equation_ids = cond.EquationIdVector(self.GetComputingModelPart().ProcessInfo)
            self.Agraph.AddEntries(equation_ids)

        self.Agraph.Finalize()

        # Build FEM matrix
        rhs = KratosMultiphysics.Vector()
        lhs = KratosMultiphysics.Matrix()
        M = KratosMultiphysics.Matrix()
        C = KratosMultiphysics.Matrix()

        first_derivative = KratosMultiphysics.Vector()
        second_derivative = KratosMultiphysics.Vector()

        self.velocity_vec=np.zeros(len(self.dof_list))
        self.acceleration_vec = np.zeros(len(self.dof_list))

        self.A = KratosMultiphysics.CsrMatrix(self.Agraph)
        self.M = KratosMultiphysics.CsrMatrix(self.Agraph)
        self.C = KratosMultiphysics.CsrMatrix(self.Agraph)
        self.b = KratosMultiphysics.SystemVector(self.Agraph)

        for i in range(self.b.size()):
            self.b[i] = 0

        self.A.BeginAssemble()
        self.M.BeginAssemble()
        self.C.BeginAssemble()
        self.b.BeginAssemble()

        for elem in self.main_model_part.Elements:
            equation_ids = elem.EquationIdVector(self.GetComputingModelPart().ProcessInfo)
            elem.CalculateLocalSystem(lhs, rhs, self.GetComputingModelPart().ProcessInfo)
            elem.CalculateMassMatrix(M, self.GetComputingModelPart().ProcessInfo)
            elem.CalculateDampingMatrix(C, self.GetComputingModelPart().ProcessInfo)


            elem.GetFirstDerivativesVector(first_derivative)
            elem.GetSecondDerivativesVector(second_derivative)

            self.velocity_vec[equation_ids] = first_derivative
            self.acceleration_vec[equation_ids] = first_derivative

            self.A.Assemble(lhs, equation_ids)
            self.M.Assemble(M, equation_ids)
            self.C.Assemble(C, equation_ids)

            self.b.Assemble(rhs, equation_ids)

        for cond in self.main_model_part.Conditions:
            equation_ids = cond.EquationIdVector(self.GetComputingModelPart().ProcessInfo)
            cond.CalculateLocalSystem(lhs, rhs, self.GetComputingModelPart().ProcessInfo)
            cond.CalculateMassMatrix(M, self.GetComputingModelPart().ProcessInfo)
            cond.CalculateDampingMatrix(C, self.GetComputingModelPart().ProcessInfo)

            self.A.Assemble(lhs, equation_ids)

            if M.Size1() > 0:
                self.M.Assemble(M, equation_ids)
            if C.Size1() > 0:
                self.C.Assemble(C, equation_ids)

            self.b.Assemble(rhs, equation_ids)

        self.A.FinalizeAssemble()
        self.M.FinalizeAssemble()
        self.C.FinalizeAssemble()
        self.b.FinalizeAssemble()

        self.ApplyDirichlet()

        self.A_scipy = KratosMultiphysics.scipy_conversion_tools.to_csr(self.A)
        self.M_scipy = KratosMultiphysics.scipy_conversion_tools.to_csr(self.M)
        self.C_scipy = KratosMultiphysics.scipy_conversion_tools.to_csr(self.C)

        # self.ApplyDirichlet()

        # self.get_dof_values()

        self.add_dynamics_to_lhs()
        self.add_dynamics_to_rhs()

        self.calculate_inverse_lhs()

        KratosMultiphysics.Timer.Stop("Kratos Assemble")

        return self.A, self.M, self.C

    def build_rhs(self):

        self.b = KratosMultiphysics.SystemVector(self.Agraph)

        for i in range(self.b.size()):
            self.b[i] = 0

        first_derivative = KratosMultiphysics.Vector()
        second_derivative = KratosMultiphysics.Vector()
        #self.acceleration_vec

        rhs = KratosMultiphysics.Vector()
        for elem in self.main_model_part.Elements:
            elem.CalculateRightHandSide(rhs, self.GetComputingModelPart().ProcessInfo)
            equation_ids = elem.EquationIdVector(self.GetComputingModelPart().ProcessInfo)

            elem.GetFirstDerivativesVector(first_derivative)
            elem.GetSecondDerivativesVector(second_derivative)

            self.velocity_vec[equation_ids] = first_derivative
            self.acceleration_vec[equation_ids] = second_derivative

            self.b.Assemble(rhs, equation_ids)

        for cond in self.main_model_part.Conditions:
            cond.CalculateRightHandSide(rhs, self.GetComputingModelPart().ProcessInfo)
            equation_ids = cond.EquationIdVector(self.GetComputingModelPart().ProcessInfo)

            self.b.Assemble(rhs, equation_ids)

        # self.builder_and_solver.ApplyDirichletConditions(self.scheme, self.main_model_part)
        self.ApplyDirichlet()
        self.add_dynamics_to_rhs()

    # def get_dof_values(self):
    #
    #     solution_list = []
    #     first_derivative_list = []
    #     second_derivative_list = []
    #     pressure_list = []
    #
    #     self.solution_vec = np.array(self.dof_list.GetValues())
    #
    #     # self.displacement_vec = []
    #     self.first_derivative_vec = []
    #     self.second_derivative_vec = []
    #
    #     ndim = 2
    #     for node in self.main_model_part.Nodes:
    #         disp = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)
    #         vel = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
    #         acc = node.GetSolutionStepValue(KratosMultiphysics.ACCELERATION)
    #         pressure = node.GetSolutionStepValue(KratosMultiphysics.WATER_PRESSURE)
    #
    #         # pressure_list.append(pressure)
    #
    #         solution_list.append(disp[0])
    #         first_derivative_list.append(vel[0])
    #         second_derivative_list.append(acc[0])
    #         solution_list.append(disp[1])
    #         first_derivative_list.append(vel[1])
    #         second_derivative_list.append(acc[1])
    #         if ndim == 3:
    #             solution_list.append(disp[2])
    #             first_derivative_list.append(vel[2])
    #             second_derivative_list.append(acc[2])
    #
    #         solution_list.append(pressure)
    #         first_derivative_list.append(0)
    #         second_derivative_list.append(0)
    #
    #     # self.solution_vec = np.asarray(solution_list)
    #     self.velocity_vec = np.flip(np.asarray(first_derivative_list))
    #     self.acceleration_vec = np.flip(np.asarray(second_derivative_list))

    def add_dynamics_to_rhs(self):



        # self.get_dof_values()

        mass_contribution = self.M_scipy.dot(self.acceleration_vec)

        damp_contribution = self.C_scipy.dot(self.velocity_vec)

        self.b = self.b - mass_contribution - damp_contribution

    def add_dynamics_to_lhs(self):
        self.dt = self.ComputeDeltaTime()

        mass_contribution = 1 / (self.beta * self.dt ** 2) * self.M_scipy
        damping_contribution = self.gamma / (self.beta * self.dt) * self.C_scipy

        self.lhs = self.A_scipy + mass_contribution + damping_contribution

    def calculate_inverse_lhs(self):
        identity = sparse.identity(self.lhs.shape[0])

        self.inverse_lhs = la.splu(self.lhs, identity)

    def ApplyDirichlet(self):
        KratosMultiphysics.Timer.Start("Apply Dirchlet")
        free_dofs_vector = KratosMultiphysics.Vector(self.A.size1())
        # self.dof_list = self.builder_and_solver.GetDofSet()

        for dof in self.dof_list:
            if (dof.IsFixed()):
                free_dofs_vector[dof.EquationId] = 0.0
            else:
                free_dofs_vector[dof.EquationId] = 1.0
        self.A.ApplyHomogeneousDirichlet(free_dofs_vector, 1.0, self.b)
        self.M.ApplyHomogeneousDirichlet(free_dofs_vector, 1.0, self.b)
        self.C.ApplyHomogeneousDirichlet(free_dofs_vector, 1.0, self.b)

        KratosMultiphysics.Timer.Stop("Apply Dirichlet")

    def SolveSolutionStep(self):

        self.dt = self.ComputeDeltaTime()

        if self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] == 1:
            self.Assemble()

        coo_lhs = self.lhs.tocoo()

        K = KratosMultiphysics.CompressedMatrix(len(self.dof_list), len(self.dof_list))

        for data, row, col in zip(coo_lhs.data, coo_lhs.row, coo_lhs.col):
            K[row, col] = data

        self.build_rhs()
        self.SolveSystem()

        # force_ext_prev = np.copy(self.b)

        # Newton Raphson loop where force is updated in every iteration
        is_converged = False
        i = 0
        while not is_converged and i < self.max_iter:
            # self.main_model_part.GetProcessInfo()["NL_ITERATION_NUMBER"].SetInt(i)
            self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.NL_ITERATION_NUMBER] = i
            self.scheme.InitializeNonLinIteration(self.main_model_part, K, self.dx, self.b)
            self.convergence_criterion.InitializeNonLinearIteration(self.main_model_part, self.dof_list,
                                                                    K, self.dx, self.b)


            #todo incorporate this if required
            # is_converged = self.convergence_criterion.PreCriteria(self.main_model_part, self.dof_list,
            #                                                         K, self.dx, self.b)

            self.build_rhs()
            self.SolveSystem()

            self.scheme.Update(self.main_model_part, self.dof_list, K, self.dx, self.b)
            self.scheme.FinalizeNonLinIteration(self.main_model_part, K, self.dx, self.b)
            self.convergence_criterion.FinalizeNonLinearIteration(self.main_model_part, self.dof_list,
                                                                  K, self.dx, self.b)


            is_converged = self.convergence_criterion.PostCriteria(self.main_model_part, self.dof_list,
                                                               K, self.dx, self.b)

            i+=1

        # self.solver.CalculateReactions(self.scheme, self.main_model_part, K, self.dx, self.b)

        return is_converged

    def SolveSystem(self):
        KratosMultiphysics.Timer.Start("system solve")
        # solve system matrix

        self.dx = self.inverse_lhs.solve(self.b)

        # self.dx = la.spsolve(self.lhs, self.b)  # note that in Kratos we always solve for the increment
        KratosMultiphysics.Timer.Stop("system solve")

