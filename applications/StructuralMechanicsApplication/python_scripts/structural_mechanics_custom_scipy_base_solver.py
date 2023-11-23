# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Import scipy modules
import KratosMultiphysics.scipy_conversion_tools
from scipy.sparse.linalg import eigsh

def CreateSolver(main_model_part, custom_settings):
    return CustomScipyBaseSolver(main_model_part, custom_settings)

class CustomScipyBaseSolver(MechanicalSolver):
    """The structural mechanics custom scipy base solver.

    This class creates the mechanical solvers to provide mass and stiffness matrices as scipy matrices.

    Derived class must override the function SolveSolutionStep. In there the Mass and Stiffness matrices
    can be obtained as scipy matrices.
    The computation of the eigenvalue problem in this implementation is only an example how this solver is to be used.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super().__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[CustomScipyBaseSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "scheme_type"         : "dynamic"
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    #### Private functions ####
    def _CreateScheme(self):
        """Create the scheme for the scipy solver.

        The scheme determines the mass and stiffness matrices
        """
        scheme_type = self.settings["scheme_type"].GetString()
        if scheme_type == "dynamic":
            solution_scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else: # here e.g. a stability scheme could be added
            err_msg =  "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"dynamic\""
            raise Exception(err_msg)

        return solution_scheme

    def _CreateLinearSolver(self):
        ''' Linear solver will not be used. But eventually the solution strategy calls the solver's clear function.
        To avoid crashing linear solver is provided here'''
        return KratosMultiphysics.LinearSolver()

    def _CreateSolutionStrategy(self):
        if self.settings["builder_and_solver_settings"]["use_block_builder"].GetBool():
            warn_msg = "In case an eigenvalue problem is computed an elimantion builder shall be used to ensure boundary conditions are applied correctly!"
            KratosMultiphysics.Logger.PrintWarning("CustomScipyBaseSolver", warn_msg)

        eigen_scheme = self._GetScheme() # The scheme defines the matrices
        computing_model_part = self.GetComputingModelPart()
        builder_and_solver = self._GetBuilderAndSolver()

        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              eigen_scheme,
                                                              builder_and_solver,
                                                              False,
                                                              False,
                                                              False,
                                                              False )

    def _MassMatrixComputation(self):
        space = KratosMultiphysics.UblasSparseSpace()
        self.GetComputingModelPart().ProcessInfo.SetValue(StructuralMechanicsApplication.BUILD_LEVEL,1) #Mass Matrix
        scheme = self._GetSolutionStrategy().GetScheme()

        aux = self._GetSolutionStrategy().GetSystemMatrix()
        space.SetToZeroMatrix(aux)

        # Create dummy vectors
        b = space.CreateEmptyVectorPointer()
        space.ResizeVector( b, space.Size1(aux) )
        space.SetToZeroVector(b)

        xD = space.CreateEmptyVectorPointer()
        space.ResizeVector( xD, space.Size1(aux) )
        space.SetToZeroVector(xD)

        # Build matrix
        builder_and_solver = self._GetBuilderAndSolver()
        builder_and_solver.Build(scheme, self.GetComputingModelPart(), aux, b)
        # Apply Constraints
        builder_and_solver.ApplyConstraints(scheme, self.GetComputingModelPart(), aux, b)
        # Apply Boundary Conditions
        builder_and_solver.ApplyDirichletConditions(scheme, self.GetComputingModelPart(), aux, xD, b)
        # Convert Mass matrix to scipy
        M = KratosMultiphysics.scipy_conversion_tools.to_csr(aux)

        return M

    def _StiffnessMatrixComputation(self):
        space = KratosMultiphysics.UblasSparseSpace()
        self.GetComputingModelPart().ProcessInfo.SetValue(StructuralMechanicsApplication.BUILD_LEVEL,2) #Stiffness Matrix
        scheme = self._GetSolutionStrategy().GetScheme()

        aux = self._GetSolutionStrategy().GetSystemMatrix()
        space.SetToZeroMatrix(aux)

        # Create dummy vectors
        b = space.CreateEmptyVectorPointer()
        space.ResizeVector( b, space.Size1(aux) )
        space.SetToZeroVector(b)

        xD = space.CreateEmptyVectorPointer()
        space.ResizeVector( xD, space.Size1(aux) )
        space.SetToZeroVector(xD)

        # Build matrix
        builder_and_solver = self._GetBuilderAndSolver()
        builder_and_solver.Build(scheme, self.GetComputingModelPart(), aux, b)
        # Apply constraints
        builder_and_solver.ApplyConstraints(scheme, self.GetComputingModelPart(), aux, b)
        # Apply boundary conditions
        builder_and_solver.ApplyDirichletConditions(scheme, self.GetComputingModelPart(), aux, xD, b)
        # Convert stiffness matrix to scipy
        K = KratosMultiphysics.scipy_conversion_tools.to_csr(aux)

        return K

    def _AssignVariables(self, eigenvalues, eigenvectors):
        num_eigenvalues = eigenvalues.size
        # Store eigenvalues in process info
        eigenvalue_vector = self.GetComputingModelPart().ProcessInfo.GetValue(StructuralMechanicsApplication.EIGENVALUE_VECTOR)
        eigenvalue_vector.Resize(num_eigenvalues)
        for i in range(num_eigenvalues):
            eigenvalue_vector[i] = eigenvalues[i]
        self.GetComputingModelPart().ProcessInfo.SetValue(StructuralMechanicsApplication.EIGENVALUE_VECTOR, eigenvalue_vector)

        # Store eigenvectors in nodes
        for node in self.GetComputingModelPart().Nodes:
            node_eigenvectors = node.GetValue(StructuralMechanicsApplication.EIGENVECTOR_MATRIX)
            # if self.settings["rotation_dofs"].GetBool() == True:
            #     dofs = [node.GetDof(KratosMultiphysics.ROTATION_X),
            #             node.GetDof(KratosMultiphysics.ROTATION_Y),
            #             node.GetDof(KratosMultiphysics.ROTATION_Z),
            #             node.GetDof(KratosMultiphysics.DISPLACEMENT_X),
            #             node.GetDof(KratosMultiphysics.DISPLACEMENT_Y),
            #             node.GetDof(KratosMultiphysics.DISPLACEMENT_Z)]

            #     node_eigenvectors.Resize(num_eigenvalues, 6 )
            # else:
            #     dofs = [node.GetDof(KratosMultiphysics.DISPLACEMENT_X),
            #             node.GetDof(KratosMultiphysics.DISPLACEMENT_Y),
            #             node.GetDof(KratosMultiphysics.DISPLACEMENT_Z)]
            #     node_eigenvectors.Resize(num_eigenvalues, 3 )

            # # Fill the eigenvector matrix
            # for i in range(num_eigenvalues):
            #     j = -1
            #     for dof in dofs:
            #         j = j + 1
            #         if dof.IsFixed():
            #             node_eigenvectors[i,j] = 0.0
            #         else:
            #             node_eigenvectors[i,j] = eigenvectors[dof.EquationId,i]
            # node.SetValue(StructuralMechanicsApplication.EIGENVECTOR_MATRIX, node_eigenvectors)

    def SolveSolutionStep(self):
        """This method must be overriden in derived class.

        The computation of the egenvalue problem is only an example how this solver is to be used.
        """
        ## Obtain scipy matrices
        M = self._MassMatrixComputation()
        K = self._StiffnessMatrixComputation()

        #mp = self.GetComputingModelPart()
        #mp.Clear()
        ## Compute eigenvalues and eigenvectors
        tolerance = 1e-3
        iteration = M.size*1000

        tol = 1e-10
        import scipy
        import numpy as np
        import os
        import sys, slepc4py
        slepc4py.init(sys.argv)
        from petsc4py import PETSc
        from slepc4py import SLEPc

        opts = PETSc.Options()

        B = PETSc.Mat().createAIJ(size=M.shape,
                                        csr=(M.indptr, M.indices,
                                        M.data),comm=PETSc.COMM_WORLD)
        B.assemble()
        A = PETSc.Mat().createAIJ(size=K.shape,
                                        csr=(K.indptr, K.indices,
                                        K.data),comm=PETSc.COMM_WORLD)
        A.assemble()
        kk = 1
        E = SLEPc.EPS().create(PETSc.COMM_WORLD)
        E.setOperators(A,B)
        E.setDimensions(kk,PETSc.DECIDE)     # set number of  eigenvalues to compute
        st = E.getST()
        st.setType(slepc4py.SLEPc.ST.Type.SINVERT)

        KK = st.getKSP()
        KK.setType(PETSc.KSP.Type.PREONLY)
        PC = KK.getPC()
        PC.setType(PETSc.PC.Type.LU)
        #PC.setFactorSolverType("mumps")
        PC.setFactorShift(PETSc.Mat.FactorShiftType.NONZERO)
        #PC.setFromOptions()

        E.setType(SLEPc.EPS.Type.KRYLOVSCHUR)
        E.setWhichEigenpairs(E.Which.TARGET_REAL)
        E.setProblemType(SLEPc.EPS.ProblemType.GNHEP)    # generalized non-Hermitian eigenvalue problem
        E.setTarget(10000)
        E.setTolerances(tol)   # set tolerance
        PC.setFromOptions()

        KK.setPC(PC)
        KK.setFromOptions()

        # if sigma is not None:
        #     E.setTarget(sigma)     # set the desired eigenvalue
        # E.setTarget(sigma)
        # if v0 is not None:
        #     E.setInitialSpace(v0)  # set the initial vector

        # st = E.getST()
        #st.setType('sinvert')

        st.setKSP(KK)
        st.setFromOptions()

        E.setST(st)
        E.setFromOptions()
        E.solve()

        xr, tmp = A.getVecs()
        xi, tmp = A.getVecs()

        print("Convergence: ", E.getConverged() )
        print("Iterations: ", E.getIterationNumber() )

        vals = []
        vecs = []
        for i in range(E.getConverged()):
            val = E.getEigenpair(i, xr, xi)
            #vecr, veci = E.getEigenvector(i, xr, xi)
            vals.append(val)
            #vecs.append(complex(vecr, veci))

        print("Scepc : ", vals)
        # vals = asarray(vals)
        # vecs = asarray(vecs).T

                                       # K_ = K.todense()
        # M_ = M.todense()

        iteration=min(M.size*1000, np.iinfo(np.int32).max)
        tolerance = 1e-6
        iteration = 1000
        # print("Size;;;;;;;;;;;;;;;;;;;;;;;;;; ", M.shape[0] )
        # print("eigenvalue Start: ")
        # #n_dofs = M.shape[0]
        # os.environ["OMP_NUM_THREADS"] = "1"
        # vals, vecs = scipy.sparse.linalg.eigsh(K, k=2, M=M, which='LM',  maxiter=iteration, tol=tolerance, return_eigenvectors=True)
        # vals, vecs = scipy.sparse.linalg.eigs(K, 3, M, which='LM', tol=tolerance, maxiter = iteration)
        # #vals, vecs = scipy.linalg.eigsh(K_, M_, subset_by_index=[n_dofs-2, n_dofs-1], turbo=False)

        # # vals, vecs = scipy.linalg.eig(K_, M_)
        #print("Scipy: ", vals)
        # vals = np.sort(vals)
        # new_vals = []
        # for i, val in enumerate(vals):
        #     if i > 6 and i < 20:
        #         new_vals.append(val)
        # print("eigenvalue End: ")
        ## Assign results to Kratos variables
        #self._AssignVariables(np.array(new_vals),vecs)
        self._AssignVariables(np.array(vals),np.array(vecs))
        return True #converged
