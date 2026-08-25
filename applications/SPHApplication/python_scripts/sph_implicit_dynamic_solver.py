import numpy as np

# Importing the Kratos Library
import KratosMultiphysics
# Import applications
import KratosMultiphysics.SPHApplication as SPH
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
# Import base class file
from KratosMultiphysics.SPHApplication.sph_solver import SPHSolver

def CreateSolver(model, custom_settings):
    return ImplicitSPHSolver(model, custom_settings)

class ImplicitSPHSolver(SPHSolver):
    
    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitSPHSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "time_integration_method" : "implicit",
            "scheme_type"             : "newmark",
            "damp_factor_m"           : 0.0,
            "newmark_beta"            : 0.25,
            "rayleigh_alpha"          : 0.0,
            "rayleigh_beta"           : 0.0
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults
    
    def AddVariables(self):
        super().AddVariables()
        self._add_dynamic_variables()
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitSPHSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        super().AddDofs()
        self._add_dynamic_dofs()
        KratosMultiphysics.Logger.PrintInfo("::[ImplicitSPHSolver]:: ", "DOF's ADDED")

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()

        # Some pre-processes may affect the system of equations, we rebuild the equation ids
        process_info = self.main_model_part.ProcessInfo
        if process_info[KratosMultiphysics.STEP] == 1 and process_info[StructuralMechanicsApplication.RESET_EQUATION_IDS]:
            # Resetting the global equations ids
            self._GetBuilderAndSolver().SetUpSystem(self.GetComputingModelPart())
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.ExposeSystemMatrix()

    def _CreateScheme(self):
        scheme_type = self.settings["scheme_type"].GetString()

        # Setting the Rayleigh damping parameters
        process_info = self.main_model_part.ProcessInfo
        process_info[StructuralMechanicsApplication.RAYLEIGH_ALPHA] = self.settings["rayleigh_alpha"].GetDouble()
        process_info[StructuralMechanicsApplication.RAYLEIGH_BETA] = self.settings["rayleigh_beta"].GetDouble()

        # Setting the time integration schemes
        if scheme_type in ("newmark", "bossak"):
            scheme_settings = KratosMultiphysics.Parameters("""{
                "damp_factor_m" : 0.0,
                "newmark_beta" : 0.0,
                "projection_variables_list" : []
            }""")
            scheme_settings["damp_factor_m"].SetDouble(0.0 if scheme_type == "newmark" else self.settings["damp_factor_m"].GetDouble())
            scheme_settings["newmark_beta"].SetDouble(self.settings["newmark_beta"].GetDouble())
            sph_scheme = StructuralMechanicsApplication.StructuralMechanicsBossakScheme(scheme_settings)
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"newmark\", \"bossak\", \"pseudo_static\", \"backward_euler\", \"bdf1\", \"bdf2\", \"bdf3\", \"bdf4\", \"bdf5\", \"relaxation\""
            raise Exception(err_msg)
        return sph_scheme
    
    ## DEBUG METHODS
    
    def ExposeSystemMatrix(self):
        A = self._GetSolutionStrategy().GetSystemMatrix() 
        vec = self._GetSolutionStrategy().GetSystemVector() 
        sol = self._GetSolutionStrategy().GetSolutionVector()

        n = A.Size1()
        m = A.Size2()
        B = np.zeros((n,m))
        vec1 = np.zeros(n)
        sol1 = np.zeros(n)

        for i in range(n):
            vec1[i] = vec[i]
            sol1[i] = sol[i]
            for j in range(m):
                B[i ,j] = A[i, j]
        
        self.system_matrix = A
        self.right_hand_side = vec
        self.solution = sol

