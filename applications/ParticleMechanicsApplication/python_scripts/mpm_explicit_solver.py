from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from KratosMultiphysics.ParticleMechanicsApplication.mpm_solver import MPMSolver

def CreateSolver(model, custom_settings):
    return MPMExplicitSolver(model, custom_settings)

class MPMExplicitSolver(MPMSolver):

    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(MPMExplicitSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[MPMExplicitSolver]:: ", "Construction is finished.")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "time_integration_method"   : "explicit",
            "scheme_type"   : "central_difference",
            "stress_update" : "usf",
            "is_fix_explicit_mp_on_grid_edge" : false
        }""")
        this_defaults.AddMissingParameters(super(MPMExplicitSolver, cls).GetDefaultSettings())
        return this_defaults


    def AddVariables(self):
        super(MPMExplicitSolver, self).AddVariables()
        self._AddDynamicVariables(self.grid_model_part)
        grid_model_part = self.GetGridModelPart()

        # Adding explicit variables
        grid_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
        grid_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.RESIDUAL_VECTOR)

        KratosMultiphysics.Logger.PrintInfo("::[MPMExplicitSolver]:: ", "Variables are all added.")

    ### Protected functions ###

    def _CreateSolutionScheme(self):
        grid_model_part = self.GetGridModelPart()
        domain_size = self._GetDomainSize()
        block_size  = domain_size
        if (self.settings["pressure_dofs"].GetBool()):
            block_size += 1

        # Check whether compressibility is considered
        is_compressible = self.settings["compressible"].GetBool()
        grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_COMPRESSIBLE, is_compressible)

        # Check if we are fixing MPs that lie directly on the edge of grid elements
        is_fix_explicit_mp_on_grid_edge = self.settings["is_fix_explicit_mp_on_grid_edge"].GetBool()
        grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_FIX_EXPLICIT_MP_ON_GRID_EDGE, is_fix_explicit_mp_on_grid_edge)

        # Setting the time integration schemes
        scheme_type = self.settings["scheme_type"].GetString()

        if scheme_type == "forward_euler":
            stress_update_option = 10
            stress_update = self.settings["stress_update"].GetString() #0 = USF, 1 = USL, 2 = MUSL
            if stress_update == "usf":
                stress_update_option = 0
            elif stress_update == "usl":
                stress_update_option = 1
            elif stress_update == "musl":
                stress_update_option = 2
            else:
                err_msg = "The requested stress update \"" + stress_update + "\" is not available!\n"
                err_msg += "Available options are: \"usf\", \"usl\",\"musl\""
            grid_model_part.ProcessInfo.SetValue(KratosParticle.EXPLICIT_STRESS_UPDATE_OPTION, stress_update_option)
            grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT_CENTRAL_DIFFERENCE, False)
        elif scheme_type == "central_difference":
            grid_model_part.ProcessInfo.SetValue(KratosParticle.EXPLICIT_STRESS_UPDATE_OPTION, 0)
            grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT_CENTRAL_DIFFERENCE, True)
        else:
            err_msg = "The requested scheme type \"" + scheme_type + "\" is not available!\n"
            err_msg += "Available options are: \"forward_euler\", \"central_difference\""
            raise Exception(err_msg)

        return KratosParticle.MPMExplicitScheme( grid_model_part)

    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
                grid_model_part = self.GetGridModelPart();
                grid_model_part.ProcessInfo.SetValue(KratosParticle.IS_EXPLICIT, True)
                solution_strategy = self._CreateLinearStrategy()
        else:
            err_msg =  "The requested explicit analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available explicit options are: \"linear\""
            raise Exception(err_msg)
        return solution_strategy


    def _CreateLinearStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        solution_scheme = self._GetSolutionScheme()
        reform_dofs_at_each_step = False ## hard-coded, but can be changed upon implementation
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()
        move_mesh_flag = False ## hard-coded
        return KratosParticle.MPMExplicitStrategy(computing_model_part,
                                                      solution_scheme,
                                                      self.settings["compute_reactions"].GetBool(),
                                                      reform_dofs_at_each_step,
                                                      move_mesh_flag)

    def _IsDynamic(self):
        return True