# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(main_model_part, custom_settings):
    return CuttingPatternMechanicalSolver(main_model_part, custom_settings)

class CuttingPatternMechanicalSolver(MechanicalSolver):
    """The structural mechanics cutting pattern solver.

    This class creates the mechanical solver for cutting pattern generation.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super().__init__(main_model_part, custom_settings)
        self.settings["initial_flattening_settings"].ValidateAndAssignDefaults(self.GetDefaultParameters()["initial_flattening_settings"])

        KratosMultiphysics.Logger.PrintInfo("::[CuttingPatternMechanicalSolver]:: ", "Construction finished")


    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "printing_format"             : "all",
            "write_cutting_pattern_geometry_file"    : true,
            "cuttingPattern_model_part_name" : "",
            "use_relaxation" : false,
            "initial_flattening_settings": {
                "projection_type"  : "planar_mean_normal",
                "global_direction" : [0.0, 0.0, 1.0],
                "echo_level"       : 0
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults



    def Finalize(self):
        super().Finalize()
        if (self.settings["write_cutting_pattern_geometry_file"].GetBool()):
            StructuralMechanicsApplication.CuttingPatternStrategy.WriteCuttingPatternMdpa(self.GetComputingModelPart())

    def SetUseRelaxation(self, use_relaxation):
        """Switches cutting-pattern elements between the least-squares stress-flattening
        system (default) and the standard elastic equilibrium (relaxation) system."""
        self._GetSolutionStrategy().SetUseRelaxation(use_relaxation)

    def _CreateScheme(self):
        return StructuralMechanicsApplication.CuttingPatternScheme()

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetScheme()
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()


        # in some cases not all elements need to be reset by the cutting pattern strategy
        cuttingpattern_model_part = self.GetComputingModelPart()
        if len(self.settings["cuttingPattern_model_part_name"].GetString())>0:
            cuttingpattern_model_part = computing_model_part.GetSubModelPart(self.settings["cuttingPattern_model_part_name"].GetString())
        return StructuralMechanicsApplication.CuttingPatternStrategy(
                                                                computing_model_part,
                                                                mechanical_scheme,
                                                                mechanical_convergence_criterion,
                                                                builder_and_solver,
                                                                cuttingpattern_model_part,
                                                                self.settings["write_cutting_pattern_geometry_file"].GetBool(),
                                                                self.settings["printing_format"].GetString(),
                                                                self.settings["initial_flattening_settings"],
                                                                self.settings["max_iteration"].GetInt(),
                                                                self.settings["compute_reactions"].GetBool(),
                                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                self.settings["move_mesh_flag"].GetBool(),
                                                                self.settings["use_relaxation"].GetBool())