# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(main_model_part, custom_settings):
    return FormfindingMechanicalSolver(main_model_part, custom_settings)

class FormfindingMechanicalSolver(MechanicalSolver):
    """The structural mechanics formfinding solver.

    This class creates the mechanical solver for formfinding.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Construct the base solver.
        super().__init__(main_model_part, custom_settings)
        custom_settings["projection_settings"].ValidateAndAssignDefaults(self.GetDefaultParameters()["projection_settings"])

        KratosMultiphysics.Logger.PrintInfo("::[FormfindingMechanicalSolver]:: ", "Construction finished")


    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "printing_format"             : "all",
            "write_formfound_geometry_file"    : true,
            "formfinding_model_part_name" : "",
            "projection_settings": {
                "model_part_name"  : "Structure",
                "echo_level"       : 0,
                "projection_type"  : "planar",
                "global_direction" : [1,0,0],
                "variable_name"    : "PLEASE_SPECIFY",
                "visualize_in_vtk" : false,
                "method_specific_settings" : { },
                "check_local_space_dimension" : false
            }
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults



    def Finalize(self):
        super().Finalize()
        if (self.settings["write_formfound_geometry_file"].GetBool()):
            StructuralMechanicsApplication.FormfindingStrategy.WriteFormFoundMdpa(self.GetComputingModelPart())

    def _CreateScheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _CreateSolutionStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetScheme()
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()


        # in some cases not all elements need to be reset by the formfinding strategy
        formfinding_model_part = self.GetComputingModelPart()
        if len(self.settings["formfinding_model_part_name"].GetString())>0:
            formfinding_model_part = computing_model_part.GetSubModelPart(self.settings["formfinding_model_part_name"].GetString())
        return StructuralMechanicsApplication.FormfindingStrategy(
                                                                computing_model_part,
                                                                mechanical_scheme,
                                                                mechanical_convergence_criterion,
                                                                builder_and_solver,
                                                                formfinding_model_part,
                                                                self.settings["write_formfound_geometry_file"].GetBool(),
                                                                self.settings["printing_format"].GetString(),
                                                                self.settings["projection_settings"],
                                                                self.settings["max_iteration"].GetInt(),
                                                                self.settings["compute_reactions"].GetBool(),
                                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                self.settings["move_mesh_flag"].GetBool())
