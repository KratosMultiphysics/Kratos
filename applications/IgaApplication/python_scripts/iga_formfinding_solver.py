from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.IgaApplication as IgaApplication

# Import base class file
import structural_mechanics_solver


def CreateSolver(main_model_part, custom_settings):
    return FormfindingMechanicalSolver(main_model_part, custom_settings)


class FormfindingMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    """The structural mechanics formfinding solver.

    This class creates the mechanical solver for formdinding.

    Public member variables:
    formfinding_settings -- settings for the formfinding solver.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, main_model_part, custom_settings):
        # Set defaults and validate custom settings.
        self.formfinding_settings = KratosMultiphysics.Parameters("""
        {
            "print_formfinding_iterations": false
        }
        """)
        self.validate_and_transfer_matching_settings(custom_settings, self.formfinding_settings)
        # Validate the remaining settings in the base class.

        # Construct the base solver.
        super(FormfindingMechanicalSolver, self).__init__(main_model_part, custom_settings)
        self.print_on_rank_zero("::[FormfindingMechanicalSolver]:: ", "Construction finished")

    def ImportModelPart(self):
        """This function imports the ModelPart"""
        self._set_nurbs_brep_modeler()

        with open(self.settings["model_import_settings"]["physics_filename"].GetString(),'r') as physics_file:
                physics_parameters = KratosMultiphysics.Parameters( physics_file.read())

        physics_file = open(self.settings["model_import_settings"]["physics_filename"].GetString(),'r')
        physics_parameters = KratosMultiphysics.Parameters( physics_file.read())
        self.nurbs_brep_modeler.ImportModelPart(self.main_model_part, physics_parameters)

    def _set_nurbs_brep_modeler(self):
        """Prepare the nurbs brep modeler and read in the necessary data. """
        # This function prepares the nurbs brep modeler and reads in the rough geometry data,
        # which can be used for surface descriptions and integration domains.

        self.nurbs_brep_modeler = IgaApplication.NurbsBrepModeler(self.main_model_part)

        if self.settings["model_import_settings"]["input_type"].GetString() == "json":
            with open(self.settings["model_import_settings"]["input_filename"].GetString() + ".json",'r') as geometry_file:
                geometry_parameters = KratosMultiphysics.Parameters( geometry_file.read())
            self.geometry_reader = IgaApplication.BrepJsonIO()

            self.nurbs_brep_modeler.ImportGeometry(self.geometry_reader, geometry_parameters)

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_mechanical_solution_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return StructuralMechanicsApplication.FormfindingUpdatedReferenceStrategy(
                                                                computing_model_part,
                                                                mechanical_scheme,
                                                                linear_solver,
                                                                mechanical_convergence_criterion,
                                                                builder_and_solver,
                                                                self.settings["max_iteration"].GetInt(),
                                                                self.settings["compute_reactions"].GetBool(),
                                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                self.settings["move_mesh_flag"].GetBool(),
                                                                self.formfinding_settings["print_formfinding_iterations"].GetBool(),
                                                                self.settings["line_search"].GetBool())
