from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_implicit_dynamic_solver


def CreateSolver(model, custom_settings):
    return IgaImplicitDynamicSolver(model, custom_settings)


class IgaImplicitDynamicSolver(structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver):
    """The iga implicit dynamic solver derived from the ImplicitMechanicalSolver.

    This class creates the mechanical solvers for implicit dynamic analysis.
    It currently supports Newmark, Bossak and dynamic relaxation schemes.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_implicit_dynamic_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        super(IgaImplicitDynamicSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[IgaImplicitDynamicSolver]:: ", "Construction finished")

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