from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.IgaApplication as IgaApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import base class file
import structural_mechanics_static_solver

def CreateSolver(model, custom_settings):
    return IgaStaticSolver(model, custom_settings)


class IgaStaticSolver(structural_mechanics_static_solver.StaticMechanicalSolver):
    """The iga static solver.

    This class is derived from the StaticMechanicalSolver but enhanced with 
    the import modifications that Iga needs.
    """

    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings in the base class.
        # Construct the base solver.
        super(IgaStaticSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[IgaStaticSolver]:: ", "Construction finished")

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
