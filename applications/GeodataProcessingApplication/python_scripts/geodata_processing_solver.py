from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.GeodataProcessingApplication as GeodataProcessingApplication
# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver
# Other imports
from KratosMultiphysics import auxiliary_solver_utilities
import KratosMultiphysics.kratos_utilities as kratos_utils

# Other imports
from importlib import import_module


class GeoDataProcessingSolver(PythonSolver):
    """The base class geo data processing solvers.

    IMPORTANT : This solver is not going to run any simulation. This
    will only output a .mdpa file and a JSON file which are to be input
    for the FluidDynamicsApplication.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _create_solution_scheme which
    constructs and returns a solution scheme. Depending on the type of
    solver, derived classes may also need to override the following functions:

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, model, custom_settings):
        super(GeoDataProcessingSolver, self).__init__(model, custom_settings)

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "terrain_file_input_parameters" : {
                "file_type":"stl",
                "file_name":"PLEASE_SPECIFY"
            },
            "buildings_file_input_parameters" : {
                "file_type":"obj",
                "file_name":"PLEASE_SPECIFY"
            },
            "extrusion_height" : 0.0,
            "domain_diameter"  : 0.0,
            "num_terrain_refinement_steps" : 1,
            "num_building_refinement_steps" : 1,
            "output_modelpart_name" : "mechanical_solver",
            "inlet_sub_modelpart_name" : "inlet",
            "outlet_sub_modelpart_name" : "outlet",
            "terrain_sub_modelpart_name" : "ground",
            "top_sub_modelpart_name" : "top",
            "domain_size" : 3
        }""")
        this_defaults.AddMissingParameters(super(MechanicalSolver, cls).GetDefaultSettings())
        return this_defaults