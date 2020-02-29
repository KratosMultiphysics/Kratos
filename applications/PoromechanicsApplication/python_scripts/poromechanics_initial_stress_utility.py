from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.PoromechanicsApplication as KratosPoro

from importlib import import_module

class InitialStressUtility(object):

    def __init__(self, model, parameters):

        self.model = model
        self.parameters = parameters

        self.mode = self.parameters["problem_data"]["initial_stress_utility_settings"]["mode"].GetString()

        domain_size = self.parameters["solver_settings"]["domain_size"].GetInt()
        if domain_size == 2:
            self.initial_stress_utility = KratosPoro.InitialStress2DUtilities()
        else:
            self.initial_stress_utility = KratosPoro.InitialStress3DUtilities()

        self.current_model_part = self.model.GetModelPart(self.parameters["solver_settings"]["model_part_name"].GetString())

    def Load(self):

        # Read the initial model part with initial stresses and perform a mapping to transfer them to the current model part
        if self.mode == 'load':

            initial_model_part_name = 'InitialPorousModelPart'

            # Create initial solver (and initial_model_part)
            initial_solver_settings = self.parameters["solver_settings"]
            initial_solver_settings["model_part_name"].SetString(initial_model_part_name)
            initial_solver_settings["model_import_settings"]["input_filename"].SetString(self.parameters["problem_data"]["initial_stress_utility_settings"]["initial_input_filename"].GetString())
            python_module_name = "KratosMultiphysics.PoromechanicsApplication"
            full_module_name = python_module_name + "." + initial_solver_settings["solver_type"].GetString()
            solver_module = import_module(full_module_name)
            initial_solver = solver_module.CreateSolver(self.model, initial_solver_settings)

            initial_solver.AddVariables()
            initial_solver.ImportModelPart()
            # initial_solver.PrepareModelPart()
            initial_solver.AddDofs()

            # Mapping between initial and current model parts
            initial_model_part = self.model.GetModelPart(initial_model_part_name)

            self.initial_stress_utility.TransferInitialStresses(initial_model_part,self.current_model_part)

            # Delete initial_model_part
            self.model.DeleteModelPart(initial_model_part_name)

    def Save(self):

        # Write an mdpa file containing the initial stress tensor as NodalData
        if self.mode == 'save':

            initial_stress_parameters = KratosMultiphysics.Parameters("{}")
            initial_stress_parameters.AddValue("initial_input_filename",self.parameters["problem_data"]["initial_stress_utility_settings"]["initial_input_filename"])

            self.current_model_part.ProcessInfo.SetValue(KratosPoro.NODAL_SMOOTHING, True)

            self.initial_stress_utility.SaveInitialStresses(initial_stress_parameters,self.current_model_part)
