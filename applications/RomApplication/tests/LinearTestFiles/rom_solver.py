from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver
import KratosMultiphysics.StructuralMechanicsApplication as SMA

import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory

# Other imports
import os

class ROMSolver(PythonSolver):

    def __init__(self, model, settings):
        model.CreateModelPart("ThermalModelPart")
        super(ROMSolver, self).__init__(model,settings)

    def AddVariables(self):
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KratosMultiphysics.AMBIENT_TEMPERATURE)
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX) 
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KratosMultiphysics.EMISSIVITY) 
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KratosMultiphysics.CONVECTION_COEFFICIENT)   
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.GetComputingModelPart().AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.TEMPERATURE, self.GetComputingModelPart())

    def ImportModelPart(self):
        self._ImportModelPart(self.GetComputingModelPart(), self.settings["model_import_settings"])

    def GetMinimumBufferSize(self):
        1

    def Initialize(self):   
        ##########################################################################
        convection_diffusion_settings =KratosMultiphysics.ConvectionDiffusionSettings()
        convection_diffusion_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)
        convection_diffusion_settings.SetTransferCoefficientVariable(KratosMultiphysics.CONVECTION_COEFFICIENT)
        convection_diffusion_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)    
        convection_diffusion_settings.SetVolumeSourceVariable(KratosMultiphysics.HEAT_FLUX)
        convection_diffusion_settings.SetSurfaceSourceVariable(KratosMultiphysics.FACE_HEAT_FLUX)   
        self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)           
        ##########################################################################    
        
        
        computing_model_part = self.GetComputingModelPart()
        convergence_criterion = self.get_convergence_criterion()
        convection_diffusion_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        rom_parameters=self.settings["rom_settings"]
        builder_and_solver = KratosMultiphysics.ROMBuilderAndSolver(linear_solver, rom_parameters)
        self.strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                              convection_diffusion_scheme,
                                                              linear_solver,
                                                              convergence_criterion,
                                                              builder_and_solver,
                                                              self.settings["max_iteration"].GetInt(),
                                                              self.settings["compute_reactions"].GetBool(),
                                                              False,
                                                              False)
        self.strategy.Initialize()
       
       
    def InitializeSolutionStep(self):
        self.strategy.InitializeSolutionStep()

    def Predict(self):
        self.strategy.Predict()
        print("Predict done")

        

    def SolveSolutionStep(self):
        is_converged = self.strategy.SolveSolutionStep()
        print("SolveSolutionStep done")
        return is_converged

    def FinalizeSolutionStep(self):
        print("Finalize Solution Step done")
        self.strategy.FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        print("Advance in time done")
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.GetComputingModelPart().ProcessInfo[KratosMultiphysics.STEP] += 1
        self.GetComputingModelPart().CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        print("Compute delta t done")
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def GetComputingModelPart(self):
        print("GetComputingModelPart done")
        return self.model.GetModelPart(self.settings["model_part_name"].GetString())
        
        
    #Residual Criterion
    #def _create_convergence_criterion(self):
    #    import base_convergence_criteria_factory as convergence_criteria_factory
    #    convergence_criterion = convergence_criteria_factory.ConvergenceCriteriaFactory(self.settings)
    #    return convergence_criterion.convergence_criterion

    #Displacement Criterion
    def _create_convergence_criterion(self):
        import convergence_criteria_factory 
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self.settings)
        return convergence_criterion.mechanical_convergence_criterion

        
    def get_convergence_criterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._create_convergence_criterion()
        return self._convergence_criterion