from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
from KratosMultiphysics import ExternalSolversApplication
from KratosMultiphysics import SolidMechanicsApplication
from KratosMultiphysics import StructuralMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

import solid_mechanics_solver

def CreateSolver(main_model_part, custom_settings):
    return EigenSolver(main_model_part, custom_settings)

class EigenSolver(solid_mechanics_solver.MechanicalSolver):

    def __init__(self, main_model_part, custom_settings): 
        
        self.main_model_part = main_model_part
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "structural_mechanics_eigensolver",
            "echo_level": 0,
            "buffer_size": 1,
            "solution_type": "Dynamic",
            "analysis_type": "Linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "eigensolver_settings":{
                "solver_type": "FEAST"
            },
            "problem_domain_sub_model_part_list": ["solid_model_part"],
            "processes_sub_model_part_list": [""]
        }
        """)
        
        ##overwrite the default settings with user-provided parameters 
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        # eigensolver_settings are validated/assigned in the linear_solver
        print("Construction of Eigensolver finished")

    def Initialize(self):
        self.compute_model_part = self.GetComputeModelPart()

        self.eigensolver_settings = self.settings["eigensolver_settings"] 
        if self.eigensolver_settings["solver_type"].GetString() == "FEAST":
            self.linear_solver = ExternalSolversApplication.FEASTSolver(self.eigensolver_settings)
        else:
            raise Exception("solver_type is not yet implemented.")

        if self.settings["solution_type"].GetString() == "Dynamic":
            self.scheme = StructuralMechanicsApplication.EigensolverDynamicScheme()
        else:
            raise Exception("solution_type is not yet implemented.")

        self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = StructuralMechanicsApplication.EigensolverStrategy(
            self.compute_model_part,
            self.scheme,
            self.builder_and_solver)

    def AddVariables(self):
        
        solid_mechanics_solver.MechanicalSolver.AddVariables(self)
   
        print("::[Structural EigenSolver]:: Variables ADDED")

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
            
        if self.settings["rotation_dofs"].GetBool():
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X)
                node.AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y)
                node.AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z)
                
        if self.settings["pressure_dofs"].GetBool():                
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.PRESSURE, KratosSolid.PRESSURE_REACTION)

        print("::[Structural EigenSolver]:: DOF's ADDED")

    def GetComputeModelPart(self):
        
        return self.main_model_part
    
    def Solve(self):
            
        self.solver.Solve()

    def SetEchoLevel(self, level):

        self.solver.SetEchoLevel(level)

