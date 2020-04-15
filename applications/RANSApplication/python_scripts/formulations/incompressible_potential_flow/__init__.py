from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import formulations
from .incompressible_potential_flow_velocity_formulation import IncompressiblePotentialFlowVelocityFormulation
from .incompressible_potential_flow_pressure_formulation import IncompressiblePotentialFlowPressureFormulation


class IncompressiblePotentialFlowFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(IncompressiblePotentialFlowFormulation, self).__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "incompressible_potential_flow",
            "velocity_potential_flow_solver_settings":{},
            "pressure_potential_flow_solver_settings":{}
        }''')
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.AddFormulation(IncompressiblePotentialFlowVelocityFormulation(model_part, settings["velocity_potential_flow_solver_settings"]))
        self.AddFormulation(IncompressiblePotentialFlowPressureFormulation(model_part, settings["pressure_potential_flow_solver_settings"]))

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.VELOCITY_POTENTIAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.PRESSURE_POTENTIAL)

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.VELOCITY_POTENTIAL, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.PRESSURE_POTENTIAL, self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step dofs.")

    def GetMinimumBufferSize(self):
        return 1

    def IsPeriodic(self):
        return False

    def SetIsPeriodic(self, value):
        super(IncompressiblePotentialFlowFormulation, self).SetIsPeriodic(False)
        if (value):
            raise Exception("Periodic conditions are not supported by incompressible potential flow solver.")