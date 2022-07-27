import math
import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

import KratosMultiphysics.FluidDynamicsApplication as KratosFluid


def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")
    return Apply2FOutletProcess(Model, settings["Parameters"])


class Apply2FOutletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"    : ""
        }
        """)

        settings.ValidateAndAssignDefaults(default_settings)

        # Set a Kratos parameters suitable for the core processes to set the PRESSURE
        pres_settings = settings.Clone()

        # Check the core processes input data
        if (pres_settings["model_part_name"].GetString() == ""):
            raise Exception(
                "Empty outlet pressure model part name. Set a valid model part name.")

        # Set the OUTLET flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[pres_settings["model_part_name"].GetString(
        )]

    # def ExecuteInitialize(self):

    #     model = self.outlet_model_part.GetModel()
    #     comp_mp = model.GetModelPart("FluidModelPart.fluid_computational_model_part")
    #     master_node = self.outlet_model_part.GetParentModelPart().Nodes[2657]
    #     mpc_id = 1
    #     for node in self.outlet_model_part.Nodes:
    #         print(node.Id)
    #         comp_mp.CreateNewMasterSlaveConstraint(
    #             "LinearMasterSlaveConstraint", mpc_id, master_node, KratosMultiphysics.VELOCITY_X, node, KratosMultiphysics.VELOCITY_X, 1.0, 0.0)
    #         mpc_id += 1


    def ExecuteInitializeSolutionStep(self):
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)

        for element in self.outlet_model_part.GetParentModelPart().Elements:
            for node in element.GetNodes():
                if node.Is(KratosMultiphysics.OUTLET):
                    element.Set(KratosMultiphysics.OUTLET, True)
                    break

        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)




    # def ExecuteInitializeSolutionStep(self):
    #     # Call the base process ExecuteInitializeSolutionStep()
    #     self.aux_pressure_process.ExecuteInitializeSolutionStep()
    #     self.aux_external_pressure_process.ExecuteInitializeSolutionStep()

    #     # If considered, add the hydrostatic component to the outlet pressure
    #     if (self.hydrostatic_outlet):
    #         self._AddOutletHydrostaticComponent()

    #     # Compute the outlet average velocity
    #     self._ComputeOutletCharacteristicVelocity()


    # def ExecuteFinalizeSolutionStep(self):
    #     # Call the base process ExecuteFinalizeSolutionStep()
    #     self.aux_pressure_process.ExecuteFinalizeSolutionStep()
    #     self.aux_external_pressure_process.ExecuteFinalizeSolutionStep()



