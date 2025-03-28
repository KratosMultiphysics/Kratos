import math
import KratosMultiphysics
from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFluidTopologyOptimizationOutletProcess(Model, settings["Parameters"])


class ApplyFluidTopologyOptimizationOutletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"            : 0,
            "model_part_name"    : "",
            "constrained"        : false,
            "pressure_value"         : {},
            "adjoint_pressure_value" : {},
            "interval"           : [0.0,"End"],
            "outlet_inflow_contribution" : false,
            "outlet_inflow_contribution_characteristic_velocity_calculation" : "outlet_average"
        }
        """)

        # Trick: allows "value" to be a double, a string or a table value (otherwise the ValidateAndAssignDefaults might fail)
        if(settings.Has("pressure_value")):
            if(settings["pressure_value"].IsString()):
                default_settings["pressure_value"].SetString("0.0")
            elif settings["pressure_value"].IsNumber():
                default_settings["pressure_value"].SetDouble(0.0)
        else:
            err_msg = "Provided settings have no 'pressure_value'. This needs to be provided."
            raise Exception(err_msg)
        if(settings.Has("adjoint_pressure_value")):
            if(settings["adjoint_pressure_value"].IsString()):
                default_settings["adjoint_pressure_value"].SetString("0.0")
            elif settings["adjoint_pressure_value"].IsNumber():
                default_settings["adjoint_pressure_value"].SetDouble(0.0)
        else:
            err_msg = "Provided settings have no 'adjoint_pressure_value'. This needs to be provided."
            raise Exception(err_msg)
        # Check the core processes input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty outlet pressure model part name. Set a valid model part name.")
        
        settings.ValidateAndAssignDefaults(default_settings)

        mesh_id = settings["mesh_id"].GetInt()
        model_part_name = settings["model_part_name"].GetString()
        constrained = settings["constrained"].GetBool()
        if (settings["pressure_value"].IsString()):
            pressure_value = settings["pressure_value"].GetString()
        else:
            pressure_value = str(settings["pressure_value"].GetDouble())
        if (settings["adjoint_pressure_value"].IsString()):
            adjoint_pressure_value = settings["adjoint_pressure_value"].GetString()
        else:
            adjoint_pressure_value = str(settings["adjoint_pressure_value"].GetDouble())
        interval = settings["interval"].PrettyPrintJsonString()
        outlet_inflow_contribution = settings["outlet_inflow_contribution"].GetBool()
        outlet_inflow_contribution_characteristic_velocity_calculation = settings["outlet_inflow_contribution_characteristic_velocity_calculation"].GetString()
        
        physics_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"            : """ + str(mesh_id) + """,
            "model_part_name"    : \"""" + model_part_name + """\",
            "variable_name"      : "PRESSURE",
            "constrained"        : """ + str(constrained).lower() + """,
            "value"              : \"""" + pressure_value + """\",
            "interval"           : """ + interval + """,
            "outlet_inflow_contribution" : """ + str(outlet_inflow_contribution).lower() + """,
            "outlet_inflow_contribution_characteristic_velocity_calculation" : \"""" + outlet_inflow_contribution_characteristic_velocity_calculation + """\"
        }
        """)
        
        adjoint_settings = KratosMultiphysics.Parameters("""
        {
            "mesh_id"            : """ + str(mesh_id) + """,
            "model_part_name"    : \"""" + model_part_name + """\",
            "variable_name"      : "PRESSURE_ADJ",
            "constrained"        : """ + str(constrained).lower() + """,
            "value"              : \"""" + adjoint_pressure_value + """\",
            "interval"           : """ + interval + """,
            "outlet_inflow_contribution" : """ + str(outlet_inflow_contribution).lower() + """,
            "outlet_inflow_contribution_characteristic_velocity_calculation" : \"""" + outlet_inflow_contribution_characteristic_velocity_calculation + """\"
        }
        """)

        pres_settings, ext_pres_settings = self.__DefinePressureSettings(physics_settings)
        adj_pres_settings, adj_ext_pres_settings = self.__DefinePressureSettings(adjoint_settings, pressure_variable="PRESSURE_ADJ", external_pressure_variable="EXTERNAL_PRESSURE_ADJ")

        # Set the OUTLET flag in the outlet model part nodes and conditions
        self.outlet_model_part = Model[model_part_name]
        for node in self.outlet_model_part.Nodes:
            node.Set(KratosMultiphysics.OUTLET, True)
        for condition in self.outlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.OUTLET, True)

        # Outlet inflow contribution
        self.__outlet_inflow_contribution = outlet_inflow_contribution
        if self.__outlet_inflow_contribution:
            self.outlet_model_part.ProcessInfo[KratosCFD.OUTLET_INFLOW_CONTRIBUTION_SWITCH] = True
            self.__outlet_char_velocity = None
            self.__outlet_char_velocity_type = outlet_inflow_contribution_characteristic_velocity_calculation
            if (not (self.__outlet_char_velocity_type in ["outlet_average", "outlet_max"])):
                raise Exception(f"Wrong 'outlet_inflow_contribution_characteristic_velocity_calculation'. Provided value is '{self.__outlet_char_velocity_type }'.")
        else:
            self.outlet_model_part.ProcessInfo[KratosCFD.OUTLET_INFLOW_CONTRIBUTION_SWITCH] = False

        # Construct the physics problem base process AssignValueProcess
        self.aux_pressure_process = AssignScalarVariableProcess(Model, pres_settings)
        self.aux_external_pressure_process = AssignScalarVariableProcess(Model, ext_pres_settings)
        # Construct the adjoint problem base process AssignValueProcess
        self.aux_adjoint_pressure_process = AssignScalarVariableProcess(Model, adj_pres_settings)
        self.aux_adjoint_external_pressure_process = AssignScalarVariableProcess(Model, adj_ext_pres_settings)


    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        self.aux_pressure_process.ExecuteInitializeSolutionStep()
        self.aux_external_pressure_process.ExecuteInitializeSolutionStep()
        self.aux_adjoint_pressure_process.ExecuteInitializeSolutionStep()
        self.aux_adjoint_external_pressure_process.ExecuteInitializeSolutionStep()

        # Set the characteristic outlet velocity in the outlet model part ProcessInfo
        if self.__outlet_inflow_contribution:
            # If required, update the outlet characteristic velocity
            if (self.outlet_model_part.ProcessInfo[KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE] == 1):
                velocity_variable=KratosMultiphysics.VELOCITY
            elif (self.outlet_model_part.ProcessInfo[KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE] == 2):
                velocity_variable=KratosMultiphysics.VELOCITY_ADJ
            else:
                raise Exception("Evaluating Outlet Inflow Contribution for FLuid Topology Optimization Outlet Process in the wrong 'FLUID_TOP_OPT_PROBLEM_STAGE'.")
            if self.__outlet_char_velocity_type == "outlet_average":
                self.__outlet_char_velocity = self.__ComputeOutletAverageVelocity(velocity_variable)
            else:
                self.__outlet_char_velocity = self.__ComputeOutletMaxVelocity(velocity_variable)
            # Save value to be used in the boundary terms integration
            self.outlet_model_part.ProcessInfo[KratosCFD.CHARACTERISTIC_VELOCITY] = self.__outlet_char_velocity


    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        self.aux_pressure_process.ExecuteFinalizeSolutionStep()
        self.aux_external_pressure_process.ExecuteFinalizeSolutionStep()


    def __ComputeOutletAverageVelocity(self, velocity_variable):
        # Compute the outlet average velocity
        outlet_avg_vel_norm = 0.0
        for node in self.outlet_model_part.GetCommunicator().LocalMesh().Nodes:
            vnode = node.GetSolutionStepValue(velocity_variable)
            outlet_avg_vel_norm += math.sqrt(vnode[0]*vnode[0] + vnode[1]*vnode[1] + vnode[2]*vnode[2])
        comm = self.outlet_model_part.GetCommunicator().GetDataCommunicator()
        outlet_avg_vel_norm = comm.SumAll(outlet_avg_vel_norm)
        tot_len = len(self.outlet_model_part.GetCommunicator().LocalMesh().Nodes)   # Partition outlet model part number of nodes
        tot_len = comm.SumAll(tot_len)                                              # Get the total outlet model part nodes
        outlet_avg_vel_norm /= tot_len
        # Check the outlet average velocity minimum value and return
        min_outlet_avg_vel_norm = 1.0e-12
        return outlet_avg_vel_norm if outlet_avg_vel_norm >= min_outlet_avg_vel_norm else min_outlet_avg_vel_norm


    def __ComputeOutletMaxVelocity(self, velocity_variable):
        # Compute the outlet max velocity
        outlet_max_vel_norm = 0.0
        for node in self.outlet_model_part.GetCommunicator().LocalMesh().Nodes:
            v_node = node.GetSolutionStepValue(velocity_variable)
            aux_val = v_node[0]*v_node[0] + v_node[1]*v_node[1] + v_node[2]*v_node[2]
            if aux_val > outlet_max_vel_norm:
                outlet_max_vel_norm = aux_val
        outlet_max_vel_norm = self.outlet_model_part.GetCommunicator().GetDataCommunicator().MaxAll(outlet_max_vel_norm)
        outlet_max_vel_norm = math.sqrt(outlet_max_vel_norm)
        # Check the outlet max velocity minimum value and return
        min_outlet_vel_norm = 1.0e-12
        return outlet_max_vel_norm if outlet_max_vel_norm >= min_outlet_vel_norm else min_outlet_vel_norm


    def __DefinePressureSettings(self, pressure_settings, pressure_variable="PRESSURE", external_pressure_variable="EXTERNAL_PRESSURE"):
        # Set a Kratos parameters suitable for the core processes to set the PRESSURE
        pres_settings = pressure_settings.Clone()
        pres_settings.RemoveValue("hydrostatic_outlet")
        pres_settings.RemoveValue("h_top")
        pres_settings.RemoveValue("outlet_inflow_contribution")
        pres_settings.RemoveValue("outlet_inflow_contribution_characteristic_velocity_value")
        pres_settings.RemoveValue("outlet_inflow_contribution_characteristic_velocity_calculation")
        # Create a copy of the PRESSURE settings to set the EXTERNAL_PRESSURE
        ext_pres_settings = pres_settings.Clone()
        ext_pres_settings["constrained"].SetBool(False)
        ext_pres_settings["variable_name"].SetString(external_pressure_variable)
        # Check the core processes input data
        if (pres_settings["model_part_name"].GetString() == ""):
            raise Exception("Empty outlet pressure model part name. Set a valid model part name.")
        elif (ext_pres_settings["model_part_name"].GetString() == ""):
            raise Exception("Empty outlet external pressure model part name. Set a valid model part name.")
        elif (pres_settings["variable_name"].GetString() != pressure_variable):
            raise Exception("Outlet pressure settings variable_name is not " + pressure_variable + ".")
        elif (ext_pres_settings["variable_name"].GetString() != external_pressure_variable):
            raise Exception("Outlet external pressure settings variable_name is not " + external_pressure_variable + ".")
        elif (pres_settings["value"].IsString()):
            if (pres_settings["value"].GetString == ""):
                raise Exception("Outlet pressure function sting is empty.")
        elif (ext_pres_settings["value"].IsString()):
            if (ext_pres_settings["value"].GetString == ""):
                raise Exception("Outlet external pressure function sting is emp+ty.")  
        return pres_settings, ext_pres_settings