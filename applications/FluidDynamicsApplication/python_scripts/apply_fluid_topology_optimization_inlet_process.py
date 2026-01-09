import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

from KratosMultiphysics.read_csv_table_utility import ReadCsvTableUtility
from KratosMultiphysics import assign_vector_by_direction_process

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyFluidTopologyOptimizationInletProcess(Model, settings["Parameters"])


class ApplyFluidTopologyOptimizationInletProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = self.GetDefaultParameters()

        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        else:
            if (not settings["velocity_modulus"].IsString()):
                raise Exception("Not string_type velocity modulus in fluid topology optimization inlet.")
            elif (settings["velocity_modulus"].GetString() == ""):
                raise Exception("Fluid topology optimization Inlet velocity modulus is empty.")
            if (not settings["adjoint_velocity_modulus"].IsString()):
                raise Exception("Not string_type adjoint velocity modulus in fluid topology optimization inlet.")
            elif (settings["adjoint_velocity_modulus"].GetString() == ""):
                raise Exception("Fluid topology optimization Inlet adjoint velocity modulus is empty.")
            else:
                if (settings.Has("direction")):
                    if (not settings["direction"].IsString()):
                        raise Exception("Not string_type direction in fluid topology optimization inlet.")
                else:
                    raise Exception("Empty direction string in fluid topology optimization inlet. Set a valid model part name.")
        
        self.use_time_dependent_velocity_scale_factor = False
        if settings.Has("time_dependent_velocity_scale_factor"):
            self.use_time_dependent_velocity_scale_factor = True
            if settings["time_dependent_velocity_scale_factor"].IsNumber():
                raise Exception("Imposed numeric_type for time_dependent_velocity_scale_factor in fluid topology optimization inlet.")
            elif settings["time_dependent_velocity_scale_factor"].IsString():
                raise Exception("Imposed string_type for time_dependent_velocity_scale_factor in fluid topology optimization inlet.")

        settings.ValidateAndAssignDefaults(default_settings)

        model_part_name = settings["model_part_name"].GetString()
        velocity_modulus = settings["velocity_modulus"].GetString()
        adjoint_velocity_modulus = settings["adjoint_velocity_modulus"].GetString()
        direction = settings["direction"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()
        if self.use_time_dependent_velocity_scale_factor:
            self.model_part = Model.GetModelPart(settings["model_part_name"].GetString())
            self.time_dependent_velocity_scale_factor_csv_table = ReadCsvTableUtility(settings["time_dependent_velocity_scale_factor"]).Read(self.model_part)

        velocity_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : \"""" + model_part_name + """\",
                "variable_name"   : "VELOCITY",
                "modulus"         : \"""" + velocity_modulus + """\",
                "constrained"     : true,
                "direction"       : \"""" + direction + """\",
                "interval"        : """ + interval + """
            }
            """)

        adjoint_velocity_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : \"""" + model_part_name + """\",
                "variable_name"   : "VELOCITY_ADJ",
                "modulus"         : \"""" + adjoint_velocity_modulus + """\",
                "constrained"     : true,
                "direction"       : \"""" + direction + """\",
                "interval"        : """ + interval + """
            }
            """)
            
        # Set the INLET flag in the inlet model part nodes and conditions
        self.inlet_model_part = Model[settings["model_part_name"].GetString()]
        for node in self.inlet_model_part.Nodes:
            node.Set(KratosMultiphysics.INLET, True)
        for condition in self.inlet_model_part.Conditions:
            condition.Set(KratosMultiphysics.INLET, True)

        # Construct the base process AssignVectorByDirectionProcess
        self.aux_process         = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, velocity_settings)
        self.aux_process_adjoint = assign_vector_by_direction_process.AssignVectorByDirectionProcess(Model, adjoint_velocity_settings)

        if self.use_time_dependent_velocity_scale_factor:
            self.base_velocity_modulus_function_str_list = []
            for process in self.aux_process.aux_processes:
                self.base_velocity_modulus_function_str_list.append(process.aux_function.FunctionBody())
    
    @staticmethod
    def GetDefaultParameters():
        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name" : "",
                "velocity_modulus"         : "0.0",
                "adjoint_velocity_modulus" : "0.0",
                "direction"       : "automatic_inwards_normal",
                "interval"        : [0.0,"End"],
                "time_dependent_velocity_scale_factor" : {}
            }
            """)
        return default_settings

    def SetTimeDependentVelocityScaleFactor(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        current_time_velocity_scale_factor = self.time_dependent_velocity_scale_factor_csv_table.GetValue(current_time)
        count_process = 0
        for process in self.aux_process.aux_processes:
            current_time_velocity_modulus_function_str = "(" + self.base_velocity_modulus_function_str_list[count_process] + ")*(" + str(current_time_velocity_scale_factor) + ")"
            process.value_is_function = True
            process.function_string = current_time_velocity_modulus_function_str
            process.aux_function = KratosMultiphysics.GenericFunctionUtility(process.function_string)
            if process.aux_function.DependsOnSpace():
                process.cpp_apply_function_utility = KratosMultiphysics.ApplyFunctionToNodesUtility(process.mesh.Nodes, process.aux_function )
            count_process += 1

    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        if self.IsPhysicsStage():
            if self.use_time_dependent_velocity_scale_factor:
                self.SetTimeDependentVelocityScaleFactor()
                
            self.aux_process.ExecuteInitializeSolutionStep()
        elif self.IsAdjointStage():
            self.aux_process_adjoint.ExecuteInitializeSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("'ExecuteInitializeSolutionStep' for ApplyFluidTopologyOptimizationInletProcess called during a non valid stage.")

    def ExecuteFinalizeSolutionStep(self):
        # Call the base process ExecuteFinalizeSolutionStep()
        if self.IsPhysicsStage():
            self.aux_process.ExecuteFinalizeSolutionStep()
        elif self.IsAdjointStage():
            self.aux_process_adjoint.ExecuteFinalizeSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("'ExecuteFinalizeSolutionStep' for ApplyFluidTopologyOptimizationInletProcess called during a non valid stage.")

    def IsPhysicsStage(self):
        top_opt_stage = self.inlet_model_part.ProcessInfo.GetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE)
        if top_opt_stage == 1:
            return True
        else:
            return False
        
    def IsAdjointStage(self):
        top_opt_stage = self.inlet_model_part.ProcessInfo.GetValue(KratosCFD.FLUID_TOP_OPT_PROBLEM_STAGE)
        if top_opt_stage == 2:
            return True
        else:
            return False
