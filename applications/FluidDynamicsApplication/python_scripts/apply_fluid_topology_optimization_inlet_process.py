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

        default_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : 0,
                "model_part_name" : "",
                "velocity_modulus"         : "0.0",
                "adjoint_velocity_modulus" : "0.0",
                "direction"       : "automatic_inwards_normal",
                "interval"        : [0.0,"End"],
                "use_time_dependent_velocity_peak_from_csv_file" : false,
                "time_dependent_velocity_peak_from_csv_file_settings" : {}
            }
            """)

        if (settings.Has("direction")):
            if (not settings["direction"].IsString()):
                raise Exception("Not string_type direction in fluid topology optimization inlet.")
        # Check input data
        if (settings["model_part_name"].GetString() == ""):
            raise Exception("Empty inlet model part name string. Set a valid model part name.")
        else:
            if (settings["velocity_modulus"].GetString() == ""):
                raise Exception("Fluid topology optimization Inlet velocity modulus is empty.")
            elif (settings["adjoint_velocity_modulus"].GetString() == ""):
                raise Exception("Fluid topology optimization Inlet adjoint velocity modulus is empty.")
        settings.ValidateAndAssignDefaults(default_settings)

        mesh_id = settings["mesh_id"].GetInt()
        model_part_name = settings["model_part_name"].GetString()
        velocity_modulus = settings["velocity_modulus"].GetString()
        adjoint_velocity_modulus = settings["adjoint_velocity_modulus"].GetString()
        direction = settings["direction"].GetString()
        interval = settings["interval"].PrettyPrintJsonString()

        velocity_settings = KratosMultiphysics.Parameters("""
            {
                "mesh_id"         : """ + str(mesh_id) +  """,
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
                "mesh_id"         : """ + str(mesh_id) + """,
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

        self.HandleCsvTableImport(Model, settings)

    def HandleCsvTableImport(self, Model, settings):
        if settings["use_time_dependent_velocity_peak_from_csv_file"].GetBool():
            self.use_time_dependent_velocity_peak_from_csv_file = True
            model_part = Model[settings["model_part_name"].GetString()]
            self.csv_table_utility = ReadCsvTableUtility(settings["time_dependent_velocity_peak_from_csv_file_settings"]).Read(model_part)
        else:
            self.use_time_dependent_velocity_peak_from_csv_file = False

    def SetTimeDependentVelocityPeakFromCsvTable(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

    def ExecuteInitializeSolutionStep(self):
        # Call the base process ExecuteInitializeSolutionStep()
        if self.IsPhysicsStage():
            if self.use_time_dependent_velocity_peak_from_csv_file:
                self.SetTimeDependentVelocityPeakFromCsvTable()
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
