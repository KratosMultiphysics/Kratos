import KratosMultiphysics
from KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMPMParticleDirichletConditionProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplyMPMParticleDirichletConditionProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
            {
                "model_part_name"           : "PLEASE_SPECIFY_MODEL_PART_NAME",
                "material_points_per_condition"   : 0,
                "imposition_type"           : "penalty",
                "penalty_factor"            : 0,
                "variable_name"             : "DISPLACEMENT",
                "constrained"               : "fixed",
                "value"                     : [0.0, "0*t", 0.0],
                "interval"                  : [0.0, 1e30],
                "option"                    : "",
                "is_equal_distributed"      : false,
                "local_axes"                : {}
            }  """ )

         # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        context_string = type(self).__name__
        old_name = 'particles_per_condition'
        new_name = 'material_points_per_condition'
        if DeprecationManager.HasDeprecatedVariable(context_string, settings, old_name, new_name):
            DeprecationManager.ReplaceDeprecatedVariableName(settings, old_name, new_name)

        settings.ValidateAndAssignDefaults(default_parameters)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.model_part_name = settings["model_part_name"].GetString()
        self.model = Model
        self.material_points_per_condition = settings["material_points_per_condition"].GetInt()
        self.imposition_type = settings["imposition_type"].GetString()
        self.is_neumann_boundary = False
        self.option = settings["option"].GetString()

        #is_equal_distributed = false (material point conditions at Gauss Point Positions)
        #is_equal_distributed = true (material point conditions equally distributed; also at nodes; only 2D)
        self.is_equal_distributed = settings["is_equal_distributed"].GetBool()

        """
        Set boundary_condition_type:
        1. penalty
        2. lagrange (WIP)
        3. fixdof (WIP)
        """

        # set type of boundary
        if (self.imposition_type == "penalty" or self.imposition_type == "Penalty"):
            self.penalty_factor = settings["penalty_factor"].GetDouble()
            self.boundary_condition_type = 1
        else:
            err_msg =  "The requested type of Dirichlet boundary imposition: \"" + self.imposition_type + "\" is not available!\n"
            err_msg += "Available option is: \"penalty\"."
            raise Exception(err_msg)

        # check constraint
        self.constrained = settings["constrained"].GetString()
        self.is_slip_boundary = False
        self.is_contact_boundary = False
        if (self.constrained == "fixed"):
            pass
        elif (self.constrained == "contact"):
            self.is_contact_boundary = True
        elif (self.constrained == "slip"):
            self.is_slip_boundary = True
        elif (self.constrained == "contact_slip"):
            self.is_contact_boundary = True
            self.is_slip_boundary = True
        else:
            err_msg =  "The requested type of constrain: \"" + self.constrained + "\" is not available!\n"
            err_msg += "Available options are: \"fixed\", \"contact\" and \"slip\"."
            raise Exception(err_msg)

        # get variable imposed and check
        variable_name = settings["variable_name"].GetString()
        variable_name_list = ["DISPLACEMENT","VELOCITY","ACCELERATION"]
        if(variable_name in variable_name_list):
            self.variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
        else:
            err_msg =  "The given variable \"" + variable_name + "\" is not available to be imposed with this process.\n"
            err_msg += "Available options are: " + ", ".join(variable_name_list)
            raise Exception(err_msg)

        self.value_is_numeric = [False, False, False]
        self.value = KratosMultiphysics.Vector(3)
        self.aux_function = ["0.0","0.0","0.0"]
        self.name = ["0.0","0.0","0.0"]
        # Loop over components X, Y and Z
        for i, variable in enumerate(["_X", "_Y", "_Z"]):
            self.name[i] = settings["variable_name"].GetString() + variable
            if settings["value"][i].IsNumber():
                self.value_is_numeric[i] = True
                self.value[i] = settings["value"][i].GetDouble()
            else:
                self.function_string = settings["value"][i].GetString()
                self.aux_function[i] = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])


        # Compute the normal on the nodes of interest
        KratosMultiphysics.NormalCalculationUtils().CalculateOnSimplex(self.model_part, self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE])
        self.modified_normal = False
        if self.option == "flip_normal":
            self.modified_normal = True

        # Set Flag BOUNDARY and variables MATERIAL_POINTS_PER_CONDITION
        if self.material_points_per_condition >= 0:
            KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.BOUNDARY, True, self.model_part.Nodes)

            for condition in self.model_part.Conditions:
                condition.Set(KratosMultiphysics.BOUNDARY, True)
                condition.Set(KratosMultiphysics.SLIP, self.is_slip_boundary)
                condition.Set(KratosMultiphysics.CONTACT, self.is_contact_boundary)
                condition.Set(KratosMultiphysics.MODIFIED, self.modified_normal)
                condition.SetValue(KratosMPM.MATERIAL_POINTS_PER_CONDITION, self.material_points_per_condition)
                condition.SetValue(KratosMPM.IS_EQUAL_DISTRIBUTED, self.is_equal_distributed)
                condition.SetValue(KratosMPM.MPC_IS_NEUMANN, self.is_neumann_boundary)
                condition.SetValue(KratosMPM.MPC_BOUNDARY_CONDITION_TYPE, self.boundary_condition_type)

                ### Set necessary essential BC variables
                if self.boundary_condition_type==1:
                    condition.SetValue(KratosMPM.PENALTY_FACTOR, self.penalty_factor)
        else:
            err_msg = '\n::[ApplyMPMParticleDirichletConditionProcess]:: W-A-R-N-I-N-G: You have specified invalid "material_points_per_condition", '
            err_msg += 'or assigned negative values. \nPlease assign: "material_points_per_condition" > 0 or = 0 (for automatic value)!\n'
            raise Exception(err_msg)

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed in before initialize the solution step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        # Get updated model_part
        if (self.model_part_name.startswith('Background_Grid.')):
            self.model_part_name = self.model_part_name.replace('Background_Grid.','')
        mpm_material_model_part_name = "MPM_Material." + self.model_part_name
        self.model_part = self.model[mpm_material_model_part_name]
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """

        for mpc in self.model_part.Conditions:
            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            mpc_coord = mpc.CalculateOnIntegrationPoints(KratosMPM.MPC_COORD,self.model_part.ProcessInfo)[0]

            if self.interval.IsInInterval(current_time):

                self.step_is_active = True

                # Loop over components X, Y and Z
                for i in range(3):
                    self.variable = self.name[i]
                    if  not self.value_is_numeric[i]:
                        if self.aux_function[i].DependsOnSpace() == False: #depends on time only
                            self.value[i] = self.aux_function[i].CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
                        else: #most general case - space varying function (possibly also time varying)
                            self.value[i] = self.aux_function[i].CallFunction(mpc_coord[0],mpc_coord[1],mpc_coord[2],current_time,0.0,0.0,0.0)


                mpc.SetValuesOnIntegrationPoints(KratosMPM.MPC_IMPOSED_DISPLACEMENT,[self.value],self.model_part.ProcessInfo)
