from optparse import Values
import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

def Factory(settings: KratosMultiphysics.Parameters, model: KratosMultiphysics.Model) -> KratosMultiphysics.Process:
    if not isinstance(model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AssignInitialVelocityToMaterialPointProcess(model, settings["Parameters"])

class AssignInitialVelocityToMaterialPointProcess(KratosMultiphysics.Process):
    def __init__(self, Model: KratosMultiphysics.Model, settings: KratosMultiphysics.Parameters) -> None:
        KratosMultiphysics.Process.__init__(self)

        # "modulus" must be a number
        if settings.Has("modulus") and not settings["modulus"].IsNumber():
            raise Exception('Parameter "modulus" must be a number')

        # "direction" must be a vector with three components
        if settings.Has("direction"):
            if not settings["direction"].IsVector() or settings["direction"].GetVector().Size() != 3:
                raise Exception('Parameter "direction" must be a vector of length 3')

        # "variable_name" is not requires: velocity is assigned to variable "MP_VELOCITY"
        if settings.Has("variable_name"):
            settings.RemoveValue("variable_name")
            wrn_msg  = 'Parameter "variable_name" has been removed and will be ignored: '
            wrn_msg += 'initial velocity is assigned to "MP_VELOCITY".'
            KratosMultiphysics.Logger.PrintWarning("AssignInitialVelocityToMaterialPointProcess", wrn_msg)

        # This parameter was not used and is removed
        if settings.Has("constrained"):
            settings.RemoveValue("constrained")
            wrn_msg = 'Parameter "constrained" has been removed and will be ignored.'
            KratosMultiphysics.Logger.PrintWarning("AssignInitialVelocityToMaterialPointProcess", wrn_msg)
    
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.interval = KratosMultiphysics.IntervalUtility(settings)

        # Get model_part
        model_part_name = settings["model_part_name"].GetString()
        # The actual initial velocity application occurs after the submodelpart is
        # transferred from the Initial_MPM_Material to the MPM_Material in the MaterialPointGeneratorUtility.
        # Therefore we change the prefix from Initial_MPM_Material to MPM_Material.
        if model_part_name.startswith('Initial_MPM_Material.'):
            model_part_name = model_part_name.replace('Initial_','')
        elif not model_part_name.startswith('MPM_Material.'):
            model_part_name = "MPM_Material." + model_part_name
        
        # self.model_part = Model[model_part_name]
        # FIXME: Model[model_part_name] has to be called after Initialize since we are currently
        #        generating MP there (MP generation should be moved to Modelers)
        self.model = Model
        self.model_part_name = model_part_name
        
        
        self.value_is_numeric = [False, False, False]
        self.numeric_value = KratosMultiphysics.Vector(3)
        self.aux_function = ["","",""]
        
        # Loop over components X, Y and Z
        for i, variable in enumerate(["_X", "_Y", "_Z"]): # FIXME: should we make this process more general? (assign flexible variables)
            if settings["value"][i].IsNumber():
                self.value_is_numeric[i] = True
                self.numeric_value[i] = settings["value"][i].GetDouble()
            elif settings["value"][i].IsString():
                self.function_string = settings["value"][i].GetString()
                self.aux_function[i] = KratosMultiphysics.GenericFunctionUtility(self.function_string, settings["local_axes"])
            else:
                raise Exception('Component {} of parameter "value" must be a number or a string'.format(i))


    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters("""
            {
                "model_part_name"      : "please_specify_model_part_name",
                "value"                : [0.0, 0.0, 0.0],
                "interval"             : [0.0, 1e30],
                "local_axes"           : {}
            }
            """)
                    # "origin" : [0.0, 0.0, 0.0]
                    # "axes"   : [ [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]

    def ExecuteBeforeSolutionLoop(self):
        # Assign velocity to MP after solver.Initialize() - only apply once at the beginning!
        # the model part is identified here, AFTER it has been transferred to the MPM_material part!
        
        # FIXME: Model[model_part_name] has to be called after Initialize since we are currently
        #        generating MP there (MP generation should be moved to Modelers)
        self.model_part = self.model[self.model_part_name]
        self.ExecuteInitializeSolutionStep()

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        
        if not self.model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
            if self.interval.IsInInterval(current_time):
                for mp in self.model_part.Elements:
                    mp_coordinate = mp.CalculateOnIntegrationPoints(KratosMPM.MP_COORD,self.model_part.ProcessInfo)[0]

                    # Evaluate values to be assigned to material point
                    current_values = self._EvaluateValues(mp_coordinate, current_time)

                    mp.SetValuesOnIntegrationPoints(KratosMPM.MP_VELOCITY,[current_values],self.model_part.ProcessInfo)

    def _EvaluateValues(self, mp_coordinate, current_time):
        current_values = KratosMultiphysics.Vector(3)
        # Loop over components X, Y and Z
        for i in range(3):
            # Evaluate velocity value
            if self.value_is_numeric[i]:
                current_values[i] = self.numeric_value[i]
            elif not self.value_is_numeric[i]:
                # use local origin (and local axes) if given
                if self.aux_function[i].UseLocalSystem() == True: 
                    if self.aux_function[i].DependsOnSpace() == True:
                        current_values[i] = self.aux_function[i].RotateAndCallFunction(mp_coordinate[0],mp_coordinate[1],mp_coordinate[2],current_time)
                    else: #depends on time only
                        current_values[i] = self.aux_function[i].RotateAndCallFunction(0.0,0.0,0.0,current_time)
                        
                elif self.aux_function[i].UseLocalSystem() == False:
                    if self.aux_function[i].DependsOnSpace() == True:
                        current_values[i] = self.aux_function[i].CallFunction(mp_coordinate[0],mp_coordinate[1],mp_coordinate[2],current_time,0.0,0.0,0.0)
                    else: #depends on time only
                        current_values[i] = self.aux_function[i].CallFunction(0.0,0.0,0.0,current_time,0.0,0.0,0.0)
                
        return current_values