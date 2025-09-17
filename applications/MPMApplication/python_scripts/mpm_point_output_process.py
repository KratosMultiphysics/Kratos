# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MPMApplication as KratosMPM

# other imports
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def Factory(settings, Model):
    if not isinstance(Model, KratosMultiphysics.Model):
        raise Exception("expected input shall be a Model object")
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MPMPointOutputProcess(Model, settings["Parameters"])

class MPMPointOutputProcess(KratosMultiphysics.OutputProcess):
    """This process writes results from a material point to a file.
    The output can be requested for material point elements and material point conditions.
    """

    def __init__(self, model, params):
        super().__init__()

        self.model = model
        self.interval = KratosMultiphysics.IntervalUtility(params)
        self.params = params
        self.params.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.output_file = None

        self.format = self.params["print_format"].GetString()
        self.search_tolerance = self.params["search_tolerance"].GetDouble()

        # Getting the ModelPart from the Model
        model_part_name = self.params["model_part_name"].GetString()
        if not model_part_name:
            raise Exception('No "model_part_name" was specified!')
        self.model_part = self.model[model_part_name]

        # Retrieving the position of the entity
        point_position = self.params["position"].GetVector()
        if point_position.Size() != 3:
            raise Exception('The position has to be provided with 3 coordinates!')
        self.point = KratosMultiphysics.Point(point_position[0], point_position[1], point_position[2])

        # Retrieving the output variables
        output_var_names = self.params["output_variables"]
        variable_names = [ output_var_names[i].GetString() for i in range(output_var_names.size()) ]
        if not variable_names:
            raise Exception('No variables specified for output!')

        self.output_variables = [ KratosMultiphysics.KratosGlobals.GetVariable(var) for var in variable_names ]

        # Validate types of variables
        for var in self.output_variables:
            if IsDoubleVariable(var) or IsVectorVariable(var) or IsArrayVariable(var):
                continue
            else:
                err_msg  = f"Type of variable {var.Name()} is not valid\n"
                err_msg +=  "It can only be double, component or array3d!"
                raise Exception(err_msg)

    @staticmethod
    def GetDefaultParameters():
        return KratosMultiphysics.Parameters('''{
            "model_part_name"      : "",
            "entity_type"          : "element",
            "interval"             : [0.0, 1e30],
            "position"             : [],
            "output_variables"     : [],
            "search_tolerance"     : 1e-6,
            "print_format"         : "",
            "output_file_settings" : {}
        }''')

    def ExecuteBeforeSolutionLoop(self):
        entity_type = self.params["entity_type"].GetString()

        # Search material point element
        if entity_type == "element":
            entity_id = KratosMPM.BruteForceMaterialPointLocator(self.model_part).FindElement(self.point, self.search_tolerance)
            if entity_id > -1:
                self.material_point = self.model_part.Elements[entity_id]
                entity_coord = self.material_point.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, self.model_part.ProcessInfo)[0]
                self.coordinates = entity_coord

        # Search material point condition
        elif entity_type == "condition":
            entity_id = KratosMPM.BruteForceMaterialPointLocator(self.model_part).FindCondition(self.point, self.search_tolerance)
            if entity_id > -1:
                self.material_point = self.model_part.Conditions[entity_id]
                entity_coord = self.material_point.CalculateOnIntegrationPoints(KratosMPM.MP_COORD, self.model_part.ProcessInfo)[0]
                self.coordinates = entity_coord

        else:
            err_msg  =  'Invalid "entity_type" specified, it can only be:\n'
            err_msg += f'"element", "condition". Found: {entity_type}.'
            raise Exception(err_msg)

        if entity_id == -1:
            warn_msg  = f'No "{entity_type}" was found for input {self.point} and '
            warn_msg += f'tolerance {self.search_tolerance:.12g}. No output is written!'
            KratosMultiphysics.Logger.PrintWarning("MPMPointOutputProcess", warn_msg)

        else:
            file_handler_params = KratosMultiphysics.Parameters(self.params["output_file_settings"])
            file_header = GetFileHeader(entity_type, entity_id, self.point, self.coordinates, self.output_variables)
            self.output_file = TimeBasedAsciiFileWriterUtility(
                self.model_part, file_handler_params, file_header).file

        if self.output_file is not None:
            self.output_file.write(self._GetVariablesNameHeader())

    def IsOutputStep(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        return self.interval.IsInInterval(time)

    def PrintOutput(self):
        time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]
        if self.output_file is not None:
            out = str(time)
            for var in self.output_variables:
                value = self.material_point.CalculateOnIntegrationPoints(var, self.model_part.ProcessInfo)[0]

                if IsArrayVariable(var) or IsVectorVariable(var):
                    out += " " + " ".join( format(v,self.format) for v in value )
                elif IsDoubleVariable(var):
                    out += " " + format(value,self.format)
                else:
                    err_msg  = f"Type of variable {var.Name()} is not valid\n"
                    err_msg +=  "It can only be double, component or array3d!\n"
                    err_msg += f"Found: {type(var)}"
                    raise Exception(err_msg)

            out += "\n"
            self.output_file.write(out)
            self.output_file.flush()

    def ExecuteFinalize(self):
        if self.output_file is not None:
            self.output_file.close()

    def _GetVariablesNameHeader(self):
        header = "# time"
        for var in self.output_variables:
            if IsArrayVariable(var):
                header += " {0}_X {0}_Y {0}_Z".format(var.Name())
            elif IsDoubleVariable(var):
                header += " " + var.Name()
            elif IsVectorVariable(var):
                tmp = self.material_point.CalculateOnIntegrationPoints(var, self.model_part.ProcessInfo)[0]
                header += " " + " ".join( f"{var.Name()}_{i}" for i in range(1,tmp.Size()+1) )
        return header + "\n"

def GetFileHeader(entity_type, entity_id, point, mp_coord, output_variables):
    header  = f"# Material point {entity_type} with Id #{entity_id} "
    header += f"and coordinates ({mp_coord[0]:.12g}, {mp_coord[1]:.12g}, "
    header += f"{mp_coord[2]:.12g}) is the closest to the selected point "
    header += f"having coordinates ({point.X:.12g}, {point.Y:.12g}, "
    header += f"{point.Z:.12g}).\n"
    return header

def IsArrayVariable(var):
    return isinstance(var,KratosMultiphysics.Array1DVariable3)

def IsVectorVariable(var):
    return isinstance(var,KratosMultiphysics.VectorVariable)

def IsDoubleVariable(var):
    return isinstance(var,KratosMultiphysics.DoubleVariable)
