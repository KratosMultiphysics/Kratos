# Importing the Kratos Library
import KratosMultiphysics as KM

# CoSimulation imports
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

# Other imports
import numpy as np

class BaseCouplingInterfaceData:
    """This class serves as interface to the data structure (Model and ModelPart)
    that holds the data used during CoSimulation
    """
    def __init__(self, custom_settings, name="default", solver_name="default_solver"):
        custom_settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.settings = custom_settings
        self.name = name
        self.solver_name = solver_name
        self.model_part_name = self.settings["model_part_name"].GetString()

        # checking names
        if self.name =="" or "." in self.name or " " in self.name:
            self._RaiseException('The name cannot be empty, contain whitespaces or "."!')
        if self.model_part_name == "":
            self._RaiseException('No "model_part_name" was specified!')

        # variable used to identify data
        variable_name = self.settings["variable_name"].GetString()
        if variable_name == "":
            self._RaiseException('No "variable_name" was specified!')
        if not KM.KratosGlobals.HasVariable(variable_name):
            # TODO here maybe we could construct a new var if necessary (maybe clashes with delayed app-import ...?)
            self._RaiseException('Variable "{}" does not exist!'.format(variable_name))

        self.variable_type = KM.KratosGlobals.GetVariableType(variable_name)

        admissible_scalar_variable_types = ["Bool", "Integer", "Unsigned Integer", "Double"]
        admissible_vector_variable_types = ["Array"]

        if not self.variable_type in admissible_scalar_variable_types and not self.variable_type in admissible_vector_variable_types:
            self._RaiseException('The input for "variable" "{}" is of variable type "{}" which is not allowed, only the following variable types are allowed:\n{}, {}'.format(variable_name, self.variable_type, ", ".join(admissible_scalar_variable_types), ", ".join(admissible_vector_variable_types)))

        self.variable = KM.KratosGlobals.GetVariable(variable_name)

        self.dtype = GetNumpyDataType(self.variable_type) # required for numpy array creation

        self.is_scalar_variable = self.variable_type in admissible_scalar_variable_types

        # location of data on ModelPart
        self.location = self.settings["location"].GetString()
        admissible_locations = ["node_historical", "node_non_historical", "element", "condition", "model_part"]
        if not self.location in admissible_locations:
            self._RaiseException('"{}" is not allowed as "location", only the following options are possible:\n{}'.format(self.location, ", ".join(admissible_locations)))

        # dimensionality of the data
        self.dimension = self.settings["dimension"].GetInt()

    @staticmethod
    def GetDefaultParameters():
        return KM.Parameters("""{
            "model_part_name" : "",
            "variable_name"   : "",
            "location"        : "node_historical",
            "dimension"       : -1
        }""")

    def IsDefinedOnThisRank(self):
        return False

    def _RaiseException(self, err_msg):
        raise Exception('CouplingInterfaceData "{}" for ModelPart "{}" of SolverWrapper "{}":\n{}'.format(self.name, self.model_part_name, self.solver_name, err_msg))


class CouplingInterfaceData(BaseCouplingInterfaceData):
    """This class serves as interface to the data structure (Model and ModelPart)
    that holds the data used during CoSimulation
    """
    def __init__(self, custom_settings, model, name="default", solver_name="default_solver"):
        super().__init__(custom_settings, name, solver_name)

        self.model = model
        self.model_part = self.model[self.model_part_name]

        # dimensionality of the data
        if self.is_scalar_variable:
            if self.dimension != -1:
                self._RaiseException('"dimension" cannot be specifed for scalar variables!')
            self.dimension = 1 # needed in other places, e.g. for "Size"
        else:
            if self.dimension < 1:
                self._RaiseException('"dimension" has to be specifed for vector variables!')
            else:
                if self.variable_type == "Array" and self.dimension not in [1,2,3]:
                    self._RaiseException('"dimension" can only be 1,2,3 when using variables of type "Array"')
                if not KM.DOMAIN_SIZE in self.model_part.ProcessInfo:
                    cs_tools.cs_print_warning('CouplingInterfaceData', 'No "DOMAIN_SIZE" was specified for ModelPart "{}"'.format(self.model_part_name))
                else:
                    domain_size = self.model_part.ProcessInfo[KM.DOMAIN_SIZE]
                    if domain_size != self.dimension:
                        cs_tools.cs_print_warning('CouplingInterfaceData', '"DOMAIN_SIZE" ({}) of ModelPart "{}" does not match dimension ({})'.format(domain_size, self.model_part_name, self.dimension))

        if self.location == "node_historical" and not self.model_part.HasNodalSolutionStepVariable(self.variable):
            self._RaiseException('"{}" is missing as SolutionStepVariable in ModelPart "{}"'.format(self.variable.Name(), self.model_part_name))

    def __str__(self):
        self_str =  'CouplingInterfaceData:\n'
        self_str += '\tName: "{}"\n'.format(self.name)
        self_str += '\tSolverWrapper: "{}"\n'.format(self.solver_name)
        self_str += '\tModelPart: "{}"\n'.format(self.model_part_name)
        self_str += '\tIsDistributed: {}\n'.format(self.IsDistributed())
        self_str += '\tVariable: "{}"'.format(self.variable.Name())
        if self.is_scalar_variable:
            self_str += ' (Scalar)'
        else:
            self_str += ' (Vector with dimension: {})'.format(self.dimension)
        self_str += '\n\tLocation: "{}"\n'.format(self.location)
        self_str += '\tSize: {}\n'.format(self.Size())

        return self_str

    def PrintInfo(self):
        print(self)

    def GetModelPart(self):
        return self.model_part

    def IsDistributed(self):
        return self.GetModelPart().IsDistributed()

    def IsDefinedOnThisRank(self):
        return self.GetModelPart().GetCommunicator().GetDataCommunicator().IsDefinedOnThisRank()

    def Size(self):
        if self.location == "model_part":
            return 1 * self.dimension
        else:
            return len(self.__GetDataContainer()) * self.dimension

    def GetBufferSize(self):
        # only historical nodal data can store multiple steps!
        if self.location == "node_historical":
            return self.GetModelPart().GetBufferSize()
        else:
            return 1

    def GetHistoricalVariableDict(self):
        # this method returns the historical variable associated to a ModelPart
        # it is intended to be used before the Mesh is read such that the historical variables
        # can be allocated beforehand. This is the reason why the name of the ModelPart is
        # retrieved from the settings and not from the ModelPart itself.
        if self.location == "node_historical":
            return {self.model_part_name : self.variable}
        else:
            return {}

    def GetData(self, solution_step_index=0):
        self.__CheckBufferSize(solution_step_index)

        if self.location == "node_historical":
            data = self.__GetDataFromContainer(self.__GetDataContainer(), GetSolutionStepValue, solution_step_index)
        elif self.location in ["node_non_historical", "element", "condition"]:
            data = self.__GetDataFromContainer(self.__GetDataContainer(), GetValue)
        elif self.location == "model_part":
            var_val = self.GetModelPart()[self.variable]
            if self.is_scalar_variable:
                data = [var_val]
            else:
                data = [var_val[i] for i in range(self.dimension)]

        return np.asarray(data, dtype=self.dtype)

    def SetData(self, new_data, solution_step_index=0):
        self.__CheckBufferSize(solution_step_index)

        # checking size of data
        if len(new_data) != self.Size():
            self._RaiseException("The sizes of the data are not matching, got: {}, expected: {}".format(len(new_data), self.Size()))

        if self.location == "node_historical":
            self.__SetDataOnContainer(self.__GetDataContainer(), SetSolutionStepValue, new_data, solution_step_index)
        elif self.location in ["node_non_historical", "element", "condition"]:
            self.__SetDataOnContainer(self.__GetDataContainer(), SetValue, new_data)
        elif self.location == "model_part":
            if self.is_scalar_variable:
                self.GetModelPart()[self.variable] = new_data[0]
            else:
                if self.variable_type == "Array":
                    vec_value = [0.0, 0.0, 0.0] # Array values require three entries
                    vec_value[:self.dimension] = new_data[:self.dimension] # apply "padding"
                    self.GetModelPart()[self.variable] = vec_value
                else:
                    self.GetModelPart()[self.variable] = new_data

    def __GetDataFromContainer(self, container, fct_ptr, *args):
        if self.is_scalar_variable:
            return [fct_ptr(entity, self.variable, *args) for entity in container]
        else:
            data = []
            for entity in container:
                vals = fct_ptr(entity, self.variable, *args)
                for i in range(self.dimension):
                    data.append(vals[i])
            return data

    def __SetDataOnContainer(self, container, fct_ptr, data, *args):
        if self.is_scalar_variable:
            [fct_ptr(entity, self.variable, *args, value) for entity, value in zip(container, data)]
        else:
            if self.variable_type == "Array":
                vec_value = [0.0, 0.0, 0.0] # Array values require three entries
                for i_entity, entity in enumerate(container):
                    slice_start = i_entity*self.dimension
                    slice_end = slice_start + self.dimension
                    vec_value[:self.dimension] = data[slice_start:slice_end] # apply "padding"
                    fct_ptr(entity, self.variable, *args, vec_value)
            else:
                for i_entity, entity in enumerate(container):
                    slice_start = i_entity*self.dimension
                    slice_end = slice_start + self.dimension
                    fct_ptr(entity, self.variable, *args, data[slice_start:slice_end])

    def __GetDataContainer(self):
        if self.location == "node_historical":
            return self.GetModelPart().GetCommunicator().LocalMesh().Nodes
        elif self.location == "node_non_historical":
            return self.GetModelPart().GetCommunicator().LocalMesh().Nodes
        elif self.location == "element":
            return self.GetModelPart().GetCommunicator().LocalMesh().Elements
        elif self.location == "condition":
            return self.GetModelPart().GetCommunicator().LocalMesh().Conditions

    def __CheckBufferSize(self, solution_step_index):
        if solution_step_index+1 > self.GetBufferSize():
            if self.location == "node_historical":
                self._RaiseException("The buffer-size is not large enough current buffer size: {} | requested solution_step_index: {}!".format(self.GetBufferSize(), solution_step_index+1))
            else:
                self._RaiseException("accessing data from previous steps is only possible with historical nodal data!")


def GetValue(entity, variable):
    return entity.GetValue(variable)

def GetSolutionStepValue(entity, variable, solution_step_index):
    return entity.GetSolutionStepValue(variable, solution_step_index)

def SetValue(entity, variable, value):
    return entity.SetValue(variable, value)

def SetSolutionStepValue(entity, variable, solution_step_index, value):
    return entity.SetSolutionStepValue(variable, solution_step_index, value)

def GetNumpyDataType(variable_type):
    # https://docs.scipy.org/doc/numpy/user/basics.types.html
    dtype_map = {
        "Bool" : np.bool,
        "Integer" : np.intc,
        "Unsigned Integer" : np.uintc,
        "Double" : np.double,
        "Component" : np.double,
        "Array" : np.double,
    }

    return dtype_map[variable_type]
