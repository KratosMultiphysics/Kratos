import numpy
from typing import Optional
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.filtering.filter import Filter
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"ExternalResponseFunctionControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"ExternalResponseFunctionControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return ExternalResponseFunctionControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class ExternalResponseFunctionControl(Control):
    def __init__(self, control_name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(control_name)

        default_parameters = Kratos.Parameters("""{
            "controlled_model_part_names": [],
            "design_variable"            : "CUSTOM_DESIGN_VARIABLE",
            "use_filtering"              : false,
            "filter_settings"            : {},
            "output_all_fields"          : false
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.optimization_problem = optimization_problem

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", parameters["controlled_model_part_names"].GetStringArray(), False)
        self.model_part: 'Optional[Kratos.ModelPart]' = None

        self.design_variable: SupportedSensitivityFieldVariableTypes = Kratos.KratosGlobals.GetVariable(parameters["design_variable"].GetString())
        self.output_all_fields = parameters["output_all_fields"].GetBool()

        if parameters["use_filtering"].GetBool():
            self.filter: 'Optional[Filter]' = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), self.design_variable, Kratos.Globals.DataLocation.NodeNonHistorical, parameters["filter_settings"])
        else:
            self.filter = None


    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

        ComponentDataView(self, self.optimization_problem).SetDataBuffer(1)

        # initialize the filter
        if self.filter:
            self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
            self.filter.Initialize()

        self.physical_phi_field = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, self.design_variable)
        self.physical_phi_field.CollectData()

        if self.filter:
            # take the physical and control field the same
            self.control_phi_field = self.filter.UnfilterField(self.physical_phi_field)
        else:
            self.control_phi_field = self.physical_phi_field.Clone()

        # now update the physical thickness field
        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        if self.filter:
            self.filter.Check()

    def Finalize(self) -> None:
        if self.filter:
            self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [self.design_variable]

    def GetEmptyField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        ta = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, self.design_variable)
        ta.data[:] = 0.0
        return ta

    def GetControlField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        return self.control_phi_field.Clone()

    def MapGradient(self, physical_gradient_variable_tensor_adaptor_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleTensorAdaptor]') -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        with TimeLogger("ExternalResponseFunctionControl::MapGradient", None, "Finished",False):
            keys = physical_gradient_variable_tensor_adaptor_map.keys()
            if len(keys) != 1:
                raise RuntimeError(f"Provided more than required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
            if self.design_variable not in keys:
                raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. {self.design_variable.Name()} not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

            physical_gradient = physical_gradient_variable_tensor_adaptor_map[self.design_variable]
            if physical_gradient.GetContainer() != self.model_part.Nodes:
                raise RuntimeError(f"Gradients for the required container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {physical_gradient.GetModelPart().FullName()} ]")

            if self.filter:
                # now filter the field
                filtered_gradient = self.filter.BackwardFilterIntegratedField(physical_gradient)
            else:
                filtered_gradient = physical_gradient.Clone()

            return filtered_gradient

    def Update(self, new_control_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> bool:
        if new_control_field.GetContainer() != self.model_part.Nodes:
            raise RuntimeError(f"Updates for the required container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()} ]")

        update = new_control_field.Clone()
        update.data[:] -= self.control_phi_field.data
        if numpy.linalg.norm(update.data) > 1e-15:
            with TimeLogger(self.__class__.__name__, f"Updating {self.GetName()}...", f"Finished updating of {self.GetName()}.",False):
                # update the control thickness field
                self.control_phi_field.data[:] = new_control_field.data
                # now update the physical field
                self._UpdateAndOutputFields(update)

                if self.filter:
                    self.filter.Update()
                return True
        return False

    def _UpdateAndOutputFields(self, update: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> None:
        if self.filter:
            # filter the control field
            filtered_update = self.filter.ForwardFilterField(update)
        else:
            filtered_update = update.Clone()

        self.physical_phi_field.data[:] += filtered_update.data
        Kratos.TensorAdaptors.VariableTensorAdaptor(self.physical_phi_field, self.design_variable, copy=False).StoreData()

        # add the values of the control to the optimization problem.
        component_data_view = ComponentDataView(self, self.optimization_problem)
        for i, v in enumerate(self.physical_phi_field.data):
            component_data_view.GetBufferedData().SetValue(f"{self.model_part.FullName()}_point_{i+1}", float(v))

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue(f"{self.GetName()}_{self.design_variable.Name()}", self.physical_phi_field.Clone(), overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue(f"{self.GetName()}_{self.design_variable.Name()}_update", update.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"{self.GetName()}_{self.design_variable.Name()}_filtered_update", filtered_update.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"{self.GetName()}_{self.design_variable.Name()}_control_phi", self.control_phi_field.Clone(), overwrite=True)
