import typing
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory
from KratosMultiphysics.OptimizationApplication.utilities.opt_projection import CreateProjection

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"SimpControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"SimpControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return SimpControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class Materials:
    def __init__(self, parameters: 'list[Kratos.Parameters]') -> None:
        default_parameters = Kratos.Parameters("""{
            "density"      : 1.0,
            "young_modulus": 1.0
        }""")

        self.__list_of_densities: 'list[float]' = []
        self.__list_of_young_modulus: 'list[float]' = []
        self.__phi: 'list[float]' = []

        data: 'list[tuple[float, float]]' = []
        for params in parameters:
            params.ValidateAndAssignDefaults(default_parameters)
            data.append((params["density"].GetDouble(), params["young_modulus"].GetDouble()))

        # sort in ascending order of densities
        data = sorted(data, key=lambda x: x[0])

        # now check whether the young modulus is also arranged in the ascending order
        for i, (_, young_modulus) in enumerate(data[:-1]):
            if young_modulus > data[i+1][1]:
                raise RuntimeError("Young modulus and densities are not in the ascending order.")

        for i, (density, young_modulus) in enumerate(data):
            self.__phi.append(i)
            self.__list_of_densities.append(density)
            self.__list_of_young_modulus.append(young_modulus)

    def GetDensities(self) -> 'list[float]':
        return self.__list_of_densities

    def GetYoungModulus(self) -> 'list[float]':
        return self.__list_of_young_modulus

    def GetPhi(self) -> 'list[float]':
        return self.__phi

class SimpControl(Control):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> None:
        super().__init__(name)

        self.parameters = parameters
        self.model = model
        self.optimization_problem = optimization_problem

        default_settings = Kratos.Parameters("""{
            "controlled_model_part_names"       : [""],
            "output_all_fields"                 : true,
            "echo_level"                        : 0,
            "consider_recursive_property_update": false,
            "density_projection_settings"       : {},
            "young_modulus_projection_settings" : {},
            "filter_settings"                   : {},
            "list_of_materials"                 : [
                {
                    "density"      : 1.0,
                    "young_modulus": 1.0
                }
            ]
        }""")

        parameters.ValidateAndAssignDefaults(default_settings)

        self.output_all_fields = parameters["output_all_fields"].GetBool()
        self.echo_level = parameters["echo_level"].GetInt()

        self.controlled_physical_variables = [Kratos.DENSITY, Kratos.YOUNG_MODULUS]

        self.materials = Materials(parameters["list_of_materials"].values())

        self.density_projection = CreateProjection(parameters["density_projection_settings"], self.optimization_problem)
        self.young_modulus_projection = CreateProjection(parameters["young_modulus_projection_settings"], self.optimization_problem)

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SimpControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

        self.consider_recursive_property_update = parameters["consider_recursive_property_update"].GetBool()

        # filtering settings
        self.filter = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), Kratos.DENSITY, Kratos.Globals.DataLocation.Element, parameters["filter_settings"])

    def Initialize(self) -> None:
        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        # init model parts
        self.model_part = self.model_part_operation.GetModelPart()

        # Creating element specific properties for Youngs Modulus and Density
        if not KratosOA.OptAppModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements, self.consider_recursive_property_update)
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
        self.filter.Initialize()

        # initialize the projections
        self.density_projection.SetProjectionSpaces(self.materials.GetPhi(), self.materials.GetDensities())
        self.young_modulus_projection.SetProjectionSpaces(self.materials.GetPhi(), self.materials.GetYoungModulus())

        # check if the density or Young's modulus is defined
        element: Kratos.Element
        for element in self.model_part.Elements:
            break
        is_density_defined = element.Properties.Has(Kratos.DENSITY)
        is_youngs_modulus_defined = element.Properties.Has(Kratos.YOUNG_MODULUS)
        if is_density_defined:
            if is_youngs_modulus_defined:
                Kratos.Logger.PrintWarning(self.__class__.__name__, f"Elements of {self.model_part.FullName()} defines both DENSITY and YOUNG_MODULUS. Using DENSITY for initial field calculation and ignoring YOUNG_MODULUS.")
            density = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(self.model_part.Elements, Kratos.DENSITY)
            density.CollectData()
            self.simp_physical_phi = self.density_projection.ProjectBackward(density)
        elif is_youngs_modulus_defined:
            young_modulus = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(self.model_part.Elements, Kratos.YOUNG_MODULUS)
            young_modulus.CollectData()
            self.simp_physical_phi = self.young_modulus_projection.ProjectBackward(young_modulus)
        else:
            raise RuntimeError(f"Elements of {self.model_part.FullName()} does not define either DENSITY or YOUNG_MODULUS.")

        # get the control field
        self.control_phi = self.filter.UnfilterField(self.simp_physical_phi)

        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        self.filter.Check()

    def Finalize(self) -> None:
        self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return self.controlled_physical_variables

    def GetEmptyField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        field = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Elements, Kratos.DENSITY)
        field.data[:] = 0.0
        return Kratos.TensorAdaptors.DoubleTensorAdaptor(field, copy=False)

    def GetControlField(self) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        return Kratos.TensorAdaptors.DoubleTensorAdaptor(self.control_phi) # returning a copy

    def MapGradient(self, physical_gradient_variable_tensor_adaptor_map: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleTensorAdaptor]') -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        keys = physical_gradient_variable_tensor_adaptor_map.keys()
        if len(keys) != 2:
            raise RuntimeError(f"Not provided required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if  Kratos.DENSITY not in keys or Kratos.YOUNG_MODULUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. DENSITY or YOUNG_MODULUS not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        # first calculate the density partial sensitivity of the response function
        d_j_d_density = physical_gradient_variable_tensor_adaptor_map[Kratos.DENSITY]

        # second calculate the young modulus partial sensitivity of the response function
        d_j_d_youngs = physical_gradient_variable_tensor_adaptor_map[Kratos.YOUNG_MODULUS]

        # now compute response function total sensitivity w.r.t. phi
        d_j_d_phi = Kratos.TensorAdaptors.DoubleTensorAdaptor(d_j_d_density)
        d_j_d_phi.data = d_j_d_density.data * self.d_density_d_phi.data + d_j_d_youngs.data * self.d_young_modulus_d_phi.data

        return self.filter.BackwardFilterIntegratedField(d_j_d_phi)

    def Update(self, control_field: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> bool:
        if control_field.GetContainer() != self.model_part.Elements:
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()} ]")

        update = Kratos.TensorAdaptors.DoubleTensorAdaptor(control_field)
        update.data = control_field.data - self.control_phi.data
        if numpy.linalg.norm(update.data) > 1e-15:
            self.control_phi.data[:] = control_field.data[:]
            self._UpdateAndOutputFields(update)
            self.filter.Update()
            self.density_projection.Update()
            self.young_modulus_projection.Update()
            return True

        self.density_projection.Update()
        self.young_modulus_projection.Update()
        return False

    def _UpdateAndOutputFields(self, update: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> None:
        # filter the control field
        physical_phi_update = self.filter.ForwardFilterField(update)
        self.simp_physical_phi.data += physical_phi_update.data

        density = self.density_projection.ProjectForward(self.simp_physical_phi)
        KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(density, Kratos.DENSITY, copy=False).StoreData()
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(density.GetContainer(), Kratos.DENSITY)
        self.un_buffered_data.SetValue("DENSITY", density, overwrite=True)

        youngs_modulus = self.young_modulus_projection.ProjectForward(self.simp_physical_phi)
        KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(youngs_modulus, Kratos.YOUNG_MODULUS, copy=False).StoreData()
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(youngs_modulus.GetContainer(), Kratos.YOUNG_MODULUS)
        self.un_buffered_data.SetValue("YOUNG_MODULUS", youngs_modulus, overwrite=True)

        # now calculate the total sensitivities of density w.r.t. phi
        self.d_density_d_phi = self.density_projection.ForwardProjectionGradient(self.simp_physical_phi)

        # now calculate the total sensitivities of young modulus w.r.t. simp_physical_phi
        self.d_young_modulus_d_phi = self.young_modulus_projection.ForwardProjectionGradient(self.simp_physical_phi)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue(f"{self.GetName()}", self.simp_physical_phi.Clone(),overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue(f"{self.GetName()}_update", physical_phi_update, overwrite=True)
            un_buffered_data.SetValue(f"dDENSITY_d{self.GetName()}", self.d_density_d_phi.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"dYOUNG_MODULUS_d{self.GetName()}", self.d_young_modulus_d_phi.Clone(), overwrite=True)
