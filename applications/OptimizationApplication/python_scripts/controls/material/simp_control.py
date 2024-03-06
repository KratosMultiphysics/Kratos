import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.filtering.filter import Factory as FilterFactory

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

        self.__list_of_densities: 'list[float]' = [0.0]
        self.__list_of_young_modulus: 'list[float]' = [0.0]
        self.__phi: 'list[float]' = [0.0]

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
            self.__phi.append(i+1)
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
            "controlled_model_part_names": [""],
            "simp_power_fac"   : 3,
            "list_of_materials": [
                {
                    "density"      : 1.0,
                    "young_modulus": 1.0
                }
            ],
            "beta_settings": {
                "initial_value": 5,
                "max_value"    : 30,
                "adaptive"     : true,
                "increase_fac" : 1.05,
                "update_period": 50
            },
            "filter_settings"            : {}
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.simp_power_fac = parameters["simp_power_fac"].GetInt()

        self.controlled_physical_variables = [Kratos.DENSITY, Kratos.YOUNG_MODULUS]

        self.materials = Materials(parameters["list_of_materials"].values())

        # beta settings
        beta_settings = parameters["beta_settings"]
        beta_settings.ValidateAndAssignDefaults(default_settings["beta_settings"])
        self.beta = beta_settings["initial_value"].GetDouble()
        self.beta_max = beta_settings["max_value"].GetDouble()
        self.beta_adaptive = beta_settings["adaptive"].GetBool()
        self.beta_increase_frac = beta_settings["increase_fac"].GetDouble()
        self.beta_update_period = beta_settings["update_period"].GetInt()
        self.beta_computed_step = 1

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SimpControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

        # filtering settings
        self.filter = FilterFactory(self.model, self.model_part_operation.GetModelPartFullName(), Kratos.DENSITY, Kratos.Globals.DataLocation.Element, parameters["filter_settings"])

    def Initialize(self) -> None:
        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        # init model parts
        self.model_part = self.model_part_operation.GetModelPart()

        # Creating element specific properties for Youngs Modulus and Density
        if not KratosOA.ModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
            KratosOA.ModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        # calculate phi from existing density
        density = Kratos.Expression.ElementExpression(self.model_part)
        KratosOA.PropertiesVariableExpressionIO.Read(density, Kratos.DENSITY)
        self.phi = KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(density, self.materials.GetPhi(), self.materials.GetDensities(), self.beta, 1)
        self.un_buffered_data.SetValue("phi", self.phi, overwrite=True)
        self.__ApplyDensityAndYoungsModulus()

        self.filter.SetComponentDataView(ComponentDataView(self, self.optimization_problem))
        self.filter.Initialize()

    def Check(self) -> None:
        self.filter.Check()

    def Finalize(self) -> None:
        self.filter.Finalize()

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return self.controlled_physical_variables

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.phi.Clone()

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 2:
            raise RuntimeError(f"Not provided required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if  Kratos.DENSITY not in keys or Kratos.YOUNG_MODULUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. DENSITY or YOUNG_MODULUS not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        # first calculate the density partial sensitivity of the response function
        d_j_d_density = physical_gradient_variable_container_expression_map[Kratos.DENSITY]

        # second calculate the young modulus partial sensitivity of the response function
        d_j_d_youngs = physical_gradient_variable_container_expression_map[Kratos.YOUNG_MODULUS]

        # now calculate the total sensitivities of density w.r.t. phi
        d_density_d_phi = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(self.phi, self.materials.GetPhi(), self.materials.GetDensities(), self.beta, 1)

        # now calculate the total sensitivities of young modulus w.r.t. phi
        d_young_modulus_d_phi = KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(self.phi, self.materials.GetPhi(), self.materials.GetYoungModulus(), self.beta, self.simp_power_fac)

        # now compute response function total sensitivity w.r.t. phi
        d_j_d_phi = d_j_d_density * d_density_d_phi + d_j_d_youngs * d_young_modulus_d_phi

        return self.filter.FilterIntegratedField(d_j_d_phi)

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
            raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        if Kratos.Expression.Utils.NormL2(self.phi - control_field) > 1e-15:
            self.phi = self.filter.FilterField(control_field)
            self.un_buffered_data.SetValue("filtered_phi", self.phi.Clone(), overwrite=True)

            self.__UpdateBeta()

            self.__ApplyDensityAndYoungsModulus()

    def __UpdateBeta(self) -> None:
        if self.beta_adaptive:
            step = self.optimization_problem.GetStep()
            if step % self.beta_update_period == 0 and self.beta_computed_step != step:
                self.beta_computed_step = step
                self.beta = min(self.beta * self.beta_increase_frac, self.beta_max)
                Kratos.Logger.PrintInfo(f"::{self.GetName()}::", f"Increased beta to {self.beta}.")

    def __ApplyDensityAndYoungsModulus(self) -> None:
        density =  KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(self.phi, self.materials.GetPhi(), self.materials.GetDensities(), self.beta, 1)
        KratosOA.PropertiesVariableExpressionIO.Write(density, Kratos.DENSITY)
        self.un_buffered_data.SetValue("DENSITY", density.Clone(), overwrite=True)

        youngs_modulus =  KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(self.phi, self.materials.GetPhi(), self.materials.GetYoungModulus(), self.beta, self.simp_power_fac)
        KratosOA.PropertiesVariableExpressionIO.Write(youngs_modulus, Kratos.YOUNG_MODULUS)
        self.un_buffered_data.SetValue("YOUNG_MODULUS", youngs_modulus.Clone(), overwrite=True)
