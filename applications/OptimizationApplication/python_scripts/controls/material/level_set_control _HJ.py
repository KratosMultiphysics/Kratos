import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import IsSameContainerExpression
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes, SupportedSensitivityFieldVariableTypes
from numpy import exp


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"LevelSetControlHJ instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"LevelSetControlHJ instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return LevelSetControlHJ(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class LevelSetControlHJ(Control):
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
            "density"                           : 1.0,
            "young_modulus"                     : 1.0,
            "k"                                 : 1000.0,
            "initial_phi_value"                 : 0.0
        }""")


        #fill in missing info from defaults
        parameters.ValidateAndAssignDefaults(default_settings)

        self.output_all_fields = parameters["output_all_fields"].GetBool()
        self.echo_level = parameters["echo_level"].GetInt()
        self.k = parameters["k"].GetDouble()
        self.initial_phi = parameters["initial_phi_value"].GetDouble()
        self.controlled_physical_variables = [Kratos.DENSITY, Kratos.YOUNG_MODULUS]

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SimpControl. [ control name = \"{self.GetName()}\"]")
        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
        # self.model_part: 'typing.Optional[Kratos.ModelPart]' = None

        self.consider_recursive_property_update = parameters["consider_recursive_property_update"].GetBool()

    def Initialize(self) -> None:
        self.un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()

        # init model parts
        self.model_part = self.model_part_operation.GetModelPart()

        # # Creating element specific properties for Youngs Modulus and Density
        if not KratosOA.OptAppModelPartUtils.CheckModelPartStatus(self.model_part, "element_specific_properties_created"):
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements, self.consider_recursive_property_update)
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, "element_specific_properties_created")

        # check if the density or Young's modulus is defined
        element: Kratos.Element
        for element in self.model_part.Elements:
            break
        is_density_defined = element.Properties.Has(Kratos.DENSITY)
        is_youngs_modulus_defined = element.Properties.Has(Kratos.YOUNG_MODULUS)
        if is_density_defined:
            if is_youngs_modulus_defined:
                Kratos.Logger.PrintWarning(self.__class__.__name__, f"Elements of {self.model_part.FullName()} defines both DENSITY and YOUNG_MODULUS. Using DENSITY for initial field calculation and ignoring YOUNG_MODULUS.")
            density = Kratos.Expression.ElementExpression(self.model_part)
            KratosOA.PropertiesVariableExpressionIO.Read(density, Kratos.DENSITY)
        elif is_youngs_modulus_defined:
            young_modulus = Kratos.Expression.ElementExpression(self.model_part)
            KratosOA.PropertiesVariableExpressionIO.Read(young_modulus, Kratos.YOUNG_MODULUS)
        else:
            raise RuntimeError(f"Elements of {self.model_part.FullName()} does not define either DENSITY or YOUNG_MODULUS.")

        # # get the control field
        self.control_phi = self.GetEmptyField()
        Kratos.Expression.LiteralExpressionIO.SetData(self.control_phi, self.initial_phi)

        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        #self.filter.Check()
        pass

    def Finalize(self) -> None:
        #self.filter.Finalize()
        pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.DENSITY, Kratos.YOUNG_MODULUS]

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_phi.Clone()

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        # Get physical variables
        keys = physical_gradient_variable_container_expression_map.keys()
        if len(keys) != 2:
            raise RuntimeError(f"Not provided required gradient fields for control \"{self.GetName()}\". Following are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))
        if  Kratos.DENSITY not in keys or Kratos.YOUNG_MODULUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. DENSITY or YOUNG_MODULUS not found. Followings are the variables:\n\t" + "\n\t".join([k.Name() for k in keys]))

        # first calculate the density partial sensitivity of the response function
        d_j_d_density = physical_gradient_variable_container_expression_map[Kratos.DENSITY]

        # second calculate the young modulus partial sensitivity of the response function
        d_j_d_youngs = physical_gradient_variable_container_expression_map[Kratos.YOUNG_MODULUS]

        # now compute response function total sensitivity w.r.t. phi
        d_j_d_phi = d_j_d_density * self.d_density_d_phi + d_j_d_youngs * self.d_young_modulus_d_phi

        return d_j_d_phi

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        # I will get an already updated LSF and I need to use the heavyside to map density and E fields

        if not IsSameContainerExpression(control_field, self.GetEmptyField()):
             raise RuntimeError(f"Updates for the required element container not found for control \"{self.GetName()}\". [ required model part name: {self.model_part.FullName()}, given model part name: {control_field.GetModelPart().FullName()} ]")

        update = control_field - self.control_phi
        if Kratos.Expression.Utils.NormL2(update) > 1e-15:
            self.control_phi = control_field
            self._UpdateAndOutputFields(update)
            return True

        # # level-set update: d_phi = ( -v * grad(phi) ) * d_t
        # # get norm of grad(phi)
        # phi_grad_norm = self._ComputeLevelSetGradientNorm(self.control_phi)

        # # get design velocity v, dj/dE missing
        # design_velocity = - self.d_young_modulus_d_phi

        # # get update value
        # update = design_velocity * phi_grad_norm

        return False

    def _UpdateAndOutputFields(self, update: ContainerExpressionTypes) -> None:
        heaviside = self._ComputeHeaviside()
        # update physical variable fields
        # density = H(phi)*density_0
        density = heaviside * self.parameters["density"].GetDouble()
        KratosOA.PropertiesVariableExpressionIO.Write(density, Kratos.DENSITY)
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(density.GetContainer(), Kratos.DENSITY)
        self.un_buffered_data.SetValue("DENSITY", density.Clone(), overwrite=True)

        # E = H(phi)*E_0
        youngs_modulus = heaviside * self.parameters["young_modulus"].GetDouble()
        KratosOA.PropertiesVariableExpressionIO.Write(youngs_modulus, Kratos.YOUNG_MODULUS)
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(youngs_modulus.GetContainer(), Kratos.YOUNG_MODULUS)
        self.un_buffered_data.SetValue("YOUNG_MODULUS", youngs_modulus.Clone(), overwrite=True)

        # Calculate the Dirac Delta of phi (dH/d_phi)
        dirac_delta = self._ComputeHeavisideGradient()
        # now calculate the total sensitivities of density and E w.r.t. phi
        self.d_density_d_phi = dirac_delta * self.parameters["density"].GetDouble()
        self.d_young_modulus_d_phi = dirac_delta * self.parameters["young_modulus"].GetDouble()

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue(f"{self.GetName()}", self.control_phi.Clone(),overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue(f"{self.GetName()}_update", update.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"dDENSITY_d{self.GetName()}", self.d_density_d_phi.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"dYOUNG_MODULUS_d{self.GetName()}", self.d_young_modulus_d_phi.Clone(), overwrite=True)

    def _ComputeHeaviside(self) -> ContainerExpressionTypes:
        e = self.control_phi.Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(e, exp(1))
        return Kratos.Expression.Utils.Pow(Kratos.Expression.Utils.Pow(e, self.control_phi * self.k * (-2)) + 1, -1)


    def _ComputeHeavisideGradient(self)  -> ContainerExpressionTypes:
        # Calculate Dirac Delta function (derivative of Heaviside function w.r.t. phi)
        e = self.control_phi.Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(e, exp(1))
        return Kratos.Expression.Utils.Pow(Kratos.Expression.Utils.Pow(e, self.control_phi * self.k * 2) + 1, -2) * Kratos.Expression.Utils.Pow(e, self.control_phi * self.k * 2) * self.k * 2


    def _ComputeLevelSetGradientNorm(self, LSF: ContainerExpressionTypes) -> ContainerExpressionTypes:
        # Compute gradient (vector expression)
        LSF_grad = KratosOA.GradientComputationUtility().ComputeNodalGradient(self.model_part, LSF)

        # Compute norm of the gradient
        LFS_grad_norm = KratosOA.ContainerExpression[Kratos.ModelPart](self.model_part)
        LFS_grad_norm.EvaluateNormOf(LSF_grad)

        # Return the norm expression
        return LFS_grad_norm

    @staticmethod
    def ComputeHeavisideValue(phi: ContainerExpressionTypes, k: float) -> ContainerExpressionTypes:  
        e = phi.Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(e, exp(1))
        return Kratos.Expression.Utils.Pow(Kratos.Expression.Utils.Pow(e, phi * k * (-2)) + 1, -1)
    @staticmethod
    def ComputeHevisideGradientValue(phi: ContainerExpressionTypes, k: float) -> ContainerExpressionTypes:
        e = phi.Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(e, exp(1))
        return Kratos.Expression.Utils.Pow(Kratos.Expression.Utils.Pow(e, phi * k * 2) + 1, -2) * Kratos.Expression.Utils.Pow(e, phi * k * 2) * k * 2

    
if __name__ == "__main__":
    import numpy
    # Prep for data
    model = Kratos.Model()
    model_part = model.CreateModelPart("test")
    k = 0.1
    for i in range(10):
        model_part.CreateNewNode(i + 1, 0, 0, 0).SetValue(Kratos.DENSITY, -2 + i / 2.5)

    phi = Kratos.Expression.NodalExpression(model_part)
    Kratos.Expression.VariableExpressionIO.Read(phi, Kratos.DENSITY, is_historical=False)

    # Calculate heaviside with fuction
    haviside_phi = LevelSetControlHJ.ComputeHeavisideValue(phi, k)

    print(f"Heaviside function: \n{haviside_phi.Evaluate()}")

    # manual calculation of heaviside
    numpy_phi = phi.Evaluate()
    numpy_heaviside = 1 / (1 + numpy.exp(numpy_phi * k * (-2)))

    print(f"Numpy Heaviside: \n{numpy_heaviside}")

    # Calculate Heaviside gradient with function
    heaviside_gradient_phi = LevelSetControlHJ.ComputeHevisideGradientValue(phi, k)

    print(f"Heaviside gradient function: \n{heaviside_gradient_phi.Evaluate()}")

    # Calculate Heaviside gradient manually with finite differences
    perturbation = 1e-5
    FD_heaviside_gradient = (LevelSetControlHJ.ComputeHeavisideValue(phi + perturbation, k) - LevelSetControlHJ.ComputeHeavisideValue(phi, k))/perturbation

    print(f"Heaviside gradient FD: \n{FD_heaviside_gradient.Evaluate()}")

    # Calculate Heaviside gradient with numpy
    numpy_heaviside_gradient = numpy.exp(numpy_phi * k * 2) * k * 2 / (1 + numpy.exp(numpy_phi * k * 2))**2

    print(f"Heaviside gradient numpy: \n{numpy_heaviside_gradient}")