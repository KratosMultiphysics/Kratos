import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from master.applications.OptimizationApplication.python_scripts.utilities.helper_utilities import IsSameContainerExpression
from numpy import tanh

from applications.OptimizationApplication.python_scripts.utilities.component_data_view import ComponentDataView
from applications.OptimizationApplication.python_scripts.utilities.union_utilities import ContainerExpressionTypes, SupportedSensitivityFieldVariableTypes


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> Control:
    if not parameters.Has("name"):
        raise RuntimeError(f"LevelSetControl instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"LevelSetControl instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return LevelSetControl(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class LevelSetControl(Control):
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
            "k"                                 : 1000.0                                 
                }
            ]
        }""")
        
        #            "density_projection_settings"       : {},
        #            "young_modulus_projection_settings" : {},
        #            "filter_settings"                   : {},

        #fill in missing info from defaults
        parameters.ValidateAndAssignDefaults(default_settings)

        self.output_all_fields = parameters["output_all_fields"].GetBool()
        self.echo_level = parameters["echo_level"].GetInt()
        self.k = parameters["k"].GetDouble()

        self.controlled_physical_variables = [Kratos.DENSITY, Kratos.YOUNG_MODULUS]

        controlled_model_names_parts = parameters["controlled_model_part_names"].GetStringArray()
        if len(controlled_model_names_parts) == 0:
            raise RuntimeError(f"No model parts are provided for SimpControl. [ control name = \"{self.GetName()}\"]")
        # self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"control_{self.GetName()}", controlled_model_names_parts, False)
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
        if not(is_density_defined & is_youngs_modulus_defined):
            raise RuntimeError(f"Elements of {self.model_part.FullName()} does not define either DENSITY or YOUNG_MODULUS.")

        # # get the control field
        # self.control_phi = self.filter.UnfilterField(self.simp_physical_phi)
        self.control_phi = self.GetEmptyField()

        self._UpdateAndOutputFields(self.GetEmptyField())

    def Check(self) -> None:
        #self.filter.Check()
        pass

    def Finalize(self) -> None:
        #self.filter.Finalize()
        pass

    def GetPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return self.controlled_physical_variables

    def GetEmptyField(self) -> ContainerExpressionTypes:
        field = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.LiteralExpressionIO.SetData(field, 0.0)
        return field

    def GetControlField(self) -> ContainerExpressionTypes:
        return self.control_phi.Clone()

    def MapGradient(self, physical_gradient_variable_container_expression_map: 'dict[SupportedSensitivityFieldVariableTypes, ContainerExpressionTypes]') -> ContainerExpressionTypes:
        # get the physical variable
        keys = physical_gradient_variable_container_expression_map.keys()
        
        # check if one gradient is provided and is it density/Young's modulus
        if len(keys) != 1:
            raise RuntimeError(f"Sensitivity w.r.t. DENSITY or YOUNG_MODULUS expected for control \"{self.GetName()}\". \"{len(keys)}\" variables provided.")
        phys_control_var = next(iter(keys))
        if  Kratos.DENSITY not in keys and Kratos.YOUNG_MODULUS not in keys:
            raise RuntimeError(f"The required gradient for control \"{self.GetName()}\" w.r.t. DENSITY or YOUNG_MODULUS not found. Following is the variable:\n \"{phys_control_var.Name()}\"")

        # first access the physical partial sensitivity of the response function w.r.t the control variable (rho or E)
        d_j_d_variable = physical_gradient_variable_container_expression_map[phys_control_var]

        # now compute response function total sensitivity w.r.t. phi
        ### or use d_density_d_phi as a filter?
        d_variable_d_phi = self._ComputePhysicalVariableGradient(phys_control_var)
        d_j_d_phi = d_j_d_variable * d_variable_d_phi

        return d_j_d_phi

    def Update(self, control_field: ContainerExpressionTypes) -> bool:
        # I think here I will get an already updated LSF and I need to use the heavyside to map density and E

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
        
        # rho = H(phi)*rho_0
        density = self._ComputePhysicalVariableField(Kratos.DENSITY)
        KratosOA.PropertiesVariableExpressionIO.Write(density, Kratos.DENSITY)
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(density.GetContainer(), Kratos.DENSITY)
        self.un_buffered_data.SetValue("DENSITY", density.Clone(), overwrite=True)

        youngs_modulus = self._ComputePhysicalVariableField(Kratos.YOUNG_MODULUS)
        KratosOA.PropertiesVariableExpressionIO.Write(youngs_modulus, Kratos.YOUNG_MODULUS)
        if self.consider_recursive_property_update:
            KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(youngs_modulus.GetContainer(), Kratos.YOUNG_MODULUS)
        self.un_buffered_data.SetValue("YOUNG_MODULUS", youngs_modulus.Clone(), overwrite=True)

        # now calculate the total sensitivities of density w.r.t. phi
        self.d_density_d_phi = self._ComputePhysicalVariableGradient(Kratos.DENSITY)

        # now calculate the total sensitivities of young modulus w.r.t. simp_physical_phi
        self.d_young_modulus_d_phi = self._ComputePhysicalVariableGradient(Kratos.YOUNG_MODULUS)

        # now output the fields
        un_buffered_data = ComponentDataView(self, self.optimization_problem).GetUnBufferedData()
        un_buffered_data.SetValue(f"{self.GetName()}", self.control_phi.Clone(),overwrite=True)
        if self.output_all_fields:
            un_buffered_data.SetValue(f"dDENSITY_d{self.GetName()}", self.d_density_d_phi.Clone(), overwrite=True)
            un_buffered_data.SetValue(f"dYOUNG_MODULUS_d{self.GetName()}", self.d_young_modulus_d_phi.Clone(), overwrite=True)


    def _ComputePhysicalVariableField(self, physical_variable) -> ContainerExpressionTypes:
        # Get variable value
        base_value = self.parameters[physical_variable.Name().lower()].GetDouble()
        # Compute Heaviside
        heaviside = 0.5 * (1.0 + tanh(self.k * self.control_phi))
        return base_value * heaviside


    def _ComputeLevelSetGradientNorm(self, LSF: ContainerExpressionTypes) -> ContainerExpressionTypes:
        # Compute gradient (vector expression)
        LSF_grad = KratosOA.GradientComputationUtility().ComputeNodalGradient(self.model_part, LSF)

        # Compute norm of the gradient
        LFS_grad_norm = KratosOA.ContainerExpression[Kratos.ModelPart](self.model_part)
        LFS_grad_norm.EvaluateNormOf(LSF_grad)

        # Return the norm expression
        return LFS_grad_norm
    

    def _ComputePhysicalVariableGradient(self, physical_variable)  -> ContainerExpressionTypes:
        
        # calculate Dirac Delta function (derivative of Heaviside function w.r.t. phi)
        d_H_d_phi = self.k / 2 *(1 - pow(tanh(self.k * self.control_phi), 2))

        # calculate variable sensitivity w.r.t. phi
        d_var_d_phi = d_H_d_phi * physical_variable

        # return the sensitivity
        return d_var_d_phi