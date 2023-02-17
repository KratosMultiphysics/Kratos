import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ResponseFunctionImplementor
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import WriteCollectiveVariableDataHolderToOptmizationInfo

class SimpControl(Control):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_names": [""],
            "simp_power_fac"  : 3,
            "youngs_modules"  : [],
            "densities"       : [],
            "beta_settings"   : {
                "initial_value": 25,
                "max_value"    : 25,
                "adaptive"     : false,
                "increase_fac" : 1.5,
                "update_period": 20
            }
        }""")
        parameters.RecursivelyValidateAndAssignDefaults(default_settings)

        self.young_modulus = parameters["youngs_modules"].GetVector()
        self.densities =  parameters["densities"].GetVector()
        self.simp_power_fac = parameters["simp_power_fac"].GetInt()
        self.beta = parameters["beta_settings"]["initial_value"].GetDouble()
        self.beta_adaptive = parameters["beta_settings"]["adaptive"].GetBool()
        self.beta_update_period = parameters["beta_settings"]["update_period"].GetInt()
        self.beta_increase_fac = parameters["beta_settings"]["increase_fac"].GetDouble()
        self.beta_max_value = parameters["beta_settings"]["max_value"].GetDouble()

        self.optimization_info = optimization_info
        self.model_parts = [model[model_part_name] for model_part_name in parameters["model_part_names"].GetStringArray()]

        self.phi_collective_container_data = KratosOA.CollectiveVariableDataHolder()

    def ExecuteInitialize(self):
        # creates element specific properties
        if "model_parts_with_element_specific_properties" not in self.optimization_info.GetSolutionStepData(0).keys():
            self.optimization_info["model_parts_with_element_specific_properties"] = []

        for model_part in self.model_parts:
            if not f"{model_part.FullName()}.Elements" in self.optimization_info["model_parts_with_element_specific_properties"]:
                KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model_part, model_part.Elements)
                self.optimization_info["model_parts_with_element_specific_properties"].append(f"{model_part.FullName()}.Elements")

                Kratos.Logger.PrintInfo(self.__class__.__name__, f"Created entity pecific properties for {model_part.FullName()}.")

            # calculate phi from existing density
            phi = KratosOA.ElementPropertiesContainerVariableDataHolder(model_part)
            phi.ReadDataFromContainerVariable(Kratos.DENSITY)
            KratosOA.SimpUtils.ProjectContainerVariableDataHolderBackward(phi, phi, self.densities, self.beta, 1.0)
            self.phi_collective_container_data.AddVariableDataHolder(phi)

            # to be consistent, predict the density from phi
            density = KratosOA.ElementPropertiesContainerVariableDataHolder(model_part)
            KratosOA.SimpUtils.ProjectContainerVariableDataHolderForward(density, phi, self.densities, self.beta, 1)
            density.AssignDataToContainerVariable(Kratos.DENSITY)

            # Now project the young modulus
            young_modulus = KratosOA.ElementPropertiesContainerVariableDataHolder(model_part)
            KratosOA.SimpUtils.ProjectContainerVariableDataHolderForward(young_modulus, phi, self.young_modulus, self.beta, self.simp_power_fac)
            young_modulus.AssignDataToContainerVariable(Kratos.YOUNG_MODULUS)

    def ExecuteInitializeSolutionStep(self):
        # write phi for post processing
        WriteCollectiveVariableDataHolderToOptmizationInfo(
            self.optimization_info,
            self.phi_collective_container_data,
            f"problem_data/control_data/<model_part_name>/{self.GetName()}/phi")

        for phi in self.phi_collective_container_data.GetVariableDataHolders():
            key_refix =  f"problem_data/control_data/{phi.GetModelPart().FullName()}/{self.GetName()}"

            # write density for post processing
            density = KratosOA.ElementPropertiesContainerVariableDataHolder(phi.GetModelPart())
            density.ReadDataFromContainerVariable(Kratos.DENSITY)
            self.optimization_info.SetValue(f"{key_refix}/DENSITY", density)

            # write young modulus for post processing
            yound_modulus = KratosOA.ElementPropertiesContainerVariableDataHolder(phi.GetModelPart())
            yound_modulus.ReadDataFromContainerVariable(Kratos.YOUNG_MODULUS)
            self.optimization_info.SetValue(f"{key_refix}/YOUNG_MODULUS", yound_modulus)

        step = self.optimization_info["step"]
        if self.beta_adaptive:
            if step % self.beta_update_period == 0 and self.beta < self.beta_max_value:
             self.beta *= self.beta_increase_fac
             if self.beta > self.beta_max_value:
                 self.beta = self.beta_max_value

    def CalculateSensitivity(self, response_function: ResponseFunctionImplementor, output_sensitivities: KratosOA.CollectiveVariableDataHolder):
        # clear the container
        output_sensitivities.ClearVariableDataHolders()

        for phi in self.phi_collective_container_data.GetVariableDataHolders():
            # raise RuntimeError(len(self.phi_collective_container_data.GetVariableDataHolders()))
            # first calculate the density partial sensitivity of the response function
            d_j_d_rho = KratosOA.ElementPropertiesContainerVariableDataHolder(phi.GetModelPart())
            response_function.GetStandardizedSensitivity(KratosOA.DENSITY_SENSITIVITY, d_j_d_rho)

            # second calculate the young modulus partial sensitivity of the response function
            d_j_d_e = KratosOA.ElementPropertiesContainerVariableDataHolder(phi.GetModelPart())
            response_function.GetStandardizedSensitivity(KratosOA.YOUNG_MODULUS_SENSITIVITY, d_j_d_e)

            # now calculate the total sensitivities of density w.r.t. phi
            d_rho_d_phi = KratosOA.ElementPropertiesContainerVariableDataHolder(phi.GetModelPart())
            KratosOA.SimpUtils.ProjectContainerVariableDataHolderDerivative(d_rho_d_phi, phi, self.densities, self.beta, 1.0)

            # now calculate the total sensitivities of young modulus w.r.t. phi
            d_e_d_phi = KratosOA.ElementPropertiesContainerVariableDataHolder(phi.GetModelPart())
            KratosOA.SimpUtils.ProjectContainerVariableDataHolderDerivative(d_e_d_phi, phi, self.young_modulus, self.beta, self.simp_power_fac)

            # now compute response function total sensitivity w.r.t. phi
            d_j_d_phi = d_j_d_rho * d_rho_d_phi + d_j_d_e * d_e_d_phi

            # add containers for post processing
            key_prefix = f"problem_data/response_data/{phi.GetModelPart().FullName()}/{response_function.GetName()}/sensitivities"
            self.optimization_info.SetValue(f"{key_prefix}/DENSITY_SENSITIVITY/{self.GetName()}/projection", d_rho_d_phi)
            self.optimization_info.SetValue(f"{key_prefix}/YOUNG_MODULUS_SENSITIVITY/{self.GetName()}/projection", d_e_d_phi)

            output_sensitivities.AddVariableDataHolder(d_j_d_phi)

    def UpdateControl(self, update: KratosOA.CollectiveVariableDataHolder):
        if not self.phi_collective_container_data.IsCompatibleWith(update):
            raise RuntimeError(f"Unsupported update found for SIMP phi [ Given update = {update}, phi = {self.phi_collective_container_data} ].")

        for phi, phi_update in zip(self.phi_collective_container_data.GetVariableDataHolders(), update.GetVariableDataHolders()):
            model_part = phi.GetModelPart()

            # update phi with update
            phi += phi_update

            # now compute and assign new density
            rho = KratosOA.ElementPropertiesContainerVariableDataHolder(model_part)
            KratosOA.SimpUtils.ProjectContainerVariableDataHolderForward(rho, phi, self.densities, self.beta, 1)
            rho.AssignDataToContainerVariable(Kratos.DENSITY)

            # now compute and assign new young modulus
            e = KratosOA.ElementPropertiesContainerVariableDataHolder(model_part)
            KratosOA.SimpUtils.ProjectContainerVariableDataHolderForward(e, phi, self.young_modulus, self.beta, self.simp_power_fac)
            e.AssignDataToContainerVariable(Kratos.YOUNG_MODULUS)
