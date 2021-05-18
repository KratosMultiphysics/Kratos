# Importing the Kratos Library
import KratosMultiphysics
# from KratosMultiphysics.python_solver import PythonSolver

import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication.fluid_solver import FluidSolver

def CreateSolver(model, custom_settings):
    return AdjointFluidSolver(model, custom_settings)

class AdjointFluidSolver(FluidSolver):

    def __init__(self, model, settings):
        super(AdjointFluidSolver,self).__init__(model, settings)

        # Overwrite the default buffer size in base FluidSolver
        # TODO: CHECK WHY THE DEFAULT BUFFER SIZE IS 2 IF THE DEFAULT TIME SCHEME IS BOSSAK
        self.min_buffer_size = 2

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_SCALAR_1, self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Adjoint fluid solver DOFs added correctly.")

    def InitializeSolutionStep(self):
        self._GetSolutionStrategy().InitializeSolutionStep()
        self.GetResponseFunction().InitializeSolutionStep()
        self.GetSensitivityBuilder().InitializeSolutionStep()

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        return self._GetSolutionStrategy().SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()
        self.GetResponseFunction().FinalizeSolutionStep()

        self.GetSensitivityBuilder().UpdateSensitivities()
        self.GetSensitivityBuilder().FinalizeSolutionStep()

    def Check(self):
        self._GetSolutionStrategy().Check()

    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self._GetElementNumNodes()
        cond_num_nodes = self._GetConditionNumNodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## If there are no elements and/or conditions, default to triangles/tetra meshes to avoid breaking the ReplaceElementsAndConditionsProcess
        ## This only affects the input name (if there are no elements or conditions to replace, nothing is replaced).
        if elem_num_nodes == 0:
            elem_num_nodes = domain_size + 1
        if cond_num_nodes == 0:
            cond_num_nodes = domain_size

        ## Complete the element name
        # TODO: EXPORT THE ADJOINT FOLLOWING THE CONVENTION. ONCE THIS IS DONE WE CAN USE THE FUNCTION IN THE BASE CLASS
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    def _ComputeDeltaTime(self):
        if self.settings["time_stepping"]["automatic_time_step"].GetBool():
            raise Exception("Automatic time stepping is not supported by adjoint fluid solver.")

        delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        return delta_time

    def _SetNodalProperties(self):
        # Get density and dynamic viscostity from the properties of the first element
        for el in self.main_model_part.Elements:
            rho = el.Properties.GetValue(KratosMultiphysics.DENSITY)
            if rho <= 0.0:
                raise Exception("DENSITY set to {0} in Properties {1}, positive number expected.".format(rho,el.Properties.Id))
            dyn_viscosity = el.Properties.GetValue(KratosMultiphysics.DYNAMIC_VISCOSITY)
            if dyn_viscosity <= 0.0:
                raise Exception("DYNAMIC_VISCOSITY set to {0} in Properties {1}, positive number expected.".format(dyn_viscosity,el.Properties.Id))
            kin_viscosity = dyn_viscosity / rho
            break
        else:
            raise Exception("No fluid elements found in the main model part.")
        # Transfer the obtained properties to the nodes
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.DENSITY, rho, self.main_model_part.Nodes)
        KratosMultiphysics.VariableUtils().SetVariable(KratosMultiphysics.VISCOSITY, kin_viscosity, self.main_model_part.Nodes)

    def GetResponseFunction(self):
        if not hasattr(self, '_response_function'):
            self._response_function = self.__CreateResponseFunction()
        return self._response_function

    def __CreateResponseFunction(self):
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        response_type = self.settings["response_function_settings"]["response_type"].GetString()
        if response_type == "drag":
            if domain_size == 2:
                response_function = KratosCFD.DragResponseFunction2D(
                    self.settings["response_function_settings"]["custom_settings"],
                    self.main_model_part)
            elif domain_size == 3:
                response_function = KratosCFD.DragResponseFunction3D(
                    self.settings["response_function_settings"]["custom_settings"],
                    self.main_model_part)
            else:
                raise Exception("Invalid DOMAIN_SIZE: " + str(domain_size))
        elif response_type == "norm_square":
            response_function = KratosCFD.VelocityPressureNormSquareResponseFunction(
                self.settings["response_function_settings"]["custom_settings"],
                self.model)
        else:
            raise Exception("Invalid response_type: " + response_type + ". Available response functions: \'drag\'.")
        return response_function

    def GetSensitivityBuilder(self):
        if not hasattr(self, '_sensitivity_builder'):
            self._sensitivity_builder = self.__CreateSensitivityBuilder()
        return self._sensitivity_builder

    def __CreateSensitivityBuilder(self):
        response_function = self.GetResponseFunction()
        time_scheme_settings = self.settings["scheme_settings"]
        time_scheme_type = time_scheme_settings["scheme_type"].GetString()

        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        # the schemes used in fluid supports SLIP conditions which rotates element/condition
        # matrices based on nodal NORMAL. Hence, the consistent adjoints also need to
        # rotate adjoint element/condition matrices accordingly and to compute derivatives
        # of rotation matrices as well. Therefore, following schemes are used.
        if (time_scheme_type == "steady"):
            self.sensitivity_builder_scheme = KratosCFD.SimpleSteadySensitivityBuilderScheme(domain_size, domain_size + 1)
        elif (time_scheme_type == "bossak"):
            self.sensitivity_builder_scheme = KratosCFD.VelocityBossakSensitivityBuilderScheme(time_scheme_settings["alpha_bossak"].GetDouble(), domain_size, domain_size + 1)

        sensitivity_builder = KratosMultiphysics.SensitivityBuilder(
            self.settings["sensitivity_settings"],
            self.main_model_part,
            response_function,
            self.sensitivity_builder_scheme
            )
        return sensitivity_builder
