# Importing the Kratos Library
import KratosMultiphysics

# Import applications modules
from KratosMultiphysics.FluidDynamicsApplication import python_solvers_wrapper_fluid
from KratosMultiphysics.ConvectionDiffusionApplication import python_solvers_wrapper_convection_diffusion

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings):
    return CoupledFluidThermalSolver(main_model_part, custom_settings)

class CoupledFluidThermalSolver(PythonSolver):

    @classmethod
    def GetDefaultParameters(cls):

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "ThermallyCoupled",
            "domain_size" : -1,
            "echo_level": 0,
            "fluid_solver_settings": {
                "solver_type": "navier_stokes_solver_vmsmonolithic",
                "model_import_settings": {
                    "input_type": "mdpa",
                    "input_filename": "unknown_name"
                }
            },
            "thermal_solver_settings": {
                "solver_type": "transient",
                "analysis_type": "linear",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                    "materials_filename": "ThermalMaterials.json"
                }
            }
        }
        """)

        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings):
        ## Cal base class constructor
        super().__init__(model, custom_settings)

        ## Get domain size
        self.domain_size = self.settings["domain_size"].GetInt()

        ## Create subdomain solvers
        self.fluid_solver = python_solvers_wrapper_fluid.CreateSolverByParameters(self.model, self.settings["fluid_solver_settings"],"OpenMP")
        self.thermal_solver = python_solvers_wrapper_convection_diffusion.CreateSolverByParameters(self.model,self.settings["thermal_solver_settings"],"OpenMP")

    def AddVariables(self):
        # Import the fluid and thermal solver variables. Then merge them to have them in both fluid and thermal solvers.
        self.fluid_solver.AddVariables()
        self.thermal_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.thermal_solver.main_model_part)

    def ImportModelPart(self):
        # Call the fluid solver to import the model part from the mdpa
        self.fluid_solver.ImportModelPart()

        # Save the convection diffusion settings
        convection_diffusion_settings = self.thermal_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

        # Here the fluid model part is cloned to be thermal model part so that the nodes are shared
        element_name, condition_name = self.__GetElementAndConditionNames()
        modeler = KratosMultiphysics.ConnectivityPreserveModeler()
        modeler.GenerateModelPart(
            self.fluid_solver.main_model_part,
            self.thermal_solver.main_model_part,
            element_name,
            condition_name)

        # Set the saved convection diffusion settings to the new thermal model part
        self.thermal_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)

    def PrepareModelPart(self):
        self.fluid_solver.PrepareModelPart()
        self.thermal_solver.PrepareModelPart()

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.thermal_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_thermal = self.thermal_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_thermal)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.thermal_solver.Initialize()

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.thermal_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.thermal_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.thermal_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        return new_time

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        self.thermal_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.thermal_solver.Predict()

    def SolveSolutionStep(self):
        fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        thermal_is_converged = self.thermal_solver.SolveSolutionStep()

        return (fluid_is_converged and thermal_is_converged)

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.thermal_solver.FinalizeSolutionStep()

    def __GetElementAndConditionNames(self):
        ''' Auxiliary function to get the element and condition names for the connectivity preserve modeler call

        This function returns the element and condition names from the domain size and number of nodes.
        Note that throughout all the substitution process a unique element type and condition is assumed.
        Also note that the connectivity preserve modeler call will create standard base elements as these are to
        be substituted by the corresponding ones in the PrepareModelPart call of the thermal solver.
        '''

        ## Get and check domain size
        fluid_model_part = self.fluid_solver.main_model_part
        domain_size = fluid_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        if domain_size not in [2,3]:
            raise Exception("DOMAIN_SIZE is not set in ProcessInfo container.")

        ## Elements
        ## Get the number of nodes from the fluid mesh elements (if there are no elements simplicial are assumed)
        num_nodes_elements = 0
        if (len(fluid_model_part.Elements) > 0):
            for elem in fluid_model_part.Elements:
                num_nodes_elements = len(elem.GetNodes())
                break
        num_nodes_elements = fluid_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_elements)
        if not num_nodes_elements:
            num_nodes_elements = domain_size + 1

        element_name = f"Element{domain_size}D{num_nodes_elements}N"

        ## Conditions
        ## Get the number of nodes from the fluid mesh conditions (if there are no elements simplicial are assumed)
        num_nodes_conditions = 0
        if (len(fluid_model_part.Conditions) > 0):
            for cond in fluid_model_part.Conditions:
                num_nodes_conditions = len(cond.GetNodes())
                break
        num_nodes_conditions = fluid_model_part.GetCommunicator().GetDataCommunicator().MaxAll(num_nodes_conditions)
        if not num_nodes_conditions:
            num_nodes_conditions = domain_size

        aux_condition_name = "LineCondition" if domain_size == 2 else "SurfaceCondition"
        condition_name = f"{aux_condition_name}{domain_size}D{num_nodes_conditions}N"

        return element_name, condition_name