# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.ConvectionDiffusionApplication as KratosCD

# Import applications modules
from KratosMultiphysics.FluidDynamicsApplication import fluid_topology_optimization_solver
from KratosMultiphysics.ConvectionDiffusionApplication import transport_topology_optimization_solver

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(main_model_part, custom_settings, isAdjointSolver = False):
    solver_settings = custom_settings["solver_settings"]
    return FluidTransportTopologyOptimizationSolver(main_model_part, solver_settings, is_adjoint_solver=isAdjointSolver)

class FluidTransportTopologyOptimizationSolver(PythonSolver):

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
            "transport_solver_settings": {
                "solver_type": "transient",
                "analysis_type": "linear",
                "model_import_settings": {
                    "input_type": "use_input_model_part"
                },
                "material_import_settings": {
                    "materials_filename": "TransportMaterials.json"
                }
            }
        }
        """)
        default_settings.AddMissingParameters(super().GetDefaultParameters())
        return default_settings

    def __init__(self, model, custom_settings, is_adjoint_solver=False):
        ## Cal base class constructor
        super().__init__(model, custom_settings)
        ## Define if the adjoint problem
        self._SetAdjointSolver(is_adjoint_solver)
        ## Get domain size
        self._SetDomainSize()
        ## Create subdomain solvers
        self._CreateFluidAndTransportSolvers()

    def _SetDomainSize(self):
        self.domain_size = self.settings["domain_size"].GetInt()
        if not self.settings["fluid_solver_settings"]["solver_settings"].Has("domain_size"):
            self.settings["fluid_solver_settings"]["solver_settings"].AddEmptyValue("domain_size")
        self.settings["fluid_solver_settings"]["solver_settings"]["domain_size"].SetInt(self.domain_size)
        if not self.settings["transport_solver_settings"]["solver_settings"].Has("domain_size"):
            self.settings["transport_solver_settings"]["solver_settings"].AddEmptyValue("domain_size")
        self.settings["transport_solver_settings"]["solver_settings"]["domain_size"].SetInt(self.domain_size)

    def _CreateFluidAndTransportSolvers(self):
        self._CreateFluidSolver()
        self._CreateTransportSolver()

    def _CreateFluidSolver(self):
        self.fluid_solver = fluid_topology_optimization_solver.CreateSolver(self.model, self.settings["fluid_solver_settings"], isAdjointSolver=self.IsAdjoint())

    def _CreateTransportSolver(self):
        self.transport_solver = transport_topology_optimization_solver.CreateSolver(self.model,self.settings["transport_solver_settings"], isAdjointSolver=self.IsAdjoint()) 
    
    def _SetAdjointSolver(self, is_adjoint_solver):
        self.is_adjoint = is_adjoint_solver
        
    def AddVariables(self):
        # Import the fluid and transport solver variables. Then merge them to have them in both fluid and transport solvers.
        self.fluid_solver.AddVariables()
        self.transport_solver.AddVariables()
        KratosMultiphysics.MergeVariableListsUtility().Merge(self.fluid_solver.main_model_part, self.transport_solver.main_model_part)

    def _DefineTransportProperties(self, decay, convection_coefficient):
        self._GetTransportSolver()._DefineProperties(decay, convection_coefficient)

    def ImportModelPart(self, model_parts=None):
        if (self.IsPhysics()):
            # Call the fluid solver to import the model part from the mdpa
            self.fluid_solver.ImportModelPart()
            # Save the convection diffusion settings
            convection_diffusion_settings = self.transport_solver.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS)

            # Here the fluid model part is cloned to be transport model part so that the nodes are shared
            element_name, condition_name = self.__GetElementAndConditionNames()
            modeler = KratosMultiphysics.ConnectivityPreserveModeler()
            modeler.GenerateModelPart(
                self.fluid_solver.main_model_part,
                self.transport_solver.main_model_part,
                element_name,
                condition_name)

            # Set the saved convection diffusion settings to the new transport model part
            self.transport_solver.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, convection_diffusion_settings)
        else:
            fluid_mp = model_parts[0]
            transport_mp = model_parts[1]
            self.fluid_solver.main_model_part = fluid_mp
            self.transport_solver.main_model_part = transport_mp

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        default_settings = self.GetDefaultParameters()
        # print(default_settings.PrettyPrintJsonString())
        # print(self.settings.PrettyPrintJsonString())
        self.settings.ValidateAndAssignDefaults(default_settings)

    
    def PrepareModelPart(self):
        self.fluid_solver.PrepareModelPart()
        self.transport_solver.PrepareModelPart()

    def AddDofs(self):
        self.fluid_solver.AddDofs()
        self.transport_solver.AddDofs()

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.fluid_solver.GetComputingModelPart()
    
    def GetMainModelPart(self):
        return self.fluid_solver.GetMainModelPart()

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        return self.fluid_solver._ComputeDeltaTime()

    def GetMinimumBufferSize(self):
        buffer_size_fluid = self.fluid_solver.GetMinimumBufferSize()
        buffer_size_transport = self.transport_solver.GetMinimumBufferSize()
        return max(buffer_size_fluid, buffer_size_transport)

    def Initialize(self):
        self.fluid_solver.Initialize()
        self.transport_solver.Initialize()

    def Clear(self):
        (self.fluid_solver).Clear()
        (self.transport_solver).Clear()

    def Check(self):
        (self.fluid_solver).Check()
        (self.transport_solver).Check()

    def SetEchoLevel(self, level):
        (self.fluid_solver).SetEchoLevel(level)
        (self.transport_solver).SetEchoLevel(level)
        (self.transport_solver).SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        #NOTE: the cloning is done ONLY ONCE since the nodes are shared
        new_time = self.fluid_solver.AdvanceInTime(current_time)
        if (self.IsPhysics()):
            self.transport_solver.main_model_part.ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_T_STEP] += 1
        else:
            self.transport_solver.main_model_part.ProcessInfo[KratosCD.TRANSPORT_TOP_OPT_ADJ_T_STEP] += 1
        return new_time

    def InitializeSolutionStep(self):
        self.fluid_solver.InitializeSolutionStep()
        self.transport_solver.InitializeSolutionStep()

    def Predict(self):
        self.fluid_solver.Predict()
        self.transport_solver.Predict()

    def SolveSolutionStep(self):
        if (self.IsPhysics()):
            fluid_is_converged = self.fluid_solver.SolveSolutionStep()
            transport_is_converged = self.transport_solver.SolveSolutionStep()
        elif (self.IsAdjoint()):
            transport_is_converged = self.transport_solver.SolveSolutionStep()
            fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        else:
            KratosMultiphysics.Logger.PrintError("WRONG TOPOLOGY OPTIMIZATION STAGE", "Running 'SolveSolutionStep' during the wrong topopology optimization stage.")
            
        return (fluid_is_converged and transport_is_converged)

    def FinalizeSolutionStep(self):
        self.fluid_solver.FinalizeSolutionStep()
        self.transport_solver.FinalizeSolutionStep()

    def __GetElementAndConditionNames(self):
        ''' Auxiliary function to get the element and condition names for the connectivity preserve modeler call
        This function returns the element and condition names from the domain size and number of nodes.
        Note that throughout all the substitution process a unique element type and condition is assumed.
        Also note that the connectivity preserve modeler call will create standard base elements as these are to
        be substituted by the corresponding ones in the PrepareModelPart call of the transport solver.
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

    def _GetFluidSolver(self):
        return self.fluid_solver
    
    def _GetTransportSolver(self):
        return self.transport_solver
    
    def IsAdjoint(self):
        return self.is_adjoint
    
    def IsPhysics(self):
        return not self.IsAdjoint()
    
    def _SetTopologyOptimizationStage(self, top_opt_stage):
        self.GetComputingModelPart().ProcessInfo.SetValue(KratosMultiphysics.TOP_OPT_PROBLEM_STAGE, top_opt_stage)
        self._SetFluidTopologyOptimizationStage(top_opt_stage)
        self._SetTransportTopologyOptimizationStage(top_opt_stage)

    def _GetTopologyOptimizationStage(self):
        return self.GetComputingModelPart().ProcessInfo.GetValue(KratosMultiphysics.TOP_OPT_PROBLEM_STAGE)    
    
    def _SetFluidTopologyOptimizationStage(self, top_opt_stage):
        self._GetFluidSolver()._SetTopologyOptimizationStage(top_opt_stage)

    def _GetFluidTopologyOptimizationStage(self):
        self._GetFluidSolver()._GetTopologyOptimizationStage()

    def _SetTransportTopologyOptimizationStage(self, top_opt_stage):
        self._GetTransportSolver()._SetTopologyOptimizationStage(top_opt_stage)

    def _GetFluidTopologyOptimizationStage(self):
        self._GetTransportSolver()._GetTopologyOptimizationStage()

    
