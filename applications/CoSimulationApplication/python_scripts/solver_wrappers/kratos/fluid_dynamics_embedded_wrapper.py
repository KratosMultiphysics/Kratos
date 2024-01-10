# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KMC
import KratosMultiphysics.MeshMovingApplication as KMM
import KratosMultiphysics.FSIApplication as KMFSI

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing FluidDynamics
if not CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    raise ImportError("The FluidDynamicsApplication is not available!")
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

# Importing FluidDynamics
if not CheckIfApplicationsAvailable("FSIApplication"):
    raise ImportError("The FSIApplication is not available!")
from KratosMultiphysics.FSIApplication.fsi_coupling_interface import FSICouplingInterface

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import data_communicator_utilities

def Create(settings, model, solver_name):
    return FluidDynamicsEmbeddedWrapper(settings, model, solver_name)

class FluidDynamicsEmbeddedWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the FluidDynamicsApplication of Kratos. NOTE: only 2D! for continuous level set! weak coupling!"""

    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)
        self.interface_name = "FSICouplingInterfaceFluid"
        self.structure_interface_submodelpart_name = "Structure.StructureInterface2D_StructureInterface"
        self.previous_displacement = [0.0, 0.0, 0.0]

    def Initialize(self):
        # Create interface model part for embedded fluid solver (distance calculation and FM_ALE)
        self.__CreateInterfaceModelPart()

        super().Initialize()

        # Set domain size of interface model part to the same as fluid
        self._interface_model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.__GetFluidComputingModelPart().ProcessInfo[KM.DOMAIN_SIZE]

        KM.CalculateNodalAreaProcess(self.__GetFluidComputingModelPart(), self.__GetFluidComputingModelPart().ProcessInfo[KM.DOMAIN_SIZE]).Execute()
        self.__UpdateLevelSet(True)

        # Create fsi utility for pressure interpolation from fluid background mesh to interface
        self._partitioned_fsi_utility = KMFSI.PartitionedFSIUtilitiesArray2D()

    def Predict(self):
        # Update level set before predicting the embedded fluid 
        # NOTE should be done after the structure has been predicted and the interface has been moved!
        self.__UpdateLevelSet()
        super().Predict()

    def SolveSolutionStep(self):
        # Move interface model part before solving the embedded fluid and update level set correspondingly
        # NOTE should be done after the structure has been solved and its movement has been transferred to the interface model part
        self.__MoveInterfaceModelPart()
        self.__UpdateLevelSet()
        super().SolveSolutionStep()
        # Calculate resultants of the fluid for the interface model part
        self.__CalculateResultantsOfInterfaceModelPart()

    def _CreateAnalysisStage(self):
        return FluidDynamicsAnalysis(self.model, self.project_parameters)

    def _GetDataCommunicator(self):
        if not KM.IsDistributedRun():
            return KM.ParallelEnvironment.GetDataCommunicator("Serial")

        # now we know that Kratos runs in MPI
        parallel_type = self.project_parameters["problem_data"]["parallel_type"].GetString()

        # first check if the solver uses MPI
        if parallel_type != "MPI":
            return data_communicator_utilities.GetRankZeroDataCommunicator()

        # now we know that the solver uses MPI, only question left is whether to use all ranks or a subset
        if self.project_parameters["solver_settings"]["solver_type"].GetString() == "ale_fluid":
            model_import_settings = self.project_parameters["solver_settings"]["fluid_solver_settings"]["model_import_settings"]
        else:
            model_import_settings = self.project_parameters["solver_settings"]["model_import_settings"]
        self._CheckDataCommunicatorIsConsistentlyDefined(model_import_settings, self.settings["mpi_settings"])

        return super()._GetDataCommunicator()
    
    def __UpdateLevelSet(self, create_distance_process=False):
        # Update level set distances
        if create_distance_process:
            self._distance_to_skin_process = self.__CreateDistanceToSkinProcess()
            self._distance_to_skin_process.Execute()
        else:
            self._distance_to_skin_process.Execute()

        # Extend level set
        if not hasattr(self, '_parallel_distance_calculator'):
            self._parallel_distance_calculator = self.__CreateParallelDistanceCalculator()
            self._parallel_distance_calculator.Execute()
        else:
            self._parallel_distance_calculator.Execute()

    def __CreateDistanceToSkinProcess(self):
        raycasting_relative_tolerance = 1.0e-10
        return KM.CalculateDistanceToSkinProcess2D(
            self.__GetFluidComputingModelPart(),
            self._interface_model_part,
            raycasting_relative_tolerance)
    
    def __CreateParallelDistanceCalculator(self):
        parallel_redistance_settings = KM.Parameters("""{
            "max_levels" : 2,
            "max_distance": 1e12
        }""")
        return KM.ParallelDistanceCalculationProcess2D(
            self.__GetFluidComputingModelPart(),
            parallel_redistance_settings)
    
    def __GetFluidComputingModelPart(self):
        return self._analysis_stage._GetSolver().GetComputingModelPart()
    
    def __CreateInterfaceModelPart(self):
        # Create structure model part for calculating the embedded distances, the FM-ALE utility and data mapping
        structure_model_part = self.model.CreateModelPart(self.interface_name)
        structure_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        #structure_model_part.AddNodalSolutionStepVariable(KM.VELOCITY)  #TODO use for prediction?
        structure_model_part.AddNodalSolutionStepVariable(KM.POSITIVE_FACE_PRESSURE)
        # Add variables needed for rigid body solver  #TODO necessary to add them?
        #structure_model_part.AddSolutionStepVariable(KMC.GLOBAL_DISPLACEMENT)
        #structure_model_part.AddSolutionStepVariable(KMC.RESULTANT_FORCE)
        #structure_model_part.AddSolutionStepVariable(KMC.RESULTANT_MOMENT)
        # Add variable needed for convergence accelerator  
        structure_model_part.AddNodalSolutionStepVariable(KM.REACTION)
        structure_model_part.AddNodalSolutionStepVariable(KM.REACTION_MOMENT)
        #TODO calculate REACTION, POSITIVE_FACE_PRESSURE, REACTION_MOMENT

        # Read the square geometry
        KM.ModelPartIO('square', KM.ModelPartIO.READ | KM.ModelPartIO.SKIP_TIMER).ReadModelPart(structure_model_part)
        self._interface_model_part = structure_model_part
        
    def __MoveInterfaceModelPart(self):
        model_part = self._interface_model_part

        # # Get global displacements
        # displacement = model_part[KMC.GLOBAL_DISPLACEMENT]

        # # Rotation is not considered
        # angle = 0
        # axis = [1,0,0]

        # # Reference point
        # reference_point = [0.0, 0.0, 0.0]

        # # Apply displacement to mesh --> calculates MESH_DISPLACEMENT
        # KMM.MoveModelPart(
        #     model_part,
        #     axis,                 # rotation axis
        #     angle,                # rotation angle
        #     reference_point, # one point of the rotation axis
        #     displacement)         # translation

        # NOTE Only meant for weak coupling, otherwise this is done multiple times per timestep
        # Advance interface model part in time
        t = self.__GetFluidComputingModelPart().ProcessInfo[KM.TIME]
        model_part.CloneTimeStep(t)
        model_part.ProcessInfo[KM.STEP] = self.__GetFluidComputingModelPart().ProcessInfo[KM.STEP]

        # Move the interface model part  # NOTE the time difference is hard-coded here
        #displacement = [0.25*t, 0.0, 0.0]  
        displacement = model_part[KMC.GLOBAL_DISPLACEMENT]  
        dt = 0.05
        velocity = [(u_n-u_nn)/dt for u_n, u_nn in zip(displacement, self.previous_displacement)]
        self.previous_displacement = displacement
        for node in model_part.Nodes:
            node.X = node.X0 + displacement[0]
            node.Y = node.Y0 + displacement[1]
            node.SetSolutionStepValue(KM.DISPLACEMENT, displacement)
            #node.SetSolutionStepValue(KM.VELOCITY, velocity)

        # NOTE embedded solver and FM_ALE take care of EMBEDDED_VELOCITY
            
    def __CalculateResultantsOfInterfaceModelPart(self):
        model_part = self._interface_model_part
        reference_point = [0.0, 0.0, 0.0]

        # Pressure interpolation from fluid background mesh to interface
        #TODO Move to Co-Sim?!
        self._partitioned_fsi_utility.EmbeddedPressureToPositiveFacePressureInterpolator(
            self.__GetFluidComputingModelPart(),
            model_part)
        
        # Calculate REACTION and REACTION_MOMENT from POSITIVE_FACE_PRESSURE
        #TODO Move to Co-Sim?!
        self._partitioned_fsi_utility.CalculateTractionFromPressureValues(
            model_part,
            KM.POSITIVE_FACE_PRESSURE,
            KM.REACTION, 
            True
        )

        # Calculate the total reaction
        reaction_force, reaction_moment = KM.ForceAndTorqueUtils.ComputeEquivalentForceAndTorque(
            model_part,
            reference_point,
            KM.REACTION,
            KM.REACTION_MOMENT
        )

        # Sign is flipped to go from reaction to action (force)
        force = -1*reaction_force
        moment= -1*reaction_moment
        
        # Save the forces in the model part
        model_part[KMC.RESULTANT_FORCE] = force
        model_part[KMC.RESULTANT_MOMENT] = moment

