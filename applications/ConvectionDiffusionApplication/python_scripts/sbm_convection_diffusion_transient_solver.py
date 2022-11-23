import KratosMultiphysics
import math
import numpy as np

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_transient_solver
from KratosMultiphysics.ConvectionDiffusionApplication.Functions import*

# from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

# Need for DistanceModificationProcess
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def CreateSolver(main_model_part, custom_settings):
    return SBMConvectionDiffusionTransientSolver(main_model_part, custom_settings)


class SBMConvectionDiffusionTransientSolver(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver):
    print('ci siamo')
    skin_model_part = Import_Structural_model_part('Line')
    # skin_model_part = Import_Structural_model_part('circle')

    def __init__(self, main_model_part, custom_settings):
        super().__init__(main_model_part, custom_settings)
        # * Here I don't have access to the model_part yet
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        print(self.skin_model_part)


    def Initialize(self):
        super().Initialize()
        time_iter = 0
        main_model_part = self.GetComputingModelPart()
        skin_model_part = self.skin_model_part
        # print(self.skin_model_part)
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(main_model_part, skin_model_part).Execute()
        # Find the surrogate boundary nodes
        surrogate_sub_model_part, tot_sur_nodes = Find_surrogate_nodes(main_model_part,time_iter)
        self.surrogate_sub_model_part = surrogate_sub_model_part
        # Total number of skin elements
        tot_skin_el = len(self.skin_model_part.Conditions)
        print('Number of skin elements: ', tot_skin_el)
        # Find the closest skin element for each surr node
        closest_element = Find_closest_skin_element(main_model_part,skin_model_part,tot_sur_nodes,tot_skin_el)
        # Find the projection onto the skin elements for each surr node
        projection_surr_nodes = Find_projections(main_model_part,skin_model_part,tot_sur_nodes,tot_skin_el,closest_element)
        # Then we create a sub_model part with just the elements & the nodes "outside" the surrogate boundary
        sub_model_part_fluid = Create_sub_model_part_fluid(main_model_part,time_iter)
        self.sub_model_part_fluid = sub_model_part_fluid
        print('Creato il sub_model_part per il calcolo del gradiente')

        for elem in main_model_part.Elements :
            if elem.IsNot(MARKER):
                elem.Set(ACTIVE,False)


        self.projection_surr_nodes = projection_surr_nodes
        self.closest_element = closest_element
        self.tot_sur_nodes = tot_sur_nodes
        # Set the BC at the skin mesh________________________________________________________________________________________________
        for node in skin_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.25*(9-node.X**2-node.Y**2-2*math.log(3) + math.log(node.X**2+node.Y**2)) + 0.25 *math.sin(node.X) * math.sinh(node.Y))
            node.Fix(KratosMultiphysics.TEMPERATURE)
        
        # Set the VELOCITY VECTOR FIELD_______________________________________________________________________________________________
        for node in main_model_part.Nodes:
            #node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [0.0, 0.0, 0.0])
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [5.0, 0.0, 0.0])
            #node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [-1/8*node.Y**2 + 5/8 *node.Y -9/32, 0.0, 0.0])
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)
       
        

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        main_model_part = self.GetComputingModelPart()
        # Compute the gradient with the function ComputeNodalGradientProcess
         # Set a fake scalar field -> VELOCITY_X
        # Compute the gradient using the sub_model_part_fluid
        KratosMultiphysics.ComputeNodalGradientProcess(
        self.sub_model_part_fluid,
        KratosMultiphysics.TEMPERATURE,
        KratosMultiphysics.TEMPERATURE_GRADIENT,
        KratosMultiphysics.NODAL_AREA).Execute()

        # for node in self.GetComputingModelPart().Nodes :
        #      if node.Is(VISITED) :
        #         grad = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)
        #         print('gradiente : ', grad)


        # Compute the Dirichlet BC at the surr nodes
        surr_BC = Dirichlet_BC (main_model_part,self.skin_model_part,self.tot_sur_nodes,self.closest_element,self.projection_surr_nodes)
        
        # Impose the embedded BC
        i = 0 
        for node in self.surrogate_sub_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,surr_BC[i])
            node.Fix(KratosMultiphysics.TEMPERATURE)
            i = i+1


        
        # parameters_to_intpose = KratosMultiphysics.Parameters("""
        #         {
        #             "model_part_name" : "ThermalModelPart.surrogate_sub_model_part",
        #             "variable_name"   : "TEMPERATURE",
        #             "interval"        : [0.0, "End"],
        #             "constrained"     : true,
        #             "value"           : 0.0
        #         }
        #         """
        #         )
        # process_applyying_bc = AssignScalarVariableProcess(self.model,parameters_to_intpose)
        # process_applyying_bc.Execute()

    # def FinilizeSolutionStep(self) :
    #     main_model_part = self.GetComputingModelPart()
    #     for node in main_model_part.Nodes :
    #         if node.IsNot(MARKER) :
    #             node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
    #             node.Fix(KratosMultiphysics.TEMPERATURE)
    #     super().FinilizeSolutionStepSolutionStep()


    def Finalize(self):
        super().Finalize()
        # This is compiles after solving all the solution steps
        file_tre = open("Surr_B.txt", "w")
        for node in self.surrogate_sub_model_part.Nodes :
            file_tre.write(str(node.X))
            file_tre.write('  ')
            file_tre.write(str(node.Y))
            file_tre.write('\n')
        file_tre.close()
        main_model_part = self.GetComputingModelPart()

      

        
        # Check the error
        # Compute_error(main_model_part)














































    #____________________________________________________________________________________________________________________________
    #____________________________________________________________________________________________________________________________
    # @classmethod
    # def __GetContinuousDistanceModificationDefaultSettings(cls):
    #     return KratosMultiphysics.Parameters(r'''{
    #         "model_part_name": "",
    #         "distance_threshold": 1e-3,
    #         "continuous_distance": true,
    #         "check_at_each_time_step": true,
    #         "avoid_almost_empty_elements": true,
    #         "deactivate_full_negative_elements": true
    #     }''')

    # def __CreateDistanceModificationProcess(self):
    #     # Set the distance modification settings according to the level set type
    #     # Note that the distance modification process is applied to the volume model part
    #     distance_modification_settings = self.settings["distance_modification_settings"]
    #     distance_modification_settings.ValidateAndAssignDefaults(self.__GetContinuousDistanceModificationDefaultSettings(self.level_set_type))
    #     aux_full_volume_part_name = self.settings["model_part_name"].GetString() + "." + self.settings["volume_model_part_name"].GetString()
    #     distance_modification_settings["model_part_name"].SetString(aux_full_volume_part_name)
    #     return KratosCFD.DistanceModificationProcess(self.model, distance_modification_settings)
    
    # def GetDistanceModificationProcess(self):
    #     if not hasattr(self, '_distance_modification_process'):
    #         self._distance_modification_process = self.__CreateDistanceModificationProcess()
    #     return self._distance_modification_process
    


    #### Private functions ####
    def _Prova_Function(self):
        print('Prova function')
        return
