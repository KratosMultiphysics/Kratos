import time
import KratosMultiphysics
import numpy as np
import math

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
    return SBMMovingConvectionDiffusionTransientSolver(main_model_part, custom_settings)


class SBMMovingConvectionDiffusionTransientSolver(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver):
    print('Ci siamo, moving object solver ready!')
    # skin_model_part = Import_Structural_model_part('circle_small')
    skin_model_part = Import_Structural_model_part('bar')

    def __init__(self, main_model_part, custom_settings):
        super().__init__(main_model_part, custom_settings)
        # * Here I don't have access to the model_part yet
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)


    def Initialize(self):
        super().Initialize()
        main_model_part = self.GetComputingModelPart()
        
        # Set the VELOCITY VECTOR FIELD_______________________________________________________________________________________________
        for node in main_model_part.Nodes:
            #node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [0.0, 0.0, 0.0])
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [5.0, 0.0, 0.0])
            #node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [-1/8*node.Y**2 + 5/8 *node.Y -9/32, 0.0, 0.0])
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)
        self.current_time = 0
        self.delta_time = 0.01
        self.iter = 1
        # self.object_velocity = [0.2, 0.0, 0.0]
        self.object_velocity = [1.0, 1.0, 0.0]
        
        

    def InitializeSolutionStep(self):
        main_model_part = self.GetComputingModelPart()
        current_model = KratosMultiphysics.Model()
        current_skin_model_part = current_model.CreateModelPart("CurrentSkinModelPart")
        current_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        for node in self.skin_model_part.Nodes :
            # # Translation **
            # current_skin_model_part.CreateNewNode(node.Id, node.X + self.current_time * self.object_velocity[0], node.Y + self.current_time * self.object_velocity[1] , 0.0)
            # Rotation **
            r = math.sqrt((node.X-2.1)**2 + (node.Y-1.25)**2)
            if node.X != 2.1 :
                if node.X > 2.1: 
                    teta = math.atan((node.Y-1.25)/(node.X-2.1))
                else :
                    teta = math.atan((node.Y-1.25)/(node.X-2.1)) + math.pi
            else :
                teta = math.pi/2
            omega = 10
            current_skin_model_part.CreateNewNode(node.Id, 2.1 + r * math.cos(teta+omega*self.current_time), 1.25 + r * math.sin(teta+omega*self.current_time) , 0.0)
        for cond in self.skin_model_part.Conditions :
            # node1 = cond.GetNodes()[0].Id
            property = self.skin_model_part.GetProperties()[0]
            current_skin_model_part.CreateNewCondition("LineCondition2D2N", cond.Id, [cond.GetNodes()[0].Id, cond.GetNodes()[1].Id],property)

        self.current_time = self.current_time + self.delta_time

        start_time = time.time()
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(main_model_part, current_skin_model_part).Execute()
        print("--> %s seconds for CalculateDistanceToSkinProcess2D" % (time.time() - start_time))
        # Find the surrogate boundary nodes    
        surrogate_sub_model_part, tot_sur_nodes = Find_surrogate_nodes(main_model_part,self.iter)
        self.surrogate_sub_model_part = surrogate_sub_model_part
        # Total number of skin elements
        tot_skin_el = len(current_skin_model_part.Conditions)
        print('Number of skin elements: ', tot_skin_el)
        # Find the closest skin element for each surr node
        closest_element = Find_closest_skin_element(main_model_part,current_skin_model_part,tot_sur_nodes,tot_skin_el)
        # Find the projection onto the skin elements for each surr node
        projection_surr_nodes = Find_projections(main_model_part,current_skin_model_part,tot_sur_nodes,tot_skin_el,closest_element)
        # Then we create a sub_model part with just the elements & the nodes "outside" the surrogate boundary
        sub_model_part_fluid = Create_sub_model_part_fluid(main_model_part,self.iter)
        self.sub_model_part_fluid = sub_model_part_fluid
        print('Creato il sub_model_part per il calcolo del gradiente')
        for elem in main_model_part.Elements :
            if elem.IsNot(MARKER):
                elem.Set(ACTIVE,False)


        # Set the BC at the skin mesh________________________________________________________________________________________________
        for node in current_skin_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.25*(9-node.X**2-node.Y**2-2*math.log(3) + math.log(node.X**2+node.Y**2)) + 0.25 *math.sin(node.X) * math.sinh(node.Y))
            node.Fix(KratosMultiphysics.TEMPERATURE)

        # Compute the gradient with the function ComputeNodalGradientProcess
        # Compute the gradient using the sub_model_part_fluid
        KratosMultiphysics.ComputeNodalGradientProcess(
        self.sub_model_part_fluid,
        KratosMultiphysics.TEMPERATURE,
        KratosMultiphysics.TEMPERATURE_GRADIENT,
        KratosMultiphysics.NODAL_AREA).Execute()

        # Compute the Dirichlet BC at the surr nodes
        surr_BC = Dirichlet_BC (main_model_part,current_skin_model_part,tot_sur_nodes,closest_element,projection_surr_nodes)
        
        # Impose the embedded BC
        i = 0 
        for node in surrogate_sub_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,surr_BC[i])
            node.Fix(KratosMultiphysics.TEMPERATURE)
            i = i+1

        
        super().InitializeSolutionStep()
        for node in main_model_part.Nodes :
            if  node.Is(VISITED):
                node.Set(VISITED, False)
            if node.Is(INTERFACE) :
                node.Set(INTERFACE, False)
        for elem in main_model_part.Elements :
            if elem.Is(MARKER):
                elem.Set(MARKER,False)
            if elem.IsNot(ACTIVE) :
                elem.Set(ACTIVE,True)
        name_surr_file = "Surr_B" + "_" + str(self.iter) + ".txt"
        # This is compiles after solving all the solution steps
        file_tre = open(name_surr_file, "w")
        for node in surrogate_sub_model_part.Nodes :
            file_tre.write(str(node.X))
            file_tre.write('  ')
            file_tre.write(str(node.Y))
            file_tre.write('\n')
        file_tre.close()
        self.iter = self.iter+1



    def FinalizeSolutionStep(self) :
        main_model_part = self.GetComputingModelPart()
        for node in main_model_part.Nodes:
            node.Free(KratosMultiphysics.TEMPERATURE)
        super().FinalizeSolutionStep()

    # def Finalize(self):
    #     super().Finalize()
        




        # main_model_part = self.GetComputingModelPart()
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
