import KratosMultiphysics
import math
import numpy as np

# Import applications
import KratosMultiphysics.FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication.stokes_solver_monolithic import StokesSolverMonolithic

import KratosMultiphysics.ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication.Functions import*

from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

# Need for DistanceModificationProcess
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def CreateSolver(main_model_part, settings):
    return SBMStabilizedStokesFormulation(main_model_part, settings)



class SBMStabilizedStokesFormulation(StokesSolverMonolithic):
    print('ci siamo')
    skin_model_part = Import_Structural_model_part('circle_small')
    def __init__(self, main_model_part, settings):
        super().__init__(main_model_part, settings)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of StokesSolverMonolithic finished.")

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_X_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY_Y_GRADIENT)
        print(self.skin_model_part)
  

    def Initialize(self):
        super().Initialize()
        main_model_part = self.GetComputingModelPart()
        skin_model_part = self.skin_model_part
        iter = 0 
        # print(self.skin_model_part)
        KratosMultiphysics.CalculateDistanceToSkinProcess2D(main_model_part, skin_model_part).Execute()
        # Find the surrogate boundary nodes
        surrogate_sub_model_part, tot_sur_nodes = Find_surrogate_nodes(main_model_part,iter)
        self.surrogate_sub_model_part = surrogate_sub_model_part
        # Total number of skin elements
        tot_skin_el = len(self.skin_model_part.Conditions)
        print('Number of skin elements: ', tot_skin_el)
        # Find the closest skin element for each surr node
        closest_element = Find_closest_skin_element(main_model_part,skin_model_part,tot_sur_nodes,tot_skin_el)
        # Find the projection onto the skin elements for each surr node
        projection_surr_nodes = Find_projections(main_model_part,skin_model_part,tot_sur_nodes,tot_skin_el,closest_element)
        # Then we create a sub_model part with just the elements & the nodes "outside" the surrogate boundary
        sub_model_part_fluid = Create_sub_model_part_fluid(main_model_part,iter)
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
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X, 0)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y, 0)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.25*(9-node.X**2-node.Y**2-2*math.log(3) + math.log(node.X**2+node.Y**2)) + 0.25 *math.sin(node.X) * math.sinh(node.Y))
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)



    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        main_model_part = self.GetComputingModelPart()
        # Compute the gradient using the sub_model_part_fluid
        KratosMultiphysics.ComputeNodalGradientProcess(
        self.sub_model_part_fluid,
        KratosMultiphysics.VELOCITY_X,
        KratosMultiphysics.VELOCITY_X_GRADIENT,
        KratosMultiphysics.NODAL_AREA).Execute()
        # Same for Velocity_y
        KratosMultiphysics.ComputeNodalGradientProcess(
        self.sub_model_part_fluid,
        KratosMultiphysics.VELOCITY_Y,
        KratosMultiphysics.VELOCITY_Y_GRADIENT,
        KratosMultiphysics.NODAL_AREA).Execute()


        # Compute the Dirichlet BC at the surr nodes
        surr_BC_x, surr_BC_y = Dirichlet_BC_CFD (main_model_part,self.skin_model_part,self.tot_sur_nodes,self.closest_element,self.projection_surr_nodes)
        # Impose the embedded BC
        i = 0 
        for node in self.surrogate_sub_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_X,surr_BC_x[i])
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY_Y,surr_BC_y[i])
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            i = i+1
    


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