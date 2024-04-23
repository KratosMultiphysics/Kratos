# from inspect import Parameter
import KratosMultiphysics
import math
import numpy as np
import KratosMultiphysics.KratosUnittest as KratosUnittest
import pdb

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver
from KratosMultiphysics.ConvectionDiffusionApplication.Functions import*

from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

# Need for DistanceModificationProcess
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def CreateSolver(main_model_part, custom_settings):
    return SBMMovingConvectionDiffusionStationarySolver(main_model_part, custom_settings)


class SBMMovingConvectionDiffusionStationarySolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    # skin_model_part = Import_Structural_model_part('Structural_circle1x1')
    # skin_model_part = Import_Structural_model_part('Structural_rectangle_vertical')
    # skin_model_part = Import_Structural_model_part('Drawings')
    skin_model_part = Import_Structural_model_part('Structural_halfvolve')
    # skin_model_part = Import_Structural_model_part('Structural_sinusoidal')


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
        
        # self.time_step = 0.000625
        # self.time_step = 0.00125
        self.time_step =0.001
        self.time = 0.0
        self.object_velocity = [0.0, 2.00, 0.0]

        self.iter = 0
        main_model_part = self.GetComputingModelPart()
        for node in main_model_part.Nodes :
            # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [40.0, 5.0, 0.0])
            # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [-20.0, -20.0, 0.0])
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [10.0, 0.0, 0.0])
            node.Fix(KratosMultiphysics.VELOCITY_X)
            node.Fix(KratosMultiphysics.VELOCITY_Y)
            node.Fix(KratosMultiphysics.VELOCITY_Z)
            if node.Is(BOUNDARY) : 
                node.Set(BOUNDARY, False)
        skin_model_part = self.skin_model_part
        # Total number of skin elements
        tot_skin_el = len(self.skin_model_part.Conditions)
        print('Number of skin elements: ', tot_skin_el)

        KratosMultiphysics.CalculateDistanceToSkinProcess2D(main_model_part, skin_model_part).Execute()
        # Find the surrogate boundary nodes
        a = KratosMultiphysics.FindSurrogateNodesProcess2D(main_model_part, skin_model_part)
        a.Execute()
        closest_element = a.FindClosestElement(main_model_part, skin_model_part)

        surrogate_sub_model_part = main_model_part.CreateSubModelPart("surrogate_sub_model_part")
        tot_sur_nodes = 0 
        for node in main_model_part.Nodes :
            if node.Is(BOUNDARY):
                surrogate_sub_model_part.AddNode(node,0)
                tot_sur_nodes = tot_sur_nodes + 1

        # Find the projection onto the skin elements for each surr node
        projection_surr_nodes = Find_projections(main_model_part,skin_model_part,tot_sur_nodes,closest_element)
        # Then we create a sub_model part with just the elements & the nodes "outside" the surrogate boundary
        sub_model_part_fluid = Create_sub_model_part_fluid(main_model_part,iter)



        # Set the BC at the skin mesh________________________________________________________________________________________________
        for node in skin_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.0)
            node.Fix(KratosMultiphysics.TEMPERATURE)

        # Calculate the required neighbours
        nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(main_model_part)
        nodal_neighbours_process.Execute()

        ## Find the so called "INTERFACE" elements
        for elem in sub_model_part_fluid.Elements :
            if elem.Is(BOUNDARY) :
                count_surr = 0
                # Let's count how many surrogate nodes the element has
                for node in elem.GetNodes() :
                    if node.Is(BOUNDARY) :
                        count_surr = count_surr + 1
                if count_surr > 1 :  # two or three nodes are surrogate nodes
                    elem.Set(INTERFACE, True) # FUNDAMENTAL

        elemental_neighbours_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(main_model_part)
        elemental_neighbours_process.Execute()


        ## Compute the gradint coefficients for each of the surrogate node
        result = ComputeGradientCoefficients (sub_model_part_fluid, self.model, surrogate_sub_model_part)
        

        ## Compute the T matrix for imposition of sbm condition
        i = 0
        for node in surrogate_sub_model_part.Nodes :
            # Create the T matrices: 0 , 1 , ... , len(boundary_sub_model_part.Nodes)-1
            nameT = "T_" + str(i)
            globals()[nameT] = Compute_T_matrix (result, node, projection_surr_nodes, i)
            i = i + 1

        # Create the CreateMasterSlaveConstraints
        j = 1
        for node in surrogate_sub_model_part.Nodes :
            name = "T_" + str(j-1)
            T = globals()[name]
            if node.Y != 0.0 and node.Y != 1.0:
                Impose_MPC_Globally (main_model_part, result, self.skin_model_part, closest_element, projection_surr_nodes, T, node, j)
            j = j + 1
        
    def InitializeSolutionStep(self):
        self.time = self.time + self.time_step
        self.iter = self.iter + 1
        main_model_part = self.GetComputingModelPart()
        for node in main_model_part.Nodes :
            if  node.IsNot(ACTIVE):
                node.Set(ACTIVE, True)
            if  node.Is(BOUNDARY):
                node.Set(BOUNDARY, False)
            if  node.Is(SLAVE):
                node.Set(SLAVE, False)
            if node.Is(INTERFACE) :
                node.Set(INTERFACE, False)
            if node.X == 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0.0)
                node.Fix(KratosMultiphysics.TEMPERATURE)
            # if node.Y == 0.0:
            #     node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,0.0)
            #     node.Fix(KratosMultiphysics.TEMPERATURE)
        for elem in main_model_part.Elements :
            if elem.IsNot(ACTIVE) :
                elem.Set(ACTIVE,True)
            if elem.Is(BOUNDARY):
                elem.Set(BOUNDARY,False)
            if elem.Is(INTERFACE):
                elem.Set(INTERFACE,False)
            if elem.Is(MARKER):
                elem.Set(MARKER,False)
            if elem.Is(VISITED) :
                elem.Set(VISITED,False)
    
        ## Remove MPC and regreate them
        # print(main_model_part.GetRootModelPart())
        for constraint in main_model_part.MasterSlaveConstraints :
            constraint.Set(KratosMultiphysics.TO_ERASE)
        main_model_part.RemoveMasterSlaveConstraintsFromAllLevels(KratosMultiphysics.TO_ERASE)

        current_model = KratosMultiphysics.Model()
        current_skin_model_part = current_model.CreateModelPart("CurrentSkinModelPart")
        current_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        CenterRotation = [0.3566987298, 0.39820508]
        for node in self.skin_model_part.Nodes :
            ## Translation **
            # current_skin_model_part.CreateNewNode(node.Id, node.X + self.time * self.object_velocity[0], -0.7 +node.Y + self.time * self.object_velocity[1] , 0.0)
            
            ##Rotation ** 
            deltaTheta = -math.pi/3/95*(self.iter-1)
            # deltaTheta = -0.0125*(self.iter-1)
            theta = math.atan((node.Y-CenterRotation[1]) / (node.X-CenterRotation[0]))
            radius = math.sqrt((node.X-CenterRotation[0])**2+(node.Y-CenterRotation[1])**2)
            if node.X >= CenterRotation[0] :
                current_skin_model_part.CreateNewNode(node.Id, CenterRotation[0]+radius * math.cos(theta+deltaTheta),  CenterRotation[1]+radius * math.sin(theta+deltaTheta), 0.0)
            else:
                current_skin_model_part.CreateNewNode(node.Id, CenterRotation[0]-radius * math.cos(theta+deltaTheta),  CenterRotation[1]-radius * math.sin(theta+deltaTheta), 0.0)
            
            # ## Deformation **
            # r = 0.1526368769
            # aspectRatios = 1+(self.iter-1)/4
            # a = math.sqrt(r**2 * aspectRatios)
            # b = r**2 / a
            # if node.X != 0.5 :
            #     theta = math.atan((node.Y-0.5) / (node.X-0.5))
            # else :
            #     if node.Y > 0.5 :
            #         theta = math.pi/2
            #     else :
            #         theta = -math.pi/2
            # x = (a*b)/math.sqrt(b**2*(math.cos(theta))**2+a**2*(math.sin(theta))**2) * math.cos(theta)
            # y = (a*b)/math.sqrt(b**2*(math.cos(theta))**2+a**2*(math.sin(theta))**2) * math.sin(theta)
            # if node.X >= 0.5:
            #     current_skin_model_part.CreateNewNode(node.Id, 0.5+x,  0.5+y, 0.0)
            # else :
            #     current_skin_model_part.CreateNewNode(node.Id, 0.5-x,  0.5-y, 0.0)

            ## Sinusoidal boundary
            # omega = 30
            # Amplitude = 2*0.2*(self.iter/4/8-1/4/8)/5
            # if node.X < 1.0 :
            #     current_skin_model_part.CreateNewNode(node.Id, node.X + Amplitude*math.sin(omega*node.Y+self.iter/4/4), node.Y , 0.0)
            # else :
            #     current_skin_model_part.CreateNewNode(node.Id, node.X, node.Y , 0.0)


        for cond in self.skin_model_part.Conditions :
            property = self.skin_model_part.GetProperties()[0]
            current_skin_model_part.CreateNewCondition("LineCondition2D2N", cond.Id, [cond.GetNodes()[0].Id, cond.GetNodes()[1].Id],property)
        
        # self.skin_model_part = current_skin_model_part

        KratosMultiphysics.CalculateDistanceToSkinProcess2D(main_model_part, current_skin_model_part).Execute()
        # Find the surrogate boundary nodes
        a = KratosMultiphysics.FindSurrogateNodesProcess2D(main_model_part, current_skin_model_part)
        a.Execute()
        closest_element = a.FindClosestElement(main_model_part, current_skin_model_part)

        name_surrogate_sub_model_part = "surrogate_sub_model_part" + "_" + str(self.iter)
        surrogate_sub_model_part = main_model_part.CreateSubModelPart(name_surrogate_sub_model_part)
        tot_sur_nodes = 0 
        for node in main_model_part.Nodes :
            if node.Is(BOUNDARY):
                surrogate_sub_model_part.AddNode(node,0)
                tot_sur_nodes = tot_sur_nodes + 1
        # Find the projection onto the skin elements for each surr node
        projection_surr_nodes = Find_projections(main_model_part,current_skin_model_part,tot_sur_nodes,closest_element)
        # Then we create a sub_model part with just the elements & the nodes "outside" the surrogate boundary
        sub_model_part_fluid = Create_sub_model_part_fluid(main_model_part,self.iter)
        print('Creato il sub_model_part per il calcolo del gradiente')

        # Set the BC at the skin mesh________________________________________________________________________________________________
        for node in current_skin_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.0)
            node.Fix(KratosMultiphysics.TEMPERATURE)
        
        # Calculate the required neighbours
        nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(main_model_part)
        nodal_neighbours_process.Execute()
        ## Find the so called "INTERFACE" elements
        for elem in sub_model_part_fluid.Elements :
            if elem.Is(BOUNDARY) :
                count_surr = 0
                # Let's count how fzmany surrogate nodes the element has
                for node in elem.GetNodes() :
                    if node.Is(BOUNDARY) :
                        count_surr = count_surr + 1
                if count_surr > 1 :  # two or three nodes are surrogate nodes
                    elem.Set(INTERFACE, True) # FUNDAMENTAL
        elemental_neighbours_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(main_model_part)
        elemental_neighbours_process.Execute()
        ## Compute the gradint coefficients for each of the surrogate node
        result = ComputeGradientCoefficients (sub_model_part_fluid, self.model, surrogate_sub_model_part)
        ## Compute the T matrix for imposition of sbm condition
        i = 0
        for node in surrogate_sub_model_part.Nodes :
            # Create the T matrices: 0 , 1 , ... , len(boundary_sub_model_part.Nodes)-1
            nameT = "T_" + str(i)
            globals()[nameT] = Compute_T_matrix (result, node, projection_surr_nodes, i)
            i = i + 1

        ## Create the CreateMasterSlaveConstraints
        j = 1
        for node in surrogate_sub_model_part.Nodes :
            name = "T_" + str(j-1)
            T = globals()[name]
            if node.Y != 0.0 and node.Y != 1.0 :
                Impose_MPC_Globally (main_model_part, result, current_skin_model_part, closest_element, projection_surr_nodes, T, node, j)
            j = j + 1
        
        if math.fmod(self.iter,1) == 0 :
            name_surr_file = "SURR_FOLDER/Surr_B" + "_" + str(self.iter) + ".txt"
            # This is compiles after solving all the solution steps
            file_tre = open(name_surr_file, "w")
            for node in surrogate_sub_model_part.Nodes :
                file_tre.write(str(node.X))
                file_tre.write('  ')
                file_tre.write(str(node.Y))
                file_tre.write('\n')
            file_tre.close()

        for node in main_model_part.Nodes :
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0 :
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
                node.Fix(KratosMultiphysics.TEMPERATURE)

        self.current_skin_model_part = current_skin_model_part
        # Initialize solution step (Alla fine!)
        super().InitializeSolutionStep()
        

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        # Unfix the distance<0 nodes
        for node in self.main_model_part.Nodes :
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0 :
                node.Free(KratosMultiphysics.TEMPERATURE)
        
        ## Sebastian check for MPC
        self.CheckIfMPCsAreAppliedCorrectly()
        
        ## White elements
        name_white_elem_file = "WHITE_ELEM_FOLDER/elem_white" + "_" + str(self.iter) + ".txt"
        file_white = open(name_white_elem_file, "w")
        for elem in self.main_model_part.Elements :
            if elem.IsNot(ACTIVE) :
                file_white.write(str(elem.Id))
                file_white.write('\n')
        file_white.close()

        ## Skin_node e Skin_elem per i plot
        name_skin_node_file = "SKIN_NODE_FOLDER/skin_node" + "_" + str(self.iter) + ".txt"
        file_skin_node = open(name_skin_node_file, "w")
        for node in self.current_skin_model_part.Nodes :
            file_skin_node.write(str(node.Id))
            file_skin_node.write('  ')
            file_skin_node.write(str(node.X))
            file_skin_node.write('  ')
            file_skin_node.write(str(node.Y))
            file_skin_node.write('  ')
            file_skin_node.write(str(node.Z))
            file_skin_node.write('\n')
        file_skin_node.close()
        
        ## Skin_node e Skin_elem per i plot
        name_skin_elem_file = "SKIN_ELEM_FOLDER/skin_elem" + "_" + str(self.iter) + ".txt"
        file_skin_elem = open(name_skin_elem_file, "w")
        for cond in self.current_skin_model_part.Conditions :
            file_skin_elem.write(str(cond.Id))
            file_skin_elem.write('  ')
            file_skin_elem.write(str(0))
            file_skin_elem.write('  ')
            file_skin_elem.write(str(cond.GetNodes()[0].Id))
            file_skin_elem.write('  ')
            file_skin_elem.write(str(cond.GetNodes()[1].Id))
            file_skin_elem.write('\n')
        file_skin_elem.close()

        # Clear
        self._convection_diffusion_solution_strategy.Clear()
        

    def Finalize(self):
        super().Finalize()

        # # main_model_part = self.GetComputingModelPart()
        main_model_part = self.main_model_part
        
        # # Check the error if the exact solution is known
        # self.errorL2, self.errorH1, self.max_err = self.Compute_error_transient(self.main_model_part, self.sub_model_part_fluid, self.time)

        # print("L_inf err = ", self.max_err)


        # ## Sebastian 2.0 --> save TEMPERATURE
        # temperature = []
        # for node in self.main_model_part.Nodes:
        #     temperature.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        # velocity_numpy = np.array(temperature)
        # print(velocity_numpy.shape)
        # np.save("SBtemperature.npy", velocity_numpy)






    def CheckIfMPCsAreAppliedCorrectly(self):
        main_model_part = self.GetComputingModelPart()
        T = KratosMultiphysics.Matrix(0,0)  # T Matrix
        c = KratosMultiphysics.Vector(0)    # Contant vector
        for constraint in main_model_part.MasterSlaveConstraints:
            constraint.CalculateLocalSystem(T,c,self.main_model_part.ProcessInfo)
            master_dofs_vector = constraint.GetMasterDofsVector()
            counter = 0
            master_solution = 0
            for master_dof in master_dofs_vector:
                master_dof_vel_x = master_dof.GetSolutionStepValue()
                master_solution += master_dof_vel_x*T[0,counter]
                counter += 1
            master_solution = master_solution + c[0]
            slave_dof = constraint.GetSlaveDofsVector()[0]
            slave_dof_solution = slave_dof.GetSolutionStepValue()
            try:
                relative_error = 100*abs(master_solution-slave_dof_solution)/abs(slave_dof_solution)
                if relative_error>1e-5:
                    print("----------------")
                    print(slave_dof.Id())
                    print(T)
                    print('Relative error : ', relative_error)
                    print('master solution : ', master_solution)
                    print('slave solution : ', slave_dof_solution)
            except:
                continue

