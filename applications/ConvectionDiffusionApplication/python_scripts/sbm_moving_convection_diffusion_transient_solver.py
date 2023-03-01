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
    return SBMMovingConvectionDiffusionTransientSolver(main_model_part, custom_settings)


class SBMMovingConvectionDiffusionTransientSolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    print('ci siamo transient solver')
    # skin_model_part = Import_Structural_model_part('Structural_circle1x1_huge')
    # skin_model_part = Import_Structural_model_part('Structural_rectangle0.1')
    skin_model_part = Import_Structural_model_part('Structural_rectangle_vertical')


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
        
        self.time_step = 0.00005
        self.time = 0.0
        self.object_velocity = [0.0, 8.00, 0.0]

        self.iter = 0
        main_model_part = self.GetComputingModelPart()
        for node in main_model_part.Nodes :
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [0.0, 0.0, 0.0])
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
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
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
            if node.X != 0.4 and node.X != 0.6:
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
        for node in self.skin_model_part.Nodes :
            ## Translation **
            current_skin_model_part.CreateNewNode(node.Id, node.X + self.time * self.object_velocity[0], node.Y + self.time * self.object_velocity[1] , 0.0)
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
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0+self.time*0)
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
            if node.X != 0.4 and node.X != 0.6 :
                Impose_MPC_Globally (main_model_part, result, current_skin_model_part, closest_element, projection_surr_nodes, T, node, j)
            j = j + 1
        
        if math.fmod(self.iter,100) == 0 :
            name_surr_file = "Surr_B" + "_" + str(self.iter) + ".txt"
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
        
        # Initialize solution step (Alla fine!)
        super().InitializeSolutionStep()
        

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        
        self._convection_diffusion_solution_strategy.Clear()
        

        ## Sebastian check for MPC
        self.CheckIfMPCsAreAppliedCorrectly()


        

    def Finalize(self):
        super().Finalize()

        # # main_model_part = self.GetComputingModelPart()
        main_model_part = self.main_model_part
        
        # Check the error if the exact solution is known
        self.errorL2, self.errorH1, self.max_err = self.Compute_error_transient(self.main_model_part, self.sub_model_part_fluid, self.time)

        print("L_inf err = ", self.max_err)


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


    def Compute_error_transient(self, main_model_part, sub_model_part_fluid, time) :
        file_due = open("error.txt", "w")
        L2_err = 0
        L2_err_area = 0
        L2_grad_err = 0
        L2_grad_err_area = 0
        KratosMultiphysics.ComputeNodalGradientProcess(
        sub_model_part_fluid,
        KratosMultiphysics.TEMPERATURE,
        KratosMultiphysics.TEMPERATURE_GRADIENT,
        KratosMultiphysics.NODAL_AREA).Execute()
        total_number_fluid_nodes = 0
        total_area = 0
        max_err = 0
        for node in main_model_part.Nodes :
            exact_grad = KratosMultiphysics.Array3()
            ## sin(x)*cos(y)
            # exact = time**2 * math.sin(2*node.X) * math.cos(2*node.Y)
            # exact_grad[0] = 2* time**2 *math.cos(2*node.X) * math.cos(2*node.Y)
            # exact_grad[1] = - 2* time**2 * math.sin(2*node.X) * math.sin(2*node.Y)
            ## (1-0.y)
            exact = (time)**2 *(node.Y-0.1)* math.sin(2*node.X) * math.cos(2*node.Y)
            exact_grad[0] = 2* time**2 *(node.Y-0.1)*math.cos(2*node.X) * math.cos(2*node.Y)
            exact_grad[1] = time**2 * math.sin(2*node.X) * (math.cos(2*node.Y)-2*(1-node.Y)*math.sin(2*node.Y))
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) > 0 :
                if abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact) > max_err :
                    max_err = abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact)
                total_number_fluid_nodes = total_number_fluid_nodes + 1
                nodal_area = node.GetValue(NODAL_AREA)
                L2_err = L2_err + (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact)**2
                L2_err_area = L2_err_area + nodal_area * (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact)**2
                L2_grad_err = L2_grad_err + (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0]-exact_grad[0])**2 + (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1]-exact_grad[1])**2 
                L2_grad_err_area = L2_grad_err_area + nodal_area * ((node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0]-exact_grad[0])**2 + (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1]-exact_grad[1])**2 )
                total_area = total_area + nodal_area
                file_due.write(str(abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)-exact)))
                file_due.write('       ')
                file_due.write(str(abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0]-exact_grad[0])))
                file_due.write('       ')
                file_due.write(str(abs(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1]-exact_grad[1])))
                file_due.write('       ')
                file_due.write(str(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0]))
                file_due.write('       ')
                file_due.write(str(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1]))
            else :
                file_due.write('0.0')
                file_due.write('     ')
                file_due.write('0.0')
                file_due.write('     ')
                file_due.write('0.0')
                file_due.write('     ')
                file_due.write('0.0')
                file_due.write('     ')
                file_due.write('0.0')
            file_due.write('\n')
            
        if total_number_fluid_nodes == 0 :
            total_number_fluid_nodes = len(sub_model_part_fluid.Nodes)
        L2_err = math.sqrt(L2_err  / total_number_fluid_nodes )
        L2_err_area = math.sqrt(L2_err_area / total_area) 
        H1_err_area = L2_err_area + math.sqrt(L2_grad_err_area / total_area ) 
        # print('Errore in norma L2 (equal areas): ', L2_err)
        print('Errore in norma L2 : ', L2_err_area)
        print('Errore in norma H1 : ', H1_err_area)
        file_due.close
        return L2_err_area, H1_err_area, max_err















































































































































































# import time
# import KratosMultiphysics
# import numpy as np
# import math

# # Import applications
# import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
# if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
#     import KratosMultiphysics.mpi as KratosMPI
#     import KratosMultiphysics.MetisApplication as KratosMetis
#     import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# # Import base class file
# from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_transient_solver
# from KratosMultiphysics.ConvectionDiffusionApplication.Functions import*

# # from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

# # Need for DistanceModificationProcess
# import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# def CreateSolver(main_model_part, custom_settings):
#     return SBMMovingConvectionDiffusionTransientSolver(main_model_part, custom_settings)


# class SBMMovingConvectionDiffusionTransientSolver(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver):
#     print('Ci siamo, moving object solver ready!')
#     skin_model_part = Import_Structural_model_part('Structural_circle1x1_huge')

#     def __init__(self, main_model_part, custom_settings):
#         super().__init__(main_model_part, custom_settings)
#         # * Here I don't have access to the model_part yet
#         KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

#     def AddVariables(self):
#         super().AddVariables()
#         self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
#         self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

#     def Initialize(self):
#         super().Initialize()
#         self.current_time = 0
#         self.delta_time = 0.01
#         self.iter = 1

#         main_model_part = self.GetComputingModelPart()
        
#         # Set the VELOCITY VECTOR FIELD_______________________________________________________________________________________________
#         for node in main_model_part.Nodes:
#             #node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [0.0, 0.0, 0.0])
#             node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [0.0, 0.0, 0.0])
#             #node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [-1/8*node.Y**2 + 5/8 *node.Y -9/32, 0.0, 0.0])
#             node.Fix(KratosMultiphysics.VELOCITY_X)
#             node.Fix(KratosMultiphysics.VELOCITY_Y)
#             node.Fix(KratosMultiphysics.VELOCITY_Z)
#             if node.Is(BOUNDARY) : 
#                 node.Set(BOUNDARY, False)

#         # self.object_velocity = [0.2, 0.0, 0.0]
#         # self.object_velocity = [0.0, 0.0, 0.0]

#         # Total number of skin elements
#         self.tot_skin_el = len(self.skin_model_part.Conditions)
#         print('Number of skin elements: ', self.tot_skin_el)
        

#     def InitializeSolutionStep(self):

#         super().InitializeSolutionStep()
#         # self.current_time = self.current_time + self.delta_time

#         main_model_part = self.GetComputingModelPart()
#         # for node in main_model_part.Nodes :
#         #     if  node.IsNot(ACTIVE):
#         #         node.Set(ACTIVE, True)
#         #     if  node.Is(BOUNDARY):
#         #         node.Set(BOUNDARY, False)
#         #     if  node.Is(SLAVE):
#         #         node.Set(SLAVE, False)
#         #     if node.Is(INTERFACE) :
#         #         node.Set(INTERFACE, False)
#         # for elem in main_model_part.Elements :
#         #     if elem.IsNot(ACTIVE) :
#         #         elem.Set(ACTIVE,True)
#         #     if elem.Is(BOUNDARY):
#         #         elem.Set(BOUNDARY,False)
#         #     if elem.Is(INTERFACE):
#         #         elem.Set(INTERFACE,False)
#         #     if elem.Is(MARKER):
#         #         elem.Set(MARKER,False)
#         #     if elem.Is(VISITED) :
#         #         elem.Set(VISITED,False)


#         # current_model = KratosMultiphysics.Model()
#         # current_skin_model_part = current_model.CreateModelPart("CurrentSkinModelPart")
#         # current_skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

#         # for node in self.skin_model_part.Nodes :
#         #     ## Translation **
#         #     current_skin_model_part.CreateNewNode(node.Id, node.X + self.current_time * self.object_velocity[0], node.Y + self.current_time * self.object_velocity[1] , 0.0)
#         #     # # Rotation **
#         #     # r = math.sqrt((node.X-2.1)**2 + (node.Y-1.25)**2)
#         #     # if node.X != 2.1 :
#         #     #     if node.X > 2.1: 
#         #     #         teta = math.atan((node.Y-1.25)/(node.X-2.1))
#         #     #     else :
#         #     #         teta = math.atan((node.Y-1.25)/(node.X-2.1)) + math.pi
#         #     # else :
#         #     #     teta = math.pi/2
#         #     # omega = 10
#         #     # current_skin_model_part.CreateNewNode(node.Id, 2.1 + r * math.cos(teta+omega*self.current_time), 1.25 + r * math.sin(teta+omega*self.current_time) , 0.0)
#         # for cond in self.skin_model_part.Conditions :
#         #     # node1 = cond.GetNodes()[0].Id
#         #     property = self.skin_model_part.GetProperties()[0]
#         #     current_skin_model_part.CreateNewCondition("LineCondition2D2N", cond.Id, [cond.GetNodes()[0].Id, cond.GetNodes()[1].Id],property)


#         current_skin_model_part = self.skin_model_part

#         start_time = time.time()
#         KratosMultiphysics.CalculateDistanceToSkinProcess2D(main_model_part, current_skin_model_part).Execute()
#         print("--> %s seconds for CalculateDistanceToSkinProcess2D" % (time.time() - start_time))
#         # Find the surrogate boundary nodes
#         a = KratosMultiphysics.FindSurrogateNodesProcess2D(main_model_part, current_skin_model_part)
#         a.Execute()
#         closest_element = a.FindClosestElement(main_model_part, current_skin_model_part)

#         name_surrogate_sub_model_part = "surrogate_sub_model_part" + "_" + str(self.iter)
#         surrogate_sub_model_part = main_model_part.CreateSubModelPart(name_surrogate_sub_model_part)
#         tot_sur_nodes = 0 
#         for node in main_model_part.Nodes :
#             if node.Is(BOUNDARY):
#                 surrogate_sub_model_part.AddNode(node,0)
#                 tot_sur_nodes = tot_sur_nodes + 1

#         projection_surr_nodes = Find_projections(main_model_part,current_skin_model_part,tot_sur_nodes,closest_element)
#         sub_model_part_fluid = Create_sub_model_part_fluid(main_model_part,self.iter)

#         self.sub_model_part_fluid = sub_model_part_fluid
#         print('Creato il sub_model_part per il calcolo del gradiente')
#         self.tot_sur_nodes = tot_sur_nodes

#         # Set the BC at the skin mesh________________________________________________________________________________________________
#         for node in current_skin_model_part.Nodes:
#             node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
#             node.Fix(KratosMultiphysics.TEMPERATURE)

#         # Calculate the required neighbours
#         nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(main_model_part)
#         nodal_neighbours_process.Execute()

#         ## Find the so called "INTERFACE" elements
#         for elem in sub_model_part_fluid.Elements :
#             if elem.Is(BOUNDARY) :
#                 count_surr = 0
#                 # Let's count how many surrogate nodes the element has
#                 for node in elem.GetNodes() :
#                     if node.Is(BOUNDARY) :
#                         count_surr = count_surr + 1
#                 if count_surr > 1 :  # two or three nodes are surrogate nodes
#                     elem.Set(INTERFACE, True) # FUNDAMENTAL

#         elemental_neighbours_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(main_model_part)
#         elemental_neighbours_process.Execute()

#         ## Compute the gradient coefficients for each of the surrogate node
#         result = ComputeGradientCoefficients (sub_model_part_fluid, self.model, surrogate_sub_model_part)
        
#         ## Compute the T matrix for imposition of sbm condition
#         i = 0
#         for node in surrogate_sub_model_part.Nodes :
#             # Create the T matrices: 0 , 1 , ... , len(boundary_sub_model_part.Nodes)-1
#             nameT = "T_" + str(i)
#             globals()[nameT] = Compute_T_matrix (result, node, projection_surr_nodes, i)
#             i = i + 1

#         ## Remove MPC
#         for constraint in self.main_model_part.MasterSlaveConstraints :
#             constraint.Set(KratosMultiphysics.TO_ERASE)
#         self.main_model_part.RemoveMasterSlaveConstraintsFromAllLevels(KratosMultiphysics.TO_ERASE)

#         # Create the CreateMasterSlaveConstraints
#         j = 1
#         for node in surrogate_sub_model_part.Nodes :
#             name = "T_" + str(j-1)
#             T = globals()[name]
#             Impose_MPC_Globally (main_model_part, result, current_skin_model_part, closest_element, projection_surr_nodes, T, node, j)
#             j = j + 1

#         name_surr_file = "Surr_B" + "_" + str(self.iter) + ".txt"
#         # This is compiles after solving all the solution steps
#         file_tre = open(name_surr_file, "w")
#         for node in surrogate_sub_model_part.Nodes :
#             file_tre.write(str(node.X))
#             file_tre.write('  ')
#             file_tre.write(str(node.Y))
#             file_tre.write('\n')
#         file_tre.close()
#         self.iter = self.iter+1




#     def FinalizeSolutionStep(self) :
#         main_model_part = self.GetComputingModelPart()
#         # for node in main_model_part.Nodes:
#         #     node.Free(KratosMultiphysics.TEMPERATURE)
#         super().FinalizeSolutionStep()
