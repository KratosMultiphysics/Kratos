import KratosMultiphysics
import numpy as np

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver
from KratosMultiphysics.ConvectionDiffusionApplication.Functions import*

def CreateSolver(main_model_part, custom_settings):
    return SBMConvectionDiffusionStationarySolver(main_model_part, custom_settings)


class SBMConvectionDiffusionStationarySolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    print('ci siamo')
    skin_model_part = Import_Structural_model_part('Structural_COMPLEX')
    #  "value"           : "0.25*(9-x**2-y**2-2*ln(3) + ln(x**2+y**2)) +  0.25 * sin(x) * sinh(y)"
    #  "value"           : "cos(x)*cos(y) - sin(x)*sin(y) + 2*sin(x)*cos(y)",
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

        iter = 0
        main_model_part = self.GetComputingModelPart()
        for node in main_model_part.Nodes :
            # node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [5.0, 10.0, 0.0])
            # node.Fix(KratosMultiphysics.VELOCITY_X)
            # node.Fix(KratosMultiphysics.VELOCITY_Y)
            # node.Fix(KratosMultiphysics.VELOCITY_Z)
            if node.Is(KratosMultiphysics.BOUNDARY) : 
                # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.25*(9-node.X**2-node.Y**2-2*math.log(3) + math.log((node.X)**2+(node.Y)**2)) + 0.25 *math.sin(node.X) * math.sinh(node.Y))
                # node.Fix(KratosMultiphysics.TEMPERATURE)
                node.Set(KratosMultiphysics.BOUNDARY, False)
        skin_model_part = self.skin_model_part
        # Total number of skin elements
        tot_skin_el = len(self.skin_model_part.Conditions)
        print('Number of skin elements: ', tot_skin_el)

        KratosMultiphysics.CalculateDistanceToSkinProcess2D(main_model_part, skin_model_part).Execute()

        # Find the surrogate boundary nodes
        a = KratosMultiphysics.FindSurrogateNodesProcess2D(main_model_part, skin_model_part)
        a.Execute()
        self.closest_element = a.FindClosestElement(main_model_part, skin_model_part)

        self.surrogate_sub_model_part = main_model_part.CreateSubModelPart("surrogate_sub_model_part")
        tot_sur_nodes = 0 
        for node in main_model_part.Nodes :
            if node.Is(KratosMultiphysics.BOUNDARY):
                self.surrogate_sub_model_part.AddNode(node,0)
                tot_sur_nodes = tot_sur_nodes + 1
        
        #### self.surrogate_sub_model_part, tot_sur_nodes = FindSurrogateNodes(main_model_part,iter)
        #### self.closest_element = FindClosestSkinElement(main_model_part,skin_model_part,tot_sur_nodes,tot_skin_el)
        
        # Find the projection onto the skin elements for each surr node
        projection_surr_nodes = Find_projections(main_model_part,skin_model_part,tot_sur_nodes,self.closest_element)
        # Then we create a sub_model part with just the elements & the nodes "outside" the surrogate boundary
        sub_model_part_fluid = Create_sub_model_part_fluid(main_model_part,iter)
        self.sub_model_part_fluid = sub_model_part_fluid
        print('Creato il sub_model_part per il calcolo del gradiente')
        self.projection_surr_nodes = projection_surr_nodes
        # self.closest_element = closest_element
        self.tot_sur_nodes = tot_sur_nodes


        # Set the BC at the skin mesh________________________________________________________________________________________________
        for node in skin_model_part.Nodes:
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, node.X + node.Y)
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.25*(9 - ((node.X)**2 + (node.Y)**2) ) ) # --> Paraboloide
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.25*(9-node.X**2-node.Y**2-2*math.log(3) + math.log((node.X)**2+(node.Y)**2)) + 0.25 *math.sin(node.X) * math.sinh(node.Y))
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, math.sin(node.X) * math.cos(node.Y))
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, math.sin(node.X)*math.cos(node.Y)+math.log(1+node.X**2+node.Y**2))
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, (1-node.X))
            node.Fix(KratosMultiphysics.TEMPERATURE)

        # Calculate the required neighbours
        nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(main_model_part)
        nodal_neighbours_process.Execute()

        ## Find the so called "INTERFACE" elements
        for elem in sub_model_part_fluid.Elements :
            if elem.Is(KratosMultiphysics.BOUNDARY) :
                count_surr = 0
                # Let's count how many surrogate nodes the element has
                for node in elem.GetNodes() :
                    if node.Is(KratosMultiphysics.BOUNDARY) :
                        count_surr = count_surr + 1
                if count_surr > 1 :  # two or three nodes are surrogate nodes
                    elem.Set(KratosMultiphysics.INTERFACE, True) # FUNDAMENTAL

        elemental_neighbours_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(main_model_part)
        elemental_neighbours_process.Execute()


        ## Compute the gradint coefficients for each of the surrogate node
        self.result = ComputeGradientCoefficients (sub_model_part_fluid, self.model, self.surrogate_sub_model_part)
        
        ## Compute the T matrix for imposition of sbm condition
        i = 0
        for node in self.surrogate_sub_model_part.Nodes :
            # Create the T matrices: 0 , 1 , ... , len(boundary_sub_model_part.Nodes)-1
            nameT = "T_" + str(i)
            globals()[nameT] = Compute_T_matrix (self.result, node, projection_surr_nodes, i)
            i = i + 1

          
        # Create the CreateMasterSlaveConstraints
        j = 1
        for node in self.surrogate_sub_model_part.Nodes :
            name = "T_" + str(j-1)
            T = globals()[name]
            if node.X != 1.0 :
                Impose_MPC_Globally (main_model_part, self.result, self.skin_model_part, self.closest_element, self.projection_surr_nodes, T, node, j)
            j = j + 1
        
    def InitializeSolutionStep(self):

        # ## Compute the gradient with the function ComputeNodalGradientProcess using the sub_model_part_fluid
        # KratosMultiphysics.ComputeNodalGradientProcess(
        # self.sub_model_part_fluid,
        # KratosMultiphysics.TEMPERATURE,
        # KratosMultiphysics.TEMPERATURE_GRADIENT,
        # KratosMultiphysics.NODAL_AREA).Execute()

        ## Compute the Dirichlet BC at the surrogate nodes
        # surr_BC = Dirichlet_BC (main_model_part,self.skin_model_part,self.tot_sur_nodes,self.closest_element,self.projection_surr_nodes)
        
        super().InitializeSolutionStep()
        
        ## Create the embedded BC con matrix T
        # embedded_BC = [0 for _ in range(len(self.surrogate_sub_model_part.Nodes))]
        # i = 0
        # for node in self.surrogate_sub_model_part.Nodes:
        #     # get the neighbourhoods nodes of node --> key
        #     my_result = self.result[node.Id]
        #     j = 0
        #     scalar_product = 0
        #     # gradient2 = [0, 0]
        #     for key, value in my_result.items() :
        #         contributing_node = self.sub_model_part_fluid.GetNode(key)
        #         name = "T_" + str(i)
        #         T = globals()[name]
        #         scalar_product = scalar_product + contributing_node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) * T[j]
        #         j = j +1 

        #     dirichlet_projection = Interpolation(self.skin_model_part,self.closest_element,self.projection_surr_nodes, i, node)
        #     # Get the BC at the surrogate node
        #     embedded_BC[i] = dirichlet_projection - scalar_product
        #     i = i+1
        # self.embedded_BC = embedded_BC
        # ## Impose the embedded BC
        # i = 0 
        # for node in self.surrogate_sub_model_part.Nodes:
        #     # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,surr_BC[i])
        #     # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,embedded_BC[i])
        #     # node.Fix(KratosMultiphysics.TEMPERATURE)
        #     i = i + 1
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        ## Sebastian check for MPC
        self.CheckIfMPCsAreAppliedCorrectly()
        

    def Finalize(self):
        super().Finalize()
        ## This is compiles after solving all the solution steps
        file_tre = open("Surr_B.txt", "w")
        for node in self.surrogate_sub_model_part.Nodes :
            file_tre.write(str(node.X))
            file_tre.write('  ')
            file_tre.write(str(node.Y))
            file_tre.write('\n')
        file_tre.close()

        # # # main_model_part = self.GetComputingModelPart()
        # main_model_part = self.main_model_part
        
        # Check the error if the exact solution is known
        self.errorL2, self.errorH1, self.max_err = Compute_error(self.main_model_part, self.sub_model_part_fluid)
        print(self.max_err)


        ## Sebastian 2.0 --> save TEMPERATURE
        temperature = []
        for node in self.main_model_part.Nodes:
            temperature.append(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        velocity_numpy = np.array(temperature)
        print(velocity_numpy.shape)
        np.save("SBtemperature.npy", velocity_numpy)

        file_white = open("elem_white.txt", "w")
        for elem in self.main_model_part.Elements :
            if elem.IsNot(KratosMultiphysics.ACTIVE) :
                file_white.write(str(elem.Id))
                file_white.write('\n')
        file_white.close()









































































              
        # with open('Values.txt') as f:
        #     i = 0
        #     for node in main_model_part.Nodes :
        #         # print(float(f.readlines(2)[0]))
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, float(f.readlines(2)[0]))
        #         # print(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE))
        #         i = i+1


        # file_21 = open("Values.txt", "w")
        # for node in main_model_part.Nodes :
        #         file_21.write(str(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)))
        #         file_21.write(' \n ')
        # file_21.close()



   


        # # Sebastian -> Create CreateNewMasterSlaveConstraint
        # j = 1
        # for node in self.surrogate_sub_model_part.Nodes :
        #     name = "T_" + str(j-1)
        #     T = globals()[name]
        #     Impose_MPC_Globally (main_model_part, self.result, self.skin_model_part, self.closest_element, self.projection_surr_nodes, T, node, j)
        #     j = j + 1

        




        ## Sebastian check for MPC
        # self.CheckIfMPCsAreAppliedCorrectly()





        # main_model_part = self.GetComputingModelPart()
        # ## Create the embedded BC con matrix T
        # embedded_BC = [0 for _ in range(len(self.surrogate_sub_model_part.Nodes))]
        # i = 0
        # for node in self.surrogate_sub_model_part.Nodes:
        #     # get the neighbourhoods nodes of node --> key
        #     my_result = self.result[node.Id]
        #     j = 0
        #     scalar_product = 0
        #     # gradient2 = [0, 0]
        #     for key, value in my_result.items() :
        #         contributing_node = self.sub_model_part_fluid.GetNode(key)
        #         name = "T_" + str(i)
        #         T = globals()[name]
        #         scalar_product = scalar_product + contributing_node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) * T[j]
        #         j = j +1 
        #     dirichlet_projection = Interpolation(self.skin_model_part,self.closest_element,self.projection_surr_nodes, i, node)
        #     embedded_BC[i] = dirichlet_projection - scalar_product
        #     i = i+1
        # # Compute the gradient with the function ComputeNodalGradientProcess using the sub_model_part_fluid
        # KratosMultiphysics.ComputeNodalGradientProcess(
        # self.sub_model_part_fluid,
        # KratosMultiphysics.TEMPERATURE,
        # KratosMultiphysics.TEMPERATURE_GRADIENT,
        # KratosMultiphysics.NODAL_AREA).Execute()
        # ## Compute the Dirichlet BC at the surrogate nodes
        # surr_BC = Dirichlet_BC (main_model_part,self.skin_model_part,self.tot_sur_nodes,self.closest_element,self.projection_surr_nodes)
        # i = 0
        # for node in self.surrogate_sub_model_part.Nodes:
        #     exact = 0.25*(9-node.X**2-node.Y**2-2*math.log(3) + math.log(node.X**2+node.Y**2)) + 0.25 *math.sin(node.X) * math.sinh(node.Y) 
        #     # print(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE),embedded_BC[i],surr_BC[i])
        #     i = i + 1






























    
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
                if relative_error>1e-8:
                    print("----------------")
                    print(slave_dof.Id())
                    print(T)
                    print('Relative error : ', relative_error)
                    print('master solution : ', master_solution)
                    print('slave solution : ', slave_dof_solution)
            except:
                continue




















































    #### Private functions ####
    def _Prova_Function(self):
        print('Prova function')
        return








# # CHECK
        # # Discretize if an element is boundary or not: an element is boundary if one of its node is a surrogate boundary node
        # # element_name = "SBMLaplacianElement2D3N"
        # element_name = "LaplacianElement2D3N"
        # condition_name = "SBMLaplacianCondition"
        # number_of_conditions = len(main_model_part.Conditions)
        # count_id_condition = 0
        # for elem in main_model_part.Elements :
        #     if elem.Is(BOUNDARY) :
        #         count_id_condition = count_id_condition + 1
        #         ## Remove the element
        #         # main_model_part.RemoveElement(elem)
        #         # Initialize the list containing the nodes involved
        #         list_nodes_involved = []
        #         count_surr = 0
        #         # Let's count how many surrogate nodes the element has
        #         list_surr_nodes = []
        #         for node in elem.GetNodes() :
        #             if node.Is(BOUNDARY) :
        #                 list_surr_nodes.append(node.Id)
        #                 count_surr = count_surr + 1
        #         my_result = self.result[list_surr_nodes[0]]
        #         for key, value in my_result.items() :
        #             list_nodes_involved.append(key)
        #         if count_surr > 1 :  # two or three nodes are surrogate nodes
        #             # Need to add some additional nodes
        #             my_result = self.result[list_surr_nodes[1]]
        #             for key, value in my_result.items() :
        #                 different = 0
        #                 for i in range(len(list_nodes_involved)) : 
        #                     if key != list_nodes_involved[i] :
        #                         different = different + 1
        #                     else :
        #                         break
        #                 if different == len(list_nodes_involved) :
        #                     list_nodes_involved.append(key)
        #             if count_surr == 3 :
        #                 print('Warning!! --> There are elements with 3 nodes that are surrogate nodes')
        #                 exit()
        #         # Create a new element
        #         # main_model_part.CreateNewElement(  element_name, elem.Id, [elem.GetNodes()[0].Id, \
        #         #     elem.GetNodes()[1].Id, elem.GetNodes()[2].Id], surrogate_sub_model_part.GetProperties()[1])
        #         # main_model_part.CreateNewElement(  element_name, elem.Id, list_nodes_involved, surrogate_sub_model_part.GetProperties()[1])
        #         main_model_part.CreateNewCondition(  condition_name, number_of_conditions + count_id_condition, list_nodes_involved, surrogate_sub_model_part.GetProperties()[1])
        #         # print(list_nodes_involved)



        # pdb.set_trace()