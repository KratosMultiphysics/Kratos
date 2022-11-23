from inspect import Parameter
import KratosMultiphysics
import math
import numpy as np
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import applications
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().IsDistributed():
    import KratosMultiphysics.mpi as KratosMPI
    import KratosMultiphysics.MetisApplication as KratosMetis
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver
# from KratosMultiphysics.ConvectionDiffusionApplication.ConvectionDiffusionSolver import convection_diffusion_stationary_solver
from KratosMultiphysics.ConvectionDiffusionApplication.Functions import*

from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

# Need for DistanceModificationProcess
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def CreateSolver(main_model_part, custom_settings):
    return SBMConvectionDiffusionStationarySolver(main_model_part, custom_settings)


class SBMConvectionDiffusionStationarySolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    print('ci siamo')
    skin_model_part = Import_Structural_model_part('Structural')

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
        iter = 0
        main_model_part = self.GetComputingModelPart()
        # main_model_part, skin_model_part = Import_Structural_SUB_model_part(main_model_part, 'Structural')
        skin_model_part = self.skin_model_part
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


        self.projection_surr_nodes = projection_surr_nodes
        self.closest_element = closest_element
        self.tot_sur_nodes = tot_sur_nodes
        # Set the BC at the skin mesh________________________________________________________________________________________________
        for node in skin_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, node.X + node.Y)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 0.25*(9-node.X**2-node.Y**2-2*math.log(3) + math.log(node.X**2+node.Y**2)) + 0.25 *math.sin(node.X) * math.sinh(node.Y))
            node.Fix(KratosMultiphysics.TEMPERATURE)
        
        boundary_sub_model_part = main_model_part.CreateSubModelPart("boundary_sub_model_part")
        for node in main_model_part.Nodes :
            if node.Is(BOUNDARY) :
                node.Set(BOUNDARY, False)
        for node in main_model_part.Nodes :
            if node.Is(VISITED) :
                boundary_sub_model_part.AddNode(node,0)
                node.Set(BOUNDARY, True)
                # node.Set(SLAVE,True)
        # for elem in main_model_part.Elements :
        #     if elem.Is(MARKER):
        #         boundary_sub_model_part.AddElement(elem,0)
    
        
        # Calculate the required neighbours
        nodal_neighbours_process = KratosMultiphysics.FindGlobalNodalNeighboursProcess(main_model_part)
        nodal_neighbours_process.Execute()
        elemental_neighbours_process = KratosMultiphysics.GenericFindElementalNeighboursProcess(main_model_part)
        elemental_neighbours_process.Execute()

        Parameters = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "ThermalModelPart",
            "boundary_sub_model_part_name" : "boundary_sub_model_part",
            "sbm_interface_condition_name" : "LineCondition2D2N",
            "conforming_basis" : true,
            "extension_operator_type" : "MLS",
            "mls_extension_operator_order" : 1,
            "levelset_variable_name" : "DISTANCE"
        }
        """)

        for elem in sub_model_part_fluid.Elements :
            if elem.IsNot(ACTIVE) :
                elem.Set(ACTIVE, True)
        
        # CALCOLO DEI COEFFICIENTI PER IL CALCOLO DEL GRADIENTE
        a = KratosMultiphysics.ShiftedBoundaryMeshlessInterfaceUtility(self.model,Parameters)
        # print(dir(a))
        result = a.SetSurrogateBoundaryNodalGradientWeights()
        self.result = result
        print('Number of results : ', len(result))
        print('Number of surr nodes : ', len(boundary_sub_model_part.Nodes))

        

        i = 0
        for node in boundary_sub_model_part.Nodes :
            my_result = result[node.Id]
            j = 0
            # Get the length of the matrix T_tilde 
            name = "T_tilde_" + str(i)
            globals()[name] = [[0 for _ in range(len(my_result))] for _ in range(2)]
            for key, value in my_result.items(): 
                globals()[name][0][j] = value[0]
                globals()[name][1][j] = value[1]       
                j = j + 1
            T_tilde = globals()[name]
            # Get the distance vector using the projections
            d_vector = [[0 for _ in range(1)] for _ in range(2)]
            d_vector[0] = projection_surr_nodes[i][0] - node.X
            d_vector[1] = projection_surr_nodes[i][1] - node.Y
            # Create the T matrices: 0 , 1 , ... , len(boundary_sub_model_part.Nodes)-1
            nameT = "T_" + str(i)
            # MATRIX MULTIPLICATION
            globals()[nameT] = np.matmul( np.transpose( T_tilde) , d_vector )
            i = i + 1
        print('Created the %d matrices T = [num_neighs x 1]' % len(boundary_sub_model_part.Nodes))
        # # Set the VELOCITY VECTOR FIELD_______________________________________________________________________________________________
        # for node in main_model_part.Nodes:
        #     node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [1.0, 0.0, 0.0])
        #     node.Fix(KratosMultiphysics.VELOCITY_X)
        #     node.Fix(KratosMultiphysics.VELOCITY_Y)
        #     node.Fix(KratosMultiphysics.VELOCITY_Z)


        # Sebastian -> Create CreateNewMasterSlaveConstraint
        j = 1
        for node in self.surrogate_sub_model_part.Nodes :
            my_result = self.result[node.Id]
            DofMasterVector = []
            CoeffVector = KratosMultiphysics.Vector(len(my_result)-1)
            ConstantVector = 0.0*KratosMultiphysics.Vector(len(my_result)-1)
            # Interpolate the value at the projection 
            dirichlet_projection = Interpolation(self.skin_model_part,self.closest_element,self.projection_surr_nodes, j-1, node)
            # print(dirichlet_projection)
            name = "T_" + str(j-1)
            T = globals()[name]
            # print(T)
            i = 0
            k = 0
            for key, value in my_result.items() :
                # Need to find the "node" term and bring it to the left-hand-side
                if node.Id != key :
                    node_master = main_model_part.GetNode(key)
                    DofMasterVector.append(node_master.GetDof(KratosMultiphysics.TEMPERATURE))
                    ConstantVector[i] = dirichlet_projection
                    CoeffVector[i] = - T[k]
                    i = i + 1
                    k = k + 1
                else :
                    Coeff_Slave = - T[k]
                    k = k + 1
            CoeffVector[:] = CoeffVector[:] / (1-Coeff_Slave)
            ConstantVector[:] = ConstantVector[:] / (1-Coeff_Slave)
            if j == 14 :
                self.DofMasterVector = DofMasterVector
                self.dirichlet_projection = dirichlet_projection
                self.CoeffVector = CoeffVector
                self.ConstantVector = ConstantVector
            CoeffMatrix = KratosMultiphysics.Matrix(np.array(CoeffVector).reshape(1,-1))
            DofSlaveVector = [node.GetDof(KratosMultiphysics.TEMPERATURE)]
            # Create the constraint
            if j == 2 or j == 4 :
                main_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", j, DofMasterVector, DofSlaveVector, CoeffMatrix, ConstantVector)
                # print(CoeffVector)
            j = j + 1



        # CHECK
        # Discretize if an element is boundary or not: an element is boundary if one of its node is a surrogate boundary node
        # new_elem_name = "LaplacianElement2D3N"
        # new_cond_name = "ThermalFace2D2N"
        # ## Set the element and condition names in the Json parameters
        # self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        # self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        # self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)
        # ## Call the replace elements and conditions process
        # KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
        element_name = "SBMLaplacianElement2D3N"
        # element_name = "WeaklyCompressibleNavierStokes2D3N"
        count = 0
        for elem in main_model_part.Elements :
            count = count + 1
            if elem.Is(BOUNDARY) :
                # Remove the element
                main_model_part.RemoveElement(elem)
                print(' \n elemento : ', elem.Id)
                # Initialize the list containing the nodes involved
                list_nodes_involved = []
                count_surr = 0
                # Let's count how many surrogate nodes the element has
                list_surr_nodes = []
                for node in elem.GetNodes() :
                    # list_nodes_involved.append(node.Id)
                    if node.Is(VISITED) :
                        list_surr_nodes.append(node.Id)
                        count_surr = count_surr + 1
                print('Number of surrogate nodes : ',count_surr)
                my_result = self.result[list_surr_nodes[0]]
                for key, value in my_result.items() :
                    list_nodes_involved.append(key)
                if count_surr > 1 :  # two or three nodes are surrogate nodes
                    # Need to add some additional nodes
                    my_result = self.result[list_surr_nodes[1]]
                    for key, value in my_result.items() :
                        different = 0
                        for i in range(len(list_nodes_involved)) : 
                            if key != list_nodes_involved[i] :
                                different = different + 1
                            else :
                                break
                        if different == len(list_nodes_involved) :
                            list_nodes_involved.append(key)
                    if count_surr == 3 :
                        print('Warning!! --> There are elements with 3 nodes that are surrogate nodes')
                        exit()
                # Create a new element
                main_model_part.CreateNewElement(  element_name, elem.Id, [elem.GetNodes()[0].Id, \
                    elem.GetNodes()[1].Id, elem.GetNodes()[2].Id], surrogate_sub_model_part.GetProperties()[1])
                print(list_nodes_involved)
        # exit()

    


    def InitializeSolutionStep(self):
        
        main_model_part = self.GetComputingModelPart()
        # Compute the gradient with the function ComputeNodalGradientProcess
         # Set a fake scalar field -> VELOCITY_X
        # Compute the gradient using the sub_model_part_fluid
        # KratosMultiphysics.ComputeNodalGradientProcess(
        # self.sub_model_part_fluid,
        # KratosMultiphysics.TEMPERATURE,
        # KratosMultiphysics.TEMPERATURE_GRADIENT,
        # KratosMultiphysics.NODAL_AREA).Execute()

        # Compute the Dirichlet BC at the surrogate nodes
        # surr_BC = Dirichlet_BC (main_model_part,self.skin_model_part,self.tot_sur_nodes,self.closest_element,self.projection_surr_nodes)
        
        
        
        super().InitializeSolutionStep()
        


        # Create the embedded BC
        embedded_BC = [0 for _ in range(len(self.surrogate_sub_model_part.Nodes))]
        i = 0
        for node in self.surrogate_sub_model_part.Nodes:
            # get the neighbourhoods nodes of node --> key
            my_result = self.result[node.Id]
            j = 0
            scalar_product = 0
            gradient2 = [0, 0]
            for key, value in my_result.items() :
                contributing_node = self.sub_model_part_fluid.GetNode(key)
                name = "T_" + str(i)
                T = globals()[name]
                scalar_product = scalar_product + contributing_node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) * T[j]
                j = j +1 
                # Check the gradient
                gradient2[0] = gradient2[0] + value[0] * contributing_node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
                gradient2[1] = gradient2[1] + value[1] * contributing_node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            # Get the Dirichlet value at the projection point  
            velocity = Interpolation(self.skin_model_part,self.closest_element,self.projection_surr_nodes, i, node)
            # if node.Id == 1457 :
            #     print()
            #     print(gradient2[0]- node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[0])
            #     print(gradient2[1]- node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT)[1])
            # Get the BC at the surrogate node
            embedded_BC[i] = velocity-scalar_product
            # embedded_BC[i] = surr_BC[i]
            i = i+1
        
        # Impose the embedded BC
        i = 0 
        for node in self.surrogate_sub_model_part.Nodes:
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE,embedded_BC[i])
            # node.Fix(KratosMultiphysics.TEMPERATURE)
            i = i + 1
        self.embedded_BC = embedded_BC
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        main_model_part = self.GetComputingModelPart()
        self.CheckIfMPCsAreAppliedCorrectly()
        # exit()
        # for constraint in main_model_part.MasterSlaveConstraints:
        #     # print(constraint)
        #     # print(self.CoeffVector)
        #     if constraint.Id== 14:
        #         master_dofs = constraint.GetMasterDofsVector()
        #         sum_of_master_sol = 0
        #         i = 0
        #         for dof in master_dofs:
        #             sum_of_master_sol += dof.GetSolutionStepValue() * self.CoeffVector[i]
        #             i += 1
        #         print()
        #         print("Master's dof solution: ",sum_of_master_sol + self.ConstantVector[0])
        #         slave_dof = constraint.GetSlaveDofsVector()
        #         for dof in slave_dof:
        #             slave_sol = dof.GetSolutionStepValue()
        #         print("Slave dof solution: ",slave_sol)
        #         exit()
        # i = 0 
        # for node in self.surrogate_sub_model_part.Nodes:
        #     if node.Id == 1457 :
        #         print(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE), self.embedded_BC[i])
        #         i = i + 1
        # # exit()
        # # Compute the gradient using the sub_model_part_fluid
        # KratosMultiphysics.ComputeNodalGradientProcess(
        # main_model_part,
        # KratosMultiphysics.TEMPERATURE,
        # KratosMultiphysics.TEMPERATURE_GRADIENT,
        # KratosMultiphysics.NODAL_AREA).Execute()
        # # for node in self.surrogate_sub_model_part.Nodes :
        # #     print(node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE_GRADIENT))


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
        
        # Check the error if the exact solution is known
        Compute_error(main_model_part)



    
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
                if relative_error>1e-4:
                    print("----------------")
                    print(slave_dof.Id())
                    print(T)
                    print('Relative error : ', relative_error)
                    print('master solution : ', master_solution)
                    print('slave solution : ', slave_dof_solution)
            except:
                continue

















































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
