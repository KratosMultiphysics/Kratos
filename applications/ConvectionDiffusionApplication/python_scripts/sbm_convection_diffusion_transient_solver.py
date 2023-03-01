from inspect import Parameter
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
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_transient_solver
from KratosMultiphysics.ConvectionDiffusionApplication.Functions import*

from KratosMultiphysics.assign_scalar_variable_process import AssignScalarVariableProcess

# Need for DistanceModificationProcess
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def CreateSolver(main_model_part, custom_settings):
    return SBMConvectionDiffusionTransientSolver(main_model_part, custom_settings)


class SBMConvectionDiffusionTransientSolver(convection_diffusion_transient_solver.ConvectionDiffusionTransientSolver):
    print('ci siamo transient solver')
    skin_model_part = Import_Structural_model_part('Structural_circle1x1_huge')
    # skin_model_part = Import_Structural_model_part('Structural_rectangle0.1')
    # skin_model_part = Import_Structural_model_part('Structural_storto_huge')


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
        
        self.time_step = 0.05
        self.time = 0.0

        iter = 0
        main_model_part = self.GetComputingModelPart()
        for node in main_model_part.Nodes :
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY,0, [1.0, 1.0, 0.0])
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
        self.closest_element = a.FindClosestElement(main_model_part, skin_model_part)

        self.surrogate_sub_model_part = main_model_part.CreateSubModelPart("surrogate_sub_model_part")
        tot_sur_nodes = 0 
        for node in main_model_part.Nodes :
            if node.Is(BOUNDARY):
                self.surrogate_sub_model_part.AddNode(node,0)
                tot_sur_nodes = tot_sur_nodes + 1

        # Find the projection onto the skin elements for each surr node
        projection_surr_nodes = Find_projections(main_model_part,skin_model_part,tot_sur_nodes,self.closest_element)
        # Then we create a sub_model part with just the elements & the nodes "outside" the surrogate boundary
        sub_model_part_fluid = Create_sub_model_part_fluid(main_model_part,iter)
        self.sub_model_part_fluid = sub_model_part_fluid
        print('Creato il sub_model_part per il calcolo del gradiente')
        self.projection_surr_nodes = projection_surr_nodes
        self.tot_sur_nodes = tot_sur_nodes


        # Set the BC at the skin mesh________________________________________________________________________________________________
        for node in skin_model_part.Nodes:
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, (self.time)**2 * math.sin(2*node.X) * math.cos(2*node.Y))
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, (self.time)**2 * math.sin(8*node.X) * math.cos(8*node.Y))
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, (self.time)**2 * math.sin(2*node.X) * math.cos(2*node.Y))
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
            if node.X != 1.0 and node.X != 0.0:
                Impose_MPC_Globally (main_model_part, self.result, self.skin_model_part, self.closest_element, self.projection_surr_nodes, T, node, j)
            j = j + 1
        
    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.time = self.time + self.time_step

        ## pdb.set_trace()
        ## modify the skin_model_part exact solution (when the exact solution depends on time) 
        for node in self.skin_model_part.Nodes:
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, (self.time)**2* (node.Y-0.1) * math.sin(8*node.X) * math.cos(8*node.Y))
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, (self.time)**2 *(node.Y-0.1)* math.sin(2*node.X) * math.cos(2*node.Y))
            node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, (self.time)**2 * math.sin(8*node.X) * math.cos(8*node.Y))
            node.Fix(KratosMultiphysics.TEMPERATURE)

        # Remove MPC and regreate them
        for constraint in self.main_model_part.MasterSlaveConstraints :
            constraint.Set(KratosMultiphysics.TO_ERASE)
        self.main_model_part.RemoveMasterSlaveConstraintsFromAllLevels(KratosMultiphysics.TO_ERASE)

        # Create the CreateMasterSlaveConstraints
        j = 1
        for node in self.surrogate_sub_model_part.Nodes :
            name = "T_" + str(j-1)
            T = globals()[name]
            if node.X != 1.0 and node.X != 0.0:
                Impose_MPC_Globally (self.main_model_part, self.result, self.skin_model_part, self.closest_element, self.projection_surr_nodes, T, node, j)
            # node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, (self.time)**2*(node.Y-1/10)* math.sin(8*node.X) * math.cos(8*node.Y))
            # node.Fix(KratosMultiphysics.TEMPERATURE)
            j = j + 1

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        ## Sebastian check for MPC
        self.CheckIfMPCsAreAppliedCorrectly()


        

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

        # # main_model_part = self.GetComputingModelPart()
        main_model_part = self.main_model_part
        
        # Check the error if the exact solution is known
        self.errorL2, self.errorH1, self.max_err = self.Compute_error_transient(self.main_model_part, self.sub_model_part_fluid, self.time)

        print("L_inf err = ", self.max_err)

        file_white = open("elem_white.txt", "w")
        for elem in self.main_model_part.Elements :
            if elem.IsNot(ACTIVE) :
                file_white.write(str(elem.Id))
                file_white.write('\n')
        file_white.close()


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
                if relative_error>1e-8:
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
            exact = time**2 * math.sin(8*node.X) * math.cos(8*node.Y)
            exact_grad[0] = 2* time**2 *math.cos(2*node.X) * math.cos(2*node.Y)
            exact_grad[1] = - 2* time**2 * math.sin(2*node.X) * math.sin(2*node.Y)
            ## (1-0.y)
            # exact = (time)**2 *(node.Y-0.1)* math.sin(8*node.X) * math.cos(8*node.Y)
            # exact_grad[0] = 2* time**2 *(node.Y-0.1)*math.cos(2*node.X) * math.cos(2*node.Y)
            # exact_grad[1] = time**2 * math.sin(2*node.X) * (math.cos(2*node.Y)-2*(1-node.Y)*math.sin(2*node.Y))
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