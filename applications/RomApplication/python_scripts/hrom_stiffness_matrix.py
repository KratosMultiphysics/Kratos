import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from math import cos, sin 
import json

def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return HROMStiffnessMatrix(model, parameters["Parameters"])

## All the processes python should be derived from "Process"
class HROMStiffnessMatrix(KratosMultiphysics.Process):
    """
    This process allows the user to obtain the Stiffness Matrices of two Master nodes related to two specific Slave surfaces.\n
    Example:\n
    Master = MasterStiffnessMatrixProcess(model, parameters)\n
    Master.Run()\n
    Public member variables:
    model -- the container of the different model parts.
    parameters -- Kratos parameters containing solver parameters.
    """

    def __init__(self, model, parameters):
        """ The default constructor of the class
        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing solver parameters.
        """

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""{
                "save_matrix": true,
                "model_part_list": ["Surface_name_1","Surface_name_2"],
                "list_of_master_coordinates": [[0.0,0.0,0.0],[1.0,0.0,0.0]],
                "eps_perturbation": 1e-5
            }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.parameters = parameters
        self.model_part = model["Structure"]
        self.eps_perturbation = self.parameters["eps_perturbation"].GetDouble()
        self.Value()
        self.inf_rot = 1e-5 #Infinitesimal rotation 
        self.master_coor = self.parameters["list_of_master_coordinates"].GetMatrix()
        self.slave_surface_name = self.parameters["model_part_list"].GetStringArray()
        self.number_master_nodes = len(self.slave_surface_name)
        self.dim = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        self.save_matrix = self.parameters["save_matrix"].GetBool()


    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.ChangeVectorValues()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.MasterStiffnessVector()

    def Value(self): #Still need to make cheaper and decide the final criterion.
        '''
        Gives an infinitesimal displacement for the 3D model part
        '''
        minimum = 1e30
        maximum = -1e30
        for node in self.model_part.Nodes:
            if min(node.X,node.Y,node.Z) < minimum:
                minimum = min(node.X,node.Y,node.Z)
            if max(node.X,node.Y,node.Z) > maximum:
                maximum = max(node.X,node.Y,node.Z)
        self.value = self.eps_perturbation*(maximum - minimum) 

    def ChangeVectorValues(self):
        """
        Similar to AddVector/ScalarValueProcess, specific for MasterStiffnessMatrixProcess.\n
        Produces an infinitesimal displacement/rotation on each direction.\n
        It is needed to provide "DoF", which stands for the loop on DoFs and n to recognize which surface Stiffness Matrix to obtain.
        """
        current_step = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)-1
        # self.n = 0 obtaining the stiffness matrix for the first slave surface; self.n = 1 obtaining the stiffness matrix for the second slave surface.
        self.n = int(current_step/(2*self.dim))

        self.DoF = current_step-(2*self.dim*self.n) #Local surface DoF

        slave_submodel_part = self.model_part.GetSubModelPart(self.slave_surface_name[self.n])#Slave surface.
        
        if self.DoF<3: # Displacement DoF(0-2)
            #Setting the infinitesimal displacement/rotation boundary condition 
            for node in slave_submodel_part.Nodes: #Slave surface
                disp = KratosMultiphysics.Vector(3,0)
                disp[self.DoF] = self.value
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,disp)
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)
        else:# Rotation DoF(3-5)
            R = self.ComputeRotationMatrix()
            v = KratosMultiphysics.Vector(3)
            for node in slave_submodel_part.Nodes: #Slave surface
                v[0] = node.X0 - self.master_coor[self.n,0]
                v[1] = node.Y0 - self.master_coor[self.n,1]
                v[2] = node.Z0 - self.master_coor[self.n,2]
                disp = R*v-v
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,disp)
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        for i in [x for x in range(self.number_master_nodes) if x != self.n]:
            fixed_submodel_part = self.model_part.GetSubModelPart(self.slave_surface_name[i])#Fixed surface.
            #Setting the infinitesimal fixed boundary condition 
            for node in fixed_submodel_part.Nodes: #Fixed surface
                disp = KratosMultiphysics.Vector(3,0)
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,disp)
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        

    def ComputeRotationMatrix(self):
        """
        This function computes the rotation matrix for the current DoF.
        """
        R = KratosMultiphysics.Matrix(3,3) 
        u = KratosMultiphysics.Vector(3,0)
        u[self.DoF-3] = 1 #Rx=[1,0,0], Ry=[0,1,0], Rz=[0,0,1]
        R[0,0] = cos(self.inf_rot)+u[0]**2*(1-cos(self.inf_rot)) 
        R[0,1] = u[0]*u[1]*(1-cos(self.inf_rot))-u[2]*sin(self.inf_rot)
        R[0,2] = u[0]*u[2]*(1-cos(self.inf_rot))+u[1]*sin(self.inf_rot)
        R[1,0] = u[1]*u[0]*(1-cos(self.inf_rot))+u[2]*sin(self.inf_rot)
        R[1,1] = cos(self.inf_rot)+u[1]**2*(1-cos(self.inf_rot))
        R[1,2] = u[1]*u[2]*(1-cos(self.inf_rot))-u[0]*sin(self.inf_rot)
        R[2,0] = u[2]*u[0]*(1-cos(self.inf_rot))-u[1]*sin(self.inf_rot)
        R[2,1] = u[2]*u[1]*(1-cos(self.inf_rot))+u[0]*sin(self.inf_rot)
        R[2,2] = cos(self.inf_rot)+u[2]**2*(1-cos(self.inf_rot))
        return R

    def MasterStiffnessVector(self):
        """
        Obtains the Stiffness column Vector of the current DoF of the n surface.
        """
        current_step = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.STEP)-1

        if self.DoF==0 and self.n==0: #Only initialize once the master stiffness matrix for each slave surface
            matrix_size = self.number_master_nodes*2*self.dim
            self.master_stiffness = KratosMultiphysics.Matrix(matrix_size,matrix_size,0)#Initialize Master Stiffness Matrix

        resultant = KratosMultiphysics.Matrix(self.number_master_nodes,2*self.dim,0)
        for k in range(self.number_master_nodes):
            surface_model_part = self.model_part.GetSubModelPart(self.slave_surface_name[k]) #Read model part 
            #Initializaing results Matrices an Vectors
            surface_num_nodes = surface_model_part.NumberOfNodes()
            surface_undeformed_coordinates = KratosMultiphysics.Matrix(surface_num_nodes,self.dim)
            surface_reaction = KratosMultiphysics.Matrix(surface_num_nodes,self.dim)
            counter=0
            for node in surface_model_part.Nodes:
                undeformed_coordinates = (node.X0,node.Y0,node.Z0) 
                for i in range(self.dim):
                    surface_undeformed_coordinates[counter,i] = undeformed_coordinates[i]
                    surface_reaction[counter,i] = node.GetSolutionStepValue(KratosMultiphysics.REACTION)[i] #Reactions of the surface
                counter+=1

            r = KratosMultiphysics.Matrix(surface_num_nodes,self.dim) # Initialize Position vector r

            for i in range(surface_num_nodes):
                for j in range(self.dim):
                    r[i,j]=surface_undeformed_coordinates[i,j]-self.master_coor[k,j] # Assign position vector r

            #Calculating the moments about the master node
            moments = self.Cross_product(r,surface_reaction,surface_num_nodes)

            for i in range(surface_num_nodes):
                for j in range(self.dim):
                    resultant[k,j]+=surface_reaction[i,j] #Resultant forces assembling
                    resultant[k,j+self.dim]+=moments[i,j] #Resultant moments assembling

        if self.DoF<3:
            constraint = self.value
        else:
            constraint = self.inf_rot

        #------ Stiffness calculation
        for k in range(self.number_master_nodes):
            for i in range(2*self.dim):
                self.master_stiffness[(2*k*self.dim)+i,current_step] = resultant[k,i]/constraint
    
        if self.DoF==5 and (self.n==self.number_master_nodes-1):
            if self.save_matrix==True:    
                KratosMultiphysics.Logger.Print("\n",self.master_stiffness, label="Master Stiffness: ")
                self.json_parameters=KratosMultiphysics.Parameters()
                self.json_parameters.AddEmptyValue("Stiffness Matrix")
                self.json_parameters["Stiffness Matrix"].SetMatrix(self.master_stiffness)
                self.json_parameters.AddEmptyValue("Displacement")
                self.json_parameters["Displacement"].SetDouble(self.value)
                self.json_parameters.AddEmptyValue("Rotation")
                self.json_parameters["Rotation"].SetDouble(self.inf_rot)
                self.json_parameters.AddEmptyValue("Interface Frame Origins")
                self.json_parameters["Interface Frame Origins"].SetMatrix(self.master_coor)
                with open('Stiffness_Matrix.json', 'w') as parameter_output_file:
                    parameter_output_file.write(self.json_parameters.PrettyPrintJsonString())


    def Cross_product(self,a,b,n):
        """
        Obtains the dot product (used to obtain the moments).
        """
        c = KratosMultiphysics.Matrix(n,self.dim)
        for i in range(n):
            c[i,0] = a[i,1]*b[i,2]-a[i,2]*b[i,1] #Moments x
            c[i,1] = a[i,2]*b[i,0]-a[i,0]*b[i,2] #Moments y
            c[i,2] = a[i,0]*b[i,1]-a[i,1]*b[i,0] #Moments z
        return c