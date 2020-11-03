import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis


class MasterStiffnessMatrixProcess(KratosMultiphysics.Process):
    """
    This process allows the user to obtain the Stiffness Matrices of two Master nodes related to two specific Slave surfaces.\n
    Example:\n
    Master = MasterStiffnessMatrixProcess(model, parameters)\n
    Master.Run()\n
    """
    #KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING) #In case we want to avoid printing

    def __init__(self, model, parameters):
        """
        It is needed to provide the parameters, with a "list_values" of the infinitesimal displacement/rotations, "list_of_master_coordinates" as a Matrix, and "processes_sub_model_part_list" which contains the names of the two slaves srufaces\n.
        """
        KratosMultiphysics.Process.__init__(self)
        self.parameters = parameters
        self.value = self.parameters["list_of_values"].GetVector()
        self.master_coor = self.parameters["list_of_master_coordinates"].GetMatrix()
        self.slave_surface_name = self.parameters["solver_settings"]["processes_sub_model_part_list"].GetStringArray()
        self.number_master_nodes = len(self.slave_surface_name)
        self.dim = self.parameters["solver_settings"]["domain_size"].GetInt()
        self.ValidateData()
            
    def ValidateData(self):
        """
        Validates the Data for this specific process.
        """
        if (type(self.parameters) != KratosMultiphysics.Parameters):
            raise Exception("Expected input -parameters- should be a Parameters object, encapsulating a json string.")
        if (type(self.value) != KratosMultiphysics.Vector):
            raise Exception("Expected input -value- should be a list.")
        if (type(self.master_coor) != KratosMultiphysics.Matrix):
            raise Exception("Expected input -master_coor- should be a list.")
        if (type(self.slave_surface_name) != list):
            raise Exception("Expected input -slave_surface_name- should be a list.")
        if (len(self.value)!=2 or self.master_coor.Size1()!=2 or self.master_coor.Size2()!=3 or len(self.slave_surface_name)!=2):
            raise Exception("Lists parameters should provide 2 surfaces information (length 2) for a 3D coordinate system.")

    def ChangeVectorValues(self,DoF,n):
        """
        Similar to AddVector/ScalarValueProcess, specific for MasterStiffnessMatrixProcess.\n
        Produces an infinitesimal displacement/rotation on each direction.\n
        It is needed to provide "DoF", which stands for the loop on DoFs and n to recognize which surface Stiffness Matrix to obtain.
        """
        
        if n==1 and DoF==0:
            self.slave_surface_name.reverse()
        slave_submodel_part_name = "Structure."+self.slave_surface_name[0]#Slave surface, the one we want to obtain the stiffness from.
        fixed_submodel_part_name = "Structure."+self.slave_surface_name[1]#Fixed surface.
        
        for i in range(self.parameters["output_processes"]["gid_output"].size()):
            gid_output = self.parameters["output_processes"]["gid_output"][i]
            gid_output_model_part_name = gid_output["Parameters"]["model_part_name"].GetString()
            if gid_output_model_part_name=="Structure.computing_domain":# Erasing ".computing_domain" from gid_output
                gid_output["Parameters"]["model_part_name"].SetString("Structure")
        for i in range(self.parameters["output_processes"]["vtk_output"].size()):
            vtk_output = self.parameters["output_processes"]["vtk_output"][i]
            vtk_output_model_part_name = vtk_output["Parameters"]["model_part_name"].GetString()
            if vtk_output_model_part_name=="Structure.computing_domain":# Erasing ".computing_domain" from vtk_output
                vtk_output["Parameters"]["model_part_name"].SetString("Structure")
                
        #Iterating through the entire list of constraints 
        for j in range(self.parameters["processes"]["constraints_process_list"].size()):
            constraints_process_list = self.parameters["processes"]["constraints_process_list"][j] 
            current_submodel_part_name = constraints_process_list["Parameters"]["model_part_name"].GetString()
            #------------------- I am working to change this part --------------------------
            #------------------- I am working to change this part --------------------------
            if current_submodel_part_name==slave_submodel_part_name: # Forcing to change the displacement vector only in slave surface
                if DoF==0:
                    constraints_process_list["Parameters"]["value"].SetVector([self.value[n], 0.0, 0.0]) #Setting the new vector of displacements.
                elif DoF==1:
                    constraints_process_list["Parameters"]["value"].SetVector([0.0, self.value[n], 0.0]) #Setting the new vector of displacements.
                elif DoF==2:
                    constraints_process_list["Parameters"]["value"].SetVector([0.0, 0.0, self.value[n]]) #Setting the new vector of displacements.
                elif DoF==3:
                    constraints_process_list["Parameters"]["value"][0].SetDouble(0.0) #Since for the rotation we give double values and strings. We have to Set them separately
                    constraints_process_list["Parameters"]["value"][1].SetString("((y-"+str(self.master_coor[n,1])+")*cos("+str(self.value[n])+")-(z-"+str(self.master_coor[n,2])+")*sin("+str(self.value[n])+"))+("+str(self.master_coor[n,1])+"-y)")
                    constraints_process_list["Parameters"]["value"][2].SetString("((y-"+str(self.master_coor[n,1])+")*sin("+str(self.value[n])+")+(z-"+str(self.master_coor[n,2])+")*cos("+str(self.value[n])+"))+("+str(self.master_coor[n,2])+"-z)")
                elif DoF==4:
                    constraints_process_list["Parameters"]["value"][0].SetString("((x-"+str(self.master_coor[n,0])+")*cos("+str(self.value[n])+")+(z-"+str(self.master_coor[n,2])+")*sin("+str(self.value[n])+"))+("+str(self.master_coor[n,0])+"-x)")#Since for the rotation we give double values and strings. We have to Set them separately
                    constraints_process_list["Parameters"]["value"][1].SetDouble(0.0)
                    constraints_process_list["Parameters"]["value"][2].SetString("(-(x-"+str(self.master_coor[n,0])+")*sin("+str(self.value[n])+")+(z-"+str(self.master_coor[n,2])+")*cos("+str(self.value[n])+"))+("+str(self.master_coor[n,2])+"-z)")
                elif DoF==5:
                    constraints_process_list["Parameters"]["value"][0].SetString("((x-"+str(self.master_coor[n,0])+")*cos("+str(self.value[n])+")-(y-"+str(self.master_coor[n,1])+")*sin("+str(self.value[n])+"))+("+str(self.master_coor[n,0])+"-x)")#Since for the rotation we give double values and strings. We have to Set them separately
                    constraints_process_list["Parameters"]["value"][1].SetString("((x-"+str(self.master_coor[n,0])+")*sin("+str(self.value[n])+")+(y-"+str(self.master_coor[n,1])+")*cos("+str(self.value[n])+"))+("+str(self.master_coor[n,1])+"-y)")
                    constraints_process_list["Parameters"]["value"][2].SetDouble(0.0) 
            elif current_submodel_part_name==fixed_submodel_part_name: # Enforcing to change the displacement vector only in fixed surface
                constraints_process_list["Parameters"]["value"].SetVector([0.0, 0.0, 0.0]) #Setting the fixed vector of displacements.
            #------------------- I am working to change this part --------------------------
            #------------------- I am working to change this part --------------------------

    def MasterStiffnessVector(self,DoF,n):
        """
        Obtains the Stiffness column Vector of the current DoF of the n surface.
        """
        slave_model_part = self.model["Structure."+self.slave_surface_name[0]] #Read model part (Slave surface)

        if DoF==0: #Only initialize once the master stiffness matrix for each slave surface
            self.master_stiffness = KratosMultiphysics.Matrix(2*self.dim,2*self.dim,0)#Initialize Master Stiffness Matrix

        #Initializaing results Matrices an Vectors
        slave_num_nodes = slave_model_part.NumberOfNodes()
        slave_deformed_coordinates = KratosMultiphysics.Matrix(slave_num_nodes,self.dim)
        slave_displacement = KratosMultiphysics.Matrix(slave_num_nodes,self.dim)
        slave_reaction = KratosMultiphysics.Matrix(slave_num_nodes,self.dim)
        counter=0
        for node in slave_model_part.Nodes:
            deformed_coordinates = (node.X,node.Y,node.Z) #Deformed coordinates of the slave surface
            for i in range(self.dim):
                slave_deformed_coordinates[counter,i] = deformed_coordinates[i]
                slave_displacement[counter,i] = node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)[i] #Displacements of the slave surface
                slave_reaction[counter,i] = node.GetSolutionStepValue(KratosMultiphysics.REACTION)[i] #Reactions of the slave surface
            counter+=1
        
        slave_undeformed_coordinates = slave_deformed_coordinates-slave_displacement #Undeformed coordinates

        r = KratosMultiphysics.Matrix(slave_num_nodes,self.dim) # Initialize Position vector r

        for i in range(slave_num_nodes):
            for j in range(self.dim):
                r[i,j]=slave_undeformed_coordinates[i,j]-self.master_coor[n,j] # Assign position vector r

        #Calculating the moments about the master node
        moments = KratosMultiphysics.Matrix(slave_num_nodes,self.dim)

        #Cross product (Obtaining moments)
        for i in range(slave_num_nodes):
            moments[i,0] = r[i,1]*slave_reaction[i,2]-r[i,2]*slave_reaction[i,1] #Moments x
            moments[i,1] = r[i,2]*slave_reaction[i,0]-r[i,0]*slave_reaction[i,2] #Moments y
            moments[i,2] = r[i,0]*slave_reaction[i,1]-r[i,1]*slave_reaction[i,0] #Moments z
        
        resultant = KratosMultiphysics.Vector(2*self.dim,0)
        for i in range(slave_num_nodes):
            for j in range(self.dim):
                resultant[j]+=slave_reaction[i,j] #Resultant forces assembling
                resultant[j+self.dim]+=moments[i,j] #Resultant moments assembling

        #------ Stiffness calculation
        for i in range(self.dim*2):
            self.master_stiffness[i,DoF] = resultant[i]/self.value[n]


    def RunSimulation(self):
        """
        Initialize a new model.\n
        Instantiation of StructuralMechanicsAnalysis.\n
        Runs the simulation.\n
        """
        simulation = StructuralMechanicsAnalysis(self.model, self.parameters)
        simulation.Run()

    def Run(self):
        """
        Runs the different simulations needed to compute the stiffness matrix of two master nodes related to two slave surfaces.\n
        ChangeVectorValues(w): Changes the vector values for the displacement/rotation of the current DoF and slave surface (n).\n
        RunSimulation(): Instance a new model and runs a simulation (StructuralMechanicsAnalysis) with the new parameters.\n
        MasterStiffnessVector(): Obtains the stiffness vector for the current DoF and slave surface(n).
        """
        for n in range(self.number_master_nodes): #Loop in both slave surfaces.
            for DoF in range(6):
                self.model = KratosMultiphysics.Model()
                self.ChangeVectorValues(DoF,n) 
                self.RunSimulation()
                self.MasterStiffnessVector(DoF,n)
            KratosMultiphysics.Logger.Print(self.slave_surface_name[0],"\n",self.master_stiffness, label="Master Stiffness Matrix") 
            