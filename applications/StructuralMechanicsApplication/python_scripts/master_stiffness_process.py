import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import numpy as np 

class MasterStiffnessMatrixProcess(KratosMultiphysics.Process):
    """
    This process allows the user to obtain an Stiffness Matrix of a Master node related to an specific Slave surface.\n
    Example:\n
    Master = MasterStiffnessMatrixProcess(parameters, 0.00001, [[0.85,0.1,0.05]], "Slave")\n
    Master.Run()
    """
    def __init__(self, parameters, Value, Master_coor, slave_surface_name):
        """
        It is needed to provide the model, paramenters, value of the infinitesimal displacement/rotation, 
        Master node's coordinate and Slave's surface name.
        """
        KratosMultiphysics.Process.__init__(self)
        self.parameters = parameters
        self.Value = Value
        self.Master_coor = np.array(Master_coor)
        self.Stiffness = np.empty((6,1)) #Initialize the Stiffness column vector
        self.Slave_name = slave_surface_name
        self.Master_Stiffness=np.empty((6,6)) #Initialize the Master Stiffness matrix


    def ChangeVectorValues(self,DoF):
        """
        Similar to AddVector/ScalarValueProcess, specific for MasterStiffnessMatrixProcess.\n
        Produces an infinitesimal displacement/rotation on each direction.\n
        It is needed to provide "DoF", which stands for the loop on DoFs.
        """
        for i in range(self.parameters["output_processes"]["gid_output"].size()):
            gid_output = self.parameters["output_processes"]["gid_output"][i]
            if gid_output["Parameters"]["model_part_name"].GetString()=="Structure.computing_domain":# Erasing ".computing_domain" from gid_output
                gid_output["Parameters"]["model_part_name"].SetString("Structure")
        for i in range(self.parameters["output_processes"]["vtk_output"].size()):
            self.vtk_output = self.parameters["output_processes"]["vtk_output"][i]
            if self.vtk_output["Parameters"]["model_part_name"].GetString()=="Structure.computing_domain":# Erasing ".computing_domain" from vtk_output
                self.vtk_output["Parameters"]["model_part_name"].SetString("Structure")
                
        #Iterating through the entire list of constraints 
        for j in range(self.parameters["processes"]["constraints_process_list"].size()):
            constraints_process_list = self.parameters["processes"]["constraints_process_list"][j] 
            if constraints_process_list["Parameters"]["model_part_name"].GetString()=="Structure.DISPLACEMENT_"+self.Slave_name: # Enforcing to change the displacement vector only in slave surface
                if DoF==0:
                    constraints_process_list["Parameters"]["value"].SetVector([self.Value, 0.0, 0.0]) #Setting the new vector of displacements.
                elif DoF==1:
                    constraints_process_list["Parameters"]["value"].SetVector([0.0, self.Value, 0.0]) #Setting the new vector of displacements.
                elif DoF==2:
                    constraints_process_list["Parameters"]["value"].SetVector([0.0, 0.0, self.Value]) #Setting the new vector of displacements.
                elif DoF==3:
                    constraints_process_list["Parameters"]["value"][0].SetDouble(0.0) #Since for the rotation we give double values and strings. We have to Set them separately
                    constraints_process_list["Parameters"]["value"][1].SetString("((y-"+str(self.Master_coor[0,1])+")*cos("+str(self.Value)+")-(z-"+str(self.Master_coor[0,2])+")*sin("+str(self.Value)+"))+("+str(self.Master_coor[0,1])+"-y)")
                    constraints_process_list["Parameters"]["value"][2].SetString("((y-"+str(self.Master_coor[0,1])+")*sin("+str(self.Value)+")+(z-"+str(self.Master_coor[0,2])+")*cos("+str(self.Value)+"))+("+str(self.Master_coor[0,2])+"-z)")
                elif DoF==4:
                    constraints_process_list["Parameters"]["value"][0].SetString("((x-"+str(self.Master_coor[0,0])+")*cos("+str(self.Value)+")+(z-"+str(self.Master_coor[0,2])+")*sin("+str(self.Value)+"))+("+str(self.Master_coor[0,0])+"-x)")#Since for the rotation we give double values and strings. We have to Set them separately
                    constraints_process_list["Parameters"]["value"][1].SetDouble(0.0)
                    constraints_process_list["Parameters"]["value"][2].SetString("(-(x-"+str(self.Master_coor[0,0])+")*sin("+str(self.Value)+")+(z-"+str(self.Master_coor[0,2])+")*cos("+str(self.Value)+"))+("+str(self.Master_coor[0,2])+"-z)")
                elif DoF==5:
                    constraints_process_list["Parameters"]["value"][0].SetString("((x-"+str(self.Master_coor[0,0])+")*cos("+str(self.Value)+")-(y-"+str(self.Master_coor[0,1])+")*sin("+str(self.Value)+"))+("+str(self.Master_coor[0,0])+"-x)")#Since for the rotation we give double values and strings. We have to Set them separately
                    constraints_process_list["Parameters"]["value"][1].SetString("((x-"+str(self.Master_coor[0,0])+")*sin("+str(self.Value)+")+(y-"+str(self.Master_coor[0,1])+")*cos("+str(self.Value)+"))+("+str(self.Master_coor[0,1])+"-y)")
                    constraints_process_list["Parameters"]["value"][2].SetDouble(0.0) 

    def MasterStiffnessVector(self):
        """
        Obtains the Stiffness column Vector of the current DoF.
        """

        slave_model_part = self.model["Structure.DISPLACEMENT_"+self.Slave_name] #Read model part (Slave surface)

        slave_deformed_coordinates = [] #Initialize lists for results
        slave_displacement = []
        slave_reaction = []
        for node in slave_model_part.Nodes:
            slave_deformed_coordinates.append([node.X,node.Y,node.Z]) #Deformed coordinates of the slave surface
            slave_displacement.append(node.GetSolutionStepValue(KratosMultiphysics.DISPLACEMENT)) #Displacements of the slave surface
            slave_reaction.append(node.GetSolutionStepValue(KratosMultiphysics.REACTION)) #Reactions of the slave surface
        slave_undeformed_coordinates = np.array(slave_deformed_coordinates)-np.array(slave_displacement)# Convert to numpy arrays
        slave_reaction = np.array(slave_reaction)

        dim = self.parameters["solver_settings"]["domain_size"].GetInt() #Dimension 

        r = np.zeros((len(slave_undeformed_coordinates),dim)) # Initialize Position vector r 

        for i in range(len(slave_undeformed_coordinates)):
            r[i,:] = slave_undeformed_coordinates[i,:]-self.Master_coor # Assign position vector r

        Moments=np.cross(r,slave_reaction) #Calculating the moments about the master node

        Resultant=[]
        for i in range(dim):
            Resultant.append(sum(slave_reaction[:,i])) #Resultant forces assembling
        for i in range(dim):
            Resultant.append(sum(Moments[:,i])) #Resultant moments assembling

        #------ Stiffness calculation
        Stiffness=[]#Only for 3D cases, change to dim+1 for 2D cases and dim for 1D cases.
        for i in range(dim*2):
            Stiffness.append(Resultant[i]/self.Value)#Obtaining the Stiffness column for the degree of fredoom

        self.Stiffness=Stiffness

    def RunSimulation(self):
        """
        Initialize a new model.\n
        Instantiation of StructuralMechanicsAnalysis.\n
        Runs the simulation.\n
        """
        self.model = KratosMultiphysics.Model() 
        simulation = StructuralMechanicsAnalysis(self.model, self.parameters)
        simulation.Run()

    def Run(self):
        """
        Runs the different simulations needed to compute the stiffness matrix of a master node related to an slave surface.\n
        ChangeVectorValues(w): Changes the vector values for the displacement/rotation of the current DoF.\n
        RunSimulation(): Instance a new model and runs a simulation (StructuralMechanicsAnalysis) with the new parameters.\n
        MasterStiffnessVector(): Obtains the stiffness vector for the current DoF (w)
        """
        for DoF in range(6):
            self.ChangeVectorValues(DoF)
            self.RunSimulation()
            self.MasterStiffnessVector()
            self.Master_Stiffness[:,DoF] = self.Stiffness #Assembling of Master Stiffness Matrix
        
        print("\n \n The Stiffness Matrix related to the master node is: \n", self.Master_Stiffness)