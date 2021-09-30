import KratosMultiphysics
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from math import cos, sin 
import json

def Factory(parameters, model):
    if not isinstance(parameters, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return UpdateBoundaryCondition(model, parameters["Parameters"])

## All the processes python should be derived from "Process"
class UpdateBoundaryCondition(KratosMultiphysics.Process):
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
                "model_part_list": ["Surface_name_1","Surface_name_2"],
                "list_of_master_coordinates": [[0.0,0.0,0.0],[1.0,0.0,0.0]],
                "boundary_condition": [0.0,0.0,0.0],
                "time_step_range": [0,1]
            }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.parameters = parameters
        self.model_part = model["Structure"]
        self.hrom_model_part = self.model_part.GetSubModelPart("COMPUTE_HROM")
        self.master_coor = self.parameters["list_of_master_coordinates"].GetMatrix()
        self.slave_surface_name = self.parameters["model_part_list"].GetStringArray()
        self.number_master_nodes = len(self.slave_surface_name)
        self.dim = self.model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)
        


    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step
        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.UpdateBCs()

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step
        Keyword arguments:
        self -- It signifies an instance of a class.
        """

    def UpdateBCs(self):
        """
        Similar to AddVector/ScalarValueProcess, specific for MasterStiffnessMatrixProcess.\n
        Produces an infinitesimal displacement/rotation on each direction.\n
        It is needed to provide "DoF", which stands for the loop on DoFs and n to recognize which surface Stiffness Matrix to obtain.
        """
        self.boundary_condition = self.parameters["boundary_condition"].GetVector()
        slave_submodel_part = self.hrom_model_part.GetSubModelPart(self.slave_surface_name[0])#Slave surface
        disp = KratosMultiphysics.Vector(3,0)
        for node in slave_submodel_part.Nodes: #Slave surface
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,disp)
            node.Fix(KratosMultiphysics.DISPLACEMENT_X)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
            node.Fix(KratosMultiphysics.DISPLACEMENT_Z)

        for i in range(self.number_master_nodes-1):
            slave_submodel_part = self.hrom_model_part.GetSubModelPart(self.slave_surface_name[i+1])#Slave surface.
            disp = KratosMultiphysics.Vector(self.dim)
            for j in range(self.dim):
                disp[j] = self.boundary_condition[6*i+j]
            Rx = self.ComputeRotationMatrix(0,self.boundary_condition[6*i+3])
            Ry = self.ComputeRotationMatrix(1,self.boundary_condition[6*i+4])
            Rz = self.ComputeRotationMatrix(2,self.boundary_condition[6*i+5])
            v = KratosMultiphysics.Vector(3)
            for node in slave_submodel_part.Nodes: #Slave surface
                v[0] = node.X0 - self.master_coor[i+1,0]
                v[1] = node.Y0 - self.master_coor[i+1,1]
                v[2] = node.Z0 - self.master_coor[i+1,2]
                rot_x = Rx*v-v
                rot_y = Ry*v-v
                rot_z = Rz*v-v
                node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT,disp+rot_x+rot_y+rot_z)
                node.Fix(KratosMultiphysics.DISPLACEMENT_X)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Y)
                node.Fix(KratosMultiphysics.DISPLACEMENT_Z)


    def ComputeRotationMatrix(self,DoF,inf_rot):
        """
        This function computes the rotation matrix for the current DoF.
        """
        R = KratosMultiphysics.Matrix(3,3) 
        u = KratosMultiphysics.Vector(3,0)
        u[DoF] = 1 #Rx=[1,0,0], Ry=[0,1,0], Rz=[0,0,1]
        R[0,0] = cos(inf_rot)+u[0]**2*(1-cos(inf_rot)) 
        R[0,1] = u[0]*u[1]*(1-cos(inf_rot))-u[2]*sin(inf_rot)
        R[0,2] = u[0]*u[2]*(1-cos(inf_rot))+u[1]*sin(inf_rot)
        R[1,0] = u[1]*u[0]*(1-cos(inf_rot))+u[2]*sin(inf_rot)
        R[1,1] = cos(inf_rot)+u[1]**2*(1-cos(inf_rot))
        R[1,2] = u[1]*u[2]*(1-cos(inf_rot))-u[0]*sin(inf_rot)
        R[2,0] = u[2]*u[0]*(1-cos(inf_rot))-u[1]*sin(inf_rot)
        R[2,1] = u[2]*u[1]*(1-cos(inf_rot))+u[0]*sin(inf_rot)
        R[2,2] = cos(inf_rot)+u[2]**2*(1-cos(inf_rot))
        return R