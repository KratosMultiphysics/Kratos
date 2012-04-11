#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
CheckForPreviousImport()

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);


def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
        
    print "variables for the mesh solver added correctly"


class MeshSolver:
    
    def __init__(self,model_part,domain_size,reform_dof_at_every_step):

        self.model_part = model_part
        self.domain_size = domain_size
        self.reform_dof_at_every_step = reform_dof_at_every_step
        
        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        
        #assignation of parameters to be used
        self.time_order = 2
        
        #definition of the solvers
        pILUPrecond = ILU0Preconditioner()
        self.linear_solver =  BICGSTABSolver(1e-5, 300,pILUPrecond)
##        pDiagPrecond = DiagonalPreconditioner()
##        self.linear_solver = CGSolver(1e-3, 300, pDiagPrecond)

    def Initialize(self):
        (self.neighbour_search).Execute()

        if(self.domain_size == 2):
            self.solver = BallVertexMeshMoving2D()
        else:
            self.solver = BallVertexMeshMoving3D()

        self.move_mesh_utilities = MoveMeshUtilities();

        if(self.reform_dof_at_every_step == False):
            self.solver.ConstructSystem(self.model_part);
            
                 
   
    def Solve(self):
        if(self.reform_dof_at_every_step == True):
            (self.neighbour_search).Execute()

            self.solver.ConstructSystem(self.model_part);

            self.solver.BuildAndSolveSystem(self.model_part,self.linear_solver);

            self.solver.ClearSystem()
        else:
            self.solver.BuildAndSolveSystem(self.model_part,self.linear_solver);

        ##move mesh and calculate mesh velocity
        self.move_mesh_utilities.BDF_MoveMesh(self.time_order, self.model_part)


    def MoveNodes(self):
        print "**************************  *********************"
        self.move_mesh_utilities.BDF_MoveMesh(self.time_order, self.model_part)
