#importing the Kratos Library
from Kratos import *
from KratosConvectionDiffusionApplication import *

def AddVariables(model_part,settings):
    print "user should include the variables as needed over the model_part"
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable());
    model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable());

def AddDofs(model_part,settings):
    for node in model_part.Nodes:
        node.AddDof(settings.GetUnknownVariable());
    print "user should include the variables as needed by nodes"


class PureConvectionSolver:
    
    def __init__(self,model_part,domain_size,my_settings):

        #neighbour search
	self.settings = my_settings
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.time_order = 2;
	self.prediction_order = 1;
        self.ReformDofAtEachIteration = True; 
        self.reform_convection_matrix = True;

        ##Variable to be convected:
        ## #1 = TEPERATURE
        ## #2 = DISTANCE
        self.scalar_var_convected = 1
        
##        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
##        self.linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
##        pILUPrecond = ILU0Preconditioner()
        self.linear_solver =  BICGSTABSolver(1e-9, 5000,pDiagPrecond)

        ##pure convection tool
        if(self.domain_size == 2):
            self.convection_solver = PureConvectionCrankNUtilities2D();
        else: 
            self.convection_solver = PureConvectionCrankNUtilities3D();

        self.first_initialize_performed = False



    def Initialize(self):
        print "Finishing Initialize"  
                
   
##    def Solve(self):
##        print self.scalar_var_convected
##        if(self.ReformDofAtEachIteration == True):
##            (self.neighbour_search).Execute()        
##        if(self.scalar_var_convected == 1):
##            # construct system -- could be done once if the mesh does not change
##            self.convection_solver.ConstructSystem(self.model_part,TEMPERATURE,VELOCITY,MESH_VELOCITY);
##
##            #calculate projections
##            self.convection_solver.CalculateProjection(self.model_part,TEMPERATURE,NODAL_AREA,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
##                
##            #perform convection step
##            self.convection_solver.ConvectScalarVar(self.model_part,self.linear_solver,TEMPERATURE,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.time_order);
##        elif(self.scalar_var_convected == 2):
##            # construct system -- could be done once if the mesh does not change
##            self.convection_solver.ConstructSystem(self.model_part,DISTANCE,VELOCITY,MESH_VELOCITY);
##
##            #calculate projections
##            self.convection_solver.CalculateProjection(self.model_part,DISTANCE,NODAL_AREA,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);
##                
##            #perform convection step
##            self.convection_solver.ConvectScalarVar(self.model_part,self.linear_solver, DISTANCE,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.time_order);
##        else:
##            print "ERROR: choose the correct variable to convect #1:TEMPERATURE #2:DISTANCE"
##        #free memory
##        self.convection_solver.ClearSystem() 

    def Solve(self, scalar_variable):
	(self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)
        if(self.ReformDofAtEachIteration == True or self.first_initialize_performed == False):
            (self.neighbour_search).Execute()
        
            # construct system -- could be done once if the mesh does not change
            self.convection_solver.ConstructSystem(self.model_part,scalar_variable,VELOCITY,MESH_VELOCITY);

            self.first_initialize_performed = True

        #calculate projections
        self.convection_solver.CalculateProjection(self.model_part,scalar_variable,NODAL_AREA,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ);

        print "97"
        #perform convection step
        self.convection_solver.ConvectScalarVar(self.model_part,self.linear_solver,scalar_variable,VELOCITY,MESH_VELOCITY,TEMP_CONV_PROJ,self.time_order);

        print "101"
        #free memory
        if(self.ReformDofAtEachIteration == True):
            self.convection_solver.ClearSystem()
