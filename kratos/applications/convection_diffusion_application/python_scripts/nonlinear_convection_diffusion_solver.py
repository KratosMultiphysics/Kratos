#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionApplication import *
CheckForPreviousImport()


def AddVariables(model_part,settings):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetProjectionVariable());
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable());

    model_part.AddNodalSolutionStepVariable(ENTHALPY);
    model_part.AddNodalSolutionStepVariable(LATENT_HEAT);
    model_part.AddNodalSolutionStepVariable(MELT_TEMPERATURE_1);
    model_part.AddNodalSolutionStepVariable(MELT_TEMPERATURE_2);

def AddDofs(model_part,settings):
    for node in model_part.Nodes:
        node.AddDof(settings.GetUnknownVariable());

    print "variables for the convection diffusion solver added correctly"


class ConvectionDiffusionSolver:
    
    def __init__(self,model_part,domain_size,my_settings):

        #neighbour search
	self.settings = my_settings
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.time_order = 1;
	self.prediction_order = 1;
        self.ReformDofAtEachIteration = False;
        self.max_iter = 15;
        self.toll = 1e-9;

        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)


    def Initialize(self):
        (self.neighbour_search).Execute()
        self.model_part.ProcessInfo
	(self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)
        self.solver = ResidualBasedConvectionDiffusionStrategyNonLinear			           (self.model_part,self.linear_solver,self.ReformDofAtEachIteration,self.time_order,self.max_iter,self.toll)   
        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
                 
   
    def Solve(self):
        if(self.ReformDofAtEachIteration == True):
            (self.neighbour_search).Execute()        
        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)
        (self.solver).Solve()

        if(self.ReformDofAtEachIteration == True):
            (self.solver).Clear()      
       

