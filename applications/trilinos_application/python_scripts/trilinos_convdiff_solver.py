# -*- coding: utf-8 -*-
#importing the Kratos Library
from Kratos import *
from KratosConvectionDiffusionApplication import *
from KratosTrilinosApplication import *

def AddVariables(model_part,settings ):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    #model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    model_part.AddNodalSolutionStepVariable(settings.GetMeshVelocityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetUnknownVariable());
    model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(SPECIFIC_HEAT);
    model_part.AddNodalSolutionStepVariable(settings.GetVolumeSourceVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDensityVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetDiffusionVariable());
    model_part.AddNodalSolutionStepVariable(settings.GetSurfaceSourceVariable());
    model_part.AddNodalSolutionStepVariable(CONVECTION_COEFFICIENT);

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

def AddDofs(model_part,settings):
    for node in model_part.Nodes:
        node.AddDof(settings.GetUnknownVariable());

    print "variables for the convection diffusion solver added correctly"


class ConvectionDiffusionSolver:
    
    def __init__(self,model_part,domain_size,my_settings):

	self.settings = my_settings

        self.model_part = model_part
        self.domain_size = domain_size

        #assignation of parameters to be used
        self.time_order = 2;
        self.prediction_order = 2;
        self.ReformDofAtEachIteration = False; 

        self.echo_level = 0

        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver","AZ_gmres");
        aztec_parameters.set("AZ_kspace",100);
        aztec_parameters.set("AZ_output","AZ_none");
        preconditioner_type = "ILU"
        preconditioner_parameters = ParameterList()
        overlap_level = 0
        nit_max = 1000
        linear_tol = 1e-9
        self.linear_solver =  AztecSolver(aztec_parameters,preconditioner_type,preconditioner_parameters,linear_tol,nit_max,overlap_level);

    def Initialize(self):
        self.model_part.ProcessInfo
	(self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)

        self.Comm = CreateCommunicator()

        self.solver = TrilinosConvectionDiffusionStrategy(self.Comm,self.model_part,self.linear_solver,self.ReformDofAtEachIteration,self.time_order,self.prediction_order)
#        (self.solver).SetEchoLevel(self.echo_level)
        print "finished initialization of the fluid strategy"
                 
   
    def Solve(self):     
        (self.model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,self.settings)

        (self.solver).Solve()

