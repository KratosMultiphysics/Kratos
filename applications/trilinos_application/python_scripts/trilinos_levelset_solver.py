# -*- coding: utf-8 -*-
#importing the Kratos Library
from Kratos import *
from KratosConvectionDiffusionApplication import *
from KratosTrilinosApplication import *
from KratosMeshingApplication import *
import trilinos_convdiff_solver
import trilinos_monolithic_solver_eulerian
import mpi

#settings for the convection solver
distance_settings = ConvectionDiffusionSettings()
distance_settings.SetDensityVariable(DENSITY_AIR)
distance_settings.SetDiffusionVariable(CONDUCTIVITY)
distance_settings.SetUnknownVariable(DISTANCE)
distance_settings.SetVolumeSourceVariable(HEAT_FLUX)
distance_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
distance_settings.SetMeshVelocityVariable(MESH_VELOCITY)

def AddVariables(model_part ):
    print distance_settings
    #variables needed for the fluid solver
    trilinos_monolithic_solver_eulerian.AddVariables(model_part)
    model_part.AddNodalSolutionStepVariable(DISTANCE)

    #variables needed for the distance solver
    trilinos_convdiff_solver.AddVariables(model_part,distance_settings)

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

def AddDofs(model_part):
    trilinos_monolithic_solver_eulerian.AddDofs(model_part)
    trilinos_convdiff_solver.AddDofs(model_part,distance_settings)

    print "variables for the convection diffusion solver added correctly"


class TrilinosLevelSetSolver:
    
    def __init__(self,model_part,domain_size):
        self.model_part = model_part
        self.domain_size = domain_size

        #construct the model part for the convection solver
        if(self.domain_size == 2):
            conv_elem = "ConvDiff2D"
            conv_cond = "Condition2D"
        else:
            conv_elem = "ConvDiff3D"
            conv_cond = "Condition3D"
        self.convection_model_part = ModelPart("convection_model_part")
        self.conv_generator = ConnectivityPreserveModeler()
        (self.conv_generator).GenerateModelPart(self.model_part,self.convection_model_part,conv_elem,conv_cond)
        
        
        (ParallelFillCommunicator(self.convection_model_part)).Execute();
        print "ln52"

        #constructing the convection solver for the distance
        self.convection_solver = trilinos_convdiff_solver.ConvectionDiffusionSolver(self.convection_model_part,self.domain_size,distance_settings)

        #constructing the fluid solver
        self.fluid_solver = trilinos_monolithic_solver_eulerian.MonolithicSolver(self.model_part,self.domain_size)
        self.vel_criteria = 1e-3
        self.press_criteria = 1e-4
        self.vel_abs_criteria = 1e-9
        self.press_abs_criteria = 1e-9
        self.fluid_solver.ReformDofSetAtEachStep = False

        #select the densities of the two fluids
        self.rho1 = 1000.0 #applied on the negative part of the domain
        self.rho2 = 100.0 #applied to the positive part of the domain

        if(self.domain_size == 2):
            self.redistance_utils = ParallelDistanceCalculator2D()
        else:
            self.redistance_utils = ParallelDistanceCalculator3D()

        self.max_levels = 10

        self.max_edge_size = self.redistance_utils.FindMaximumEdgeSize(self.convection_model_part)
        self.max_distance = self.max_edge_size * 3.0;

        #assigning the fluid properties
        conductivity = 0.0;
        density = 1.0;
        specific_heat = 1.0;
        for node in model_part.Nodes:
            node.SetSolutionStepValue(CONDUCTIVITY,0,conductivity);
            node.SetSolutionStepValue(DENSITY_AIR,0,density);
            node.SetSolutionStepValue(SPECIFIC_HEAT,0,specific_heat);

        self.max_ns_iterations = 20
        self.dynamic_tau = 0.001







    def Initialize(self):
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
        self.fluid_solver.vel_criteria = self.vel_criteria
        self.fluid_solver.press_criteria = self.press_criteria
        self.fluid_solver.vel_abs_criteria = self.vel_abs_criteria
        self.fluid_solver.press_abs_criteria = self.press_abs_criteria
        self.fluid_solver.max_iter = self.max_ns_iterations

        self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA,self.max_levels,self.max_distance)
        self.convection_solver.Initialize()
        self.fluid_solver.Initialize()
                 
   
    def Solve(self):
        self.convection_model_part.ProcessInfo = self.model_part.ProcessInfo
        (self.convection_model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS,distance_settings)
        if(mpi.rank == 0):
            print "line 106"


        
        #redistance if required
        self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA,self.max_levels,self.max_distance)
        if(mpi.rank == 0):
            print "line 113"

        #convect distance
        (self.convection_solver).Solve()
        if(mpi.rank == 0):
            print "line 118"

        #apply density
        for node in self.model_part.Nodes:
            dist = node.GetSolutionStepValue(DISTANCE)
            if(dist <= 0):
                node.SetSolutionStepValue(DENSITY,0,self.rho1)
            else:
                node.SetSolutionStepValue(DENSITY,0,self.rho2)
                
        #solve fluid
        (self.fluid_solver).Solve()
        if(mpi.rank == 0):
            print "line 131"


