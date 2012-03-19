# -*- coding: utf-8 -*-
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
import edgebased_eulerian_solver

def AddVariables(fluid_model_part,particle_model_part):
    edgebased_eulerian_solver.AddVariables(fluid_model_part)
    edgebased_eulerian_solver.AddVariables(particle_model_part)
    fluid_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    fluid_model_part.AddNodalSolutionStepVariable(VISCOSITY)
    fluid_model_part.AddNodalSolutionStepVariable(DENSITY)
    fluid_model_part.AddNodalSolutionStepVariable(NODAL_H)
    fluid_model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    fluid_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    fluid_model_part.AddNodalSolutionStepVariable(FORCE)
    fluid_model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    
    particle_model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    particle_model_part.AddNodalSolutionStepVariable(NODAL_H)
    particle_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
    particle_model_part.AddNodalSolutionStepVariable(FORCE)

def AddDofs(fluid_model_part,particle_model_part):
    edgebased_eulerian_solver.AddDofs(fluid_model_part)
    #edgebased_eulerian_solver.AddDofs(particle_model_part)

    print "dofs for the xivas solver added correctly"
    

def ReadRestartFile(FileName,nodes):
   NODES = nodes
   aaa = open(FileName)
   for line in aaa:
       exec(line)
          
import math

class XIVASSolver:
    
    def __init__(self,fluid_model_part,particle_model_part,domain_size,body_force,viscosity,density):
        self.fluid_model_part = fluid_model_part
        self.particle_model_part = particle_model_part
        self.domain_size = domain_size

        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(self.fluid_model_part,number_of_avg_elems,number_of_avg_nodes)


        #assignation of parameters to be used
	self.assume_constant_pressure =  False
	self.stabdt_pressure_factor = 1.0
	self.use_mass_correction = False
        self.echo_level = 0

        #definition of the solvers
        pDiagPrecond = DiagonalPreconditioner()
        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)

        self.compute_reactions=True
        
	self.particle_utils = LagrangianUtils2D()

	self.node_locator = BinBasedFastPointLocator2D(self.fluid_model_part)
	
	self.body_force = body_force
	self.bf = Array3()
	self.bf[0] = body_force[0]; self.bf[1] = body_force[1]; self.bf[2] = body_force[2];
	self.viscosity = viscosity
	self.density = density
	self.substeps = 2.0
	self.restart_with_eulerian_vel = False
	
	self.max_particles_in_element = 7 #7
	self.min_particles_in_element = 3 #1
	
	self.perform_cheap_correction_step = True
	#self.implicit_viscous_correction = True
	
	#assign nodal H
	aux = CalculateNodalAreaProcess(fluid_model_part,domain_size)
	aux.Execute()
	for node in fluid_model_part.Nodes:
	    a = node.GetSolutionStepValue(NODAL_AREA);
	    h = math.sqrt(2.0*a);
	    if(h == 0):
	      print "node ", node.Id, " has zero h"
	      raise "error, node found with 0 h"
	    node.SetSolutionStepValue(NODAL_H,0,h)
        


    def Initialize(self):
        (self.neighbour_search).Execute()

        self.domain_size = int(self.domain_size)
        
	self.fluid_solver = edgebased_eulerian_solver.EdgeBasedLevelSetSolver(self.fluid_model_part,self.domain_size,self.body_force,self.viscosity,self.density)
	self.fluid_solver.assume_constant_pressure =  True
	self.fluid_solver.stabdt_pressure_factor = 0.01
	self.fluid_solver.use_mass_correction = False

	self.fluid_solver.pressure_linear_solver=self.pressure_linear_solver
	self.fluid_solver.Initialize() 
	
        #(self.fluid_solver).SetEchoLevel(self.echo_level)
        hmin = (self.fluid_solver.fluid_solver).ComputeMinimum_Havg()

        print "minimum nodal havg found on the mesh = ",hmin
        #self.node_locator.UpdateSearchDatabase()  
        self.node_locator.UpdateSearchDatabaseAssignedSize(hmin)
        
        
        self.particle_utils.Reseed(self.fluid_model_part,self.particle_model_part)


        print "finished initialization of the xivas solver"
        
    def EstimateTimeStep(self,safety_factor,max_Dt):
	return self.fluid_solver.EstimateTimeStep(safety_factor,max_Dt)        
   
    def Solve(self):
      
	Dt = self.particle_model_part.ProcessInfo[DELTA_TIME]
	Dt_check = self.fluid_model_part.ProcessInfo[DELTA_TIME]
	if(Dt != Dt_check):
	  raise "error, time step for particle_model_part is not appropriately cloned (not syncronized with fluid_model_part)"
        
	(self.fluid_solver.fluid_solver).ComputeViscousForces()
	print "ccc"
	#self.particle_utils.StreamlineMove(self.bf,self.density,Dt,self.substeps,self.fluid_model_part,self.particle_model_part,self.restart_with_eulerian_vel,self.node_locator)
	self.particle_utils.StreamlineMove(self.bf,self.density,Dt,self.fluid_model_part,self.particle_model_part,self.restart_with_eulerian_vel,self.node_locator)
		
	self.particle_utils.TransferToEulerianMeshShapeBased(self.fluid_model_part,self.particle_model_part,self.node_locator)
	
	#if(self.implicit_viscous_correction == True):
	  #(self.fluid_solver.fluid_solver).ViscosityCorrectionStep()

        (self.fluid_solver.fluid_solver).ComputePressureStabilization()

	(self.fluid_solver.fluid_solver).SolveStep2(self.pressure_linear_solver);
        
        (self.fluid_solver.fluid_solver).SolveStep3();

	(self.fluid_solver.fluid_solver).ComputeViscousForces()
	
	if(self.perform_cheap_correction_step == True):
	  self.particle_utils.StreamlineCorrect(self.density,Dt,self.fluid_model_part,self.particle_model_part,self.node_locator)
	  
	else:
	  self.particle_utils.StreamlineMove(self.bf,self.density,Dt,self.substeps,self.fluid_model_part,self.particle_model_part,self.restart_with_eulerian_vel,self.node_locator)
	  
	self.particle_utils.ReseedEmptyElements(self.fluid_model_part,self.particle_model_part,self.node_locator,self.min_particles_in_element,self.max_particles_in_element)
        
        if(self.compute_reactions == True):
	    exclude_convection_terms = True
            (self.fluid_solver.fluid_solver).ComputeReactions(exclude_convection_terms)

    def Clear(self):
        (self.fluid_solver.fluid_solver).Clear()
         

    def WriteRestartFile(self,FileName):
        restart_file = open(FileName+".mdpa",'w')
        import new_restart_utilities
        new_restart_utilities.PrintProperties(restart_file)
        new_restart_utilities.PrintNodes(self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintElements("Fluid3D",self.model_part.Elements,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_X,"VELOCITY_X",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Y,"VELOCITY_Y",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VELOCITY_Z,"VELOCITY_Z",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(PRESSURE,"PRESSURE",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(VISCOSITY,"VISCOSITY",self.model_part.Nodes,restart_file)
        new_restart_utilities.PrintRestart_ScalarVariable(DENSITY,"DENSITY",self.model_part.Nodes,restart_file)
        restart_file.close() 




        
        

