#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.FreezingSoilApplication import *
CheckForPreviousImport()

import sys

def AddVariables(model_part):	
	########################################	
	# from KratosMengmengApplication	
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_NULL);
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_EINS);
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_DT);
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_EINS_DT);
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_NULL_DT);  
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_ACCELERATION);
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_NULL_ACCELERATION);
	model_part.AddNodalSolutionStepVariable(TEMPERATURE_EINS_ACCELERATION);		
					
	model_part.AddNodalSolutionStepVariable(FACE_WATER_FLUX);
	model_part.AddNodalSolutionStepVariable(FACE_LOAD_PRESSURE);
		
	model_part.AddNodalSolutionStepVariable(LINEAR_STRAIN); 
	model_part.AddNodalSolutionStepVariable(EFFECTIVE_STRESS);  
	model_part.AddNodalSolutionStepVariable(TOTAL_STRESS);  
	model_part.AddNodalSolutionStepVariable(SUCTION); 

	model_part.AddNodalSolutionStepVariable(KRATOS_WATCH_FLAG); 
	model_part.AddNodalSolutionStepVariable(ASSIGN_PRESTRESS_FLAG);  
	model_part.AddNodalSolutionStepVariable(PLASTIC_FLAG);  
	
	model_part.AddNodalSolutionStepVariable(PRESTRESS);

	model_part.AddNodalSolutionStepVariable(ICE_MASS); 
	model_part.AddNodalSolutionStepVariable(WATER_MASS); 
	model_part.AddNodalSolutionStepVariable(ICE_PRESSURE); 
	model_part.AddNodalSolutionStepVariable(ICE_SATURATION); 
	model_part.AddNodalSolutionStepVariable(ICE_VOLUME_FRACTION);	
		
	model_part.AddNodalSolutionStepVariable(WATER_FLOW); 
	model_part.AddNodalSolutionStepVariable(ICE_FLOW); 
	model_part.AddNodalSolutionStepVariable(HEAT_FLOW); 
	
	model_part.AddNodalSolutionStepVariable(STRESS_TENSOR); 
	model_part.AddNodalSolutionStepVariable(STRAIN_TENSOR); 
	
	model_part.AddNodalSolutionStepVariable(SCALE_U); 
	model_part.AddNodalSolutionStepVariable(SCALE_O); 
	
	model_part.AddNodalSolutionStepVariable(MECH_DISSIPATION); 
	
	model_part.AddNodalSolutionStepVariable(ICE_DENSITY); 
	model_part.AddNodalSolutionStepVariable(WATER_DENSITY); 
	
	model_part.AddNodalSolutionStepVariable(PLASTICITY_INDICATOR); 
	model_part.AddNodalSolutionStepVariable(INSITU_STRESS_SCALE); 
	model_part.AddNodalSolutionStepVariable(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR); 
	
	model_part.AddNodalSolutionStepVariable(PRECONSOLIDATION); 
	model_part.AddNodalSolutionStepVariable(EQUIVALENT_VOLUMETRIC_STRAIN); 
	model_part.AddNodalSolutionStepVariable(EQUIVALENT_DEVIATORIC_STRAIN); 
	model_part.AddNodalSolutionStepVariable(EQUIVALENT_VOLUMETRIC_STRESS); 
	model_part.AddNodalSolutionStepVariable(EQUIVALENT_DEVIATORIC_STRESS);  
	model_part.AddNodalSolutionStepVariable(LOG_EQUIVALENT_VOLUMETRIC_STRESS);  
	
	model_part.AddNodalSolutionStepVariable(ELEMENT_PARAMETERS);  
	
	
	########################################
	# from KRATOS
	model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
	model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL);
	model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS);
	model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT);
	model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT);
	model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT);  
	model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);  
	model_part.AddNodalSolutionStepVariable(ACCELERATION);
	model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL);
	model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS);
			
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);	
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL);
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS);		
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_DT);		
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL_DT);		
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS_DT);	
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_ACCELERATION);	
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_NULL_ACCELERATION);
	model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_EINS_ACCELERATION);
	
	model_part.AddNodalSolutionStepVariable(TEMPERATURE);
	model_part.AddNodalSolutionStepVariable(POROSITY);
	
	model_part.AddNodalSolutionStepVariable(DENSITY_WATER);   
	model_part.AddNodalSolutionStepVariable(DENSITY);   
	model_part.AddNodalSolutionStepVariable(VISCOSITY); 
	model_part.AddNodalSolutionStepVariable(DISTANCE); 
	
	model_part.AddNodalSolutionStepVariable(SATURATION); 
	model_part.AddNodalSolutionStepVariable(FACE_LOAD);     
	
	model_part.AddNodalSolutionStepVariable(YIELD_STRESS);  
	    
	model_part.AddNodalSolutionStepVariable(PRESSURE);     
	
	model_part.AddNodalSolutionStepVariable(ERROR_RATIO); 
	model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE); 
	
	model_part.AddNodalSolutionStepVariable(FACE_HEAT_FLUX); 
	
	model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE); 
	model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE); 
	
	model_part.AddNodalSolutionStepVariable(FORCE); 
	model_part.AddNodalSolutionStepVariable(ROTATION); 
				
	model_part.AddNodalSolutionStepVariable(AIR_PRESSURE);	
	model_part.AddNodalSolutionStepVariable(FRICTION_COEFFICIENT);	
	model_part.AddNodalSolutionStepVariable(BULK_MODULUS);	
	#model_part.AddNodalSolutionStepVariable(SHEAR_MODULUS);	
	
	model_part.AddNodalSolutionStepVariable(CONVECTION_COEFFICIENT);  
	model_part.AddNodalSolutionStepVariable(AMBIENT_TEMPERATURE);  
	model_part.AddNodalSolutionStepVariable(EMISSIVITY);  
	
	
	model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT);  
	
	model_part.AddNodalSolutionStepVariable(SCALE);  
	
print "variables for the dynamic structural solution added correctly"
	
def AddDofs(model_part):
	for node in model_part.Nodes:
		#adding dofs
		node.AddDof(DISPLACEMENT_X);
		node.AddDof(DISPLACEMENT_Y);
		node.AddDof(DISPLACEMENT_Z);
		node.AddDof(WATER_PRESSURE);
		node.AddDof(TEMPERATURE);
	print "dofs for the dynamic structural solution added correctly"


class MengmengSolver:
        #######################################################################
        
        
	def __init__(self,model_part, domain_size, abs_tol, rel_tol):
	
		self.model_part = model_part
		self.echo_level = 0
		
		self.damp_factor = 1.0
		#self.toll = 1e-12
		#self.absolute_tol = 1e-14
		self.toll = rel_tol
		self.absolute_tol = abs_tol
	
		#definition of the solvers
		self.structure_linear_solver =  SkylineLUFactorizationSolver()
		
		#definition of the convergence criteria
		#self.conv_criteria = MultiPhaseFlowCriteria(1e-6, 1e-6)
	
	#######################################################################
	def Initialize(self):
	
		self.time_scheme = GeneralizedAlphaMengmeng(self.damp_factor)
	
		#definition of the convergence criteria
		self.conv_criteria = MengmengCriteria(self.toll,self.absolute_tol)
		#definition of BuilderAndSolver
		#builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
	
		#creating the solution strategy
		self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,30,False,True,False)
		print "self.echo_level = " , self.echo_level
		(self.solver).SetEchoLevel(self.echo_level)
	
		print "finished initialization of the mengmeng strategy"
	
			
	#######################################################################   
	def Solve(self):
		(self.solver).Solve()
		
	#######################################################################   
	def SetEchoLevel(self,level):
		(self.solver).SetEchoLevel(level)
	
