from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
#from KratosMultiphysics import *

import KratosMultiphysics.IGAStructuralMechanicsApplication # import *

#importing the Kratos Library
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication 

import ImportModelPart


# check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(model_part, custom_settings):
	return IGAStructuralMechanicsSolver(model_part, custom_settings)

class IGAStructuralMechanicsSolver:
	def __init__(self, model_part, custom_settings):
		#TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
		self.model_part = model_part 

		##settings string in json format
		default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "iga_structural_mechanics_solver",
            "model_import_settings": {
                "input_type": "txt",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "echo_level": 1,
            "compute_reactions": false,
            "reform_dofs_at_each_iteration": true,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-7,
            "linear_solver_settings"        : {
				"solver_type"   : "BiConjugate_gradient_stabilized",
				"max_iteration" : 500,
				"tolerance"     : 1e-9,
				"scaling"       : true,
				"verbosity"     : 0
            },
			"problem_domain_sub_model_part_list" : [],
			"processes_sub_model_part_list"      : [],
			"rotation_dofs"                      : "false"
        }""")   

        ##overwrite the default settings with user-provided parameters
		self.settings = custom_settings
		self.settings.ValidateAndAssignDefaults(default_settings)

        #construct the linear solver
		import linear_solver_factory
		self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])


	def GetMinimumBufferSize(self):
		return 3;

	def AddVariables(self):
		# Add displacements
		self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER)
		## Add dynamic variables
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
		## Add reactions for the displacements
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
		## Add nodal force variables
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FORCE)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCE)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
		## Add specific variables for the problem conditions
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SolidMechanicsApplication.LINE_LOAD)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SolidMechanicsApplication.SURFACE_LOAD)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IGAStructuralMechanicsApplication.INTEGRATION_WEIGHT)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_VALUES)
		#self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IGAStructuralMechanicsApplication.SHAPE_FUNCTION_LOCAL_DERIVATIVES)

		#t = True

		#if (t):#self.settings["rotation_dofs"].GetBool():
		#	# Add specific variables for the problem (rotation dofs)
		#	self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
		#	self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
		#	self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.StructuralMechanicsApplication.POINT_TORQUE)
		#	self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
		#	self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

		self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER)
		print("Added Variables: DISPLACEMENT and VECTOR_LAGRANGE_MULTIPLIER")

		
	def ImportModelPart(self, projectparameters):
		#print ("\n------------------Import Elements--------------------")
		
		if(self.settings["model_import_settings"]["input_type"].GetString() == "txt"):
            #here it would be the place to import restart data if required
			ModelPartIO = ImportModelPart.Factory(self.settings["model_import_settings"]["input_filename"].GetString(), self.settings, projectparameters)
			ModelPartIO.ReadModelPart(self.model_part)

		else:
			raise Exception("Other input options are not implemented.")

		current_buffer_size = self.model_part.GetBufferSize()
		if(self.GetMinimumBufferSize() > current_buffer_size):
			self.model_part.SetBufferSize( self.GetMinimumBufferSize() )

		print ("Model reading finished.")
	
	def ImportModelPartNurbsBrep(self, model_part_nurbs_brep):
		#print ("\n------------------Import Elements--------------------")
		
		if(self.settings["model_import_settings"]["input_type"].GetString() == "txt"):
            #here it would be the place to import restart data if required
			ModelPartIO = ImportModelPart.Factory(model_part_nurbs_brep, self.settings)
			ModelPartIO.ReadModelPart(self.model_part)

		else:
			raise Exception("Other input options are not implemented.")

		current_buffer_size = self.model_part.GetBufferSize()
		if(self.GetMinimumBufferSize() > current_buffer_size):
			self.model_part.SetBufferSize( self.GetMinimumBufferSize() )

		print ("Model reading finished.")
	
		
	def AddDofs(self):
		for node in self.model_part.Nodes:
			# adding dofs
			print(node.Id)
			node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X);
			node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y);
			node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z);

			node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_X);
			node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Y);
			node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Z);
		print("Added DOFs: DISPLACEMENT and VECTOR_LAGRANGE_MULTIPLIER")

	def Initialize(self):
		self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

		builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver) 

		# definition of the type of the solver
		self.laplace_solver = KratosMultiphysics.SkylineLUFactorizationSolver() 

		# definition of the convergence criteria (tolerance)
		self.conv_criteria = KratosMultiphysics.DisplacementCriteria(self.settings["relative_tolerance"].GetDouble(), self.settings["absolute_tolerance"].GetDouble())
		CalculateReactionFlag = False
		MoveMeshFlag = True
		self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria, builder_and_solver,  self.settings["maximum_iterations"].GetInt(), self.settings["compute_reactions"].GetBool(), self.settings["reform_dofs_at_each_iteration"].GetBool(), MoveMeshFlag)
		print ("Initialization finished\n")

	def Solve(self):
		#self.solver.Clear() 
		self.solver.Solve()

		
	def SetEchoLevel(self, level):
		(self.solver).SetEchoLevel(level)

	def GetComputeModelPart(self):
		return self.model_part

	def Clear(self):
		self.solver.Clear()
	
	def Check(self):
		self.solver.Check()
