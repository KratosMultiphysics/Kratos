#importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.LagrangianMPMApplication
#import KratosMultiphysics.SolidMechanicsApplication

def CreateSolver(main_model_part, custom_settings):
    return LagrangianMPMSolver(main_model_part, custom_settings)

class LagrangianMPMSolver:

    ##constructor. the constructor shall only take care of storing the settings
    ##and the pointer to the main_model part. This is needed since at the point of constructing the
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings):

        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "lagrangian_mpm_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "maximum_iterations": 10,
            "echo_level": 1,
            "compute_reactions": false,
            "reform_dofs_at_each_step": true,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-7,
            "linear_solver_settings"        : {
                "solver_type" : "Super LU",
                "scaling": true
            },
            "processes_sub_model_part_list" : [""]
        }""")

	##overwrite the default settings with user-provided parameters
	self.settings = custom_settings
	self.settings.ValidateAndAssignDefaults(default_settings)

	#self.EraseNodes = NodeEraseProcess(self.main_model_part)

	#construct the linear solver
	import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
	self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

	print("Construction of NavierStokesSolver_VMSMonolithic finished")


    def GetMinimumBufferSize(self):
	return 2;

    def AddVariables(self):
	self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
	self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
	self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
	self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
	self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.LAGRANGE_DISPLACEMENT)
	#self.main_model_part.AddNodalSolutionStepVariable(NODAL_STRESS);
	#self.main_model_part.AddNodalSolutionStepVariable(NODAL_AREA);
	#self.main_model_part.AddNodalSolutionStepVariable(NODAL_STRAIN);
	self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE);

	print("variables for the lagrangian MPM solver added correctly")

    def ImportModelPart(self):

	if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
	    #here it would be the place to import restart data if required

	    KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)



	    print("model part successfully read")

	    #generate the MPM particles

	    property_id = 1
	    KratosMultiphysics.LagrangianMPMApplication.CreateMLSParticleGauss(self.main_model_part, "UpdatedLagrangianMPMElement","Condition",property_id).Execute()

	else:
	    raise Exception("Other input options are not yet implemented.")

	current_buffer_size = self.main_model_part.GetBufferSize()
	if(self.GetMinimumBufferSize() > current_buffer_size):
	    self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

	print ("Model reading finished.")


    def AddDofs(self):

	for node in self.main_model_part.Nodes:
	    # adding dofs
	    node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X)
	    node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y)
	    node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z)
	    node.AddDof(KratosMultiphysics.VELOCITY_X)
	    node.AddDof(KratosMultiphysics.VELOCITY_Y)
	    node.AddDof(KratosMultiphysics.VELOCITY_Z)
	    node.AddDof(KratosMultiphysics.ACCELERATION_X)
	    node.AddDof(KratosMultiphysics.ACCELERATION_Y)
	    node.AddDof(KratosMultiphysics.ACCELERATION_Z)
	    node.AddDof(KratosMultiphysics.LAGRANGE_DISPLACEMENT_X)
	    node.AddDof(KratosMultiphysics.LAGRANGE_DISPLACEMENT_Y)


	print("DOFs for the VMS fluid solver added correctly.")








    def Initialize(self):

	self.computing_model_part = self.GetComputingModelPart()



	# creating the solution strategy
	self.conv_criteria = KratosMultiphysics.DisplacementCriteria(self.settings["relative_tolerance"].GetDouble(),
						    self.settings["absolute_tolerance"].GetDouble())
	print("echo level", self.settings["echo_level"].GetInt())
	self.conv_criteria.SetEchoLevel(self.settings["echo_level"].GetInt())
	print(self.conv_criteria)


	self.damp_factor_m = -0.01;
	#self.damp_factor_m = 0.0;
	self.dynamic_factor = 1.0;

	#self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
	self.time_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(self.damp_factor_m)
	#self.time_scheme = KratosMultiphysics.ResidualBasedNewmarkDisplacementScheme()


	builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
	#builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)

	MoveMeshFlag = False
	self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.main_model_part,
									    self.time_scheme,
									    self.linear_solver,
									    self.conv_criteria,
									    builder_and_solver,
									    self.settings["maximum_iterations"].GetInt(),
									    self.settings["compute_reactions"].GetBool(),
									    self.settings["reform_dofs_at_each_step"].GetBool(),
									    MoveMeshFlag)


	(self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
	(self.solver).Initialize()

	#self.solver.Check() #TODO: define a correct checking in the element
	#self.Check()
	print ("lagrangian MPM solver initialization finished.")

    def GetComputingModelPart(self):
	return self.main_model_part

    def GetOutputVariables(self):
	pass

    def ComputeDeltaTime(self):
	pass

    def SaveRestart(self):
	pass #one should write the restart file here

    def Solve(self):

	#KratosMultiphysics.LagrangianMPMApplication.NodeAndElementEraseProcess(self.main_model_part).Execute()


	#(self.PfemUtils).MarkOuterNodes(self.box_corner1, self.box_corner2, (self.model_part).Nodes);
	#(self.EraseNodes).Execute();

	#for node in self.main_model_part.Nodes:
	    #node.X = node.X0 + node.GetSolutionStepValue(DISPLACEMENT_X)
	    #node.Y = node.Y0 + node.GetSolutionStepValue(DISPLACEMENT_Y)

	#recompute the cloud of nodes - this refle
	KratosMultiphysics.LagrangianMPMApplication.RecomputeNeighboursProcess(self.main_model_part).Execute()
	#KratosMultiphysics.LagrangianMPMApplication.Gauss_Coordinates_Update_Process(self.main_model_part).Execute()

	#for elem in self.main_model_part.Elements:
	    #for node in elem.GetNodes():
		#print(node.Id)
	    #print(" ")            for node in self.main_model_part.Nodes:


	#compute the shape functions on the geometry. Note that the node position is still the one at the end of the previous step, so it is 0!!
	KratosMultiphysics.LagrangianMPMApplication.ComputeMLSShapeFunctionsProcess(self.main_model_part).Execute()

	self._CheckShapeFunctions()
	#KratosMultiphysics.LagrangianMPMApplication.NodeAndElementEraseProcess(self.main_model_part).Execute()

	#ensure that the system graph is to be reformed
	self.solver.Clear()


	self.solver.Solve()




    def SetEchoLevel(self, level):
	self.solver.SetEchoLevel(level)

    def Clear(self):
	self.solver.Clear()

    def Check(self):
	self.solver.Check()

    def _CheckShapeFunctions(self):
	for elem in self.main_model_part.Elements:
	    nnodes = len(elem.GetNodes())
	    DN = elem.GetValue(KratosMultiphysics.LagrangianMPMApplication.SHAPE_FUNCTIONS_DERIVATIVES)
	    N = elem.GetValue(KratosMultiphysics.LagrangianMPMApplication.SHAPE_FUNCTIONS)
	    dnx = 0.0
	    dny = 0.0
	    N_tot = 0.0

	    for i in range(len(elem.GetNodes())):
		N_tot += N[i]
		dnx += DN[i,0]
		dny += DN[i,1]
	    if(abs(N_tot - 1.0) > 1e-6):
		print(" N_tot ", N_tot)
	    if(abs(dnx + dny) > 1e-6):
		print(" DN sum  ", dnx+dny)


