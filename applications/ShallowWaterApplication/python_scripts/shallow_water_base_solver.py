from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(model_part, custom_settings):
    return ShallowWaterBaseSolver(model_part, custom_settings)

class ShallowWaterBaseSolver(object):
    def __init__(self,model_part, custom_settings):  # Constructor of the class 
        self.model_part  = model_part
        self.domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type"                  : "shallow_water_base_solver",
            "model_import_settings"        : {
                "input_type"          : "mdpa",
                "input_filename"      : "unknown_name"
            },
            "echo_level"                   : 0,
            "convergence_echo_level"       : 1,
            "solver_echo_level"            : 0,
            "buffer_size"                  : 2,
            "relative_tolerance"           : 1e-6,
            "absolute_tolerance"           : 1e-9,
            "maximum_iterations"           : 20,
            "compute_reactions"            : false,
            "reform_dofs_at_each_step"     : false,
            "move_mesh_flag"               : false,
            "volume_model_part_name"       : "volume_model_part",
            "linear_solver_settings"       : {
                    "solver_type"     : "SkylineLUFactorizationSolver"
            },
            "time_stepping"                : {
                "automatic_time_step" : false,
                "time_step"           : 0.01
            }
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def AddVariables(self):
        # Basic variables
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.HEIGHT);
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY);
        # Physic problem parameters
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.GRAVITY);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.BATHYMETRY);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.RAIN);
        # Auxiliar variables
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)

    def AddDofs(self):
        pass

    def GetMinimumBufferSize(self):
        return 1

    def GetComputingModelPart(self):
        return self.model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString())

    def ImportModelPart(self):
        ## Read model part
        self._ModelPartReading()
        ## Replace elements and conditions, check the input reading
        self._ExecuteAfterReading()
        ## Set buffer size
        self._SetBufferSize()

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # If needed, create the estimate time step utility
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            self.EstimateDeltaTimeUtility = KratosFluid.EstimateDtUtility2D(self.model_part,
                                                                            self.settings["time_stepping"])

        # Creating the solution strategy for the mesh stage
        self.conv_criteria = KratosMultiphysics.DisplacementCriteria(self.settings["relative_tolerance"].GetDouble(),
                                                                     self.settings["absolute_tolerance"].GetDouble())
        (self.conv_criteria).SetEchoLevel(self.settings["convergence_echo_level"].GetInt())

        #~ self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticSchemeSlip(self.domain_size,   # DomainSize
                                                                                             self.domain_size+1) # BlockSize

        builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.model_part,
                                                                            self.time_scheme,
                                                                            self.linear_solver,
                                                                            self.conv_criteria,
                                                                            builder_and_solver,
                                                                            self.settings["maximum_iterations"].GetInt(),
                                                                            self.settings["compute_reactions"].GetBool(),
                                                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                            self.settings["move_mesh_flag"].GetBool())

        (self.solver).SetEchoLevel(self.settings["solver_echo_level"].GetInt())
        (self.solver).Check()

        (self.solver).Initialize()

        print ("Mesh stage solver initialization finished.")

        # Initialize dry/wet state utility
        #~ self.drybedutility = DryBedUtility(self.model_part)

    def Solve(self):
        # Solve equations on mesh
        (self.solver).Solve()

    def Clear(self):
        (self.solver).Clear()

    def ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        return delta_time

    #### Specific internal functions ####

    def _ModelPartReading(self):
        ## Model part reading
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            ## Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.model_part)
        else:
            raise Exception("Other input options are not yet implemented.")

    def _ExecuteAfterReading(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])

    def _SetBufferSize(self):
        ## Set the buffer size
        self.model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        if(self.GetMinimumBufferSize() > self.model_part.GetBufferSize() ):
            self.model_part.SetBufferSize(self.GetMinimumBufferSize())

    def _AddParticleVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.PROJECTED_HEIGHT);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.DELTA_HEIGHT);
        # Specific variables to convect particles
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YP);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.MEAN_SIZE);

    def _AddPrimitiveDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(KratosMultiphysics.VELOCITY_X);
            node.AddDof(KratosMultiphysics.VELOCITY_Y);
            node.AddDof(KratosShallow.HEIGHT);
        print ("variables for the SWE solver added correctly")

    def _AddConservedDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(KratosMultiphysics.MOMENTUM_X);
            node.AddDof(KratosMultiphysics.MOMENTUM_Y);
            node.AddDof(KratosShallow.HEIGHT);
        print ("variables for the SWE solver added correctly")

    def _pfem2_init(self,model_part):
        # For the pfem2
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        (self.neighbour_search).Execute()
        self.neighbour_elements_search= KratosMultiphysics.FindElementalNeighboursProcess(model_part,self.domain_size,number_of_avg_elems)
        (self.neighbour_elements_search).Execute()

    def ExecuteParticlesUtilitiesBeforeSolve(self):
        # Move particles
        (self.moveparticles).CalculateVelOverElemSize();
        (self.moveparticles).MoveParticles();
        pre_minimum_number_of_particles=self.domain_size;
        (self.moveparticles).PreReseed(pre_minimum_number_of_particles);    
        (self.moveparticles).TransferLagrangianToEulerian();
        (self.VariableUtils).CopyScalarVar(KratosShallow.PROJECTED_HEIGHT,KratosShallow.HEIGHT,self.model_part.Nodes)
        (self.VariableUtils).CopyVectorVar(KratosShallow.PROJECTED_VELOCITY,KratosMultiphysics.VELOCITY,self.model_part.Nodes)
        (self.moveparticles).ResetBoundaryConditions()
        #(self.moveparticles).CopyScalarVarToPreviousTimeStep(KratosShallow.PROJECTED_HEIGHT,self.model_part.Nodes)
        #(self.moveparticles).CopyVectorVarToPreviousTimeStep(KratosShallow.PROJECTED_VELOCITY,self.model_part.Nodes)

    def ExecuteParticlesUtilitiesAfetrSolve(self):
        # Update particles
        (self.moveparticles).CalculateDeltaVariables();
        (self.moveparticles).CorrectParticlesWithoutMovingUsingDeltaVariables();
        post_minimum_number_of_particles=self.domain_size*2;
        (self.moveparticles).PostReseed(post_minimum_number_of_particles);
