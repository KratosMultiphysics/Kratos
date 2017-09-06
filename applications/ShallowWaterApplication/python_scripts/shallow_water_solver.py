from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication as KratosShallow

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(model_part, custom_settings):
    
    return ShallowWaterSolver(model_part, custom_settings)

class ShallowWaterSolver:
    def __init__(self,model_part, custom_settings):  # Constructor of the class 
        self.model_part = model_part
        self.settings = custom_settings
        
        #~ self.domain_size = self.settings["domain_size"].GetInt()
        self.domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        #~ self.domain_size = self.model_part.ProcessInfo[DOMAIN_SIZE]
        self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        # Definition of the mesh stage solver
        self.swe_linear_solver =  KratosMultiphysics.SkylineLUFactorizationSolver()  # We set the type of solver that we want 
        
        # definition of the convergence criteria
        self.conv_criteria = KratosMultiphysics.DisplacementCriteria(1e-6,1e-9)  # Tolerance for the solver 
        
        # For the pfem2
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)
        (self.neighbour_search).Execute()
        self.neighbour_elements_search= KratosMultiphysics.FindElementalNeighboursProcess(model_part,self.domain_size,number_of_avg_elems)
        (self.neighbour_elements_search).Execute()

    def AddVariables(self):  #this way er only need one command to add all the variables to our problem 
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.PROJECTED_VELOCITY);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.DELTA_VELOCITY);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.HEIGHT);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.PROJECTED_HEIGHT);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.DELTA_HEIGHT);
        # Physic problem parameters
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.BATHYMETRY);
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.GRAVITY);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.RAIN);
        # Specific variables to convect particles
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YP);
        self.model_part.AddNodalSolutionStepVariable(KratosShallow.MEAN_SIZE);

    def AddDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(KratosMultiphysics.VELOCITY_X);
            node.AddDof(KratosMultiphysics.VELOCITY_Y);
            node.AddDof(KratosShallow.HEIGHT);
        print ("variables for the SWE solver added correctly")

    def GetMinimumBufferSize(self):
        return 1

    def GetComputingModelPart(self):
        return self.model_part.GetSubModelPart("Parts_water")

    def ImportModelPart(self):
        ## Read model part
        self._ModelPartReading()
        ## Replace elements and conditions, check the input reading and set KINEMATIC_VISCOSITY and DENSITY
        self._ExecuteAfterReading()
        ## Set buffer size
        self._SetBufferSize()

    def Initialize(self):
        # Creating the solution strategy for the mesh stage
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = False
        import strategy_python
        self.mesh_solver = strategy_python.SolvingStrategyPython(self.model_part,
                                            self.time_scheme, self.swe_linear_solver,
                                            self.conv_criteria, CalculateReactionFlag,
                                            ReformDofSetAtEachStep, MoveMeshFlag)

        # Creating the solution strategy for the particle stage
        self.VariableUtils = KratosMultiphysics.VariableUtils()
        maximum_number_of_particles= 8*self.domain_size
        self.moveparticles = KratosShallow.MoveShallowWaterParticleUtility(self.model_part,maximum_number_of_particles)  
        self.moveparticles.MountBin()
        
        # Initialize dry/wet state utility
        #~ self.drybedutility = DryBedUtility(self.model_part)

    def Solve(self):
        # Check dry bed
        #~ (self.drybedutility).CheckConservedVariables(self.model_part.Nodes)
        
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
        
        # Solve equations on mesh
        (self.mesh_solver).Solve()
        
        # Update particles
        (self.moveparticles).CalculateDeltaVariables();        
        (self.moveparticles).CorrectParticlesWithoutMovingUsingDeltaVariables();
        post_minimum_number_of_particles=self.domain_size*2;
        (self.moveparticles).PostReseed(post_minimum_number_of_particles); 

    def SetEchoLevel(self,level):
        (self.mesh_solver).SetEchoLevel(level)

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

