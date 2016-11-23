from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# MPI
import KratosMultiphysics.mpi as mpi
import KratosMultiphysics.TrilinosApplication as TrilinosApplication
import KratosMultiphysics.MetisApplication as MetisApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import solid_mechanics_solver

def CreateSolver(main_model_part, custom_settings):
    return MechanicalSolver(main_model_part, custom_settings)

class MechanicalSolver(solid_mechanics_solver.MechanicalSolver):
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "structural_mechanics_solver_MPI",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "analysis_type": "Non-Linear",
            "time_integration_method": "Implicit",
            "scheme_type": "Newmark",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": null,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "implex": false,
            "compute_reactions": false,
            "compute_contact_forces": false,
            "block_builder": true,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type" : "Klu",
                "scaling": false
            },
            "bodies_list": [],
            "problem_domain_sub_model_part_list": ["solid"],
            "processes_sub_model_part_list": [""]
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solver
        import new_trilinos_linear_solver_factory # TODO: Is new_trilinos_linear_solver_factory or trilinos_linear_solver_factory?
        self.linear_solver = new_trilinos_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        # Creating the partitions
        self.CreatingPartitions()
        
    def AddVariables(self):
        solid_mechanics_solver.MechanicalSolver.AddVariables(self)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
        
    def CreatingPartitions(self):
        # Creating the partitions
        self.domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
                
        input_file_name = self.settings["model_import_settings"]["input_filename"].GetString()
        
        number_of_partitions = mpi.mpi.size  # we set it equal to the number of processors
        
        if mpi.mpi.size > 1:
            model_part_io = KratosMultiphysics.ModelPartIO(input_file_name)
            if mpi.mpi.rank == 0:
                verbosity = 1
                sync_conditions = True # Make sure that the condition goes to the same partition as the element is is a face of
                partitioner = MetisApplication.MetisDivideHeterogeneousInputProcess( model_part_io, number_of_partitions, self.domain_size, verbosity, sync_conditions)
                partitioner.Execute()

            mpi.mpi.world.barrier()

            self.settings["model_import_settings"]["input_filename"].SetString(input_file_name + "_" + str(mpi.mpi.rank))

            MPICommSetup = MetisApplication.SetMPICommunicatorProcess(self.main_model_part)
            MPICommSetup.Execute()

            self.Comm = TrilinosApplication.CreateCommunicator()
        else:
            raise NameError('Your number of mpi.size is 1')
        
        print("Construction of MPI MechanicalSolver finished")
    
    def ImportModelPart(self):

        if(self.domain_size == 2): # TODO: Add to the input
            self.guess_row_size = 15
        else:
            self.guess_row_size = 45

        solid_mechanics_solver.MechanicalSolver.ImportModelPart(self)
        
    def _GetConvergenceCriterion(self): 
        # Creation of an auxiliar Kratos parameters object to store the convergence settings
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("component_wise",self.settings["component_wise"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])
        
        # Construction of the class convergence_criterion
        import convergence_criteria_factory_MPI
        convergence_criterion = convergence_criteria_factory_MPI.convergence_criterion(conv_params)
        
        return convergence_criterion.mechanical_convergence_criterion
    
    def _GetBuilderAndSolver(self, component_wise, block_builder):
        if(block_builder):
            # To keep matrix blocks in builder
            builder_and_solver = "standard"
        else:
            builder_and_solver = "residual"
    
        return builder_and_solver
    
            
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search, implex):
        import trilinos_strategy_python
        self.mechanical_solver = trilinos_strategy_python.SolvingStrategyPython(builder_and_solver, self.main_model_part, mechanical_scheme, self.linear_solver, mechanical_convergence_criterion, compute_reactions, reform_step_dofs, move_mesh_flag, self.Comm, self.guess_row_size)