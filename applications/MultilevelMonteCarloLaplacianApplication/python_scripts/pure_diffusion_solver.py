from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.MyMultilevelMonteCarloLaplacianApplication as Poisson
from python_solver import PythonSolver
# import sys
# import pprint

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(model, custom_settings):
    return PureDiffusionSolver(model, custom_settings)


class PureDiffusionSolver(PythonSolver):


    def __init__(self, model, custom_settings):  # Constructor of the class

        ## Overwrite the default settings with user-provided parameters
        settings = self._ValidateSettings(custom_settings)
        
        super(PureDiffusionSolver,self).__init__(model, settings)
        
        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True
        ## Set the element and condition names for the replace settings
        self.element_name = "MyLaplacianElement"
        # self.condition_name = "PointSourceCondition"
        self.condition_name = None
        
        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(Poisson.SOLUTION)
        self.main_model_part.AddNodalSolutionStepVariable(Poisson.FORCING)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)

        
    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(Poisson.SOLUTION,self.main_model_part)
        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("PureDiffusionSolver", "PureDiffusion solver DOFs added correctly.")

    
    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

        
    def PrepareModelPart(self):
        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            self.model.CreateModelPart(self.settings["model_part_name"].GetString())
            
        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("PureDiffusionSolverSolver", "Model reading finished.")

            
    def GetMinimumBufferSize(self):
        return 1
            
            
    def Initialize(self):
        # Creating the solution strategy
        self.conv_criteria = KratosMultiphysics.DisplacementCriteria(self.settings["relative_tolerance"].GetDouble(),
                                                                     self.settings["absolute_tolerance"].GetDouble())
        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())
        self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        # self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.main_model_part,
        #                                                              self.time_scheme,
        #                                                              self.linear_solver,
        #                                                              self.builder_and_solver,
        #                                                              self.settings["compute_reactions"].GetBool(),
        #                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
        #                                                              self.settings["calculate_norm_dx"].GetBool(),
        #                                                              self.settings["move_mesh_flag"].GetBool())

        self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.main_model_part,
                                                                     self.time_scheme,
                                                                     self.linear_solver,
                                                                     self.builder_and_solver,
                                                                     False,
                                                                     False,
                                                                     False,
                                                                     False)


        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())
        (self.solver).Initialize()
        (self.solver).Check()
        # print ("Pure diffusion solver initialization finished")

        
    def Solve(self):
        # Solve equations on mesh
        print("!!!WARNING: YOU SHOULD NOT USE solver.Solve()!!!")
        (self.solver).Solve()

        
    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.solver.SolveSolutionStep()
            if not is_converged and self._IsPrintingRank():
                msg  = "Pure Diffusion solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
                msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
                KratosMultiphysics.Logger.PrintWarning("PureDiffusionSolver",msg)

        
    def Finalize(self):
        # now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file
        # since in the .json file the output configuration was not defined, we define it now
        gid_mode = KratosMultiphysics.GiDPostMode.GiD_PostAscii  #we import the python file that includes the commands that we need
        multifile = KratosMultiphysics.MultiFileFlag.SingleFile
        deformed_mesh_flag = KratosMultiphysics.WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = KratosMultiphysics.WriteConditionsFlag.WriteElementsOnly
        gid_io = KratosMultiphysics.GidIO("MultilevelMonteCarloLaplacian", gid_mode,multifile, deformed_mesh_flag, write_conditions)

        # we create a mesh for the postprocess
        mesh_name = 0.0
        gid_io.InitializeMesh( mesh_name )
        gid_io.WriteMesh((self.main_model_part).GetMesh())
        gid_io.FinalizeMesh()

        # and we print the results
        # select which result to print on the nodes
        gid_io.InitializeResults(mesh_name,(self.main_model_part).GetMesh())
        gid_io.WriteNodalResults(Poisson.SOLUTION,self.main_model_part.Nodes,0,0)
        gid_io.FinalizeResults()


    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        return new_time

    
    def _ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

        return delta_time

    
    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.GetMinimumBufferSize()

    
    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.InitializeSolutionStep()
            

    def Predict(self):
        if self._TimeBufferIsInitialized():
            self.solver.Predict()


    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            (self.solver).FinalizeSolutionStep()

            
    def Check(self):
        (self.solver).Check()

        
    def Clear(self):
        (self.solver).Clear()

        
    def GetComputingModelPart(self):
        if not self.main_model_part.HasSubModelPart("Parts_Domain"):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.main_model_part.GetSubModelPart("Parts_Domain")


    def _ValidateSettings(self, settings):
        ## Settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
        "model_import_settings"              : {
            "input_type"     : "mdpa",
            "input_filename" : "../tests/FirstTestCoarserMesh/FirstTestCoarserMesh"
        },
	    "model_part_name" : "MLMCLaplacianModelPart",
	    "buffer_size"     : 2,
	    "domain_size"     : 2,
        "echo_level"                         : 1,
        "maximum_iterations"                 : 20,
        "relative_tolerance"                 : 1e-6,
        "absolute_tolerance"                 : 1e-9,
	    "compute_reactions"          	     : false,
	    "reform_dofs_at_each_step"    	     : false,
	    "calculate_norm_dx"           	     : true,
	    "move_mesh_flag"             	     : false,
        "problem_domain_sub_model_part_list" : ["Parts_Domain"],
        "processes_sub_model_part_list"      : ["Subpart_Boundary"],
	    "linear_solver_settings"       : {
			"solver_type"     : "SkylineLUFactorizationSolver"
	    },
        "time_stepping" : {
        "automatic_time_step" : false,
        "time_step" : 1.1
        }
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        return settings

    
    def _IsPrintingRank(self):
        return self._is_printing_rank


