---
title: Creating the Python Solver file
keywords: 
tags: [Tutorial-Creating-the-Python-Solver-file.md]
sidebar: kratos_for_developers
summary: 
---

The purpose of this file is to act as an interface between your problem's script (also in python) file and the Kratos functions. We will create a class and functions already customized for our type of problem, so that the calls required in the .py from the problem become simpler. The file should be in the `/MyLaplacianApplication/python_scripts` folder.

## pure_diffusion_solver.py  

```python

# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MyLaplacianApplication as Poisson

def CreateSolver(Model, custom_settings):
    return PureDiffusionSolver(Model, custom_settings)

class PureDiffusionSolver(object):
    def __init__(self, Model, custom_settings):  # Constructor of the class
        ## Settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name"              : "Main",
            "model_import_settings"        : {
                "input_type"          : "mdpa",
                "input_filename"      : "unknown_name"
            },
            "echo_level"                   : 0,
            "buffer_size"                  : 2,
            "relative_tolerance"           : 1e-6,
            "absolute_tolerance"           : 1e-9,
            "maximum_iterations"           : 20,
            "compute_reactions"            : false,
            "reform_dofs_at_each_step"     : false,
            "calculate_norm_dx"            : true,
            "move_mesh_flag"               : false,
            "linear_solver_settings"       : {
                    "solver_type"     : "skyline_lu_factorization"
            }
        }""")

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # We get the model part
        self.model_part  = Model[self.settings["model_part_name"].GetString()]
        self.domain_size = self.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## Construct the linear solver
        import KratosMultiphysics.python_linear_solver_factory as python_linear_solver_factory
        self.linear_solver = python_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE);
        self.model_part.AddNodalSolutionStepVariable(Poisson.POINT_HEAT_SOURCE);

    def AddDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(KratosMultiphysics.TEMPERATURE);
        
        print ("variables for the Poisson solver added correctly")

    def GetMinimumBufferSize(self):
        return 1

    def ImportModelPart(self):
        ## Model part reading
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            ## Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.model_part)
        else:
            raise Exception("Other input options are not yet implemented")
        ## Set the buffer size
        self.model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        if(self.GetMinimumBufferSize() > self.model_part.GetBufferSize() ):
            self.model_part.SetBufferSize(self.GetMinimumBufferSize())

    def Initialize(self):

        # Creating the solution strategy
        self.conv_criteria = KratosMultiphysics.DisplacementCriteria(self.settings["relative_tolerance"].GetDouble(),
                                                                     self.settings["absolute_tolerance"].GetDouble())
        (self.conv_criteria).SetEchoLevel(self.settings["echo_level"].GetInt())

        self.time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

        self.builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.model_part,
                                                                     self.time_scheme,
                                                                     self.linear_solver,
                                                                     self.builder_and_solver,
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["calculate_norm_dx"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())


        (self.solver).SetEchoLevel(self.settings["echo_level"].GetInt())

        (self.solver).Initialize()

        print ("Pure diffusion solver initialization finished")

    def Solve(self):
        # Solve equations on mesh
        (self.solver).Solve()
```
