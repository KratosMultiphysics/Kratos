---
title: How to create unitary tests
keywords: 
tags: [How-to-create-unitary-tests.md]
sidebar: kratos_for_developers
summary: 
---

## Overview

_Kratos Multiphysics_ has a mechanism to automatically test your code called `KratosUnittest`. If you are familiar with Python, this module is basically an extension of [Unittest](https://docs.python.org/2/library/unittest.html)  and you can expect to find every functionality that exists there in `KratosUnittest` as well.

In order for your application to be robust, it is recommended to add Unittests to it. Here we present some guidelines on how to do it. 

## Python

### How to run tests 

Kratos unittest are executed using the `run_tests.py` script located in the `kratos/kratos/python_scripts` folder.

Usage is the following: 

```console
python run_tests.py
```

you can specify the following options: 

```console
-l,--level:        Select the suit. Values: "All", "Nightly", "Small"(default)
-v,--verbosity:    Select the verbosity level of the output. Values: 0, 1 (default) , 2
-a,--applications: List of applications to run separated by ":". For example "-a KratosKore:IncompressibleFluidApplication"
                   All applications compiled are run by default.
```

### Basic structure

Tests are defined in a python script. To keep tests organized we recommend you to create the tests in your "application/tests/" directory.

Tests are organized in suites. Suites are a collection of tests that will be run together as a package. In _Kratos_, we define three basic suites:

* `All`: which should contain all the tests
* `Nightly`: which should contain a set of tests that could be executed in less than 10 min
* `Small`: which should contain a set of tests that could be executed in less than 1 min 


All applications should implement those packages as they can be automatically run, for example in the `Nightly` runs. In order to add tests you should create at least a couple of files, one to define the tests and one to define the suites: 

```console
"application/tests/test_NAME_OF_OUR_APPLICATION.py"
```

for example: 

```console
kratos/applications/example_application/tests/test_example_application.py
```

### General structure

This file will define the suites to run the tests. We define three different levels of suites in kratos:

* **All**: which should contain all the tests.
* **Nightly**: which should contain a set of tests that could be executed in less than 10 min.
* **Small**: which should contain a set of tests that could be executed in less than 1 min.

In order to add a test to some of these suites, one can do it as shown in this example: 

```python
from __future__ import print_function, absolute_import, division

# import Kratos
from KratosMultiphysics import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suites. For example
from test_my_app_example_tests_1 import TestCase1 as TestCase1
from test_my_app_example_tests_2 import TestCase2 as TestCase2

def AssembleTestSuites():
    ''' 
    Populates the test suites to run.
    
    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"
    
    Return
    ------
    
    suites: A dictionary of suites
        The set of suites with its test_cases added.
    ''' 
    
    # Get the already defined suites
    suites = KratosUnittest.KratosSuites
    
    # Get the small suite and populate it with some tests from TestCase1 and TestCase2
    smallSuite = suites['small']
    smallSuite.addTest(TestCase1('test_example_small_boo_1'))
    smallSuite.addTest(TestCase2('test_example_small_foo_1'))
    
    # Get the small suite and populate it with some tests from TestCase1 and TestCase2
    nightSuite = suites['nightly']
    smallSuite.addTest(TestCase1('test_example_nightly_boo_1'))
    smallSuite.addTest(TestCase1('test_example_nightly_boo_2'))
    smallSuite.addTest(TestCase2('test_example_nightly_foo_1'))
    
    # Get the small suite and populate it with all tests from TestCase1 and TestCase2
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TestCase1,
            TestCase2
        ])
    )
    
    # Return the suites
    return suites
    
# The main function executes the tests
if __name__ == '__main__':
    KratosUnittest.runTests(AssembleTestSuites())
```

### Structure for examples

Examples are cases that test your code/application, but are more complex than the unitary tests, take longer and involve many functionalities. A typical example is a Finite Element calculation of few time steps and a few elements. In order to compute some examples (the examples that used to be located in the folder `test_examples`) some additional considerations must to be taken into account. These considerations can be observed in the following example that is located in the StructuralMechanicsApplication/tests. In order to compute a dynamic test (_Bossak_ and _Newmark_) and a patch test (with bending and membrane behaviours), the following files has been considered, first the main file, which will be the one launched. 

```python
# import Kratos
import KratosMultiphysics 
import KratosMultiphysics.ExternalSolversApplication 
import KratosMultiphysics.SolidMechanicsApplication 
import KratosMultiphysics.StructuralMechanicsApplication 

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests o test_classes to create the suits
## SMALL TESTS
from SmallTests import SimpleMeshMovingTest as TSimpleMeshMovingTest
from SmallTests import DynamicBossakTests as TDynamicBossakTests
from SmallTests import DynamicNewmarkTests as TDynamicNewmarkTests
from SmallTests import SprismMembranePatchTests as TSprismMembranePatchTests
from SmallTests import SprismBendingPatchTests as TSprismBendingPatchTests
from SmallTests import ShellQ4ThickBendingRollUpTests as TShellQ4ThickBendingRollUpTests
from SmallTests import ShellQ4ThickDrillingRollUpTests as TShellQ4ThickDrillingRollUpTests
from SmallTests import ShellT3ThinBendingRollUpTests as TShellT3ThinBendingRollUpTests
from SmallTests import ShellT3ThinDrillingRollUpTests as TShellT3ThinDrillingRollUpTests
from SmallTests import EigenQ4Thick2x2PlateTests as TEigenQ4Thick2x2PlateTests
from SmallTests import EigenTL3D8NCubeTests as TEigenTL3D8NCubeTests

## NIGTHLY TESTS
from NightlyTests import ShellT3IsotropicScordelisTests as TShellT3IsotropicScordelisTests

## VALIDATION TESTS
from ValidationTests import SprismPanTests as TSprismPanTests

def AssambleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should pupulate the suites:
    "small", "nighlty" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''
    suites = KratosUnittest.KratosSuites

    # Create a test suit with the selected tests (Small tests):
    smallSuite = suites['small']
    smallSuite.addTest(TSimpleMeshMovingTest('test_execution'))
    smallSuite.addTest(TDynamicBossakTests('test_execution'))
    smallSuite.addTest(TDynamicNewmarkTests('test_execution'))
    smallSuite.addTest(TSprismMembranePatchTests('test_execution'))
    smallSuite.addTest(TSprismBendingPatchTests('test_execution'))
    smallSuite.addTest(TShellQ4ThickBendingRollUpTests('test_execution'))
    smallSuite.addTest(TShellQ4ThickDrillingRollUpTests('test_execution'))
    smallSuite.addTest(TShellT3ThinBendingRollUpTests('test_execution'))
    smallSuite.addTest(TShellT3ThinDrillingRollUpTests('test_execution'))
    smallSuite.addTest(TEigenQ4Thick2x2PlateTests('test_execution'))
    smallSuite.addTest(TEigenTL3D8NCubeTests('test_execution'))

    # Create a test suit with the selected tests plus all small tests
    nightSuite = suites['nightly']
    nightSuite.addTests(smallSuite)
    nightSuite.addTest(TShellT3IsotropicScordelisTests('test_execution'))
    
    # For very long tests that should not be in nighly and you can use to validate 
    validationSuite = suites['validation']
    validationSuite.addTest(TSprismPanTests('test_execution'))

    # Create a test suit that contains all the tests:
    allSuite = suites['all']
    allSuite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            TSimpleMeshMovingTest,
            TDynamicBossakTests,
            TDynamicNewmarkTests,
            TSprismMembranePatchTests,
            TSprismBendingPatchTests,
            TShellQ4ThickBendingRollUpTests,
            TShellQ4ThickDrillingRollUpTests,
            TShellT3ThinBendingRollUpTests,
            TShellT3ThinDrillingRollUpTests,
            TShellT3IsotropicScordelisTests
        ])
    )
    
    if( hasattr(KratosMultiphysics.ExternalSolversApplication,  "FEASTSolver") ):
        allSuite.addTests(
            KratosUnittest.TestLoader().loadTestsFromTestCases([
                TEigenQ4Thick2x2PlateTests,
                TEigenTL3D8NCubeTests
            ])
        )
    else:
        print("FEASTSolver solver is not included in the compilation of the External Solvers Application")

    return suites

if __name__ == '__main__':
    KratosUnittest.runTests(AssambleTestSuites())
```

The following can be added as commentary to clarify what the script is doing:

* The tests are added using the string `test_nametest` from the imported class.
* Nightly should consider at least all the small tests, in the same way all should consider all the nightly tests. 

The following corresponds to the script that contains the test examples: 

```python
import os

# Import Kratos
from KratosMultiphysics import *

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
import Kratos_Execute_Solid_Test as Execute_Test

# This utility will control the execution scope in case we need to acces files or we depend
# on specific relative locations of the files.

# TODO: Should we move this to KratosUnittest?
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

class StructuralMechanichsTestFactory(KratosUnittest.TestCase):

    def setUp(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Initialize GiD  I/O
            parameter_file = open(self.file_name + "_parameters.json", 'r')
            ProjectParameters = Parameters(parameter_file.read())

            # Creating the model part
            self.test = Execute_Test.Kratos_Execute_Test(ProjectParameters)

    def test_execution(self):
        # Within this location context:
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Solve()

    def tearDown(self):
        pass


class SimpleMeshMovingTest(StructuralMechanichsTestFactory):
    file_name = "mesh_moving_test/simple_mesh_moving_test"


class DynamicBossakTests(StructuralMechanichsTestFactory):
    file_name = "dynamic_test/dynamic_bossak_test"


class DynamicNewmarkTests(StructuralMechanichsTestFactory):
    file_name = "dynamic_test/dynamic_newmark_test"


class SprismMembranePatchTests(StructuralMechanichsTestFactory):
    file_name = "sprism_test/patch_membrane_test"


class SprismBendingPatchTests(StructuralMechanichsTestFactory):
    file_name = "sprism_test/patch_bending_test"


class ShellQ4ThickBendingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__BendingRollUp_test"


class ShellQ4ThickDrillingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_Q4_Thick__DrillingRollUp_test"


class ShellT3ThinBendingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__BendingRollUp_test"


class ShellT3ThinDrillingRollUpTests(StructuralMechanichsTestFactory):
    file_name = "shell_test/Shell_T3_Thin__DrillingRollUp_test"


class EigenQ4Thick2x2PlateTests(StructuralMechanichsTestFactory):
    file_name = "eigen_test/Eigen_Q4_Thick_2x2_Plate_test"


class EigenTL3D8NCubeTests(StructuralMechanichsTestFactory):
    file_name = "eigen_test/Eigen_TL_3D8N_Cube_test"

```

This file calls the main file, which is very similar to the one used by the GiD's problem type: 

```python
import os

from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

import process_factory

class Kratos_Execute_Test:

    def __init__(self, ProjectParameters):

        self.ProjectParameters = ProjectParameters

        self.main_model_part = ModelPart(self.ProjectParameters["problem_data"]["model_part_name"].GetString())
        self.main_model_part.ProcessInfo.SetValue(DOMAIN_SIZE, self.ProjectParameters["problem_data"]["domain_size"].GetInt())

        self.Model = {self.ProjectParameters["problem_data"]["model_part_name"].GetString(): self.main_model_part}

        # Construct the solver (main setting methods are located in the solver_module)
        solver_module = __import__(self.ProjectParameters["solver_settings"]["solver_type"].GetString())
        self.solver = solver_module.CreateSolver(self.main_model_part, self.ProjectParameters["solver_settings"])

        # Add variables (always before importing the model part) (it must be integrated in the ImportModelPart)
        # If we integrate it in the model part we cannot use combined solvers
        self.solver.AddVariables()

        # Read model_part (note: the buffer_size is set here) (restart can be read here)
        self.solver.ImportModelPart()

        # Add dofs (always after importing the model part) (it must be integrated in the ImportModelPart)
        # If we integrate it in the model part we cannot use combined solvers
        self.solver.AddDofs()

        # Build sub_model_parts or submeshes (rearrange parts for the application of custom processes)
        # #Get the list of the submodel part in the object Model
        for i in range(self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"].size()):
            part_name = self.ProjectParameters["solver_settings"]["processes_sub_model_part_list"][i].GetString()
            self.Model.update({part_name: self.main_model_part.GetSubModelPart(part_name)})

        # Obtain the list of the processes to be applied
        self.list_of_processes = process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["constraints_process_list"])
        self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["loads_process_list"])
        self.list_of_processes += process_factory.KratosProcessFactory(self.Model).ConstructListOfProcesses(self.ProjectParameters["list_other_processes"])

        for process in self.list_of_processes:
            process.ExecuteInitialize()

        # ### START SOLUTION ####

        self.computing_model_part = self.solver.GetComputingModelPart()

        # ### Output settings start ####
        self.problem_path = os.getcwd()
        self.problem_name = self.ProjectParameters["problem_data"]["problem_name"].GetString()

        # ### Output settings start ####
        self.output_post = ProjectParameters.Has("output_configuration")
        if (self.output_post == True):
            from gid_output_process import GiDOutputProcess
            output_settings = ProjectParameters["output_configuration"]
            self.gid_output = GiDOutputProcess(self.computing_model_part,
                                               self.problem_name,
                                               output_settings)
            self.gid_output.ExecuteInitialize()
            
        # Sets strategies, builders, linear solvers, schemes and solving info, and fills the buffer
        self.solver.Initialize()
        self.solver.SetEchoLevel(0) # Avoid to print anything 
        
        if (self.output_post == True):
            self.gid_output.ExecuteBeforeSolutionLoop()

    def Solve(self):
        for process in self.list_of_processes:
            process.ExecuteBeforeSolutionLoop()

        # #Stepping and time settings (get from process info or solving info)
        # Delta time
        delta_time = self.ProjectParameters["problem_data"]["time_step"].GetDouble()
        # Start step
        self.main_model_part.ProcessInfo[TIME_STEPS] = 0
        # Start time
        time = self.ProjectParameters["problem_data"]["start_time"].GetDouble()
        # End time
        end_time = self.ProjectParameters["problem_data"]["end_time"].GetDouble()

        # Solving the problem (time integration)
        while(time <= end_time):
            time = time + delta_time
            self.main_model_part.ProcessInfo[TIME_STEPS] += 1
            self.main_model_part.CloneTimeStep(time)

            for process in self.list_of_processes:
                process.ExecuteInitializeSolutionStep()
                
            if (self.output_post == True):
                self.gid_output.ExecuteInitializeSolutionStep()
                        
            self.solver.Solve()
            
            if (self.output_post == True):
                self.gid_output.ExecuteFinalizeSolutionStep()

            for process in self.list_of_processes:
                process.ExecuteFinalizeSolutionStep()

            for process in self.list_of_processes:
                process.ExecuteBeforeOutputStep()

            for process in self.list_of_processes:
                process.ExecuteAfterOutputStep()

            if (self.output_post == True):
                if self.gid_output.IsOutputStep():
                    self.gid_output.PrintOutput()

        if (self.output_post == True):
            self.gid_output.ExecuteFinalize()

        for process in self.list_of_processes:
            process.ExecuteFinalize()
```

The main is always the same, just the processes from the project parameters inside the `.json` are changed, look for example the `.json` that corresponds to the Bossak dynamic test: 

```json
{
    "problem_data"             : {
        "problem_name"    : "dynamic_test",
        "model_part_name" : "Structure",
        "domain_size"     : 2,
        "time_step"       : 0.001,
        "start_time"      : 0.001,
        "end_time"        : 0.10,
        "echo_level"      : 0
    },
    "solver_settings"          : {
        "solver_type"                        : "solid_mechanics_implicit_dynamic_solver",
        "echo_level"                         : 0,
        "solution_type"                      : "Dynamic",
        "time_integration_method"            : "Implicit",
        "scheme_type"                        : "Bossak",
        "model_import_settings"              : {
            "input_type"     : "mdpa",
            "input_filename" : "dynamic_test/dynamic_test"
        },
        "line_search"                        : false,
        "convergence_criterion"              : "Residual_criterion",
        "displacement_relative_tolerance"    : 0.0001,
        "displacement_absolute_tolerance"    : 1e-9,
        "residual_relative_tolerance"        : 0.0001,
        "residual_absolute_tolerance"        : 1e-9,
        "max_iteration"                      : 10,
        "linear_solver_settings"             : {
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
        },
        "problem_domain_sub_model_part_list" : ["Parts_Parts_Auto1"],
        "processes_sub_model_part_list"      : ["DISPLACEMENT_Displacement_Auto1"]
    },
    "constraints_process_list" : [
    {
        "python_module"   : "impose_vector_value_by_components_process",
        "kratos_module" : "KratosMultiphysics",
        "help"                  : "This process fixes the selected components of a given vector variable",
        "process_name"          : "ImposeVectorValueByComponentsProcess",
        "Parameters"            : {
            "mesh_id"         : 0,
            "model_part_name" : "DISPLACEMENT_Displacement_Auto1",
            "variable_name"   : "DISPLACEMENT",
            "is_fixed_x"      : false,
            "is_fixed_y"      : true,
            "is_fixed_z"      : true,
            "value"           : [0.1, 0.0, 0.0]
        }
    }
    ],
    "loads_process_list" : [],
    "list_other_processes" :[
    {
        "python_module"   : "from_json_check_result_process",
        "kratos_module" : "KratosMultiphysics",
        "help"                  : "",
        "process_name"          : "FromJsonCheckResultProcess",
        "Parameters"            : {
            "check_variables" : ["DISPLACEMENT_X"],
            "input_file_name" : "dynamic_test/dynamic_bossak_test_results.json",
            "model_part_name"  : "DISPLACEMENT_Displacement_Auto1",
            "time_frequency"   : 0.01
        }
    }
    ],
    "print_output_process" : [
    {
        "python_module"   : "json_output_process",
        "kratos_module" : "KratosMultiphysics",
        "help"                  : "",
        "process_name"          : "JsonOutputProcess",
        "Parameters"            : {
            "output_variables" : ["DISPLACEMENT_X","VELOCITY_X","ACCELERATION_X"],
            "output_file_name" : "dynamic_test/dynamic_bossak_test_results.json",
            "model_part_name"  : "DISPLACEMENT_Displacement_Auto1",
            "time_frequency"   : 0.01
        }
    }
    ],
    "check_json_results_process" : [
    {
        "python_module"   : "from_json_check_result_process",
        "kratos_module" : "KratosMultiphysics",
        "help"                  : "",
        "process_name"          : "FromJsonCheckResultProcess",
        "Parameters"            : {
            "check_variables" : ["DISPLACEMENT_X","VELOCITY_X","ACCELERATION_X"],
            "input_file_name" : "dynamic_test/dynamic_bossak_test_results.json",
            "model_part_name"  : "DISPLACEMENT_Displacement_Auto1",
            "time_frequency"   : 0.01
        }
    }
    ],
    "check_analytic_results_process" : [
    {
        "python_module"   : "from_analytic_check_result_process",
        "kratos_module" : "KratosMultiphysics",
        "help"                  : "",
        "process_name"          : "FromAnalyticCheckResultProcess",
        "Parameters"            : {
            "variable_name"     : "DISPLACEMENT_X",
            "mesh_id"           : 0,
            "f(x,y,z,t)="       : "cos(10.0*t)",
            "model_part_name"   : "DISPLACEMENT_Displacement_Auto1",
            "time_frequency"    : 0.01
        }
    }
    ],    
    "apply_custom_function_process" : [],
    "restart_options"          : {
        "SaveRestart"      : false,
        "RestartFrequency" : 0,
        "LoadRestart"      : false,
        "Restart_Step"     : 0
    },
    "constraints_data"         : {
        "incremental_load"         : false,
        "incremental_displacement" : false
    }
}

```

In the previous example the tests can be checked both analitycally with `check_analytic_results_process` or numerically with `check_json_results_process` (for that check you need to generate first the `.json` results file with the `print_output_process`).

For some specific problems, where the common processes are not enough to define your test, you can define a local process, like the following one, where the displacement is imposed in the boundary and after solving the displacements they are checked in all the domain with 'AssertAlmostEqual': 

```python
# Importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *

CheckForPreviousImport()

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

def Factory(settings, Model):
    if (type(settings) != Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLocalProcess(Model, settings["Parameters"])

class ApplyLocalProcess(Process, KratosUnittest.TestCase):

    def __init__(self,model_part,params):

        self.model_part = model_part[params["model_part_name"].GetString()]
        self.params = params
        
    def ExecuteInitialize(self):
        # Find neighbours if required
        sprism_neighbour_search = SprismNeighbours(self.model_part)
        sprism_neighbour_search.Execute()
        #pass
        
    def ExecuteBeforeSolutionLoop(self):
        # Add BC
        for node in self.model_part.Nodes:
            if (node.X >2.40000e-01 -1.0e-5) | (node.Y > 1.20000e-01 -1.0e-5) | (node.X < 1.0e-5) | (node.Y < 1.0e-5):
                node.Fix(DISPLACEMENT_X)
                node.SetSolutionStepValue(DISPLACEMENT_X, 0, -1.0e-7 * (node.Z - 0.0005) * (node.X + node.Y / 2))
                node.Fix(DISPLACEMENT_Y)
                node.SetSolutionStepValue(DISPLACEMENT_Y, 0, -1.0e-7 * (node.Z - 0.0005) * (node.Y + node.X / 2))
                node.Fix(DISPLACEMENT_Z)
                node.SetSolutionStepValue(DISPLACEMENT_Z, 0, 0.5 * 1.0e-7 * (node.X ** 2 + node.X * node.Y + node.Y ** 2))

    
    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        for node in self.model_part.Nodes:
            value = node.GetSolutionStepValue(DISPLACEMENT_X,0)
            self.assertAlmostEqual(value, -1.0e-7 * (node.Z - 0.0005) * (node.X + node.Y / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Y,0)
            self.assertAlmostEqual(value,-1.0e-7 * (node.Z - 0.0005) * (node.Y + node.X / 2))
            value = node.GetSolutionStepValue(DISPLACEMENT_Z,0)
            self.assertAlmostEqual(value, 0.5 * 1.0e-7 * (node.X ** 2 + node.X * node.Y + node.Y **2 ))
            
    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteFinalize(self):
        pass
```

This specifics results can be found in the StructuralMechanicsApplication tests folder.

Important additional comments:

* The solver, constitutive laws, etc..., everything is created from the main file, that should be **ALWAYS** the same, and just change the processes, or even include your local processes, but never change the main.
 * You can create a json with data from previous simulations for future benchmarking using the process json_output_process.py and check it with from_json_check_result_process.py.
* In order to be recognized automatically, all tests need to have the name "test_" at the beginning of the definition (`test_MembranePatch` and `test_BendingPatch` in this example).
* The `tearDown` doesn't do anything, it just closes files, clears memory...
    `imports` should be at the beginning of the file to avoid conflicts.
* **NO PRINTS**: please don't print anything on screen by calling the `print()` function. It is not necessary, the unittest already prints all the necessary information. 

### Some common commands in Unittest

The [link](https://docs.python.org/2/library/unittest.html)  contents the commands that can be used with the original Unittest from Python, not all the commands from the python Unittest are available in Kratos right now, but they will be available in a near future. The following are the most common and essential commands in Kratos Unittest:

* `assertTrue/assertFalse`:They check if the two inputs are true or false respectively.
* `assertEqual`: Check if the two inputs are exactly equal, if the inputs are doubles it is recommended to use the next assert, in order to avoid precision errors.
* `assertAlmostEqual`: Check if the solution is almost equal in the two inputs, with a certain precision (check the unittest python wiki for more details [link](https://docs.python.org/2/library/unittest.html#unittest.TestCase.assertAlmostEqual)) 

### Skipping tests with external dependencies

In general, tests for classes and functions defined in the Kratos Core should only use core objects and tests for classes and functions defined in an application should only use objects from the core or that application. However, in some cases it might be necessary to import something (an element, an utility,...) that belongs to a different application. When that happens, we have to account for the fact that our end user might not be compiling that extra application, causing the test to fail.

One way to solve this problem is to skip tests that require different applications if that application could not be found at runtime. This section shows how to do so. For the following example, assume that we are defining a new test for the FSI Application, which will test a new tool to do a coupled FSI simulation using features from the Structural Mechanics Application and the Fluid Dynamics Application. 

When we define our test case, we use a `KratosMultiphysics.kratos_utilities.CheckIfApplicationsAvailable` to import other applications. Note that this is only necessary if functionalities from the other applications are directly used. If they are implicitly used (i.e. though json-input) then it is not necessary to import the other applications.

```python
import KratosMultiphysics
from KratosMultiphysics import kratos_utilities

if kratos_utilities.CheckIfApplicationsAvailable("StructuralMechanicsApplication", "FluidDynamicsApplication")
    from KratosMultiphysics import FluidDynamicsApplication
    from KratosMultiphysics import StructuralMechanicsApplication
from KratosMultiphysics import FSIApplication
```

This will allow the script to run even if the applications are not available. Note that code execution will still fail if we try to use anything from a missing application.

Next, when we define a test that needs something from the external applications, we specify that it should be skipped if either the Structural Mechanics Application or the Fluid Dynamics Application are not available. This can be done using the `skipIfApplicationsNotAvailable` decorator:

```python
@KratosUnittest.skipIfApplicationsNotAvailable("StructuralMechanicsApplication", "FluidDynamicsApplication")
class FSIProblemEmulatorTest(FSIProblemEmulatorTestFactory):
    # Define the test as normal
```

The same can be done with using `KratosUnittest.skipTestIfApplicationsNotAvailable`:

```python
class FSIProblemEmulatorTest(FSIProblemEmulatorTestFactory):
    def test_FSI_problem(self):
        self.skipTestIfApplicationsNotAvailable("StructuralMechanicsApplication", "FluidDynamicsApplication")
        # Define the test as normal
```

Now, if we try to run the tests but we have not compiled one of the extra applications, the test will be skipped. If we set the verbosity level to `2`, the test script will tell us why it was skipped.

```
$ python run_tests.py -v 2 -a 'FSIApplication'
...
test_execution (SmallTests.FSIProblemEmulatorTest) ... skipped 'Required Applications are missing: "StructuralMechanicsApplication", "FluidDynamicsApplication"'
```

Note that, if the test uses auxiliary python script files, the `import` statement for the extra applications should be placed in a `try`/`except` block in these files too, otherwise it will produce a runtime error.

## C++

### How to run tests 

First remember to add the following to your `configure.sh`:

```console
-DKRATOS_BUILD_TESTING=ON                                                                                   \
```

For the tests defined in C++ it is necessary to call the `Tester` from `KratosMultiPhysics`. Usage is the following: 

```python
from KratosMultiphysics import *

# If you want to run tests defined in an application, import it here.

Tester.SetVerbosity(Tester.Verbosity.TESTS_OUTPUTS) # Set the verbosity level

Tester.RunAllTestCases() #Test all cases

Tester.RunTestSuite("Suite_name_here") #Test a whole suite

Tester.RunTestCases("Test_name_here") #Test a specific case
```

But can consult the functioning of the C++ cases looking at the [code](https://github.com/KratosMultiphysics/Kratos/blob/master/kratos/python/add_testing_to_python.cpp). In general the main arguments are:

*  `SetVerbosity()` : To set the verbosity level. You can choose between:
    * `QUITE`: No message will be shown
    * `PROGRESS`: Only the progress of the tests run is shown 
    * `TESTS_LIST`: Shows the list of tests
    * `FAILED_TESTS_OUTPUTS`: It returns the tests outputs (only in the failing ones)
    * `TESTS_OUTPUTS`: It returns all the tests outputs
* `ProfileAllTestCases()`:
* `ProfileTestSuite()`:
* `ProfileAllTestCases()`:
* `NumberOfFailedTestCases()`: it returns the number of testcases failed.
* `ResetAllTestCasesResults()`: 
* `ListOfAllTestCases()`: It lists all the avalaible test cases

### MPI

Please notice that when executing the tests with MPI, all following severity levels: `QUITE`, `PROGRESS`, `TESTS_LIST`, `FAILED_TESTS_OUTPUTS`, will omit the output from all processes except the master, and only `TESTS_OUTPUTS` will provide you with the output of all processes.

In order to ease the lecture of the outputs which will be mixed up in the terminal, we recommend running the tests with:

`mpirun -np [N] -output-filename [filename_folder]` 

So a separate `stdout` and `stderr` are created for every rank of the execution in the `filename_folder` specified.

### General structure

In order to add a test to some of these suites, one can do it as shown in this example: 

```C++
// System includes

// External includes

// Project includes
#include "testing/testing.h"

namespace Kratos {
namespace Testing {

	KRATOS_TEST_CASE_IN_SUITE(TestSuite, KratosCoreFastSuite)
	{
	  Tester::CreateTestSuite("MyTestSuite")
	  {
	  }
	}
}
} // namespace Kratos.
```

This was the simplest example, tests with more complexity can be created, as the following one:

```C++
// System includes
#include <limits>

// External includes


// Project includes
#include "testing/testing.h"

// Utility includes
#include "utilities/math_utils.h"

namespace Kratos 
{
    namespace Testing 
    {
        constexpr double EPSILON = std::numeric_limits<double>::epsilon();
        constexpr double TOLERANCE = 1e-6;
        
        /// Tests
      
        /** Checks if the area of the triangle is calculated correctly using Heron equation.
         * Checks if the area of the triangle is calculated correctly using Heron equation.
         */
        
        KRATOS_TEST_CASE_IN_SUITE(MathUtilsHeronTest, KratosCoreMathUtilsFastSuite) 
        {
            const double area = MathUtils<double>::Heron<false>(std::sqrt(2.0), 1.0, 1.0);

            KRATOS_EXPECT_NEAR(area, 0.5, TOLERANCE);
        }
    } // namespace Testing
}  // namespace Kratos.
```

Remember to add the following to the CmakeLists.txt of 

```CMake
if(${KRATOS_BUILD_TESTING} MATCHES ON)
	file(GLOB_RECURSE KRATOS_NAME_YOUR_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)
```

You need to add your `KRATOS_NAME_YOUR_APPLICATION_TESTING_SOURCES` to your library, for example:

```CMake
add_library(KratosYourApplicationApplication SHARED ${KRATOS_NAME_YOUR_APPLICATION_SOURCES} ${KRATOS_NAME_YOUR_APPLICATION_TESTING_SOURCES})
```

It is recommended to put the tests files in a folder called `cpp_tests` inside the tests folder in order to follow the same structure as the **Core** and the **CoreApplications**.

> **Please Note**
>
> To run tests (or suites) defined by an application, you need to import that application in the python script used to launch the tests.

### Some common commands

The following are the most common and essential commands in Kratos Unittest:

* `KRATOS_EXPECT_TRUE`/`KRATOS_EXPECT_FALSE`:They check if the input is true or false
* `KRATOS_CHECK_EXCEPTION_RAISED`: Check if an exception is raised
* `KRATOS_EXPECT_STREQ`/`KRATOS_CHECK_C_STRING__NOT_EQUAL`/`KRATOS_EXPECT_HAS_SUBSTRING`: Check if the string is the same, or not, or a substring
* `KRATOS_EXPECT_EQ`: Check if the two inputs are exactly equal, if the inputs are doubles it is recommended to use the next assert, in order to avoid precision errors.
* `KRATOS_EXPECT_NEAR`: Check if the solution is almost equal in the two inputs, with a certain precision 
* `KRATOS_EXPECT_LT`/`KRATOS_EXPECT_GT`/`KRATOS_EXPECT_LE`/`KRATOS_EXPECT_GE`:  It checks if the value is <, <=, > or >= than the other 
