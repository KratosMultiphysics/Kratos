from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing Kratos
import KratosMultiphysics

# Other imports
import os

class PythonSolver(object):
    """The base class for the Python Solvers in the applications
    Changes to this BaseClass have to be discussed first!
    """
    def __init__(self, model, settings):
        """The constructor of the PythonSolver-Object.

        It is intended to be called from the constructor
        of deriving classes:
        super(DerivedSolver, self).__init__(settings)

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- The Model to be used
        settings -- The solver settings used
        """
        if (type(model) != KratosMultiphysics.Model):
            raise Exception("Input is expected to be provided as a Kratos Model object")

        if (type(settings) != KratosMultiphysics.Parameters):
            raise Exception("Input is expected to be provided as a Kratos Parameters object")

        self.model = model
        self.settings = settings

        self.echo_level = 0 # default to zero
        if self.settings.Has("echo_level"):
            self.echo_level = self.settings["echo_level"].GetInt()

    def AddVariables(self):
        """This function add the Variables needed by this PythonSolver to the the ModelPart
        It has to be called BEFORE the ModelPart is read!
        """
        pass

    def AddDofs(self):
        """This function add the Dofs needed by this PythonSolver to the the ModelPart
        It has to be called AFTER the ModelPart is read!
        """
        pass

    def ImportModelPart(self):
        """This function reads the ModelPart
        """
        raise Exception("This function has to be implemented in the derived class")

    def PrepareModelPart(self):
        """This function prepares the ModelPart for being used by the PythonSolver
        """
        pass

    def GetMinimumBufferSize(self):
        """This function returns the minimum buffer size needed for this PythonSolver
        """
        raise Exception('Please implement "GetMinimumBufferSize" in your derived solver')

    def ExportModelPart(self):
        """This function exports the ModelPart to and mdpa-file
        """
        raise Exception("This function has to be implemented in the derived class")

    def AdvanceInTime(self, current_time):
        """This function advances the PythonSolver in time
        Usage: It is designed to be called once per solution step, before performing the solution
        """
        pass

    def Initialize(self):
        """This function initializes the PythonSolver
        Usage: It is designed to be called ONCE, BEFORE the execution of the solution-loop
        """
        pass

    def Finalize(self):
        """This function finalizes the PythonSolver
        Usage: It is designed to be called ONCE, AFTER the execution of the solution-loop
        """
        pass

    def Predict(self):
        """This function performs all the required operations that should be executed
        (for each step) ONCE, AFTER initializing the solution step.
        """
        pass

    def InitializeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) BEFORE solving the solution step.
        """
        pass

    def FinalizeSolutionStep(self):
        """This function performs all the required operations that should be executed
        (for each step) AFTER solving the solution step.
        """
        pass

    def SolveSolutionStep(self):
        """This function solves the current step.
        It can be called multiple times within one solution step
        """
        pass

    def Check(self):
        """This function checks the PythonSolver. It usually calls the "Check" function of a solving strategy
        """
        pass

    def Clear(self):
        """This function clears the PythonSolver
        """
        pass

    def Solve(self):
        warning_msg  = 'Using "Solve" is deprecated and will be removed in the future!\n'
        warning_msg += 'Use the separate calls to "Initialize", "InitializeSolutionStep", "Predict", '
        warning_msg += '"SolveSolutionStep" and "FinalizeSolutionStep"'
        self.print_warning_on_rank_zero("::[PythonSolver]::", warning_msg)
        self.Initialize()
        self.Predict()
        self.InitializeSolutionStep()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    def GetComputingModelPart(self):
        raise Exception("This function has to be implemented in the derived class")

    def _ImportModelPart(self, model_part, model_part_import_settings):
        """This function imports the ModelPart
        """
        self.print_on_rank_zero("::[PythonSolver]::", "Reading model part.")
        input_type = model_part_import_settings["input_type"].GetString()

        if (input_type == "mdpa"):
            problem_path = os.getcwd()
            input_filename = model_part_import_settings["input_filename"].GetString()
            import_flags = KratosMultiphysics.ModelPartIO.READ
            if model_part_import_settings.Has("ignore_variables_not_in_solution_step_data"):
                if model_part_import_settings["ignore_variables_not_in_solution_step_data"].GetBool():
                    import_flags = KratosMultiphysics.ModelPartIO.IGNORE_VARIABLES_ERROR|KratosMultiphysics.ModelPartIO.READ

            # Import model part from mdpa file.
            self.print_on_rank_zero("::[PythonSolver]::", "Reading model part from file: " + os.path.join(problem_path, input_filename) + ".mdpa")
            KratosMultiphysics.ModelPartIO(input_filename, import_flags).ReadModelPart(model_part)
            if (model_part_import_settings.Has("reorder") and model_part_import_settings["reorder"].GetBool()):
                tmp = KratosMultiphysics.Parameters("{}")
                KratosMultiphysics.ReorderAndOptimizeModelPartProcess(model_part, tmp).Execute()
            self.print_on_rank_zero("::[PythonSolver]::", "Finished reading model part from mdpa file.")
        elif (input_type == "rest"):
            self.print_on_rank_zero("::[PythonSolver]::", "Loading model part from restart file.")
            from restart_utility import RestartUtility
            RestartUtility(model_part, self._GetRestartSettings(model_part_import_settings)).LoadRestart()
            self.print_on_rank_zero("::[PythonSolver]::", "Finished loading model part from restart file.")
        elif(input_type == "use_input_model_part"):
            pass
        else:
            raise Exception("Other model part input options are not yet implemented.")
        self.print_on_rank_zero("ModelPart", model_part)
        self.print_on_rank_zero("::[PythonSolver]:: ", "Finished reading model part.")

    def _GetRestartSettings(self, model_part_import_settings):
        restart_settings = model_part_import_settings.Clone()
        restart_settings.RemoveValue("input_type")
        if not restart_settings.Has("restart_load_file_label"):
            raise Exception('"restart_load_file_label" must be specified when starting from a restart-file!')
        if model_part_import_settings.Has("echo_level"):
            restart_settings.AddValue("echo_level", model_part_import_settings["echo_level"])

        return restart_settings

    #### Auxiliar functions ####

    def print_on_rank_zero(self, *args):
        # This function will be overridden in the trilinos-solvers
        KratosMultiphysics.Logger.PrintInfo(" ".join(map(str,args)))

    def print_warning_on_rank_zero(self, *args):
        # This function will be overridden in the trilinos-solvers
        KratosMultiphysics.Logger.PrintWarning(" ".join(map(str,args)))

    def validate_and_transfer_matching_settings(self, origin_settings, destination_settings):
        """Transfer matching settings from origin to destination.

        If a name in origin matches a name in destination, then the setting is
        validated against the destination.

        The typical use is for validating and extracting settings in derived classes:

        class A:
            def __init__(self, model_part, a_settings):
                default_a_settings = Parameters('''{
                    ...
                }''')
                a_settings.ValidateAndAssignDefaults(default_a_settings)
        class B(A):
            def __init__(self, model_part, custom_settings):
                b_settings = Parameters('''{
                    ...
                }''') # Here the settings contain default values.
                self.validate_and_transfer_matching_settings(custom_settings, b_settings)
                super().__init__(model_part, custom_settings)
        """
        for name, dest_value in destination_settings.items():
            if origin_settings.Has(name): # Validate and transfer value.
                orig_value = origin_settings[name]
                if dest_value.IsDouble() and orig_value.IsDouble():
                    destination_settings[name].SetDouble(origin_settings[name].GetDouble())
                elif dest_value.IsInt() and orig_value.IsInt():
                    destination_settings[name].SetInt(origin_settings[name].GetInt())
                elif dest_value.IsBool() and orig_value.IsBool():
                    destination_settings[name].SetBool(origin_settings[name].GetBool())
                elif dest_value.IsString() and orig_value.IsString():
                    destination_settings[name].SetString(origin_settings[name].GetString())
                elif dest_value.IsArray() and orig_value.IsArray():
                    if dest_value.size() != orig_value.size():
                        raise Exception('len("' + name + '") != ' + str(dest_value.size()))
                    for i in range(dest_value.size()):
                        if dest_value[i].IsDouble() and orig_value[i].IsDouble():
                            dest_value[i].SetDouble(orig_value[i].GetDouble())
                        elif dest_value[i].IsInt() and orig_value[i].IsInt():
                            dest_value[i].SetInt(orig_value[i].GetInt())
                        elif dest_value[i].IsBool() and orig_value[i].IsBool():
                            dest_value[i].SetBool(orig_value[i].GetBool())
                        elif dest_value[i].IsString() and orig_value[i].IsString():
                            dest_value[i].SetString(orig_value[i].GetString())
                        elif dest_value[i].IsSubParameter() and orig_value[i].IsSubParameter():
                            self.validate_and_transfer_matching_settings(orig_value[i], dest_value[i])
                            if len(orig_value[i].items()) != 0:
                                raise Exception('Json settings not found in default settings: ' + orig_value[i].PrettyPrintJsonString())
                        else:
                            raise Exception('Unsupported parameter type.')
                elif dest_value.IsSubParameter() and orig_value.IsSubParameter():
                    self.validate_and_transfer_matching_settings(orig_value, dest_value)
                    if len(orig_value.items()) != 0:
                        raise Exception('Json settings not found in default settings: ' + orig_value.PrettyPrintJsonString())
                else:
                    raise Exception('Unsupported parameter type.')
                origin_settings.RemoveValue(name)
