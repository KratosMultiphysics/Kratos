from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesBaseSolver(main_model_part, custom_settings)

class NavierStokesBaseSolver(object):

    def __init__(self, main_model_part, custom_settings):
        ## Set the element and condition names for the replace settings
        self.element_name = None
        self.condition_name = None
        self.min_buffer_size = 3

        KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Construction of NavierStokesBaseSolver finished.")

    def AddVariables(self):
        raise Exception("Trying to add Navier-Stokes base solver variables. Implement the AddVariables() method in the specific derived solver.")

    def ImportModelPart(self):
        ## Read model part
        self._model_part_reading()
        ## Replace default elements and conditions
        self._replace_elements_and_conditions()
        ## Executes the check and prepare model process
        self._execute_check_and_prepare()
        ## Set buffer size
        self._set_buffer_size()

        KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Model reading finished.")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Model export finished.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Fluid solver DOFs added correctly.")

    def AdaptMesh(self):
        pass

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        # Automatic time step computation according to user defined CFL number
        if (self.settings["time_stepping"]["automatic_time_step"].GetBool()):
            delta_time = self.EstimateDeltaTimeUtility.EstimateDt()
        # User-defined delta time
        else:
            delta_time = self.settings["time_stepping"]["time_step"].GetDouble()

        return delta_time

    def Initialize(self):
        raise Exception("Calling the Navier-Stokes base solver. Please implement the custom Initialize() method of your solver.")

    def SaveRestart(self):
        pass #one should write the restart file here

    def Clear(self):
        (self.solver).Clear()

    def Check(self):
        (self.solver).Check()

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def InitializeSolutionStep(self):
        (self.solver).InitializeSolutionStep()

    def Predict(self):
        (self.solver).Predict()

    def SolveSolutionStep(self):
        is_converged = (self.solver).SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        (self.solver).FinalizeSolutionStep()

    def Solve(self):
        # We always have one extra old step (step 0, read from input)
        if self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.min_buffer_size:
            self.solver.Solve()

    def _model_part_reading(self):
        ## Model part reading
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            ## Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            if(self.settings["reorder"].GetBool()):
                tmp = KratosMultiphysics.Parameters("{}")
                KratosMultiphysics.ReorderAndOptimizeModelPartProcess(self.main_model_part, tmp).Execute()
        else:
            raise Exception("Other input options are not implemented yet.")

    def _execute_check_and_prepare(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

        import check_and_prepare_model_process_fluid
        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()


    def _set_buffer_size(self):
        current_buffer_size = self.main_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.main_model_part.SetBufferSize(self.min_buffer_size)

    def _get_automatic_time_stepping_utility(self):
        if (self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.computing_model_part,
                                                                     self.settings["time_stepping"])
        else:
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.computing_model_part,
                                                                     self.settings["time_stepping"])

        return EstimateDeltaTimeUtility

    def _replace_elements_and_conditions(self):
        ## Get element number of nodes and domain size
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        elem_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
        cond_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(domain_size) + "D" + str(elem_num_nodes) + "N"
        else:
            raise Exception("There is no element name. Implement the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(domain_size) + "D" + str(cond_num_nodes) + "N"
        else:
            raise Exception("There is no condition name. Implement the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{}""")
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
