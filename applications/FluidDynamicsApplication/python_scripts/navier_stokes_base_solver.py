from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

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
        self.step = 0

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

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

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Model reading finished.")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Model export finished.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)

        if self._IsPrintingRank():
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

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

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

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            self.solver.InitializeSolutionStep()

    def Predict(self):
        if self._TimeBufferIsInitialized():
            self.solver.Predict()

    def SolveSolutionStep(self):
        if self._TimeBufferIsInitialized():
            is_converged = self.solver.SolveSolutionStep()
            if not is_converged and self._IsPrintingRank():
                msg  = "Navier-Stokes solver did not converge for iteration " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
                msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
                KratosMultiphysics.Logger.PrintWarning("NavierStokesBaseSolver",msg)

    def FinalizeSolutionStep(self):
        if self._TimeBufferIsInitialized():
            (self.solver).FinalizeSolutionStep()

    def Solve(self):
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

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
        ## Get number of nodes and domain size
        elem_num_nodes = self._get_element_num_nodes()
        cond_num_nodes = self._get_condition_num_nodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## If there are no elements and/or conditions, default to triangles/tetra meshes to avoid breaking the ReplaceElementsAndConditionsProcess
        ## This only affects the input name (if there are no elements or conditions to replace, nothing is replaced).
        if elem_num_nodes == 0:
            elem_num_nodes = domain_size + 1
        if cond_num_nodes == 0:
            cond_num_nodes = domain_size

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D" + str(int(elem_num_nodes)) + "N"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        #self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""{}""")
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    def _get_element_num_nodes(self):
        if self.main_model_part.NumberOfElements() != 0:
            if sys.version_info[0] >= 3: # python3 syntax
                element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
            else: # python2 syntax
                element_num_nodes = len(self.main_model_part.Elements.__iter__().next().GetNodes())
        else:
            element_num_nodes = 0

        element_num_nodes = self.main_model_part.GetCommunicator().MaxAll(element_num_nodes)
        return element_num_nodes

    def _get_condition_num_nodes(self):
        if self.main_model_part.NumberOfConditions() != 0:
            if sys.version_info[0] >= 3: # python3 syntax
                condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())
            else: # python2 syntax
                condition_num_nodes = len(self.main_model_part.Conditions.__iter__().next().GetNodes())
        else:
            condition_num_nodes = 0

        condition_num_nodes = self.main_model_part.GetCommunicator().MaxAll(condition_num_nodes)
        return condition_num_nodes

    def _IsPrintingRank(self):
        return self._is_printing_rank

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.GetMinimumBufferSize()
