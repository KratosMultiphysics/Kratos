from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics
from python_sovler import PythonSolver

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("FluidDynamicsApplication")

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

def CreateSolver(main_model_part, custom_settings):
    return NavierStokesBaseSolver(main_model_part, custom_settings)

class NavierStokesBaseSolver(PythonSolver):

    def __init__(self, model, custom_settings):
        super(NavierStokesBaseSolver,self).__init___(model, custom_settings)

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = None
        self.condition_name = None
        self.min_buffer_size = 3

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = ModelPart(model_part_name)
            self.model.AddModelPart(self.main_model_part)

    def AddVariables(self):
        raise Exception("Trying to call NavierStokesBaseSolver.AddVariables(). Implement the AddVariables() method in the specific derived solver.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Fluid solver DOFs added correctly.")

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart()

    def PrepareModelPart(self):
        ## Replace default elements and conditions
        self._replace_elements_and_conditions()
        ## Executes the check and prepare model process
        self._execute_check_and_prepare()
        ## Set buffer size
        self.main_model_part.SetBufferSize(self.min_buffer_size)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Model reading finished.")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("NavierStokesBaseSolver", "Model export finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Initialize(self):
        raise Exception("Calling the Navier-Stokes base solver. Please implement the custom Initialize() method of your solver.")

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

    def Check(self):
        (self.solver).Check()

    def Clear(self):
        (self.solver).Clear()

    def Solve(self):
        message = "".join(
            "Calling NavierStokesBaseSolver.Solve() method, which is deprecated\n",
            "Please call the individual methods instead:\n",
            "solver.InitializeSolutionStep()\n",
            "solver.Predict()\n",
            "solver.SolveSolutionStep()\n",
            "solver.FinalizeSolutionStep()\n"
        )
        KratosMultiphysics.Logger.PrintWarning("NavierStokesBaseSolver",message)
        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    ## Fluid-specific additions (?)

    def AdaptMesh(self):
        pass

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

    def SaveRestart(self):
        pass #one should write the restart file here


    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def _get_automatic_time_stepping_utility(self):
        if (self.computing_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility2D(self.computing_model_part,
                                                                     self.settings["time_stepping"])
        else:
            EstimateDeltaTimeUtility = KratosCFD.EstimateDtUtility3D(self.computing_model_part,
                                                                     self.settings["time_stepping"])

        return EstimateDeltaTimeUtility

    ## The following are required by derived classes

    def _TimeBufferIsInitialized(self):
        # We always have one extra old step (step 0, read from input)
        return self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] + 1 >= self.GetMinimumBufferSize()

    def _IsPrintingRank(self):
        return self._is_printing_rank

    def _GetDefaultSettings(self):
        raise Exception("Please define the default solver settings in the derived solver class")

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

    def _execute_check_and_prepare(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

        import check_and_prepare_model_process_fluid
        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()