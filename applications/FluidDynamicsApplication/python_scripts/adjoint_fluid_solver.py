from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics
from python_solver import PythonSolver

import KratosFluidDynamicsApplication as KratosCFD

def CreateSolver(model, custom_settings):
    return AdjointFluidSolver(model, custom_settings)

class AdjointFluidSolver(PythonSolver):

    def __init__(self, model, custom_settings):
        settings = self._ValidateSettings(custom_settings)

        super(AdjointFluidSolver,self).__init__(model, settings)

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        ## Set the element and condition names for the replace settings
        ## These should be defined in derived classes
        self.element_name = None
        self.condition_name = None
        self.min_buffer_size = 2

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    def AddVariables(self):
        raise Exception("Trying to call AdjointFluidSolver.AddVariables(). Implement the AddVariables() method in the specific derived solver.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_VECTOR_1_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosCFD.ADJOINT_FLUID_SCALAR_1, self.main_model_part)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Adjoint fluid solver DOFs added correctly.")

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Replace default elements and conditions
            self._ReplaceElementsAndConditions()
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.min_buffer_size)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model reading finished.")

    def ExportModelPart(self):
        ## Model part writing
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Model export finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Initialize(self):
        raise Exception("Calling AdjointFluidSolver.Initialize() base method. Please implement a custom Initialize() method for your solver.")

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()
        self.response_function.InitializeSolutionStep()


    def Predict(self):
        self.solver.Predict()

    def SolveSolutionStep(self):
        self.solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        (self.solver).FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()
        self.sensitivity_builder.UpdateSensitivities()

    def Check(self):
        (self.solver).Check()

    def Clear(self):
        (self.solver).Clear()

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    ## FluidSolver specific methods.

    def _IsPrintingRank(self):
        return self._is_printing_rank

    def _ValidateSettings(self, settings):
        raise Exception("Please define the _ValidateSettings() method in your derived solver class to validate the Kratos::Parameters configuration.")
        # Suggested implementation:
        #settings.ValidateAndAssignDefaults(KratosMultiphysics.Parameters(r'{}'))
        #return settings

    def _ReplaceElementsAndConditions(self):
        ## Get number of nodes and domain size
        elem_num_nodes = self._GetElementNumNodes()
        cond_num_nodes = self._GetConditionNumNodes()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        ## If there are no elements and/or conditions, default to triangles/tetra meshes to avoid breaking the ReplaceElementsAndConditionsProcess
        ## This only affects the input name (if there are no elements or conditions to replace, nothing is replaced).
        if elem_num_nodes == 0:
            elem_num_nodes = domain_size + 1
        if cond_num_nodes == 0:
            cond_num_nodes = domain_size

        ## Complete the element name
        if (self.element_name is not None):
            new_elem_name = self.element_name + str(int(domain_size)) + "D"
        else:
            raise Exception("There is no element name. Define the self.element_name string variable in your derived solver.")

        ## Complete the condition name
        if (self.condition_name is not None):
            new_cond_name = self.condition_name + str(int(domain_size)) + "D" + str(int(cond_num_nodes)) + "N"
        else:
            raise Exception("There is no condition name. Define the self.condition_name string variable in your derived solver.")

        ## Set the element and condition names in the Json parameters
        self.settings.AddValue("element_replace_settings", KratosMultiphysics.Parameters("""{}"""))
        self.settings["element_replace_settings"].AddEmptyValue("element_name").SetString(new_elem_name)
        self.settings["element_replace_settings"].AddEmptyValue("condition_name").SetString(new_cond_name)

        ## Call the replace elements and conditions process
        KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()

    def _GetElementNumNodes(self):
        if self.main_model_part.NumberOfElements() != 0:
            if sys.version_info[0] >= 3: # python3 syntax
                element_num_nodes = len(self.main_model_part.Elements.__iter__().__next__().GetNodes())
            else: # python2 syntax
                element_num_nodes = len(self.main_model_part.Elements.__iter__().next().GetNodes())
        else:
            element_num_nodes = 0

        element_num_nodes = self.main_model_part.GetCommunicator().MaxAll(element_num_nodes)
        return element_num_nodes

    def _GetConditionNumNodes(self):
        if self.main_model_part.NumberOfConditions() != 0:
            if sys.version_info[0] >= 3: # python3 syntax
                condition_num_nodes = len(self.main_model_part.Conditions.__iter__().__next__().GetNodes())
            else: # python2 syntax
                condition_num_nodes = len(self.main_model_part.Conditions.__iter__().next().GetNodes())
        else:
            condition_num_nodes = 0

        condition_num_nodes = self.main_model_part.GetCommunicator().MaxAll(condition_num_nodes)
        return condition_num_nodes

    def _ExecuteCheckAndPrepare(self):
        ## Check that the input read has the shape we like
        prepare_model_part_settings = KratosMultiphysics.Parameters("{}")
        prepare_model_part_settings.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
        prepare_model_part_settings.AddValue("skin_parts",self.settings["skin_parts"])

        import check_and_prepare_model_process_fluid
        check_and_prepare_model_process_fluid.CheckAndPrepareModelProcess(self.main_model_part, prepare_model_part_settings).Execute()

        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )

    def _ComputeDeltaTime(self):
        if self.settings["time_stepping"]["automatic_time_step"].GetBool():
            raise Exception("Automatic time stepping is not supported by adjoint fluid solver.")

        delta_time = self.settings["time_stepping"]["time_step"].GetDouble()
        return delta_time
