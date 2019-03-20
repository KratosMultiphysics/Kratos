from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys

# Importing the Kratos Library
import KratosMultiphysics

import KratosMultiphysics.CompressiblePotentialFlowApplication as KCPFApp

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver


def CreateSolver(model, custom_settings):
    return PotentialFlowSolver(model, custom_settings)

class PotentialFlowSolver(PythonSolver):
    def __init__(self, model, custom_settings):

        super(PotentialFlowSolver, self).__init__(model, custom_settings)

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        # Default settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "potential_flow_solver",
            "domain_size": -1,
            "model_part_name": "",
            "echo_level": 1,
            "relative_tolerance": 1e-5,
            "absolute_tolerance": 1e-9,
            "maximum_iterations": 1,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "calculate_solution_norm": false,
            "volume_model_part_name": "volume_model_part",
            "skin_parts":[],
            "no_skin_parts": [],
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "linear_solver_settings": {
                "solver_type": "amgcl",
                "max_iteration": 400,
                "gmres_krylov_space_dimension": 100,
                "smoother_type":"ilu0",
                "coarsening_type":"ruge_stuben",
                "coarse_enough" : 5000,
                "krylov_type": "lgmres",
                "tolerance": 1e-9,
                "verbosity": 3,
                "scaling": false
            }
        }""")

        # Overwrite the default settings with user-provided parameters.
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        model_part_name = self.settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        # Set the element and condition names for the replace settings
        self.element_name = "IncompressiblePotentialFlowElement"
        self.condition_name = "PotentialWallCondition"

        # Set the main model part
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # Construct the linear solvers
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def AddVariables(self):
        # Degrees of freedom
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.VELOCITY_POTENTIAL)
        self.main_model_part.AddNodalSolutionStepVariable(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL)

        # Kratos variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.VELOCITY_POTENTIAL, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KCPFApp.AUXILIARY_VELOCITY_POTENTIAL, self.main_model_part)

    def Initialize(self):
        move_mesh_flag = False
        time_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        self.incompressible_solution_stratety = KratosMultiphysics.ResidualBasedLinearStrategy(
            self.GetComputingModelPart(),
            time_scheme,
            self.linear_solver,
            self.settings["compute_reactions"].GetBool(),
            self.settings["reform_dofs_at_each_step"].GetBool(),
            self.settings["calculate_solution_norm"].GetBool(),
            move_mesh_flag)

        (self.incompressible_solution_stratety).SetEchoLevel(self.settings["echo_level"].GetInt())
        (self.incompressible_solution_stratety).Initialize()

    def Check(self):
        self.incompressible_solution_stratety.Check()

    def ImportModelPart(self):
        # We can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Replace default elements and conditions
            self._ReplaceElementsAndConditions()
            ## Executes the check and prepare model process
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self.main_model_part.SetBufferSize(self.GetMinimumBufferSize())

        if self._IsPrintingRank():
            KratosMultiphysics.Logger.PrintInfo("PotentialFlowSolver", "Model reading finished.")

    def GetMinimumBufferSize(self):
        return 1

    def GetComputingModelPart(self):
        if not self.main_model_part.HasSubModelPart("fluid_computational_model_part"):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.main_model_part.GetSubModelPart("fluid_computational_model_part")

    def InitializeSolutionStep(self):
        self.incompressible_solution_stratety.InitializeSolutionStep()

    def Predict(self):
        self.incompressible_solution_stratety.Predict()

    def SolveSolutionStep(self):
        self.incompressible_solution_stratety.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.incompressible_solution_stratety.FinalizeSolutionStep()

    def SetEchoLevel(self, level):
        self.incompressible_solution_stratety.SetEchoLevel(level)

    def Clear(self):
        self.incompressible_solution_stratety.Clear()

    def AdvanceInTime(self, current_time):
        raise Exception("AdvanceInTime is not implemented. Potential Flow simulations are steady state.")

    def _IsPrintingRank(self):
        return self._is_printing_rank

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
        ## TODO: Implement it in such a way that the check_and_prepare_model_process_fluid can be used
        volume_model_part_name = self.settings["volume_model_part_name"].GetString()
        volume_model_part = self.main_model_part.GetSubModelPart(volume_model_part_name)

        skin_parts = []
        skin_name_list = self.settings["skin_parts"]
        for i in range(skin_name_list.size()):
            skin_parts.append(self.main_model_part.GetSubModelPart(skin_name_list[i].GetString()))

        # Construct a model part which contains both the skin and the volume
        # Temporarily we call it "fluid_computational_model_part"
        self.main_model_part.CreateSubModelPart("fluid_computational_model_part")
        fluid_computational_model_part = self.main_model_part.GetSubModelPart("fluid_computational_model_part")
        fluid_computational_model_part.ProcessInfo = self.main_model_part.ProcessInfo

        for node in volume_model_part.Nodes:
            fluid_computational_model_part.AddNode(node, 0)
        for elem in volume_model_part.Elements:
            fluid_computational_model_part.AddElement(elem, 0)

        ##TODO: Do some gymnastics to have this done fast. - create an ordered list to be added
        list_of_ids = set()
        for part in skin_parts:
            for cond in part.Conditions:
                list_of_ids.add(cond.Id)

        fluid_computational_model_part.AddConditions(list(list_of_ids))

        # Verify the orientation of the skin
        tmoc = KratosMultiphysics.TetrahedralMeshOrientationCheck
        throw_errors = False
        flags = tmoc.NOT_COMPUTE_NODAL_NORMALS | tmoc.COMPUTE_CONDITION_NORMALS | tmoc.ASSIGN_NEIGHBOUR_ELEMENTS_TO_CONDITIONS 
        KratosMultiphysics.TetrahedralMeshOrientationCheck(self.main_model_part, throw_errors, flags).Execute()
