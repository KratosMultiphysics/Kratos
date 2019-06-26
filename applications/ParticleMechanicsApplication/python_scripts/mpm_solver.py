from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import structural mechanics applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
# Import base class file from structural mechanics application
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(model, custom_settings):
    return MpmSolver(model, custom_settings)

class MpmSolver(MechanicalSolver):
    """The MPM solver is based on the structural mechanics solver.

    Structural mechanics application has to be compiled to run this solver.
    """
    ### Solver constructor
    def __init__(self, model, custom_settings):
        super(MpmSolver, self).__init__(model, custom_settings)

        # Add model part containers
        self._add_model_part_containers()

        self.pressure_dofs = self.settings["pressure_dofs"].GetBool()
        self.axis_symmetric_flag = self.settings["axis_symmetric_flag"].GetBool()

        KratosMultiphysics.Logger.PrintInfo("::[MpmSolver]:: ", "Construction finished.")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "grid_model_import_settings": {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name_Grid"
            },
            "axis_symmetric_flag"                : false,
            "pressure_dofs"                      : false,
            "element_search_settings"            : {
                "max_number_of_results"          : 1000,
                "searching_tolerance"            : 1.0E-5
            }
        }""")
        this_defaults.AddMissingParameters(super(MpmSolver, cls).GetDefaultSettings())
        return this_defaults

    ### Solver public functions
    def AddVariables(self):
        # Add variables to different model parts
        self._add_variables_to_model_part(self.grid_model_part)
        #self._add_variables_to_model_part(self.main_model_part)
        #self._add_variables_to_model_part(self.initial_material_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[MpmSolver]:: ", "Variables are added.")

    def ImportModelPart(self):
        # Read model part
        self._model_part_reading()

                # Identify geometry type
        self._identify_geometry_type()

        # Set default solver_settings parameters
        if self.geometry_element == "Triangle":
            if (self.domain_size == 2):
                if (self.pressure_dofs):
                    self.new_element = KratosParticle.CreateUpdatedLagragianUP2D3N()
                else:
                    if (self.axis_symmetric_flag):
                        self.new_element = KratosParticle.CreateUpdatedLagragianAxis2D3N()
                    else:
                        self.new_element = KratosParticle.CreateUpdatedLagragian2D3N()
            else:
                if (self.pressure_dofs):
                    raise Exception("Element for mixed U-P formulation in 3D for Tetrahedral Element is not yet implemented.")
                else:
                    self.new_element = KratosParticle.CreateUpdatedLagragian3D4N()
        elif self.geometry_element == "Quadrilateral":
            if (self.domain_size == 2):
                if (self.pressure_dofs):
                    raise Exception("Element for mixed U-P formulation in 2D for Quadrilateral Element is not yet implemented.")
                else:
                    if (self.axis_symmetric_flag):
                        self.new_element = KratosParticle.CreateUpdatedLagragianAxis2D4N()
                    else:
                        self.new_element = KratosParticle.CreateUpdatedLagragian2D4N()
            else:
                if (self.pressure_dofs):
                    raise Exception("Element for mixed U-P formulation in 3D for Hexahedral Element is not yet implemented.")
                else:
                    self.new_element = KratosParticle.CreateUpdatedLagragian3D8N()

        KratosParticle.CreateMaterialPointElement(
            self.grid_model_part, self.main_model_part, self.initial_material_model_part, self.new_element, False)
        KratosParticle.CreateMaterialPointCondition(
            self.grid_model_part, self.main_model_part, self.initial_material_model_part)

        for prop in self.main_model_part.Properties:
            print(prop)

        print(self.grid_model_part)
        print(self.main_model_part)
        print(self.initial_material_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[MpmSolver]:: ","Models are imported.")

    def _identify_geometry_type(self):
        for mpm in self.grid_model_part.Elements:
            num_nodes = len(mpm.GetNodes())
            break

        self.domain_size = self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)

        if (self.domain_size == 2 and num_nodes == 3) or (self.domain_size == 3 and num_nodes == 4):
            self.geometry_element = "Triangle"
        elif (self.domain_size == 2 and num_nodes == 4) or (self.domain_size == 3 and num_nodes == 8):
            self.geometry_element = "Quadrilateral"

    def AddDofs(self):
        # Add dofs to different model parts
        self._add_dofs_to_model_part(self.grid_model_part)
        #self._add_dofs_to_model_part(self.main_model_part)
        #self._add_dofs_to_model_part(self.initial_material_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[MpmSolver]:: ","DOFs are added.")

    def SearchElement(self):
        KratosParticle.SearchElement2D(self.grid_model_part, self.main_model_part, 1000, 1.0e-5)

    def InitializeSolutionStep(self):
        self.SearchElement()
        print(self.main_model_part)
        self.get_mechanical_solution_strategy().Initialize()
        self.get_mechanical_solution_strategy().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.get_mechanical_solution_strategy().FinalizeSolutionStep()

        self.get_mechanical_solution_strategy().Clear()

    ### Solver private functions
    def _add_model_part_containers(self):

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size not in [2,3]:
            err_msg  = "The input \"domain_size\" is wrong!"
            err_msg += "Available options are: \"2\" or \"3\""
            raise Exception(err_msg)

        ### In MPM three model parts are needed
        ## Material model part definition
        material_model_part_name = self.settings["model_part_name"].GetString()
        if not self.model.HasModelPart(material_model_part_name):
            self.main_model_part = self.model.CreateModelPart(material_model_part_name) # Equivalent to model_part3 in the old format
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Initial material model part definition
        initial_material_model_part_name = "Initial_" + material_model_part_name
        if not self.model.HasModelPart(initial_material_model_part_name):
            self.initial_material_model_part = self.model.CreateModelPart(initial_material_model_part_name) #Equivalent to model_part2 in the old format
            self.initial_material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Grid model part definition
        if not self.model.HasModelPart("Background_Grid"):
            self.grid_model_part = self.model.CreateModelPart("Background_Grid") #Equivalent to model_part1 in the old format
            self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    def _add_variables_to_model_part(self, model_part):
        # Add displacements
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)

        # Add dynamic variables
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)

        # Add reactions for the displacements
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # Add specific variables for the problem conditions
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MOMENTUM)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_INERTIA)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)

        # Add variables for arbitrary slope with slip
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.IS_STRUCTURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)

        # Add variables for specific cases
        if self.settings["pressure_dofs"].GetBool():
            # add specific variables for the problem (pressure dofs)
            model_part.AddNodalSolutionStepVariable(KratosParticle.PRESSURE_REACTION)
            model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MPRESSURE)

    def _model_part_reading(self):
        # reading the model part of the background grid
        if(self.settings["grid_model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.grid_model_part, self.settings["grid_model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

        # reading the model part of the material point
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.initial_material_model_part, self.settings["model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

    def _execute_after_reading(self):
        # Specific active node and element check for particle MPM solver
        for node in self.grid_model_part.Nodes:
            if (node.Is(KratosMultiphysics.ACTIVE)):
                KratosMultiphysics.Logger.PrintInfo("::[MpmSolver]:: ","WARNING: This grid node have been set active: ", node.Id)

        # Setting active initial elements
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.initial_material_model_part.Elements)

        # Specify domain size
        self.domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        """Prepare computing model part and import constitutive laws. """
        # Auxiliary parameters object for the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",self.settings["model_part_name"])
        params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        # Assign mesh entities from domain and process sub model parts to the computing model part.
        import check_and_prepare_model_process_structural
        check_and_prepare_model_process_structural.CheckAndPrepareModelProcess(self.model, params).Execute()

        # Read material property
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[MpmSolver]:: ","Constitutive law was successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintWarning("::[MpmSolver]:: ","Constitutive law was not imported.")

        # Clone property of model_part2 to model_part3
        self.main_model_part.Properties = self.initial_material_model_part.Properties

    def _add_dofs_to_model_part(self, model_part):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, model_part)

        if self.settings["pressure_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosParticle.PRESSURE_REACTION, model_part)
