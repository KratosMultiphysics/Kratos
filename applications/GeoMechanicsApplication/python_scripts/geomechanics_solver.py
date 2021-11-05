# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication

def CreateSolver(model, custom_settings):
    return GeoMechanicalSolver(model, custom_settings)

class GeoMechanicalSolver(PythonSolver):
    """The base class for geomechanics solvers.

    This class provides shared functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes can override the functions
    """
    def __init__(self, model, custom_settings):

        super().__init__(model, custom_settings)

        # # Overwrite the default settings with user-provided parameters.
        # self.settings.ValidateAndAssignDefaults(default_settings)
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # This will be changed once the Model is fully supported!
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
            self.solver_imports_model_part = False
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')

            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.solver_imports_model_part = True

        self.min_buffer_size = 2

        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Construction finished")

        # Set if the analysis is restarted
        if self.settings["model_import_settings"]["input_type"].GetString() == "rest":
            KratosMultiphysics.Logger.PrintInfo("geomechanics_solver", "is a restarted model")
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        else:
            KratosMultiphysics.Logger.PrintInfo("geomechanics_solver", "is not a restarted model")
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def AddVariables(self):
        # this can safely be called also for restarts, it is internally checked if the variables exist already
        # Variables for 1-phase types of calculations:
        # Add displacements.
        self._add_displacement_variables()

        # Add rotational variables
        self._add_rotational_variables()

        # Add dynamic variables
        self._add_dynamic_variables()

        # Variables for 2-phase types of calculations:
        ## Fluid Variables
        self._add_water_variables()

        ## smoothing variables
        self._add_smoothing_variables()

        # Add variables that the user defined in the ProjectParameters
        if (self.settings.Has("auxiliary_variables_list")):
            for i in range(self.settings["auxiliary_variables_list"].size()):
                variable_name = self.settings["auxiliary_variables_list"][i].GetString()
                variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        pass

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def ImportModelPart(self):
        """This function imports the ModelPart
        """
        if self.solver_imports_model_part:
            self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])

    def PrepareModelPart(self):
        """This function prepares the ModelPart for being used by the PythonSolver
        """
        pass

    def KeepAdvancingSolutionLoop(self, end_time):
        current_time_corrected = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        return current_time_corrected < end_time

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        pass

    def InitializeSolutionStep(self):
        pass

    def Predict(self):
        pass

    def SolveSolutionStep(self):
        pass

    def FinalizeSolutionStep(self):
        pass

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        current_time_corrected = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        new_time = current_time_corrected + dt
        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def ComputeDeltaTime(self):
        return self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()


    #### Specific internal functions ####

    def import_constitutive_laws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: importing constitutive law", materials_filename)
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    #### Private functions ####

    def _add_dynamic_variables(self):
        # For being consistent for Serial and Trilinos
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)

    def _add_displacement_variables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.TOTAL_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_CONTACT_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TANGENTIAL_CONTACT_STRESS)

    def _add_rotational_variables(self):
        if (self.settings.Has("rotation_dofs")):
            if self.settings["rotation_dofs"].GetBool():
                # Add specific variables for the problem (rotation dofs).
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
                self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

    def _add_water_variables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.WATER_PRESSURE)
        # Add reactions for the water pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.DT_WATER_PRESSURE)
        # Add variables for the water conditions
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NORMAL_FLUID_FLUX)
        # Add variables for the water conditions
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.HYDRAULIC_DISCHARGE)

    def _add_smoothing_variables(self):
        if (self.settings.Has("nodal_smoothing")):
            if(self.settings["nodal_smoothing"].GetBool() == True):
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_CAUCHY_STRESS_TENSOR)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_DAMAGE_VARIABLE)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_AREA)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_WIDTH)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_DAMAGE)

    def _add_dynamic_dofs(self):
        # For being consistent for Serial and Trilinos
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Z,self.main_model_part)
        if(self.settings["rotation_dofs"].GetBool()):
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Z,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Z,self.main_model_part)

