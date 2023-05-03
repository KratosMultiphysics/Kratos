# Importing Kratos
import KratosMultiphysics

# Importing the base class
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis

# Importing the formfinding solvers
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_formfinding_solver import FormfindingMechanicalSolver
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_static_solver import StaticMechanicalSolver

# Importing third-party libraries
import json

class StructuralMechanicsFormfindingAndAnalysis(StructuralMechanicsAnalysis):

    def __init__(self, model, project_parameters):

        '''
        NOTES:
        - project_parameters assumed to be default "run-no write" from GiD for
        a formfinding analysis. That is why the CreateAnalysisSolver, for 
        example, needs to remove values. This can be refined at a later stage.
        - The classes _Formfinding, _Analysis work "in the background". 
        The StructuralMechanicsFormfindingAnalysis is in itself, no simulation,
        it has attributes of the various simulations (nested simulations)

        TO-DO:
        [ ] Design check evtl. invoke formfinding again with higher prestress
        '''

        # Control Parameters
        self.stage = 'formfinding - start'

        # Empty simluations (PyLint req. within __init__)
        self.formfinding_simulation = None
        self.analysis_simulation = None

        # Project Parameters
        self.formfinding_parameters = project_parameters.Clone()
        self.analysis_parameters = project_parameters.Clone()

        # Convergence Criterion

        '''
        In the event that the User desires different convergence_criterion 
        for formfinding and analysis, this can be achieved by creating the 
        following parameters in the ProjectParameters.json file:

                formfinding_convergence_criterion
                analysis_convergence_criterion

        This code block then sets the parameters accordingly and removes values 
        such that validation does not fail.
        '''
        try:
            project_parameters["solver_settings"]["formfinding_convergence_criterion"].GetString()
            project_parameters["solver_settings"]["analysis_convergence_criterion"].GetString()
        except:
            KratosMultiphysics.Logger.PrintInfo("::[FormfindingAndAnalysis]:: ", "Using standard convergence criterion for both formfinding and analysis")
        else:
            KratosMultiphysics.Logger.PrintInfo("::[FormfindingAndAnalysis]:: ", "Different convergence criterion detected, applying appropriately to solver_settings")
            
            self.formfinding_parameters["solver_settings"]["convergence_criterion"].SetString(project_parameters["solver_settings"]["formfinding_convergence_criterion"].GetString())
            self.analysis_parameters["solver_settings"]["convergence_criterion"].SetString(project_parameters["solver_settings"]["analysis_convergence_criterion"].GetString())
            
            self.formfinding_parameters["solver_settings"].RemoveValue("formfinding_convergence_criterion")
            self.formfinding_parameters["solver_settings"].RemoveValue("analysis_convergence_criterion")
            
            self.analysis_parameters["solver_settings"].RemoveValue("formfinding_convergence_criterion")
            self.analysis_parameters["solver_settings"].RemoveValue("analysis_convergence_criterion")

        # Models
        self.input_model = model
        self.formfound_model = None
        self.analysed_model = None

        # Properties / Materials and loads
        self.materials = {}
        self.loads = []

    def Run(self):

        # Preprocessing
        self.ManageSimulationParameters()

        # Formfinding
        self.formfinding_simulation = Formfinding(self.input_model, self.formfinding_parameters)
        self.CheckProjectionSettings()
        self.LoadsAndMaterialsManager()
        self.formfinding_simulation.Initialize()
        self.formfinding_simulation.RunSolutionLoop()
        self.RefreshSimulation()

        # Analysis
        self.analysis_simulation = Analysis(self.formfound_model, self.analysis_parameters)
        self.LoadsAndMaterialsManager()
        self.analysis_simulation.Initialize()
        self.analysis_simulation.RunSolutionLoop()
        self.RefreshSimulation()

        self.formfinding_simulation.Finalize()
        self.analysis_simulation.Finalize()
    
    def ManageSimulationParameters(self):

        '''
        It is assumed, that the export from GiD is of that for formfinding. As
        such, parameters need to be adjusted / removed for analysis. If the 
        export does not align to this assumption, a NotImplementedError is
        raised.
        '''

        KratosMultiphysics.Logger.PrintInfo("::[FormfindingAndAnalysis]:: ", "Adjusting parameter files")

        if self.analysis_parameters["solver_settings"]["solver_type"].GetString() != "formfinding":
            raise NotImplementedError("Required manipulations are based off of an export off a model with Kratos Solver Type 'Formfinding'.")

        # The analysis uses result from formfinding
        self.analysis_parameters["solver_settings"]["model_import_settings"]["input_type"].SetString("use_input_model_part")
        self.analysis_parameters["solver_settings"]["model_import_settings"].RemoveValue("input_filename")

        # Some project_parameters for formdinding not needed for anaylsis and cause validation errors
        self.analysis_parameters["solver_settings"].RemoveValue('formfinding_model_part_name')
        self.analysis_parameters["solver_settings"].RemoveValue('printing_format')
        self.analysis_parameters["solver_settings"].RemoveValue('write_formfound_geometry_file')
        self.analysis_parameters["solver_settings"].RemoveValue('projection_settings')
        self.analysis_parameters["solver_settings"].RemoveValue('printing_format')

        # Time stepping - the default time_step_table was always manually deleted from ProjectParameters.json and replaced with time_step = 1.1
        try:
            self.formfinding_parameters["solver_settings"]["time_stepping"]["time_step"].GetDouble()

        except:
            KratosMultiphysics.Logger.PrintInfo("::[FormfindingAndAnalysis]:: ", "'time_step' manually set to 1.1. 'time_step_table' parameter removed.")
            self.formfinding_parameters["solver_settings"]["time_stepping"].RemoveValue("time_step_table")
            self.formfinding_parameters["solver_settings"]["time_stepping"].AddDouble("time_step", 1.1)
            self.analysis_parameters["solver_settings"]["time_stepping"].RemoveValue("time_step_table")
            self.analysis_parameters["solver_settings"]["time_stepping"].AddDouble("time_step", 1.1)
        
        # list_other_processes default from GiD causes ModuleNotFoundError, remove
        self.formfinding_parameters["processes"].RemoveValue("list_other_processes")
        self.analysis_parameters["processes"].RemoveValue("list_other_processes")

    def CheckProjectionSettings(self):

        '''
        With the export from GiD, the projection_settings in the json 
        ProjectParameters file is always manually updated. Avoiding this step 
        by setting the required values which the ValidateAndAssignDefaults() 
        doesn't do.
        '''

        global_echo_level = self.formfinding_simulation.project_parameters["solver_settings"]["echo_level"].GetDouble()
        self.formfinding_simulation.project_parameters["solver_settings"]["projection_settings"]["echo_level"].SetDouble(global_echo_level)
        
        if self.formfinding_simulation.project_parameters["solver_settings"]["projection_settings"]["variable_name"].GetString() != "LOCAL_PRESTRESS_AXIS_1":
            self.formfinding_simulation.project_parameters["solver_settings"]["projection_settings"]["variable_name"].SetString("LOCAL_PRESTRESS_AXIS_1")
            KratosMultiphysics.Logger.PrintInfo("::[FormfindingUpdate]:: ", "Projection setting adjusted to LOCAL_PRESTRESS_AXIS_1")

    def LoadsAndMaterialsManager(self):

        '''
        Invokes a series of sub-processes. Refer to sub-processes for details.
        '''

        self._ManageLoadsForSimulation()
        self._ManageMaterialsForSimulation()

    def _ManageLoadsForSimulation(self):

        '''
        For formfinding this means:
            - Remove all loads

        For analysis this means:
            - Re-instate all loads

        This is done by manipulating the loads_process_list of the Kratos 
        parameters object.
        '''

        if self.stage == "formfinding - start":

            KratosMultiphysics.Logger.PrintInfo("::[FormfindingUpdate]:: ", "Removing external loads for formfinding")
            
            for load in self.formfinding_parameters["processes"]["loads_process_list"]:
                self.loads.append(load["Parameters"]["modulus"].GetDouble())
                load["Parameters"]["modulus"].SetDouble(0.0)
        
        if self.stage == "analysis - start":

            KratosMultiphysics.Logger.PrintInfo("::[AnalysisUpdate]:: ", "Re-instating external loads for analysis")

            for index, load in enumerate(self.analysis_parameters["processes"]["loads_process_list"]):
                load["Parameters"]["modulus"].SetDouble(self.loads[index])

    def _ManageMaterialsForSimulation(self):

        '''
        For formfinding this means:
            - Setting all density and E to zero (or mathematically zero)

        For analysis this means:
            - Re-instating all originally defined values

        This is done by manipulating the json file. This is simpler than editing
        the Parameters object - as for this, one would have to create extra 
        functions to link the property id with the respective groups.
        '''

        material_file = self.formfinding_parameters["solver_settings"]["material_import_settings"]["materials_filename"].GetString()
        mp_name = self.formfinding_parameters["solver_settings"]["model_part_name"].GetString()

        jsonFile = open(material_file, "r")
        data = json.load(jsonFile)
        jsonFile.close()

        if self.stage == 'formfinding - start':

            KratosMultiphysics.Logger.PrintInfo("::[FormfindingUpdate]:: ", "Setting density and E_Mod to zero for formfinding")

            self._StructuralMaterialsUpdate(mp_name, data["properties"], "Membrane")
            self._StructuralMaterialsUpdate(mp_name, data["properties"], "Cable", {"DENSITY": 0, "YOUNG_MODULUS": 1e-18})

        if self.stage == 'analysis - start':

            KratosMultiphysics.Logger.PrintInfo("::[AnalysisUpdate]:: ", "Re-instating density and E_Mod to zero for analysis")

            self._StructuralMaterialsUpdate(mp_name, data["properties"], "Membrane", {"DENSITY": self.materials["Membrane"]["DENSITY"], "YOUNG_MODULUS": self.materials["Membrane"]["YOUNG_MODULUS"]})
            self._StructuralMaterialsUpdate(mp_name, data["properties"], "Cable", {"DENSITY": self.materials["Cable"]["DENSITY"], "YOUNG_MODULUS": self.materials["Cable"]["YOUNG_MODULUS"]})
        
        jsonFile = open(material_file, "w+")
        jsonFile.write(json.dumps(data))
        jsonFile.close()
    
    def _StructuralMaterialsUpdate(self,
                                   mp_name : str = "Structure",
                                   materials : list = None,
                                   group : str = None,
                                   new_variables : dict = {"DENSITY": 0.0, "YOUNG_MODULUS": 0}):

        '''
        Since _ManageMaterialsForSimulation() repeats the below code block several
        times, this was wrapped in a seperate function. Easier to fix a bug 
        once, than four times!
        '''

        if group == "Membrane":
            group_name = mp_name + ".Parts_Membrane"
        if group == "Cable":
            group_name = mp_name + ".Parts_Cable"
        
        for material in materials:

            if material.get("model_part_name").startswith(group_name):

                if self.stage == 'formfinding - start':
                    
                    density = material.get("Material").get("Variables").get("DENSITY")
                    young_modulus = material.get("Material").get("Variables").get("YOUNG_MODULUS")

                    variables = dict(zip(["DENSITY", "YOUNG_MODULUS"], [density, young_modulus]))
                    self.materials[group] = variables

                material.get("Material").get("Variables").update(new_variables)

    def RefreshSimulation(self):

        '''
        Invokes a series of sub-processes.
        '''

        self._ManageModels()
        self._SwitchStages()

    def _ManageModels(self):

        '''
        Stores the results of various stages as parameters. The model objects
        are stored.
        '''

        if self.stage == 'formfinding - start':
            self.stage = "formfinding - end"
            self.formfound_model = self.formfinding_simulation.model
        if self.stage == 'analysis - start':
            self.stage = "analysis - end"
            self.analysed_model = self.analysis_simulation.model

    def _SwitchStages(self):

        '''
        Possibly seen as uneccesary, however, envisaged future usage for 
        Bemessungen / Nachweise can re-trigger a formfinding process?
        '''

        if self.stage == 'formfinding - end':
            KratosMultiphysics.Logger.PrintInfo("Switching simulation from formfinding to analysis")
            self.stage = 'analysis - start'
        if self.stage == 'analysis - end':
            self.stage = 'completed'

class Formfinding(StructuralMechanicsAnalysis):

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

    def _CreateSolver(self):
        return FormfindingMechanicalSolver(self.model, self.project_parameters['solver_settings'])

class Analysis(StructuralMechanicsAnalysis):

    def __init__(self, model, project_parameters):
        super().__init__(model, project_parameters)

    def _CreateSolver(self):
        return StaticMechanicalSolver(self.model, self.project_parameters['solver_settings'])
