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
        '''

        # Control Parameters
        self.stage = 'formfinding - start'

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
        except RuntimeError:
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
        self.model_part_name = project_parameters["solver_settings"]["model_part_name"].GetString()
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

        # Design
        self.design_check = Design(self.analysis_simulation, self.model_part_name, self.materials)
        self.design_check.Initialize()
        self.design_check.MembraneDesignCheck()

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

        jsonFile = open(material_file, "r")
        data = json.load(jsonFile)
        jsonFile.close()

        if self.stage == 'formfinding - start':

            KratosMultiphysics.Logger.PrintInfo("::[FormfindingUpdate]:: ", "Setting density and E_Mod to zero for formfinding")

            self._StructuralMaterialsUpdate(self.model_part_name, data["properties"], "Membrane")
            self._StructuralMaterialsUpdate(self.model_part_name, data["properties"], "Cable", {"DENSITY": 0, "YOUNG_MODULUS": 1e-18})

        if self.stage == 'analysis - start':

            KratosMultiphysics.Logger.PrintInfo("::[AnalysisUpdate]:: ", "Re-instating density and E_Mod to zero for analysis")

            self._StructuralMaterialsUpdate(self.model_part_name, data["properties"], "Membrane", {"DENSITY": self.materials["Membrane"]["DENSITY"], "YOUNG_MODULUS": self.materials["Membrane"]["YOUNG_MODULUS"]})
            self._StructuralMaterialsUpdate(self.model_part_name, data["properties"], "Cable", {"DENSITY": self.materials["Cable"]["DENSITY"], "YOUNG_MODULUS": self.materials["Cable"]["YOUNG_MODULUS"]})
        
        jsonFile = open(material_file, "w+")
        jsonFile.write(json.dumps(data))
        jsonFile.close()
    
    def _StructuralMaterialsUpdate(self,
                                   mp_name : str = "Structure",
                                   materials : list = None,
                                   group : str = None,
                                   new_variables : dict = {"DENSITY": 0.0, "YOUNG_MODULUS": 0}):

        '''
        Workflow:
        (1) self.materials stores original material parameters state
        (2) Material parameters within data gets manipulated depending on state
        (3) Both reference same dictionary object -> deepcopy required to break
        the link!

        Notes:
        Since _ManageMaterialsForSimulation() repeats the below code block several
        times, this was wrapped in a seperate function.
        '''

        # Break link to dictionary object (for manipulation purposes)
        from copy import deepcopy

        # Find element type
        if group == "Membrane":
            group_name = mp_name + ".Parts_Membrane"
        if group == "Cable":
            group_name = mp_name + ".Parts_Cable"

        # Iterate through defined materials                
        for material in materials:

            # If material belongs to input group_name
            if material.get("model_part_name").startswith(group_name):

                # Check if material has already been stored (provision for multiple iterations)
                try:
                    self.materials[group]

                # If it hasn't been stored, store it (Note: deepcopy, self.materials != material (depends on design phase))
                except KeyError:
                    self.materials[group] = deepcopy(material.get("Material").get("Variables"))

                # Update material (note, not self.materials) according to argument
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

class Design():

    def __init__(self, analysis_simulation, model_part_name, materials):

        # Analysis model
        self.model_part = analysis_simulation.model.GetModelPart(model_part_name)

        # ModelPart to be seperated into membrane and cables
        self.membrane_elements = [] # List of dictionaries
        self.cable_elements = [] # List of dictionaries

        # Materials
        self.materials = materials

        # Design parameters
        jsonFile = open("DesignParameters.json", "r")
        self.design_parameters = json.load(jsonFile)
        jsonFile.close()

    def Initialize(self):

        '''
        Workflow:
        (1) Seperate into membrane and cable elements
        (2) Generate dictionary of desired information
        (3) Append to attributes self.membrane_elements and self.cable_elements

        Notes:
        - Seperation needed because design checks differ per element type
        - Not relying on group model parts created by the User. Reliance on
        the groups may cause bugs (e.g. cable and membrane elements in one 
        group).
        
        Questions:
        ? Surely it gives a better way to seperate by ELEMENT_TYPE
        '''

        for index, element in enumerate(self.model_part.Elements):

            if "Triangle3D3" in element.GetSpecifications()["compatible_geometries"].GetStringArray():
                dict_to_append = self._PoppulateElementDictionary(index + 1, element)
                self.membrane_elements.append(dict_to_append)
            if "Line3D2" in element.GetSpecifications()["compatible_geometries"].GetStringArray():
                self._PoppulateElementDictionary(index + 1, element)
                self.membrane_elements.append(dict_to_append)
    
    def MembraneDesignCheck(self):

        self._CalculateDesignMembraneStress()
        self._CalculateMembraneResistance()
        self._CalculateMembraneUtilizations()
    
    def _PoppulateElementDictionary(self, index, obj):
            
            '''
            Notes:
            - setattr() and __setattr__() not permitted to an existing Kratos.Element
            - Workaround is to append a dictionary object as an attribute to the 
            Design() class object.
            - This function poppulates each dictionary instance 
            to a standard convention. This way, the membrane_element and 
            cable_elements can be poppulated with one line of code respectively.

            Workflow:
            (1) Get the stress representations (cauchy and principal)
            (2) Fill in template dictionary item
            '''

            stress_representations = self._CalculateStresses(obj)

            element_dict = {"index" : index,
                            "obj" : obj,
                            "cauchy_stress_vector" : {"Sxx"  : stress_representations[0][0],
                                                      "Syy"  : stress_representations[0][1],
                                                      "Sxy"  : stress_representations[0][2]},
                            "principal_stresses" :   {"Smax" : stress_representations[1][0],
                                                      "Smin" : stress_representations[1][1]}}

            return element_dict

    def _CalculateStresses(self, element):

        '''
        Notes:
        - Principal stresses conservatively used for design check

        Workflow:
        (1) Get cauchy_stress_vector from Kratos 
        (2) Calculate principal stresses
        (3) Return list of lists [[Sxx, Syy, Sxy], [Smax, Smin]]

        Questions:
        - Need to calculate principal stresses? Maybe == STRESSES from Kratos
        '''

        # Get CAUCHY_STRESS_VECTOR
        cauchy_stress_vector = element.CalculateOnIntegrationPoints(KratosMultiphysics.CAUCHY_STRESS_VECTOR, self.model_part.ProcessInfo)[0]

        # Return is vector of "Sxx", "Syy" and "Sxy"
        Sxx = cauchy_stress_vector[0]
        Syy = cauchy_stress_vector[1]
        Sxy = cauchy_stress_vector[2]

        # Calculate principal stresses using above
        Smax = ((Sxx + Syy) / 2) + (((Sxx - Syy) / 2)**2 + Sxy**2)**0.5
        Smin = ((Sxx + Syy) / 2) - (((Sxx - Syy) / 2)**2 + Sxy**2)**0.5
        principal_stresses = [Smax, Smin]

        # Return a nested list of two stress representations
        return [cauchy_stress_vector, principal_stresses]

    def _CalculateDesignMembraneStress(self):

        '''
        Notes:
        - Convert stress to a kN/m value (over the thickness)

        Workflow:
        (1) Take maximum principal stress (see _CalculateStresses())
        (2) Multiply by 0.001 to convert Pa to kN/m²
        (3) Multiply by thickness to get kN/m value
        (4) Add fEd to dictionary
        '''

        for element in self.membrane_elements:

            fEd = element["principal_stresses"]["Smax"] * 0.001 * self.materials["Membrane"]["THICKNESS"]
            element["fEd"] = fEd

    def _CalculateMembraneResistance(self):

        '''
        Notes:
        - Calculate membrane resistance § prCEN/TS 19102:2021 (E), 8.2.1
        - Constrained to Type II or Type III with warp == weft

        Workflow:
        (1) Check assumption regarding type
        (2) Characteristic tensile strength read from .json file (User-Input)
        (3) Partial factor(s) read from .json file (User-Input)
        (4) Design condition read from .json file (User-Input)
            (4.1) Strings as per prCEN/TS 19102:2021 (E), Table 8.2
            (4.2) Strings exactly as per table (replace spaces with "_")
        (5) Perform product of modification factors according to design condition
        (6) Calculate design tensile strength and write to class parameter
        '''

        # Assumption check
        if self.design_parameters["fabric_material"]["type"] not in ["II", "III"]:
            raise NotImplementedError("Design checks restricted to membranes of type II and III.")

        # Characteristic tensile strength of a membrane in a uniaxial short term tensile test at T = 23°C [prCEN/TS 19102:2021 (E), 3.2.2]
        fk23 = self.design_parameters["fabric_material"]["tensile_strength"]
        
        # The resistance of the material (γM0) is in this case the γM [prCEN/TS 19102:2021 (E), 8.1(3) and 8.2.1(3)]
        gamma_M = self.design_parameters["partial_factors"]["gamma_M"]

        # Design condition and respective modification factors [prCEN/TS 19102:2021 (E), Table 8.2]
        design_condition = self.design_parameters["general"]["design_condition"]
        modification_factors = self.design_parameters["modification_factors"]
        if design_condition == "prestress":
            k = modification_factors["k_biax"] * modification_factors["k_age"] * modification_factors["k_durP"] * modification_factors["k_temp"] * modification_factors["k_size"]
        elif design_condition == "prestress_temporarily_increased":
            k = modification_factors["k_biax"] * modification_factors["k_age"] * modification_factors["k_durL"] * modification_factors["k_temp"] * modification_factors["k_size"]
        elif design_condition == "snow_>_1000_m_altitude":
            k = modification_factors["k_biax"] * modification_factors["k_age"] * modification_factors["k_durL"] * modification_factors["k_size"]
        elif design_condition == "snow_<=_1000_m_altitude":
            k = modification_factors["k_biax"] * modification_factors["k_age"] * modification_factors["k_durM"] * modification_factors["k_size"]
        elif design_condition == "wind":
            k = modification_factors["k_biax"] * modification_factors["k_age"] * modification_factors["k_size"]
        elif design_condition == "wind_at_elevated_temperature":
            k = modification_factors["k_biax"] * modification_factors["k_age"] * modification_factors["k_temp"] * modification_factors["k_size"]
        else:
            available_design_conditions = ["prestress", "prestress_temporarily_increased", "snow_>_1000_m_altitude", "snow_<=_1000_m_altitude", "wind", "wind_at_elevated_temperature"]
            err_msg = "Available options are: {}".format(', '.join(available_design_conditions))
            raise Exception(err_msg)
    
        # Design tensile strength of fabrice material [prCEN/TS 19102:2021 (E), Formula (8.2)]
        self.design_tensile_strength = fk23 / (gamma_M * k)

    def _CalculateMembraneUtilizations(self):

        '''
        Notes:
        - Check utilization for each memrane element

        Workflow:
        (1) Iterate over list of membrane elements
        (2) Calculate utilization § prCEN/TS 19102:2021 (E), Formula (8.1)
        (3) Write utilization to element parameter
        '''

        for element in self.membrane_elements:
            element["utilization"] = element["fEd"] / self.design_tensile_strength
