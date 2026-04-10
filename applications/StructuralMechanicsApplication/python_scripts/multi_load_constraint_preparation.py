import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.process_factory import KratosProcessFactory
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.scipy_conversion_tools #just for debugging

class MultiLoadConstraintPreparation(AnalysisStage):

    def __init__(self, model, project_parameters):
        self.model = model
        self.settings = project_parameters["solver_settings"]
        self.scheme_settings = project_parameters["future_scheme_settings"]
        self.model_part_name = self.settings["model_part_name"].GetString()

        if self.model_part_name == "":
            raise Exception('Please specify a model_part name!')
        
        #check if model already has model part, otherwise create it and set domain size
        if self.model.HasModelPart(self.model_part_name):
            self.main_model_part = self.model[self.model_part_name]
        else:
            self.main_model_part = self.model.CreateModelPart(self.model_part_name)
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        super().__init__(model, project_parameters)
        if self.project_parameters.Has("process_combinations"):
            if self.project_parameters["process_combinations"].Has("load_processes"):
                self.load_processes = self.project_parameters["process_combinations"]["load_processes"]

    def Initialize(self):
        super().Initialize()
        self.__BuildLHS()
        #self.ApplyBoundaryConditions() -> this is now done in buildrhss
        self.__BuildRHSs()
        
    def __BuildRHSs(self):

        #stores right hand sides
        self.rhss = {}

        for load_definition in self.load_processes.values():
            self.__ResetRHS()
            process_id = load_definition["process_id"].GetString() #get the loadcase id
            process = self.__CreateProcess(load_definition) 
            process.ExecuteInitialize()
            process.ExecuteBeforeSolutionLoop()
            process.ExecuteInitializeSolutionStep() #usually done in ApplyBoundaryConditions()

            rhs = self.__BuildRHS()
            self.rhss[process_id] = rhs.copy()
        KratosMultiphysics.Logger.PrintInfo("::[PreparationStage]:: ", "RHSs built")

    def __ResetRHS(self):
        """Bug was: 
                [6](0,0,0,0,0,30000)
                [6](0,0,0,0,0,30000)
                lc3: SystemVector


                [6](40000,0,0,0,0,30000)
                [6](40000,0,0,0,0,30000)
                lc4: SystemVector
                
                instead of
                [6](0,0,0,0,0,30000)
                [6](0,0,0,0,0,30000)
                lc3: SystemVector


                [6](40000,0,0,0,0,0)
                [6](40000,0,0,0,0,0)
                lc4: SystemVector

                So reset of conditions is needed. Maybe there is a better way (i.e. Kratos already has that functionality)?
                """
        variables = [
            StructuralMechanicsApplication.POINT_LOAD,
            StructuralMechanicsApplication.LINE_LOAD,
            StructuralMechanicsApplication.SURFACE_LOAD,
            StructuralMechanicsApplication.POINT_MOMENT,
            KratosMultiphysics.POSITIVE_FACE_PRESSURE,
            KratosMultiphysics.NEGATIVE_FACE_PRESSURE,
        ]
        zero_vector = KratosMultiphysics.Vector(3)
        zero_vector[0] = 0.0
        zero_vector[1] = 0.0
        zero_vector[2] = 0.0

        for condition in self.main_model_part.Conditions:
            for variable in variables:
                if condition.Has(variable):
                    condition.SetValue(variable, zero_vector)

    def __BuildRHS(self):
        linear_system = self.strategy_data.GetLinearSystem()
        rhs = linear_system.GetVector(KratosMultiphysics.Future.DenseVectorTag.RHS)

        rhs.SetValue(0.0)
        self.scheme.Build(rhs)
        return rhs.copy()

    def __CreateProcess(self, process_definition):
        factory = KratosProcessFactory(self.model) #instantiate factory
        # the function "ConstructListOfProcesses" expects a Parameters array, so the process_list is stores in one here
        process_list = KratosMultiphysics.Parameters("[]") 
        process_list.Append(process_definition) 
        return factory.ConstructListOfProcesses(process_list)[0] # function returns a list, but we create processes one by one so it only holds one entry


    def RunSolutionLoop(self):
        """Changed it because only matrix generation is done"""
        self.InitializeSolutionStep()
        self.SolveSolutionStep()
    
    def Check(self):
        pass

    def __BuildLHS(self):
        """Construct the lhs"""
        #get linear system container stored inside the strategy data
        linear_system = self.strategy_data.GetLinearSystem()
        #Get the matrix stored under the tag LHS
        #Can store MassMatrix, DampingMatrix etc.
        lhs = linear_system.GetMatrix(KratosMultiphysics.Future.SparseMatrixTag.LHS)
        lhs.SetValue(0.0)
        self.scheme.Build(lhs)
        self.lhs = lhs
        KratosMultiphysics.Logger.PrintInfo("::[PreparationStage]:: ", "LHS built")
        #KratosMultiphysics.Logger.PrintInfo(lhs)
        print(KratosMultiphysics.scipy_conversion_tools.to_csr(lhs).todense())

    def _InitializeInternals(self):
        self.__InitializeScheme()
    
    def __InitializeScheme(self):
        self.strategy_data = KratosMultiphysics.Future.ImplicitStrategyData() #stores all arrays required to setup linear system
        self.scheme = KratosMultiphysics.Future.StaticScheme(self.main_model_part, self.scheme_settings) # scheme settings from json file
        self.scheme.Initialize(self.strategy_data) #initialize the scheme

    def GetComputingModelPart(self):
        """Is used in Initialize()"""
        return self.main_model_part

    def _AddDofs(self):
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append(["DISPLACEMENT_X", "REACTION_X"])
        dofs_and_reactions_to_add.append(["DISPLACEMENT_Y", "REACTION_Y"])
        dofs_and_reactions_to_add.append(["DISPLACEMENT_Z", "REACTION_Z"])
        if self.settings["rotation_dofs"].GetBool():
            dofs_and_reactions_to_add.append(["ROTATION_X", "REACTION_MOMENT_X"])
            dofs_and_reactions_to_add.append(["ROTATION_Y", "REACTION_MOMENT_Y"])
            dofs_and_reactions_to_add.append(["ROTATION_Z", "REACTION_MOMENT_Z"])

        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[PreparationStage]:: ", "DOF's ADDED")

    def _AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.SURFACE_LOAD)
        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs).
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)

    def _PrepareModelPart(self):
        self.__ReadMaterials()

    def __ReadMaterials(self):
        """Reads materials. Was done by solver before in 'PrepareModelPart'"""
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if materials_filename != "":
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}}""")
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            KratosMultiphysics.Logger.PrintInfo("::[PreparationStage]:: ", "Materials successfully imported")
        else:
            raise Exception("Please specify a 'materials_filename'!")
        
    def GetFinalData(self):
        """Returns the final data dictionary.

        The main purpose of this function is to retrieve any data (in a key-value format) from outside the stage.
        Note that even though it can be called at any point, it is intended to be called at the end of the stage run.
        """
        return {"lhs": self.lhs,
                "rhss": self.rhss,
                "model_part_name": self.main_model_part.FullName()}

if __name__ == "__main__":
    from sys import argv

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = MultiLoadConstraintPreparation(model, parameters)
    simulation.Run()