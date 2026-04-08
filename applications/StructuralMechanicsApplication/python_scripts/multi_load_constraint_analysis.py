import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage

class MultiLoadConstraintAnalysis(AnalysisStage):

    def __init__(self, model, project_parameters):
        self.model = model
        self.settings = project_parameters["solver_settings"]
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

    def _GetListOfProcesses(self):
        """Otherwise error, because Initialize() of the baseclass loops through self._list_of_processes and calls a method for each process.
        I would have to add variables for this to work and I think this should not be done right now."""
        return []
        
    def RunSolutionLoop(self):
        """Changed it because only matrix generation is done"""
        self.InitializeSolutionStep()
        self.SolveSolutionStep()

    def InitializeSolutionStep(self):
        #not necessary, but why not
        self.PrintAnalysisStageProgressInformation()

    def SolveSolutionStep(self):
        self.BuildLHS()

    def BuildLHS(self):
        """Construct the lhs"""
        linear_system = self.strategy_data.GetLinearSystem()
        lhs = linear_system.GetMatrix(KratosMultiphysics.Future.SparseMatrixTag.LHS)

        lhs.SetValue(0.0)
        KratosMultiphysics.Logger.PrintInfo(lhs)
        self.scheme.Build(lhs)
        KratosMultiphysics.Logger.PrintInfo(lhs)

    def _InitializeInternals(self):
        self.InitializeScheme()
    
    def InitializeScheme(self):
        self.strategy_data = KratosMultiphysics.Future.ImplicitStrategyData()
        scheme_settings = self.project_parameters["future_scheme_settings"]
        self.scheme = KratosMultiphysics.Future.StaticScheme(self.main_model_part, scheme_settings)
        self.scheme.Initialize(self.strategy_data)

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
        KratosMultiphysics.Logger.PrintInfo("::[Null]:: ", "DOF's ADDED")

    def _PrepareModelPart(self):
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] = 0
        self.main_model_part.ProcessInfo[KratosMultiphysics.TIME] = self.project_parameters["problem_data"]["start_time"].GetDouble()
        self.read_materials()

    def read_materials(self):
        """Reads materials. Was done by solver before in 'PrepareModelPart'"""
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if materials_filename != "":
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}}""")
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            KratosMultiphysics.Logger.PrintInfo("::[Null]:: ", "Materials successfully imported")
        else:
            raise Exception("Please specify a 'materials_filename'!")

if __name__ == "__main__":
    from sys import argv

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = MultiLoadConstraintAnalysis(model, parameters)
    simulation.Run()