class ResponseFunctionBase(object):
    def __init__(self, identifier, project_parameters):
        self.identifier = identifier
        self.project_parameters = project_parameters
    def Initialize(self):
        pass
    def CalculateValue(self):
        raise NotImplementedError("CalculateValue needs to be implemented by the base class")
    def GetValue(self):
        raise NotImplementedError("GetValue needs to be implemented by the base class")
    def CalculateGradient(self):
        raise NotImplementedError("CalculateGradient needs to be implemented by the base class")
    def GetGradient(self):
        raise NotImplementedError("GetGradient needs to be implemented by the base class")
    def Finalize(self):
        pass

class AdjointResponseFunction(ResponseFunctionBase):
    def __init__(self, identifier, project_parameters, model_part = None):
        # TODO use **kwargs for optional model_part, primal_analysis, adjoint_analysis
        super(AdjointResponseFunction, self).__init__(identifier, project_parameters)

        # Create the primal solver
        ProjectParametersPrimal = Parameters( open(project_parameters["primal_settings"].GetString(),'r').read() )
        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal, model_part)

        # Create the adjoint solver
        ProjectParametersAdjoint = Parameters( open(project_parameters["adjoint_settings"].GetString(),'r').read() )
        ProjectParametersAdjoint["solver_settings"].AddValue("response_function_settings", project_parameters)
        self.adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        # TODO find out why it is not possible to use the same model_part
    def Initialize(self):
        self.primal_analysis.Initialize()
        self.adjoint_analysis.Initialize()

        # TODO should be created here, not in solver!
        self.response_function_utility = self.adjoint_analysis.GetSolver().response_function
    def CalculateValue(self):

        for node in self.adjoint_analysis.GetModelPart().Nodes:
            ref_node = self.primal_analysis.GetModelPart().Nodes[node.Id]
            node.X0 = ref_node.X0
            node.Y0 = ref_node.Y0
            node.Z0 = ref_node.Z0
            node.X = ref_node.X
            node.Y = ref_node.Y
            node.Z = ref_node.Z

        print("\n> Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        self.primal_analysis.InitializeTimeStep()
        self.primal_analysis.SolveTimeStep()
        self.primal_analysis.FinalizeTimeStep()
        print("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")
        startTime = timer.time()
        self.value = self.response_function_utility.CalculateValue(self.primal_analysis.GetModelPart())
        print("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")
    def GetValue(self):
        return self.value
    def CalculateGradient(self):
        print("\n> Starting adjoint analysis for response:", self.identifier)
        startTime = timer.time()
        self.adjoint_analysis.InitializeTimeStep()
        self.adjoint_analysis.SolveTimeStep()
        self.adjoint_analysis.FinalizeTimeStep()
        print("> Time needed for solving the adjoint analysis = ",round(timer.time() - startTime,2),"s")
    def GetGradient(self):
        gradient = {}
        for node in self.adjoint_analysis.GetModelPart().Nodes:
            grad = node.GetSolutionStepValue(SHAPE_SENSITIVITY)
            gradient[node.Id] = grad
        return gradient

        MeshControllerUtilities( self.primal_analysis.GetModelPart() ).SetDeformationVariablesToZero()
        MeshControllerUtilities( self.primal_analysis.GetModelPart() ).SetDeformationVariablesToZero()
        # TODO reset also ADJOINT_DISPLACEMENT and ADJOINT_ROTATION
    def Finalize(self):
        self.primal_analysis.Finalize()
        self.adjoint_analysis.Finalize()

class SimpleResponseFunctionWrapper(ResponseFunctionBase):
    def __init__(self, identifier, project_parameters, response_function_utility, model_part = None):
        super(SimpleResponseFunctionWrapper, self).__init__(identifier, project_parameters)

        # Create the primal solver
        ProjectParametersPrimal = Parameters( open(project_parameters["primal_settings"].GetString(),'r').read() )
        self.primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal, model_part)

        self.response_function_utility = response_function_utility
    def Initialize(self):
        self.primal_analysis.Initialize()
        self.response_function_utility.Initialize()
    def CalculateValue(self):
        print("\n> Starting primal analysis for response:", self.identifier)
        startTime = timer.time()
        self.primal_analysis.InitializeTimeStep()
        self.primal_analysis.SolveTimeStep()
        self.primal_analysis.FinalizeTimeStep()
        print("> Time needed for solving the primal analysis = ",round(timer.time() - startTime,2),"s")
        startTime = timer.time()
        self.response_function_utility.CalculateValue()
        print("> Time needed for calculating the response value = ",round(timer.time() - startTime,2),"s")
    def GetValue(self):
        return self.response_function_utility.GetValue()
    def CalculateGradient(self):
        self.response_function_utility.CalculateGradient()
    def GetGradient(self):
        return self.response_function_utility.GetGradient()
    def Finalize(self):
        self.primal_analysis.Finalize()

class MassResponseFunctionWrapper(ResponseFunctionBase):
    def __init__(self, identifier, project_parameters, response_function_utility):
        super(MassResponseFunctionWrapper, self).__init__(identifier, project_parameters)
        self.response_function_utility = response_function_utility
    def Initialize(self):
        self.response_function_utility.Initialize()
    def CalculateValue(self):
        self.response_function_utility.CalculateValue()
    def GetValue(self):
        return self.response_function_utility.GetValue()
    def CalculateGradient(self):
        self.response_function_utility.CalculateGradient()
    def GetGradient(self):
        return self.response_function_utility.GetGradient()
    def Finalize(self):
        pass
