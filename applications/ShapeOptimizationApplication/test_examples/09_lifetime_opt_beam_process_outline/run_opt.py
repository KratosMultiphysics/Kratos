# Import Kratos core and apps
import KratosMultiphysics as Kratos
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape
import KratosMultiphysics.StructuralMechanicsApplication as KratosCSM
import structural_response_function_factory

# Additional imports
from analyzer_base import AnalyzerBaseClass

# Read parameters
with open("optimization_parameters.json",'r') as parameter_file:
    parameters = Kratos.Parameters(parameter_file.read())

model = Kratos.Model()

# Definition of external analyzer
class CustomAnalyzer(AnalyzerBaseClass):

    def InitializeBeforeOptimizationLoop(self):
        opt_model_part = model.GetModelPart("solid_cantilever")

        material_settings = Kratos.Parameters("""{"Parameters": {"materials_filename": ""}} """)
        material_settings["Parameters"]["materials_filename"].SetString("materials.json")
        Kratos.ReadMaterialsUtility(material_settings, model)

        for elem in opt_model_part.Elements:
            elem.Initialize()

    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        # Initialization of response
        tmp_model = Kratos.Model()

        kratos_response_settings = parameters["optimization_settings"]["objectives"][0]["kratos_response_settings"]
        response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, tmp_model)

        primal_model_part = response.primal_model_part
        adjoint_model_part = response.adjoint_model_part

        adjoint_model_part.AddNodalSolutionStepVariable(KratosCSM.DFDX)
        adjoint_model_part.AddNodalSolutionStepVariable(KratosCSM.DFDU)

        response.Initialize()

        # Update Mesh from optimizer
        mdpa_from_optimizer = model.GetModelPart(parameters["optimization_settings"]["model_settings"]["model_part_name"].GetString())
        for primal_node, adjoint_node, opt_node in zip(primal_model_part.Nodes, adjoint_model_part.Nodes, mdpa_from_optimizer.Nodes):

            X_opt = opt_node.X0
            Y_opt = opt_node.Y0
            Z_opt = opt_node.Z0

            primal_node.X = X_opt
            primal_node.X0 = X_opt

            primal_node.Y = Y_opt
            primal_node.Y0 = Y_opt

            primal_node.Z = Z_opt
            primal_node.Z0 = Z_opt

            adjoint_node.X = X_opt
            adjoint_node.X0 = X_opt

            adjoint_node.Y = Y_opt
            adjoint_node.Y0 = Y_opt

            adjoint_node.Z = Z_opt
            adjoint_node.Z0 = Z_opt

        # Calculate first primal field
        response.InitializeSolutionStep()
        response.CalculateValue()

        # Now, output displacemnets
        # for node in ....

        # Calculate lifetime

        # Calculate lifetime derivatives

        # Read value and derivative

        value = 10*optimization_iteration

        # External input
        for node in adjoint_model_part.Nodes:
            node.SetSolutionStepValue(KratosCSM.DFDU,[0,0,0])
            node.SetSolutionStepValue(KratosCSM.DFDX,[0,66,0])

        # Report values
        if communicator.isRequestingValueOf("lifetime"):
            communicator.reportValue("lifetime", value)

        if communicator.isRequestingGradientOf("lifetime"):
            response.CalculateGradient()
            communicator.reportGradient("lifetime", response.GetShapeGradient())

        # kratos_response_settings = parameters["optimization_settings"]["objectives"][0]["kratos_response_settings"]
        # response = structural_response_function_factory.CreateResponseFunction(kratos_response_settings["response_type"].GetString(), kratos_response_settings, tmp_model)
        # response.Initialize()
        # response.InitializeSolutionStep()
        # response.CalculateValue()
        # response.CalculateGradient()
        # response.FinalizeSolutionStep()
        # response.Finalize()

# Create optimizer and perform optimization
import optimizer_factory
optimizer = optimizer_factory.CreateOptimizer(parameters["optimization_settings"], model, CustomAnalyzer())
optimizer.Optimize()