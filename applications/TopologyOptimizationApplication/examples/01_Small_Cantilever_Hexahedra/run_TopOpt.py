#import kratos core and applications
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as ksm
import KratosMultiphysics.TopologyOptimizationApplication as kto
import KratosMultiphysics.LinearSolversApplication as kls
from KratosMultiphysics import process_factory
import os
import OptimizationParameters as opt_parameters
from importlib import import_module
from KratosMultiphysics.gid_output_process import GiDOutputProcess
import KratosMultiphysics.KratosUnittest as KratosUnittest

parameter_file = open("/home/philipp/opt/kratosDev/applications/TopologyOptimizationApplication/examples/01_Small_Cantilever_Hexahedra/ProjectParameters.json",'r')
ProjectParameters = km.Parameters(parameter_file.read())
echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
current_model = km.Model()
dimension = 2
model_part = current_model.CreateModelPart("Structure")
solver_module = ProjectParameters["solver_settings"]["solver_type"].GetString()
#mod = 'KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_static_solver'
mod = 'KratosMultiphysics.TopologyOptimizationApplication.topology_optimization_simp_static_solver'
solver = import_module(mod).CreateSolver(current_model, ProjectParameters["solver_settings"])
solver.AddVariables()
solver.ImportModelPart()
solver.AddDofs()
model_part.GetProperties()[1].SetValue(km.YOUNG_MODULUS, 200.0e9)
model_part.GetProperties()[1].SetValue(km.POISSON_RATIO, 0.4)
model_part.GetProperties()[1].SetValue(km.DENSITY, 1.0)

cons_law = ksm.LinearElastic3DLaw()
model_part.GetProperties()[1].SetValue(km.CONSTITUTIVE_LAW, cons_law)

if(echo_level>1):
    print(main_model_part)
    for properties in main_model_part.Properties:
        print(properties)

#obtain the list of the processes to be applied (the process order of execution is important)
list_of_processes  = process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses( ProjectParameters["processes"]["constraints_process_list"] )
list_of_processes += process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses( ProjectParameters["processes"]["loads_process_list"] )
if(ProjectParameters.Has("problem_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses( ProjectParameters["processes"]["problem_process_list"] )
if(ProjectParameters.Has("output_process_list")):
    list_of_processes += process_factory.KratosProcessFactory(current_model).ConstructListOfProcesses( ProjectParameters["processes"]["output_process_list"] )
            
if(echo_level>1):
    for process in list_of_processes:
        print(process)
for process in list_of_processes:
    process.ExecuteInitialize()

computing_model_part = solver.GetComputingModelPart()
problem_path = os.getcwd()
problem_name = ProjectParameters["problem_data"]["problem_name"].GetString()

# initialize GiD  I/O (gid outputs, file_lists)
#from gid_output_process import GiDOutputProcess
output_settings = ProjectParameters["processes"]["output_configuration"]
gid_output = GiDOutputProcess(computing_model_part,
                              problem_name,
                              output_settings)

gid_output.ExecuteInitialize()
solver.Initialize()
for process in list_of_processes:
    process.ExecuteBeforeSolutionLoop()
gid_output.ExecuteBeforeSolutionLoop()

""" # Solve the problem
linear_solver = km.SkylineLUFactorizationSolver()
builder_and_solver = km.ResidualBasedBlockBuilderAndSolver(linear_solver)
scheme = km.ResidualBasedIncrementalUpdateStaticScheme()
compute_reactions = True
reform_step_dofs = True
calculate_norm_dx = False
move_mesh_flag = True

strategy = km.ResidualBasedLinearStrategy(
    model_part,
    scheme,
    builder_and_solver,
    compute_reactions,
    reform_step_dofs,
    calculate_norm_dx,
    move_mesh_flag)

strategy.Solve() """

for process in list_of_processes:
    process.ExecuteInitializeSolutionStep()
gid_output.ExecuteInitializeSolutionStep()       
# Solve the problem
linear_solver = km.SkylineLUFactorizationSolver()
builder_and_solver = km.ResidualBasedBlockBuilderAndSolver(linear_solver)
scheme = km.ResidualBasedIncrementalUpdateStaticScheme()
compute_reactions = True
reform_step_dofs = True
calculate_norm_dx = False
move_mesh_flag = True

strategy = km.ResidualBasedLinearStrategy(
    model_part,
    scheme,
    builder_and_solver,
    compute_reactions,
    reform_step_dofs,
    calculate_norm_dx,
    move_mesh_flag)

strategy.Solve()
for process in list_of_processes:
    process.ExecuteFinalizeSolutionStep()
gid_output.ExecuteFinalizeSolutionStep()
for process in list_of_processes:
    process.ExecuteFinalizeSolutionStep()
for process in list_of_processes:
    process.ExecuteBeforeOutputStep()
if(gid_output.IsOutputStep()):
    gid_output.PrintOutput()
for process in list_of_processes:
    process.ExecuteAfterOutputStep()

response_analyzer = kto.StructureResponseFunctionUtilities(model_part)
linear_solver = km.python_linear_solver_factory.ConstructSolver(ProjectParameters["solver_settings"]["linear_solver_settings"])
sensitivity_solver = kto.StructureAdjointSensitivityStrategy(model_part, linear_solver,ProjectParameters["solver_settings"]["domain_size"].GetInt())

"""     # Compute structure solution to get displacement field u
    strategy.Solve()
    # Calculate objective function value based on u and save into container
    response["strain_energy"]["func"] = response_analyzer.ComputeStrainEnergy() """
# Compute constraint function value

"""     target_volume_fraction = opt_parameters.initial_volume_fraction
    response["volume_fraction"]["func"] = response_analyzer.ComputeVolumeFraction() - target_volume_fraction
# Compute sensitivities of objective function
if(controls["strain_energy"]["calc_grad"]):         """
sensitivity_solver.ComputeStrainEnergySensitivities()
""" # Compute sensitivities of constraint function
if(controls["volume_fraction"]["calc_grad"]):  
    sensitivity_solver.ComputeVolumeFractionSensitivities() """
















""" class RunTopOPt(KratosUnittest.TestCase):

    def _apply_material_properties(self, ModelPart, Dimesion):
         # Define material properties
        ModelPart.GetProperties()[1].SetValue(km.YOUNG_MODULUS, 200.0e9)
        ModelPart.GetProperties()[1].SetValue(km.POISSON_RATIO, 0.4)
        ModelPart.GetProperties()[1].SetValue(km.DENSITY, 1.0)

        # Define body force
        g = [0,0,0]
        ModelPart.GetProperties()[1].SetValue(km.VOLUME_ACCELERATION, g)

        # Define constitutive law
        ModelPart.GetProperties()[1].SetValue(km.CONSTITUTIVE_LAW, cons_law)

    def _solve(self, ModelPart):
        # Define a linear strategy to solve the problem
        linear_solver = km.SkylineLUFactorizationSolver()
        builder_and_solver = km.ResidualBasedBlockBuilderAndSolver(linear_solver)
        scheme = km.ResidualBasedIncrementalUpdateStaticScheme()
        compute_reactions = True
        reform_step_dofs = True
        calculate_norm_dx = False
        move_mesh_flag = True

        strategy = km.ResidualBasedLinearStrategy(
            ModelPart,
            scheme,
            builder_and_solver,
            compute_reactions,
            reform_step_dofs,
            calculate_norm_dx,
            move_mesh_flag)
        strategy.SetEchoLevel(0)
        strategy.Check()

        # Solve the problem
        strategy.Solve()

    def testSmallDisplacement(self):
        parameter_file = open("/home/philipp/opt/kratosDev/applications/TopologyOptimizationApplication/examples/01_Small_Cantilever_Hexahedra/ProjectParameters.json",'r')
        ProjectParameters = km.Parameters(parameter_file.read())
        echo_level = ProjectParameters["problem_data"]["echo_level"].GetInt()
        dimension = 2
        current_model = km.Model()
        model_part = current_model.CreateModelPart("Structure")
        model_part.AddVariables()
        self.AddDofs(model_part)
        self.ImportModelPart(model_part)
        self._apply_material_properties(model_part, dimension)

        self._solve(model_part)

    def __post_process(self, main_model_part, post_type = "gid"):
        if post_type == "gid":
            self.gid_output = GiDOutputProcess(
                main_model_part,
                main_model_part.Name,
                km.Parameters(r

            self.gid_output.ExecuteInitialize()
            self.gid_output.ExecuteBeforeSolutionLoop()
            self.gid_output.ExecuteInitializeSolutionStep()
            self.gid_output.PrintOutput()
            self.gid_output.ExecuteFinalizeSolutionStep()
            self.gid_output.ExecuteFinalize()
if __name__ == '__main__':
    km.Logger.GetDefaultOutput().SetSeverity(km.Logger.Severity.WARNING)
    KratosUnittest.main() """