
import KratosMultiphysics as Kratos
from KratosMultiphysics.MeshingApplication.mmg_process import MmgProcess

def RefineMesh(
    input_mdpa_name,
    output_mdpa_name,
    mmg_parameters,
    solution_step_variables_list,
    set_variable_data,
    set_model_parts = lambda _: _):

    model = Kratos.Model()
    model_part = model.CreateModelPart("MeshingModelPart")

    for solution_step_variable in solution_step_variables_list:
        model_part.AddNodalSolutionStepVariable(solution_step_variable)

    Kratos.ModelPartIO(input_mdpa_name, Kratos.IO.READ).ReadModelPart(model_part)

    set_model_parts(model_part)
    set_variable_data(model_part)

    model_part.Set(Kratos.MODIFIED, False)
    model_part.ProcessInfo.SetValue(Kratos.DOMAIN_SIZE, 2)
    model_part.ProcessInfo.SetValue(Kratos.STEP, 1)
    model_part.ProcessInfo.SetValue(Kratos.TIME, 1.0)

    mmg_parameters["model_part_name"].SetString("MeshingModelPart")
    remeshing_process = MmgProcess(model, mmg_parameters)
    remeshing_process.ExecuteInitialize()
    remeshing_process.ExecuteInitializeSolutionStep()
    remeshing_process.ExecuteFinalizeSolutionStep()
    remeshing_process.ExecuteFinalize()

    Kratos.ModelPartIO(output_mdpa_name, Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(model_part)
    return model_part




