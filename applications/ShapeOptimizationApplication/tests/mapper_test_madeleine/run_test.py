import KratosMultiphysics as km
import KratosMultiphysics.ShapeOptimizationApplication as kso
import mapper_factory
from vtk_output_process import VtkOutputProcess
import distance_response

def import_model_part(model, file_name):
    model_part = model.CreateModelPart(file_name, 1)

    # Adding variables BEFORE reading the .mdpa
    model_part.AddNodalSolutionStepVariable(kso.DF1DX)
    model_part.AddNodalSolutionStepVariable(kso.DF1DX_MAPPED)
    model_part.AddNodalSolutionStepVariable(kso.DC1DX)
    model_part.AddNodalSolutionStepVariable(kso.DC1DX_MAPPED)
    model_part.AddNodalSolutionStepVariable(kso.SHAPE_CHANGE)
    model_part.AddNodalSolutionStepVariable(kso.SHAPE_UPDATE)

    model_part_io = km.ModelPartIO(file_name) #path to file without ".mdpa"
    model_part_io.ReadModelPart(model_part)
    return model_part

def create_io(model, model_part_name, folder_name="VTK_Output"):
    vtk_output_configuration = km.Parameters("""{
            "model_part_name"        : \""""+model_part_name+"""\",
            "folder_name"            : \""""+folder_name+"""\",
            "custom_name_prefix"     : \""""+folder_name+"""\",
            "output_sub_model_parts" : false,
            "nodal_solution_step_data_variables" : ["DF1DX", "DF1DX_MAPPED", "DC1DX", "DC1DX_MAPPED", "SHAPE_CHANGE", "SHAPE_UPDATE"]
        }""")

    vtk_output = VtkOutputProcess(model, vtk_output_configuration)
    return vtk_output

def write_output_step(vtk_output):
    vtk_output.ExecuteInitialize()
    vtk_output.ExecuteBeforeSolutionLoop()
    vtk_output.ExecuteInitializeSolutionStep()
    vtk_output.PrintOutput()
    vtk_output.ExecuteFinalizeSolutionStep()
    vtk_output.ExecuteFinalize()

def write_output(model, model_part_name, folder_name="VTK_Output"):
    vtk_output = create_io(model, model_part_name, folder_name)
    write_output_step(vtk_output)


def run(mp_name, mapper_settings, folder_name):
    model = km.Model()
    model_part = import_model_part(model, mp_name)

    #distance_response.DiscreteVolumeResponse().GetGradient(model_part, kso.DF1DX)
    #distance_response.SumOfNodalDistancesResponse().GetGradient(model_part, kso.DF1DX)
    #distance_response.SumOfNodalDistancesLinearResponse().GetGradient(model_part, kso.DF1DX)
    distance_response.DiscreteLinearResponse().GetGradient(model_part, kso.DF1DX)

    mapper = mapper_factory.CreateMapper(model_part, model_part, mapper_settings)
    mapper.InverseMap(kso.DF1DX, kso.DF1DX_MAPPED)
    mapper.Map(kso.DF1DX_MAPPED, kso.SHAPE_UPDATE)

    write_output(model, mp_name, folder_name=folder_name)

def main():

    conservative_node_sum = km.Parameters("""
    {
        "filter_function_type"       : "linear",
        "filter_radius"              : 0.4,
        "max_nodes_in_filter_radius" : 10000
    }""")

    #run("split_plate", conservative_node_sum, "s_conse_ns")
    run("diagonal_plate", conservative_node_sum, "d_conse_ns")

    consistent_node_sum = km.Parameters("""
    {
        "filter_function_type"       : "linear",
        "filter_radius"              : 0.4,
        "max_nodes_in_filter_radius" : 10000,
        "consistent_mapping"         : true
    }""")

    #run("split_plate", consistent_node_sum, "s_consi_ns")
    run("diagonal_plate", consistent_node_sum, "d_consi_ns")


    # Test improved integration
    conservative_GI3 = km.Parameters("""
    {
        "filter_function_type"       : "linear",
        "filter_radius"              : 0.4,
        "max_nodes_in_filter_radius" : 10000,
        "improved_integration"       : true,
        "integration_method"         : "gauss_integration",
        "number_of_gauss_points"     : 3
    }""")

    #run("split_plate", conservative_GI3, "s_conse_GI3")
    run("diagonal_plate", conservative_GI3, "d_conse_GI3")

    # Test improved integration
    consistent_GI3 = km.Parameters("""
    {
        "filter_function_type"       : "linear",
        "filter_radius"              : 0.4,
        "max_nodes_in_filter_radius" : 10000,
        "consistent_mapping"         : true,
        "improved_integration"       : true,
        "integration_method"         : "gauss_integration",
        "number_of_gauss_points"     : 3
    }""")

    #run("split_plate", consistent_GI3, "s_consi_GI3")
    run("diagonal_plate", consistent_GI3, "d_consi_GI3")


    # Test improved integration
    consistent_area = km.Parameters("""
    {
        "filter_function_type"       : "linear",
        "filter_radius"              : 0.4,
        "max_nodes_in_filter_radius" : 10000,
        "improved_integration"       : true,
        "integration_method"         : "area_weighted_sum"
    }""")

    #run("split_plate", consistent_area, "s_conse_a")
    run("diagonal_plate", consistent_area, "d_conse_a")

    # Test improved integration
    consistent_area = km.Parameters("""
    {
        "filter_function_type"       : "linear",
        "filter_radius"              : 0.4,
        "max_nodes_in_filter_radius" : 10000,
        "consistent_mapping"         : true,
        "improved_integration"       : true,
        "integration_method"         : "area_weighted_sum"
    }""")

    #run("split_plate", consistent_area, "s_consi_a")
    run("diagonal_plate", consistent_area, "d_consi_a")

if __name__=="__main__":
    main()