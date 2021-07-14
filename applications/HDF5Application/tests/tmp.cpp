// --- Core Includes ---
#include "testing/testing.h"
#include "containers/model.h"
#include "utilities/brute_force_point_locator.h"


// --- Application Includes ---
#include "custom_io/hdf5_container_component_io.cpp"
#include "custom_io/hdf5_container_component_io.h"
#include "custom_io/hdf5_file.h"
#include "custom_io/hdf5_vertex_container_io.h"
#include "hdf5_application_define.h"
#include "tests/test_utils.h"
#include "custom_utilities/vertex.h"
#include "custom_processes/hdf5_point_set_output_process.h"


namespace Kratos
{
namespace Testing
{


ModelPart& MakeModelPart(Model& rModel)
{
    ModelPart& r_model_part = rModel.CreateModelPart("main");
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);

    r_model_part.CreateNewNode(
        1,
        0.0, 0.0, 0.0)->GetSolutionStepValue(DISPLACEMENT) = array_1d<double,3>{{0.0, 0.0, 0.0}};
    r_model_part.CreateNewNode(
        2,
        1.0, 0.0, 0.0)->GetSolutionStepValue(DISPLACEMENT) = array_1d<double,3>{{1.0, 0.0, 0.0}};
    r_model_part.CreateNewNode(
        3,
        1.0, 1.0, 0.0)->GetSolutionStepValue(DISPLACEMENT) = array_1d<double,3>{{2.0, 0.0, 0.0}};
    r_model_part.CreateNewNode(
        4,
        0.0, 1.0, 0.0)->GetSolutionStepValue(DISPLACEMENT) = array_1d<double,3>{{3.0, 0.0, 0.0}};

    r_model_part.CreateNewElement(
        "Element2D4N",
        0,
        {1, 2, 3, 4},
        make_shared<Properties>(0));

    return r_model_part;
}


KRATOS_TEST_CASE_IN_SUITE(HDF5_TMP, KratosHDF5TestSuite)
{
    Model model;
    MakeModelPart(model);

    Parameters parameters(R"(
    {
        "model_part_name" : "main",
        "positions" : [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [1.0, 1.0, 0.0]
        ],
        "output_variables" : ["DISPLACEMENT"],
        "file_path" : "test.h5"
    })");

    HDF5::PointSetOutputProcess point_output_process(parameters, model);
    point_output_process.ExecuteInitialize();
    point_output_process.ExecuteFinalizeSolutionStep();
}


} // namespace Testing
} // namespace Kratos