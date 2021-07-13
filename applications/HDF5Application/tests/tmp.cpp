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
    ModelPart& r_model_part = MakeModelPart(model);

    try
    {
        Parameters file_parameters(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "sec2"
        })");
        
        auto p_file = HDF5::File::Pointer(new HDF5::FileSerial(file_parameters));

        Parameters write_parameters(R"({
            "prefix" : "/prefix",
            "list_of_variables" : ["DISPLACEMENT"]
        })");

        HDF5::Detail::VertexContainerType vertices;
        BruteForcePointLocator locator(r_model_part);

        vertices.push_back(HDF5::Detail::Vertex::Pointer(new HDF5::Detail::Vertex(
            array_1d<double,3>{{0.0, 0.0, 0.0}},
            r_model_part)));
        vertices.push_back(HDF5::Detail::Vertex::Pointer(new HDF5::Detail::Vertex(
            array_1d<double,3>{{0.5, 0.5, 0.0}},
            r_model_part)));
        vertices.push_back(HDF5::Detail::Vertex::Pointer(new HDF5::Detail::Vertex(
            array_1d<double,3>{{1.0, 1.0, 0.0}},
            r_model_part)));

        for (HDF5::Detail::Vertex& rVertex : vertices) {
            rVertex.Locate(locator);
        }

        try
        {
            HDF5::VertexContainerIO writer(write_parameters, p_file);
            writer.WriteCoordinates(vertices);
            writer.WriteVariables(vertices);
        }
        catch (std::exception& rException)
        {
            puts(rException.what());
            puts("nope");
            return;
        }
    }
    catch (std::exception& rException)
    {
        puts(rException.what());
        puts("something went wrong");
        return;
    }
}


} // namespace Testing
} // namespace Kratos