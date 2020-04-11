//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Suneth Warnakulasuriya, https://github.com/sunethwarna
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "testing/testing.h"

// Application includes
#include "custom_io/hdf5_element_data_value_io.h"
#include "custom_io/hdf5_file_serial.h"
#include "tests/test_utils.h"

namespace Kratos
{
namespace Testing
{
KRATOS_TEST_CASE_IN_SUITE(HDF5PointsData_ReadElementResults, KratosHDF5TestSuite)
{
    Parameters file_params(R"(
        {
            "file_name" : "test.h5",
            "file_access_mode": "exclusive",
            "file_driver": "core"
        })");
    auto p_test_file = Kratos::make_shared<HDF5::FileSerial>(file_params);

    Model this_model;
    ModelPart& r_read_model_part = this_model.CreateModelPart("test_read");
    ModelPart& r_write_model_part = this_model.CreateModelPart("test_write");
    TestModelPartFactory::CreateModelPart(r_write_model_part,
                                          {{"Element2D3N"}});
    TestModelPartFactory::CreateModelPart(r_read_model_part, {{"Element2D3N"}});

    r_read_model_part.SetBufferSize(2);
    r_write_model_part.SetBufferSize(2);

    std::vector<std::string> variables_list = {
        {"DISPLACEMENT"}, {"PRESSURE"}, {"REFINEMENT_LEVEL"}, {"GREEN_LAGRANGE_STRAIN_TENSOR"}};

    for (auto& r_element : r_write_model_part.Elements())
    {
        TestModelPartFactory::AssignDataValueContainer(
            r_element.Data(), r_element, variables_list);
    }

    Parameters io_params(R"(
        {
            "prefix": "/Step",
            "list_of_variables": ["DISPLACEMENT", "PRESSURE", "REFINEMENT_LEVEL", "GREEN_LAGRANGE_STRAIN_TENSOR"]
        })");

    HDF5::ElementDataValueIO data_io(io_params, p_test_file);
    data_io.WriteElementResults(r_write_model_part.Elements());
    data_io.ReadElementResults(r_read_model_part.Elements(),
                               r_read_model_part.GetCommunicator());

    for (auto& r_write_element : r_write_model_part.Elements())
    {
        HDF5::ElementType& r_read_element =
            r_read_model_part.Elements()[r_write_element.Id()];
        CompareDataValueContainers(r_read_element.Data(), r_read_element,
                                   r_write_element.Data(), r_write_element);
    }
}

} // namespace Testing
} // namespace Kratos.
