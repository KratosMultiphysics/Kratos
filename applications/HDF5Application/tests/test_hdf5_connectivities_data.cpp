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
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/node.h"
#include "includes/kratos_components.h"
#include "utilities/compare_elements_and_conditions_utility.h"
#include "includes/kratos_parameters.h"

// Application includes
#include "custom_io/hdf5_file_serial.h"
#include "custom_io/hdf5_connectivities_data.h"
#include "custom_utilities/factor_elements_and_conditions_utility.h"
#include "tests/test_utils.h"

namespace Kratos
{
namespace Testing
{

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_ConnectivitiesData1, KratosHDF5TestSuite)
{
    ModelPart test_model_part;
    CreateTestModelPart(test_model_part);
    KRATOS_CHECK(test_model_part.Elements().size() > 0);
    HDF5::Internals::ConnectivitiesData data;
    data.SetData(FactorElements(test_model_part.Elements()).front());
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::WriteInfo info;
    data.WriteData(test_file, "/Elements", info);
    data.Clear();
    KRATOS_WATCH(data.size());
    KRATOS_CHECK(data.size() == 0);
    data.ReadData(test_file, "/Elements", info.StartIndex, info.BlockSize);
    HDF5::ElementsContainerType new_elements;
    data.CreateEntities(test_model_part.Nodes(), test_model_part.rProperties(), new_elements);
    CompareElements(new_elements, test_model_part.Elements());
}

KRATOS_TEST_CASE_IN_SUITE(HDF5_Internals_ConnectivitiesData2, KratosHDF5TestSuite)
{
    ModelPart test_model_part;
    CreateTestModelPart(test_model_part);
    KRATOS_CHECK(test_model_part.Conditions().size() > 0);
    HDF5::Internals::ConnectivitiesData data;
    data.SetData(FactorConditions(test_model_part.Conditions()).front());
    Parameters test_params(R"(
        {
            "file_name" : "test.h5",
            "file_driver": "core"
        })");
    HDF5::FileSerial test_file(test_params);
    HDF5::WriteInfo info;
    data.WriteData(test_file, "/Conditions", info);
    data.Clear();
    KRATOS_CHECK(data.size() == 0);
    data.ReadData(test_file, "/Conditions", info.StartIndex, info.BlockSize);
    HDF5::ConditionsContainerType new_conditions;
    data.CreateEntities(test_model_part.Nodes(), test_model_part.rProperties(), new_conditions);
    CompareConditions(new_conditions, test_model_part.Conditions());
}

} // namespace Testing
} // namespace Kratos.
