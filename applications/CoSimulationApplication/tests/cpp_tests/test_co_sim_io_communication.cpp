// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes
#include <vector>
#include <thread>

// External includes

// Project includes
#include "includes/checks.h"
#include "testing/testing.h"
#include "co_simulation_io/internals/co_sim_file_communication.h"

namespace Kratos {
namespace Testing {

namespace {

typedef CoSimIO::Internals::CoSimCommunication CoSimComm;
typedef CoSimIO::Internals::CoSimFileCommunication CoSimFileComm;
typedef CoSimIO::Internals::DataContainer<double> BaseDataContainter;
typedef CoSimIO::Internals::DataContainerStdVector<double> DataContainterStdVec;
typedef CoSimIO::Internals::DataContainerRawMemory<double> DataContainterRawMem;
typedef std::vector<std::pair<std::size_t, double>> SizeValueVectorType;

void ExportDataDetail(
    CoSimComm& rCoSimComm,
    BaseDataContainter& rDataContainer,
    const std::string& rIdentifier,
    const std::vector<std::pair<std::size_t, double>>& rSizesValues)
{
    rCoSimComm.Connect();

    for (const auto& r_size_val_pair : rSizesValues) {
        const std::size_t current_size = r_size_val_pair.first;
        rDataContainer.resize_if_smaller(current_size);
        for (std::size_t i=0; i<current_size; ++i) {
            rDataContainer[i] = r_size_val_pair.second;
        }
        rCoSimComm.ExportData(rIdentifier, current_size, rDataContainer);
    }

    rCoSimComm.Disconnect();
}

void ImportDataDetail(
    CoSimComm& rCoSimComm,
    BaseDataContainter& rDataContainer,
    const std::string& rIdentifier,
    const std::vector<std::pair<std::size_t, double>>& rSizesValues)
{
    int received_size;

    rCoSimComm.Connect();

    for (const auto& r_size_val_pair : rSizesValues) {
        rCoSimComm.ImportData(rIdentifier, received_size, rDataContainer);

        const int expected_size = r_size_val_pair.first;
        const double expected_value = r_size_val_pair.second;

        KRATOS_CHECK_EQUAL(expected_size, received_size);

        // cannot use the vector check macro, bcs the size might be larger / does not work with manual memory!
        for (int i=0; i<received_size; ++i) {
            KRATOS_CHECK_DOUBLE_EQUAL(expected_value, rDataContainer[i]);
        }
    }

    rCoSimComm.Disconnect();
}

void ExportImportData_StdVector(
    CoSimComm& rCoSimCommExport,
    CoSimComm& rCoSimCommImport,
    const SizeValueVectorType& rSizesValues)
{
    std::string identifier("dummy_data");

    // Exporting (done in separate thread to avoid deadlocks)
    std::vector<double> values_export;
    Kratos::unique_ptr<BaseDataContainter> p_data_container_export(Kratos::make_unique<DataContainterStdVec>(values_export));

    std::thread export_thread(ExportDataDetail, std::ref(rCoSimCommExport), std::ref(*p_data_container_export), identifier, std::ref(rSizesValues));

    // Importing
    std::vector<double> values_import;
    Kratos::unique_ptr<BaseDataContainter> p_data_container_import(Kratos::make_unique<DataContainterStdVec>(values_import));
    ImportDataDetail(rCoSimCommImport, *p_data_container_import, identifier, rSizesValues);

    export_thread.join();
}

// void ExportImportMesh(CoSimComm& rCoSimComm)
// {

// }

// void SendReceiveControlSignal(CoSimComm& rCoSimComm)
// {

// }

} // helpers namespace

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Data_Connect_Disconnect, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings = { // only disabling prints, otherwise use default configuration
            {"echo_level",   "0"},
            {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_Data_Connect_Disconnect");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    const SizeValueVectorType sizes_values; // assigning nothing here results in connecting and disconnecting

    ExportImportData_StdVector(*exporter, *importer, sizes_values);
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Data_Once_StdVector, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings = { // only disabling prints, otherwise use default configuration
            {"echo_level",   "0"},
            {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_Data_Once_StdVector");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    const SizeValueVectorType sizes_values = {
        {15, 1.884}
    };

    ExportImportData_StdVector(*exporter, *importer, sizes_values);
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Data_Multiple_StdVector, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings = { // only disabling prints, otherwise use default configuration
            {"echo_level",   "0"},
            {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_Data_Multiple_StdVector");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    const SizeValueVectorType sizes_values = {
        {15, 1.884},
        {25, -147.884},
        {10, 11447.556},
        {13, -55.78},
        {32, 10.123}
    };

    ExportImportData_StdVector(*exporter, *importer, sizes_values);
}



KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Data_RawMemory, KratosCoSimulationFastSuite)
{
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Mesh_StdVector, KratosCoSimulationFastSuite)
{
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Mesh_RawMemory, KratosCoSimulationFastSuite)
{
}

} // namespace Testing
}  // namespace Kratos.
