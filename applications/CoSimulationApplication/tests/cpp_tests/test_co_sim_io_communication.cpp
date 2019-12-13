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
#include "testing/testing.h"
#include "co_simulation_io/impl/co_sim_file_communication.h"

namespace Kratos {
namespace Testing {

namespace {

typedef CoSimIO::Internals::CoSimCommunication CoSimComm;
typedef CoSimIO::Internals::CoSimFileCommunication CoSimFileComm;
typedef CoSimIO::Internals::DataContainer<double> BaseDataContainterDouble;
typedef CoSimIO::Internals::DataContainerStdVector<double> DataContainterStdVecDouble;
typedef CoSimIO::Internals::DataContainerRawMemory<double> DataContainterRawMemDouble;
typedef CoSimIO::Internals::DataContainer<int> BaseDataContainterInt;
typedef CoSimIO::Internals::DataContainerStdVector<int> DataContainterStdVecInt;
typedef CoSimIO::Internals::DataContainerRawMemory<int> DataContainterRawMemInt;

typedef std::vector<std::pair<std::size_t, double>> SizeValuePairsVectorType;
typedef std::vector<std::pair<std::string, CoSimIO::ControlSignal>> IdentifierSignalPairsVectorType;

typedef std::vector<std::vector<double>> NodalCoordiantesVectorType;
typedef std::vector<std::vector<int>> ElementQuantitiesVectorType;

void SendControlSignals(
    CoSimComm& rCoSimComm,
    const IdentifierSignalPairsVectorType& rIdentifiersSignals)
{
    rCoSimComm.Connect();

    for (const auto& ident_signal_pair : rIdentifiersSignals) {
        rCoSimComm.SendControlSignal(ident_signal_pair.first, ident_signal_pair.second);
    }

    rCoSimComm.Disconnect();
}

void RecvControlSignals(
    CoSimComm& rCoSimComm,
    const IdentifierSignalPairsVectorType& rIdentifiersSignals)
{
    rCoSimComm.Connect();

    std::string received_identifier;
    CoSimIO::ControlSignal received_signal;

    for (const auto& ident_signal_pair : rIdentifiersSignals) {
        received_signal = rCoSimComm.RecvControlSignal(received_identifier);
        KRATOS_CHECK_STRING_EQUAL(ident_signal_pair.first, received_identifier);
        KRATOS_CHECK_EQUAL(static_cast<int>(ident_signal_pair.second), static_cast<int>(received_signal));
    }

    rCoSimComm.Disconnect();
}

// ****************************************************************************
// ****************************************************************************

void ExportDataDetail(
    CoSimComm& rCoSimComm,
    BaseDataContainterDouble& rDataContainer,
    const std::string& rIdentifier,
    const SizeValuePairsVectorType& rSizesValues)
{
    rCoSimComm.Connect();

    for (const auto& r_size_val_pair : rSizesValues) {
        const std::size_t current_size = r_size_val_pair.first;
        rDataContainer.resize(current_size);
        for (std::size_t i=0; i<current_size; ++i) {
            rDataContainer[i] = r_size_val_pair.second;
        }
        rCoSimComm.ExportData(rIdentifier, rDataContainer);
    }

    rCoSimComm.Disconnect();
}

void ImportDataDetail(
    CoSimComm& rCoSimComm,
    BaseDataContainterDouble& rDataContainer,
    const std::string& rIdentifier,
    const SizeValuePairsVectorType& rSizesValues)
{
    rCoSimComm.Connect();

    for (const auto& r_size_val_pair : rSizesValues) {
        rCoSimComm.ImportData(rIdentifier, rDataContainer);

        const int expected_size = r_size_val_pair.first;
        const double expected_value = r_size_val_pair.second;
        const int received_size = rDataContainer.size();
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
    const SizeValuePairsVectorType& rSizesValues)
{
    std::string identifier("dummy_data");

    // Exporting (done in separate thread to avoid deadlocks)
    std::vector<double> values_export;
    Kratos::unique_ptr<BaseDataContainterDouble> p_data_container_export(Kratos::make_unique<DataContainterStdVecDouble>(values_export));

    std::thread export_thread(ExportDataDetail, std::ref(rCoSimCommExport), std::ref(*p_data_container_export), identifier, std::ref(rSizesValues));

    // Importing
    std::vector<double> values_import;
    Kratos::unique_ptr<BaseDataContainterDouble> p_data_container_import(Kratos::make_unique<DataContainterStdVecDouble>(values_import));
    ImportDataDetail(rCoSimCommImport, *p_data_container_import, identifier, rSizesValues);

    export_thread.join();
}

void ExportImportData_RawMemory(
    CoSimComm& rCoSimCommExport,
    CoSimComm& rCoSimCommImport,
    const SizeValuePairsVectorType& rSizesValues)
{
    std::string identifier("dummy_data");

    // Exporting (done in separate thread to avoid deadlocks)
    double** values_raw_export = (double**)malloc(sizeof(double*)*1);
    values_raw_export[0]= NULL;

    Kratos::unique_ptr<BaseDataContainterDouble> p_data_container_export(Kratos::make_unique<DataContainterRawMemDouble>(values_raw_export, 0));

    std::thread export_thread(ExportDataDetail, std::ref(rCoSimCommExport), std::ref(*p_data_container_export), identifier, std::ref(rSizesValues));

    // Importing
    double** values_raw_import = (double**)malloc(sizeof(double*)*1);
    values_raw_import[0]= NULL;
    Kratos::unique_ptr<BaseDataContainterDouble> p_data_container_import(Kratos::make_unique<DataContainterRawMemDouble>(values_raw_import, 0));
    ImportDataDetail(rCoSimCommImport, *p_data_container_import, identifier, rSizesValues);

    export_thread.join();

    // deallocating memory
    free(*values_raw_export);
    free(*values_raw_import);
    free(values_raw_export);
    free(values_raw_import);
}

// ****************************************************************************
// ****************************************************************************

void ExportMeshDetail(
    CoSimComm& rCoSimComm,
    BaseDataContainterDouble& rDataContainerNodalCoords,
    BaseDataContainterInt& rDataContainerElementConnectivities,
    BaseDataContainterInt& rDataContainerElementTypes,
    const std::string& rIdentifier,
    const NodalCoordiantesVectorType& rNodeCoords,
    const ElementQuantitiesVectorType& rElemConnectivities,
    const ElementQuantitiesVectorType& rElemTypes)
{
    rCoSimComm.Connect();

    // prelim checks if input is correct
    const std::size_t num_meshes = rNodeCoords.size();
    KRATOS_CHECK_EQUAL(num_meshes, rElemConnectivities.size());
    KRATOS_CHECK_EQUAL(num_meshes, rElemTypes.size());

    for (std::size_t i=0; i<num_meshes; ++i) {
        const std::size_t size_nodal_coords = rNodeCoords[i].size();
        const std::size_t num_connectivities = rElemConnectivities[i].size();
        const std::size_t num_elems = rElemTypes[i].size();

        rDataContainerNodalCoords.resize(size_nodal_coords);
        rDataContainerElementConnectivities.resize(num_connectivities);
        rDataContainerElementTypes.resize(num_elems);

        for (std::size_t j=0; j<size_nodal_coords; ++j) {
            rDataContainerNodalCoords[j] = rNodeCoords[i][j];
        }

        for (std::size_t j=0; j<num_connectivities; ++j) {
            rDataContainerElementConnectivities[j] = rElemConnectivities[i][j];
        }

        for (std::size_t j=0; j<num_elems; ++j) {
            rDataContainerElementTypes[j] = rElemTypes[i][j];
        }

        rCoSimComm.ExportMesh(rIdentifier, rDataContainerNodalCoords, rDataContainerElementConnectivities, rDataContainerElementTypes);
    }

    rCoSimComm.Disconnect();
}

void ImportMeshDetail(
    CoSimComm& rCoSimComm,
    BaseDataContainterDouble& rDataContainerNodalCoords,
    BaseDataContainterInt& rDataContainerElementConnectivities,
    BaseDataContainterInt& rDataContainerElementTypes,
    const std::string& rIdentifier,
    const NodalCoordiantesVectorType& rNodeCoords,
    const ElementQuantitiesVectorType& rElemConnectivities,
    const ElementQuantitiesVectorType& rElemTypes)
{

    rCoSimComm.Connect();

    // prelim checks if input is correct
    const std::size_t num_meshes = rNodeCoords.size();
    KRATOS_CHECK_EQUAL(num_meshes, rElemConnectivities.size());
    KRATOS_CHECK_EQUAL(num_meshes, rElemTypes.size());

    for (std::size_t i=0; i<num_meshes; ++i) {
        rCoSimComm.ImportMesh(rIdentifier, rDataContainerNodalCoords, rDataContainerElementConnectivities, rDataContainerElementTypes);

        const int received_num_nodes = rDataContainerNodalCoords.size()/3;
        const int received_num_elems = rDataContainerElementTypes.size();
        const int exp_num_nodes = rNodeCoords[i].size()/3;
        const int exp_num_elems = rElemTypes[i].size();
        KRATOS_CHECK_EQUAL(exp_num_nodes, received_num_nodes);
        KRATOS_CHECK_EQUAL(exp_num_elems, received_num_elems);

        // cannot use the vector check macro, bcs the size might be larger / does not work with manual memory!
        for (std::size_t j=0; j<rNodeCoords[i].size(); ++j) {
            KRATOS_CHECK_DOUBLE_EQUAL(rNodeCoords[i][j], rDataContainerNodalCoords[j]);
        }

        for (int j=0; j<exp_num_elems; ++j) {
            KRATOS_CHECK_EQUAL(rElemTypes[i][j], rDataContainerElementTypes[j]);
        }

        for (std::size_t j=0; j<rElemConnectivities[i].size(); ++j) {
            KRATOS_CHECK_EQUAL(rElemConnectivities[i][j], rDataContainerElementConnectivities[j]);
        }
    }

    rCoSimComm.Disconnect();
}

void ExportImportMesh_StdVector(
    CoSimComm& rCoSimCommExport,
    CoSimComm& rCoSimCommImport,
    const NodalCoordiantesVectorType& rNodeCoords,
    const ElementQuantitiesVectorType& rElemConnectivities,
    const ElementQuantitiesVectorType& rElemTypes)
{
    std::string identifier("dummy_mesh");

    // Exporting (done in separate thread to avoid deadlocks)
    std::vector<double> coords_vec_exp;
    std::vector<int> conn_vec_exp;
    std::vector<int> types_vec_exp;
    Kratos::unique_ptr<BaseDataContainterDouble> p_data_container_export_coords(Kratos::make_unique<DataContainterStdVecDouble>(coords_vec_exp));
    Kratos::unique_ptr<BaseDataContainterInt> p_data_container_export_conn(Kratos::make_unique<DataContainterStdVecInt>(conn_vec_exp));
    Kratos::unique_ptr<BaseDataContainterInt> p_data_container_export_types(Kratos::make_unique<DataContainterStdVecInt>(types_vec_exp));

    std::thread export_thread(ExportMeshDetail,
        std::ref(rCoSimCommExport),
        std::ref(*p_data_container_export_coords),
        std::ref(*p_data_container_export_conn),
        std::ref(*p_data_container_export_types),
        identifier,
        std::ref(rNodeCoords),
        std::ref(rElemConnectivities),
        std::ref(rElemTypes));

    // Importing
    std::vector<double> coords_vec_imp;
    std::vector<int> conn_vec_imp;
    std::vector<int> types_vec_imp;
    Kratos::unique_ptr<BaseDataContainterDouble> p_data_container_import_coords(Kratos::make_unique<DataContainterStdVecDouble>(coords_vec_imp));
    Kratos::unique_ptr<BaseDataContainterInt> p_data_container_import_conn(Kratos::make_unique<DataContainterStdVecInt>(conn_vec_imp));
    Kratos::unique_ptr<BaseDataContainterInt> p_data_container_import_types(Kratos::make_unique<DataContainterStdVecInt>(types_vec_imp));
    ImportMeshDetail(
        rCoSimCommImport,
        *p_data_container_import_coords,
        *p_data_container_import_conn,
        *p_data_container_import_types,
        identifier,
        rNodeCoords,
        rElemConnectivities,
        rElemTypes);

    export_thread.join();
}

void ExportImportMesh_RawMemory(
    CoSimComm& rCoSimCommExport,
    CoSimComm& rCoSimCommImport,
    const NodalCoordiantesVectorType& rNodeCoords,
    const ElementQuantitiesVectorType& rElemConnectivities,
    const ElementQuantitiesVectorType& rElemTypes)
{
    std::string identifier("dummy_mesh");

    // Exporting (done in separate thread to avoid deadlocks)
    double** coords_raw_exp = (double**)malloc(sizeof(double*)*1);
    coords_raw_exp[0]= NULL;
    int** conn_raw_exp = (int**)malloc(sizeof(int*)*1);
    conn_raw_exp[0]= NULL;
    int** types_raw_exp = (int**)malloc(sizeof(int*)*1);
    types_raw_exp[0]= NULL;

    Kratos::unique_ptr<BaseDataContainterDouble> p_data_container_export_coords(Kratos::make_unique<DataContainterRawMemDouble>(coords_raw_exp, 0));
    Kratos::unique_ptr<BaseDataContainterInt> p_data_container_export_conn(Kratos::make_unique<DataContainterRawMemInt>(conn_raw_exp, 0));
    Kratos::unique_ptr<BaseDataContainterInt> p_data_container_export_types(Kratos::make_unique<DataContainterRawMemInt>(types_raw_exp, 0));

    std::thread export_thread(ExportMeshDetail,
        std::ref(rCoSimCommExport),
        std::ref(*p_data_container_export_coords),
        std::ref(*p_data_container_export_conn),
        std::ref(*p_data_container_export_types),
        identifier,
        std::ref(rNodeCoords),
        std::ref(rElemConnectivities),
        std::ref(rElemTypes));

    // Importing
    double** coords_raw_imp = (double**)malloc(sizeof(double*)*1);
    coords_raw_imp[0]= NULL;
    int** conn_raw_imp = (int**)malloc(sizeof(int*)*1);
    conn_raw_imp[0]= NULL;
    int** types_raw_imp = (int**)malloc(sizeof(int*)*1);
    types_raw_imp[0]= NULL;

    Kratos::unique_ptr<BaseDataContainterDouble> p_data_container_import_coords(Kratos::make_unique<DataContainterRawMemDouble>(coords_raw_imp, 0));
    Kratos::unique_ptr<BaseDataContainterInt> p_data_container_import_conn(Kratos::make_unique<DataContainterRawMemInt>(conn_raw_imp, 0));
    Kratos::unique_ptr<BaseDataContainterInt> p_data_container_import_types(Kratos::make_unique<DataContainterRawMemInt>(types_raw_imp, 0));

    ImportMeshDetail(
        rCoSimCommImport,
        *p_data_container_import_coords,
        *p_data_container_import_conn,
        *p_data_container_import_types,
        identifier,
        rNodeCoords,
        rElemConnectivities,
        rElemTypes);

    export_thread.join();

    // deallocating memory
    free(*coords_raw_exp);
    free(*conn_raw_exp);
    free(*types_raw_exp);
    free(*coords_raw_imp);
    free(*conn_raw_imp);
    free(*types_raw_imp);
    free(coords_raw_exp);
    free(conn_raw_exp);
    free(types_raw_exp);
    free(coords_raw_imp);
    free(conn_raw_imp);
    free(types_raw_imp);
}

} // helpers namespace

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Data_Connect_Disconnect, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings { // only disabling prints, otherwise use default configuration
        {"echo_level",   "0"},
        {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_Data_Connect_Disconnect");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    const SizeValuePairsVectorType sizes_values; // assigning nothing here results in connecting and disconnecting

    ExportImportData_StdVector(*exporter, *importer, sizes_values);
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_ControlSignal_Once, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings { // only disabling prints, otherwise use default configuration
        {"echo_level",   "0"},
        {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_ControlSignal_Once");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    const IdentifierSignalPairsVectorType identifier_signals_pairs {
        {"ab258cc", CoSimIO::ControlSignal::Dummy}
    };

    std::thread export_thread(SendControlSignals, std::ref(*exporter), std::ref(identifier_signals_pairs));

    RecvControlSignals(*importer, identifier_signals_pairs);

    export_thread.join();
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_ControlSignal_Multiple, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings { // only disabling prints, otherwise use default configuration
        {"echo_level",   "0"},
        {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_ControlSignal_Multiple");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    const IdentifierSignalPairsVectorType identifier_signals_pairs {
        {"abcc", CoSimIO::ControlSignal::Dummy},
        {"frt1", CoSimIO::ControlSignal::ConvergenceAchieved},
        {"bbhh", CoSimIO::ControlSignal::AdvanceInTime},
        {"oott", CoSimIO::ControlSignal::SolveSolutionStep},
        {"lkuz", CoSimIO::ControlSignal::ExportData}
    };

    std::thread export_thread(SendControlSignals, std::ref(*exporter), std::ref(identifier_signals_pairs));

    RecvControlSignals(*importer, identifier_signals_pairs);

    export_thread.join();
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Data_Once, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings { // only disabling prints, otherwise use default configuration
        {"echo_level",   "0"},
        {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_Data_Once");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    const SizeValuePairsVectorType sizes_values {
        {15, 1.884}
    };

    ExportImportData_StdVector(*exporter, *importer, sizes_values);
    ExportImportData_RawMemory(*exporter, *importer, sizes_values);
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Data_Multiple, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings { // only disabling prints, otherwise use default configuration
        {"echo_level",   "0"},
        {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_Data_Multiple");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    const SizeValuePairsVectorType sizes_values {
        {15, 1.884},
        {25, -147.884},
        {10, 11447.556},
        {13, -55.78},
        {32, 10.123}
    };

    ExportImportData_StdVector(*exporter, *importer, sizes_values);
    ExportImportData_RawMemory(*exporter, *importer, sizes_values);
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Mesh_Once, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings { // only disabling prints, otherwise use default configuration
        {"echo_level",   "0"},
        {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_Mesh_Once");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    /*         6
       7 x-----x-----x 8
         |     |     |
         |  1  |  2  |
         |     |     |
       3 x-----x-----x 5
         |     |2    |
         |  3  |  4  |
         |     |     |
       0 x-----x-----x 4
               1
    */

    NodalCoordiantesVectorType node_coords {{
        -1.5, -1.8, 0.0, // 0
         0.0, -1.8, 0.0, // 1
         0.0,  0.0, 0.0, // 2
        -1.5,  0.0, 0.0, // 3
         1.5, -1.8, 0.0, // 4
         1.5,  0.0, 0.0, // 5
         0.0,  1.8, 0.0, // 6
        -1.5,  1.8, 0.0, // 7
         1.5,  1.8, 0.0  // 8
    }};

    ElementQuantitiesVectorType elem_connectivities {{
        7, 3, 2, 6, // 1
        6, 2, 5, 8, // 2
        3, 0, 1, 2, // 3
        2, 1, 4, 5  // 4
    }};

    ElementQuantitiesVectorType elem_types {{9,9,9,9}}; // VTK_QUAD

    std::vector<int> connectivities_offsets {1}; // required otherwise import fails

    ExportImportMesh_StdVector(*exporter, *importer, node_coords, elem_connectivities, elem_types);
    ExportImportMesh_RawMemory(*exporter, *importer, node_coords, elem_connectivities, elem_types);
}

KRATOS_TEST_CASE_IN_SUITE(FileCommunication_Mesh_Multiple, KratosCoSimulationFastSuite)
{
    CoSimIO::SettingsType settings { // only disabling prints, otherwise use default configuration
        {"echo_level",   "0"},
        {"print_timing", "0"}
    };
    std::string connection_name("FileCommunication_Mesh_Multiple");

    Kratos::unique_ptr<CoSimComm> exporter(Kratos::make_unique<CoSimFileComm>(connection_name, settings, true));
    Kratos::unique_ptr<CoSimComm> importer(Kratos::make_unique<CoSimFileComm>(connection_name, settings, false));

    /*    -- Mesh 1 --
               6
       7 x-----x-----x 8
         |     |     |
         |  1  |  2  |
         |     |     |
       3 x-----x-----x 5
         |     |2    |
         |  3  |  4  |
         |     |     |
       0 x-----x-----x 4
               1
    */

    /*    -- Mesh 2 --
        0      2      3
        x------x------x
         \     |     /|\
          \  1 |  2 / | \
           \   |   /  |  \
            \  |  /   |   \
             \ | /  3 |  4 \
              \|/     |     \
               x------x-----x
               1      4     5
    */

    NodalCoordiantesVectorType node_coords {{
        -1.5, -1.8, 0.0, // 0
         0.0, -1.8, 0.0, // 1
         0.0,  0.0, 0.0, // 2
        -1.5,  0.0, 0.0, // 3
         1.5, -1.8, 0.0, // 4
         1.5,  0.0, 0.0, // 5
         0.0,  1.8, 0.0, // 6
        -1.5,  1.8, 0.0, // 7
         1.5,  1.8, 0.0  // 8
    },{
        0.0, 2.5, 1.0, // 0
        2.0, 0.0, 1.5, // 1
        2.0, 2.5, 1.5, // 2
        4.0, 2.5, 1.7, // 3
        4.0, 0.0, 1.7, // 4
        6.0, 0.0, 1.8 //  5
    },{
        -1.5, -1.8, 0.0, // 0
         0.0, -1.8, 0.0, // 1
         0.0,  0.0, 0.0, // 2
        -1.5,  0.0, 0.0, // 3
    }};

    ElementQuantitiesVectorType elem_connectivities {{
        7, 3, 2, 6, // 1
        6, 2, 5, 8, // 2
        3, 0, 1, 2, // 3
        2, 1, 4, 5  // 4
    },{
        0, 1, 2, // 1
        1, 3, 2, // 2
        1, 4, 3, // 3
        3, 4, 5, // 4
    },{ }}; // no elements!

    ElementQuantitiesVectorType elem_types {
        {9,9,9,9}, // VTK_QUAD
        {5,5,5,5}, // VTK_TRIANGLE
        {} // no elements!
    };

    ExportImportMesh_StdVector(*exporter, *importer, node_coords, elem_connectivities, elem_types);
    ExportImportMesh_RawMemory(*exporter, *importer, node_coords, elem_connectivities, elem_types);
}

} // namespace Testing
}  // namespace Kratos.
