// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef KRATOS_CO_SIM_FILE_COMM_H_INCLUDED
#define KRATOS_CO_SIM_FILE_COMM_H_INCLUDED

// System includes
#include <chrono>
#include <thread>
#include <iomanip>
#include <algorithm>

// std::filesystem is needed for file communication in a folder
// std::filesystem is part of C++17 and not supported by every compiler. Here we check if it is available.
#if defined(__cplusplus) && __cplusplus >= 201703L
    #if defined(__has_include) && __has_include(<filesystem>) // has_include is C++17, hence has to be checked in a separate line
        #define KRATOS_CO_SIM_IO_FILESYSTEM_AVAILABLE
        #include <filesystem>
        namespace fs = std::filesystem;
    #elif __has_include(<experimental/filesystem>)
        #define KRATOS_CO_SIM_IO_FILESYSTEM_AVAILABLE
        #include <experimental/filesystem>
        namespace fs = std::experimental::filesystem;
    #endif
#endif


// Project includes
#include "co_sim_communication.h"

namespace CoSimIO {
namespace Internals {

namespace { // helpers namespace

static double ElapsedSeconds(const std::chrono::steady_clock::time_point& rStartTime)
{
    using namespace std::chrono;
    return duration_cast<duration<double>>(steady_clock::now() - rStartTime).count();
}

static bool FileExists(const std::string& rFileName)
{
    std::ifstream infile(rFileName);
    return infile.good(); // no need to close manually
}

static void RemoveFile(const std::string& rFileName)
{
    if (std::remove(rFileName.c_str()) != 0) {
        KRATOS_CO_SIM_INFO("CoSimIO") << "Warning: \"" << rFileName << "\" could not be deleted!" << std::endl;
    }
}

template <typename T>
static void CheckStream(const T& rStream, const std::string& rFileName)
{
    KRATOS_CO_SIM_ERROR_IF_NOT(rStream.is_open()) << rFileName << " could not be opened!" << std::endl;
}

static int GetNumNodesForVtkCellType(const int VtkCellType)
{
    const std::unordered_map<int, int> vtk_cell_type_map {
        { /*Point3D,          */ 1 ,  1},
        { /*Line3D2,          */ 3 ,  2},
        { /*Triangle3D3,      */ 5 ,  3},
        { /*Quadrilateral3D4, */ 9 ,  4},
        { /*Tetrahedra3D4,    */ 10 , 4},
        { /*Hexahedra3D8,     */ 12 , 8},
        { /*Prism3D6,         */ 13 , 6},
        { /*Line3D3,          */ 21 , 3},
        { /*Triangle3D6,      */ 22 , 6},
        { /*Quadrilateral3D8, */ 23 , 7},
        { /*Tetrahedra3D10,   */ 24,  10}
    };

    if (vtk_cell_type_map.count(VtkCellType) > 0) {
        return vtk_cell_type_map.at(VtkCellType);
    } else {
        KRATOS_CO_SIM_ERROR << "Unsupported cell type: " << VtkCellType << std::endl;
        return 0;
    }
}

} // helpers namespace


class CoSimFileCommunication : public CoSimCommunication
{
public:
    explicit CoSimFileCommunication(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster)
        : CoSimCommunication(rName, rSettings, IsConnectionMaster)
    {
        const SettingsType default_settings {
            {"use_folder_for_communication" , "0"}
        };
        Internals::AddMissingSettings(default_settings, mrSettings);

        mCommFolder = ".CoSimIOFileComm_"+rName;
        mCommInFolder = (mrSettings.at("use_folder_for_communication") == "1");

        #ifndef KRATOS_CO_SIM_IO_FILESYSTEM_AVAILABLE
        KRATOS_CO_SIM_ERROR_IF(mCommInFolder) << "Communication is a folder can only be used if std::filesystem (C++17) is available" << std::endl;
        #endif

        if (mCommInFolder && GetIsConnectionMaster()) {
            #ifdef KRATOS_CO_SIM_IO_FILESYSTEM_AVAILABLE
            // delete and recreate directory to remove potential leftovers
            fs::remove_all(mCommFolder);
            fs::create_directory(mCommFolder);
            #endif
        }
    }

    ~CoSimFileCommunication() override
    {
        if (GetIsConnected()) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
            Disconnect();
        }
    }

private:

    std::string mCommFolder = "";
    bool mCommInFolder = false;

    bool ConnectDetail() override
    {
        return true; // nothing needed here for file-based communication (maybe do sth here?)
    }

    bool DisconnectDetail() override
    {
        return true; // nothing needed here for file-based communication (maybe do sth here?)
    }

     void ImportDataImpl(
        const std::string& rIdentifier,
        CoSimIO::Internals::DataContainer<double>& rData) override
    {
        const std::string file_name(GetFullPath("CoSimIO_data_" + GetConnectionName() + "_" + rIdentifier + ".dat"));

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive array \"" << rIdentifier << "\" in file \"" << file_name << "\" ..." << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << " ImportDataImpl " << std::endl;
        WaitForFile(file_name);

        const auto start_time(std::chrono::steady_clock::now());

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        input_file >> std::setprecision(14); // TODO maybe this should be configurable

        int size_read;
        input_file >> size_read; // the first number in the file is the size of the array

        rData.resize(size_read);

        for (int i=0; i<size_read; ++i) {
            input_file >> rData[i];
        }

        RemoveFile(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving array with size: " << size_read << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetPrintTiming()) << "Receiving Array \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ExportDataImpl(
        const std::string& rIdentifier,
        const CoSimIO::Internals::DataContainer<double>& rData) override
    {
        const std::string file_name(GetFullPath("CoSimIO_data_" + GetConnectionName() + "_" + rIdentifier + ".dat"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        const int size = rData.size();
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to send array \"" << rIdentifier << "\" with size: " << size << " in file \"" << file_name << "\" ..." << std::endl;

        const auto start_time(std::chrono::steady_clock::now());

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << std::scientific << std::setprecision(14); // TODO maybe this should be configurable

        output_file << size << "\n";

        for (int i=0; i<size-1; ++i) {
            output_file << rData[i] << " ";
        }
        // TODO check if size == 0!
        output_file << rData[size-1]; // outside to not have trailing whitespace

        output_file.close();
        MakeFileVisible(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished sending array" << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetPrintTiming()) << "Sending Array \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ImportMeshImpl(
        const std::string& rIdentifier,
        CoSimIO::Internals::DataContainer<double>& rNodalCoordinates,
        CoSimIO::Internals::DataContainer<int>& rElementConnectivities,
        CoSimIO::Internals::DataContainer<int>& rElementTypes) override
    {
        const std::string file_name(GetFullPath("CoSimIO_mesh_" + GetConnectionName() + "_" + rIdentifier + ".vtk"));

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive mesh \"" << rIdentifier << "\" in file \"" << file_name << "\" ..." << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << " ImportMeshImpl " << std::endl;
        WaitForFile(file_name);

        const auto start_time(std::chrono::steady_clock::now());

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        // reading file
        std::string current_line;
        bool nodes_read = false;
        bool cells_read = false;

        while (std::getline(input_file, current_line)) {
            // reading nodes
            if (current_line.find("POINTS") != std::string::npos) {
                KRATOS_CO_SIM_ERROR_IF(nodes_read) << "The nodes were read already!" << std::endl;
                KRATOS_CO_SIM_ERROR_IF(cells_read) << "The cells were read already!" << std::endl;
                nodes_read = true;

                int num_nodes;
                current_line = current_line.substr(current_line.find("POINTS") + 7); // removing "POINTS"
                std::istringstream line_stream(current_line);
                line_stream >> num_nodes;

                KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Mesh contains " << num_nodes << " Nodes" << std::endl;

                rNodalCoordinates.resize(3*num_nodes);

                for (int i=0; i<num_nodes*3; ++i) {
                    input_file >> rNodalCoordinates[i];
                }
            }

            // reading cells
            if (current_line.find("CELLS") != std::string::npos) {
                KRATOS_CO_SIM_ERROR_IF_NOT(nodes_read) << "The nodes were not yet read!" << std::endl;
                KRATOS_CO_SIM_ERROR_IF(cells_read) << "The cells were read already!" << std::endl;
                cells_read = true;

                int num_nodes_per_cell, num_cells, elem_conn, cell_list_size;
                current_line = current_line.substr(current_line.find("CELLS") + 6); // removing "CELLS"
                std::istringstream line_stream(current_line);
                line_stream >> num_cells;
                line_stream >> cell_list_size;

                rElementConnectivities.resize(cell_list_size);
                rElementTypes.resize(num_cells);

                KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Mesh contains " << num_cells << " Elements" << std::endl;

                int counter=0;
                for (int i=0; i<num_cells; ++i) {
                    input_file >> num_nodes_per_cell;
                    for (int j=0; j<num_nodes_per_cell; ++j) {
                        input_file >> elem_conn;
                        rElementConnectivities[counter++] = elem_conn;
                    }
                }
            }

            // reading cell types
            if (current_line.find("CELL_TYPES") != std::string::npos) {
                KRATOS_CO_SIM_ERROR_IF_NOT(nodes_read) << "The nodes were not yet read!" << std::endl;
                KRATOS_CO_SIM_ERROR_IF_NOT(cells_read) << "The cells were not yet read!" << std::endl;

                for (std::size_t i=0; i<rElementTypes.size(); ++i) { // rElementTypes was resized to correct size above
                    input_file >> rElementTypes[i];
                }

            }
        }

        RemoveFile(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving mesh" << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetPrintTiming()) << "Receiving Mesh \"" << file_name << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ExportMeshImpl(
        const std::string& rIdentifier,
        CoSimIO::Internals::DataContainer<double>& rNodalCoordinates,
        CoSimIO::Internals::DataContainer<int>& rElementConnectivities,
        CoSimIO::Internals::DataContainer<int>& rElementTypes) override
    {
        const std::string file_name(GetFullPath("CoSimIO_mesh_" + GetConnectionName() + "_" + rIdentifier + ".vtk"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        const int num_nodes = rNodalCoordinates.size()/3;
        const int num_elems = rElementTypes.size();

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to send mesh \"" << rIdentifier << "\" with " << num_nodes << " Nodes | " << num_elems << " Elements in file \"" << file_name << "\" ..." << std::endl;

        const auto start_time(std::chrono::steady_clock::now());

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << std::scientific << std::setprecision(7); // TODO maybe this should be configurable

        // write file header
        output_file << "# vtk DataFile Version 4.0\n";
        output_file << "vtk output\n";
        output_file << "ASCII\n";
        output_file << "DATASET UNSTRUCTURED_GRID\n\n";

        // write nodes
        output_file << "POINTS " << num_nodes << " float\n";
        for (int i=0; i<num_nodes; ++i) {
            output_file << rNodalCoordinates[i*3] << " " << rNodalCoordinates[i*3+1] << " " << rNodalCoordinates[i*3+2] << "\n";
        }
        output_file << "\n";

        // get connectivity information
        int cell_list_size = 0;
        int counter = 0;
        int connectivities_offset = std::numeric_limits<int>::max(); //in paraview the connectivities start from 0, hence we have to check beforehand what is the connectivities offset
        for (int i=0; i<num_elems; ++i) {
            const int num_nodes_cell = GetNumNodesForVtkCellType(rElementTypes[i]);
            cell_list_size += num_nodes_cell + 1; // +1 for size of connectivity
            for (int j=0; j<num_nodes_cell; ++j) {
                connectivities_offset = std::min(connectivities_offset, rElementConnectivities[counter++]);
            }
        }

        KRATOS_CO_SIM_ERROR_IF(num_elems > 0 && connectivities_offset != 0) << "Connectivities have an offset of " << connectivities_offset << " which is not allowed!" << std::endl;

        // write cells connectivity
        counter = 0;
        output_file << "CELLS " << num_elems << " " << cell_list_size << "\n";
        for (int i=0; i<num_elems; ++i) {
            const int num_nodes_cell = GetNumNodesForVtkCellType(rElementTypes[i]);
            output_file << num_nodes_cell << " ";
            for (int j=0; j<num_nodes_cell; ++j) {
                output_file << (rElementConnectivities[counter++]-connectivities_offset);
                if (j<num_nodes_cell-1) output_file << " "; // not adding a whitespace after last number
            }
            output_file << "\n";
        }

        output_file << "\n";

        // write cell types
        output_file << "CELL_TYPES " << num_elems << "\n";
        for (int i=0; i<num_elems; ++i) {
            output_file << rElementTypes[i] << "\n";
        }

        output_file.close();
        MakeFileVisible(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished sending mesh" << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetPrintTiming()) << "Sending Mesh \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void SendControlSignalDetail(const std::string& rIdentifier, const CoSimIO::ControlSignal Signal) override
    {
        std::cout << "co_sim_file_communication SendControlSignalDetail rIdentifier" << rIdentifier << std::endl;
        std::cout << "co_sim_file_communication SendControlSignalDetail Signal" << static_cast<int>(Signal) << std::endl;
        const std::string file_name(GetFullPath("CoSimIO_control_signal_" + GetConnectionName() + ".dat"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to send control signal in file \"" << file_name << "\" ..." << std::endl;

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << static_cast<int>(Signal) << " " << rIdentifier;

        output_file.close();
        MakeFileVisible(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished sending control signal" << std::endl;
    }

    CoSimIO::ControlSignal RecvControlSignalDetail(std::string& rIdentifier) override
    {
        std::cout << "co_sim_file_communication RecvControlSignalDetail rIdentifier" << rIdentifier << std::endl;
        const std::string file_name(GetFullPath("CoSimIO_control_signal_" + GetConnectionName() + ".dat"));

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive control signal in file \"" << file_name << "\" ..." << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << " RecvControlSignalDetail " << std::endl;
        WaitForFile(file_name);

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        int control_signal;
        input_file >> control_signal;
        std::cout << "control-signal: " << control_signal << std::endl;
        input_file >> rIdentifier;

        RemoveFile(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving control signal" << std::endl;

        return static_cast<CoSimIO::ControlSignal>(control_signal);
    }

    std::string GetTempFileName(const std::string& rFileName)
    {
        if (mCommInFolder) {
            // TODO check this
            return std::string(rFileName).insert(mCommFolder.length()+1, ".");
        } else {
            return "." + rFileName;
        }
    }

    std::string GetFullPath(const std::string& rFileName)
    {
        if (mCommInFolder) {
            // TODO check this
            return mCommFolder + "/" + rFileName;  // using portable separator "/"
        } else {
            return rFileName;
        }
    }

    void WaitForFile(const std::string& rFileName)
    {
        std::cout << "CoSimIO WaitForFile" << std::endl;
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << " WaitForFile " << std::endl;
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Waiting for file: \"" << rFileName << "\"" << std::endl;
        while(!FileExists(rFileName)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50)); // wait 0.05s before next check
            KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>2) << "    Waiting" << std::endl;
        }
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Found file: \"" << rFileName << "\"" << std::endl;
    }

    void WaitUntilFileIsRemoved(const std::string& rFileName)
    {
        if (FileExists(rFileName)) { // only issue the wating message if the file exists initially
            KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Waiting for file: \"" << rFileName << "\" to be removed" << std::endl;
            while(FileExists(rFileName)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(50)); // wait 0.05s before next check
                KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>2) << "    Waiting" << std::endl;
            }
            KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "File: \"" << rFileName << "\" was removed" << std::endl;
        }
    }

    void MakeFileVisible(const std::string& rFinalFileName)
    {
        if (std::rename(GetTempFileName(rFinalFileName).c_str(), rFinalFileName.c_str()) != 0) {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Warning: \"" << rFinalFileName << "\" could not be made visible!" << std::endl;
        }
    }

};

} // namespace Internals
} // namespace CoSimIO

#endif /* KRATOS_CO_SIM_FILE_COMM_H_INCLUDED */
