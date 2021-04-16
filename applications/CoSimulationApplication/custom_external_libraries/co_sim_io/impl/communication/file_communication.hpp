//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_FILE_COMMUNICATION_H_INCLUDED
#define CO_SIM_IO_FILE_COMMUNICATION_H_INCLUDED

// System includes
#include <chrono>
#include <thread>
#include <iomanip>
#include <algorithm>

// std::filesystem is needed for file communication in a folder
// std::filesystem is part of C++17 and not supported by every compiler. Here we check if it is available.
#if defined(__cplusplus) && __cplusplus >= 201703L
    #if defined(__has_include) && __has_include(<filesystem>) // has_include is C++17, hence has to be checked in a separate line
        #define CO_SIM_IO_FILESYSTEM_AVAILABLE
        #include <filesystem>
        namespace fs = std::filesystem;
    #elif __has_include(<experimental/filesystem>)
        #define CO_SIM_IO_FILESYSTEM_AVAILABLE
        #include <experimental/filesystem>
        namespace fs = std::experimental::filesystem;
    #endif
#endif


// Project includes
#include "communication.hpp"

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
        CO_SIM_IO_INFO("CoSimIO") << "Warning: \"" << rFileName << "\" could not be deleted!" << std::endl;
    }
}

template <typename T>
static void CheckStream(const T& rStream, const std::string& rFileName)
{
    CO_SIM_IO_ERROR_IF_NOT(rStream.is_open()) << rFileName << " could not be opened!" << std::endl;
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
        CO_SIM_IO_ERROR << "Unsupported cell type: " << VtkCellType << std::endl;
        return 0;
    }
}

} // helpers namespace


class FileCommunication : public Communication
{
public:
    explicit FileCommunication(const std::string& rName, const Info& I_Settings, const bool IsConnectionMaster)
        : Communication(rName, I_Settings, IsConnectionMaster)
    {
        if (I_Settings.Has("use_folder_for_communication")) {
            mCommInFolder = I_Settings.Get<bool>("use_folder_for_communication");
        }

        mCommFolder = ".CoSimIOFileComm_"+rName;

        #ifndef CO_SIM_IO_FILESYSTEM_AVAILABLE
        CO_SIM_IO_ERROR_IF(mCommInFolder) << "Communication is a folder can only be used if std::filesystem (C++17) is available" << std::endl;
        #endif

        if (mCommInFolder && GetIsConnectionMaster()) {
            #ifdef CO_SIM_IO_FILESYSTEM_AVAILABLE
            // delete and recreate directory to remove potential leftovers
            fs::remove_all(mCommFolder);
            fs::create_directory(mCommFolder);
            #endif
        }
    }

    ~FileCommunication() override
    {
        if (GetIsConnected()) {
            CO_SIM_IO_INFO("CoSimIO") << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
            Disconnect();
        }
    }

private:

    std::string mCommFolder = "";
    bool mCommInFolder = false;

    bool ConnectDetail() override
    {
        return true; // nothing needed here for file-based communication (maybe do sth here?)
        // master could write a file that gets deleted by slave to aknowledge connection... Probably not a bad idea! => slave returns once it found and deleted file, master waits for deletion of file
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

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive array \"" << rIdentifier << "\" in file \"" << file_name << "\" ..." << std::endl;

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

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving array with size: " << size_read << std::endl;

        CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Receiving Array \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ExportDataImpl(
        const std::string& rIdentifier,
        const CoSimIO::Internals::DataContainer<double>& rData) override
    {
        const std::string file_name(GetFullPath("CoSimIO_data_" + GetConnectionName() + "_" + rIdentifier + ".dat"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        const std::size_t size = rData.size();
        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to send array \"" << rIdentifier << "\" with size: " << size << " in file \"" << file_name << "\" ..." << std::endl;

        const auto start_time(std::chrono::steady_clock::now());

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << std::scientific << std::setprecision(14); // TODO maybe this should be configurable

        output_file << size << "\n";

        for (std::size_t i=0; i<size-1; ++i) {
            output_file << rData[i] << " ";
        }
        // TODO check if size == 0!
        output_file << rData[size-1]; // outside to not have trailing whitespace

        output_file.close();
        MakeFileVisible(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished sending array" << std::endl;

        CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Sending Array \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ImportMeshImpl(
        const std::string& rIdentifier,
        CoSimIO::Internals::DataContainer<double>& rNodalCoordinates,
        CoSimIO::Internals::DataContainer<int>& rElementConnectivities,
        CoSimIO::Internals::DataContainer<int>& rElementTypes) override
    {
        const std::string file_name(GetFullPath("CoSimIO_mesh_" + GetConnectionName() + "_" + rIdentifier + ".vtk"));

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive mesh \"" << rIdentifier << "\" in file \"" << file_name << "\" ..." << std::endl;

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
                CO_SIM_IO_ERROR_IF(nodes_read) << "The nodes were read already!" << std::endl;
                CO_SIM_IO_ERROR_IF(cells_read) << "The cells were read already!" << std::endl;
                nodes_read = true;

                int num_nodes;
                current_line = current_line.substr(current_line.find("POINTS") + 7); // removing "POINTS"
                std::istringstream line_stream(current_line);
                line_stream >> num_nodes;

                CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Mesh contains " << num_nodes << " Nodes" << std::endl;

                rNodalCoordinates.resize(3*num_nodes);

                for (int i=0; i<num_nodes*3; ++i) {
                    input_file >> rNodalCoordinates[i];
                }
            }

            // reading cells
            if (current_line.find("CELLS") != std::string::npos) {
                CO_SIM_IO_ERROR_IF_NOT(nodes_read) << "The nodes were not yet read!" << std::endl;
                CO_SIM_IO_ERROR_IF(cells_read) << "The cells were read already!" << std::endl;
                cells_read = true;

                int num_nodes_per_cell, num_cells, elem_conn, cell_list_size;
                current_line = current_line.substr(current_line.find("CELLS") + 6); // removing "CELLS"
                std::istringstream line_stream(current_line);
                line_stream >> num_cells;
                line_stream >> cell_list_size;

                rElementConnectivities.resize(cell_list_size-num_cells); // the first in number in each line is the number of connectivities, which is not needed bcs it can be derived form the elements-type
                rElementTypes.resize(num_cells);

                CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Mesh contains " << num_cells << " Elements" << std::endl;

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
                CO_SIM_IO_ERROR_IF_NOT(nodes_read) << "The nodes were not yet read!" << std::endl;
                CO_SIM_IO_ERROR_IF_NOT(cells_read) << "The cells were not yet read!" << std::endl;

                for (std::size_t i=0; i<rElementTypes.size(); ++i) { // rElementTypes was resized to correct size above
                    input_file >> rElementTypes[i];
                }

            }
        }

        RemoveFile(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving mesh" << std::endl;

        CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Receiving Mesh \"" << file_name << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ExportMeshImpl(
        const std::string& rIdentifier,
        const CoSimIO::Internals::DataContainer<double>& rNodalCoordinates,
        const CoSimIO::Internals::DataContainer<int>& rElementConnectivities,
        const CoSimIO::Internals::DataContainer<int>& rElementTypes) override
    {
        const std::string file_name(GetFullPath("CoSimIO_mesh_" + GetConnectionName() + "_" + rIdentifier + ".vtk"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        const std::size_t num_nodes = rNodalCoordinates.size()/3;
        const std::size_t num_elems = rElementTypes.size();

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to send mesh \"" << rIdentifier << "\" with " << num_nodes << " Nodes | " << num_elems << " Elements in file \"" << file_name << "\" ..." << std::endl;

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
        for (std::size_t i=0; i<num_nodes; ++i) {
            output_file << rNodalCoordinates[i*3] << " " << rNodalCoordinates[i*3+1] << " " << rNodalCoordinates[i*3+2] << "\n";
        }
        output_file << "\n";

        // get connectivity information
        std::size_t cell_list_size = 0;
        std::size_t counter = 0;
        int connectivities_offset = std::numeric_limits<int>::max(); //in paraview the connectivities start from 0, hence we have to check beforehand what is the connectivities offset
        for (std::size_t i=0; i<num_elems; ++i) {
            const std::size_t num_nodes_cell = GetNumNodesForVtkCellType(rElementTypes[i]);
            cell_list_size += num_nodes_cell + 1; // +1 for size of connectivity
            for (std::size_t j=0; j<num_nodes_cell; ++j) {
                connectivities_offset = std::min(connectivities_offset, rElementConnectivities[counter++]);
            }
        }

        CO_SIM_IO_ERROR_IF(num_elems > 0 && connectivities_offset != 0) << "Connectivities have an offset of " << connectivities_offset << " which is not allowed!" << std::endl;

        // write cells connectivity
        counter = 0;
        output_file << "CELLS " << num_elems << " " << cell_list_size << "\n";
        for (std::size_t i=0; i<num_elems; ++i) {
            const std::size_t num_nodes_cell = GetNumNodesForVtkCellType(rElementTypes[i]);
            output_file << num_nodes_cell << " ";
            for (std::size_t j=0; j<num_nodes_cell; ++j) {
                output_file << (rElementConnectivities[counter++]-connectivities_offset);
                if (j<num_nodes_cell-1) output_file << " "; // not adding a whitespace after last number
            }
            output_file << "\n";
        }

        output_file << "\n";

        // write cell types
        output_file << "CELL_TYPES " << num_elems << "\n";
        for (std::size_t i=0; i<num_elems; ++i) {
            output_file << rElementTypes[i] << "\n";
        }

        output_file.close();
        MakeFileVisible(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished sending mesh" << std::endl;

        CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Sending Mesh \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void SendControlSignalDetail(const std::string& rIdentifier, const CoSimIO::ControlSignal Signal) override
    {
        const std::string file_name(GetFullPath("CoSimIO_control_signal_" + GetConnectionName() + ".dat"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to send control signal in file \"" << file_name << "\" ..." << std::endl;

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << static_cast<int>(Signal) << " " << rIdentifier;

        output_file.close();
        MakeFileVisible(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished sending control signal" << std::endl;
    }

    CoSimIO::ControlSignal RecvControlSignalDetail(std::string& rIdentifier) override
    {
        const std::string file_name(GetFullPath("CoSimIO_control_signal_" + GetConnectionName() + ".dat"));

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive control signal in file \"" << file_name << "\" ..." << std::endl;

        WaitForFile(file_name);

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        int control_signal;
        input_file >> control_signal;
        input_file >> rIdentifier;

        RemoveFile(file_name);

        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving control signal" << std::endl;

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
        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Waiting for file: \"" << rFileName << "\"" << std::endl;
        while(!FileExists(rFileName)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(50)); // wait 0.05s before next check
            CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>2) << "    Waiting" << std::endl;
        }
        CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Found file: \"" << rFileName << "\"" << std::endl;
    }

    void WaitUntilFileIsRemoved(const std::string& rFileName)
    {
        if (FileExists(rFileName)) { // only issue the wating message if the file exists initially
            CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Waiting for file: \"" << rFileName << "\" to be removed" << std::endl;
            while(FileExists(rFileName)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(50)); // wait 0.05s before next check
                CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>2) << "    Waiting" << std::endl;
            }
            CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>0) << "File: \"" << rFileName << "\" was removed" << std::endl;
        }
    }

    void MakeFileVisible(const std::string& rFinalFileName)
    {
        if (std::rename(GetTempFileName(rFinalFileName).c_str(), rFinalFileName.c_str()) != 0) {
            CO_SIM_IO_INFO("CoSimIO") << "Warning: \"" << rFinalFileName << "\" could not be made visible!" << std::endl;
        }
    }

};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_FILE_COMMUNICATION_H_INCLUDED