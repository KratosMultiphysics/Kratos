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

#ifndef KRATOS_CO_SIM_FILE_COMM_H_INCLUDED
#define KRATOS_CO_SIM_FILE_COMM_H_INCLUDED

// System includes
#include <chrono>
#include <thread>
#include <iomanip>

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
    const std::unordered_map<int, int> vtk_cell_type_map = {
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
        const SettingsType default_settings = {
            {"use_folder_for_communication" , "0"}
        };
        Internals::AddMissingSettings(default_settings, mrSettings);

        mCommFolder = ".CoSimIOFileComm_"+rName;
        // mCommInFolder = (mrSettings.at("use_folder_for_communication") == "1"); // this is not yet supported

        if (mCommInFolder && GetIsConnectionMaster()) {
            // delete and recreate folder
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
        int& rSize,
        CoSimIO::Internals::DataContainer<double>& rData) override
    {
        const std::string file_name(GetFullPath("CoSimIO_data_" + GetName() + "_" + rIdentifier + ".dat"));

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive array \"" << rIdentifier << "\" in file \"" << file_name << "\" ..." << std::endl;

        WaitForFile(file_name);

        const auto start_time(std::chrono::steady_clock::now());

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        input_file >> std::setprecision(14); // TODO maybe this should be configurable

        int size_read;
        input_file >> size_read; // the first number in the file is the size of the array
        rSize = size_read;

        // here if incoming size is larger than allocated size => resize
        rData.resize_if_smaller(size_read); // TODO do we want to go this way???

        for (int i=0; i<size_read; ++i) {
            input_file >> rData[i];
        }

        RemoveFile(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving array" << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetPrintTiming()) << "Receiving Array \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ExportDataImpl(
        const std::string& rIdentifier,
        const int Size,
        const CoSimIO::Internals::DataContainer<double>& rData) override
    {
        const std::string file_name(GetFullPath("CoSimIO_data_" + GetName() + "_" + rIdentifier + ".dat"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to send array \"" << rIdentifier << "\" with size: " << Size << " in file \"" << file_name << "\" ..." << std::endl;

        const auto start_time(std::chrono::steady_clock::now());

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << std::scientific << std::setprecision(14); // TODO maybe this should be configurable

        output_file << Size << "\n";

        for (int i=0; i<Size-1; ++i) {
            output_file << rData[i] << " ";
        }
        // TODO check if size == 0!
        output_file << rData[Size-1]; // outside to not have trailing whitespace

        output_file.close();
        MakeFileVisible(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished sending array" << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetPrintTiming()) << "Sending Array \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ImportMeshImpl(
        const std::string& rIdentifier,
        int& rNumberOfNodes,
        int& rNumberOfElements,
        CoSimIO::Internals::DataContainer<double>& rNodalCoordinates,
        CoSimIO::Internals::DataContainer<int>& rElementConnectivities,
        CoSimIO::Internals::DataContainer<int>& rElementTypes) override
    {
        const std::string file_name(GetFullPath("CoSimIO_mesh_" + GetName() + "_" + rIdentifier + ".vtk"));

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive mesh \"" << rIdentifier << "\" in file \"" << file_name << "\" ..." << std::endl;

        WaitForFile(file_name);

        const auto start_time(std::chrono::steady_clock::now());

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        // reading file
        std::string current_line;
        bool nodes_read = false;

        while (std::getline(input_file, current_line)) {
            // reading nodes
            if (current_line.find("POINTS") != std::string::npos) {
                KRATOS_CO_SIM_ERROR_IF(nodes_read) << "The nodes were read already!" << std::endl;
                nodes_read = true;

                int num_nodes;
                current_line = current_line.substr(current_line.find("POINTS") + 7); // removing "POINTS"
                std::istringstream line_stream(current_line);
                line_stream >> num_nodes;
                rNumberOfNodes = num_nodes;

                KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Mesh contains " << num_nodes << " Nodes" << std::endl;

                rNodalCoordinates.resize_if_smaller(3*num_nodes);

                for (int i=0; i<num_nodes*3; ++i) {
                    input_file >> rNodalCoordinates[i];
                }
            }

            // reading cells
            if (current_line.find("CELLS") != std::string::npos) {
                KRATOS_CO_SIM_ERROR_IF_NOT(nodes_read) << "The nodes were not yet read!" << std::endl;

                int /*num_nodes_per_cell,*/ num_cells, /*node_id,*/ cell_list_size;
                current_line = current_line.substr(current_line.find("CELLS") + 6); // removing "CELLS"
                std::istringstream line_stream(current_line);
                line_stream >> num_cells;
                line_stream >> cell_list_size;
                rNumberOfElements = num_cells;

                KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Mesh contains " << num_cells << " Cells" << std::endl;

                // int counter=0;
                // for (int i=0; i<num_cells; ++i) {
                //     input_file >> num_nodes_per_cell;
                //     (*numNodesPerElem)[i] = num_nodes_per_cell;
                //     for (int j=0; j<num_nodes_per_cell; ++j) {
                //         input_file >> node_id;
                //         (*elem)[counter++] = node_id+1; // Node Ids have an offset of 1 from Kratos to VTK
                //     }
                // }
                break; // no further information reading required => CELL_TYPES are not used here
            }
        }

        RemoveFile(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving mesh" << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetPrintTiming()) << "Receiving Mesh \"" << file_name << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void ExportMeshImpl(
        const std::string& rIdentifier,
        const int NumberOfNodes,
        const int NumberOfElements,
        CoSimIO::Internals::DataContainer<double>& rNodalCoordinates,
        CoSimIO::Internals::DataContainer<int>& rElementConnectivities,
        CoSimIO::Internals::DataContainer<int>& rElementTypes) override
    {
        const std::string file_name(GetFullPath("CoSimIO_mesh_" + GetName() + "_" + rIdentifier + ".vtk"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to send mesh \"" << rIdentifier << "\" with " << NumberOfNodes << " Nodes | " << NumberOfElements << " Cells in file \"" << file_name << "\" ..." << std::endl;

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
        output_file << "POINTS " << NumberOfNodes << " float\n";
        for (int i=0; i<NumberOfNodes; ++i) {
            output_file << rNodalCoordinates[i*3] << " " << rNodalCoordinates[i*3+1] << " " << rNodalCoordinates[i*3+2] << "\n";
        }
        output_file << "\n";

        // write cells connectivity
        int cell_list_size = 0;
        for (int i=0; i<NumberOfElements; ++i) {
            cell_list_size += GetNumNodesForVtkCellType(rElementTypes[i]) + 1;
        }

        // write cells connectivity
        int counter=0;
        output_file << "CELLS " << NumberOfElements << " " << cell_list_size << "\n";
        for (int i=0; i<NumberOfElements; ++i) {
            const int num_nodes_cell = GetNumNodesForVtkCellType(rElementTypes[i]);
            output_file << num_nodes_cell << " ";
            for (int j=0; j<num_nodes_cell; ++j) {
                output_file << rElementConnectivities[counter++];
                if (j<num_nodes_cell-1) output_file << " "; // not adding a whitespace after last number
            }
            output_file << "\n";
        }

        output_file << "\n";

        // write cell types
        output_file << "CELL_TYPES " << NumberOfElements << "\n";
        for (int i=0; i<NumberOfElements; ++i) {
            output_file << rElementTypes[i] << "\n";
        }

        output_file.close();
        MakeFileVisible(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished sending mesh" << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetPrintTiming()) << "Sending Mesh \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    void SendControlSignalDetail(Internals::ControlSignal Signal, const std::string& rIdentifier) override
    {
        const std::string file_name(GetFullPath("CoSimIO_control_signal_" + GetName() + ".dat"));

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

    Internals::ControlSignal RecvControlSignalDetail(std::string& rIdentifier) override
    {
        const std::string file_name(GetFullPath("CoSimIO_control_signal_" + GetName() + ".dat"));

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to receive control signal in file \"" << file_name << "\" ..." << std::endl;

        WaitForFile(file_name);

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        int control_signal;
        input_file >> control_signal;
        input_file >> rIdentifier;

        RemoveFile(file_name);

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished receiving control signal" << std::endl;

        return static_cast<Internals::ControlSignal>(control_signal);
    }

    std::string GetTempFileName(const std::string& rFileName)
    {
        // return std::string(rFileName).insert(CommDir.length()+1, ".");
        return "." + rFileName;
    }

    std::string GetFullPath(const std::string& rFileName)
    {
        // return CommDir + "/" + rFileName; // TODO check if this work in Win
        return rFileName;
    }

    void WaitForFile(const std::string& rFileName)
    {
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Waiting for file: \"" << rFileName << "\"" << std::endl;
        while(!FileExists(rFileName)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(500)); // wait 0.5s before next check
            KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>2) << "    Waiting" << std::endl;
        }
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Found file: \"" << rFileName << "\"" << std::endl;
    }

    void WaitUntilFileIsRemoved(const std::string& rFileName)
    {
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Waiting for file: \"" << rFileName << "\" to be removed" << std::endl;
        while(FileExists(rFileName)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(500)); // wait 0.5s before next check
            KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>2) << "    Waiting" << std::endl;
        }
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "File: \"" << rFileName << "\" was removed" << std::endl;
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