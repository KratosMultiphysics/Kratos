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
#include "co_sim_comm.h"

namespace CoSim {

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
        CS_LOG << "Warning: \"" << rFileName << "\" could not be deleted!" << std::endl;
    }
}

template <typename T>
static void CheckStream(const T& rStream, const std::string& rFileName)
{
    if (!rStream.is_open()) {
        std::stringstream err_msg;
        err_msg << rFileName << " could not be opened!";
        throw std::runtime_error(err_msg.str());
    }
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
        throw std::runtime_error("Unsupported cell type: "+std::to_string(VtkCellType));
    }
}

} // helpers namespace


class FileComm : public CoSimComm
{
public:
    explicit FileComm(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster)
        : CoSimComm(rName, rSettings, IsConnectionMaster)
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

    ~FileComm() override
    {
        if (GetIsConnected()) {
            CS_LOG << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
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

    bool ImportDetail(DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        const std::string file_name(GetFullPath("CoSimIO_mesh_" + GetName() + "_" + rIdentifier + ".vtk"));

        CS_LOG_IF(GetEchoLevel()>1) << "Attempting to receive mesh \"" << rIdentifier << "\" in file \"" << file_name << "\" ..." << std::endl;

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
                if (nodes_read) throw std::runtime_error("The nodes were read already!");
                nodes_read = true;

                int num_nodes;
                current_line = current_line.substr(current_line.find("POINTS") + 7); // removing "POINTS"
                std::istringstream line_stream(current_line);
                line_stream >> num_nodes;

                CS_LOG_IF(GetEchoLevel()>1) << "Mesh contains " << num_nodes << " Nodes" << std::endl;

                rDataContainer.node_coords.resize(3*num_nodes);

                for (int i=0; i<num_nodes*3; ++i) {
                    input_file >> rDataContainer.node_coords[i];
                }
            }

            // reading cells
            if (current_line.find("CELLS") != std::string::npos) {
                if (!nodes_read) throw std::runtime_error("The nodes were not yet read!");

                int num_nodes_per_cell, num_cells, node_id, cell_list_size;
                current_line = current_line.substr(current_line.find("CELLS") + 6); // removing "CELLS"
                std::istringstream line_stream(current_line);
                line_stream >> num_cells;
                line_stream >> cell_list_size;

                CS_LOG_IF(GetEchoLevel()>1) << "Mesh contains " << num_cells << " Cells" << std::endl;

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

        CS_LOG_IF(GetEchoLevel()>1) << "Finished receiving mesh" << std::endl;

        CS_LOG_IF(GetPrintTiming()) << "Receiving Mesh \"" << file_name << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;

        return true;
    }

    bool ExportDetail(const DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        const std::string file_name(GetFullPath("CoSimIO_mesh_" + GetName() + "_" + rIdentifier + ".vtk"));

        const int num_nodes = rDataContainer.node_coords.size()/3;
        const int num_cells = rDataContainer.cell_types.size();

        CS_LOG_IF(GetEchoLevel()>1) << "Attempting to send mesh \"" << rIdentifier << "\" with " << num_nodes << " Nodes | " << num_cells << " Cells in file \"" << file_name << "\" ..." << std::endl;

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
            output_file << rDataContainer.node_coords[i*3] << " " << rDataContainer.node_coords[i*3+1] << " " << rDataContainer.node_coords[i*3+2] << "\n";
        }
        output_file << "\n";

        // write cells connectivity
        int counter=0;
        output_file << "CELLS " << num_cells << " " << rDataContainer.connectivities.size() << "\n";
        for (int i=0; i<num_cells; ++i) {
            const int num_nodes_cell = GetNumNodesForVtkCellType(rDataContainer.cell_types[i]);
            output_file << num_nodes_cell << " ";
            for (int j=0; j<num_nodes_cell; ++j) {
                output_file << rDataContainer.connectivities[counter++];
                if (j<num_nodes_cell-1) output_file << " "; // not adding a whitespace after last number
            }
            output_file << "\n";
        }

        output_file << "\n";

        // write cell types
        output_file << "CELL_TYPES " << num_cells << "\n";
        for (int i=0; i<num_cells; ++i) {
            output_file << rDataContainer.cell_types[i] << "\n";
        }

        output_file.close();
        MakeFileVisible(file_name);

        CS_LOG_IF(GetEchoLevel()>1) << "Finished sending mesh" << std::endl;

        CS_LOG_IF(GetPrintTiming()) << "Sending Mesh \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;

        return true;
    }

    bool ImportDetail(DataContainers::Data& rDataContainer, const std::string& rIdentifier) override
    {
        ReceiveArray(rIdentifier, rDataContainer.data);
        return true;
    }

    bool ExportDetail(const DataContainers::Data& rDataContainer, const std::string& rIdentifier) override
    {
        SendArray(rIdentifier, rDataContainer.data);
        return true;
    }

    void SendControlSignalDetail(const int rSignal, const std::string& rIdentifier) override
    {
        const std::string file_name(GetFullPath("CoSimIO_control_signal_" + GetName() + ".dat"));

        WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

        CS_LOG_IF(GetEchoLevel()>1) << "Attempting to send control signal in file \"" << file_name << "\" ..." << std::endl;

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << rSignal << " " << rIdentifier;

        output_file.close();
        MakeFileVisible(file_name);

        CS_LOG_IF(GetEchoLevel()>1) << "Finished sending control signal" << std::endl;
    }

    int RecvControlSignalDetail(std::string& rIdentifier) override
    {
        const std::string file_name(GetFullPath("CoSimIO_control_signal_" + GetName() + ".dat"));

        CS_LOG_IF(GetEchoLevel()>1) << "Attempting to receive control signal in file \"" << file_name << "\" ..." << std::endl;

        WaitForFile(file_name);

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        int control_signal;
        input_file >> control_signal;
        input_file >> rIdentifier;

        RemoveFile(file_name);

        CS_LOG_IF(GetEchoLevel()>1) << "Finished receiving control signal" << std::endl;

        return control_signal;
    }


    template<typename T>
    void SendArray(const std::string& rIdentifier, const std::vector<T>& rArray)
    {
        const std::string file_name(GetFullPath("CoSimIO_data_" + GetName() + "_" + rIdentifier + ".dat"));

        const int size = rArray.size();

        CS_LOG_IF(GetEchoLevel()>1) << "Attempting to send array \"" << rIdentifier << "\" with size: " << size << " in file \"" << file_name << "\" ..." << std::endl;

        const auto start_time(std::chrono::steady_clock::now());

        std::ofstream output_file;
        output_file.open(GetTempFileName(file_name));
        CheckStream(output_file, file_name);

        output_file << std::scientific << std::setprecision(14); // TODO maybe this should be configurable

        output_file << size << "\n";

        for (int i=0; i<size-1; ++i) {
            output_file << rArray[i] << " ";
        }
        // TODO check if size == 0!
        output_file << rArray[size-1]; // outside to not have trailing whitespace

        output_file.close();
        MakeFileVisible(file_name);

        CS_LOG_IF(GetEchoLevel()>1) << "Finished sending array" << std::endl;

        CS_LOG_IF(GetPrintTiming()) << "Sending Array \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }

    template<typename T>
    void ReceiveArray(const std::string& rIdentifier, std::vector<T>& rArray)
    {
        const std::string file_name(GetFullPath("CoSimIO_data_" + GetName() + "_" + rIdentifier + ".dat"));

        CS_LOG_IF(GetEchoLevel()>1) << "Attempting to receive array \"" << rIdentifier << "\" in file \"" << file_name << "\" ..." << std::endl;

        WaitForFile(file_name);

        const auto start_time(std::chrono::steady_clock::now());

        std::ifstream input_file(file_name);
        CheckStream(input_file, file_name);

        input_file >> std::setprecision(14); // TODO maybe this should be configurable

        int size_read;
        input_file >> size_read; // the first number in the file is the size of the array

        rArray.resize(size_read);

        for (int i=0; i<size_read; ++i) {
            input_file >> rArray[i];
        }

        RemoveFile(file_name);

        CS_LOG_IF(GetEchoLevel()>1) << "Finished receiving array" << std::endl;

        CS_LOG_IF(GetPrintTiming()) << "Receiving Array \"" << rIdentifier << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
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
        CS_LOG_IF(GetEchoLevel()>0) << "Waiting for file: \"" << rFileName << "\"" << std::endl;
        while(!FileExists(rFileName)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(500)); // wait 0.5s before next check
            CS_LOG_IF(GetEchoLevel()>2) << "    Waiting" << std::endl;
        }
        CS_LOG_IF(GetEchoLevel()>0) << "Found file: \"" << rFileName << "\"" << std::endl;
    }

    void WaitUntilFileIsRemoved(const std::string& rFileName)
    {
        CS_LOG_IF(GetEchoLevel()>0) << "Waiting for file: \"" << rFileName << "\" to be removed" << std::endl;
        while(FileExists(rFileName)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(500)); // wait 0.5s before next check
            CS_LOG_IF(GetEchoLevel()>2) << "    Waiting" << std::endl;
        }
        CS_LOG_IF(GetEchoLevel()>0) << "File: \"" << rFileName << "\" was removed" << std::endl;
    }

    void MakeFileVisible(const std::string& rFinalFileName)
    {
        if (std::rename(GetTempFileName(rFinalFileName).c_str(), rFinalFileName.c_str()) != 0) {
            CS_LOG << "Warning: \"" << rFinalFileName << "\" could not be made visible!" << std::endl;
        }
    }

};

} // namespace CoSim

#endif /* KRATOS_CO_SIM_FILE_COMM_H_INCLUDED */