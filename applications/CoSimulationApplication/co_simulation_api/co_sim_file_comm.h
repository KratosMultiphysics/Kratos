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

#define LOG(req_level, level) if(level>=req_level) std::cout << "[CoSimIO] "

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

static std::string GetTempFileName(const std::string& rFileName)
{
    // return std::string(rFileName).insert(CommDir.length()+1, ".");
    return "." + rFileName;
}

static std::string GetFullPath(const std::string& rFileName)
{
    // return CommDir + "/" + rFileName; // TODO check if this work in Win
}

static void WaitForFile(const std::string& rFileName, const int EchoLevel)
{
    LOG(1, EchoLevel) << "Waiting for file: \"" << rFileName << "\"" << std::endl;
    while(!FileExists(rFileName)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500)); // wait 0.5s before next check
        LOG(3, EchoLevel) << "    Waiting" << std::endl;
    }
    LOG(1, EchoLevel) << "Found file: \"" << rFileName << "\"" << std::endl;
}

static void RemoveFile(const std::string& rFileName, const int EchoLevel)
{
    if (std::remove(rFileName.c_str()) != 0) {
        LOG(0, EchoLevel) << "Warning: \"" << rFileName << "\" could not be deleted!" << std::endl;
    }
}

static void MakeFileVisible(const std::string& rFinalFileName, const int EchoLevel)
{
    if (std::rename(GetTempFileName(rFinalFileName).c_str(), rFinalFileName.c_str()) != 0) {
        LOG(0, EchoLevel) << "Warning: \"" << rFinalFileName << "\" could not be made visible!" << std::endl;
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

template<typename T>
void SendArray(const std::string& rFileName, const std::vector<T>& rArray, const int EchoLevel)
{
    const int size = rArray.size();

    LOG(0, EchoLevel) << "Attempting to send array \"" << rFileName << "\" with size: " << size << " ..." << std::endl;

    const auto start_time(std::chrono::steady_clock::now());

    std::ofstream output_file;
    output_file.open(GetTempFileName(rFileName));
    CheckStream(output_file, rFileName);

    output_file << std::scientific << std::setprecision(14); // TODO maybe this should be configurable

    output_file << size << "\n";

    for (int i=0; i<size-1; ++i) {
        output_file << rArray[i] << " ";
    }
    // TODO check if size == 0!
    output_file << rArray[size-1]; // outside to not have trailing whitespace

    output_file.close();
    MakeFileVisible(rFileName, EchoLevel);

    LOG(2, EchoLevel) << "Finished sending array" << std::endl;

    // if (PrintTiming) {
    //     // EMPIRE_API_LOG(0) << "Sending Array \"" << rFileName << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    // }
}

template<typename T>
void ReceiveArray(const std::string& rFileName, std::vector<T>& rArray, const int EchoLevel)
{
    LOG(2, EchoLevel) << "Attempting to receive array \"" << rFileName << "\" ..." << std::endl;

    WaitForFile(rFileName, EchoLevel);

    const auto start_time(std::chrono::steady_clock::now());

    std::ifstream input_file(rFileName);
    CheckStream(input_file, rFileName);

    input_file >> std::setprecision(14); // TODO maybe this should be configurable

    int size_read;
    input_file >> size_read; // the first number in the file is the size of the array

    rArray.resize(size_read);

    for (int i=0; i<size_read; ++i) {
        input_file >> rArray[i];
    }

    RemoveFile(rFileName, EchoLevel);

    LOG(2, EchoLevel) << "Finished receiving array" << std::endl;

    // if (PrintTiming) {
    //     // EMPIRE_API_LOG(0) << "Receiving Array \"" << rFileName << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    // }
}

static int GetVtkCellType(const int NumberOfNodes)
{
    if (NumberOfNodes == 1) {
        return 1;
    } else if (NumberOfNodes == 2) {
        return 3;
    } else if (NumberOfNodes == 3) {
        return 5;
    } else if (NumberOfNodes == 4) {
        return 9;
    } else {
        std::stringstream err_msg;
        err_msg << "Unsupported number of nodes/element: " << NumberOfNodes;
        throw std::runtime_error(err_msg.str());
    }
}

} // helpers namespace


class FileComm : public CoSimComm
{
public:
    explicit FileComm(const std::string& rName, SettingsType& rSettings)
        : CoSimComm(rName, rSettings)
    {
        const SettingsType default_settings = {
            {"communication_folder_name_suffix", ""},
            {"use_folder_for_communication" , "0"}
        };
        Tools::AddMissingSettings(default_settings, CoSimComm::mrSettings);

        mCommFolderSuffix = CoSimComm::mrSettings.at("communication_folder_name_suffix");
        mCommInFolder = (CoSimComm::mrSettings.at("use_folder_for_communication") == "1");
    }

private:

    std::string mCommFolderSuffix = "";
    bool mCommInFolder = false;

    bool ConnectDetail() override
    {
        return true; // nothing needed here for file-based communication
    }

    bool DisconnectDetail() override
    {
        return true; // nothing needed here for file-based communication
    }

    bool ImportDetail(DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        // EMPIRE_API_LOG(2) << "Attempting to receive mesh \"" << std::string(name) << "\" ..." << std::endl;

        // const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_mesh_" + std::string(name) + ".vtk"));

        // EMPIRE_API_helpers::WaitForFile(file_name);

        // const auto start_time(std::chrono::steady_clock::now());

        // std::ifstream input_file(file_name);
        // EMPIRE_API_helpers::CheckStream(input_file, file_name);

        // // reading file
        // std::string current_line;
        // bool nodes_read = false;

        // while (std::getline(input_file, current_line)) {
        //     // reading nodes
        //     if (current_line.find("POINTS") != std::string::npos) {
        //         if (nodes_read) throw std::runtime_error("The nodes were read already!");
        //         nodes_read = true;

        //         EMPIRE_API_helpers::ReadNumberAfterKeyword("POINTS", current_line, *numNodes);

        //         EMPIRE_API_LOG(2) << "Mesh contains " << *numNodes << " Nodes" << std::endl;

        //         // allocating memory for nodes
        //         // note that this has to be deleted by the client!
        //         EMPIRE_API_helpers::AllocateMemory(nodes, (*numNodes) * 3); // *nodes = new double[(*numNodes) * 3];
        //         EMPIRE_API_helpers::AllocateMemory(nodeIDs, *numNodes); // *nodeIDs = new int[*numNodes];

        //         for (int i=0; i<*numNodes*3; ++i) {
        //             input_file >> (*nodes)[i];
        //         }

        //         for (int i=0; i<*numNodes; ++i) {
        //             (*nodeIDs)[i] = i+1; // Node Ids have an offset of 1 from Kratos to VTK
        //         }
        //     }

        //     // reading elements
        //     if (current_line.find("CELLS") != std::string::npos) {
        //         if (!nodes_read) throw std::runtime_error("The nodes were not yet read!");

        //         int num_nodes_per_elem, node_id, cell_list_size;
        //         current_line = current_line.substr(current_line.find("CELLS") + 6); // removing "CELLS"
        //         std::istringstream line_stream(current_line);
        //         line_stream >> *numElems;
        //         line_stream >> cell_list_size;

        //         EMPIRE_API_LOG(2) << "Mesh contains " << *numElems << " Elements" << std::endl;

        //         // allocating memory for elements
        //         // note that this has to be deleted by the client!
        //         EMPIRE_API_helpers::AllocateMemory(numNodesPerElem, *numElems); // *numNodesPerElem = new int[*numElems];
        //         EMPIRE_API_helpers::AllocateMemory(elem, cell_list_size-*numElems); // *elem = new int[cell_list_size-*numElems];

        //         int counter=0;
        //         for (int i=0; i<*numElems; ++i) {
        //             input_file >> num_nodes_per_elem;
        //             (*numNodesPerElem)[i] = num_nodes_per_elem;
        //             for (int j=0; j<num_nodes_per_elem; ++j) {
        //                 input_file >> node_id;
        //                 (*elem)[counter++] = node_id+1; // Node Ids have an offset of 1 from Kratos to VTK
        //             }
        //         }
        //         break; // no further information reading required => CELL_TYPES are not used here
        //     }
        // }

        // EMPIRE_API_helpers::RemoveFile(file_name);

        // EMPIRE_API_LOG(2) << "Finished receiving mesh" << std::endl;

        // if (EMPIRE_API_helpers::PrintTiming) {
        //     EMPIRE_API_LOG(0) << "Receiving Mesh \"" << file_name << "\" took: " << EMPIRE_API_helpers::ElapsedSeconds(start_time) << " [sec]" << std::endl;
        // }

        return true;
    }

    bool ExportDetail(const DataContainers::Mesh& rDataContainer, const std::string& rIdentifier) override
    {
        // EMPIRE_API_LOG(2) << "Attempting to send mesh \"" << std::string(name) << "\" with " << numNodes << " Nodes | " << numElems << " Elements ..." << std::endl;

        // const auto start_time(std::chrono::steady_clock::now());

        // const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_mesh_" + std::string(name) + ".vtk"));

        // std::ofstream output_file;
        // output_file.open(EMPIRE_API_helpers::GetTempFileName(file_name));
        // EMPIRE_API_helpers::CheckStream(output_file, file_name);

        // output_file << std::scientific << std::setprecision(7); // TODO maybe this should be configurable

        // // write file header
        // output_file << "# vtk DataFile Version 4.0\n";
        // output_file << "vtk output\n";
        // output_file << "ASCII\n";
        // output_file << "DATASET UNSTRUCTURED_GRID\n\n";

        // // write nodes
        // int vtk_id = 0;
        // std::unordered_map<int, int> node_vtk_id_map;
        // output_file << "POINTS " << numNodes << " float\n";
        // for (int i=0; i<numNodes; ++i) {
        //     output_file << nodes[i*3] << " " << nodes[i*3+1] << " " << nodes[i*3+2] << "\n";
        //     node_vtk_id_map[nodeIDs[i]] = vtk_id++;
        // }
        // output_file << "\n";

        // // write cells connectivity
        // int cell_list_size = 0;
        // for (int i=0; i<numElems; ++i) {
        //     cell_list_size += numNodesPerElem[i] + 1;
        // }

        // int counter=0;
        // output_file << "CELLS " << numElems << " " << cell_list_size << "\n";
        // for (int i=0; i<numElems; ++i) {
        //     const int num_nodes_elem = numNodesPerElem[i];
        //     output_file << num_nodes_elem << " ";
        //     for (int j=0; j<num_nodes_elem; ++j) {
        //         output_file << node_vtk_id_map.at(elems[counter++]);
        //         if (j<num_nodes_elem-1) output_file << " "; // not adding a whitespace after last number
        //     }
        //     output_file << "\n";
        // }

        // output_file << "\n";

        // // write cell types
        // output_file << "CELL_TYPES " << numElems << "\n";
        // for (int i=0; i<numElems; ++i) {
        //     output_file << EMPIRE_API_helpers::GetVtkCellType(numNodesPerElem[i]) << "\n";
        // }

        // output_file.close();
        // EMPIRE_API_helpers::MakeFileVisible(file_name);

        // EMPIRE_API_LOG(2) << "Finished sending mesh" << std::endl;

        // if (EMPIRE_API_helpers::PrintTiming) {
        //     EMPIRE_API_LOG(0) << "Sending Mesh \"" << file_name << "\" took: " << EMPIRE_API_helpers::ElapsedSeconds(start_time) << " [sec]" << std::endl;
        // }

        return true;
    }

    bool ImportDetail(DataContainers::Data& rDataContainer, const std::string& rIdentifier) override
    {
        ReceiveArray(rIdentifier, rDataContainer.data, mEchoLevel);
        return true;
    }

    bool ExportDetail(const DataContainers::Data& rDataContainer, const std::string& rIdentifier) override
    {
        SendArray(rIdentifier, rDataContainer.data, mEchoLevel);
        return true;
    }

    bool ImportDetail(int& rDataContainer, const std::string& rIdentifier) override
    {
        std::vector<int> data_vec;
        ReceiveArray(rIdentifier, data_vec, mEchoLevel);
        rDataContainer = data_vec[0];
        return true;
    }

    bool ExportDetail(const int& rDataContainer, const std::string& rIdentifier) override
    {
        std::vector<int> data_vec = {rDataContainer};
        SendArray(rIdentifier, data_vec, mEchoLevel);
        return true;
    }

};

} // namespace CoSim

#endif /* KRATOS_CO_SIM_FILE_COMM_H_INCLUDED */