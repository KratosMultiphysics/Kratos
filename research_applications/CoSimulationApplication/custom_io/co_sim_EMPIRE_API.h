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

#ifndef KRATOS_CO_SIM_EMPIRE_API_H_INCLUDED
#define KRATOS_CO_SIM_EMPIRE_API_H_INCLUDED

/*
This file has the same API as EMPIRE (defined in "EMPIRE_API.h"), hence it can be included instead of EMPIRE
It is used FileIO for data-exchange, in VTK-format
Note:
- This file cannot have Kratos-includes, because it is also included in other codes!
- This file is intended to be header-only, such that other codes do not have to link against a library
- Requires c++11 to compile
- All functions have to be static, otherwise linking errors occur if it is included in multiple files
*/

// System includes
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <chrono>
#include <thread>
#include <unordered_map>
#include <chrono>
#include <vector>

namespace EMPIRE_API_helpers {

// Some options that can be configured
static const std::string CommDir = ".EmpireIO";
static int PrintTiming = 0;
static int EchoLevel = 1;

#define EMPIRE_API_LOG(level) if(EMPIRE_API_helpers::EchoLevel>=level) std::cout << "[EMPIRE_API] "

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
    return std::string(rFileName).insert(CommDir.length()+1, ".");
}

static std::string GetFullPath(const std::string& rFileName)
{
    return CommDir + "/" + rFileName; // TODO check if this work in Win
}

static void WaitForFile(const std::string& rFileName)
{
    EMPIRE_API_LOG(1) << "Waiting for file: \"" << rFileName << "\"" << std::endl;
    while(!FileExists(rFileName)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500)); // wait 0.5s before next check
        EMPIRE_API_LOG(3) << "    Waiting" << std::endl;
    }
    EMPIRE_API_LOG(1) << "Found file: \"" << rFileName << "\"" << std::endl;
}

static void RemoveFile(const std::string& rFileName)
{
    if (std::remove(rFileName.c_str()) != 0) {
        EMPIRE_API_LOG(0) << "Warning: \"" << rFileName << "\" could not be deleted!" << std::endl;
    }
}

static void MakeFileVisible(const std::string& rFinalFileName)
{
    if (std::rename(GetTempFileName(rFinalFileName).c_str(), rFinalFileName.c_str()) != 0) {
        EMPIRE_API_LOG(0) << "Warning: \"" << rFinalFileName << "\" could not be made visible!" << std::endl;
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

static void SendArray(const std::string& rFileName, const int sizeOfArray, const double *data)
{
    EMPIRE_API_LOG(2) << "Attempting to send array \"" << rFileName << "\" with size: " << sizeOfArray << " ..." << std::endl;

    const auto start_time(std::chrono::steady_clock::now());

    std::ofstream output_file;
    output_file.open(GetTempFileName(rFileName));
    CheckStream(output_file, rFileName);

    output_file << std::scientific << std::setprecision(14); // TODO maybe this should be configurable

    output_file << sizeOfArray << "\n";

    for (int i=0; i<sizeOfArray-1; ++i) {
        output_file << data[i] << " ";
    }
    output_file << data[sizeOfArray-1]; // outside to not have trailing whitespace

    output_file.close();
    MakeFileVisible(rFileName);

    EMPIRE_API_LOG(2) << "Finished sending array" << std::endl;

    if (PrintTiming) {
        EMPIRE_API_LOG(0) << "Sending Array \"" << rFileName << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }
}

static void ReceiveArray(const std::string& rFileName, const int sizeOfArray, double *data)
{
    EMPIRE_API_LOG(2) << "Attempting to receive array \"" << rFileName << "\" with size: " << sizeOfArray << " ..." << std::endl;

    WaitForFile(rFileName);

    const auto start_time(std::chrono::steady_clock::now());

    std::ifstream input_file(rFileName);
    CheckStream(input_file, rFileName);

    input_file >> std::setprecision(14); // TODO maybe this should be configurable

    int size_read;
    input_file >> size_read; // the first number in the file is the size of the array

    if (size_read != sizeOfArray) {
        std::stringstream err_msg;
        err_msg << "The received size for array \"" << rFileName << "\" is different from what is expected:";
        err_msg << "\n    Expected size: " << sizeOfArray;
        err_msg << "\n    Received size: " << size_read;
        throw std::runtime_error(err_msg.str());
    }

    for (int i=0; i<sizeOfArray; ++i) {
        input_file >> data[i];
    }

    input_file.close();

    RemoveFile(rFileName);

    EMPIRE_API_LOG(2) << "Finished receiving array" << std::endl;

    if (PrintTiming) {
        EMPIRE_API_LOG(0) << "Receiving Array \"" << rFileName << "\" took: " << ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }
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

static bool StringEndsWith(const std::string& rString, const std::string& rSuffix)
{
    return (rString.size() >= rSuffix.size() && rString.compare(rString.size() - rSuffix.size(), rSuffix.size(), rSuffix) == 0);
}

static void ReadNumberAfterKeyword(const std::string& rKeyWord, const std::string& rCurrentLine, int& rNumber)
{
    std::string current_line_substr = rCurrentLine.substr(rCurrentLine.find(rKeyWord) + rKeyWord.length());
    std::istringstream line_stream(current_line_substr);
    line_stream >> rNumber;
}

template <typename TDataType>
static void AllocateMemory(TDataType** ppContainer, const std::size_t Size)
{
    *ppContainer = new TDataType[Size];
}

template <typename TDataType>
static void AllocateMemory(std::vector<TDataType>* pContainer, const std::size_t Size)
{
    if (pContainer->size() != Size) {
        pContainer->resize(Size);
    }
}

} // namespace EMPIRE_API_helpers

/***********************************************************************************************
 * \brief Establishes the necessary connection with the Emperor
 ***********/
static void EMPIRE_API_Connect(const char* inputFileName)
{
    const std::string file_name(inputFileName);

    if (EMPIRE_API_helpers::StringEndsWith(file_name, ".xml")) {
        EMPIRE_API_LOG(0) << "Called \"EMPIRE_API_Connect\" with an xml-file, which is no longer supported.\nPlease pass a config file instead" << std::endl;
    } else {
        if (EMPIRE_API_helpers::FileExists(file_name)) {
            // reading config-file
            std::ifstream input_file(file_name);
            EMPIRE_API_helpers::CheckStream(input_file, file_name);

            // reading file
            std::string current_line;

            while (std::getline(input_file, current_line)) {
                if (current_line.find("EchoLevel") != std::string::npos) {
                    EMPIRE_API_helpers::ReadNumberAfterKeyword("EchoLevel", current_line, EMPIRE_API_helpers::EchoLevel);
                }
                if (current_line.find("PrintTiming") != std::string::npos) {
                    EMPIRE_API_helpers::ReadNumberAfterKeyword("PrintTiming", current_line, EMPIRE_API_helpers::PrintTiming);
                }
            }
        } else {
            EMPIRE_API_LOG(0) << "Config-file \"" << file_name << "\" not found!" << std::endl;
        }
    }

    EMPIRE_API_LOG(0) << "Configuration:\n"
                      << "    EchoLevel: " << EMPIRE_API_helpers::EchoLevel << "\n"
                      << "    PrintTiming: " << EMPIRE_API_helpers::PrintTiming << std::endl;
}

/***********************************************************************************************
 * \brief Get user defined text by the element name in the XML input file
 * \param[in] elementName name of the XML element
 * \return user defined text
 ***********/
static char *EMPIRE_API_getUserDefinedText(char *elementName)
{
    EMPIRE_API_LOG(0) << "Called \"EMPIRE_API_getUserDefinedText\" with \"" << elementName << "\" which is no longer supported and can be removed" << std::endl;
    return const_cast<char*>("");
}

/***********************************************************************************************
 * \brief Send the mesh to the server
 * \param[in] name name of the mesh
 * \param[in] numNodes number of nodes
 * \param[in] numElems number of elements
 * \param[in] nodes coordinates of all nodes
 * \param[in] nodeIDs IDs of all nodes
 * \param[in] numNodesPerElem number of nodes per element
 * \param[in] elems connectivity table of all elements
 ***********/
static void EMPIRE_API_sendMesh(const char *name, const int numNodes, const int numElems, const double *nodes, const int *nodeIDs, const int *numNodesPerElem, const int *elems)
{
    EMPIRE_API_LOG(2) << "Attempting to send mesh \"" << std::string(name) << "\" with " << numNodes << " Nodes | " << numElems << " Elements ..." << std::endl;

    const auto start_time(std::chrono::steady_clock::now());

    const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_mesh_" + std::string(name) + ".vtk"));

    std::ofstream output_file;
    output_file.open(EMPIRE_API_helpers::GetTempFileName(file_name));
    EMPIRE_API_helpers::CheckStream(output_file, file_name);

    output_file << std::scientific << std::setprecision(7); // TODO maybe this should be configurable

    // write file header
    output_file << "# vtk DataFile Version 4.0\n";
    output_file << "vtk output\n";
    output_file << "ASCII\n";
    output_file << "DATASET UNSTRUCTURED_GRID\n\n";

    // write nodes
    int vtk_id = 0;
    std::unordered_map<int, int> node_vtk_id_map;
    output_file << "POINTS " << numNodes << " float\n";
    for (int i=0; i<numNodes; ++i) {
        output_file << nodes[i*3] << " " << nodes[i*3+1] << " " << nodes[i*3+2] << "\n";
        node_vtk_id_map[nodeIDs[i]] = vtk_id++;
    }
    output_file << "\n";

    // write cells connectivity
    int cell_list_size = 0;
    for (int i=0; i<numElems; ++i) {
        cell_list_size += numNodesPerElem[i] + 1;
    }

    int counter=0;
    output_file << "CELLS " << numElems << " " << cell_list_size << "\n";
    for (int i=0; i<numElems; ++i) {
        const int num_nodes_elem = numNodesPerElem[i];
        output_file << num_nodes_elem << " ";
        for (int j=0; j<num_nodes_elem; ++j) {
            output_file << node_vtk_id_map.at(elems[counter++]);
            if (j<num_nodes_elem-1) output_file << " "; // not adding a whitespace after last number
        }
        output_file << "\n";
    }

    output_file << "\n";

    // write cell types
    output_file << "CELL_TYPES " << numElems << "\n";
    for (int i=0; i<numElems; ++i) {
        output_file << EMPIRE_API_helpers::GetVtkCellType(numNodesPerElem[i]) << "\n";
    }

    output_file.close();
    EMPIRE_API_helpers::MakeFileVisible(file_name);

    EMPIRE_API_LOG(2) << "Finished sending mesh" << std::endl;

    if (EMPIRE_API_helpers::PrintTiming) {
        EMPIRE_API_LOG(0) << "Sending Mesh \"" << file_name << "\" took: " << EMPIRE_API_helpers::ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }
}

/***********************************************************************************************
 * \brief Recieve mesh from the server
 * \param[in] name name of the mesh
 * \param[in] numNodes number of nodes
 * \param[in] numElems number of elements
 * \param[in] nodes coordinates of all nodes
 * \param[in] nodeIDs IDs of all nodes
 * \param[in] numNodesPerElem number of nodes per element
 * \param[in] elems connectivity table of all elements
 ***********/
template <typename TDouble, typename TInt>
static void EMPIRE_API_recvMesh(const char *name, int *numNodes, int *numElems, TDouble* nodes, TInt* nodeIDs, TInt* numNodesPerElem, TInt* elem)
{
    EMPIRE_API_LOG(2) << "Attempting to receive mesh \"" << std::string(name) << "\" ..." << std::endl;

    const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_mesh_" + std::string(name) + ".vtk"));

    EMPIRE_API_helpers::WaitForFile(file_name);

    const auto start_time(std::chrono::steady_clock::now());

    std::ifstream input_file(file_name);
    EMPIRE_API_helpers::CheckStream(input_file, file_name);

    // reading file
    std::string current_line;
    bool nodes_read = false;

    while (std::getline(input_file, current_line)) {
        // reading nodes
        if (current_line.find("POINTS") != std::string::npos) {
            if (nodes_read) throw std::runtime_error("The nodes were read already!");
            nodes_read = true;

            EMPIRE_API_helpers::ReadNumberAfterKeyword("POINTS", current_line, *numNodes);

            EMPIRE_API_LOG(2) << "Mesh contains " << *numNodes << " Nodes" << std::endl;

            // allocating memory for nodes
            // note that this has to be deleted by the client!
            EMPIRE_API_helpers::AllocateMemory(nodes, (*numNodes) * 3); // *nodes = new double[(*numNodes) * 3];
            EMPIRE_API_helpers::AllocateMemory(nodeIDs, *numNodes); // *nodeIDs = new int[*numNodes];

            for (int i=0; i<*numNodes*3; ++i) {
                input_file >> (*nodes)[i];
            }

            for (int i=0; i<*numNodes; ++i) {
                (*nodeIDs)[i] = i+1; // Node Ids have an offset of 1 from Kratos to VTK
            }
        }

        // reading elements
        if (current_line.find("CELLS") != std::string::npos) {
            if (!nodes_read) throw std::runtime_error("The nodes were not yet read!");

            int num_nodes_per_elem, node_id, cell_list_size;
            current_line = current_line.substr(current_line.find("CELLS") + 6); // removing "CELLS"
            std::istringstream line_stream(current_line);
            line_stream >> *numElems;
            line_stream >> cell_list_size;

            EMPIRE_API_LOG(2) << "Mesh contains " << *numElems << " Elements" << std::endl;

            // allocating memory for elements
            // note that this has to be deleted by the client!
            EMPIRE_API_helpers::AllocateMemory(numNodesPerElem, *numElems); // *numNodesPerElem = new int[*numElems];
            EMPIRE_API_helpers::AllocateMemory(elem, cell_list_size-*numElems); // *elem = new int[cell_list_size-*numElems];

            int counter=0;
            for (int i=0; i<*numElems; ++i) {
                input_file >> num_nodes_per_elem;
                (*numNodesPerElem)[i] = num_nodes_per_elem;
                for (int j=0; j<num_nodes_per_elem; ++j) {
                    input_file >> node_id;
                    (*elem)[counter++] = node_id+1; // Node Ids have an offset of 1 from Kratos to VTK
                }
            }
            break; // no further information reading required => CELL_TYPES are not used here
        }
    }

    input_file.close();

    EMPIRE_API_helpers::RemoveFile(file_name);

    EMPIRE_API_LOG(2) << "Finished receiving mesh" << std::endl;

    if (EMPIRE_API_helpers::PrintTiming) {
        EMPIRE_API_LOG(0) << "Receiving Mesh \"" << file_name << "\" took: " << EMPIRE_API_helpers::ElapsedSeconds(start_time) << " [sec]" << std::endl;
    }
}

/***********************************************************************************************
 * \brief Send data field to the server
 * \param[in] name name of the field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[in] dataField the data field to be sent
 ***********/
static void EMPIRE_API_sendDataField(const char *name, const int sizeOfArray, const double *dataField)
{
    const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_datafield_" + std::string(name) + ".dat"));

    EMPIRE_API_helpers::SendArray(file_name, sizeOfArray, dataField);
}

/***********************************************************************************************
 * \brief Receive data field from the server
 * \param[in] name name of the field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[out] dataField the data field to be received
 ***********/
static void EMPIRE_API_recvDataField(const char *name, const int sizeOfArray, double *dataField)
{
    const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_datafield_" + std::string(name) + ".dat"));

    EMPIRE_API_helpers::ReceiveArray(file_name, sizeOfArray, dataField);
}

/***********************************************************************************************
 * \brief Send signal to the server
 * \param[in] name name of the signal
 * \param[in] sizeOfArray size of the array (signal)
 * \param[in] signal the signal
 ***********/
static void EMPIRE_API_sendSignal_double(const char *name, const int sizeOfArray, const double *signal)
{
    const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_signal_" + std::string(name) + ".dat"));

    EMPIRE_API_helpers::SendArray(file_name, sizeOfArray, signal);
}

/***********************************************************************************************
 * \brief Receive signal from the server
 * \param[in] name name of the signal
 * \param[in] sizeOfArray size of the array (signal)
 * \param[in] signal the signal
 ***********/
static void EMPIRE_API_recvSignal_double(const char *name, const int sizeOfArray, double *signal)
{
    const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_signal_" + std::string(name) + ".dat"));

    EMPIRE_API_helpers::ReceiveArray(file_name, sizeOfArray, signal);
}

/***********************************************************************************************
 * \brief Receive the convergence signal of an loop
 * \return 1 means convergence, 0 means non-convergence
 ***********/
static int EMPIRE_API_recvConvergenceSignal(const std::string& rFileNameExtension="default")
{
    const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_convergence_signal_" + std::string(rFileNameExtension) + ".dat"));

    EMPIRE_API_helpers::WaitForFile(file_name);

    std::ifstream input_file(file_name);
    EMPIRE_API_helpers::CheckStream(input_file, file_name);

    int signal;
    input_file >> signal;

    if (!(signal == 0 || signal == 1)) {
        std::stringstream err_msg;
        err_msg << "Read an invalid convergence signal: " << signal << ", can only be 0 for non-convergence or 1 for convergence";
        throw std::runtime_error(err_msg.str());
    }

    input_file.close();

    EMPIRE_API_helpers::RemoveFile(file_name);

    return signal;
}

/***********************************************************************************************
 * \brief Send the convergence signal of an loop
 * \param[in] signal 1 means convergence, 0 means non-convergence
 ***********/
static void EMPIRE_API_sendConvergenceSignal(const int signal, const std::string& rFileNameExtension="default")
{
    const std::string file_name(EMPIRE_API_helpers::GetFullPath("EMPIRE_convergence_signal_" + std::string(rFileNameExtension) + ".dat"));

    if (!(signal == 0 || signal == 1)) {
        std::stringstream err_msg;
        err_msg << "Input can only be 0 for non-convergence or 1 for convergence, called with: " << signal;
        throw std::runtime_error(err_msg.str());
    }

    std::ofstream output_file;
    output_file.open(EMPIRE_API_helpers::GetTempFileName(file_name));
    EMPIRE_API_helpers::CheckStream(output_file, file_name);

    output_file << signal;

    output_file.close();
    EMPIRE_API_helpers::MakeFileVisible(file_name);
}

/***********************************************************************************************
 * \brief Performs disconnection and finalization operations to the Emperor
 ***********/
static void EMPIRE_API_Disconnect()
{
    EMPIRE_API_LOG(0) << "Called \"EMPIRE_API_Disconnect\" which is no longer necessary and can be removed" << std::endl;
}

// defining aliases to the functions to avoid the "unused-functions" compiler-warning if not all functions are used
const auto _alias_EMPIRE_API_Connect = EMPIRE_API_Connect;
const auto _alias_EMPIRE_API_Disconnect = EMPIRE_API_Disconnect;
const auto _alias_EMPIRE_API_getUserDefinedText = EMPIRE_API_getUserDefinedText;
const auto _alias_EMPIRE_API_sendMesh = EMPIRE_API_sendMesh;
const auto _alias_EMPIRE_API_sendDataField = EMPIRE_API_sendDataField;
const auto _alias_EMPIRE_API_recvDataField = EMPIRE_API_recvDataField;
const auto _alias_EMPIRE_API_sendSignal_double = EMPIRE_API_sendSignal_double;
const auto _alias_EMPIRE_API_recvSignal_double = EMPIRE_API_recvSignal_double;
const auto _alias_EMPIRE_API_EMPIRE_API_sendConvergenceSignal = EMPIRE_API_sendConvergenceSignal;
const auto _alias_EMPIRE_API_EMPIRE_API_recvConvergenceSignal = EMPIRE_API_recvConvergenceSignal;

#endif /* KRATOS_CO_SIM_EMPIRE_API_H_INCLUDED */
