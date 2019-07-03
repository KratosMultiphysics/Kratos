// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Philipp Bucher
//

#ifndef KRATOS_CO_SIM_EMPIRE_API_H_INCLUDED
#define KRATOS_CO_SIM_EMPIRE_API_H_INCLUDED

/*
This file has the same API as EMPIRE (defined in "EMPIRE_API.h"), hence it can be included instead of EMPIRE
It used FileIO for data-exchange, in VTK-format
Note:
- This file cannot have Kratos-includes, because it is also included in other codes!
- This file is intended to be header-only, such that other codes do not have to link against a library
- Requires c++11 to compile

*/

// System includes
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stdio.h>
#include <chrono>
#include <thread>

//#define VTK_USE_BINARY // comment this for using binary for the files

namespace CoSimEMPIRE_API {

namespace helpers {

const std::string ConvergenceSignalFileName = "EMPIRE_convergence_signal.dat";
const std::string TempFilePreString = ".";
const int EchoLevel = 0;

bool FileExists(const std::string& rFileName)
{
    std::ifstream infile(rFileName);
    return infile.good(); // no need to close manually
}

std::string GetTempFileName(const std::string& rFileName)
{
    // wrapped in a function such that it could be changed easily (e.g. if files are in a folder)
    return TempFilePreString + rFileName;
}

void WaitForFile(const std::string& rFileName)
{
    if (EchoLevel>0) std::cout << "Waiting for file: \"" << rFileName << "\"" << std::endl;
    while(!FileExists(rFileName)) {
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        if (EchoLevel>1) std::cout << "Waiting" << std::endl;
    }
    if (EchoLevel>0) std::cout << "Found file: \"" << rFileName << "\"" << std::endl;
}

void RemoveFile(const std::string& rFileName)
{
    if (std::remove(rFileName.c_str()) != 0) {
        std::cout << "Warning: \"" << rFileName << "\" could not be deleted!" << std::endl;
    }
}

void MakeFileVisible(const std::string& rFinalFileName)
{
    if (std::rename(GetTempFileName(rFinalFileName).c_str(), rFinalFileName.c_str()) != 0) {
        std::cout << "Warning: \"" << rFinalFileName << "\" could not be made visible!" << std::endl;
    }
}

template <typename T>
void CheckStream(const T& rStream, const std::string& rFileName)
{
    if (!rStream.is_open()) {
        std::stringstream err_msg;
        err_msg << rFileName << " could not be opened!";
        throw std::runtime_error(err_msg.str());
    }
}

} // namespace helpers

/***********************************************************************************************
 * \brief Establishes the necessary connection with the Emperor
 ***********/
void EMPIRE_API_Connect(char* inputFileName)
{
    std::cout << "Called \"EMPIRE_API_Connect\" which is no longer necessary and can be removed" << std::endl;
}

/***********************************************************************************************
 * \brief Get user defined text by the element name in the XML input file
 * \param[in] elementName name of the XML element
 * \return user defined text
 ***********/
char *EMPIRE_API_getUserDefinedText(char *elementName)
{
    std::cout << "Called \"EMPIRE_API_getUserDefinedText\" with \"" << elementName << "\" which is no longer working and can be removed" << std::endl;
    return ""; // TODO this gives a warning, find better solution
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
void EMPIRE_API_sendMesh(char *name, int numNodes, int numElems, double *nodes, int *nodeIDs, int *numNodesPerElem, int *elems)
{

    // rename file after writing such that it becomes visible
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
void EMPIRE_API_recvMesh(char *name, int *numNodes, int *numElems, double **nodes, int **nodeIDs, int **numNodesPerElem, int **elem)
{
    // wait for file

    // delete file after reading? // TODO

}

/***********************************************************************************************
 * \brief Send data field to the server
 * \param[in] name name of the field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[in] dataField the data field to be sent
 ***********/
void EMPIRE_API_sendDataField(char *name, int sizeOfArray, double *dataField)
{

    // rename file after writing such that it becomes visible
}

/***********************************************************************************************
 * \brief Receive data field from the server
 * \param[in] name name of the field
 * \param[in] sizeOfArray size of the array (data field)
 * \param[out] dataField the data field to be received
 ***********/
void EMPIRE_API_recvDataField(char *name, int sizeOfArray, double *dataField)
{
    // wait for file

    // delete file after reading? // TODO

}

/***********************************************************************************************
 * \brief Send signal to the server
 * \param[in] name name of the signal
 * \param[in] sizeOfArray size of the array (signal)
 * \param[in] signal the signal
 ***********/
void EMPIRE_API_sendSignal_double(char *name, int sizeOfArray, double *signal)
{
    const std::string file_name("EMPIRE_signal_" + std::string(name));

    std::ofstream output_file;
    output_file.open(helpers::GetTempFileName(file_name));
    helpers::CheckStream(output_file, file_name);

    // TODO write size in first line?

    for (int i=0; i<sizeOfArray-1; ++i) {
        output_file << signal[i] << " ";
    }
    output_file << signal[sizeOfArray-1]; // outside to not have trailing whitespace

    output_file.close();

    helpers::MakeFileVisible(file_name);
}

/***********************************************************************************************
 * \brief Receive signal from the server
 * \param[in] name name of the signal
 * \param[in] sizeOfArray size of the array (signal)
 * \param[in] signal the signal
 ***********/
void EMPIRE_API_recvSignal_double(char *name, int sizeOfArray, double *signal)
{
    // wait for file

    // delete file after reading? // TODO

}

/***********************************************************************************************
 * \brief Receive the convergence signal of an loop
 * \return 1 means convergence, 0 means non-convergence
 ***********/
int EMPIRE_API_recvConvergenceSignal()
{
    helpers::WaitForFile(helpers::ConvergenceSignalFileName);

    std::ifstream input_file(helpers::ConvergenceSignalFileName);
    helpers::CheckStream(input_file, helpers::ConvergenceSignalFileName);

    int signal;
    input_file >> signal;

    if (!(signal == 0 || signal == 1)) {
        std::stringstream err_msg;
        err_msg << "Read an invalid convergence signal: " << signal << ", can only be 0 for non-convergence or 1 for convergence";
        throw std::runtime_error(err_msg.str());
    }

    helpers::RemoveFile(helpers::ConvergenceSignalFileName);

    return signal;
}

/***********************************************************************************************
 * \brief Send the convergence signal of an loop
 * \param[in] signal 1 means convergence, 0 means non-convergence
 ***********/
void EMPIRE_API_sendConvergenceSignal(int signal)
{
    if (!(signal == 0 || signal == 1)) {
        std::stringstream err_msg;
        err_msg << "Input can only be 0 for non-convergence or 1 for convergence, called with: " << signal;
        throw std::runtime_error(err_msg.str());
    }

    std::ofstream output_file;
    output_file.open(helpers::GetTempFileName(helpers::ConvergenceSignalFileName));
    helpers::CheckStream(output_file, helpers::ConvergenceSignalFileName);

    output_file << signal;
    output_file.close();

    helpers::MakeFileVisible(helpers::ConvergenceSignalFileName);
}

/***********************************************************************************************
 * \brief Performs disconnection and finalization operations to the Emperor
 ***********/
void EMPIRE_API_Disconnect()
{
    std::cout << "Called \"EMPIRE_API_Disconnect\" which is no longer necessary and can be removed" << std::endl;
}

} // namespace CoSimEMPIRE_API

#endif /* KRATOS_CO_SIM_EMPIRE_API_H_INCLUDED */
