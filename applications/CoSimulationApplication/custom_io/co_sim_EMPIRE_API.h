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

*/

// System includes
#include <iostream>
#include <fstream>
#include <stdexcept>

//#define VTK_USE_BINARY // comment this for using binary for the files

#ifdef __cplusplus // TODO check if this is needed
extern "C" { // Define extern C if C++ compiler is used
#endif

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

    // rename file after writing such that it becomes visible
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
    // wait for file

    // delete file after reading? // TODO

}

/***********************************************************************************************
 * \brief Send the convergence signal of an loop
 * \param[in] signal 1 means convergence, 0 means non-convergence
 ***********/
void EMPIRE_API_sendConvergenceSignal(int signal)
{
    if (!(signal == 0 || signal == 1)) {
        std::stringstream err_msg;
        err_msg << "Input can only be 0 non-convergence or 1 convergence";
        err_msg << ", called with: " << signal;
        throw std::runtime_error(err_msg.str());
    }

    std::ofstream output_file;
    output_file.open(".EMPIRE_convergence_signal.dat");
    output_file << signal;
    output_file.close();

    // rename file after writing such that it becomes visible
    std::rename(".EMPIRE_convergence_signal.dat", "EMPIRE_convergence_signal.dat");
}

/***********************************************************************************************
 * \brief Performs disconnection and finalization operations to the Emperor
 ***********/
void EMPIRE_API_Disconnect()
{
    std::cout << "Called \"EMPIRE_API_Disconnect\" which is no longer necessary and can be removed" << std::endl;
}

#ifdef __cplusplus
}
#endif

#endif /* KRATOS_CO_SIM_EMPIRE_API_H_INCLUDED */
