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

#ifndef KRATOS_CO_SIM_C_IO_H_INCLUDED
#define KRATOS_CO_SIM_C_IO_H_INCLUDED

// C Interface for CoSimulation
// see "co_sim_io.h"

#ifdef __cplusplus
extern "C" { // Define extern C if C++ compiler is used
#endif

void CoSimIO_Connect(
    const char* pConnectionName,
    const char* pSettingsFileName);

void CoSimIO_Disconnect(
    const char* pConnectionName);

void CoSimIO_ImportData(
    const char* pConnectionName,
    const char* pIdentifier,
    int* pSize,
    double** ppData);

void CoSimIO_ExportData(
    const char* pConnectionName,
    const char* pIdentifier,
    int Size,
    double** pData);

void CoSimIO_ImportMesh(
    const char* pConnectionName,
    const char* pIdentifier,
    int* pNumberOfNodes,
    int* pNumberOfElements,
    double** ppNodalCoordinates,
    int** ppElementConnectivities,
    int** ppElementTypes);

void CoSimIO_ExportMesh(
    const char* pConnectionName,
    const char* pIdentifier,
    int NumberOfNodes,
    int NumberOfElements,
    double** pNodalCoordinates,
    int** pElementConnectivities,
    int** pElementTypes);

void CoSimIO_RegisterAdvanceInTime(
    const char* pConnectionName,
    double (*pFunctionPointer)(double));

void CoSimIO_RegisterSolvingFunction(
    const char* pConnectionName,
    const char* pFunctionName,
    void (*pFunctionPointer)());

void CoSimIO_RegisterDataExchangeFunction(
    const char* pConnectionName,
    const char* pFunctionName,
    void (*pFunctionPointer)(const char*, const char*));

void CoSimIO_Run(const char* pConnectionName);

int CoSimIO_IsConverged(const char* pConnectionName);

// The following functions are intended to only be used from Fortran
// This is because memory used in the CoSimIO should be allocated and freed in C
void _AllocateMemoryInt(const int Size, int** ppData);
void _AllocateMemoryDouble(const int Size, double** ppData);
void _FreeMemory(void** ppData);

#ifdef __cplusplus
}
#endif

#endif // KRATOS_CO_SIM_C_IO_H_INCLUDED
