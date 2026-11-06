/*   ______     _____ _           ________
    / ____/___ / ___/(_)___ ___  /  _/ __ |
   / /   / __ \\__ \/ / __ `__ \ / // / / /
  / /___/ /_/ /__/ / / / / / / // // /_/ /
  \____/\____/____/_/_/ /_/ /_/___/\____/
  Kratos CoSimulationApplication

  License:         BSD License, see license.txt

  Main authors:    Philipp Bucher (https://github.com/philbucher)
*/

#ifndef CO_SIM_IO_C_INCLUDED
#define CO_SIM_IO_C_INCLUDED

/* C Interface for CoSimulation
   see "co_sim_io.hpp"
*/

#include <stdio.h>

#include "co_sim_io_c_info.h"
#include "co_sim_io_c_model_part.h"

#ifdef __cplusplus
extern "C" { /* Define extern C if C++ compiler is used */
#endif

enum CoSimIO_ConnectionStatus
{
    CoSimIO_NotConnected,
    CoSimIO_Connected,
    CoSimIO_Disconnected,
    CoSimIO_ConnectionError,
    CoSimIO_DisconnectionError
};

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_Hello(
    void);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_Connect(
    const CoSimIO_Info I_Settings);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_Disconnect(
    const CoSimIO_Info I_Info);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_ImportData(
    const CoSimIO_Info I_Info,
    int* O_Size,
    double** O_Data);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_ExportData(
    const CoSimIO_Info I_Info,
    const int I_Size,
    const double* I_Data);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_ImportMesh(
    const CoSimIO_Info I_Info,
    CoSimIO_ModelPart O_ModelPart);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_ExportMesh(
    const CoSimIO_Info I_Info,
    const CoSimIO_ModelPart I_ModelPart);


void CoSimIO_PrintInfo(FILE *Stream,
    const CoSimIO_Info I_Info);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_ImportInfo(
    const CoSimIO_Info I_Info);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_ExportInfo(
    const CoSimIO_Info I_Info);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_Register(
    const CoSimIO_Info I_Info,
    CoSimIO_Info (*I_FunctionPointer)(const CoSimIO_Info I_Info));

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_Run(const CoSimIO_Info I_Info);

void* CoSimIO_Malloc(size_t size);

void CoSimIO_Free (void* ptr);

#ifdef __cplusplus
}
#endif

#endif /* CO_SIM_IO_C_INCLUDED */
