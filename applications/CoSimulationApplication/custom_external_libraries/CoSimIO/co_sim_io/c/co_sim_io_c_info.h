/*   ______     _____ _           ________
    / ____/___ / ___/(_)___ ___  /  _/ __ |
   / /   / __ \\__ \/ / __ `__ \ / // / / /
  / /___/ /_/ /__/ / / / / / / // // /_/ /
  \____/\____/____/_/_/ /_/ /_/___/\____/
  Kratos CoSimulationApplication

  License:         BSD License, see license.txt

  Main authors:    Philipp Bucher (https://github.com/philbucher)
*/

#ifndef CO_SIM_IO_C_INFO_INCLUDED
#define CO_SIM_IO_C_INFO_INCLUDED


#include "define_c.h"

typedef struct CoSimIO_Info
{
    void* PtrCppInfo;
} CoSimIO_Info;


CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_CreateInfo(void);

CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_CopyInfo(const CoSimIO_Info I_Info);

int CoSimIO_FreeInfo(CoSimIO_Info I_Info);

int CoSimIO_Info_Has(const CoSimIO_Info I_Info, const char* I_Key);
void CoSimIO_Info_Erase(const CoSimIO_Info I_Info, const char* I_Key);
void CoSimIO_Info_Clear(const CoSimIO_Info I_Info);
int CoSimIO_Info_Size(const CoSimIO_Info I_Info);

int CoSimIO_Info_GetInt(const CoSimIO_Info I_Info, const char* I_Key);
double CoSimIO_Info_GetDouble(const CoSimIO_Info I_Info, const char* I_Key);
int CoSimIO_Info_GetBool(const CoSimIO_Info I_Info, const char* I_Key);
const char* CoSimIO_Info_GetString(const CoSimIO_Info I_Info, const char* I_Key);
CO_SIM_IO_NODISCARD CoSimIO_Info CoSimIO_Info_GetInfo(const CoSimIO_Info I_Info, const char* I_Key);

void CoSimIO_Info_SetInt(CoSimIO_Info I_Info, const char* I_Key, const int I_Value);
void CoSimIO_Info_SetDouble(CoSimIO_Info I_Info, const char* I_Key, const double I_Value);
void CoSimIO_Info_SetBool(CoSimIO_Info I_Info, const char* I_Key, const int I_Value);
void CoSimIO_Info_SetString(CoSimIO_Info I_Info, const char* I_Key, const char* I_Value);
void CoSimIO_Info_SetInfo(CoSimIO_Info I_Info, const char* I_Key, CoSimIO_Info I_Value);

#endif /* CO_SIM_IO_C_INFO_INCLUDED */
