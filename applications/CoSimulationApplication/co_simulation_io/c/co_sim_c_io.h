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

#ifndef KRATOS_CO_SIM_C_IO_H_INCLUDED
#define KRATOS_CO_SIM_C_IO_H_INCLUDED


#ifdef __cplusplus
extern "C" { // Define extern C if C++ compiler is used
#endif

static void Connect(const char* pName);

static void Disconnect(const char* pName);

static void ImportData(const char* pName);
static void ExportData(const char* pName);

static void ImportMesh(const char* pName);
static void ExportMesh(const char* pName);

static void ImportGeometry(const char* pName);
static void ExportGeometry(const char* pName);

#ifdef __cplusplus
}
#endif

#endif // KRATOS_CO_SIM_C_IO_H_INCLUDED
