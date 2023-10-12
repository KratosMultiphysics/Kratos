// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//
#include "apply_normal_load_table_process.h"

namespace Kratos
{

void ApplyNormalLoadTableProcess::Execute()
{
    Process::Execute();
}

void ApplyNormalLoadTableProcess::ExecuteInitialize()
{
    for (const auto& process : mProcesses) {
        process->ExecuteInitialize();
    }}

void ApplyNormalLoadTableProcess::ExecuteInitializeSolutionStep()
{
    for (const auto& process : mProcesses) {
        process->ExecuteInitializeSolutionStep();
    }}
}
