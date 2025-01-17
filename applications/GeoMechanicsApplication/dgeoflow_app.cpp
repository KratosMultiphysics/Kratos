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

#include <iostream>
#include <string>

#include "custom_workflows/dgeoflow.h"

using namespace Kratos;


int main(int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "Too few arguments: provide at least the working directory and a project parameters file" << std::endl;
    return 1;
  }

  auto workingDirectory = argv[1];
  auto projectFile      = argv[2];

  auto dummy_log_callback = [](const char*){};
  auto dummy_report_progress = [](double){};
  auto dummy_report_textual_progress = [](const char*){};
  auto never_cancel = [](){ return false; };

  KratosExecute execute;
  const int status = execute.ExecuteFlowAnalysis(workingDirectory, projectFile, 4.0, 6.0, 0.05,
                                                 "PorousDomain.0", dummy_log_callback, dummy_report_progress,dummy_report_textual_progress, never_cancel);

  return 0;
}