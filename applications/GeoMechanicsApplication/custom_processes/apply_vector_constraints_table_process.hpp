// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra,
//                   Anne van de Graaf
//
#pragma once

#include "processes/process.h"

namespace Kratos
{

class ModelPart;
class Parameters;


class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyVectorConstraintsTableProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyVectorConstraintsTableProcess);

    ApplyVectorConstraintsTableProcess(ModelPart&        rModelPart,
                                       const Parameters& rSettings);

    ~ApplyVectorConstraintsTableProcess() override;

    ApplyVectorConstraintsTableProcess(const ApplyVectorConstraintsTableProcess&) = delete;
    ApplyVectorConstraintsTableProcess& operator=(const ApplyVectorConstraintsTableProcess&) = delete;

    using ProcessUniquePointer = std::unique_ptr<Process>;

    void ExecuteInitialize() override;
    void ExecuteInitializeSolutionStep() override;

private:
    static std::vector<Parameters> CreateParametersForActiveComponents(const Parameters& rSettings);
    static std::vector<char> ActiveComponents(const Parameters& rSettings);
    static Parameters CreateParametersForComponent(const Parameters& rSettings, char component);
    static std::size_t ComponentToIndex(char component);
    ProcessUniquePointer MakeProcessFor(const Parameters& rParameters) const;

    ModelPart& mrModelPart;
    std::vector<ProcessUniquePointer> mProcesses;
};

}
