// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Lorenzo Gracia,
//                   Aron Noordam,
//                   Vahid Galavi,
//                   Marjan Fathian
//
#pragma once

#include "processes/process.h"

namespace Kratos
{


class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyExcavationProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyExcavationProcess);

    ApplyExcavationProcess(ModelPart&        rModelPart,
                           const Parameters& rSettings);

    ~ApplyExcavationProcess() override;

    ApplyExcavationProcess(const ApplyExcavationProcess&) = delete;
    ApplyExcavationProcess& operator=(const ApplyExcavationProcess&) = delete;

    void ExecuteInitialize() override;

private:
    ModelPart& mrModelPart;
    bool mDeactivateSoilPart;
};

}