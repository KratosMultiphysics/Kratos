// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#pragma once

// System includes
#include <cmath>
#include <iostream>

// Project includes
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyK0ProcedureProcess : public Process
{
  public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyK0ProcedureProcess);

    ApplyK0ProcedureProcess(ModelPart&  model_part,
                           const Parameters& );

    ~ApplyK0ProcedureProcess() override = default;

    void ExecuteFinalizeSolutionStep() override;

    std::string Info() const override;

  private:
      ModelPart& mrModelPart;

      void CalculateK0Stresses(Element& rElement);

};

}