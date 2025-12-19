// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

#include "includes/table.h"

#include "custom_processes/apply_constant_boundary_phreatic_line_pressure_process.h"

namespace Kratos
{
class ModelPart;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyBoundaryPhreaticLinePressureTableProcess
    : public ApplyConstantBoundaryPhreaticLinePressureProcess
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyBoundaryPhreaticLinePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;

    ApplyBoundaryPhreaticLinePressureTableProcess(ModelPart& model_part, Parameters rParameters);
    ApplyBoundaryPhreaticLinePressureTableProcess(const ApplyBoundaryPhreaticLinePressureTableProcess&) = delete;
    ApplyBoundaryPhreaticLinePressureTableProcess& operator=(const ApplyBoundaryPhreaticLinePressureTableProcess&) = delete;
    ~ApplyBoundaryPhreaticLinePressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

    /// Turn back information as a string.
    std::string Info() const override;

private:
    /// Member Variables
    TableType::Pointer mpTable;
    double             mTimeUnitConverter;
};

} // namespace Kratos