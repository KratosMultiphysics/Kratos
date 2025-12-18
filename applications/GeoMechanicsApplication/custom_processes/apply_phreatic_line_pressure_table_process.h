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

#include "custom_processes/apply_constant_phreatic_line_pressure_process.h"

namespace Kratos
{
class ModelPart;
class Node;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyPhreaticLinePressureTableProcess : public ApplyConstantPhreaticLinePressureProcess
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPhreaticLinePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;

    ApplyPhreaticLinePressureTableProcess(ModelPart& model_part, Parameters rParameters);
    ApplyPhreaticLinePressureTableProcess(const ApplyPhreaticLinePressureTableProcess&) = delete;
    ApplyPhreaticLinePressureTableProcess& operator=(const ApplyPhreaticLinePressureTableProcess&) = delete;
    ~ApplyPhreaticLinePressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

    /// Turn back information as a string.
    std::string Info() const override;

protected:
    double CalculatePressurewithTable(const Node& rNode, const array_1d<double, 2>& y) const;

private:
    /// Member Variables
    array_1d<TableType::Pointer, 2> mpTable;
    double                          mTimeUnitConverter;
};

} // namespace Kratos