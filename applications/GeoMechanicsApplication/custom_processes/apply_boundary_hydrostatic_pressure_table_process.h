// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "includes/table.h"

#include "custom_processes/apply_constant_boundary_hydrostatic_pressure_process.hpp"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyBoundaryHydrostaticPressureTableProcess
    : public ApplyConstantBoundaryHydrostaticPressureProcess
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyBoundaryHydrostaticPressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;

    ApplyBoundaryHydrostaticPressureTableProcess(ModelPart& model_part, Parameters rParameters);
    ApplyBoundaryHydrostaticPressureTableProcess(const ApplyBoundaryHydrostaticPressureTableProcess&) = delete;
    ApplyBoundaryHydrostaticPressureTableProcess& operator=(const ApplyBoundaryHydrostaticPressureTableProcess&) = delete;
    ~ApplyBoundaryHydrostaticPressureTableProcess() override = default;

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