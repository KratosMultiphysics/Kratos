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

#include "custom_processes/apply_constant_phreatic_surface_pressure_process.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
class ModelPart;

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyPhreaticSurfacePressureTableProcess : public ApplyConstantPhreaticSurfacePressureProcess
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyPhreaticSurfacePressureTableProcess);

    /// Defining a table with double argument and result type as table type.
    using TableType = Table<double, double>;

    ApplyPhreaticSurfacePressureTableProcess(ModelPart& model_part, Parameters rParameters);
    ApplyPhreaticSurfacePressureTableProcess(const ApplyPhreaticSurfacePressureTableProcess&) = delete;
    ApplyPhreaticSurfacePressureTableProcess& operator=(const ApplyPhreaticSurfacePressureTableProcess&) = delete;
    ~ApplyPhreaticSurfacePressureTableProcess() override = default;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;

    /// Turn back information as a string.
    std::string Info() const override;

private:
    TableType::Pointer mpTable;
    double             mTimeUnitConverter;
};

} // namespace Kratos