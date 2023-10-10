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
//                   Jonathan Nuttall

#pragma once

#include "includes/table.h"

#include "custom_processes/apply_constant_phreatic_multi_line_pressure_process.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ApplyPhreaticMultiLinePressureTableProcess : public ApplyConstantPhreaticMultiLinePressureProcess
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(ApplyPhreaticMultiLinePressureTableProcess);

    ApplyPhreaticMultiLinePressureTableProcess(ModelPart& model_part, Parameters rParameters);
    ~ApplyPhreaticMultiLinePressureTableProcess() override = default;
    ApplyPhreaticMultiLinePressureTableProcess& operator=(ApplyPhreaticMultiLinePressureTableProcess const&) = delete;
    ApplyPhreaticMultiLinePressureTableProcess(ApplyPhreaticMultiLinePressureTableProcess const&) = delete;

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override;
    std::string Info() const override;
    void PrintInfo(std::ostream& rOStream) const override;

private:
    using TableType = Table<double, double>;
    std::vector<TableType::Pointer> mpTable;
    double mTimeUnitConverter;
}; // Class ApplyPhreaticMultiLinePressureTableProcess

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyPhreaticMultiLinePressureTableProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyPhreaticMultiLinePressureTableProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


} // namespace Kratos.

