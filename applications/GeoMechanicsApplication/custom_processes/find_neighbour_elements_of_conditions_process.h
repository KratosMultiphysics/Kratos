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
//                   Richard Faasse
//

#pragma once

#include "processes/process.h"

namespace Kratos
{
class ModelPart;

class KRATOS_API(GEO_MECHANICS_APPLICATION) FindNeighbourElementsOfConditionsProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(FindNeighbourElementsOfConditionsProcess);

    explicit FindNeighbourElementsOfConditionsProcess(ModelPart& rModelPart);
    FindNeighbourElementsOfConditionsProcess& operator=(const FindNeighbourElementsOfConditionsProcess&) = delete;
    FindNeighbourElementsOfConditionsProcess(const FindNeighbourElementsOfConditionsProcess&) = delete;
    ~FindNeighbourElementsOfConditionsProcess() override = default;

    void                      Execute() override;
    [[nodiscard]] std::string Info() const override;
    void                      PrintData(std::ostream& rOStream) const override;

private:
    ModelPart& mrModelPart;

    void FindNeighbouringElementsForAllBoundaryTypes();

    [[nodiscard]] bool AllConditionsHaveAtLeastOneNeighbour() const;
    [[noreturn]] void  ReportConditionsWithoutNeighboursAndThrow() const;
};

std::ostream& operator<<(std::ostream& rOStream, const FindNeighbourElementsOfConditionsProcess& rThis);

} // namespace Kratos