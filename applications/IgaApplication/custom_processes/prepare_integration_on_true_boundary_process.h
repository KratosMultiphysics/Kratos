//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli
//

#pragma once

// Project includes
#include "containers/model.h"
#include "processes/process.h"
#include "spatial_containers/bins_dynamic.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) PrepareIntegrationOnTrueBoundaryProcess
    : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PrepareIntegrationOnTrueBoundaryProcess);

    using IndexType = std::size_t;
    using SizeType = std::size_t;

    using PointType = Node;
    using PointTypePointer = Node::Pointer;
    using PointVector = std::vector<PointTypePointer>;
    using PointIterator = PointVector::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = DistanceVector::iterator;
    using DynamicBins = BinsDynamic<3, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>;

    PrepareIntegrationOnTrueBoundaryProcess(
        Model& rModel,
        Parameters ThisParameters);

    ~PrepareIntegrationOnTrueBoundaryProcess() override = default;

    void Execute() override;

    const Parameters GetDefaultParameters() const override;

    std::string Info() const override
    {
        return "PrepareIntegrationOnTrueBoundaryProcess";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "PrepareIntegrationOnTrueBoundaryProcess";
    }

    void PrintData(std::ostream& rOStream) const override
    {
    }

private:
    Model& mrModel;
    Parameters mParameters;

    double ComputeSearchRadius(const ModelPart& rAnalysisModelPart) const;

    void ResetIntegrationData(ModelPart& rSurrogateModelPart) const;

    void AppendIntegrationData(
        Condition& rCondition,
        const array_1d<double, 3>& rIntegrationPoint,
        const double IntegrationWeight,
        const array_1d<double, 3>& rNormal) const;
};

inline std::istream& operator >> (
    std::istream& rIStream,
    PrepareIntegrationOnTrueBoundaryProcess& rThis);

inline std::ostream& operator << (
    std::ostream& rOStream,
    const PrepareIntegrationOnTrueBoundaryProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos
