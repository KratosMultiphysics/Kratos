//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#pragma once

// System includes
#include <array>
#include <vector>

// Project includes
#include "includes/model_part.h"
#include "modeler/modeler.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) LocalRefinementModeler : public Modeler
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(LocalRefinementModeler);

    using SizeType = std::size_t;
    using IndexType = std::size_t;
    using RectangleType = std::array<double, 4>; // [umin, umax, vmin, vmax]

    LocalRefinementModeler() = default;

    LocalRefinementModeler(
        Model& rModel,
        const Parameters ModelParameters = Parameters());

    ~LocalRefinementModeler() override = default;

    Modeler::Pointer Create(
        Model& rModel,
        const Parameters ModelParameters) const override;

    void SetupModelPart() override;

    const Parameters GetDefaultParameters() const override;

private:
    using RectType = RectangleType;

    void GenerateRefinementRegions();

    void RunGapSbmPatchModelers();

    void ProcessSkinModelPart(
        const std::string& rSkinModelPartName,
        const std::string& rUpdatedSubModelPartName,
        const RectType& rRect,
        const IndexType RegionIndex,
        const bool ReverseBorderOrientation) const;

    void RunGapSbmPatch(
        const RectType& rRect,
        const std::string& rPatchSuffix,
        const std::string& rPatchSkinModelPartName,
        const std::string& rInnerInitialSkinModelPartName,
        const std::string& rOuterInitialSkinModelPartName,
        const Parameters* pRegionParameters) const;

    std::string BuildPatchSkinModelPartName(const std::string& rPatchSuffix) const;

    ModelPart& CreateOrResetSubModelPart(
        ModelPart& rParentModelPart,
        const std::string& rName) const;

    Model* mpModel = nullptr;
    SizeType mEchoLevel = 0;
    RectType mBaseRect{0.0, 1.0, 0.0, 1.0};
    std::vector<RectType> mRefinementRegions;
    std::vector<IndexType> mRefinementRegionSourceIndices;
};

} // namespace Kratos
