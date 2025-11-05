//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolò Antonelli

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "modeler/modeler.h"

// Application includes
#include "custom_modelers/nurbs_geometry_modeler_sbm.h"
#include "custom_modelers/iga_modeler_sbm.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) PatchSubdivisionModeler : public Modeler
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PatchSubdivisionModeler);

    using RectangleType = std::array<double, 4>; // [umin, umax, vmin, vmax]

    struct RefinementRegionData
    {
        RectangleType Rectangle;
        bool HasPolynomialOrder = false;
        std::vector<int> PolynomialOrder;
        bool HasNumberOfKnotSpans = false;
        std::vector<int> NumberOfKnotSpans;
    };

    PatchSubdivisionModeler() = default;

    PatchSubdivisionModeler(Model& rModel, const Parameters ModelParameters = Parameters());

    ~PatchSubdivisionModeler() override = default;

    Modeler::Pointer Create(Model& rModel, const Parameters ModelParameters) const override;

    void SetupModelPart() override;

    const Parameters GetDefaultParameters() const override;

    const std::vector<RectangleType>& GetSubdomains() const;

private:
    // Optional alias for one rectangular subdomain:
    // [U_MIN, U_MAX, V_MIN, V_MAX].
    using RectType = std::array<double, 4>;

    // Performs input validation, parameter normalization, logging, subdivision,
    // and writes global PARAMETER_SPACE_CORNERS. Returns the target ModelPart.
    ModelPart& InitializeSetup(
        Parameters& rParameters,
        Parameters& rGeometryBaseOut,
        Parameters& rAnalysisBaseOut,
        std::vector<std::string>& rAnalysisTargetsOut,
        std::unordered_map<std::string, std::vector<std::string>>& rCondNameToPartsOut);

    // Processes a single patch (subdomain). Builds/upgrades geometry and analysis,
    // sets per-patch values (corners, spans, orders), classifies BREP edges,
    // duplicates/filters condition entries as needed, and merges into the parent.
    void ProcessPatch(
        std::size_t iPatch,
        const RectType& rRect,
        const Parameters& rGeometryBase,
        const Parameters& rAnalysisBase,
        const std::vector<std::string>& rAnalysisTargetSubmodelParts,
        const std::unordered_map<std::string, std::vector<std::string>>& rConditionNameToModelParts,
        ModelPart& rParentModelPart,
        const std::string& rPrefix);
        
    void GenerateSubdivision();

    // Returns true if the skin model part (by name) is fully contained
    // within the rectangular patch in physical XY space derived from
    // the given rect (parameter-space) and the base geometry extents.
    bool IsSkinFullyInsidePatch(
        const std::string& rSkinModelPartName,
        const RectType& rRect,
        const Parameters& rGeometryBase) const;

    std::vector<RectangleType> mSubdomains;
    std::vector<RefinementRegionData> mRefinementRegions;
    Model* mpModel = nullptr;

    // (removed) BuildGlobalSubModelParts: unused in PatchSubdivisionModeler

    // Classify BREP edges of a patch as external/internal and append external entry to list
    void ClassifyAndAppendBrepEdges(
        const Vector& rPatchLowerUvw,
        const Vector& rPatchUpperUvw,
        const Vector& rBaseLowerUvw,
        const Vector& rBaseUpperUvw,
        ModelPart& rPatchModelPart,
        Parameters& rExtTemplate,
        Parameters& rIntTemplate,
        Parameters& rReducedList) const;

    // Build final element_condition_list for SBM patches based on body-fitted list and original list
    Parameters BuildSbmElementConditionList(
        const Parameters& rOriginalConditionList,
        const Parameters& rBodyFittedReduced) const;

};

} // namespace Kratos
