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
#include "custom_modelers/nurbs_geometry_modeler_gap_sbm.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) MultipatchModeler : public Modeler
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(MultipatchModeler);

    using RectangleType = std::array<double, 4>; // [umin, umax, vmin, vmax]

    struct RefinementRegionData
    {
        RectangleType Rectangle;
        bool HasPolynomialOrder = false;
        std::vector<int> PolynomialOrder;
        bool HasNumberOfKnotSpans = false;
        std::vector<int> NumberOfKnotSpans;
    };

    MultipatchModeler() = default;

    MultipatchModeler(Model& rModel, const Parameters ModelParameters = Parameters());

    ~MultipatchModeler() override = default;

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
        Parameters& r_parameters,
        Parameters& rGeometryBaseOut,
        Parameters& rAnalysisBaseOut,
        std::vector<std::string>& rAnalysisTargetsOut,
        std::unordered_map<std::string, std::vector<std::string>>& rCondNameToPartsOut);

    // Processes a single patch (subdomain). Builds/upgrades geometry and analysis,
    // sets per-patch values (corners, spans, orders), classifies BREP edges,
    // duplicates/filters condition entries as needed, and merges into the parent.
    void ProcessPatch_(
        std::size_t i_patch,
        const RectType& rect,
        const Parameters& geometry_base,
        const Parameters& analysis_base,
        const std::vector<std::string>& analysis_target_submodel_parts,
        const std::unordered_map<std::string, std::vector<std::string>>& condition_name_to_model_parts,
        ModelPart& r_parent_model_part,
        const std::string& prefix);

    // Specialized: process only the base domain (SBM if skin is present)
    void ProcessBasePatch_(
        const RectType& rect,
        const Parameters& geometry_base,
        const Parameters& analysis_base,
        ModelPart& r_parent_model_part,
        const std::string& prefix);

    // Specialized: process refinement patch as body-fitted, only elements (FluidDomain)
    void ProcessRefPatch(
        const RectType& rect,
        const Parameters& geometry_base,
        const Parameters& analysis_base,
        ModelPart& r_parent_model_part,
        const std::string& prefix);
        
    void GenerateSubdivision();

    // Returns true if the skin model part (by name) is fully contained
    // within the rectangular patch in physical XY space derived from
    // the given rect (parameter-space) and the base geometry extents.
    bool IsSkinFullyInsidePatch(
        const std::string& rSkinModelPartName,
        const RectType& rect,
        const Parameters& geometry_base) const;

    // Utility: ensure a model part exists and is empty of nodes/conditions/elements/geometries
    ModelPart& CreateOrResetModelPart_(const std::string& rName) const;

    // Utility: append a rectangular closed loop (4 LineCondition2D2N) to a model part using parameter-space coords
    void AppendRectangleSkinLoop(ModelPart& rModelPart, const RectType& rect) const;

    // Builds a skin model part composed by the refinement regions as closed rectangular loops
    // Returns the created/cleared ModelPart reference
    ModelPart& CreateSkinCouplingModelPartForRefinements(const std::string& rSkinModelPartName) const;

    // Builds the inner-initial skin model part by copying the surrogate outer
    // loop of the refinement patch (Patch2) as a set of degree-1 NURBS curves,
    // one curve per surrogate condition.
    // Returns the created/cleared ModelPart reference
    ModelPart& CreateSkinInnerInitialFromRefinementSurrogateOuter(const std::string& rSkinModelPartName) const;

    std::vector<RectangleType> mSubdomains;
    std::vector<RefinementRegionData> mRefinementRegions;
    RectangleType mBaseRect{0.0, 1.0, 0.0, 1.0};
    Model* mpModel = nullptr;

    /// Build (or refresh) a global "ConvectionDiffusionDomain" submodelpart under the given root
    /// by collecting the elements from every submodelpart whose name starts with "Patch".
    /// Also ensures the nodes used by those elements are members of the global submodelpart.
    void BuildGlobalSubModelParts(
        ModelPart& r_parent_model_part,
        const std::vector<std::string>& rTargetNames) const;

    struct PatchPreparationResult
    {
        std::vector<ModelPart::IndexType> GeometryIdsBefore;
        Parameters PatchGeometryForSbmCoupling;
    };

    PatchPreparationResult PreparePatchGeometryAndData(
        ModelPart& r_patch_model_part,
        Parameters& r_patch_geometry,
        const Parameters& r_geometry_base,
        const RectType& rect,
        const bool use_sbm_for_this_patch,
        const std::string& patch_full_name,
        const bool is_base_patch) const;

    PatchPreparationResult PreparePatchGeometryAndDataGapSbm(
        ModelPart& r_patch_model_part,
        Parameters& r_patch_geometry,
        const Parameters& r_geometry_base,
        const RectType& rect,
        const bool use_sbm_for_this_patch,
        const std::string& patch_full_name,
        const bool is_base_patch) const;

};

} // namespace Kratos
