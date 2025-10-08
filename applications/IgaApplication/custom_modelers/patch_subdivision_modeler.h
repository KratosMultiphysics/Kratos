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
#include <array>
#include <vector>

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
    ModelPart& InitializeSetup_(
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
        
    void GenerateSubdivision();

    // Returns true if the skin model part (by name) is fully contained
    // within the rectangular patch in physical XY space derived from
    // the given rect (parameter-space) and the base geometry extents.
    bool IsSkinFullyInsidePatch_(
        const std::string& rSkinModelPartName,
        const RectType& rect,
        const Parameters& geometry_base) const;

    std::vector<RectangleType> mSubdomains;
    std::vector<RefinementRegionData> mRefinementRegions;
    Model* mpModel = nullptr;

    /// Build (or refresh) a global "ConvectionDiffusionDomain" submodelpart under the given root
    /// by collecting the elements from every submodelpart whose name starts with "Patch".
    /// Also ensures the nodes used by those elements are members of the global submodelpart.
    void BuildGlobalSubModelParts(
        ModelPart& r_parent_model_part,
        const std::vector<std::string>& rTargetNames) const;

};

} // namespace Kratos
