//     |  /           |
//     ' /   __| _` | __|  _ \   __|
//     . \  |   (   | |   (   |\__ `
//    _|\_\_|  \__,_|\__|\___/ ____/
//                    Multi-Physics
//
//   Fluid Dynamics Application
//
//   License:         BSD License
//                    Kratos default license: kratos/license.txt
//
//   Main author:     Nicolò Antonelli
//

#include "custom_processes/create_fluid_patch_gap_coupling_conditions_process.h"
#include "custom_utilities/fluid_patch_gap_element_check.h"
#include "custom_utilities/qsvms_gap_element_utility.h"
#include "utilities/multipatch_gap_utilities.h"

namespace Kratos {
CreateFluidPatchGapCouplingConditionsProcess::CreateFluidPatchGapCouplingConditionsProcess(
    Model& rModel, Parameters Settings) : mrModel(rModel), mSettings(Settings.Clone())
{
    KRATOS_TRY

    auto defaults = GetDefaultParameters();
    // Add the optional key to the validation schema only when supplied;
    // this value is not used as a default interface location.
    if (mSettings.Has("interface_height")) {
        defaults.AddDouble("interface_height", 0.0);
    }
    mSettings.ValidateAndAssignDefaults(defaults);
    const auto formulation = mSettings["formulation"].GetString();
    KRATOS_ERROR_IF((formulation == "two_phase_stokes" || formulation == "two_phase_qsvms") &&
                    !mSettings.Has("interface_height"))
        << "The horizontal two-phase formulation '" << formulation
        << "' requires an explicit 'interface_height' setting." << std::endl;

    KRATOS_CATCH("")
}

const Parameters CreateFluidPatchGapCouplingConditionsProcess::GetDefaultParameters() const
{
    return Parameters(R"({
        "background_model_part_name":"",
        "inner_model_part_name":"",
        "interface_model_part_name":"",
        "properties_id":0,
        "tolerance":1e-10,
        "allow_transient":false,
        "formulation":"stokes",
        "level_set_contour":[]
    })");
}

void CreateFluidPatchGapCouplingConditionsProcess::Execute()
{
    KRATOS_TRY
    const std::string formulation = mSettings["formulation"].GetString();
    KRATOS_ERROR_IF(formulation != "stokes" && formulation != "qsvms" &&
                    formulation != "two_phase_stokes" && formulation != "two_phase_qsvms" &&
                    formulation != "two_phase_level_set_stokes" && formulation != "two_phase_level_set_qsvms")
        << "Unknown fluid patch gap formulation: " << formulation << std::endl;
    const bool level_set = formulation == "two_phase_level_set_stokes" || formulation == "two_phase_level_set_qsvms";
    const bool use_qsvms = formulation == "qsvms" || formulation == "two_phase_qsvms" || formulation == "two_phase_level_set_qsvms";
    const bool two_phase = formulation == "two_phase_stokes" || formulation == "two_phase_qsvms" || level_set;
    KRATOS_ERROR_IF(two_phase && mSettings["allow_transient"].GetBool() &&
                    formulation != "two_phase_level_set_qsvms")
        << "Transient two-phase gap coupling requires the experimental level-set QSVMS path." << std::endl;
    KRATOS_ERROR_IF(use_qsvms && !two_phase && !mSettings["allow_transient"].GetBool())
        << "QSVMS gap coupling requires the transient BDF assembly path." << std::endl;
    auto& r_background = mrModel.GetModelPart(mSettings["background_model_part_name"].GetString());
    auto& r_inner = mrModel.GetModelPart(mSettings["inner_model_part_name"].GetString());
    auto& r_root = r_background.GetRootModelPart();
    KRATOS_ERROR_IF(&r_root!=&r_inner.GetRootModelPart()) << "Patches must share a root ModelPart." << std::endl;
    KRATOS_ERROR_IF(&r_background==&r_root || &r_inner==&r_root || &r_background==&r_inner)
        << "Gap coupling requires two independent patch submodelparts." << std::endl;
    const std::string destination = mSettings["interface_model_part_name"].GetString();
    KRATOS_ERROR_IF(destination.empty()) << "interface_model_part_name is required." << std::endl;
    if (mrModel.HasModelPart(destination)) {
        const auto& r_part = mrModel.GetModelPart(destination);
        KRATOS_ERROR_IF(&r_part.GetRootModelPart()!=&r_root || &r_part==&r_root ||
                        r_part.NumberOfElements()!=0 || r_part.NumberOfConditions()!=0 || r_part.NumberOfSubModelParts()!=0)
            << "Gap coupling destination must be an empty submodelpart; execute only once." << std::endl;
    }
    const int properties_id = mSettings["properties_id"].GetInt();
    KRATOS_ERROR_IF(properties_id<0 || !r_root.HasProperties(properties_id))
        << "Coupling Properties must already exist." << std::endl;
    auto p_properties = r_root.pGetProperties(properties_id);
    KRATOS_ERROR_IF_NOT(p_properties->Has(PENALTY_COEFFICIENT)) << "PENALTY_COEFFICIENT is required." << std::endl;
    const double coefficient = (*p_properties)[PENALTY_COEFFICIENT];
    KRATOS_ERROR_IF(!std::isfinite(coefficient) || (coefficient<0.0 && coefficient!=-1.0))
        << "Invalid PENALTY_COEFFICIENT." << std::endl;
    const auto& r_info = r_root.GetProcessInfo();
    KRATOS_ERROR_IF_NOT(r_info.Has(BDF_COEFFICIENTS)) << "Gap Stokes requires BDF_COEFFICIENTS (zero for steady assembly)." << std::endl;
    KRATOS_ERROR_IF(r_info[BDF_COEFFICIENTS].size()!=3)
        << "Gap Stokes requires three BDF coefficients." << std::endl;
    KRATOS_ERROR_IF(!mSettings["allow_transient"].GetBool() && norm_2(r_info[BDF_COEFFICIENTS])!=0.0)
        << "Set allow_transient for BDF2 Stokes gap assembly." << std::endl;
    for (const double value:r_info[BDF_COEFFICIENTS]) {
        KRATOS_ERROR_IF_NOT(std::isfinite(value)) << "Non-finite BDF coefficient." << std::endl;
    }
    if (use_qsvms) {
        KRATOS_ERROR_IF_NOT(r_info.Has(OSS_SWITCH) && r_info[OSS_SWITCH] == 0)
            << "QSVMS gap coupling currently requires OSS_SWITCH=0 (ASGS)." << std::endl;
        KRATOS_ERROR_IF(r_inner.NumberOfNodes()==0) << "Inner patch has no nodes." << std::endl;
        if (two_phase && !mSettings["allow_transient"].GetBool()) {
            KRATOS_ERROR_IF_NOT(r_info.Has(DYNAMIC_TAU) && r_info[DYNAMIC_TAU] == 0.0)
                << "Stationary two-phase QSVMS requires DYNAMIC_TAU=0." << std::endl;
        }
    }
    std::set<ModelPart::IndexType> node_ids;
    for (auto* p_patch:{&r_background, &r_inner}) {
        for (const auto& r_node:p_patch->Nodes()) {
            KRATOS_ERROR_IF_NOT(node_ids.insert(r_node.Id()).second) << "Patches must have independent nodes." << std::endl;
            if (mSettings["allow_transient"].GetBool()) {
                KRATOS_ERROR_IF(r_node.GetBufferSize()<3) << "Transient gap Stokes requires three history slots." << std::endl;
                if (p_patch==&r_background && r_node.SolutionStepsDataHas(MESH_VELOCITY)) {
                    KRATOS_ERROR_IF(norm_2(r_node.FastGetSolutionStepValue(MESH_VELOCITY))>1e-12)
                        << "Gap extension requires a fixed Eulerian background." << std::endl;
                }
                if (use_qsvms && p_patch==&r_inner && r_node.SolutionStepsDataHas(MESH_VELOCITY)) {
                    const auto& r_mesh_velocity = r_node.FastGetSolutionStepValue(MESH_VELOCITY);
                    const auto& r_reference_velocity = r_inner.NodesBegin()->FastGetSolutionStepValue(MESH_VELOCITY);
                    KRATOS_ERROR_IF(std::abs(r_mesh_velocity[2])>1e-12 ||
                                    norm_2(r_mesh_velocity-r_reference_velocity)>1e-12)
                        << "QSVMS inner patch requires uniform planar mesh velocity." << std::endl;
                }
            }
        }
        for (const auto& r_element:p_patch->Elements()) {
            bool moving_inner = false;
            if (!use_qsvms && p_patch==&r_inner && mSettings["allow_transient"].GetBool()) {
                for (const auto& r_node:r_element.GetGeometry()) {
                    moving_inner = moving_inner || (r_node.SolutionStepsDataHas(MESH_VELOCITY) &&
                        norm_2(r_node.FastGetSolutionStepValue(MESH_VELOCITY))>1e-12);
                }
            }
            FluidPatchGapElementCheck::Check(r_element, use_qsvms, moving_inner);
            if (use_qsvms) {
                for (const auto& r_node : r_element.GetGeometry()) {
                    KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(MESH_VELOCITY))
                        << "QSVMS gap coupling requires MESH_VELOCITY on all patch nodes." << std::endl;
                    if (two_phase) {
                        KRATOS_ERROR_IF(norm_2(r_node.FastGetSolutionStepValue(MESH_VELOCITY)) > 1e-12)
                            << "Two-phase QSVMS requires fixed meshes." << std::endl;
                        KRATOS_ERROR_IF(r_node.GetBufferSize() < 3)
                            << "TimeIntegratedQSVMS requires three allocated solution slots." << std::endl;
                    }
                }
            }
            const auto& r_material = r_element.GetProperties();
            KRATOS_ERROR_IF_NOT(r_material.Has(DYNAMIC_VISCOSITY) &&
                                std::isfinite(r_material[DYNAMIC_VISCOSITY]) && r_material[DYNAMIC_VISCOSITY]>0.0 &&
                                r_material.Has(DENSITY) && std::isfinite(r_material[DENSITY]) && r_material[DENSITY]>0.0 &&
                                r_material.Has(CONSTITUTIVE_LAW) && r_material[CONSTITUTIVE_LAW]!=nullptr &&
                                r_material[CONSTITUTIVE_LAW]->Info()=="Newtonian2DLaw")
                << "Gap-SBM requires a constant Newtonian2DLaw material per parent." << std::endl;
        }
    }
    FemMultipatchGapUtilities::GapData data;
    const double tolerance = mSettings["tolerance"].GetDouble();
    if (two_phase) {
        Matrix contour;
        if (level_set) {
            const auto coordinates = mSettings["level_set_contour"];
            contour.resize(coordinates.size(), 2, false);
            for (std::size_t i = 0; i < coordinates.size(); ++i) {
                KRATOS_ERROR_IF(coordinates[i].size() != 2) << "Each level-set contour point must have two coordinates." << std::endl;
                for (std::size_t d = 0; d < 2; ++d) contour(i, d) = coordinates[i][d].GetDouble();
            }
        }
        auto water = level_set
            ? FemMultipatchGapUtilities::BuildLevelSetPhase(r_background, contour, true, tolerance)
            : FemMultipatchGapUtilities::BuildHorizontalPhase(r_background, mSettings["interface_height"].GetDouble(), true, tolerance);
        auto air = level_set
            ? FemMultipatchGapUtilities::BuildLevelSetPhase(r_inner, contour, false, tolerance)
            : FemMultipatchGapUtilities::BuildHorizontalPhase(r_inner, mSettings["interface_height"].GetDouble(), false, tolerance);
        KRATOS_ERROR_IF(water.Physical.size() != air.Physical.size())
            << "Two-phase interface requires matching horizontal partitions." << std::endl;
        for (std::size_t i = 0; i < water.Physical.size(); ++i) {
            const auto& r_water_segment = water.Physical[i];
            const auto& r_air_segment = air.Physical[i];
            KRATOS_ERROR_IF(norm_2(r_water_segment.Start - r_air_segment.End) > tolerance || norm_2(r_water_segment.End - r_air_segment.Start) > tolerance)
                << "Two-phase interface segments do not match." << std::endl;
            data.Physical.push_back({r_water_segment.pFirst, r_air_segment.pFirst, r_water_segment.Start, r_water_segment.End});
        }
        data.Cells = std::move(water.Cells);
        data.Cells.insert(data.Cells.end(), air.Cells.begin(), air.Cells.end());
        data.Connectors = std::move(water.Connectors);
        data.Connectors.insert(data.Connectors.end(), air.Connectors.begin(), air.Connectors.end());
    } else {
        data = FemMultipatchGapUtilities::Build(r_background, r_inner, tolerance);
    }
    auto& r_output = mrModel.HasModelPart(destination)?mrModel.GetModelPart(destination):mrModel.CreateModelPart(destination);
    KRATOS_ERROR_IF(&r_output.GetRootModelPart()!=&r_root || &r_output==&r_root)
        << "Gap coupling destination must belong to the patches' root." << std::endl;
    auto& r_volume = r_output.CreateSubModelPart("GapVolumes");
    auto& r_physical = r_output.CreateSubModelPart("PhysicalInterfaces");
    auto& r_connectors = r_output.CreateSubModelPart("GapInterfaces");
    ModelPart::IndexType element_id = 1, condition_id = 1;
    for (const auto& r_element:r_root.Elements()) element_id = std::max(element_id, r_element.Id()+1);
    for (const auto& r_condition:r_root.Conditions()) condition_id = std::max(condition_id, r_condition.Id()+1);
    for (const auto& r_cell:data.Cells) {
        auto p_element = use_qsvms
            ? QSVMSGapElementUtility::Create(element_id++, *r_cell.pParent, r_cell.pQuadrature)
            : r_cell.pParent->Create(element_id++, r_cell.pQuadrature, r_cell.pParent->pGetProperties());
        GlobalPointersVector<Element> parents;
        parents.push_back(GlobalPointer<Element>(r_cell.pParent.get()));
        p_element->SetValue(NEIGHBOUR_ELEMENTS, parents); p_element->Set(ACTIVE, true);
        r_volume.AddElement(p_element);
        for (const auto& r_node:p_element->GetGeometry()) r_volume.AddNode(r_root.pGetNode(r_node.Id()));
    }
    const auto create = [&](const FemMultipatchGapUtilities::Segment& rSegment, ModelPart& rPart, const bool Connector) {
        std::vector<ModelPart::IndexType> ids;
        std::set<ModelPart::IndexType> inserted;
        for (const auto* p_parent:{rSegment.pFirst.get(), rSegment.pSecond.get()}) {
            for (const auto& r_node:p_parent->GetGeometry()) if (inserted.insert(r_node.Id()).second) ids.push_back(r_node.Id());
        }
        const std::string condition_name = Connector ? "FluidPatchGapInterfaceCondition2D" :
            (two_phase ? "FluidTwoPhaseGapCouplingCondition2D8N" : "FluidPatchGapCouplingCondition2D8N");
        auto p_condition = rPart.CreateNewCondition(condition_name,
                                                    condition_id++, ids, p_properties);
        GlobalPointersVector<Element> parents;
        parents.push_back(GlobalPointer<Element>(rSegment.pFirst.get()));
        parents.push_back(GlobalPointer<Element>(rSegment.pSecond.get()));
        p_condition->SetValue(NEIGHBOUR_ELEMENTS, parents);
        Matrix points = ZeroMatrix(2, 3);
        for (std::size_t d = 0; d<3; ++d) { points(0, d) = rSegment.Start[d]; points(1, d) = rSegment.End[d]; }
        p_condition->SetValue(Connector?SURROGATE_BOUNDARY_FACE_COORDINATES:COUPLING_INTERFACE_COORDINATES, points);
        p_condition->Set(ACTIVE, true); rPart.AddNodes(ids);
    };
    for (const auto& r_segment:data.Physical) create(r_segment, r_physical, false);
    for (const auto& r_segment:data.Connectors) create(r_segment, r_connectors, true);
    KRATOS_CATCH("")
}
} // namespace Kratos
