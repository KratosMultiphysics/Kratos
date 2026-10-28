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
{ mSettings.ValidateAndAssignDefaults(GetDefaultParameters()); }

const Parameters CreateFluidPatchGapCouplingConditionsProcess::GetDefaultParameters() const
{
    return Parameters(R"({
        "background_model_part_name":"",
        "inner_model_part_name":"",
        "interface_model_part_name":"",
        "properties_id":0,
        "tolerance":1e-10,
        "allow_transient":false,
        "formulation":"stokes"
    })");
}

void CreateFluidPatchGapCouplingConditionsProcess::Execute()
{
    KRATOS_TRY
    const std::string formulation = mSettings["formulation"].GetString();
    KRATOS_ERROR_IF(formulation != "stokes" && formulation != "qsvms")
        << "Unknown fluid patch gap formulation: " << formulation << std::endl;
    const bool use_qsvms = formulation == "qsvms";
    KRATOS_ERROR_IF(use_qsvms && !mSettings["allow_transient"].GetBool())
        << "QSVMS gap coupling requires the transient BDF assembly path." << std::endl;
    auto& background = mrModel.GetModelPart(mSettings["background_model_part_name"].GetString());
    auto& inner = mrModel.GetModelPart(mSettings["inner_model_part_name"].GetString());
    auto& root = background.GetRootModelPart();
    KRATOS_ERROR_IF(&root!=&inner.GetRootModelPart()) << "Patches must share a root ModelPart." << std::endl;
    KRATOS_ERROR_IF(&background==&root || &inner==&root || &background==&inner)
        << "Gap coupling requires two independent patch submodelparts." << std::endl;
    const std::string destination = mSettings["interface_model_part_name"].GetString();
    KRATOS_ERROR_IF(destination.empty()) << "interface_model_part_name is required." << std::endl;
    if (mrModel.HasModelPart(destination)) {
        const auto& part = mrModel.GetModelPart(destination);
        KRATOS_ERROR_IF(&part.GetRootModelPart()!=&root || &part==&root ||
                        part.NumberOfElements()!=0 || part.NumberOfConditions()!=0 || part.NumberOfSubModelParts()!=0)
            << "Gap coupling destination must be an empty submodelpart; execute only once." << std::endl;
    }
    const int properties_id = mSettings["properties_id"].GetInt();
    KRATOS_ERROR_IF(properties_id<0 || !root.HasProperties(properties_id))
        << "Coupling Properties must already exist." << std::endl;
    auto p_properties = root.pGetProperties(properties_id);
    KRATOS_ERROR_IF_NOT(p_properties->Has(PENALTY_COEFFICIENT)) << "PENALTY_COEFFICIENT is required." << std::endl;
    const double coefficient = (*p_properties)[PENALTY_COEFFICIENT];
    KRATOS_ERROR_IF(!std::isfinite(coefficient) || (coefficient<0.0 && coefficient!=-1.0))
        << "Invalid PENALTY_COEFFICIENT." << std::endl;
    const auto& info = root.GetProcessInfo();
    KRATOS_ERROR_IF_NOT(info.Has(BDF_COEFFICIENTS)) << "Gap Stokes requires BDF_COEFFICIENTS (zero for steady assembly)." << std::endl;
    KRATOS_ERROR_IF(info[BDF_COEFFICIENTS].size()!=3)
        << "Gap Stokes requires three BDF coefficients." << std::endl;
    KRATOS_ERROR_IF(!mSettings["allow_transient"].GetBool() && norm_2(info[BDF_COEFFICIENTS])!=0.0)
        << "Set allow_transient for BDF2 Stokes gap assembly." << std::endl;
    for (const double value:info[BDF_COEFFICIENTS]) {
        KRATOS_ERROR_IF_NOT(std::isfinite(value)) << "Non-finite BDF coefficient." << std::endl;
    }
    if (use_qsvms) {
        KRATOS_ERROR_IF_NOT(info.Has(OSS_SWITCH) && info[OSS_SWITCH] == 0)
            << "QSVMS gap coupling currently requires OSS_SWITCH=0 (ASGS)." << std::endl;
        KRATOS_ERROR_IF(inner.NumberOfNodes()==0) << "Inner patch has no nodes." << std::endl;
    }
    std::set<ModelPart::IndexType> node_ids;
    for (auto* patch:{&background, &inner}) {
        for (const auto& node:patch->Nodes()) {
            KRATOS_ERROR_IF_NOT(node_ids.insert(node.Id()).second) << "Patches must have independent nodes." << std::endl;
            if (mSettings["allow_transient"].GetBool()) {
                KRATOS_ERROR_IF(node.GetBufferSize()<3) << "Transient gap Stokes requires three history slots." << std::endl;
                if (patch==&background && node.SolutionStepsDataHas(MESH_VELOCITY)) {
                    KRATOS_ERROR_IF(norm_2(node.FastGetSolutionStepValue(MESH_VELOCITY))>1e-12)
                        << "Gap extension requires a fixed Eulerian background." << std::endl;
                }
                if (use_qsvms && patch==&inner && node.SolutionStepsDataHas(MESH_VELOCITY)) {
                    const auto& r_mesh_velocity = node.FastGetSolutionStepValue(MESH_VELOCITY);
                    const auto& r_reference_velocity = inner.NodesBegin()->FastGetSolutionStepValue(MESH_VELOCITY);
                    KRATOS_ERROR_IF(std::abs(r_mesh_velocity[2])>1e-12 ||
                                    norm_2(r_mesh_velocity-r_reference_velocity)>1e-12)
                        << "QSVMS inner patch requires uniform planar mesh velocity." << std::endl;
                }
            }
        }
        for (const auto& element:patch->Elements()) {
            bool moving_inner = false;
            if (!use_qsvms && patch==&inner && mSettings["allow_transient"].GetBool()) {
                for (const auto& node:element.GetGeometry()) {
                    moving_inner = moving_inner || (node.SolutionStepsDataHas(MESH_VELOCITY) &&
                        norm_2(node.FastGetSolutionStepValue(MESH_VELOCITY))>1e-12);
                }
            }
            FluidPatchGapElementCheck::Check(element, use_qsvms, moving_inner);
            if (use_qsvms) {
                for (const auto& node : element.GetGeometry()) {
                    KRATOS_ERROR_IF_NOT(node.SolutionStepsDataHas(MESH_VELOCITY))
                        << "QSVMS gap coupling requires MESH_VELOCITY on all patch nodes." << std::endl;
                }
            }
            const auto& material = element.GetProperties();
            KRATOS_ERROR_IF_NOT(material.Has(DYNAMIC_VISCOSITY) &&
                                std::isfinite(material[DYNAMIC_VISCOSITY]) && material[DYNAMIC_VISCOSITY]>0.0 &&
                                material.Has(DENSITY) && std::isfinite(material[DENSITY]) && material[DENSITY]>0.0 &&
                                material.Has(CONSTITUTIVE_LAW) && material[CONSTITUTIVE_LAW]!=nullptr &&
                                material[CONSTITUTIVE_LAW]->Info()=="Newtonian2DLaw")
                << "Gap-SBM requires a constant Newtonian2DLaw material per parent." << std::endl;
        }
    }
    const auto data = FemMultipatchGapUtilities::Build(background, inner, mSettings["tolerance"].GetDouble());
    auto& output = mrModel.HasModelPart(destination)?mrModel.GetModelPart(destination):mrModel.CreateModelPart(destination);
    KRATOS_ERROR_IF(&output.GetRootModelPart()!=&root || &output==&root)
        << "Gap coupling destination must belong to the patches' root." << std::endl;
    auto& volume = output.CreateSubModelPart("GapVolumes");
    auto& physical = output.CreateSubModelPart("PhysicalInterfaces");
    auto& connectors = output.CreateSubModelPart("GapInterfaces");
    ModelPart::IndexType element_id = 1, condition_id = 1;
    for (const auto& e:root.Elements()) element_id = std::max(element_id, e.Id()+1);
    for (const auto& c:root.Conditions()) condition_id = std::max(condition_id, c.Id()+1);
    for (const auto& cell:data.Cells) {
        auto p_element = use_qsvms
            ? QSVMSGapElementUtility::Create(element_id++, *cell.pParent, cell.pQuadrature)
            : cell.pParent->Create(element_id++, cell.pQuadrature, cell.pParent->pGetProperties());
        GlobalPointersVector<Element> parents;
        parents.push_back(GlobalPointer<Element>(cell.pParent.get()));
        p_element->SetValue(NEIGHBOUR_ELEMENTS, parents); p_element->Set(ACTIVE, true);
        volume.AddElement(p_element);
        for (const auto& n:p_element->GetGeometry()) volume.AddNode(root.pGetNode(n.Id()));
    }
    const auto create = [&](const FemMultipatchGapUtilities::Segment& s, ModelPart& part, const bool Connector) {
        std::vector<ModelPart::IndexType> ids;
        std::set<ModelPart::IndexType> inserted;
        for (const auto* parent:{s.pFirst.get(), s.pSecond.get()}) {
            for (const auto& node:parent->GetGeometry()) if (inserted.insert(node.Id()).second) ids.push_back(node.Id());
        }
        auto condition = part.CreateNewCondition(Connector?"FluidPatchGapInterfaceCondition2D":"FluidPatchGapCouplingCondition2D8N",
                                               condition_id++, ids, p_properties);
        GlobalPointersVector<Element> parents;
        parents.push_back(GlobalPointer<Element>(s.pFirst.get()));
        parents.push_back(GlobalPointer<Element>(s.pSecond.get()));
        condition->SetValue(NEIGHBOUR_ELEMENTS, parents);
        Matrix points = ZeroMatrix(2, 3);
        for (std::size_t d = 0; d<3; ++d) { points(0, d) = s.Start[d]; points(1, d) = s.End[d]; }
        condition->SetValue(Connector?SURROGATE_BOUNDARY_FACE_COORDINATES:COUPLING_INTERFACE_COORDINATES, points);
        condition->Set(ACTIVE, true); part.AddNodes(ids);
    };
    for (const auto& segment:data.Physical) create(segment, physical, false);
    for (const auto& segment:data.Connectors) create(segment, connectors, true);
    KRATOS_CATCH("")
}
} // namespace Kratos
