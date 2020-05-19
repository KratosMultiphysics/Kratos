//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
//

// Application includes
#include <omp.h>
#include "custom_processes/apply_chimera_process.h"
#include "containers/model.h"
#include "utilities/builtin_timer.h"
#include "utilities/variable_utils.h"

namespace Kratos {

template <int TDim>
ApplyChimera<TDim>::ApplyChimera(ModelPart& rMainModelPart, Parameters iParameters)
    : mrMainModelPart(rMainModelPart), mParameters(iParameters)
{
    // This is only for example
    Parameters example_parameters(R"(
            {
               	"chimera_parts"   :   [
									[{
										"model_part_name":"PLEASE_SPECIFY",
                                        "search_model_part_name":"PLEASE_SPECIFY",
                                        "boundary_model_part_name":"PLEASE_SPECIFY",
										"overlap_distance":0.0
									}],
									[{
										"model_part_name":"PLEASE_SPECIFY",
                                        "search_model_part_name":"PLEASE_SPECIFY",
                                        "boundary_model_part_name":"PLEASE_SPECIFY",
										"overlap_distance":0.0
									}],
									[{
										"model_part_name":"PLEASE_SPECIFY",
                                        "search_model_part_name":"PLEASE_SPECIFY",
                                        "boundary_model_part_name":"PLEASE_SPECIFY",
										"overlap_distance":0.0
									}]
								]
            })");
    mNumberOfLevels = mParameters.size();
    KRATOS_ERROR_IF(mNumberOfLevels < 2)
        << "Chimera requires atleast one Patch !" << std::endl;

    mEchoLevel = 0;
    mReformulateEveryStep = false;
    mIsFormulated = false;
}

template <int TDim>
void ApplyChimera<TDim>::SetEchoLevel(int EchoLevel)
{
    mEchoLevel = EchoLevel;
}

template <int TDim>
void ApplyChimera<TDim>::SetReformulateEveryStep(bool Reformulate)
{
    mReformulateEveryStep = Reformulate;
}

template <int TDim>
void ApplyChimera<TDim>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;
    VariableUtils().SetFlag(ACTIVE, true, mrMainModelPart.Elements());
    // Actual execution of the functionality of this class
    if (mReformulateEveryStep || !mIsFormulated) {
        DoChimeraLoop();
        mIsFormulated = true;
    }
    KRATOS_CATCH("");
}

template <int TDim>
void ApplyChimera<TDim>::ExecuteFinalizeSolutionStep()
{
    VariableUtils().SetFlag(VISITED, false, mrMainModelPart.Nodes());
    VariableUtils().SetFlag(VISITED, false, mrMainModelPart.Elements());
    VariableUtils().SetNonHistoricalVariable(SPLIT_ELEMENT, false,
                                             mrMainModelPart.Elements());

    if (mReformulateEveryStep) {
        mrMainModelPart.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
        mIsFormulated = false;
    }
}

template <int TDim>
std::string ApplyChimera<TDim>::Info() const
{
    return "ApplyChimera";
}

template <int TDim>
void ApplyChimera<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ApplyChimera" << std::endl;
}

template <int TDim>
void ApplyChimera<TDim>::PrintData(std::ostream& rOStream) const
{
    rOStream << "ApplyChimera" << std::endl;
}

template <int TDim>
void ApplyChimera<TDim>::DoChimeraLoop()
{
    const int num_elements = static_cast<int>(mrMainModelPart.NumberOfElements());
    const auto elem_begin = mrMainModelPart.ElementsBegin();

    Parameters parameters_for_validation(R"(
                                {
                                    "model_part_name":"PLEASE_SPECIFY",
                                    "search_model_part_name":"PLEASE_SPECIFY",
                                    "boundary_model_part_name":"PLEASE_SPECIFY",
                                    "model_import_settings":{},
                                    "overlap_distance":0.0,
                                    "internal_parts_for_chimera":[]
                                }
        )");

#pragma omp parallel for
    for (int i_be = 0; i_be < num_elements; ++i_be) {
        auto i_elem = elem_begin + i_be;
        if (!i_elem->Is(VISITED)) // for multipatch
            i_elem->Set(ACTIVE, true);
    }

    const int num_nodes = static_cast<int>(mrMainModelPart.NumberOfNodes());
    const auto nodes_begin = mrMainModelPart.ElementsBegin();

#pragma omp parallel for
    for (int i_bn = 0; i_bn < num_nodes; ++i_bn) {
        auto i_node = nodes_begin + i_bn;
        i_node->Set(VISITED, false);
    }
    BuiltinTimer do_chimera_loop_time;

    int i_current_level = 0;
    for (auto& current_level : mParameters) {
        ChimeraHoleCuttingUtility::Domain domain_type =
            ChimeraHoleCuttingUtility::Domain::MAIN_BACKGROUND;
        for (auto& background_patch_param :
             current_level) { // Gives the current background
            background_patch_param.ValidateAndAssignDefaults(parameters_for_validation);
            Model& current_model = mrMainModelPart.GetModel();
            ModelPart& r_background_model_part = current_model.GetModelPart(
                background_patch_param["model_part_name"].GetString());
            // compute the outerboundary of the background to save
            if (!r_background_model_part.HasSubModelPart(mBoundaryName)) {
                auto& r_boundary_model_part =
                    r_background_model_part.CreateSubModelPart(mBoundaryName);
                BuiltinTimer extraction_time;
                if (i_current_level == 0)
                    ChimeraHoleCuttingUtility().ExtractBoundaryMesh<TDim>(
                        r_background_model_part, r_boundary_model_part);
                else
                    ChimeraHoleCuttingUtility().ExtractBoundaryMesh<TDim>(r_background_model_part, r_boundary_model_part, ChimeraHoleCuttingUtility::SideToExtract::INSIDE);
                KRATOS_INFO_IF(
                    "ApplyChimera : Extraction of background boundary took "
                    "   : ",
                    mEchoLevel > 1)
                    << extraction_time.ElapsedSeconds() << " seconds" << std::endl;
            }


            for (int i_slave_level = i_current_level + 1;
                 i_slave_level < mNumberOfLevels; ++i_slave_level) {
                for (auto& slave_patch_param :
                     mParameters[i_slave_level]) // Loop over all other slave
                                                 // patches
                {
                    slave_patch_param.ValidateAndAssignDefaults(parameters_for_validation);
                    KRATOS_INFO_IF(
                        "Formulating Chimera for the combination :: ", mEchoLevel > 0)
                        << "\n :: Background :: \n"
                        << background_patch_param << "\n :: Patch :: \n"
                        << slave_patch_param << std::endl;
                    if (i_current_level == 0) // a check to identify
                                              // computational Domain
                                              // boundary
                        domain_type = ChimeraHoleCuttingUtility::Domain::OTHER;
                    FormulateChimera(background_patch_param, slave_patch_param, domain_type);
                    mPointLocatorsMap.erase(
                        background_patch_param["model_part_name"].GetString());
                }
            }
        }
        ++i_current_level;
    }

    mNodeIdToConstraintIdsMap.clear();
    mPointLocatorsMap.clear();

    KRATOS_INFO_IF(
        "ApplyChimera : Chimera Initialization took               : ", mEchoLevel > 0)
        << do_chimera_loop_time.ElapsedSeconds() << " seconds" << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of constraints for Chimera         : ", mEchoLevel > 0)
        << mrMainModelPart.NumberOfMasterSlaveConstraints() << std::endl;
    KRATOS_INFO("End of Formulate Chimera") << std::endl;
}

template <int TDim>
void ApplyChimera<TDim>::FormulateChimera(const Parameters BackgroundParam,
                                          const Parameters PatchParameters,
                                          ChimeraHoleCuttingUtility::Domain DomainType)
{
    Model& current_model = mrMainModelPart.GetModel();
    ModelPart& r_background_model_part =
        current_model.GetModelPart(BackgroundParam["model_part_name"].GetString());
    ModelPart& r_background_boundary_model_part =
        r_background_model_part.GetSubModelPart(mBoundaryName);
    ModelPart& r_patch_model_part =
        current_model.GetModelPart(PatchParameters["model_part_name"].GetString());
    const std::string bg_search_mp_name =
        BackgroundParam["search_model_part_name"].GetString();
    auto& r_background_search_model_part =
        current_model.HasModelPart(bg_search_mp_name)
            ? current_model.GetModelPart(bg_search_mp_name)
            : r_background_model_part;

    const double overlap_bg = BackgroundParam["overlap_distance"].GetDouble();
    const double overlap_pt = PatchParameters["overlap_distance"].GetDouble();
    const double over_lap_distance = (overlap_bg > overlap_pt) ? overlap_bg : overlap_pt;

    BuiltinTimer search_creation_time;
    PointLocatorPointerType p_point_locator_on_background =
        GetPointLocator(r_background_search_model_part);
    PointLocatorPointerType p_pointer_locator_on_patch = GetPointLocator(r_patch_model_part);
    KRATOS_INFO_IF(
        "ApplyChimera : Creation of search structures took        : ", mEchoLevel > 0)
        << search_creation_time.ElapsedSeconds() << " seconds" << std::endl;
    KRATOS_ERROR_IF(over_lap_distance < 1e-12)
        << "Overlap distance should be a positive and non-zero number." << std::endl;

    ModelPart& r_hole_model_part = r_background_model_part.CreateSubModelPart(mHoleName);
    ModelPart& r_hole_boundary_model_part =
        r_hole_model_part.CreateSubModelPart(mBoundaryName);

    auto& r_patch_boundary_model_part = ExtractPatchBoundary(
        PatchParameters, r_background_boundary_model_part, DomainType);

    BuiltinTimer bg_distance_calc_time;
    ChimeraDistanceCalculationUtility<TDim>::CalculateDistance(
        r_background_model_part, r_patch_boundary_model_part);
    KRATOS_INFO_IF(
        "Distance calculation on background took                  : ", mEchoLevel > 0)
        << bg_distance_calc_time.ElapsedSeconds() << " seconds" << std::endl;

    BuiltinTimer hole_creation_time;
    ChimeraHoleCuttingUtility().CreateHoleAfterDistance<TDim>(
        r_background_model_part, r_hole_model_part, r_hole_boundary_model_part, over_lap_distance);
    KRATOS_INFO_IF(
        "ApplyChimera : Hole creation took                        : ", mEchoLevel > 0)
        << hole_creation_time.ElapsedSeconds() << " seconds" << std::endl;

    const int n_elements = static_cast<int>(r_hole_model_part.NumberOfElements());
#pragma omp parallel for
    for (int i_elem = 0; i_elem < n_elements; ++i_elem) {
        ModelPart::ElementsContainerType::iterator it_elem =
            r_hole_model_part.ElementsBegin() + i_elem;
        it_elem->Set(VISITED, true); // for multipatch
    }

    BuiltinTimer mpc_time;
    ApplyContinuityWithMpcs(r_patch_boundary_model_part, *p_point_locator_on_background);
    ApplyContinuityWithMpcs(r_hole_boundary_model_part, *p_pointer_locator_on_patch);
    KRATOS_INFO_IF(
        "ApplyChimera : Creation of MPC for chimera took          : ", mEchoLevel > 0)
        << mpc_time.ElapsedSeconds() << " seconds" << std::endl;

    r_hole_boundary_model_part.RemoveSubModelPart(mBoundaryName);
    r_background_model_part.RemoveSubModelPart(mHoleName);
    r_patch_model_part.RemoveSubModelPart(mModifiedName);
}

template <int TDim>
void ApplyChimera<TDim>::CreateConstraintIds(std::vector<int>& rIdVector,
                                             const IndexType NumberOfConstraintsRequired)
{
    IndexType max_constraint_id = 0;
    // Get current maximum constraint ID
    if (mrMainModelPart.MasterSlaveConstraints().size() != 0) {
        mrMainModelPart.MasterSlaveConstraints().Sort();
        ModelPart::MasterSlaveConstraintContainerType::iterator it =
            mrMainModelPart.MasterSlaveConstraintsEnd() - 1;
        max_constraint_id = (*it).Id();
        ++max_constraint_id;
    }

    // Now create a vector size NumberOfConstraintsRequired
    rIdVector.resize(NumberOfConstraintsRequired * (TDim + 1));
    std::iota(std::begin(rIdVector), std::end(rIdVector), max_constraint_id); // Fill with consecutive integers
}

template <int TDim>
void ApplyChimera<TDim>::ApplyContinuityWithMpcs(ModelPart& rBoundaryModelPart,
                                                 PointLocatorType& pBinLocator)
{
}

template <int TDim>
template <typename TVariableType>
void ApplyChimera<TDim>::AddMasterSlaveRelation(
    MasterSlaveConstraintContainerType& rMasterSlaveContainer,
    const LinearMasterSlaveConstraint& rCloneConstraint,
    unsigned int ConstraintId,
    NodeType& rMasterNode,
    TVariableType& rMasterVariable,
    NodeType& rSlaveNode,
    TVariableType& rSlaveVariable,
    const double Weight,
    const double Constant)
{
    rSlaveNode.Set(SLAVE);
    ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint =
        rCloneConstraint.Create(ConstraintId, rMasterNode, rMasterVariable,
                                rSlaveNode, rSlaveVariable, Weight, Constant);
    p_new_constraint->Set(TO_ERASE);
    mNodeIdToConstraintIdsMap[rSlaveNode.Id()].push_back(ConstraintId);
    rMasterSlaveContainer.push_back(p_new_constraint);
}

template <int TDim>
template <typename TVariableType>
void ApplyChimera<TDim>::ApplyContinuityWithElement(
    Geometry<NodeType>& rGeometry,
    NodeType& rBoundaryNode,
    Vector& rShapeFuncWeights,
    TVariableType& rVariable,
    unsigned int StartIndex,
    std::vector<int>& rConstraintIdVector,
    MasterSlaveConstraintContainerType& rMsContainer)
{
    const auto& r_clone_constraint = new LinearMasterSlaveConstraint();
    // Initialise the boundary nodes dofs to 0 at ever time steps
    rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) = 0.0;

    for (IndexType i = 0; i < rGeometry.size(); i++) {
        // Interpolation of rVariable
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) +=
            rGeometry[i].GetDof(rVariable).GetSolutionStepValue(0) * rShapeFuncWeights[i];
        AddMasterSlaveRelation(rMsContainer, *r_clone_constraint,
                               rConstraintIdVector[StartIndex++], rGeometry[i], rVariable,
                               rBoundaryNode, rVariable, rShapeFuncWeights[i]);
    } // end of loop over host element nodes

    // Setting the buffer 1 same buffer 0
    rBoundaryNode.FastGetSolutionStepValue(rVariable, 1) =
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 0);
}

template <int TDim>
void ApplyChimera<TDim>::ReserveMemoryForConstraintContainers(
    ModelPart& rModelPart, MasterSlaveContainerVectorType& rContainerVector)
{
#pragma omp parallel
    {
        const IndexType num_constraints_per_thread =
            (rModelPart.NumberOfNodes() * 4) / omp_get_num_threads();
#pragma omp single
        {
            rContainerVector.resize(omp_get_num_threads());
            for (auto& container : rContainerVector)
                container.reserve(num_constraints_per_thread);
        }
    }
}

template <int TDim>
int ApplyChimera<TDim>::RemoveExistingConstraintsForNode(ModelPart::NodeType& rBoundaryNode)
{
    ConstraintIdsVectorType constrainIds_for_the_node;
    int removed_counter = 0;
    constrainIds_for_the_node = mNodeIdToConstraintIdsMap[rBoundaryNode.Id()];
    for (auto const& constraint_id : constrainIds_for_the_node) {
#pragma omp critical
        {
            mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
            removed_counter++;
        }
    }
    constrainIds_for_the_node.clear();
    return removed_counter;
}

template <int TDim>
void ApplyChimera<TDim>::AddConstraintsToModelpart(ModelPart& rModelPart,
                                                   MasterSlaveContainerVectorType& rContainerVector)
{
    IndexType n_total_constraints = 0;
    for (auto& container : rContainerVector) {
        const int n_constraints = static_cast<int>(container.size());
        n_total_constraints += n_constraints;
    }

    auto& constraints = rModelPart.MasterSlaveConstraints();
    constraints.reserve(n_total_constraints);
    auto& constraints_data = constraints.GetContainer();
    for (auto& container : rContainerVector) {
        constraints_data.insert(constraints_data.end(), container.ptr_begin(),
                                container.ptr_end());
    }
    constraints.Sort();
}

template <int TDim>
void ApplyChimera<TDim>::FormulateConstraints(
    ModelPart& rBoundaryModelPart,
    PointLocatorType& rBinLocator,
    MasterSlaveContainerVectorType& rVelocityMasterSlaveContainerVector,
    MasterSlaveContainerVectorType& rPressureMasterSlaveContainerVector)
{
    std::vector<int> vector_of_non_found_nodes;
    const int n_boundary_nodes = static_cast<int>(rBoundaryModelPart.Nodes().size());
    std::vector<int> constraints_id_vector;
    int num_constraints_required = (TDim + 1) * (rBoundaryModelPart.Nodes().size());
    CreateConstraintIds(constraints_id_vector, num_constraints_required);

    IndexType found_counter = 0;
    IndexType removed_counter = 0;

    BuiltinTimer loop_over_b_nodes;
#pragma omp parallel for shared(constraints_id_vector,               \
                                rVelocityMasterSlaveContainerVector, \
                                rPressureMasterSlaveContainerVector, \
                                rBinLocator) reduction(+ : found_counter)
    for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn) {
        ModelPart::NodesContainerType::iterator i_boundary_node =
            rBoundaryModelPart.NodesBegin() + i_bn;
        NodeType& r_boundary_node = *(*(i_boundary_node.base()));
        unsigned int start_constraint_id = i_bn * (TDim + 1) * (TDim + 1);
        Element::Pointer r_host_element;
        Vector weights;
        bool is_found = SearchNode(rBinLocator, r_boundary_node, r_host_element, weights);
        if (is_found) {
            auto& ms_velocity_container =
                rVelocityMasterSlaveContainerVector[omp_get_thread_num()];
            auto& ms_pressure_container =
                rPressureMasterSlaveContainerVector[omp_get_thread_num()];
            removed_counter += RemoveExistingConstraintsForNode(r_boundary_node);
            MakeConstraints(r_boundary_node, r_host_element, weights,
                            ms_velocity_container, ms_pressure_container,
                            constraints_id_vector, start_constraint_id);
            found_counter += 1;
        }
    }

    double loop_time = loop_over_b_nodes.ElapsedSeconds();
    KRATOS_INFO_IF(
        "ApplyChimera : Loop over boundary nodes took             : ", mEchoLevel > 0)
        << loop_time << " seconds" << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of Boundary nodes                  : ", mEchoLevel > 1)
        << n_boundary_nodes << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of Boundary nodes found            : ", mEchoLevel > 1)
        << found_counter << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of Boundary nodes not found        : ", mEchoLevel > 1)
        << n_boundary_nodes - found_counter << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of constraints made                : ", mEchoLevel > 1)
        << found_counter * 9 << std::endl;
    KRATOS_INFO_IF(
        "ApplyChimera : Number of constraints removed             : ", mEchoLevel > 1)
        << removed_counter << std::endl;
}

template <int TDim>
ModelPart& ApplyChimera<TDim>::ExtractPatchBoundary(const Parameters PatchParameters,
                                                    ModelPart& rBackgroundBoundaryModelpart,
                                                    const ChimeraHoleCuttingUtility::Domain DomainType)
{
    Model& current_model = mrMainModelPart.GetModel();
    const std::string patch_boundary_mp_name =
        PatchParameters["boundary_model_part_name"].GetString();

    if (!current_model.HasModelPart(patch_boundary_mp_name)) {
        ModelPart& r_patch_model_part =
            current_model.GetModelPart(PatchParameters["model_part_name"].GetString());
        ModelPart& r_modified_patch_model_part =
            r_patch_model_part.CreateSubModelPart(mModifiedName);
        ModelPart& r_modified_patch_boundary_model_part =
            r_modified_patch_model_part.CreateSubModelPart(
                mBoundaryName + r_modified_patch_model_part.Name());
        BuiltinTimer distance_calc_time_patch;
        ChimeraDistanceCalculationUtility<TDim>::CalculateDistance(
            r_patch_model_part, rBackgroundBoundaryModelpart);
        KRATOS_INFO_IF(
            "Distance calculation on patch took                       : ", mEchoLevel > 0)
            << distance_calc_time_patch.ElapsedSeconds() << " seconds" << std::endl;
        // TODO: Below is brutforce. Check if the boundary of bg is actually
        // cutting the patch.
        BuiltinTimer rem_out_domain_time;
        ChimeraHoleCuttingUtility().RemoveOutOfDomainElements<TDim>(
            r_patch_model_part, r_modified_patch_model_part, DomainType,
            ChimeraHoleCuttingUtility::SideToExtract::INSIDE);
        KRATOS_INFO_IF(
            "ApplyChimera : Removing out of domain patch took         : ", mEchoLevel > 0)
            << rem_out_domain_time.ElapsedSeconds() << " seconds" << std::endl;

        BuiltinTimer patch_boundary_extraction_time;
        ChimeraHoleCuttingUtility().ExtractBoundaryMesh<TDim>(
            r_modified_patch_model_part, r_modified_patch_boundary_model_part);
        KRATOS_INFO_IF(
            "ApplyChimera : Extraction of patch boundary took         : ", mEchoLevel > 0)
            << patch_boundary_extraction_time.ElapsedSeconds() << " seconds"
            << std::endl;

        return r_modified_patch_boundary_model_part;
    }
    else {
        return current_model.GetModelPart(patch_boundary_mp_name);
    }
}

template <int TDim>
typename ApplyChimera<TDim>::PointLocatorPointerType ApplyChimera<TDim>::GetPointLocator(ModelPart& rModelPart)
{
    if (mPointLocatorsMap.count(rModelPart.Name()) == 0) {
        PointLocatorPointerType p_point_locator =
            Kratos::make_shared<PointLocatorType>(rModelPart);
        p_point_locator->UpdateSearchDatabase();
        mPointLocatorsMap[rModelPart.Name()] = p_point_locator;
        return p_point_locator;
    }
    else {
        auto p_point_locator = mPointLocatorsMap.at(rModelPart.Name());
        if (mReformulateEveryStep)
            p_point_locator->UpdateSearchDatabase();
        return p_point_locator;
    }
}

template <int TDim>
bool ApplyChimera<TDim>::SearchNode(PointLocatorType& rBinLocator,
                                    NodeType& rNodeToFind,
                                    Element::Pointer& prHostElement,
                                    Vector& rWeights)
{
    const int max_results = 10000;
    typename PointLocatorType::ResultContainerType results(max_results);
    typename PointLocatorType::ResultIteratorType result_begin =
        results.begin();

    bool is_found = false;
    is_found = rBinLocator.FindPointOnMesh(rNodeToFind.Coordinates(), rWeights,
                                           prHostElement, result_begin, max_results);

    bool node_coupled = false;
    if (rNodeToFind.IsDefined(VISITED))
        node_coupled = rNodeToFind.Is(VISITED);

    if (node_coupled && is_found) {
        auto& constrainIds_for_the_node = mNodeIdToConstraintIdsMap[rNodeToFind.Id()];
        for (auto const& constraint_id : constrainIds_for_the_node) {
#pragma omp critical
            {
                mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
            }
        }
        constrainIds_for_the_node.clear();
        rNodeToFind.Set(VISITED, false);
    }

    return is_found;
}

template <int TDim>
void ApplyChimera<TDim>::MakeConstraints(
    NodeType& rNodeToFind,
    Element::Pointer& rHostElement,
    Vector& rWeights,
    MasterSlaveConstraintContainerType& rVelocityMsConstraintsVector,
    MasterSlaveConstraintContainerType& rPressureMsConstraintsVector,
    std::vector<int>& rConstraintIdVector,
    const IndexType StartConstraintId)
{
    Geometry<NodeType>& r_geom = rHostElement->GetGeometry();
    int init_index = 0;
    ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, VELOCITY_X,
                               StartConstraintId + init_index,
                               rConstraintIdVector, rVelocityMsConstraintsVector);
    init_index += (TDim + 1);
    ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, VELOCITY_Y,
                               StartConstraintId + init_index,
                               rConstraintIdVector, rVelocityMsConstraintsVector);
    init_index += (TDim + 1);
    if (TDim == 3) {
        ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, VELOCITY_Z,
                                   StartConstraintId + init_index, rConstraintIdVector,
                                   rVelocityMsConstraintsVector);
        init_index += (TDim + 1);
    }
    ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, PRESSURE,
                               StartConstraintId + init_index,
                               rConstraintIdVector, rPressureMsConstraintsVector);
}

// Template declarations
template class ApplyChimera<2>;
template class ApplyChimera<3>;

} // namespace Kratos.
