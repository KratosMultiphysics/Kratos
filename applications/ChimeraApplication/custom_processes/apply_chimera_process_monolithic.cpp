// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//					 Navaneeth K Narayanan
//					 Rishith Ellath Meethal
// ==============================================================================
//

// System includes

// Application includes
#include "custom_processes/apply_chimera_process_monolithic.h"

namespace Kratos
{

template <int TDim>
ApplyChimeraProcessMonolithic<TDim>::ApplyChimeraProcessMonolithic(ModelPart &rMainModelPart, Parameters iParameters) : Process(), mrMainModelPart(rMainModelPart), mParameters(iParameters)
{
    Parameters default_parameters(R"(
            {
               	"chimera"   :   [
									[{
										"model_part_name":"GENERIC_background",
										"model_part_inside_boundary_name" :"GENERIC_domainboundary",
										"overlap_distance":0.045
									}],
									[{
										"model_part_name":"GENERIC_patch_1_1",
										"model_part_inside_boundary_name":"GENERIC_structure_1_1",
										"overlap_distance":0.045
									}],
									[{
										"model_part_name":"GENERIC_patch_2_1",
										"model_part_inside_boundary_name":"GENERIC_strcuture2_1",
										"overlap_distance":0.045
									}]
								]
            })");

    mNumberOfLevels = mParameters.size();
    for (int i = 0; i < mNumberOfLevels; i++)
        mLevelTable.push_back(mParameters[i].size());

    ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
    this->mpHoleCuttingUtility = ChimeraHoleCuttingUtility::Pointer(new ChimeraHoleCuttingUtility());
    this->mpCalculateDistanceProcess = typename CustomCalculateSignedDistanceProcess<TDim>::Pointer(new CustomCalculateSignedDistanceProcess<TDim>());
}

template <int TDim>
ApplyChimeraProcessMonolithic<TDim>::~ApplyChimeraProcessMonolithic()
{
}

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;
    // Actual execution of the functionality of this class
    for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
        it->SetValue(SPLIT_ELEMENT, false);
    DoChimeraLoop();
    KRATOS_CATCH("");
}

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::ExecuteFinalizeSolutionStep()
{
    Clear();
    //for multipatch
    for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
    {
        it->Set(VISITED, false);
        it->SetValue(SPLIT_ELEMENT, false);
    }

    mrMainModelPart.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
}

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::DoChimeraLoop() //selecting patch and background combination for chimera method
{
    for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
    {
        if (!it->Is(VISITED)) //for multipatch
            it->Set(ACTIVE, true);
    }

    for (ModelPart::NodesContainerType::iterator it = mrMainModelPart.NodesBegin(); it != mrMainModelPart.NodesEnd(); ++it)
        it->Set(VISITED, false);

    int MainDomainOrNot = 1;

    for (int BG_i = 0; BG_i < mNumberOfLevels; BG_i++) // Iteration for selecting background
    {
        for (int BG_j = 0; BG_j < mLevelTable[BG_i]; BG_j++) //TODO change the names
        {
            for (int patch_i = BG_i + 1; patch_i < mNumberOfLevels; patch_i++) // Iteration for selecting patch
            {
                for (int patch_j = 0; patch_j < mLevelTable[patch_i]; patch_j++)
                {
                    m_background_model_part_name = mParameters[BG_i][BG_j]["model_part_name"].GetString();
                    m_domain_boundary_model_part_name = mParameters[BG_i][BG_j]["model_part_inside_boundary_name"].GetString();
                    m_patch_model_part_name = mParameters[patch_i][patch_j]["model_part_name"].GetString();
                    m_patch_inside_boundary_model_part_name = mParameters[patch_i][patch_j]["model_part_inside_boundary_name"].GetString();

                    double mesh_size_1 = mParameters[BG_i][BG_j]["overlap_distance"].GetDouble();
                    double mesh_size_2 = mParameters[patch_i][patch_j]["overlap_distance"].GetDouble();

                    if (mesh_size_1 > mesh_size_2)
                        mOverlapDistance = mesh_size_1;
                    else
                        mOverlapDistance = mesh_size_2;

                    KRATOS_INFO("Formulating Chimera for the combination \n background::") << m_background_model_part_name << "  \t Patch::" << m_patch_model_part_name << std::endl;

                    MainDomainOrNot = 1;
                    if (BG_i == 0) // a check to identify computational Domain boundary
                        MainDomainOrNot = -1;

                    FormulateChimera(MainDomainOrNot);
                }
            }
        }
    }
    KRATOS_INFO("End of chimera loop") << std::endl;

    KRATOS_INFO("Total number of constraints created so far") << mrMainModelPart.NumberOfMasterSlaveConstraints() << std::endl;
}

//Apply Chimera with or without overlap
template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::FormulateChimera(int MainDomainOrNot)
{
    ModelPart &r_background_model_part = mrMainModelPart.GetSubModelPart(m_background_model_part_name);
    ModelPart &r_patch_model_part = mrMainModelPart.GetSubModelPart(m_patch_model_part_name);
    ModelPart &r_domain_boundary_model_part = mrMainModelPart.GetSubModelPart(m_domain_boundary_model_part_name);
    ModelPart &r_patch_inside_boundary_model_part = mrMainModelPart.GetSubModelPart(m_patch_inside_boundary_model_part_name);

    PointLocatorPointerType p_point_locator_on_background = PointLocatorPointerType(new PointLocatorType(r_background_model_part));
    p_point_locator_on_background->UpdateSearchDatabase();
    PointLocatorPointerType p_pointer_locator_on_patch = PointLocatorPointerType(new PointLocatorType(r_patch_model_part));
    p_pointer_locator_on_patch->UpdateSearchDatabase();

    const double eps = 1e-12;
    if (mOverlapDistance < eps)
        KRATOS_THROW_ERROR("", "Overlap distance should be a positive and non-zero number \n", "");

    if (mOverlapDistance > eps)
    {
        Model &current_model = mrMainModelPart.GetModel();
        ModelPart &r_hole_model_part = current_model.CreateModelPart("HoleModelpart");
        ModelPart &r_hole_boundary_model_part = current_model.CreateModelPart("HoleBoundaryModelPart");
        ModelPart &r_modified_patch_boundary_model_part = current_model.CreateModelPart("ModifiedPatchBoundary");
        ModelPart &r_modified_patch_model_part = current_model.CreateModelPart("ModifiedPatch");
        bool has_overlap = BoundingBoxTest(r_background_model_part, r_patch_model_part); // true if they dont overlap

        if (has_overlap)
        {
            KRATOS_INFO("Bounding boxes overlap , So finding the modified patch boundary") << std::endl;
            this->mpCalculateDistanceProcess->CalculateSignedDistance(r_patch_model_part, r_domain_boundary_model_part);
            //TODO: Below is brutforce. Check if the boundary of bg is actually cutting the patch.
            this->mpHoleCuttingUtility->RemoveOutOfDomainElements(r_patch_model_part, r_modified_patch_model_part, MainDomainOrNot);
        }

        mpHoleCuttingUtility->FindOutsideBoundaryOfModelPartGivenInside(r_modified_patch_model_part, r_patch_inside_boundary_model_part, r_modified_patch_boundary_model_part);
        this->mpCalculateDistanceProcess->CalculateSignedDistance(r_background_model_part, r_modified_patch_boundary_model_part);
        this->mpHoleCuttingUtility->CreateHoleAfterDistance(r_background_model_part, r_hole_model_part, r_hole_boundary_model_part, mOverlapDistance);

        //for multipatch
        const unsigned int n_elements = r_hole_model_part.NumberOfElements();
        for (unsigned int i_elem = 0; i_elem < n_elements; ++i_elem)
        {
            ModelPart::ElementsContainerType::iterator it_elem = r_hole_model_part.ElementsBegin() + i_elem;
            it_elem->Set(VISITED, true);
        }

        ApplyContinuityWithMpcs(r_modified_patch_boundary_model_part, p_point_locator_on_background);
        KRATOS_INFO("") << std::endl;
        KRATOS_INFO("") << std::endl;
        ApplyContinuityWithMpcs(r_hole_boundary_model_part, p_pointer_locator_on_patch);

        const unsigned int n_nodes = r_hole_model_part.NumberOfNodes();
        for (unsigned int i_node = 0; i_node < n_nodes; ++i_node)
        {
            auto it_node = r_hole_model_part.NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(VELOCITY_X) = 0.0;
            it_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
            it_node->FastGetSolutionStepValue(PRESSURE) = 0.0;
        }

        current_model.DeleteModelPart("HoleModelpart");
        current_model.DeleteModelPart("HoleBoundaryModelPart");
        current_model.DeleteModelPart("ModifiedPatchBoundary");
        current_model.DeleteModelPart("ModifiedPatch");
    }
    KRATOS_INFO("End of Formulate Chimera") << std::endl;
}

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::CreateConstraintIds(std::vector<int> &rIdVector, const IndexType NumberOfConstraintsRequired)
{
    IndexType max_constraint_id = 0;
    // Get current maximum constraint ID
    if (mrMainModelPart.MasterSlaveConstraints().size() != 0)
    {
        mrMainModelPart.MasterSlaveConstraints().Sort();
        ModelPart::MasterSlaveConstraintContainerType::iterator it = mrMainModelPart.MasterSlaveConstraintsEnd() - 1;
        max_constraint_id = (*it).Id();
        ++max_constraint_id;
    }

    // Now create a vector size NumberOfConstraintsRequired
    rIdVector.resize(NumberOfConstraintsRequired * (TDim + 1));
    std::iota(std::begin(rIdVector), std::end(rIdVector), max_constraint_id); // Fill with consecutive integers
}

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorPointerType &pBinLocator)
{
    //loop over nodes and find the triangle in which it falls, then do interpolation
    MasterSlaveContainerVectorType master_slave_container_vector;
#pragma omp parallel
    {
#pragma omp single
        {
            master_slave_container_vector.resize(omp_get_num_threads());
            for (auto &container : master_slave_container_vector)
                container.reserve(1000);
        }
    }
    std::vector<int> constraints_id_vector;

    int num_constraints_required = (TDim + 1) * (rBoundaryModelPart.Nodes().size());
    CreateConstraintIds(constraints_id_vector, num_constraints_required);

    const int max_results = 10000;
    const unsigned int n_boundary_nodes = rBoundaryModelPart.Nodes().size();
    std::size_t counter = 0;
    std::size_t removed_counter = 0;
    std::size_t not_found_counter = 0;

    for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
    {
        ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
        Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());

        mNodeIdToConstraintIdsMap[p_boundary_node->Id()].reserve(150);
    }

#pragma omp parallel for shared(constraints_id_vector, master_slave_container_vector, pBinLocator) reduction(+:not_found_counter) reduction(+:removed_counter) reduction(+ \
                                                                                                                                           : counter)
    for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
    {

        Vector shape_fun_weights;
        typename PointLocatorType::ResultContainerType results(max_results);
        auto &ms_container = master_slave_container_vector[omp_get_thread_num()]; //TODO: change this to out of loop.

        ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
        Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());
        ConstraintIdsVectorType constrainIds_for_the_node;
        unsigned int start_constraint_id = i_bn * (TDim + 1) * (TDim + 1);
        bool node_coupled = false;
        if ((p_boundary_node)->IsDefined(VISITED))
            node_coupled = (p_boundary_node)->Is(VISITED);

        typename PointLocatorType::ResultIteratorType result_begin = results.begin();
        Element::Pointer p_element;
        bool is_found = false;
        is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), shape_fun_weights, p_element, result_begin, max_results);

        if (node_coupled && is_found)
        {
            constrainIds_for_the_node = mNodeIdToConstraintIdsMap[p_boundary_node->Id()];
            for (auto const &constraint_id : constrainIds_for_the_node)
            {
                mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
                removed_counter++;
            }
            p_boundary_node->Set(VISITED, false);
        }

        if (is_found == true)
        {
            Geometry<Node<3>> &r_geom = p_element->GetGeometry();
            ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, start_constraint_id, constraints_id_vector, ms_container);
            counter += 1;
        }else{
            not_found_counter+=1;
        }
        p_boundary_node->Set(VISITED, true);
    } // end of loop over boundary nodes

    for (auto &container : master_slave_container_vector)
        mrMainModelPart.AddMasterSlaveConstraints(container.begin(), container.end());

    KRATOS_INFO("Number of boundary nodes in : ") << rBoundaryModelPart.Name() << " is coupled " << rBoundaryModelPart.NumberOfNodes() << std::endl;
    KRATOS_INFO("Number of Boundary nodes found : ") << counter<<". Number of constraints : "<<counter*9<< std::endl;
    KRATOS_INFO("Number of Boundary nodes not found  : ") << not_found_counter << std::endl;
}

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::GetBoundingBox(ModelPart &rModelPart, std::vector<double> &rLowPoint, std::vector<double> &rHighPoint)
{
    double rLowPoint0 = 1e10;
    double rLowPoint1 = 1e10;
    double rLowPoint2 = 1e10;

    double rHighPoint0 = -1e10;
    double rHighPoint1 = -1e10;
    double rHighPoint2 = -1e10;

    const unsigned int num_nodes = rModelPart.Nodes().size();
#pragma omp parallel for reduction(min                                                                                                                           \
                                   : rLowPoint0) reduction(min                                                                                                   \
                                                           : rLowPoint1) reduction(min                                                                           \
                                                                                   : rLowPoint2) reduction(max                                                   \
                                                                                                           : rHighPoint0) reduction(max                          \
                                                                                                                                    : rHighPoint1) reduction(max \
                                                                                                                                                             : rHighPoint2)
    for (unsigned int i_node = 0; i_node < num_nodes; ++i_node)
    {
        ModelPart::NodesContainerType::iterator it_node = rModelPart.NodesBegin() + i_node;
        rLowPoint0 = std::min(it_node->X(), rLowPoint0);
        rLowPoint1 = std::min(it_node->Y(), rLowPoint1);
        rLowPoint2 = std::min(it_node->Z(), rLowPoint2);

        rHighPoint0 = std::max(it_node->X(), rHighPoint0);
        rHighPoint1 = std::max(it_node->Y(), rHighPoint1);
        rHighPoint2 = std::max(it_node->Z(), rHighPoint2);
    }

    rHighPoint[0] = rHighPoint0;
    rHighPoint[1] = rHighPoint1;
    rHighPoint[2] = rHighPoint2;

    rLowPoint[0] = rLowPoint0;
    rLowPoint[1] = rLowPoint1;
    rLowPoint[2] = rLowPoint2;
}

template <int TDim>
bool ApplyChimeraProcessMonolithic<TDim>::BoundingBoxTest(ModelPart &rModelPartA, ModelPart &rModelPartB) //background A and Patch B
{
    std::vector<double> min_cornerA(3), max_cornerA(3), min_cornerB(3), max_cornerB(3);
    GetBoundingBox(rModelPartA, min_cornerA, max_cornerA);
    GetBoundingBox(rModelPartB, min_cornerB, max_cornerB);
    const int dim = rModelPartA.GetProcessInfo().GetValue(DOMAIN_SIZE);

    KRATOS_INFO("Bounding box of Background") << min_cornerA[0] << "::" << max_cornerA[0] << std::endl;
    KRATOS_INFO("Bounding box of patch") << min_cornerB[0] << "::" << max_cornerB[0] << std::endl;

    for (int i = 0; i < dim; i++)
    {
        if (min_cornerA[i] > max_cornerB[i])
            return false;
        if (max_cornerA[i] < min_cornerB[i])
            return false;
    }
    return true;
}

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::ApplyContinuityWithElement(Geometry<Node<3>> &rGeometry,
                                                                     Node<3> &rBoundaryNode,
                                                                     Vector &rShapeFuncWeights,
                                                                     unsigned int StartId,
                                                                     std::vector<int> &ConstraintIdVector,
                                                                     MasterSlaveConstraintContainerType &rMsContainer)
{
    const auto &r_clone_constraint = (LinearMasterSlaveConstraint)KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
    // Initialise the boundary nodes dofs to 0 at ever time steps
    rBoundaryNode.FastGetSolutionStepValue(VELOCITY_X, 0) = 0.0;
    rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Y, 0) = 0.0;
    if (TDim == 3)
        rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Z, 0) = 0.0;
    rBoundaryNode.FastGetSolutionStepValue(PRESSURE, 0) = 0.0;
    for (std::size_t i = 0; i < rGeometry.size(); i++)
    {
        //Interpolation of velocity
        rBoundaryNode.FastGetSolutionStepValue(VELOCITY_X, 0) += rGeometry[i].GetDof(VELOCITY_X).GetSolutionStepValue(0) * rShapeFuncWeights[i];
        rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Y, 0) += rGeometry[i].GetDof(VELOCITY_Y).GetSolutionStepValue(0) * rShapeFuncWeights[i];
        //Interpolation of pressure
        rBoundaryNode.FastGetSolutionStepValue(PRESSURE, 0) += rGeometry[i].GetDof(PRESSURE).GetSolutionStepValue(0) * rShapeFuncWeights[i];

        //Define master slave relation for velocity X and Y
        AddMasterSlaveRelation(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], VELOCITY_X, rBoundaryNode, VELOCITY_X, rShapeFuncWeights[i]);
        AddMasterSlaveRelation(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], VELOCITY_Y, rBoundaryNode, VELOCITY_Y, rShapeFuncWeights[i]);
        if (TDim == 3)
        {
            //Interpolation of velocity Z
            rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Z, 0) += rGeometry[i].GetDof(VELOCITY_Z).GetSolutionStepValue(0) * rShapeFuncWeights[i];
            AddMasterSlaveRelation(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], VELOCITY_Z, rBoundaryNode, VELOCITY_Z, rShapeFuncWeights[i]);
        }
        //Defining master slave relation for pressure
        AddMasterSlaveRelation(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], PRESSURE, rBoundaryNode, PRESSURE, rShapeFuncWeights[i]);
    } // end of loop over host element nodes

    // Setting the buffer 1 same buffer 0
    rBoundaryNode.FastGetSolutionStepValue(VELOCITY_X, 1) = rBoundaryNode.FastGetSolutionStepValue(VELOCITY_X, 0);
    rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Y, 1) = rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Y, 0);
    if (TDim == 3)
        rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Z, 1) = rBoundaryNode.FastGetSolutionStepValue(VELOCITY_Z, 0);
    rBoundaryNode.FastGetSolutionStepValue(PRESSURE, 1) = rBoundaryNode.FastGetSolutionStepValue(PRESSURE, 0);
}

template class ApplyChimeraProcessMonolithic<2>;
template class ApplyChimeraProcessMonolithic<3>;

} // namespace Kratos.

