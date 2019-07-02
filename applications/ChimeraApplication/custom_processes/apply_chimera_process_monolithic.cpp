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

template <int TDim, class TDistanceCalculatorType>
ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::ApplyChimeraProcessMonolithic(ModelPart &rMainModelPart, Parameters iParameters) : Process(), mrMainModelPart(rMainModelPart), mParameters(iParameters)
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
    mpHoleCuttingUtility = ChimeraHoleCuttingUtility::Pointer(new ChimeraHoleCuttingUtility());
    mpCalculateDistanceProcess = typename DistanceCalculatorType::Pointer(new DistanceCalculatorType());
}

template <int TDim, class TDistanceCalculatorType>
ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::~ApplyChimeraProcessMonolithic()
{
}

template <int TDim, class TDistanceCalculatorType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY;
    // Actual execution of the functionality of this class
    for (ModelPart::ElementsContainerType::iterator it = mrMainModelPart.ElementsBegin(); it != mrMainModelPart.ElementsEnd(); ++it)
        it->SetValue(SPLIT_ELEMENT, false);
    DoChimeraLoop();
    KRATOS_CATCH("");
}

template <int TDim, class TDistanceCalculatorType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::ExecuteFinalizeSolutionStep()
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

template <int TDim, class TDistanceCalculatorType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::DoChimeraLoop() //selecting patch and background combination for chimera method
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
template <int TDim, class TDistanceCalculatorType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::FormulateChimera(int MainDomainOrNot)
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
            mpCalculateDistanceProcess->CalculateSignedDistance(r_patch_model_part, r_domain_boundary_model_part);
            //TODO: Below is brutforce. Check if the boundary of bg is actually cutting the patch.
            mpHoleCuttingUtility->RemoveOutOfDomainElements(r_patch_model_part, r_modified_patch_model_part, MainDomainOrNot, false);
        }

        mpHoleCuttingUtility->FindOutsideBoundaryOfModelPartGivenInside(r_modified_patch_model_part, r_patch_inside_boundary_model_part, r_modified_patch_boundary_model_part);
        mpCalculateDistanceProcess->CalculateSignedDistance(r_background_model_part, r_modified_patch_boundary_model_part);
        mpHoleCuttingUtility->CreateHoleAfterDistance(r_background_model_part, r_hole_model_part, r_hole_boundary_model_part, mOverlapDistance);

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

        current_model.DeleteModelPart("HoleModelpart");
        current_model.DeleteModelPart("HoleBoundaryModelPart");
        current_model.DeleteModelPart("ModifiedPatchBoundary");
        current_model.DeleteModelPart("ModifiedPatch");
    }
    KRATOS_INFO("End of Formulate Chimera") << std::endl;
}

template <int TDim, class TDistanceCalculatorType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::CreateConstraintIds(std::vector<int> &rIdVector, const IndexType NumberOfConstraintsRequired)
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

template <int TDim, class TDistanceCalculatorType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorPointerType &pBinLocator)
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
            int init_index = 0;
            ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_X, start_constraint_id+init_index, constraints_id_vector, ms_container);
            init_index+=3;
            ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_Y, start_constraint_id+init_index, constraints_id_vector, ms_container);
            init_index+=3;
            if(TDim == 3){
                ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, VELOCITY_Z, start_constraint_id+init_index, constraints_id_vector, ms_container);
                init_index+=3;
            }
            ApplyContinuityWithElement(r_geom, *p_boundary_node, shape_fun_weights, PRESSURE, start_constraint_id+init_index, constraints_id_vector, ms_container);
            init_index+=3;
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

template <int TDim, class TDistanceCalculatorType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::GetBoundingBox(ModelPart &rModelPart, std::vector<double> &rLowPoint, std::vector<double> &rHighPoint)
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

template <int TDim, class TDistanceCalculatorType>
bool ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::BoundingBoxTest(ModelPart &rModelPartA, ModelPart &rModelPartB) //background A and Patch B
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

template <int TDim, class TDistanceCalculatorType>
template <typename TVariableType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::ApplyContinuityWithElement(Geometry<Node<3>> &rGeometry,
                                                                     Node<3> &rBoundaryNode,
                                                                     Vector &rShapeFuncWeights,
                                                                     TVariableType& rVariable,
                                                                     unsigned int StartId,
                                                                     std::vector<int> &ConstraintIdVector,
                                                                     MasterSlaveConstraintContainerType &rMsContainer)
{
    const auto &r_clone_constraint = (LinearMasterSlaveConstraint)KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
    // Initialise the boundary nodes dofs to 0 at ever time steps
    rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) = 0.0;
    for (std::size_t i = 0; i < rGeometry.size(); i++)
    {
        //Interpolation of rVariable
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) += rGeometry[i].GetDof(rVariable).GetSolutionStepValue(0) * rShapeFuncWeights[i];
        AddMasterSlaveRelation(rMsContainer, r_clone_constraint, ConstraintIdVector[StartId++], rGeometry[i], rVariable, rBoundaryNode, rVariable, rShapeFuncWeights[i]);
    } // end of loop over host element nodes

    // Setting the buffer 1 same buffer 0
    rBoundaryNode.FastGetSolutionStepValue(rVariable, 1) = rBoundaryNode.FastGetSolutionStepValue(rVariable, 0);
}



template <int TDim, class TDistanceCalculatorType>
template <typename TVariableType>
void ApplyChimeraProcessMonolithic<TDim, TDistanceCalculatorType>::AddMasterSlaveRelation(MasterSlaveConstraintContainerType &rMasterSlaveContainer,
                                    const LinearMasterSlaveConstraint &rCloneConstraint,
                                    unsigned int ConstraintId,
                                    Node<3> &rMasterNode,
                                    TVariableType &rMasterVariable,
                                    Node<3> &rSlaveNode,
                                    TVariableType &rSlaveVariable,
                                    const double Weight,
                                    const double Constant)
{
    rSlaveNode.Set(SLAVE);
    ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = rCloneConstraint.Create(ConstraintId, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
    p_new_constraint->Set(TO_ERASE);
    mNodeIdToConstraintIdsMap[rSlaveNode.Id()].push_back(ConstraintId);
    rMasterSlaveContainer.insert(rMasterSlaveContainer.begin(), p_new_constraint);
}

typedef CustomCalculateSignedDistanceProcess<2> DistanceCalculator2DType;
typedef CustomCalculateSignedDistanceProcess<3> DistanceCalculator3DType;

template class ApplyChimeraProcessMonolithic<2, DistanceCalculator2DType>;
template class ApplyChimeraProcessMonolithic<3, DistanceCalculator3DType>;

} // namespace Kratos.

