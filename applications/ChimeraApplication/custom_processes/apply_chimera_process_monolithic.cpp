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
    this->mpHoleCuttingProcess = CustomHoleCuttingProcess::Pointer(new CustomHoleCuttingProcess());
    this->mpCalculateDistanceProcess = typename CustomCalculateSignedDistanceProcess<TDim>::Pointer(new CustomCalculateSignedDistanceProcess<TDim>());

    mNumberOfConstraintsAdded = 0;
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
    int number_of_constraints = mrMainModelPart.MasterSlaveConstraints().size();
    KRATOS_INFO("Total number of constraints created so far") << number_of_constraints << std::endl;
}

//Apply Chimera with or without overlap

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::FormulateChimera(int MainDomainOrNot)
{
    ModelPart &rBackgroundModelPart = mrMainModelPart.GetSubModelPart(m_background_model_part_name);
    ModelPart &rPatchModelPart = mrMainModelPart.GetSubModelPart(m_patch_model_part_name);
    ModelPart &rDomainBoundaryModelPart = mrMainModelPart.GetSubModelPart(m_domain_boundary_model_part_name);
    ModelPart &rPatchInsideBoundaryModelPart = mrMainModelPart.GetSubModelPart(m_patch_inside_boundary_model_part_name);

    this->mpBinLocatorForBackground = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(rBackgroundModelPart));
    this->mpBinLocatorForPatch = BinBasedPointLocatorPointerType(new BinBasedFastPointLocator<TDim>(rPatchModelPart));
    this->mpBinLocatorForBackground->UpdateSearchDatabase();
    this->mpBinLocatorForPatch->UpdateSearchDatabase();

    const double eps = 1e-12;
    if (mOverlapDistance < eps)
        KRATOS_THROW_ERROR("", "Overlap distance should be a positive and non-zero number \n", "");

    if (mOverlapDistance > eps)
    {
        Model &current_model = mrMainModelPart.GetModel();
        ModelPart &pHoleModelPart = current_model.CreateModelPart("HoleModelpart");
        ModelPart &pHoleBoundaryModelPart = current_model.CreateModelPart("HoleBoundaryModelPart");
        ModelPart &pModifiedPatchBoundaryModelPart = current_model.CreateModelPart("ModifiedPatchBoundary");
        ModelPart &pModifiedPatchModelPart = current_model.CreateModelPart("ModifiedPatch");
        bool has_overlap = BoundingBoxTest(rBackgroundModelPart, rPatchModelPart); // true if they dont overlap

        if (has_overlap)
        {
            KRATOS_INFO("Bounding boxes overlap , So finding the modified patch boundary") << std::endl;
            this->mpCalculateDistanceProcess->CalculateSignedDistance(rPatchModelPart, rDomainBoundaryModelPart);
            this->mpHoleCuttingProcess->RemoveOutOfDomainPatchAndReturnModifiedPatch(rPatchModelPart, rPatchInsideBoundaryModelPart, pModifiedPatchModelPart, pModifiedPatchBoundaryModelPart, MainDomainOrNot);
        }
        else
        {
            KRATOS_INFO("Bounding boxes does NOT overlap , So finding outside boundary of patch using the inside boundary") << std::endl;
            FindOutsideBoundaryOfModelPartGivenInside(rPatchModelPart, rPatchInsideBoundaryModelPart, pModifiedPatchBoundaryModelPart);
        }

        this->mpCalculateDistanceProcess->CalculateSignedDistance(rBackgroundModelPart, pModifiedPatchBoundaryModelPart);
        this->mpHoleCuttingProcess->CreateHoleAfterDistance(rBackgroundModelPart, pHoleModelPart, pHoleBoundaryModelPart, mOverlapDistance);

        //for multipatch
        for (ModelPart::ElementsContainerType::iterator it = pHoleModelPart.ElementsBegin(); it != pHoleModelPart.ElementsEnd(); ++it)
            it->Set(VISITED, true);

        KRATOS_INFO("Formulate Chimera: Number of nodes in modified patch boundary : ") << pModifiedPatchBoundaryModelPart.Nodes().size() << std::endl;
        KRATOS_INFO("Formulate Chimera: Number of nodes in hole boundary : ") << pHoleBoundaryModelPart.Nodes().size() << std::endl;
        mNumberOfConstraintsAdded = 0;
        ApplyMpcConstraint(pModifiedPatchBoundaryModelPart, mpBinLocatorForBackground);
        KRATOS_INFO("Formulate Chimera: Constraints formulated for modified patch boundary ... ") << std::endl;

        mNumberOfConstraintsAdded = 0;
        ApplyMpcConstraint(pHoleBoundaryModelPart, mpBinLocatorForPatch);
        KRATOS_INFO("Formulate Chimera: Constraints formulated for hole boundary ... ") << std::endl;

        KRATOS_INFO("Patch boundary coupled with background & HoleBoundary  coupled with patch using nearest element approach") << std::endl;
        KRATOS_INFO("Formulate Chimera: Appplied MPCs") << std::endl;

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
void ApplyChimeraProcessMonolithic<TDim>::SetOverlapDistance(double distance)
{
    this->mOverlapDistance = distance;
}

template <int TDim>
void ApplyChimeraProcessMonolithic<TDim>::ApplyMpcConstraint(ModelPart &rBoundaryModelPart, BinBasedPointLocatorPointerType &pBinLocator)
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
    const auto &r_clone_constraint = (LinearMasterSlaveConstraint)KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
    {
        ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
        Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());

        mNodeIdToConstraintIdsMap[p_boundary_node->Id()].reserve(150);
    }

#pragma omp parallel for shared(constraints_id_vector, master_slave_container_vector, pBinLocator) firstprivate(removed_counter, r_clone_constraint) reduction(+ \
                                                                                                                                                               : counter)
    for (unsigned int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
    {

        Vector N;
        typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
        auto &ms_container = master_slave_container_vector[omp_get_thread_num()];

        ModelPart::NodesContainerType::iterator i_boundary_node = rBoundaryModelPart.NodesBegin() + i_bn;
        Node<3>::Pointer p_boundary_node = *(i_boundary_node.base());
        ConstraintIdsVectorType ConstrainIdsForTheNode;
        unsigned int start_constraint_id = i_bn * (TDim + 1) * (TDim + 1);
        bool NodeCoupled = false;
        if ((p_boundary_node)->IsDefined(VISITED))
            NodeCoupled = (p_boundary_node)->Is(VISITED);

        typename BinBasedFastPointLocator<TDim>::ResultIteratorType result_begin = results.begin();
        Element::Pointer pElement;
        bool is_found = false;
        is_found = pBinLocator->FindPointOnMesh(p_boundary_node->Coordinates(), N, pElement, result_begin, max_results);

        if (NodeCoupled && is_found)
        {
            ConstrainIdsForTheNode = mNodeIdToConstraintIdsMap[p_boundary_node->Id()];
            for (auto const &constraint_id : ConstrainIdsForTheNode)
            {
                mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
                removed_counter++;
            }
            p_boundary_node->Set(VISITED, false);
        }

        // Initialise the boundary nodes dofs to 0 at ever time steps
        p_boundary_node->FastGetSolutionStepValue(VELOCITY_X, 0) = 0.0;
        p_boundary_node->FastGetSolutionStepValue(VELOCITY_Y, 0) = 0.0;
        if (TDim == 3)
            p_boundary_node->FastGetSolutionStepValue(VELOCITY_Z, 0) = 0.0;
        p_boundary_node->FastGetSolutionStepValue(PRESSURE, 0) = 0.0;

        if (is_found == true)
        {
            Geometry<Node<3>> &geom = pElement->GetGeometry();
            for (std::size_t i = 0; i < geom.size(); i++)
            {

                //Interpolation of velocity
                p_boundary_node->FastGetSolutionStepValue(VELOCITY_X, 0) += geom[i].GetDof(VELOCITY_X).GetSolutionStepValue(0) * N[i];
                p_boundary_node->FastGetSolutionStepValue(VELOCITY_Y, 0) += geom[i].GetDof(VELOCITY_Y).GetSolutionStepValue(0) * N[i];
                //Interpolation of pressure
                p_boundary_node->FastGetSolutionStepValue(PRESSURE, 0) += geom[i].GetDof(PRESSURE).GetSolutionStepValue(0) * N[i];

                //Define master slave relation for velocity X and Y
                AddMasterSlaveRelation(ms_container, r_clone_constraint, constraints_id_vector[start_constraint_id++], geom[i], VELOCITY_X, *p_boundary_node, VELOCITY_X, N[i]);
                AddMasterSlaveRelation(ms_container, r_clone_constraint, constraints_id_vector[start_constraint_id++], geom[i], VELOCITY_Y, *p_boundary_node, VELOCITY_Y, N[i]);
                if (TDim == 3)
                {
                    //Interpolation of velocity Z
                    p_boundary_node->FastGetSolutionStepValue(VELOCITY_Z, 0) += geom[i].GetDof(VELOCITY_Z).GetSolutionStepValue(0) * N[i];
                    AddMasterSlaveRelation(ms_container, r_clone_constraint, constraints_id_vector[start_constraint_id++], geom[i], VELOCITY_Z, *p_boundary_node, VELOCITY_Z, N[i]);
                }
                //Defining master slave relation for pressure
                AddMasterSlaveRelation(ms_container, r_clone_constraint, constraints_id_vector[start_constraint_id++], geom[i], PRESSURE, *p_boundary_node, PRESSURE, N[i]);

                counter++;

            } // end of loop over host element nodes

            // Setting the buffer 1 same buffer 0
            p_boundary_node->FastGetSolutionStepValue(VELOCITY_X, 1) = p_boundary_node->FastGetSolutionStepValue(VELOCITY_X, 0);
            p_boundary_node->FastGetSolutionStepValue(VELOCITY_Y, 1) = p_boundary_node->FastGetSolutionStepValue(VELOCITY_Y, 0);
            if (TDim == 3)
                p_boundary_node->FastGetSolutionStepValue(VELOCITY_Z, 1) = p_boundary_node->FastGetSolutionStepValue(VELOCITY_Z, 0);
            p_boundary_node->FastGetSolutionStepValue(PRESSURE, 1) = p_boundary_node->FastGetSolutionStepValue(PRESSURE, 0);
        }
        p_boundary_node->Set(VISITED, true);
    } // end of loop over boundary nodes

    for (auto &container : master_slave_container_vector)
        mrMainModelPart.AddMasterSlaveConstraints(container.begin(), container.end());

    counter /= TDim + 1;
    KRATOS_INFO("Pressure nodes from") << rBoundaryModelPart.Name() << " is coupled" << counter << std::endl;
    KRATOS_INFO("number of constraints created for this combination") << counter * 9 << std::endl;
    KRATOS_INFO("number of constraints removed for this combination") << removed_counter << std::endl;
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
void ApplyChimeraProcessMonolithic<TDim>::FindOutsideBoundaryOfModelPartGivenInside(ModelPart &rModelPart, ModelPart &rInsideBoundary, ModelPart &rExtractedBoundaryModelPart)
{
    std::size_t n_nodes = rModelPart.ElementsBegin()->GetGeometry().size();

    if (n_nodes == 3)
        this->mpHoleCuttingProcess->ExtractOutsideBoundaryMesh(rInsideBoundary, rModelPart, rExtractedBoundaryModelPart);
    else if (n_nodes == 4)
        this->mpHoleCuttingProcess->ExtractOutsideSurfaceMesh(rInsideBoundary, rModelPart, rExtractedBoundaryModelPart);
}

template class ApplyChimeraProcessMonolithic<2>;
template class ApplyChimeraProcessMonolithic<3>;

} // namespace Kratos.
