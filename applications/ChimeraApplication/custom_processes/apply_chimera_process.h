//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Authors:        Aditya Ghantasala, https://github.com/adityaghantasala
// 					Navaneeth K Narayanan
//					Rishith Ellath Meethal
// ==============================================================================
//
#if !defined(KRATOS_APPLY_CHIMERA_H_INCLUDED)
#define KRATOS_APPLY_CHIMERA_H_INCLUDED

// System includes
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include "omp.h"

// External includes

// Project includes
#include "includes/define.h"
#include "includes/process_info.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "includes/linear_master_slave_constraint.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "input_output/vtk_output.h"
#include "utilities/binbased_fast_point_locator.h"
#include "factories/standard_linear_solver_factory.h"
#include "utilities/builtin_timer.h"
#include "utilities/variable_utils.h"
#ifdef KRATOS_USING_MPI
#include "mpi/utilities/parallel_fill_communicator.h"
#endif

// Application includes
#include "chimera_application_variables.h"
#include "custom_utilities/hole_cutting_utility.h"
#include "custom_utilities/distance_calcuation_utility.h"
#include "custom_utilities/gather_modelpart_on_all_ranks.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <int TDim>
class KRATOS_API(CHIMERA_APPLICATION) ApplyChimera : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    typedef ProcessInfo::Pointer ProcessInfoPointerType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    typedef std::size_t IndexType;
    typedef ModelPart::NodeType NodeType;
    typedef Kratos::Variable<double> VariableType;
    typedef std::vector<IndexType> ConstraintIdsVectorType;
    typedef typename ModelPart::MasterSlaveConstraintType MasterSlaveConstraintType;
    typedef typename ModelPart::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef std::vector<MasterSlaveConstraintContainerType> MasterSlaveContainerVectorType;
    typedef BinBasedFastPointLocator<TDim> PointLocatorType;
    typedef typename PointLocatorType::Pointer PointLocatorPointerType;
    typedef ModelPart::NodesContainerType NodesContainerType;

    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimera);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor
     * @param rMainModelPart The reference to the modelpart which will be used for computations later on.
     * @param iParameters The settings parameters.
     */
    explicit ApplyChimera(ModelPart &rMainModelPart, Parameters iParameters) : mrMainModelPart(rMainModelPart), mParameters(iParameters)
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
        KRATOS_ERROR_IF(mNumberOfLevels < 2) << "Chimera requires atleast two levels !. Put Background in one level and the patch in second one." << std::endl;

        ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
        mEchoLevel = 0;
        mReformulateEveryStep = false;
        mIsFormulated = false;
    }

    /// Destructor.
    virtual ~ApplyChimera()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetEchoLevel(int EchoLevel)
    {
        mEchoLevel = EchoLevel;
    }

    void SetReformulateEveryStep(bool Reformulate)
    {
        mReformulateEveryStep = Reformulate;
    }

    virtual void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        VariableUtils().SetFlag(ACTIVE, true, mrMainModelPart.Elements());
        // Actual execution of the functionality of this class
        if (mReformulateEveryStep || !mIsFormulated)
        {
            const auto &r_comm = mrMainModelPart.GetCommunicator();
            ModelPart::MeshType::Pointer p_interface_mesh = r_comm.pGhostMesh();
            if(r_comm.IsDistributed()){
                const int mpi_rank = r_comm.MyPID();
                for(const auto& node_id : mRemoteNodes)
                {
                    if(mrMainModelPart.GetNode(node_id).FastGetSolutionStepValue(PARTITION_INDEX) != mpi_rank)
                        p_interface_mesh->RemoveNode(node_id);
                }
                mRemoteNodes.clear();
            }
            DoChimeraLoop();
            mIsFormulated = true;
        }
        KRATOS_CATCH("");
    }

    virtual void ExecuteFinalizeSolutionStep() override
    {
        VariableUtils().SetFlag(VISITED, false, mrMainModelPart.Nodes());
        VariableUtils().SetFlag(VISITED, false, mrMainModelPart.Elements());
        VariableUtils().SetNonHistoricalVariable(SPLIT_ELEMENT, false, mrMainModelPart.Elements());

        if (mReformulateEveryStep)
        {
            mrMainModelPart.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
            mIsFormulated = false;
        }
    }

    virtual std::string Info() const override
    {
        return "ApplyChimera";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "ApplyChimera" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {
        rOStream << "ApplyChimera" << std::endl;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    ModelPart &mrMainModelPart;
    IndexType mNumberOfLevels;
    Parameters mParameters;
    std::unordered_map<IndexType, ConstraintIdsVectorType> mNodeIdToConstraintIdsMap;
    int mEchoLevel;
    bool mReformulateEveryStep;
    std::map<std::string, PointLocatorPointerType> mPointLocatorsMap;
    bool mIsFormulated;
    std::vector<IndexType> mRemoteNodes;

    // Modelpart names which are generated here
    std::string mModifiedName = "ChimeraModified";
    std::string mBoundaryName = "ChimeraBoundary";
    std::string mHoleName = "ChimeraHole";
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Does a loop on the background and patch combinations possible and uses FormulateChimera method.
     */
    virtual void DoChimeraLoop()
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
        for (int i_be = 0; i_be < num_elements; ++i_be)
        {
            auto i_elem = elem_begin + i_be;
            if (!i_elem->Is(VISITED)) //for multipatch
                i_elem->Set(ACTIVE, true);
        }

        const int num_nodes = static_cast<int>(mrMainModelPart.NumberOfNodes());
        const auto nodes_begin = mrMainModelPart.ElementsBegin();

#pragma omp parallel for
        for (int i_bn = 0; i_bn < num_nodes; ++i_bn)
        {
            auto i_node = nodes_begin + i_bn;
            i_node->Set(VISITED, false);
        }
        BuiltinTimer do_chimera_loop_time;

        int i_current_level = 0;
        for (auto &current_level : mParameters)
        {
            ChimeraHoleCuttingUtility::Domain domain_type = ChimeraHoleCuttingUtility::Domain::MAIN_BACKGROUND;
            for (auto &background_patch_param : current_level)
            { // Gives the current background
                background_patch_param.ValidateAndAssignDefaults(parameters_for_validation);
                Model &current_model = mrMainModelPart.GetModel();
                ModelPart &r_background_model_part = current_model.GetModelPart(background_patch_param["model_part_name"].GetString());
                // compute the outerboundary of the background to save
                if (i_current_level == 0)
                    if (!r_background_model_part.HasSubModelPart(mBoundaryName))
                    {
                        auto &r_boundary_model_part = r_background_model_part.CreateSubModelPart(mBoundaryName);
                        BuiltinTimer extraction_time;
                        ChimeraHoleCuttingUtility().ExtractBoundaryMesh<TDim>(r_background_model_part, r_boundary_model_part);
                        KRATOS_INFO_IF("ApplyChimera : Extraction of boundary mesh took          : ", mEchoLevel > 0) << extraction_time.ElapsedSeconds() << " seconds" << std::endl;
                    }

                for (IndexType i_slave_level = i_current_level + 1; i_slave_level < mNumberOfLevels; ++i_slave_level)
                {
                    for (auto &slave_patch_param : mParameters[i_slave_level]) // Loop over all other slave patches
                    {
                        slave_patch_param.ValidateAndAssignDefaults(parameters_for_validation);
                        KRATOS_INFO_IF("Formulating Chimera for the combination :: ", mEchoLevel > 0) << "\n :: Background :: \n" << background_patch_param << "\n :: Patch :: \n" << slave_patch_param << std::endl;
                        if (i_current_level == 0) // a check to identify computational Domain boundary
                            domain_type = ChimeraHoleCuttingUtility::Domain::OTHER;
                        FormulateChimera(background_patch_param, slave_patch_param, domain_type);
                        mPointLocatorsMap.erase(background_patch_param["model_part_name"].GetString());
                    }
                }
            }
            ++i_current_level;
        }

        mNodeIdToConstraintIdsMap.clear();
        mPointLocatorsMap.clear();

        const auto &r_comm = mrMainModelPart.GetCommunicator().GetDataCommunicator();
        double chimera_elased_seconds = do_chimera_loop_time.ElapsedSeconds();
        KRATOS_INFO_IF("ApplyChimera : Chrimera Initialization took              : ", mEchoLevel > 0)<< r_comm.Max(chimera_elased_seconds, 0) << " seconds" << std::endl;
        int num_local_constraints = mrMainModelPart.NumberOfMasterSlaveConstraints();

        KRATOS_INFO_IF("ApplyChimera : Total number of constraints created so far:", mEchoLevel > 0) << r_comm.Sum(num_local_constraints, 0) << std::endl;
    }

    /**
     * @brief Formulates the Chimera conditions with a given set of background and patch combination.
     * @param BackgroundParam Parameters/Settings for the background
     * @param PatchParameters Parameters/Settings for the Patch
     * @param DomainType Flag specifying if the background is the main bg or not
     */
    virtual void FormulateChimera(const Parameters BackgroundParam, const Parameters PatchParameters, ChimeraHoleCuttingUtility::Domain DomainType)
    {
        Model &current_model = mrMainModelPart.GetModel();
        const auto &r_comm = mrMainModelPart.GetCommunicator().GetDataCommunicator();
        ModelPart &r_background_model_part = current_model.GetModelPart(BackgroundParam["model_part_name"].GetString());
        ModelPart &r_background_boundary_model_part = r_background_model_part.GetSubModelPart(mBoundaryName);
        ModelPart &r_patch_model_part = current_model.GetModelPart(PatchParameters["model_part_name"].GetString());
        const std::string bg_search_mp_name = BackgroundParam["search_model_part_name"].GetString();
        auto &r_background_search_model_part = current_model.HasModelPart(bg_search_mp_name) ? current_model.GetModelPart(bg_search_mp_name) : r_background_model_part;

        const double overlap_bg = BackgroundParam["overlap_distance"].GetDouble();
        const double overlap_pt = PatchParameters["overlap_distance"].GetDouble();
        const double over_lap_distance = (overlap_bg > overlap_pt) ? overlap_bg : overlap_pt;

        BuiltinTimer search_creation_time;
        PointLocatorPointerType p_point_locator_on_background = GetPointLocator(r_background_search_model_part);
        PointLocatorPointerType p_pointer_locator_on_patch = GetPointLocator(r_patch_model_part);
        auto elapsed_search_creation_time = search_creation_time.ElapsedSeconds();
        KRATOS_INFO_IF("ApplyChimera : Creation of search structures took        : ", mEchoLevel > 0) << r_comm.Max(elapsed_search_creation_time, 0) << " seconds" << std::endl;
        KRATOS_ERROR_IF(over_lap_distance < 1e-12) << "Overlap distance should be a positive and non-zero number." << std::endl;

        ModelPart &r_hole_model_part = r_background_model_part.CreateSubModelPart(mHoleName);
        ModelPart &r_hole_boundary_model_part = r_hole_model_part.CreateSubModelPart(mBoundaryName);

        auto &r_patch_boundary_model_part = ExtractPatchBoundary(PatchParameters, r_background_boundary_model_part, DomainType);

        BuiltinTimer bg_distance_calc_time;
        DistanceCalculationUtility <TDim>::CalculateDistance(r_background_model_part, r_patch_boundary_model_part);
        KRATOS_INFO_IF("Distance calculation on background took                  : ", mEchoLevel > 0)<< bg_distance_calc_time.ElapsedSeconds()<< " seconds"<< std::endl;

        BuiltinTimer hole_creation_time;
        ChimeraHoleCuttingUtility().CreateHoleAfterDistance<TDim>(r_background_model_part, r_hole_model_part, r_hole_boundary_model_part, over_lap_distance);
        auto hole_creation_time_elapsed = hole_creation_time.ElapsedSeconds();
        KRATOS_INFO_IF("ApplyChimera : Hole creation took                        : ", mEchoLevel > 0) << r_comm.Max(hole_creation_time_elapsed, 0) << " seconds" << std::endl;

        // WriteModelPart(r_hole_model_part);
        // WriteModelPart(r_background_boundary_model_part);
        // WriteModelPart(r_background_model_part);
        // WriteModelPart(r_patch_boundary_model_part);
        // WriteModelPart(r_hole_boundary_model_part);
        // WriteModelPart(r_background_model_part);

        const int n_elements = static_cast<int>(r_hole_model_part.NumberOfElements());
#pragma omp parallel for
        for (int i_elem = 0; i_elem < n_elements; ++i_elem)
        {
            ModelPart::ElementsContainerType::iterator it_elem = r_hole_model_part.ElementsBegin() + i_elem;
            it_elem->Set(VISITED, true); //for multipatch
        }

        BuiltinTimer mpc_time;
        ApplyContinuityWithMpcs(r_patch_boundary_model_part, *p_point_locator_on_background);
        ApplyContinuityWithMpcs(r_hole_boundary_model_part, *p_pointer_locator_on_patch);
        auto mpc_time_elapsed = mpc_time.ElapsedSeconds();
        KRATOS_INFO_IF("ApplyChimera : Creation of MPC for chimera took          : ", mEchoLevel > 0) << r_comm.Max(mpc_time_elapsed, 0) << " seconds" << std::endl;

        r_hole_boundary_model_part.RemoveSubModelPart(mBoundaryName);
        r_background_model_part.RemoveSubModelPart(mHoleName);
        r_patch_model_part.RemoveSubModelPart(mModifiedName);

#ifdef KRATOS_USING_MPI
        BuiltinTimer par_fill_comm;
        ParallelFillCommunicator(mrMainModelPart).Execute();
        double par_fill_time = par_fill_comm.ElapsedSeconds();
        KRATOS_INFO_IF("SynchronizeNodes : Time taken for parallel fill comm     : ", mEchoLevel > 1) << r_comm.Max(par_fill_time, 0) << std::endl;
#endif

        KRATOS_INFO("End of Formulate Chimera") << std::endl;
    }

    /**
     * @brief Creates a vector of unique constraint ids based on how many required and how many are already present in the mrModelPart.
     * @param rIdVector The vector which is populated with unique constraint ids.
     * @param NumberOfConstraintsRequired The number of further constraints required. used for calculation of unique ids.
     */
    virtual void CreateConstraintIds(std::vector<int> &rIdVector, const IndexType NumberOfConstraintsRequired)
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

    /**
     * @brief Applies the continuity between the boundary modelpart and the background.
     * @param rBoundaryModelPart The boundary modelpart for which the continuity is to be enforced.
     * @param pBinLocator The bin based locator formulated on the background. This is used to locate nodes of rBoundaryModelPart on background.
     */
    virtual void ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorType &pBinLocator)
    {
    }

    /**
     * @brief Applies the master-slave constraint between the given master and slave nodes with corresponding variable.
     * @param rMasterSlaveContainer The container to which the constraint to be added (useful to so OpenMP loop)
     * @param rCloneConstraint The prototype of constraint which is to be added.
     * @param ConstraintId The ID of the constraint to be added.
     * @param rMasterNode The Master node of the constraint.
     * @param rMasterVariable The variable for the master node.
     * @param rSlaveNode The Slave node of the constraint.
     * @param rSlaveVariable The variable for the slave node.
     * @param Weight The weight of the Master node.
     * @param Constant The constant of the master slave relation.
     */
    template <typename TVariableType>
    void AddMasterSlaveRelation(MasterSlaveConstraintContainerType &rMasterSlaveContainer,
                                const LinearMasterSlaveConstraint &rCloneConstraint,
                                unsigned int ConstraintId,
                                NodeType &rMasterNode,
                                TVariableType &rMasterVariable,
                                NodeType &rSlaveNode,
                                TVariableType &rSlaveVariable,
                                const double Weight,
                                const double Constant = 0.0)
    {
        rSlaveNode.Set(SLAVE);
        ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = rCloneConstraint.Create(ConstraintId, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        p_new_constraint->Set(TO_ERASE);
        if(rSlaveNode.SolutionStepsDataHas(PARTITION_INDEX))
            p_new_constraint->SetValue(PARTITION_INDEX, rSlaveNode.GetSolutionStepValue(PARTITION_INDEX) );
        mNodeIdToConstraintIdsMap[rSlaveNode.Id()].push_back(ConstraintId);
        //TODO: Check if an insert(does a sort) is required here or just use a push_back and then
        //      sort once at the end.
        // rMasterSlaveContainer.insert(rMasterSlaveContainer.begin(), p_new_constraint);
        rMasterSlaveContainer.push_back(p_new_constraint);
    }

    /**
     * @brief Applies the master-slave constraint to enforce the continuity between a given geometry/element and a boundary node
     * @param rGeometry The geometry of the element
     * @param rBoundaryNode The boundary node for which the connections are to be made.
     * @param rShapeFuncWeights The shape function weights for the node in the rGeometry.
     * @param StartIndex The start Index of the constraints which are to be added.
     * @param rConstraintIdVector The vector of the constraints Ids which is accessed with StartIndex.
     * @param rMsContainer The Constraint container to which the contraints are added.
     */
    template <typename TVariableType>
    void ApplyContinuityWithElement(Geometry<NodeType> &rGeometry,
                                    NodeType &rBoundaryNode,
                                    Vector &rShapeFuncWeights,
                                    TVariableType &rVariable,
                                    unsigned int StartIndex,
                                    std::vector<int> &rConstraintIdVector,
                                    MasterSlaveConstraintContainerType &rMsContainer)
    {
        // const auto &r_clone_constraint = (LinearMasterSlaveConstraint)KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
        const auto &r_clone_constraint = new LinearMasterSlaveConstraint();
        // Initialise the boundary nodes dofs to 0 at ever time steps
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) = 0.0;

        for (std::size_t i = 0; i < rGeometry.size(); i++)
        {
            //Interpolation of rVariable
            rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) += rGeometry[i].GetDof(rVariable).GetSolutionStepValue(0) * rShapeFuncWeights[i];
            AddMasterSlaveRelation(rMsContainer, *r_clone_constraint, rConstraintIdVector[StartIndex++], rGeometry[i], rVariable, rBoundaryNode, rVariable, rShapeFuncWeights[i]);
        } // end of loop over host element nodes

        // Setting the buffer 1 same buffer 0
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 1) = rBoundaryNode.FastGetSolutionStepValue(rVariable, 0);
    }

    /**
     * @brief This function reserves the necessary memory for the contraints in each thread
     * @param rModelPart The modelpart to which the constraints are to be added.
     * @param rContainerVector The container vector which has constraints to transfer
     */
    void ReserveMemoryForConstraintContainers(ModelPart &rModelPart, MasterSlaveContainerVectorType &rContainerVector)
    {
#pragma omp parallel
        {
            const IndexType num_constraints_per_thread = (rModelPart.NumberOfNodes() * 4) / omp_get_num_threads();
#pragma omp single
            {
                rContainerVector.resize(omp_get_num_threads());
                for (auto &container : rContainerVector)
                    container.reserve(num_constraints_per_thread);
            }
        }
    }

    /**
     * @brief Given a node, this funcition finds and deletes all the existing constraints for that node
     * @param rBoundaryNode The boundary node for which the connections are to be made.
     */
    int RemoveExistingConstraintsForNode(ModelPart::NodeType &rBoundaryNode)
    {
        ConstraintIdsVectorType constrainIds_for_the_node;
        int removed_counter = 0;
        constrainIds_for_the_node = mNodeIdToConstraintIdsMap[rBoundaryNode.Id()];
        for (auto const &constraint_id : constrainIds_for_the_node)
        {
#pragma omp critical
            {
                mrMainModelPart.RemoveMasterSlaveConstraintFromAllLevels(constraint_id);
                removed_counter++;
            }
        }
        constrainIds_for_the_node.clear();
        return removed_counter;
    }

    /**
     * @brief The function transfers all the constraints in the container vector to the modelpart.
     *          IMPORTANT: The constraints are directly added to the constraints container of the
     *                     modelpart. So the parent and child modelparts have no info about these
     *                     constraints.
     * @param rModelPart The modelpart to which the constraints are to be added.
     * @param rContainerVector The container vector which has constraints to transfer
     */
    void AddConstraintsToModelpart(ModelPart &rModelPart, MasterSlaveContainerVectorType &rContainerVector)
    {
        IndexType n_total_constraints = 0;
        for (auto &container : rContainerVector)
        {
            const int n_constraints = static_cast<int>(container.size());
            n_total_constraints += n_constraints;
        }

        auto &constraints = rModelPart.MasterSlaveConstraints();
        constraints.reserve(n_total_constraints);
        auto &constraints_data = constraints.GetContainer();
        for (auto &container : rContainerVector)
        {
            constraints_data.insert(constraints_data.end(), container.ptr_begin(), container.ptr_end());
        }
        constraints.Sort();
    }

    /**
     * @brief Loops over the nodes of the given modelpart and uses the binlocater to locate them on a
     *          element and formulates respective constraints
     * @param rBoundaryModelPart The modelpart whose nodes are to be found.
     * @param rBinLocator The bin based locator formulated on the background. This is used to locate nodes on rBoundaryModelPart.
     * @param rVelocityMasterSlaveContainerVector the vector of velocity constraints vectors (one for each thread)
     * @param rPressureMasterSlaveContainerVector the vector of pressure constraints vectors (one for each thread)
     */
    void FormulateConstraints(ModelPart &rBoundaryModelPart, PointLocatorType &rBinLocator, MasterSlaveContainerVectorType &rVelocityMasterSlaveContainerVector,
                              MasterSlaveContainerVectorType &rPressureMasterSlaveContainerVector)
    {
        const DataCommunicator &r_comm = mrMainModelPart.GetCommunicator().GetDataCommunicator();
        int mpi_rank = r_comm.Rank();
        int mpi_size = r_comm.Size();
        const bool is_comm_distributed = r_comm.IsDistributed();
        std::vector<NodesContainerType> SendNodes(mpi_size);
        Model &current_model = mrMainModelPart.GetModel();
        auto& gathered_modelpart = is_comm_distributed ? current_model.CreateModelPart("GatheredBoundary") : rBoundaryModelPart;
        if(is_comm_distributed)
            GatherModelPartOnAllRanksUtility::GatherModelPartOnAllRanks(rBoundaryModelPart, gathered_modelpart);

        // WriteModelPart(gathered_modelpart);
        std::vector<int> vector_of_non_found_nodes;
        const int n_boundary_nodes = static_cast<int>(gathered_modelpart.Nodes().size());
        std::vector<int> constraints_id_vector;
        int num_constraints_required = (TDim + 1) * (gathered_modelpart.Nodes().size());
        CreateConstraintIds(constraints_id_vector, num_constraints_required);

        IndexType found_counter = 0;
        std::vector<IndexType> nodes_to_add;

        BuiltinTimer loop_over_b_nodes;
#pragma omp parallel for shared(constraints_id_vector, rVelocityMasterSlaveContainerVector, rPressureMasterSlaveContainerVector, rBinLocator)
        for (int i_bn = 0; i_bn < n_boundary_nodes; ++i_bn)
        {
            auto &ms_velocity_container = rVelocityMasterSlaveContainerVector[omp_get_thread_num()];
            auto &ms_pressure_container = rPressureMasterSlaveContainerVector[omp_get_thread_num()];

            ModelPart::NodesContainerType::iterator i_boundary_node = gathered_modelpart.NodesBegin() + i_bn;
            NodeType& r_boundary_node = *( *(i_boundary_node.base()) );
            const int node_p_index = is_comm_distributed ? r_boundary_node.GetSolutionStepValue(PARTITION_INDEX) : 0;
            unsigned int start_constraint_id = i_bn * (TDim + 1) * (TDim + 1);
            Element::Pointer r_host_element;
            Vector weights;
            bool is_found = SearchNode(rBinLocator, r_boundary_node, r_host_element, weights);
            if (is_found){
                MakeConstraints(r_boundary_node, r_host_element, weights, ms_velocity_container, ms_pressure_container,
                                                            constraints_id_vector, start_constraint_id);
                found_counter += 1;
                if(node_p_index != mpi_rank){
                    auto p_geom = r_host_element->GetGeometry(); 
                    for(auto p_node=p_geom.ptr_begin(); p_node != p_geom.ptr_end(); ++p_node)
                        SendNodes[node_p_index].push_back(*p_node);
                }
            }
        }

        if (is_comm_distributed)
        {
            SynchronizeNodes(mrMainModelPart, SendNodes);
            SynchronizeConstraints(rVelocityMasterSlaveContainerVector);
            SynchronizeConstraints(rPressureMasterSlaveContainerVector);
            current_model.DeleteModelPart("GatheredBoundary");
        }

        double loop_time = loop_over_b_nodes.ElapsedSeconds();
        KRATOS_INFO_IF("ApplyChimera : Loop over boundary nodes took             : ", mEchoLevel > 1) << r_comm.Max(loop_time, 0) << " seconds" << std::endl;
        int global_num_found = r_comm.Sum(found_counter, 0);
        KRATOS_INFO_IF("ApplyChimera : Number of Boundary nodes                  : ", mEchoLevel > 1) << n_boundary_nodes << std::endl;
        KRATOS_INFO_IF("ApplyChimera : Number of Boundary nodes found            : ", mEchoLevel > 1) << global_num_found << std::endl;
        KRATOS_INFO_IF("ApplyChimera : Number of Boundary nodes not found        : ", mEchoLevel > 1) << n_boundary_nodes - global_num_found << std::endl;
        KRATOS_INFO_IF("ApplyChimera : Number of constraints made                : ", mEchoLevel > 1) << global_num_found*9 << std::endl;
    }

    void SynchronizeConstraints(MasterSlaveContainerVectorType& rConstraintsContainterVector)
    {
        ProcessInfo &r_current_process_info = mrMainModelPart.GetProcessInfo();
        const DataCommunicator &r_comm = mrMainModelPart.GetCommunicator().GetDataCommunicator();
        auto& r_nodes = mrMainModelPart.Nodes();
        int mpi_size = r_comm.Size();
        int mpi_rank = r_comm.Rank();
        MasterSlaveContainerVectorType send_constraints(mpi_size);
        MasterSlaveContainerVectorType recv_constraints(mpi_size);
        for(auto& constraints_container : rConstraintsContainterVector){
            for(auto p_constraint=constraints_container.ptr_begin(); p_constraint!=constraints_container.ptr_end(); ++p_constraint){
                int const_rank = (*p_constraint)->GetValue(PARTITION_INDEX);
                if(const_rank != mpi_rank){
                    send_constraints[const_rank].push_back(*p_constraint);
                }
            }
        }
        mrMainModelPart.GetCommunicator().TransferObjects(send_constraints, recv_constraints);

        for(auto& constraints_container : send_constraints){
            for(auto p_constraint=constraints_container.ptr_begin(); p_constraint!=constraints_container.ptr_end(); ++p_constraint){
                rConstraintsContainterVector[0].erase((*p_constraint)->Id()); // IMPORTANT: This means no hybrid parallism
            }
        }
        typedef Element::DofsVectorType DofsVectorType;

        auto& constraints_container = rConstraintsContainterVector[0]; // This does not matter for OMP and MPI hybrid parallism
        for(auto& recv_constraints_container : recv_constraints)
        {
            for(auto p_constraint=recv_constraints_container.ptr_begin(); p_constraint!=recv_constraints_container.ptr_end(); ++p_constraint){
                DofsVectorType slave_dof_list, master_dof_list, new_master_dof_list, new_slave_dof_list;
                (*p_constraint)->GetDofList(slave_dof_list, master_dof_list, r_current_process_info);
                for(auto master_dof : master_dof_list){
                    const int master_node_id = master_dof->Id();
                    auto it_master_node = r_nodes.find(master_node_id);
                    if( it_master_node != r_nodes.end()){
                        if((*it_master_node).GetSolutionStepValue(PARTITION_INDEX) != mpi_rank){
                            auto& r_master_variable = master_dof->GetVariable();
                            auto p_new_master = (*it_master_node).pGetDof(r_master_variable);
                            new_master_dof_list.push_back(p_new_master);
                        }else{
                            new_master_dof_list.push_back(master_dof);
                        }
                    }
                }
                (*p_constraint)->SetMasterDofsVector(new_master_dof_list);

                for(auto slave_dof : slave_dof_list){
                    const int slave_node_id = slave_dof->Id();
                    auto it_slave_node = r_nodes.find(slave_node_id);
                    if( it_slave_node != r_nodes.end()){
                        if((*it_slave_node).GetSolutionStepValue(PARTITION_INDEX) != mpi_rank){
                            KRATOS_ERROR<<"A slave which do not belong to this partition "<<mpi_rank<<" is received !"<<std::endl;
                        }else{
                            auto& r_slave_variable = slave_dof->GetVariable();
                            auto p_new_slave = (*it_slave_node).pGetDof(r_slave_variable);
                            (*it_slave_node).Set(SLAVE);
                            new_slave_dof_list.push_back(p_new_slave);
                        }
                    } else {
                        KRATOS_ERROR<<"A slave which do not belong to this partition "<<mpi_rank<<" is received !"<<std::endl;
                    }
                }
                (*p_constraint)->SetSlaveDofsVector(new_slave_dof_list);


                constraints_container.push_back(*p_constraint);
            }
        }
    }


    void SynchronizeNodes(ModelPart& rModelpart, std::vector<NodesContainerType>& rSendNodes)
    {
        const DataCommunicator &r_comm = rModelpart.GetCommunicator().GetDataCommunicator();
        const int mpi_size = r_comm.Size();
        // const int mpi_rank = r_comm.Rank();
        std::vector<NodesContainerType> RecvNodes(mpi_size);
        BuiltinTimer send_recv_nodes;
        rModelpart.GetCommunicator().TransferObjects(rSendNodes, RecvNodes);
        double time_send_recv = send_recv_nodes.ElapsedSeconds();
        KRATOS_INFO_IF("SynchronizeNodes : Time taken to send recv nodes         : ", mEchoLevel > 1) << r_comm.Max(time_send_recv, 0) << std::endl;

        BuiltinTimer add_nodes;
        for (unsigned int i = 0; i < RecvNodes.size(); i++)
        {
            for (NodesContainerType::iterator it = RecvNodes[i].begin();
                 it != RecvNodes[i].end(); ++it){
                        auto p_node = *it.base();
                        mRemoteNodes.push_back(p_node->Id());
                        rModelpart.Nodes().push_back(p_node);
                 }
        }
        double time_add_nodes = add_nodes.ElapsedSeconds();
        KRATOS_INFO_IF("SynchronizeNodes : Time taken to add nodes               : ", mEchoLevel > 1) << r_comm.Max(time_add_nodes, 0)  << std::endl;

        BuiltinTimer unique_nodes;
        rModelpart.Nodes().Unique();
        double time_unique_nodes = add_nodes.ElapsedSeconds();
        KRATOS_INFO_IF("SynchronizeNodes : Time taken to unique nodes            : ", mEchoLevel > 1) << r_comm.Max(time_unique_nodes, 0) << std::endl;

        int temp = rModelpart.Nodes().size();
        KRATOS_ERROR_IF(temp != int(rModelpart.Nodes().size()))
            << "the rModelpart has repeated nodes";
        rSendNodes.clear();
        RecvNodes.clear();
    }
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Extracts the patch boundary modelpart
     * @param rBackgroundBoundaryModelpart background boundary to remove out of domain patch
     * @param PatchParameters Parameters/Settings for the Patch
     * @param DomainType Flag specifying if the background is the main bg or not
     */
    ModelPart &ExtractPatchBoundary(const Parameters PatchParameters, ModelPart &rBackgroundBoundaryModelpart, const ChimeraHoleCuttingUtility::Domain DomainType)
    {
        Model &current_model = mrMainModelPart.GetModel();
        const auto &r_comm = mrMainModelPart.GetCommunicator().GetDataCommunicator();
        const std::string patch_boundary_mp_name = PatchParameters["boundary_model_part_name"].GetString();

        if (!current_model.HasModelPart(patch_boundary_mp_name))
        {
            ModelPart &r_patch_model_part = current_model.GetModelPart(PatchParameters["model_part_name"].GetString());
            ModelPart &r_modified_patch_model_part = r_patch_model_part.CreateSubModelPart(mModifiedName);
            ModelPart &r_modified_patch_boundary_model_part = r_modified_patch_model_part.CreateSubModelPart(mBoundaryName+r_modified_patch_model_part.Name());
            BuiltinTimer distance_calc_time_patch;
            DistanceCalculationUtility <TDim>::CalculateDistance(r_patch_model_part,
                                                                                                            rBackgroundBoundaryModelpart);
            KRATOS_INFO_IF("Distance calculation on patch took                       : ", mEchoLevel > 0)<< distance_calc_time_patch.ElapsedSeconds()<< " seconds"<< std::endl;
            //TODO: Below is brutforce. Check if the boundary of bg is actually cutting the patch.
            BuiltinTimer rem_out_domain_time;
            ChimeraHoleCuttingUtility().RemoveOutOfDomainElements<TDim>(r_patch_model_part, r_modified_patch_model_part, DomainType, ChimeraHoleCuttingUtility::SideToExtract::INSIDE);
            auto rem_out_domain_time_elapsed = rem_out_domain_time.ElapsedSeconds();
            r_comm.Max(rem_out_domain_time_elapsed, 0);
            KRATOS_INFO_IF("ApplyChimera : Removing out of domain patch took         : ", mEchoLevel > 0) << rem_out_domain_time_elapsed << " seconds" << std::endl;

            BuiltinTimer patch_boundary_extraction_time;
            ChimeraHoleCuttingUtility().ExtractBoundaryMesh<TDim>(r_modified_patch_model_part, r_modified_patch_boundary_model_part);
            auto patch_boundary_extraction_time_elapsed = patch_boundary_extraction_time.ElapsedSeconds();
            r_comm.Max(patch_boundary_extraction_time_elapsed, 0);
            KRATOS_INFO_IF("ApplyChimera : Extraction of patch boundary took         : ", mEchoLevel > 0) << patch_boundary_extraction_time_elapsed << " seconds" << std::endl;

            return r_modified_patch_boundary_model_part;
        }
        else
        {
            return current_model.GetModelPart(patch_boundary_mp_name);
        }
    }

    /**
     * @brief Creates or returns an existing point locator on a given ModelPart
     * @param rModelPart The modelpart on which a point locator is to be obtained.
     */
    PointLocatorPointerType GetPointLocator(ModelPart &rModelPart)
    {
        if (mPointLocatorsMap.count(rModelPart.Name()) == 0)
        {
            PointLocatorPointerType p_point_locator = Kratos::make_shared<PointLocatorType>(rModelPart);
            p_point_locator->UpdateSearchDatabase();
            mPointLocatorsMap[rModelPart.Name()] = p_point_locator;
            return p_point_locator;
        }
        else
        {
            auto p_point_locator = mPointLocatorsMap.at(rModelPart.Name());
            if (mReformulateEveryStep)
                p_point_locator->UpdateSearchDatabase();
            return p_point_locator;
        }
    }


    /**
     * @brief Searches for a given node using given locator and adds the velocity
     *        and pressureconstraints to the respecitive containers.
     * @param rBinLocator The bin based locator formulated on the background. This is used to locate nodes on rBoundaryModelPart.
     * @param pNodeToFind The node which is to be found
     * @param[out] prHostElement The element where the node is found.
     * @param[out] rWeights the values of the shape functions at the node inside the elements.
     */
    bool SearchNode(PointLocatorType &rBinLocator, NodeType& rNodeToFind, Element::Pointer& prHostElement, Vector& rWeights)
    {
        const int max_results = 10000;
        typename PointLocatorType::ResultContainerType results(max_results);
        typename PointLocatorType::ResultIteratorType result_begin = results.begin();

        bool is_found = false;
        is_found = rBinLocator.FindPointOnMesh(rNodeToFind.Coordinates(), rWeights, prHostElement, result_begin, max_results);

        bool node_coupled = false;
        if (rNodeToFind.IsDefined(VISITED))
            node_coupled = rNodeToFind.Is(VISITED);

        if (node_coupled && is_found)
        {
            auto &constrainIds_for_the_node = mNodeIdToConstraintIdsMap[rNodeToFind.Id()];
            for (auto const &constraint_id : constrainIds_for_the_node)
            {
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

    /**
     * @brief Searches for a given node using given locator and adds the velocity
     *        and pressureconstraints to the respecitive containers.
     * @param pNodeToFind The node which is to be found
     * @param pHostElement The element where the node is found.
     * @param rWeights The weights (#Nodes on host elem) for constraint relations
     * @param rVelocityMsConstraintsVector the velocity constraints vector
     * @param rPressureMsConstraintsVector the pressure constraints vector
     * @param rConstraintIdVector the vector of constraint ids which are to be used.
     * @param StartConstraintId the start index of the constraints
     */
    void MakeConstraints(NodeType& rNodeToFind, Element::Pointer& rHostElement, Vector& rWeights,
                                      MasterSlaveConstraintContainerType &rVelocityMsConstraintsVector,
                                      MasterSlaveConstraintContainerType &rPressureMsConstraintsVector,
                                      std::vector<int> &rConstraintIdVector,
                                      const IndexType StartConstraintId)
    {
        Geometry<NodeType> &r_geom = rHostElement->GetGeometry();
        int init_index = 0;
        ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, VELOCITY_X, StartConstraintId + init_index, rConstraintIdVector, rVelocityMsConstraintsVector);
        init_index += (TDim + 1);
        ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, VELOCITY_Y, StartConstraintId + init_index, rConstraintIdVector, rVelocityMsConstraintsVector);
        init_index += (TDim + 1);
        if (TDim == 3)
        {
            ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, VELOCITY_Z, StartConstraintId + init_index, rConstraintIdVector, rVelocityMsConstraintsVector);
            init_index += (TDim + 1);
        }
        ApplyContinuityWithElement(r_geom, rNodeToFind, rWeights, PRESSURE, StartConstraintId + init_index, rConstraintIdVector, rPressureMsConstraintsVector);
        init_index += (TDim + 1);
    }


    /**
     * @brief Utility function to print out various intermediate modelparts
     * @param rModelPart ModelPart to be printed.
     */
    void WriteModelPart(ModelPart &rModelPart)
    {
        Parameters vtk_parameters(R"(
                {
                    "output_control_type"                : "step",
                    "output_frequency"                   : 1,
                    "file_format"                        : "binary",
                    "output_precision"                   : 3,
                    "output_sub_model_parts"             : false,
                    "folder_name"                        : "test_vtk_output",
                    "save_output_files_in_folder"        : false,
                    "nodal_solution_step_data_variables" : [],
                    "nodal_data_value_variables"         : [],
                    "element_flags"                      : [],
                    "nodal_flags"                        : [],
                    "element_data_value_variables"       : [],
                    "condition_data_value_variables"     : [],
                    "write_ids"                          :true
                }
                )");

        VtkOutput vtk_output(rModelPart, vtk_parameters);
        vtk_output.PrintOutput();
    }

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ApplyChimera &operator=(ApplyChimera const &rOther);

    ///@}
}; // Class ApplyChimera

} // namespace Kratos.

#endif //  KRATOS_APPLY_CHIMERA_H_INCLUDED defined
