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
#include "processes/process.h"
#include "includes/model_part.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "includes/linear_master_slave_constraint.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "input_output/vtk_output.h"
#include "utilities/binbased_fast_point_locator.h"
#include "factories/standard_linear_solver_factory.h"
#include "processes/variational_distance_calculation_process.h"
#include "utilities/parallel_levelset_distance_calculator.h"
#include "utilities/builtin_timer.h"

// Application includes
#include "chimera_application_variables.h"
#include "custom_utilities/hole_cutting_utility.h"
#include "calculate_signed_distance_to_2d_condition_skin_process.h"
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"

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

template <int TDim, class TSparseSpaceType, class TLocalSpaceType>
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
    typedef Kratos::Variable<double> VariableType;
    typedef std::vector<IndexType> ConstraintIdsVectorType;
    typedef typename ModelPart::MasterSlaveConstraintType MasterSlaveConstraintType;
    typedef typename ModelPart::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef std::vector<MasterSlaveConstraintContainerType> MasterSlaveContainerVectorType;
    typedef BinBasedFastPointLocator<TDim> PointLocatorType;
    typedef ChimeraHoleCuttingUtility<TDim> HoleCuttingUtilityType;
    typedef typename ChimeraHoleCuttingUtility<TDim>::Pointer HoleCuttingUtilityTypePointerType;
    typedef typename PointLocatorType::Pointer PointLocatorPointerType;

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
										"overlap_distance":0.0
									}],
									[{
										"model_part_name":"PLEASE_SPECIFY",
										"overlap_distance":0.0
									}],
									[{
										"model_part_name":"PLEASE_SPECIFY",
										"overlap_distance":0.0
									}]
								]
            })");

        mNumberOfLevels = mParameters.size();
        KRATOS_ERROR_IF(mNumberOfLevels < 2) << "Chimera requires atleast two levels !. Put Background in one level and the patch in second one." << std::endl;

        ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
        mpHoleCuttingUtility = Kratos::make_shared<HoleCuttingUtilityType> ();

        mEchoLevel = 0;
        mReformulateEveryStep = false;
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
        // Actual execution of the functionality of this class
        if(mReformulateEveryStep)
            DoChimeraLoop();
        KRATOS_CATCH("");
    }

    virtual void ExecuteFinalizeSolutionStep() override
    {

        //for multipatch
        const unsigned int num_elements = mrMainModelPart.NumberOfElements();
        const auto elem_begin = mrMainModelPart.ElementsBegin();

#pragma omp parallel for
        for (unsigned int i_be = 0; i_be < num_elements; ++i_be)
        {
            auto i_elem = elem_begin + i_be;
            i_elem->Set(VISITED, false);
            i_elem->SetValue(SPLIT_ELEMENT, false);
        }

        if(mReformulateEveryStep)
            mrMainModelPart.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
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
        KRATOS_INFO("ApplyChimera") << std::endl;
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

    HoleCuttingUtilityTypePointerType mpHoleCuttingUtility;
    ModelPart &mrMainModelPart;
    double mOverlapDistance;
    IndexType mNumberOfLevels;
    Parameters mParameters;
    std::unordered_map<IndexType, ConstraintIdsVectorType> mNodeIdToConstraintIdsMap;
    int mEchoLevel;
    bool mReformulateEveryStep;
    std::map<std::string, PointLocatorPointerType> mPointLocatorsMap;
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
        const unsigned int num_elements = mrMainModelPart.NumberOfElements();
        const auto elem_begin = mrMainModelPart.ElementsBegin();


#pragma omp parallel for
        for (unsigned int i_be = 0; i_be < num_elements; ++i_be)
        {
            auto i_elem = elem_begin + i_be;
            if (!i_elem->Is(VISITED)) //for multipatch
                i_elem->Set(ACTIVE, true);
        }

        const unsigned int num_nodes = mrMainModelPart.NumberOfNodes();
        const auto nodes_begin = mrMainModelPart.ElementsBegin();

#pragma omp parallel for
        for (unsigned int i_bn = 0; i_bn < num_nodes; ++i_bn)
        {
            auto i_node = nodes_begin + i_bn;
            i_node->Set(VISITED, false);
        }
        BuiltinTimer do_chimera_loop_time;

        int i_current_level = 0;
        for (const auto &current_level : mParameters)
        {
            int is_main_background = 1;
            for (const auto &background_patch_param : current_level)
            { // Gives the current background patch
                // compute the outerboundary of the background to save

                for (IndexType i_slave_level = i_current_level + 1; i_slave_level < mNumberOfLevels; ++i_slave_level)
                {
                    for (const auto &slave_patch_param : mParameters[i_slave_level]) // Loop over all other slave patches
                    {
                        if (i_current_level == 0) // a check to identify computational Domain boundary
                            is_main_background = -1;
                        ModelPart &r_background_model_part = mrMainModelPart.GetSubModelPart(background_patch_param["model_part_name"].GetString());
                        if(!r_background_model_part.HasSubModelPart("chimera_boundary_mp")){
                            auto& r_boundary_model_part = r_background_model_part.CreateSubModelPart("chimera_boundary_mp");
                            BuiltinTimer extraction_time;
                            if(i_current_level==0){
                                mpHoleCuttingUtility->ExtractBoundaryMesh(r_background_model_part, r_boundary_model_part);
                            }else{
                                mpHoleCuttingUtility->ExtractBoundaryMesh(r_background_model_part, r_boundary_model_part, true);
                            }
                            KRATOS_INFO_IF("Extraction of boundary mesh took : ", mEchoLevel > 0)<< extraction_time.ElapsedSeconds()<< " seconds"<< std::endl;
                        }
                        FormulateChimera(background_patch_param, slave_patch_param, is_main_background);
                        //KRATOS_INFO("Formulating Chimera for the combination :: \n") << "Background" << background_patch_param << "\n Patch::" << slave_patch_param << std::endl;
                    }
                }
            }
            ++i_current_level;
        }

        KRATOS_INFO_IF("Chrimera Initialization took : ", mEchoLevel > 0)
                        << do_chimera_loop_time.ElapsedSeconds()<< " seconds"<< std::endl;

        KRATOS_INFO_IF("Total number of constraints created so far", mEchoLevel > 0) << mrMainModelPart.NumberOfMasterSlaveConstraints() << std::endl;
    }

    /**
     * @brief Formulates the Chimera conditions with a given set of background and patch combination.
     * @param BackgroundParam Parameters/Settings for the background
     * @param PatchParameters Parameters/Settings for the Patch
     * @param MainDomainOrNot Flag specifying if the background is the main bg or not
     */
    virtual void FormulateChimera(const Parameters BackgroundParam, const Parameters PatchParameters, int MainDomainOrNot)
    {
        Model &current_model = mrMainModelPart.GetModel();
        ModelPart &r_background_model_part = mrMainModelPart.GetSubModelPart(BackgroundParam["model_part_name"].GetString());
        ModelPart &r_background_boundary_model_part = r_background_model_part.GetSubModelPart("chimera_boundary_mp");


        ModelPart &r_patch_model_part = mrMainModelPart.GetSubModelPart(PatchParameters["model_part_name"].GetString());

        const double overlap_bg = BackgroundParam["overlap_distance"].GetDouble();
        const double overlap_pt = PatchParameters["overlap_distance"].GetDouble();

        const double over_lap_distance = (overlap_bg > overlap_pt) ? overlap_bg : overlap_pt;

        BuiltinTimer search_creation_time;
        PointLocatorPointerType p_point_locator_on_background = GetPointLocator(r_background_model_part);
        PointLocatorPointerType p_pointer_locator_on_patch = GetPointLocator(r_patch_model_part);
        KRATOS_INFO_IF("Creation of search structures took : ", mEchoLevel > 0)<< search_creation_time.ElapsedSeconds()<< " seconds"<< std::endl;

        const double eps = 1e-12;
        KRATOS_ERROR_IF(over_lap_distance < eps) << "Overlap distance should be a positive and non-zero number." << std::endl;

        if (over_lap_distance > eps)
        {
            ModelPart &r_hole_model_part = current_model.CreateModelPart("HoleModelpart");
            ModelPart &r_hole_boundary_model_part = current_model.CreateModelPart("HoleBoundaryModelPart");
            ModelPart &r_modified_patch_boundary_model_part = current_model.CreateModelPart("ModifiedPatchBoundary");
            ModelPart &r_modified_patch_model_part = current_model.CreateModelPart("ModifiedPatch");

            BuiltinTimer distance_calc_time_patch;
            CalculateDistanceChimeraApplication(r_patch_model_part, r_background_boundary_model_part, over_lap_distance);
            KRATOS_INFO_IF("Distance calculation on patch took : ", mEchoLevel > 0)<< distance_calc_time_patch.ElapsedSeconds()<< " seconds"<< std::endl;
            //TODO: Below is brutforce. Check if the boundary of bg is actually cutting the patch.
            BuiltinTimer rem_out_domain_time;
            mpHoleCuttingUtility->RemoveOutOfDomainElements(r_patch_model_part, r_modified_patch_model_part, MainDomainOrNot, false);
            KRATOS_INFO_IF("Removing out of domain patch took : ", mEchoLevel > 0)<< rem_out_domain_time.ElapsedSeconds()<< " seconds"<< std::endl;

            BuiltinTimer patch_boundary_extraction_time;
            mpHoleCuttingUtility->ExtractBoundaryMesh(r_modified_patch_model_part, r_modified_patch_boundary_model_part);
            KRATOS_INFO_IF("Extraction of patch boundary took : ", mEchoLevel > 0)<< patch_boundary_extraction_time.ElapsedSeconds()<< " seconds"<< std::endl;

            BuiltinTimer bg_distance_calc_time;
            CalculateDistanceChimeraApplication(r_background_model_part, r_modified_patch_boundary_model_part, over_lap_distance);
            KRATOS_INFO_IF("Distance calculation on background took : ", mEchoLevel > 0)<< bg_distance_calc_time.ElapsedSeconds()<< " seconds"<< std::endl;

            BuiltinTimer hole_creation_time;
            mpHoleCuttingUtility->CreateHoleAfterDistance(r_background_model_part, r_hole_model_part, r_hole_boundary_model_part, over_lap_distance);
            KRATOS_INFO_IF("Hole creation took : ", mEchoLevel > 0)<< hole_creation_time.ElapsedSeconds()<< " seconds"<< std::endl;

            // WriteModelPart(r_hole_model_part);
            // WriteModelPart(r_modified_patch_boundary_model_part);
            // WriteModelPart(r_modified_patch_model_part);
            // WriteModelPart(r_patch_model_part);
            // WriteModelPart(r_background_boundary_model_part);
            // WriteModelPart(r_background_model_part);

            const unsigned int n_elements = r_hole_model_part.NumberOfElements();
#pragma omp parallel for
            for (IndexType i_elem = 0; i_elem < n_elements; ++i_elem)
            {
                ModelPart::ElementsContainerType::iterator it_elem = r_hole_model_part.ElementsBegin() + i_elem;
                it_elem->Set(VISITED, true); //for multipatch
            }

            BuiltinTimer mpc_time;
            ApplyContinuityWithMpcs(r_modified_patch_boundary_model_part, p_point_locator_on_background);
            ApplyContinuityWithMpcs(r_hole_boundary_model_part, p_pointer_locator_on_patch);
            KRATOS_INFO_IF("Creation of MPC for chimera took : ", mEchoLevel > 0)<< mpc_time.ElapsedSeconds()<< " seconds"<< std::endl;

            current_model.DeleteModelPart("HoleModelpart");
            current_model.DeleteModelPart("HoleBoundaryModelPart");
            current_model.DeleteModelPart("ModifiedPatchBoundary");
            current_model.DeleteModelPart("ModifiedPatch");
        }
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
    virtual void ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorPointerType &pBinLocator)
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
                                Node<3> &rMasterNode,
                                TVariableType &rMasterVariable,
                                Node<3> &rSlaveNode,
                                TVariableType &rSlaveVariable,
                                const double Weight,
                                const double Constant = 0.0)
    {
        rSlaveNode.Set(SLAVE);
        ModelPart::MasterSlaveConstraintType::Pointer p_new_constraint = rCloneConstraint.Create(ConstraintId, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        p_new_constraint->Set(TO_ERASE);
        mNodeIdToConstraintIdsMap[rSlaveNode.Id()].push_back(ConstraintId);
        rMasterSlaveContainer.insert(rMasterSlaveContainer.begin(), p_new_constraint);
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
    void ApplyContinuityWithElement(Geometry<Node<3>> &rGeometry,
                                    Node<3> &rBoundaryNode,
                                    Vector &rShapeFuncWeights,
                                    TVariableType &rVariable,
                                    unsigned int StartIndex,
                                    std::vector<int> &rConstraintIdVector,
                                    MasterSlaveConstraintContainerType &rMsContainer)
    {
        const auto &r_clone_constraint = (LinearMasterSlaveConstraint)KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
        // Initialise the boundary nodes dofs to 0 at ever time steps
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) = 0.0;
        for (std::size_t i = 0; i < rGeometry.size(); i++)
        {
            //Interpolation of rVariable
            rBoundaryNode.FastGetSolutionStepValue(rVariable, 0) += rGeometry[i].GetDof(rVariable).GetSolutionStepValue(0) * rShapeFuncWeights[i];
            AddMasterSlaveRelation(rMsContainer, r_clone_constraint, rConstraintIdVector[StartIndex++], rGeometry[i], rVariable, rBoundaryNode, rVariable, rShapeFuncWeights[i]);
        } // end of loop over host element nodes

        // Setting the buffer 1 same buffer 0
        rBoundaryNode.FastGetSolutionStepValue(rVariable, 1) = rBoundaryNode.FastGetSolutionStepValue(rVariable, 0);
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
     * @brief Creates or returns an existing point locator on a given ModelPart
     * @param rModelPart The modelpart on which a point locator is to be obtained.
     */
    PointLocatorPointerType GetPointLocator(ModelPart& rModelPart)
    {
        if(mPointLocatorsMap.count(rModelPart.Name()) == 0){
            PointLocatorPointerType p_point_locator = Kratos::make_shared<PointLocatorType>(rModelPart);
            p_point_locator->UpdateSearchDatabase();
            mPointLocatorsMap[rModelPart.Name()] = p_point_locator;
            return p_point_locator;
        } else {
            auto p_point_locator = mPointLocatorsMap.at(rModelPart.Name());
            if(mReformulateEveryStep)
                p_point_locator->UpdateSearchDatabase();
            return p_point_locator;
        }
    }

    /**
     * @brief Calculates distance on the whole of rBackgroundModelPart from rSkinModelPart
     * @param rBackgroundModelPart The background modelpart where distances are calculated.
     * @param rSkinModelPart The skin modelpart from where the distances are calculated
     */
    void CalculateDistance(ModelPart &rBackgroundModelPart, ModelPart &rSkinModelPart, const double OverlapDistance)
    {
        typedef CalculateDistanceToSkinProcess<TDim> CalculateDistanceToSkinProcessType;
        IndexType nnodes = static_cast<IndexType>(rBackgroundModelPart.NumberOfNodes());
        auto p_distance_smoother = Kratos::make_shared<ParallelDistanceCalculator<TDim>>();
        unsigned int max_level = 100;
		double max_distance = 200;

#pragma omp parallel for
        for (IndexType i_node = 0; i_node < nnodes; ++i_node)
        {
            auto it_node = rBackgroundModelPart.NodesBegin() + i_node;
            double &node_distance = it_node->FastGetSolutionStepValue(DISTANCE);
            node_distance = 0;
        }

        CalculateDistanceToSkinProcessType(rBackgroundModelPart, rSkinModelPart).Execute();

        p_distance_smoother->CalculateDistances(rBackgroundModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);
    }

    /**
     * @brief Calculates distance on the whole of rBackgroundModelPart from rSkinModelPart
     * @param rBackgroundModelPart The background modelpart where distances are calculated.
     * @param rSkinModelPart The skin modelpart from where the distances are calculated
     */
    void CalculateDistanceChimeraApplication(ModelPart &rBackgroundModelPart, ModelPart &rSkinModelPart, const double OverlapDistance)
    {
        IndexType nnodes = static_cast<IndexType>(rBackgroundModelPart.NumberOfNodes());
        auto p_distance_smoother = Kratos::make_shared<ParallelDistanceCalculator<TDim>>();
        unsigned int max_level = 100;
		double max_distance = 200;

#pragma omp parallel for
        for (IndexType i_node = 0; i_node < nnodes; ++i_node)
        {
            auto it_node = rBackgroundModelPart.NodesBegin() + i_node;
            double &node_distance = it_node->FastGetSolutionStepValue(DISTANCE);
            node_distance = 0;
        }

        if(TDim ==2)
            CalculateSignedDistanceTo2DConditionSkinProcess(rSkinModelPart, rBackgroundModelPart).Execute();
        else if(TDim==3)
            CalculateSignedDistanceTo3DConditionSkinProcess(rSkinModelPart, rBackgroundModelPart).Execute();

        p_distance_smoother->CalculateDistances(rBackgroundModelPart, DISTANCE, NODAL_AREA, max_level, max_distance);
    }

    /**
     * @brief Utility function to print out various intermediate modelparts
     * @param rModelPart ModelPart to be printed.
     */
    void WriteModelPart(ModelPart &rModelPart)
    {
        Parameters vtk_parameters(R"(
                {
                    "model_part_name"                    : "HoleModelpart",
                    "output_control_type"                : "step",
                    "output_frequency"                   : 1,
                    "file_format"                        : "ascii",
                    "output_precision"                   : 3,
                    "output_sub_model_parts"             : false,
                    "folder_name"                        : "test_vtk_output",
                    "save_output_files_in_folder"        : true,
                    "nodal_solution_step_data_variables" : ["VELOCITY","PRESSURE","DISTANCE"],
                    "nodal_data_value_variables"         : [],
                    "element_flags"                      : ["ACTIVE"],
                    "nodal_flags"                        : ["VISITED"],
                    "element_data_value_variables"       : [],
                    "condition_data_value_variables"     : []
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
