// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
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

// Application includes
#include "chimera_application_variables.h"
#include "custom_utilities/hole_cutting_utility.h"

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
    typedef ChimeraHoleCuttingUtility HoleCuttingUtilityType;
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
        Parameters default_parameters(R"(
            {
               	"chimera_parts"   :   [
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
        KRATOS_ERROR_IF(mNumberOfLevels < 2) << "Chimera requires atleast two levels !. Put Background in one level and the patch in second one." << std::endl;

        ProcessInfoPointerType info = mrMainModelPart.pGetProcessInfo();
        mpHoleCuttingUtility = ChimeraHoleCuttingUtility::Pointer(new ChimeraHoleCuttingUtility());
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

    virtual void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        // Actual execution of the functionality of this class
        DoChimeraLoop();
        KRATOS_CATCH("");
    }

    virtual void ExecuteFinalizeSolutionStep() override
    {
        //for multipatch
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

    ChimeraHoleCuttingUtility::Pointer mpHoleCuttingUtility;
    ModelPart &mrMainModelPart;
    double mOverlapDistance;
    IndexType mNumberOfLevels;
    Parameters mParameters;
    std::unordered_map<IndexType, ConstraintIdsVectorType> mNodeIdToConstraintIdsMap;
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

        int i_current_level = 0;
        for (const auto &current_level : mParameters)
        {
            int is_main_background = 1;
            for (const auto &background_patch_param : current_level)
            { // Gives the current background patch
                for (IndexType i_slave_level = i_current_level + 1; i_slave_level < mNumberOfLevels; ++i_slave_level)
                {
                    for (const auto &slave_patch_param : mParameters[i_slave_level]) // Loop over all other slave patches
                    {
                        if (i_current_level == 0) // a check to identify computational Domain boundary
                            is_main_background = -1;
                        FormulateChimera(background_patch_param, slave_patch_param, is_main_background);
                        KRATOS_INFO("Formulating Chimera for the combination :: \n") << "Background" << background_patch_param << "\n Patch::" << slave_patch_param << std::endl;
                    }
                }
            }
            ++i_current_level;
        }

        KRATOS_INFO("End of chimera loop") << std::endl;

        KRATOS_INFO("Total number of constraints created so far") << mrMainModelPart.NumberOfMasterSlaveConstraints() << std::endl;
    }

    /**
     * @brief Formulates the Chimera conditions with a given set of background and patch combination.
     * @param BackgroundParam Parameters/Settings for the background
     * @param PatchParameters Parameters/Settings for the Patch
     * @param MainDomainOrNot Flag specifying if the background is the main bg or not
     */
    virtual void FormulateChimera(const Parameters BackgroundParam, const Parameters PatchParameters, int MainDomainOrNot)
    {
        ModelPart &r_background_model_part = mrMainModelPart.GetSubModelPart(BackgroundParam["model_part_name"].GetString());
        ModelPart &r_domain_boundary_model_part = mrMainModelPart.GetSubModelPart(BackgroundParam["model_part_inside_boundary_name"].GetString());

        ModelPart &r_patch_model_part = mrMainModelPart.GetSubModelPart(PatchParameters["model_part_name"].GetString());
        ModelPart &r_patch_inside_boundary_model_part = mrMainModelPart.GetSubModelPart(PatchParameters["model_part_inside_boundary_name"].GetString());

        const double overlap_bg = BackgroundParam["overlap_distance"].GetDouble();
        const double overlap_pt = PatchParameters["overlap_distance"].GetDouble();

        const double over_lap_distance = (overlap_bg > overlap_pt) ? overlap_bg : overlap_pt;

        PointLocatorPointerType p_point_locator_on_background = Kratos::make_shared<PointLocatorType>(r_background_model_part);
        p_point_locator_on_background->UpdateSearchDatabase();
        PointLocatorPointerType p_pointer_locator_on_patch = Kratos::make_shared<PointLocatorType>(r_patch_model_part);
        p_pointer_locator_on_patch->UpdateSearchDatabase();

        const double eps = 1e-12;
        KRATOS_ERROR_IF(over_lap_distance < eps) << "Overlap distance should be a positive and non-zero number." << std::endl;

        if (over_lap_distance > eps)
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
                CalculateDistance(r_patch_model_part, r_domain_boundary_model_part, over_lap_distance);
                //TODO: Below is brutforce. Check if the boundary of bg is actually cutting the patch.
                mpHoleCuttingUtility->RemoveOutOfDomainElements(r_patch_model_part, r_modified_patch_model_part, MainDomainOrNot, false);
            }

            mpHoleCuttingUtility->FindOutsideBoundaryOfModelPartGivenInside(r_modified_patch_model_part, r_patch_inside_boundary_model_part, r_modified_patch_boundary_model_part);
            CalculateDistance(r_background_model_part, r_modified_patch_boundary_model_part, over_lap_distance);
            mpHoleCuttingUtility->CreateHoleAfterDistance(r_background_model_part, r_hole_model_part, r_hole_boundary_model_part, over_lap_distance);

            //WriteModelPart(r_hole_model_part);
            //WriteModelPart(r_modified_patch_boundary_model_part);
            //WriteModelPart(r_modified_patch_model_part);
            //WriteModelPart(r_hole_boundary_model_part);

            //for multipatch
            const unsigned int n_elements = r_hole_model_part.NumberOfElements();
#pragma omp parallel for
            for (IndexType i_elem = 0; i_elem < n_elements; ++i_elem)
            {
                ModelPart::ElementsContainerType::iterator it_elem = r_hole_model_part.ElementsBegin() + i_elem;
                it_elem->Set(VISITED, true);
            }

            ApplyContinuityWithMpcs(r_modified_patch_boundary_model_part, p_point_locator_on_background);
            ApplyContinuityWithMpcs(r_hole_boundary_model_part, p_pointer_locator_on_patch);

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
     * @brief Computes the bounding box of the modelpart given. The low and high points (brute force way)
     * @param rModelPart Modelpart for which the bounding box is to be computed.
     * @param rLowPoint The lowest point in the modelpart (returned)
     * @param rHighPoint The highest point in the modelpart (returned)
     */
    virtual void GetBoundingBox(ModelPart &rModelPart, std::vector<double> &rLowPoint, std::vector<double> &rHighPoint)
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

    /**
     * @brief Checks if two given modelparts (A and B) have bounding box overlaps
     *                  Here in Chimera A is usually the background and B is the patch
     * @param rModelPartA ModelPartA
     * @param rModelPartB ModelPartB
     * @return bool if the bounding boxes intersect or not.
     */
    bool BoundingBoxTest(ModelPart &rModelPartA, ModelPart &rModelPartB)
    {
        std::vector<double> min_cornerA(3), max_cornerA(3), min_cornerB(3), max_cornerB(3);
        GetBoundingBox(rModelPartA, min_cornerA, max_cornerA);
        GetBoundingBox(rModelPartB, min_cornerB, max_cornerB);
        const int dim = rModelPartA.GetProcessInfo().GetValue(DOMAIN_SIZE);

        for (int i = 0; i < dim; i++)
        {
            if (min_cornerA[i] > max_cornerB[i])
                return false;
            if (max_cornerA[i] < min_cornerB[i])
                return false;
        }
        return true;
    }

    /**
     * @brief Calculates distance on the whole of rBackgroundModelPart from rSkinModelPart
     * @param rBackgroundModelPart The background modelpart where distances are calculated.
     * @param rSkinModelPart The skin modelpart from where the distances are calculated
     */
    void CalculateDistance(ModelPart &rBackgroundModelPart, ModelPart &rSkinModelPart, const double OverlapDistance)
    {
        typedef LinearSolverFactory<SparseSpaceType, LocalSparseSpaceType> LinearSolverFactoryType;
        typedef LinearSolver<SparseSpaceType, TLocalSpaceType> LinearSolverType;
        typedef VariationalDistanceCalculationProcess<TDim, TSparseSpaceType, TLocalSpaceType, LinearSolverType> VariationalDistanceCalculationProcessType;
        typedef CalculateDistanceToSkinProcess<TDim> CalculateDistanceToSkinProcessType;
        IndexType nnodes = static_cast<IndexType>(rBackgroundModelPart.NumberOfNodes());

#pragma omp parallel for
        for (IndexType i_node = 0; i_node < nnodes; ++i_node)
        {
            auto it_node = rBackgroundModelPart.NodesBegin() + i_node;
            double &node_distance = it_node->FastGetSolutionStepValue(DISTANCE);
            node_distance = 0;
        }

        CalculateDistanceToSkinProcessType(rBackgroundModelPart, rSkinModelPart).Execute();
        int num_fixed_nodes = 0;
// Setting the boundary conditions for the VariationalDistanceCalculationProcess
#pragma omp parallel for reduction(+ \
                                   : num_fixed_nodes)
        for (IndexType i_node = 0; i_node < nnodes; ++i_node)
        {
            auto it_node = rBackgroundModelPart.NodesBegin() + i_node;
            const double node_distance = it_node->FastGetSolutionStepValue(DISTANCE);
            if (std::abs(node_distance) < OverlapDistance)
            {
                num_fixed_nodes++;
                it_node->Fix(DISTANCE);
            }
        }

        if (num_fixed_nodes > 0)
        {
            Parameters amgcl_settings(R"(
                {
                    "solver_type"                    : "amgcl"
                }
                )");

            LinearSolverFactoryType const &linear_solver_factory = KratosComponents<LinearSolverFactoryType>::Get("amgcl");
            auto amgcl_solver = linear_solver_factory.Create(amgcl_settings);
            VariationalDistanceCalculationProcessType(rBackgroundModelPart, amgcl_solver).Execute();
        }
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
