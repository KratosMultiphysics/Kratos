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

#if !defined(KRATOS_CUSTOM_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED)
#define KRATOS_CUSTOM_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED

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
#include "processes/calculate_signed_distance_to_3d_condition_skin_process.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "input_output/vtk_output.h"

// Application includes
#include "chimera_application_variables.h"
#include "custom_utilities/hole_cutting_utility.h"
#include "utilities/binbased_fast_point_locator.h"
#include "custom_processes/calculate_signed_distance_to_2d_condition_skin_process.h"

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

template <int TDim, class TDistanceCalculatorType>
class ApplyChimeraProcessMonolithic : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    typedef ProcessInfo::Pointer                                                                    ProcessInfoPointerType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>>  VariableComponentType;
    typedef std::size_t                                                                             IndexType;
    typedef Kratos::Variable<double>                                                                VariableType;
    typedef std::vector<IndexType>                                                                  ConstraintIdsVectorType;
    typedef typename ModelPart::MasterSlaveConstraintType                                           MasterSlaveConstraintType;
    typedef typename ModelPart::MasterSlaveConstraintContainerType                                  MasterSlaveConstraintContainerType;
    typedef std::vector<MasterSlaveConstraintContainerType>                                         MasterSlaveContainerVectorType;

    typedef BinBasedFastPointLocator<TDim>                                                          PointLocatorType;
    typedef TDistanceCalculatorType                                                                 DistanceCalculatorType;
    typedef ChimeraHoleCuttingUtility                                                               HoleCuttingUtilityType;

    typedef typename PointLocatorType::Pointer PointLocatorPointerType;

    KRATOS_CLASS_POINTER_DEFINITION(ApplyChimeraProcessMonolithic);

    ///@}
    ///@name Life Cycle
    ///@{

    ApplyChimeraProcessMonolithic(ModelPart &rMainModelPart, Parameters iParameters);

    /// Destructor.
    virtual ~ApplyChimeraProcessMonolithic() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void ExecuteInitializeSolutionStep() override;

    virtual void ExecuteFinalizeSolutionStep() override;

    virtual std::string Info() const override
    {
        return "ApplyChimeraProcessMonolithic";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "ApplyChimeraProcessMonolithic" << std::endl;
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const override
    {
        KRATOS_INFO("\nNumber of slave nodes :: ") << std::endl;
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
    virtual void DoChimeraLoop();

    /**
     * @brief Formulates the Chimera conditions with a given set of background and patch combination.
     * @param BackgroundParam Parameters/Settings for the background
     * @param PatchParameters Parameters/Settings for the Patch
     * @param MainDomainOrNot Flag specifying if the background is the main bg or not
     */
    virtual void FormulateChimera(const Parameters BackgroundParam, const Parameters PatchParameters, int MainDomainOrNot);

    /**
     * @brief Creates a vector of unique constraint ids based on how many required and how many are already present in the mrModelPart.
     * @param rIdVector The vector which is populated with unique constraint ids.
     * @param NumberOfConstraintsRequired The number of further constraints required. used for calculation of unique ids.
     */
    virtual void CreateConstraintIds(std::vector<int> &rIdVector, const IndexType NumberOfConstraintsRequired);

    /**
     * @brief Applies the continuity between the boundary modelpart and the background.
     * @param rBoundaryModelPart The boundary modelpart for which the continuity is to be enforced.
     * @param pBinLocator The bin based locator formulated on the background. This is used to locate nodes on rBoundaryModelPart.
     */
    virtual void ApplyContinuityWithMpcs(ModelPart &rBoundaryModelPart, PointLocatorPointerType &pBinLocator);

    /**
     * @brief Computes the bounding box of the modelpart given. The low and high points (brute force way)
     * @param rModelPart Modelpart for which the bounding box is to be computed.
     * @param rLowPoint The lowest point in the modelpart (returned)
     * @param rHighPoint The highest point in the modelpart (returned)
     */
    virtual void GetBoundingBox(ModelPart &rModelPart, std::vector<double> &rLowPoint, std::vector<double> &rHighPoint);

    /**
     * @brief Checks if two given modelparts (A and B) have bounding box overlaps
     *                  Here in Chimera A is usually the background and B is the patch
     * @param rModelPartA ModelPartA
     * @param rModelPartB ModelPartB
     * @return bool if the bounding boxes intersect or not.
     */
    virtual bool BoundingBoxTest(ModelPart &rModelPartA, ModelPart &rModelPartB);

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
    inline void AddMasterSlaveRelation(MasterSlaveConstraintContainerType &rMasterSlaveContainer,
                                       const LinearMasterSlaveConstraint &rCloneConstraint,
                                       unsigned int ConstraintId,
                                       Node<3> &rMasterNode,
                                       TVariableType &rMasterVariable,
                                       Node<3> &rSlaveNode,
                                       TVariableType &rSlaveVariable,
                                       const double Weight,
                                       const double Constant = 0.0);


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
                                    TVariableType& rVariable,
                                    unsigned int StartIndex,
                                    std::vector<int> &rConstraintIdVector,
                                    MasterSlaveConstraintContainerType &rMsContainer);


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

    void WriteModelPart(ModelPart& rModelPart)
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
                    "nodal_solution_step_data_variables" : ["VELOCITY","PRESSURE"],
                    "nodal_data_value_variables"         : [],
                    "element_flags"                      : ["ACTIVE"],
                    "element_data_value_variables"       : [],
                    "condition_data_value_variables"     : []
                }
                )");

        VtkOutput vtk_output (rModelPart, vtk_parameters);
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
    ApplyChimeraProcessMonolithic &operator=(ApplyChimeraProcessMonolithic const &rOther);

    ///@}
}; // Class ApplyChimeraProcessMonolithic

} // namespace Kratos.

#endif //  KRATOS_CUSTOM_APPLY_CHIMERA_MONOLITHIC_H_INCLUDED defined
