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


#if !defined(KRATOS_APPLY_CHIMERA_H_INCLUDED)
#define KRATOS_APPLY_CHIMERA_H_INCLUDED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/linear_master_slave_constraint.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "utilities/binbased_fast_point_locator.h"

// Application includes
#include "chimera_application_variables.h"
#include "custom_utilities/distance_calcuation_utility.h"
#include "custom_utilities/hole_cutting_utility.h"

namespace Kratos {

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

/**
 * @class ApplyChimera
 *
 * @ingroup ChimeraApplication
 *
 * @brief This class contains methods applies the continuity between the patch and background using linear master-slave constraints.
 * @details This serves as a base class for monolithic and fractional step processes which choose who and where the constraints created are stored. for example velocity and pressure constraints are to be saved seperately for fractional step. 
 *
*/

template <int TDim>
class KRATOS_API(CHIMERA_APPLICATION) ApplyChimera : public Process {
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    typedef ModelPart::IndexType IndexType;
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
     * @param rMainModelPart The reference to the modelpart which will be used
     * for computations later on.
     * @param iParameters The settings parameters.
     */
    explicit ApplyChimera(ModelPart& rMainModelPart, Parameters iParameters);

    /// Destructor.
    virtual ~ApplyChimera()=default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void SetEchoLevel(int EchoLevel);

    void SetReformulateEveryStep(bool Reformulate);

    virtual void ExecuteInitializeSolutionStep() override;

    virtual void ExecuteFinalizeSolutionStep() override;

    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override;

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
    ModelPart& mrMainModelPart;
    int mNumberOfLevels;
    Parameters mParameters;
    std::unordered_map<IndexType, ConstraintIdsVectorType> mNodeIdToConstraintIdsMap;
    int mEchoLevel;
    bool mReformulateEveryStep;
    std::map<std::string, PointLocatorPointerType> mPointLocatorsMap;
    bool mIsFormulated;

    // Modelpart names which are generated here
    const std::string mModifiedName = "ChimeraModified";
    const std::string mBoundaryName = "ChimeraBoundary";
    const std::string mHoleName = "ChimeraHole";
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Does a loop on the background and patch combinations possible and
     * uses FormulateChimera method.
     */
    virtual void DoChimeraLoop();

    /**
     * @brief Formulates the Chimera conditions with a given set of background
     * and patch combination.
     * @param BackgroundParam Parameters/Settings for the background
     * @param PatchParameters Parameters/Settings for the Patch
     * @param DomainType Flag specifying if the background is the main bg or not
     */
    virtual void FormulateChimera(const Parameters BackgroundParam,
                                  const Parameters PatchParameters,
                                  ChimeraHoleCuttingUtility::Domain DomainType);

    /**
     * @brief Creates a vector of unique constraint ids based on how many
     * required and how many are already present in the mrModelPart.
     * @param rIdVector The vector which is populated with unique constraint
     * ids.
     * @param NumberOfConstraintsRequired The number of further constraints
     * required. used for calculation of unique ids.
     */
    void CreateConstraintIds(std::vector<int>& rIdVector,
                                     const IndexType NumberOfConstraintsRequired);
    /**
     * @brief Applies the continuity between the boundary modelpart and the
     * background.
     * @param rBoundaryModelPart The boundary modelpart for which the continuity
     * is to be enforced.
     * @param pBinLocator The bin based locator formulated on the background.
     * This is used to locate nodes of rBoundaryModelPart on background.
     */
    virtual void ApplyContinuityWithMpcs(ModelPart& rBoundaryModelPart,
                                         PointLocatorType& pBinLocator);

    /**
     * @brief Applies the master-slave constraint between the given master and
     * slave nodes with corresponding variable.
     * @param rMasterSlaveContainer The container to which the constraint to be
     * added (useful to so OpenMP loop)
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
    void AddMasterSlaveRelation(MasterSlaveConstraintContainerType& rMasterSlaveContainer,
                                const LinearMasterSlaveConstraint& rCloneConstraint,
                                unsigned int ConstraintId,
                                NodeType& rMasterNode,
                                TVariableType& rMasterVariable,
                                NodeType& rSlaveNode,
                                TVariableType& rSlaveVariable,
                                const double Weight,
                                const double Constant = 0.0);

    /**
     * @brief Applies the master-slave constraint to enforce the continuity
     * between a given geometry/element and a boundary node
     * @param rGeometry The geometry of the element
     * @param rBoundaryNode The boundary node for which the connections are to
     * be made.
     * @param rShapeFuncWeights The shape function weights for the node in the
     * rGeometry.
     * @param StartIndex The start Index of the constraints which are to be
     * added.
     * @param rConstraintIdVector The vector of the constraints Ids which is
     * accessed with StartIndex.
     * @param rMsContainer The Constraint container to which the contraints are
     * added.
     */
    template <typename TVariableType>
    void ApplyContinuityWithElement(Geometry<NodeType>& rGeometry,
                                    NodeType& rBoundaryNode,
                                    Vector& rShapeFuncWeights,
                                    TVariableType& rVariable,
                                    unsigned int StartIndex,
                                    std::vector<int>& rConstraintIdVector,
                                    MasterSlaveConstraintContainerType& rMsContainer);

    /**
     * @brief This function reserves the necessary memory for the contraints in
     * each thread
     * @param rModelPart The modelpart to which the constraints are to be added.
     * @param rContainerVector The container vector which has constraints to
     * transfer
     */
    void ReserveMemoryForConstraintContainers(ModelPart& rModelPart,
                                              MasterSlaveContainerVectorType& rContainerVector);

    /**
     * @brief Given a node, this funcition finds and deletes all the existing
     * constraints for that node
     * @param rBoundaryNode The boundary node for which the connections are to
     * be made.
     */
    int RemoveExistingConstraintsForNode(ModelPart::NodeType& rBoundaryNode);

    /**
     * @brief The function transfers all the constraints in the container vector
     * to the modelpart.
     *          IMPORTANT: The constraints are directly added to the constraints
     * container of the
     *                     modelpart. So the parent and child modelparts have no
     * info about these
     *                     constraints.
     * @param rModelPart The modelpart to which the constraints are to be added.
     * @param rContainerVector The container vector which has constraints to
     * transfer
     */
    void AddConstraintsToModelpart(ModelPart& rModelPart,
                                   MasterSlaveContainerVectorType& rContainerVector);

    /**
     * @brief Loops over the nodes of the given modelpart and uses the
     * binlocater to locate them on a
     *          element and formulates respective constraints
     * @param rBoundaryModelPart The modelpart whose nodes are to be found.
     * @param rBinLocator The bin based locator formulated on the background.
     * This is used to locate nodes on rBoundaryModelPart.
     * @param rVelocityMasterSlaveContainerVector the vector of velocity
     * constraints vectors (one for each thread)
     * @param rPressureMasterSlaveContainerVector the vector of pressure
     * constraints vectors (one for each thread)
     */
    void FormulateConstraints(ModelPart& rBoundaryModelPart,
                              PointLocatorType& rBinLocator,
                              MasterSlaveContainerVectorType& rVelocityMasterSlaveContainerVector,
                              MasterSlaveContainerVectorType& rPressureMasterSlaveContainerVector);
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
     * @param rBackgroundBoundaryModelpart background boundary to remove out of
     * domain patch
     * @param PatchParameters Parameters/Settings for the Patch
     * @param DomainType Flag specifying if the background is the main bg or not
     */
    ModelPart& ExtractPatchBoundary(const Parameters PatchParameters,
                                    ModelPart& rBackgroundBoundaryModelpart,
                                    const ChimeraHoleCuttingUtility::Domain DomainType);

    /**
     * @brief Creates or returns an existing point locator on a given ModelPart
     * @param rModelPart The modelpart on which a point locator is to be
     * obtained.
     */
    PointLocatorPointerType GetPointLocator(ModelPart& rModelPart);

    /**
     * @brief Searches for a given node using given locator and adds the
     * velocity
     *        and pressureconstraints to the respective containers.
     * @param rBinLocator The bin based locator formulated on the background.
     * This is used to locate nodes on rBoundaryModelPart.
     * @param pNodeToFind The node which is to be found
     * @param[out] prHostElement The element where the node is found.
     * @param[out] rWeights the values of the shape functions at the node inside
     * the elements.
     */
    bool SearchNode(PointLocatorType& rBinLocator,
                    NodeType& rNodeToFind,
                    Element::Pointer& prHostElement,
                    Vector& rWeights);

    /**
     * @brief Creates the constraints and adds them respective containers.
     * @param pNodeToFind The node which is to be found
     * @param pHostElement The element where the node is found.
     * @param rWeights The weights (#Nodes on host elem) for constraint
     * relations
     * @param rVelocityMsConstraintsVector the velocity constraints vector
     * @param rPressureMsConstraintsVector the pressure constraints vector
     * @param rConstraintIdVector the vector of constraint ids which are to be
     * used.
     * @param StartConstraintId the start index of the constraints
     */
    void MakeConstraints(NodeType& rNodeToFind,
                         Element::Pointer& rHostElement,
                         Vector& rWeights,
                         MasterSlaveConstraintContainerType& rVelocityMsConstraintsVector,
                         MasterSlaveConstraintContainerType& rPressureMsConstraintsVector,
                         std::vector<int>& rConstraintIdVector,
                         const IndexType StartConstraintId);
    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
}; // Class ApplyChimera

} // namespace Kratos.

#endif //  KRATOS_APPLY_CHIMERA_H_INCLUDED defined
