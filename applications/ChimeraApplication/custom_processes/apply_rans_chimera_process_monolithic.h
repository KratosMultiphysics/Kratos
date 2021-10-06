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
//                  Rahul Kikkeri Nagaraja
//

#if !defined(KRATOS_APPLY_RANS_CHIMERA_MONOLITHIC_H_INCLUDED)
#define KRATOS_APPLY_RANS_CHIMERA_MONOLITHIC_H_INCLUDED

// System includes

// External includes

// Project includes
#include "apply_chimera_process.h"

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
 * @class ApplyRANSChimeraProcessMonolithic
 *
 * @ingroup ChimeraApplication
 *
 * @brief This class extends ApplyChimera base class and overrides the function ApplyContinuityWithMpcs to use same container for storing pressure and velocity constraints.
 *
*/
template <int TDim>
class KRATOS_API(CHIMERA_APPLICATION) ApplyRANSChimeraProcessMonolithic
    : public ApplyChimera<TDim> {
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    KRATOS_CLASS_POINTER_DEFINITION(ApplyRANSChimeraProcessMonolithic);
    typedef ApplyChimera<TDim> BaseType;
    typedef typename BaseType::MasterSlaveConstraintContainerType MasterSlaveConstraintContainerType;
    typedef typename BaseType::PointLocatorType PointLocatorType;
    typedef typename BaseType::PointLocatorPointerType PointLocatorPointerType;
    typedef typename BaseType::MasterSlaveContainerVectorType MasterSlaveContainerVectorType;
    typedef typename BaseType::ConstraintIdsVectorType ConstraintIdsVectorType;
    typedef typename BaseType::NodeType NodeType;
    typedef std::vector<std::string> SolvingVariablesVectorType;
    typedef std::map<std::string, PointLocatorPointerType> PointLocatorsMapType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor
     * @param rMainModelPart The reference to the modelpart which will be used
     * for computations later on.
     * @param iParameters The settings parameters.
     * @param SolvingVariablesVector The vector of solving variables name used 
     * to identify solving variables which are constrained at chimera boundary.
     */
    explicit ApplyRANSChimeraProcessMonolithic(ModelPart& rMainModelPart, Parameters iParameters, 
                                               SolvingVariablesVectorType rSolvingVariablesVector);

    /// Destructor.
    ~ApplyRANSChimeraProcessMonolithic() = default;
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    /**
     * @brief Returns the map of submodelpart name and its corresponding PointLocatorPointer.
     */
    PointLocatorsMapType const& GetPointLocatorsMap() const; // removed reference to const return type to reference type

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
    SolvingVariablesVectorType mSolvingVariableNamesVector;
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief Applies the continuity between the boundary modelpart and the
     * background.
     * @param rBoundaryModelPart The boundary modelpart for which the continuity
     * is to be enforced.
     * @param rBinLocator The bin based locator formulated on the background.
     * This is used to locate nodes of rBoundaryModelPart on background.
     */
    void ApplyContinuityWithMpcs(ModelPart& rBoundaryModelPart, PointLocatorType& rBinLocator) override;

    /**
     * @brief This function reserves the necessary memory for the contraints in
     * each thread.
     * @param rModelPart The modelpart to which the constraints are to be added.
     * @param rMsConstraintContainer The vector of master-slave constraints.
     */
    void ReserveMemoryForConstraintContainers(ModelPart& rModelPart,
                                              MasterSlaveConstraintContainerType& rMsConstraintContainer);

    /**
     * @brief Loops over the nodes of the given modelpart and uses the
     * binlocater to locate them on a
     *          element and formulates respective constraints
     * @param rBoundaryModelPart The modelpart whose nodes are to be found.
     * @param rBinLocator The bin based locator formulated on the background.
     * This is used to locate nodes on rBoundaryModelPart.
     * @param rMsConstraintContainer The vector of master-slave constraints.
     */
    void FormulateConstraints(ModelPart& rBoundaryModelPart,
                              PointLocatorType& rBinLocator,
                              MasterSlaveConstraintContainerType& rMsConstraintContainer);

    /**
     * @brief Creates the constraints and adds them respective containers.
     * @param pNodeToFind The node which is to be found
     * @param pHostElement The element where the node is found.
     * @param rWeights The weights (#Nodes on host elem) for constraint
     * relations
     * @param rMsConstraintContainer The vector of master-slave constraints.
     * @param rConstraintIdVector the vector of constraint ids which are to be
     * used.
     * @param StartConstraintId the start index of the constraints
     */
    void MakeConstraints(NodeType& rNodeToFind,
                         Element::Pointer& rHostElement,
                         Vector& rWeights,
                         MasterSlaveConstraintContainerType& rMsConstraintContainer,
                         std::vector<int>& rConstraintIdVector,
                         const IndexType StartConstraintId);

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
     * @brief The function transfers all the constraints in the container to 
     * the modelpart.
     * IMPORTANT: The constraints are directly added to the constraints
     * container of the modelpart. So the parent and child modelparts have no
     * info about these constraints.
     * @param rModelPart The modelpart to which the constraints are to be added.
     * @param rContainerVector The vector which has constraints to transfer.
     */
    void AddConstraintsToModelpart(ModelPart& rModelPart,
                                   MasterSlaveConstraintContainerType& rMsConstraintContainer);

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
    ApplyRANSChimeraProcessMonolithic& operator=(ApplyRANSChimeraProcessMonolithic const& rOther);

    ///@}
}; // Class ApplyRANSChimeraProcessMonolithic

} // namespace Kratos.

#endif //  KRATOS_APPLY_RANS_CHIMERA_MONOLITHIC_H_INCLUDED defined
