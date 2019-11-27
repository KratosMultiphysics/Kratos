// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(CONTACT_MASTER_SLAVE_CONSTRAINT_H)
#define CONTACT_MASTER_SLAVE_CONSTRAINT_H

// System includes

// External includes

// Project includes
#include "includes/linear_master_slave_constraint.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef std::size_t SizeType;

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
 * @class ContactMasterSlaveConstraint
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This is a constraint for contact mechanics based in a linear kinematic MPC constriant
 * @details This constraint is based on LinearMasterSlaveConstraint. It adds additional consideration related with contact
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ContactMasterSlaveConstraint
    :  public LinearMasterSlaveConstraint
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the base constraint class
    typedef MasterSlaveConstraint BaseConstraintType;

    /// The definition of the base class, we take the rest of the definitions from the base class
    typedef LinearMasterSlaveConstraint BaseType;

    /// The index type definition
    typedef BaseType::IndexType IndexType;

    /// The DoF type definition
    typedef BaseType::DofType DofType;

    /// The DoF pointer vector type definition
    typedef BaseType::DofPointerVectorType DofPointerVectorType;

    /// The node type definition
    typedef BaseType::NodeType NodeType;

    /// The equation Id vector type definition
    typedef BaseType::EquationIdVectorType EquationIdVectorType;

    /// The matrix type definition
    typedef BaseType::MatrixType MatrixType;

    /// The vector type definition
    typedef BaseType::VectorType VectorType;

    /// The variable type definition (double)
    typedef BaseType::VariableType VariableType;

    /// The component variable type definition
    typedef BaseType::VariableComponentType VariableComponentType;

    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(ContactMasterSlaveConstraint);

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The default constructor
     * @param IndexType The Id of the new created constraint
     */
    explicit ContactMasterSlaveConstraint(IndexType Id = 0);

    /**
     * @brief Constructor by passing a vector of Master and slave dofs and corresponding Matrix and constant vector
     * @param IndexType The Id of the new created constraint
     * @param rMasterDofsVector The vector containing the DoF of the master side
     * @param rSlaveDofsVector The vector containing the DoF of the slave side
     * @param rRelationMatrix The relation matrix between the master/slave DoF
     * @param rConstantVector The vector containing the additional kinematic relationship
     */
    ContactMasterSlaveConstraint(
        IndexType Id,
        DofPointerVectorType& rMasterDofsVector,
        DofPointerVectorType& rSlaveDofsVector,
        const MatrixType& rRelationMatrix,
        const VectorType& rConstantVector
        );

    /**
     * @brief Constructor by passing a single Master and slave dofs and corresponding weight and constant for a variable component
     * @param IndexType The Id of the new created constraint
     * @param rMasterNode The node of master side
     * @param rMasterVariable The variable of the master DoF
     * @param rSlaveNode The node of slave side
     * @param rSlaveVariable The variable of the slave DoF
     * @param Weight The relation between the master/slave DoF
     * @param Constant The additional kinematic relationship
     */
    ContactMasterSlaveConstraint(
        IndexType Id,
        NodeType& rMasterNode,
        const VariableType& rMasterVariable,
        NodeType& rSlaveNode,
        const VariableType& rSlaveVariable,
        const double Weight,
        const double Constant
        );

    /**
     * @brief Constructor by passing a single Master and slave dofs and corresponding weight and constant for a variable component
     * @param IndexType The Id of the new created constraint
     * @param rMasterNode The node of master side
     * @param rMasterVariable The variable of the master DoF
     * @param rSlaveNode The node of slave side
     * @param rSlaveVariable The variable of the slave DoF
     * @param Weight The relation between the master/slave DoF
     * @param Constant The additional kinematic relationship
     */
    ContactMasterSlaveConstraint(
        IndexType Id,
        NodeType& rMasterNode,
        const VariableComponentType& rMasterVariable,
        NodeType& rSlaveNode,
        const VariableComponentType& rSlaveVariable,
        const double Weight,
        const double Constant
        );

    /// Destructor.
    ~ContactMasterSlaveConstraint() override;

    /// Copy Constructor
    ContactMasterSlaveConstraint(const ContactMasterSlaveConstraint& rOther);

    /// Assignment operator
    ContactMasterSlaveConstraint& operator=(const ContactMasterSlaveConstraint& rOther);

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method by passing a single Master and slave dofs and corresponding weight and constant for a variable component
     * @param IndexType The Id of the new created constraint
     * @param rMasterDofsVector The DoFs of master side
     * @param rSlaveDofsVector The DoFs of slave side
     * @param rRelationMatrix The relation matrix between the master/slave DoF
     * @param rConstantVector The vector containing the additional kinematic relationship
     * @return A Pointer to the new constraint
     */
    MasterSlaveConstraint::Pointer Create(
        IndexType Id,
        DofPointerVectorType& rMasterDofsVector,
        DofPointerVectorType& rSlaveDofsVector,
        const MatrixType& rRelationMatrix,
        const VectorType& rConstantVector
        ) const override;

    /**
     * @brief Create method by passing a single Master and slave dofs and corresponding weight and constant for a variable component
     * @param IndexType The Id of the new created constraint
     * @param rMasterNode The node of master side
     * @param rMasterVariable The variable of the master DoF
     * @param rSlaveNode The node of slave side
     * @param rSlaveVariable The variable of the slave DoF
     * @param Weight The relation between the master/slave DoF
     * @param Constant The additional kinematic relationship
     * @return A Pointer to the new constraint
     */
    MasterSlaveConstraint::Pointer Create(
        IndexType Id,
        NodeType& rMasterNode,
        const VariableType& rMasterVariable,
        NodeType& rSlaveNode,
        const VariableType& rSlaveVariable,
        const double Weight,
        const double Constant
        ) const override;

    /**
     * @brief Create method by passing a single Master and slave dofs and corresponding weight and constant for a variable component
     * @param IndexType The Id of the new created constraint
     * @param rMasterNode The node of master side
     * @param rMasterVariable The variable of the master DoF
     * @param rSlaveNode The node of slave side
     * @param rSlaveVariable The variable of the slave DoF
     * @param Weight The relation between the master/slave DoF
     * @param Constant The additional kinematic relationship
     * @return A Pointer to the new constraint
     */
    MasterSlaveConstraint::Pointer Create(
        IndexType Id,
        NodeType& rMasterNode,
        const VariableComponentType& rMasterVariable,
        NodeType& rSlaveNode,
        const VariableComponentType& rSlaveVariable,
        const double Weight,
        const double Constant
        ) const override;

    /**
     * @brief This is called for non-linear analysis at the end of the iteration process
     */
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Returns the string containing a detailed description of this object.
     * @return the string with informations
     */
    std::string GetInfo() const override;

    /**
     * @brief This method prints the current Constraint Id
     * @param rOStream The buffer where the information is given
     */
    void PrintInfo(std::ostream &rOStream) const override;

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer &rSerializer) const override;

    void load(Serializer &rSerializer) override;
};

///@name Input/Output funcitons
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, ContactMasterSlaveConstraint& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const ContactMasterSlaveConstraint& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;

    return rOStream;
}

///@}


} // namespace Kratos

#endif // USER_PROVIDED_CONTACT_MASTER_SLAVE_CONSTRAINT_H
