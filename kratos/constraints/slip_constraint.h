//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(SLIP_CONSTRAINT_H)
#define SLIP_CONSTRAINT_H

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/linear_master_slave_constraint.h"

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

/**
 * @class SlipConstraint
 * @ingroup KratosCore
 * @brief This constructs a constraint which imposes that v * n = 0
 * where:
 *  v is a variable
 *  n is the unit normal to the node
 * @author Riccardo Rossi
 */
class KRATOS_API(KRATOS_CORE) SlipConstraint
    :  public LinearMasterSlaveConstraint
{
public:
    ///@name Type Definitions
    ///@{

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
    KRATOS_CLASS_POINTER_DEFINITION(SlipConstraint);

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

    /**
     * @brief Constructor2D
     */    
    SlipConstraint(
        IndexType Id,
        Dof<double>* pDofX,
        Dof<double>* pDofY,
        const array_1d<double,3> NormalVector
    );

    /**
     * @brief Constructor3D
     */
    SlipConstraint(
        IndexType Id,
        Dof<double>* pDofX,
        Dof<double>* pDofY,
        Dof<double>* pDofZ,
        const array_1d<double,3> NormalVector
    );

    /// Destructor.
    ~SlipConstraint() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Determines the constrant's slave and master list of DOFs
     * @param rSlaveDofsVector The list of slave DOFs
     * @param rMasterDofsVector The list of slave DOFs
     * @param rCurrentProcessInfo The current process info instance
     */
    void SetDofList(
        const DofPointerVectorType& rSlaveDofsVector,
        const DofPointerVectorType& rMasterDofsVector,
        const ProcessInfo& rCurrentProcessInfo
    ) override ;

    /**
     * @brief This method returns the slave dof vector
     * @return The vector containing the slave dofs
     */
    void SetSlaveDofsVector(const DofPointerVectorType& rSlaveDofsVector) override ;

    /**
     * @brief This method returns the slave dof vector
     * @return The vector containing the slave dofs
     */
    void SetMasterDofsVector(const DofPointerVectorType& rMasterDofsVector) override ;


    /**
     * @brief This method allows to set the Local System in case is not computed on tunning time (internal variable)
     * @param rRelationMatrix the matrix which relates the master and slave degree of freedom
     * @param rConstant The constant vector (one entry for each slave)
     * @param rCurrentProcessInfo The current process info instance
     */
    void SetLocalSystem(
        const MatrixType& rRelationMatrix,
        const VectorType& rConstantVector,
        const ProcessInfo& rCurrentProcessInfo
    ) override ;


    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Returns the string containing a detailed description of this object.
     * @return the string with informations
     */
    std::string GetInfo() const override ;

    /**
     * @brief This method prints the current Constraint Id
     * @param rOStream The buffer where the information is given
     */
    void PrintInfo(std::ostream &rOStream) const override;

    ///@}
protected:        
    void ConstructorHelper(
        DofPointerVectorType& rAllDofs,
        array_1d<double,3>& rNormalVector
        );

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

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearMasterSlaveConstraint);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearMasterSlaveConstraint);
    }
};

///@name Input/Output funcitons
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, SlipConstraint& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const SlipConstraint& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;

    return rOStream;
}

///@}


} // namespace Kratos

#endif // SLIP_CONSTRAINT_H
