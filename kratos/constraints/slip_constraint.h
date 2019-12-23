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
 *
 */
class SlipConstraint
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
    // explicit SlipConstraint(IndexType Id = 0)
    //     : BaseType(Id)
    // {
    // }

    /**
     * @brief Constructor by passing a vector of Master and slave dofs and corresponding Matrix and constant vector
     * @param IndexType The Id of the new created constraint
     * @param rMasterDofsVector The vector containing the DoF of the master side
     * @param rSlaveDofsVector The vector containing the DoF of the slave side
     * @param rRelationMatrix The relation matrix between the master/slave DoF
     * @param rConstantVector The vector containing the additional kinematic relationship
     */
    SlipConstraint(
        IndexType Id,
        Node<3>& rNode,
        const VariableComponentType& rVarX,
        const VariableComponentType& rVarY,
        const VariableComponentType& rVarZ,
        const Variable<array_1d<double,3> >& rNormalVar
        ) : BaseType(Id),
            mrNode(rNode),
            mpNormalVar(&rNormalVar)
    {
        mRelationMatrix.resize(1,2,false);
        mRelationMatrix.clear();

        mConstantVector.resize(1,false);
        mConstantVector.clear();

        DofPointerVectorType all_dofs;
        all_dofs.reserve(3);
        all_dofs.push_back(rNode.pGetDof(rVarX));
        all_dofs.push_back(rNode.pGetDof(rVarY));
        all_dofs.push_back(rNode.pGetDof(rVarZ));

        auto n = rNode.FastGetSolutionStepValue(rNormalVar);
        n /= norm_2(n);

        array_1d<double,3> v;
        v[0] = all_dofs[0]->GetSolutionStepValue(rVarX);
        v[1] = all_dofs[1]->GetSolutionStepValue(rVarY);
        v[2] = all_dofs[2]->GetSolutionStepValue(rVarZ);

        unsigned int max_n_component_index = 0;
        double max_abs_n_component = std::abs(n[0]);
        for(unsigned int i=1; i<3; ++i)
        {
            if(std::abs(n[i]) > max_abs_n_component)
            {
                max_abs_n_component=std::abs(n[i]);
                max_n_component_index = i;
            }
        }

        mMasterDofsVector.clear();
        mSlaveDofsVector.clear();

        unsigned int counter = 0;
        unsigned int Tcounter = 0;
        for(auto& dof : all_dofs)
        {
            if(counter == max_n_component_index)
            {
                mSlaveDofsVector.push_back(dof);
            }
            else
            {
                mMasterDofsVector.push_back(dof);
                mRelationMatrix(0,Tcounter++) = -n[counter]/n[max_n_component_index];
            }

            counter++;
        }
    }

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
        ) override
    {
        KRATOS_ERROR << "Slip constraint doesn't allow SetDofList" << std::endl;
    }

    /**
     * @brief This method returns the slave dof vector
     * @return The vector containing the slave dofs
     */
    void SetSlaveDofsVector(const DofPointerVectorType& rSlaveDofsVector) override
    {
        KRATOS_ERROR << "Slip constraint doesn't allow SetSlaveDofsVector" << std::endl;
    }

    /**
     * @brief This method returns the slave dof vector
     * @return The vector containing the slave dofs
     */
    void SetMasterDofsVector(const DofPointerVectorType& rMasterDofsVector) override
    {
        KRATOS_ERROR << "Slip constraint doesn't allow SetMasterDofsVector" << std::endl;
    }


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
        ) override
    {
        KRATOS_ERROR << "Slip constraint doesn't allow SetLocalSystem" << std::endl;
    }


    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Returns the string containing a detailed description of this object.
     * @return the string with informations
     */
    std::string GetInfo() const override
    {
        return "SlipConstraint class !";
    }

    /**
     * @brief This method prints the current Constraint Id
     * @param rOStream The buffer where the information is given
     */
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << " SlipConstraint Id  : " << this->Id() << std::endl;
        rOStream << " node Id " << mrNode.Id() << std::endl;
        rOStream << " normal : " << mrNode.FastGetSolutionStepValue(*mpNormalVar) << std::endl;
        rOStream << " Number of Slaves          : " << mSlaveDofsVector.size() << std::endl;
        rOStream << " Number of Masters         : " << mMasterDofsVector.size() << std::endl;
    }

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{
    Node<3>& mrNode;
    const Variable<array_1d<double,3> >* mpNormalVar;

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

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LinearMasterSlaveConstraint);
        rSerializer.save("mrNode", mrNode);
        rSerializer.save("mpNormalVar", mpNormalVar);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LinearMasterSlaveConstraint);
        rSerializer.load("mrNode", mrNode);
        rSerializer.load("mpNormalVar", mpNormalVar);
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
