//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//  Collaborators:   Vicente Mataix
//

#pragma once

// System includes

// project includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/flags.h"
#include "containers/variable.h"
#include "includes/process_info.h"
#include "includes/indexed_object.h"

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
 * @class MasterSlaveConstraint
 * @ingroup KratosCore
 * @brief A class that implements the interface for different master-slave constraints to be applied on a system.
 * @details This is the part that is seen by the user from the python level. Objects of this class are
 * first class citizens of the modelpart.
 *
 * This class allows to add a master-slave constraint which is of the form
 *
 * SlaveDofVector = T * MasterDofVector + rConstantVector. (Processing of this is currently not implemented.)
 *
 * or
 *
 * SlaveDof = weight * MasterDof + Constant
 *
 * This class's object will provide its slave, master details and relation matrix between them.
 *
 * One can add two MasterSlaveConstraint objects with same slave but different masters and weights.
 * Consider user adds : SlaveDof = weight1 * MasterDof1 + Constant1
 *              and   : SlaveDof = weight2 * MasterDof2 + Constant2
 *
 * These are later consolidated in the builder and solver to make
 *                    : SlaveDof = weight1 * MasterDof1 + weight2 * MasterDof2 + Constant1+Constant2
 *       and then converted to :
 *                    : SlaveEqID = weight1 * MasterEqId1 + weight2 * MasterEqId2 + Constant1+Constant2
 * This unique equation is used later on to modify the equation system.
 * @author Aditya Ghantasala
 */
class KRATOS_API(KRATOS_CORE) MasterSlaveConstraint
    :  public IndexedObject, public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the base class
    typedef IndexedObject BaseType;

    /// The index type definition
    typedef std::size_t IndexType;

    /// The DoF type definition
    typedef Dof<double> DofType;

    /// The DoF pointer vector type definition
    typedef std::vector< DofType::Pointer > DofPointerVectorType;

    /// The node type definition
    typedef Node NodeType;

    /// The equation Id vector type definition
    typedef std::vector<std::size_t> EquationIdVectorType;

    /// The matrix type definition
    typedef Matrix MatrixType;

    /// The vector type definition
    typedef Vector VectorType;

    /// The variable type definition (double)
    typedef Kratos::Variable<double> VariableType;

    /// Pointer definition of MasterSlaveConstraint
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveConstraint);

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
    explicit MasterSlaveConstraint(IndexType Id = 0) : IndexedObject(Id), Flags()
    {
    }

    /// Destructor.
    virtual ~MasterSlaveConstraint() override
    {

    }

    /// Copy Constructor
    MasterSlaveConstraint(const MasterSlaveConstraint& rOther)
        : BaseType(rOther),
          mData(rOther.mData)
    {
    }

    /// Assignment operator
    MasterSlaveConstraint& operator=(const MasterSlaveConstraint& rOther)
    {
        BaseType::operator=( rOther );
        mData = rOther.mData;
        return *this;
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new constraint pointer
     * @param Id the ID of the new constraint
     * @param rMasterDofsVector the vector of master degree of freedoms.
     * @param rSlaveDofsVector the vector of slave degree of freedoms.
     * @param rRelationMatrix The matrix of weights relating the master DOFs and Slave DOFs
     * @param rConstantVector The vector of the constants, one entry for each of the slave.
     * @return A Pointer to the new constraint
     */
    virtual MasterSlaveConstraint::Pointer Create(
        IndexType Id,
        DofPointerVectorType& rMasterDofsVector,
        DofPointerVectorType& rSlaveDofsVector,
        const MatrixType& rRelationMatrix,
        const VectorType& rConstantVector
        ) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Create not implemented in MasterSlaveConstraintBaseClass" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * creates a new constraint pointer
     * @param Id the ID of the new constraint
     * @param rMasterNode Node which is the master of for this constraint.
     * @param rMasterVariable the scalar variable which is on the master node. (DOF)
     * @param rSlaveNode Node which is the slave of for this constraint.
     * @param rSlaveVariable the scalar variable which is on the slave node. (DOF)
     * @param Weight The weight with which the master and slave are related s = w*m + c
     * @param Constant The constant in the master slave relation
     * @return A Pointer to the new constraint
     */
    virtual MasterSlaveConstraint::Pointer Create(
        IndexType Id,
        NodeType& rMasterNode,
        const VariableType& rMasterVariable,
        NodeType& rSlaveNode,
        const VariableType& rSlaveVariable,
        const double Weight,
        const double Constant
        ) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Create not implemented in MasterSlaveConstraintBaseClass" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief It creates a new constraint pointer and clones the previous constraint data
     * @param NewId the ID of the new constraint
     * @return a Pointer to the new constraint
     */
    virtual Pointer Clone (IndexType NewId) const
    {
        KRATOS_TRY

        KRATOS_WARNING("MasterSlaveConstraint") << " Call base class constraint Clone " << std::endl;
        MasterSlaveConstraint::Pointer p_new_const = Kratos::make_shared<MasterSlaveConstraint>(*this);
        p_new_const->SetId(NewId);
        p_new_const->SetData(this->GetData());
        p_new_const->Set(Flags(*this));
        return p_new_const;

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the maps contents
     */
    virtual void Clear()
    {
    }

    /**
     * @brief It is called to initialize the constraint
     * @details If the constraint needs to perform any operation before any calculation is done
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * @brief It is called to finalize the constraint
     * @details If the constraint needs to perform any operation before any calculation is done
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void Finalize(const ProcessInfo& rCurrentProcessInfo)
    {
        this->Clear();
    }

    /**
     * @brief This is called in the beginning of each solution step
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * @brief This is called for non-linear analysis at the end of the iteration process
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * @brief This is called at the end of each solution step
     */
    virtual void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * @brief Determines the constrant's slave and master list of DOFs
     * @param rSlaveDofsVector The list of slave DOFs
     * @param rMasterDofsVector The list of slave DOFs
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetDofList(
        DofPointerVectorType& rSlaveDofsVector,
        DofPointerVectorType& rMasterDofsVector,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        KRATOS_ERROR << "GetDofList not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * @brief Determines the constrant's slave and master list of DOFs
     * @param rSlaveDofsVector The list of slave DOFs
     * @param rMasterDofsVector The list of slave DOFs
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void SetDofList(
        const DofPointerVectorType& rSlaveDofsVector,
        const DofPointerVectorType& rMasterDofsVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_ERROR << "SetDofList not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * @brief This determines the master equation IDs connected to this constraint
     * @param rSlaveEquationIds The vector of slave equation ids.
     * @param rMasterEquationIds The vector of master equation ids.
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void EquationIdVector(
        EquationIdVectorType& rSlaveEquationIds,
        EquationIdVectorType& rMasterEquationIds,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        if (rSlaveEquationIds.size() != 0)
            rSlaveEquationIds.resize(0);

        if (rMasterEquationIds.size() != 0)
            rMasterEquationIds.resize(0);
    }

    /**
     * @brief This method returns the slave dof vector
     * @return The vector containing the slave dofs
     */
    virtual const DofPointerVectorType& GetSlaveDofsVector() const
    {
        KRATOS_ERROR << "GetSlaveDofsVector not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * @brief This method returns the slave dof vector
     * @return The vector containing the slave dofs
     */
    virtual void SetSlaveDofsVector(const DofPointerVectorType& rSlaveDofsVector)
    {
        KRATOS_ERROR << "SetSlaveDofsVector not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * @brief This method returns the slave dof vector
     * @return The vector containing the slave dofs
     */
    virtual const DofPointerVectorType& GetMasterDofsVector() const
    {
        KRATOS_ERROR << "GetMasterDofsVector not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * @brief This method returns the slave dof vector
     * @return The vector containing the slave dofs
     */
    virtual void SetMasterDofsVector(const DofPointerVectorType& rMasterDofsVector)
    {
        KRATOS_ERROR << "SetMasterDofsVector not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * @brief This method resets the values of the slave dofs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void ResetSlaveDofs(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR << "ResetSlaveDofs not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * @brief This method directly applies the master/slave relationship
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void Apply(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR << "Apply not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * @brief This method allows to set the Local System in case is not computed on running time (internal variable)
     * @param rRelationMatrix the matrix which relates the master and slave degree of freedom
     * @param rConstant The constant vector (one entry for each slave)
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void SetLocalSystem(
        const MatrixType& rRelationMatrix,
        const VectorType& rConstantVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        KRATOS_ERROR << "SetLocalSystem not implemented in MasterSlaveConstraintBaseClass" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method allows to get the Local System in case is not computed on running time (internal variable)
     * @param rRelationMatrix the matrix which relates the master and slave degree of freedom
     * @param rConstant The constant vector (one entry for each slave)
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetLocalSystem(
        MatrixType& rRelationMatrix,
        VectorType& rConstantVector,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        KRATOS_TRY

        this->CalculateLocalSystem(rRelationMatrix, rConstantVector, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    /**
     * @brief This is called during the assembling process in order
     * @details To calculate the relation between the master and slave.
     * @param rRelationMatrix the matrix which relates the master and slave degree of freedom
     * @param rConstant The constant vector (one entry for each slave)
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateLocalSystem(
        MatrixType& rRelationMatrix,
        VectorType& rConstantVector,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        if (rRelationMatrix.size1() != 0) {
            rRelationMatrix.resize(0, 0, false);
        }

        if (rConstantVector.size() != 0) {
            rConstantVector.resize(0, false);
        }
    }

    /**
     * @brief This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rCurrentProcessInfo
     * @note This method is: MANDATORY
     */
    virtual int Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR_IF( this->Id() < 1 ) << "MasterSlaveConstraint found with Id " << this->Id() << std::endl;

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Returns the string containing a detailed description of this object.
     * @return the string with information
     */
    virtual std::string GetInfo() const
    {
        return " Constraint base class !";
    }

    /**
     * @brief This method prints the current Constraint Id
     * @param rOStream The buffer where the information is given
     */
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << " MasterSlaveConstraint Id  : " << this->Id() << std::endl;
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the data container of the constraint
     * @return The data container mData
     */
    DataValueContainer& Data()
    {
        return mData;
    }

    /**
     * @brief This method returns the data container of the constraint (constant)
     * @return The data container mData
     */
    DataValueContainer const& GetData() const
    {
        return mData;
    }

    /**
     * @brief This method sets the data container of the constraint
     * @param rThisData The data container to set on mData
     */
    void SetData(DataValueContainer const& rThisData)
    {
        mData = rThisData;
    }

    /**
     * @brief Check if the Data exists with Has(..) methods:
     * @param rThisVariable The variable to be check
     */
    template<class TDataType> 
    bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    /**
     * @brief Set Data with SetValue and the Variable to set
     * @param rThisVariable The variable to be set
     * @param rValue The value to be set
     */
    template<class TVariableType>
    void SetValue(
        const TVariableType& rThisVariable,
        typename TVariableType::Type const& rValue
        )
    {
        mData.SetValue(rThisVariable, rValue);
    }

    /**
     * @brief Get Data with GetValue and the Variable to get
     * @param rThisVariable The variable to get
     */
    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }

    /**
     * @brief Get Data with GetValue and the Variable to get
     * @param rThisVariable The variable to get
     */
    template<class TVariableType>
    typename TVariableType::Type& GetValue(const TVariableType& rThisVariable) const
    {
        return mData.GetValue(rThisVariable);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if the GeometricalObject is active
     * @return True by default, otherwise depending on the ACTIVE flag
     */
    bool IsActive() const;

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

    DataValueContainer mData; /// Pointer to the data related to this constraint

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject);
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags);
        rSerializer.save("Data", mData);
    }

    virtual void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject);
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags);
        rSerializer.load("Data", mData);
    }
};

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<MasterSlaveConstraint>;

///@name Input/Output functions
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                MasterSlaveConstraint& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                const MasterSlaveConstraint& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;

    return rOStream;
}

///@}

} // namespace Kratos
