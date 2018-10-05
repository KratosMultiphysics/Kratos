//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//
//

#if !defined(MASTER_SLAVE_CONSTRAINT_H)
#define MASTER_SLAVE_CONSTRAINT_H
// System includes

// project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "includes/process_info.h"
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
 * SlaveDofVector = T * MasterDofVector + ConstantVector. (Processing of this is currently not implemented.)
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
class MasterSlaveConstraint
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
    typedef Node<3> NodeType;

    /// The equation Id vector type definition
    typedef std::vector<std::size_t> EquationIdVectorType;

    /// The matrix type definition
    typedef Matrix MatrixType;

    /// The vector type definition
    typedef Vector VectorType;

    /// The variable type definition (double)
    typedef Kratos::Variable<double> VariableType;

    /// The component variable type definition
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;

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
     * @param MasterDofsVector the vector of master degree of freedoms.
     * @param SlaveDofsVector the vector of slave degree of freedoms.
     * @param RelationMatrix The matrix of weights relating the master DOFs and Slave DOFs
     * @param ConstantVector The vector of the constants, one entry for each of the slave.
     * @return A Pointer to the new constraint
     */
    virtual MasterSlaveConstraint::Pointer Create(
        IndexType Id,
        DofPointerVectorType& MasterDofsVector,
        DofPointerVectorType& SlaveDofsVector,
        const MatrixType& RelationMatrix,
        const VectorType& ConstantVector
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
     * @brief Creates a new constraint pointer
     * @param Id the ID of the new constraint
     * @param rMasterNode Node which is the master of for this constraint.
     * @param rMasterVariable the component of vector variable which is on the master node. (DOF)
     * @param rSlaveNode Node which is the slave of for this constraint.
     * @param rSlaveVariable the component of vector variable which is on the slave node. (DOF)
     * @param Weight The weight with which the master and slave are related s = w*m + c
     * @param Constant The constant in the master slave relation
     * @return A Pointer to the new constraint
     */
    virtual MasterSlaveConstraint::Pointer Create(
        IndexType Id,
        NodeType& rMasterNode,
        const VariableComponentType& rMasterVariable,
        NodeType& rSlaveNode,
        const VariableComponentType& rSlaveVariable,
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
        KRATOS_ERROR << "Clone not implemented in MasterSlaveConstraintBaseClass" << std::endl;
        return nullptr;
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
     */
    virtual void Initialize()
    {
    }

    /**
     * @brief It is called to finalize the constraint
     * @details If the constraint needs to perform any operation before any calculation is done
     */
    virtual void Finalize()
    {
        this->Clear();
    }

    /**
     * @brief This is called in the beginning of each solution step
     */
    virtual void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * @brief This is called for non-linear analysis at the beginning of the iteration process
     */
    virtual void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
    {
    }

    /**
     * @brief This is called for non-linear analysis at the end of the iteration process
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
     * @param rSlaveDofList The list of slave DOFs
     * @param rMasterDofList The list of slave DOFs
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetDofList(
        DofPointerVectorType& rSlaveDofList,
        DofPointerVectorType& rMasterDofList,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        KRATOS_ERROR << "Create not implemented in MasterSlaveConstraintBaseClass" << std::endl;
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
     * @brief This is called during the assembling process in order
     * @details To calculate the relation between the master and slave.
     * @param rTransformationMatrix the matrix which relates the master and slave degree of freedom
     * @param rConstant The constant vector (one entry for each slave)
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateLocalSystem(
        MatrixType& rTransformationMatrix,
        VectorType& rConstantVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        if (rTransformationMatrix.size1() != 0) {
            rTransformationMatrix.resize(0, 0, false);
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
    virtual int Check(const ProcessInfo& rCurrentProcessInfo)
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
     * @return the string with informations
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
    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    /**
     * @brief Set Data with SetValue and the Variable to set
     * @param rThisVariable The variable to be set
     * @param rValue The value to be set
     */
    template<class TVariableType> void SetValue(
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
    template<class TVariableType> typename TVariableType::Type& GetValue(const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }

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


///@name Input/Output funcitons
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

#endif // MASTER_SLAVE_CONSTRAINT_H
