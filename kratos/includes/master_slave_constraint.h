//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
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
/** * @class MasterSlaveConstraint
    * @ingroup KratosCore
    * @brief
	* A class that implements the interface for different master-slave constraints to be applied on a system.
    * This is the part that is seen by the user from the python level. Objects of this class are
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
    */
class MasterSlaveConstraint :  public IndexedObject, public Flags
{

  public:
    /// Pointer definition of MasterSlaveConstraint
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveConstraint);
    typedef IndexedObject BaseType;
    typedef std::size_t IndexType;
    typedef Dof<double> DofType;
    typedef std::vector< DofType::Pointer > DofPointerVectorType;
    typedef Node<3> NodeType;
    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef Matrix MatrixType;
    typedef Vector VectorType;
    typedef Kratos::Variable<double> VariableType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    ///@name Life Cycle
    ///@{


    /// Empty Constructor
    explicit MasterSlaveConstraint(IndexType Id = 0) : IndexedObject(Id), Flags()
    {
    }

    /// Destructor.
    virtual ~MasterSlaveConstraint() override
    {

    }

    /**
     * creates a new constraint pointer
     * @param Id the ID of the new constraint
     * @param MasterDofsVector the vector of master degree of freedoms.
     * @param SlaveDofsVector the vector of slave degree of freedoms.
     * @param RelationMatrix The matrix of weights relating the master DOFs and Slave DOFs
     * @param ConstantVector The vector of the constants, one entry for each of the slave.
     * @return a Pointer to the new constraint
     */
    virtual MasterSlaveConstraint::Pointer Create(IndexType Id,
                                                  DofPointerVectorType& MasterDofsVector,
                                                  DofPointerVectorType& SlaveDofsVector,
                                                  const MatrixType& RelationMatrix,
                                                  const VectorType& ConstantVector) const
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
     * @return a Pointer to the new constraint
     */
    virtual MasterSlaveConstraint::Pointer Create(IndexType Id,
                                                  NodeType& rMasterNode,
                                                  const VariableType& rMasterVariable,
                                                  NodeType& rSlaveNode,
                                                  const VariableType& rSlaveVariable,
                                                  const double Weight,
                                                  const double Constant) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Create not implemented in MasterSlaveConstraintBaseClass" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * creates a new constraint pointer
     * @param Id the ID of the new constraint
     * @param rMasterNode Node which is the master of for this constraint.
     * @param rMasterVariable the component of vector variable which is on the master node. (DOF)
     * @param rSlaveNode Node which is the slave of for this constraint.
     * @param rSlaveVariable the component of vector variable which is on the slave node. (DOF)
     * @param Weight The weight with which the master and slave are related s = w*m + c
     * @param Constant The constant in the master slave relation
     * @return a Pointer to the new constraint
     */
    virtual MasterSlaveConstraint::Pointer Create(IndexType Id,
                                                  NodeType& rMasterNode,
                                                  const VariableComponentType& rMasterVariable,
                                                  NodeType& rSlaveNode,
                                                  const VariableComponentType& rSlaveVariable,
                                                  const double Weight,
                                                  const double Constant) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Create not implemented in MasterSlaveConstraintBaseClass" << std::endl;

        KRATOS_CATCH("");
    }


    ///@}

    ///@name Access
    ///@{

    /**
	* Clears the maps contents
	*/
    virtual void Clear()
    {
    }

    /**
     * is called to initialize the constraint
     * if the constraint needs to perform any operation before any calculation is done
     */
    virtual void Initialize()
    {
    }

    /**
     * is called to finalize the constraint
     * if the constraint needs to perform any operation before any calculation is done
     */
    virtual void Finalize()
    {
        this->Clear();
    }

    /**
     * this is called in the beginning of each solution step
     */
    virtual void InitializeSolutionStep()
    {
    }


    /**
     * this is called for non-linear analysis at the beginning of the iteration process
     */
    virtual void InitializeNonLinearIteration()
    {
    }

    /**
     * this is called for non-linear analysis at the end of the iteration process
     */
    virtual void FinalizeNonLinearIteration()
    {
    }

    /**
     * this is called at the end of each solution step
     */
    virtual void FinalizeSolutionStep()
    {
    }

    /**
     * determines the constrant's slvae and master list of DOFs
     * @param rSlaveDofList the list of slave DOFs
     * @param rMasterDofList the list of slave DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void GetDofList(DofPointerVectorType& rSlaveDofList,
                            DofPointerVectorType& rMasterDofList,
                            ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR << "Create not implemented in MasterSlaveConstraintBaseClass" << std::endl;
    }

    /**
     * this determines the master equation IDs connected to this constraint
     * @param rSlaveEquationIds the vector of slave equation ids.
     * @param rMasterEquationIds the vector of master equation ids.
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void EquationIdVector(EquationIdVectorType& rSlaveEquationIds,
                                  EquationIdVectorType& rMasterEquationIds,
                                  ProcessInfo& rCurrentProcessInfo)
    {
        if (rSlaveEquationIds.size() != 0)
            rSlaveEquationIds.resize(0);

        if (rMasterEquationIds.size() != 0)
            rMasterEquationIds.resize(0);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the relation between the master and slave.
     * matrix and the right hand side
     * @param rTransformationMatrix the matrix which relates the master and slave degree of freedom
     * @param rConstant The constant vector (one entry for each slave)
     * @param rCurrentProcessInfo the current process info instance
     */
    virtual void CalculateLocalSystem(MatrixType& rTransformationMatrix,
                                      VectorType& rConstantVector,
                                      ProcessInfo& rCurrentProcessInfo)
    {
      if (rTransformationMatrix.size1() != 0)
      {
    	rTransformationMatrix.resize(0, 0, false);
      }

      if (rConstantVector.size() != 0)
      {
    	rConstantVector.resize(0, false);
      }
    }

    /**
	* Returns the string containing a detailed description of this object.
	* @return the string with informations
	*/
    virtual std::string GetInfo() const
    {
        return " Constraint base class !";
    }

    ///@

    ///@name Static Operations
    ///
    //@{

    ///@}
    virtual void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << " MasterSlaveConstraint Id  : " <<this->Id()<<std::endl;
    }

    /**
     * Check if the Data exists with Has(..) methods:
     */
    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return mData.Has(rThisVariable);
    }

    /**
     * Set Data with SetValue and the Variable to set:
     */
    template<class TVariableType> void SetValue(
        const TVariableType& rThisVariable,
        typename TVariableType::Type const& rValue)
    {
        mData.SetValue(rThisVariable, rValue);
    }

    /**
     * Get Data with GetValue and the Variable to get:
     */
    template<class TVariableType> typename TVariableType::Type& GetValue(
        const TVariableType& rThisVariable)
    {
        return mData.GetValue(rThisVariable);
    }


  private:
    ///@}

    /**
     * pointer to the data related to this constraint
     */
    DataValueContainer mData;
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