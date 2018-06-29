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

#if !defined(CONSTRAINT_H)
#define CONSTRAINT_H
// System includes

// project includes
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/node.h"
#include "includes/kratos_flags.h"
#include "geometries/geometry.h"
#include "containers/flags.h"
#include "containers/variable.h"
#include "containers/variable_component.h"
#include "containers/vector_component_adaptor.h"
#include "containers/variable_data.h"
#include "includes/process_info.h"
namespace Kratos
{
/** \brief Constraint * @class MasterSlaveRelation
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
    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(MasterSlaveConstraint);
    typedef IndexedObject BaseType;
    typedef std::size_t IndexType;
    typedef Dof<double> DofType;
    typedef std::vector< DofType::Pointer > DofPointerVectorType;
    typedef Node<3> NodeType;
    typedef PointerVectorSet<NodeType, IndexedObject> NodesContainerType;
    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef double ConstantType;
    typedef Matrix MatrixType;
    typedef Vector VectorType;
    typedef Kratos::Variable<double> VariableType;
    typedef Kratos::VariableComponent<Kratos::VectorComponentAdaptor<Kratos::array_1d<double, 3>>> VariableComponentType;
    ///@name Life Cycle
    ///@{

   /*
    * Empty Constructor
    */
    MasterSlaveConstraint(IndexType Id = 0):IndexedObject(Id), Flags()
    {
    }

    /*
    * Constructor by passing a vector of Master and slave dofs and corresponding Matrix and constant vector
    */
    MasterSlaveConstraint(IndexType Id, DofPointerVectorType& rMasterDofsVector,
                                        DofPointerVectorType& rSlaveDofsVector,
                                        MatrixType& rRelationMatrix,
                                        VectorType& rConstantVector):IndexedObject(Id), Flags()
    {
        mSlaveDofsVector = rSlaveDofsVector;
        mMasterDofsVector = rMasterDofsVector;
        mRelationMatrix = rRelationMatrix;
        mConstantVector = rConstantVector;
    }

    /*
    * Constructor by passing a single Master and slave dofs and corresponding weight and constant for a double variable
    */
    MasterSlaveConstraint(IndexType Id, NodeType& rMasterNode,
                                        VariableType& rMasterVariable,
                                        NodeType& rSlaveNode,
                                        VariableType& rSlaveVariable,
                                        double Weight,
                                        double Constant):IndexedObject(Id), Flags()
    {
        // Resizing the memeber variables
        mRelationMatrix.resize(1,1,false);
        mConstantVector.resize(1,false);

        // Obtaining the dofs from the variables
        mSlaveDofsVector.push_back(rSlaveNode.pGetDof(rSlaveVariable));
        mMasterDofsVector.push_back(rMasterNode.pGetDof(rMasterVariable));

        mRelationMatrix(0,0) = Weight;
        mConstantVector(0) = Constant;

        // Setting the slave flag on the node
        rSlaveNode.Set(SLAVE);
    }

    /*
    * Constructor by passing a single Master and slave dofs and corresponding weight and constant for a variable component
    */
    MasterSlaveConstraint(IndexType Id, NodeType& rMasterNode,
                                        VariableComponentType& rMasterVariable,
                                        NodeType& rSlaveNode,
                                        VariableComponentType& rSlaveVariable,
                                        double Weight,
                                        double Constant):IndexedObject(Id), Flags()
    {
        // Resizing the memeber variables
        mRelationMatrix.resize(1,1,false);
        mConstantVector.resize(1,false);

        // Obtaining the dofs from the variables
        mSlaveDofsVector.push_back(rSlaveNode.pGetDof(rSlaveVariable));
        mMasterDofsVector.push_back(rMasterNode.pGetDof(rMasterVariable));

        mRelationMatrix(0,0) = Weight;
        mConstantVector(0) = Constant;

        // Setting the slave flag on the node
        rSlaveNode.Set(SLAVE);
    }

    /// Destructor.
    virtual ~MasterSlaveConstraint() override
    {

    }

    /// Copy Constructor
    MasterSlaveConstraint(const MasterSlaveConstraint& rOther)
    {
        this->mSlaveDofsVector = rOther.mSlaveDofsVector;
        this->mMasterDofsVector = rOther.mMasterDofsVector;

        this->mRelationMatrix = rOther.mRelationMatrix;
        this->mConstantVector = rOther.mConstantVector;
    }

    /**
     * creates a new constraint pointer
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    virtual Pointer Create(IndexType Id, DofPointerVectorType& MasterDofsVector, DofPointerVectorType& SlaveDofsVector, MatrixType RelationMatrix, VectorType ConstantVector) const
    {
        KRATOS_TRY
        auto new_pointer = Kratos::make_shared<MasterSlaveConstraint>(Id, MasterDofsVector, SlaveDofsVector, RelationMatrix, ConstantVector);
        return new_pointer;
        KRATOS_CATCH("");
    }

    virtual Pointer Create(IndexType Id, NodeType& rMasterNode,
                                        VariableType& rMasterVariable,
                                        NodeType& rSlaveNode,
                                        VariableType& rSlaveVariable,
                                        double Weight,
                                        double Constant) const
    {
        KRATOS_TRY
        auto new_pointer = Kratos::make_shared<MasterSlaveConstraint>(Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        return new_pointer;
        KRATOS_CATCH("");
    }

    virtual Pointer Create(IndexType Id, NodeType& rMasterNode,
                                        VariableComponentType& rMasterVariable,
                                        NodeType& rSlaveNode,
                                        VariableComponentType& rSlaveVariable,
                                        double Weight,
                                        double Constant) const
    {
        KRATOS_TRY
        auto new_pointer = Kratos::make_shared<MasterSlaveConstraint>(Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        return new_pointer;
        KRATOS_CATCH("");
    }


    ///@}

    ///@name Access
    ///@{

    /**
	* Clears the maps contents
	*/
    void Clear()
    {
        //TODO: clear the relation matrix and the constant vector.
    }


    /**
     * is called to initialize the constraint
     * if the constraint needs to perform any operation before any calculation is done
     */
    virtual void Initialize()
    {
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
     * this determines the master equation IDs connected to this constraint
     * @param rResult the elemental equation ID vector
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
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rTransformationMatrix the elemental left hand side matrix
     * @param rConstant the elemental right hand side
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
        rOStream << " MasterSlaveConstraint Id  : " <<Id()<<std::endl;
        rOStream << " Number of Slaves          : " <<mSlaveDofsVector.size()<<std::endl;
        rOStream << " Number of Masters         : " <<mMasterDofsVector.size()<<std::endl;
    }

    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save(Serializer &rSerializer) const override
    {

    }

    virtual void load(Serializer &rSerializer) override
    {

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
    DofPointerVectorType mSlaveDofsVector;
    DofPointerVectorType mMasterDofsVector;
    MatrixType mRelationMatrix;
    VectorType mConstantVector;

    /**
     * pointer to the data related to this constraint
     */
    DataValueContainer mData;
    ///@}
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

#endif // CONSTRAINT_H_INCLUDED