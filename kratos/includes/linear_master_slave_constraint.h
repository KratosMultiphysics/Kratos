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

#if !defined(LINEAR_MASTER_SLAVE_CONSTRAINT_H)
#define LINEAR_MASTER_SLAVE_CONSTRAINT_H
// System includes

// project includes
#include "includes/define.h"
#include "includes/master_slave_constraint.h"

namespace Kratos
{
  /** @class LinearMasterSlaveConstraint
    * @ingroup KratosCore
    * @brief
    *
    * This class allows to add a master-slave constraint which is of the form
    *
    * SlaveDofVector = T * MasterDofVector + ConstantVector.
    *
    * or
    *
    * SlaveDof = weight * MasterDof + Constant
    *
    * the data T and ConstantVector (or the equivalent scalars) are not stored in the base class, since they can be eventually evaluated runtime.
    */
class LinearMasterSlaveConstraint :  public MasterSlaveConstraint
{
public:
    /// Pointer definition of DataValueContainer
    KRATOS_CLASS_POINTER_DEFINITION(LinearMasterSlaveConstraint);
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

    explicit LinearMasterSlaveConstraint(IndexType Id = 0)
        :
        MasterSlaveConstraint(Id)
    {
    }

    /*
    * Constructor by passing a vector of Master and slave dofs and corresponding Matrix and constant vector
    */
    LinearMasterSlaveConstraint(IndexType Id,
                                DofPointerVectorType& rMasterDofsVector,
                                DofPointerVectorType& rSlaveDofsVector,
                                const MatrixType& rRelationMatrix,
                                const VectorType& rConstantVector)
        : MasterSlaveConstraint(Id), mSlaveDofsVector(rSlaveDofsVector), mMasterDofsVector(rMasterDofsVector),
                                    mRelationMatrix(rRelationMatrix), mConstantVector(rConstantVector)
    {
    }

    /*
    * Constructor by passing a single Master and slave dofs and corresponding weight and constant for a variable component
    */
    LinearMasterSlaveConstraint(IndexType Id,
                                NodeType& rMasterNode,
                                const VariableType& rMasterVariable,
                                NodeType& rSlaveNode,
                                const VariableType& rSlaveVariable,
                                const double Weight,
                                const double Constant)
        : MasterSlaveConstraint(Id)
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
    LinearMasterSlaveConstraint(IndexType Id,
                                NodeType& rMasterNode,
                                const VariableComponentType& rMasterVariable,
                                NodeType& rSlaveNode,
                                const VariableComponentType& rSlaveVariable,
                                const double Weight,
                                const double Constant)
        : MasterSlaveConstraint(Id)
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
    ~LinearMasterSlaveConstraint() override
    {

    }

    /// Copy Constructor
    LinearMasterSlaveConstraint(const LinearMasterSlaveConstraint& rOther)
    {
        this->SetId(rOther.Id());
        // this->Flags = rOther.Flags;
        this->mSlaveDofsVector = rOther.mSlaveDofsVector;
        this->mMasterDofsVector = rOther.mMasterDofsVector;
        this->mRelationMatrix = rOther.mRelationMatrix;
        this->mConstantVector = rOther.mConstantVector;
    }

    /// Assignment operator
    LinearMasterSlaveConstraint& operator=(const LinearMasterSlaveConstraint& rOther)
    {
        this->SetId(rOther.Id());
        // this->Flags = rOther.Flags;
        this->mSlaveDofsVector = rOther.mSlaveDofsVector;
        this->mMasterDofsVector = rOther.mMasterDofsVector;
        this->mRelationMatrix = rOther.mRelationMatrix;
        this->mConstantVector = rOther.mConstantVector;

        return *this;
    }

    MasterSlaveConstraint::Pointer Create(IndexType Id,
                                          DofPointerVectorType& MasterDofsVector,
                                          DofPointerVectorType& SlaveDofsVector,
                                          const MatrixType& RelationMatrix,
                                          const VectorType& ConstantVector) const override
    {
        KRATOS_TRY
        auto new_pointer = Kratos::make_shared<LinearMasterSlaveConstraint>(Id, MasterDofsVector, SlaveDofsVector, RelationMatrix, ConstantVector);
        return new_pointer;
        KRATOS_CATCH("");
    }

    MasterSlaveConstraint::Pointer Create(IndexType Id,
                                          NodeType& rMasterNode,
                                          const VariableType& rMasterVariable,
                                          NodeType& rSlaveNode,
                                          const VariableType& rSlaveVariable,
                                          const double Weight,
                                          const double Constant) const override
    {
        KRATOS_TRY
        return Kratos::make_shared<LinearMasterSlaveConstraint>(Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        KRATOS_CATCH("");
    }

    MasterSlaveConstraint::Pointer Create(IndexType Id,
                                          NodeType& rMasterNode,
                                          const VariableComponentType& rMasterVariable,
                                          NodeType& rSlaveNode,
                                          const VariableComponentType& rSlaveVariable,
                                          const double Weight,
                                          const double Constant) const override
    {
        KRATOS_TRY
        return Kratos::make_shared<LinearMasterSlaveConstraint>(Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
        KRATOS_CATCH("");
    }


    ///@}

    ///@name Access
    ///@{

    virtual void GetDofList(DofPointerVectorType& rSlaveDofsVector,
                            DofPointerVectorType& rMasterDofsVector,
                            ProcessInfo& rCurrentProcessInfo) override
    {
        rSlaveDofsVector = mSlaveDofsVector;
        rMasterDofsVector = mMasterDofsVector;
    }

    void EquationIdVector(EquationIdVectorType& rSlaveEquationIds,
                          EquationIdVectorType& rMasterEquationIds,
                          ProcessInfo& rCurrentProcessInfo) override
    {
        if (rSlaveEquationIds.size() != mSlaveDofsVector.size())
            rSlaveEquationIds.resize(mSlaveDofsVector.size());

        if (rMasterEquationIds.size() != mMasterDofsVector.size())
            rMasterEquationIds.resize(mMasterDofsVector.size());

        for(unsigned int i=0; i<rSlaveEquationIds.size(); ++i)
            rSlaveEquationIds[i] = mSlaveDofsVector[i]->EquationId();

        for(unsigned int i=0; i<rMasterEquationIds.size(); ++i)
            rMasterEquationIds[i] = mMasterDofsVector[i]->EquationId();
    }

    void CalculateLocalSystem(MatrixType& rTransformationMatrix,
                              VectorType& rConstantVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
      rTransformationMatrix = mRelationMatrix;
      rConstantVector = mConstantVector;
    }

    std::string GetInfo() const override
    {
        return "Linear User Provded Master Slave Constraint class !";
    }

    ///@

    ///@name Static Operations
    ///
    //@{

    ///@}
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << " LinearMasterSlaveConstraint Id  : " << this->Id() << std::endl;
        rOStream << " Number of Slaves          : " << this->mSlaveDofsVector.size() << std::endl;
        rOStream << " Number of Masters         : " << this->mMasterDofsVector.size() << std::endl;
    }



private:
    ///@}
    DofPointerVectorType mSlaveDofsVector;
    DofPointerVectorType mMasterDofsVector;
    MatrixType mRelationMatrix;
    VectorType mConstantVector;

    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, MasterSlaveConstraint);
        rSerializer.save("SlaveDofVec", mSlaveDofsVector);
        rSerializer.save("MasterDofVec", mMasterDofsVector);
        rSerializer.save("RelationMat", mRelationMatrix);
        rSerializer.save("ConstantVec", mConstantVector);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, MasterSlaveConstraint);
        rSerializer.load("SlaveDofVec", mSlaveDofsVector);
        rSerializer.load("MasterDofVec", mMasterDofsVector);
        rSerializer.load("RelationMat", mRelationMatrix);
        rSerializer.load("ConstantVec", mConstantVector);
    }
};

///@name Input/Output funcitons
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, LinearMasterSlaveConstraint& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const LinearMasterSlaveConstraint& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;

    return rOStream;
}

///@}


} // namespace Kratos

#endif // USER_PROVIDED_LINEAR_MASTER_SLAVE_CONSTRAINT_H
