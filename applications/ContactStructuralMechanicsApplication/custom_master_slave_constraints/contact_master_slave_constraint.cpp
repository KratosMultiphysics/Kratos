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

// System includes

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "custom_master_slave_constraints/contact_master_slave_constraint.h"

namespace Kratos
{
ContactMasterSlaveConstraint::ContactMasterSlaveConstraint(IndexType Id)
    : BaseType(Id)
{
}

/***********************************************************************************/
/***********************************************************************************/

ContactMasterSlaveConstraint::ContactMasterSlaveConstraint(
    IndexType Id,
    DofPointerVectorType& rMasterDofsVector,
    DofPointerVectorType& rSlaveDofsVector,
    const MatrixType& rRelationMatrix,
    const VectorType& rConstantVector
    ) : BaseType(Id, rMasterDofsVector, rSlaveDofsVector, rRelationMatrix, rConstantVector)
{
}

/***********************************************************************************/
/***********************************************************************************/

ContactMasterSlaveConstraint::ContactMasterSlaveConstraint(
    IndexType Id,
    NodeType& rMasterNode,
    const VariableType& rMasterVariable,
    NodeType& rSlaveNode,
    const VariableType& rSlaveVariable,
    const double Weight,
    const double Constant
    ) : BaseType(Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant)
{
    KRATOS_ERROR << "ContactMasterSlaveConstraint :: Please don't use this constructor. A components variable is expected" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

ContactMasterSlaveConstraint::ContactMasterSlaveConstraint(
    IndexType Id,
    NodeType& rMasterNode,
    const VariableComponentType& rMasterVariable,
    NodeType& rSlaveNode,
    const VariableComponentType& rSlaveVariable,
    const double Weight,
    const double Constant
    )
{
    KRATOS_ERROR << "ContactMasterSlaveConstraint :: Please don't use this constructor. A components variable is expected" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

ContactMasterSlaveConstraint::~ContactMasterSlaveConstraint()
{

}

/***********************************************************************************/
/***********************************************************************************/

ContactMasterSlaveConstraint::ContactMasterSlaveConstraint(const ContactMasterSlaveConstraint& rOther)
    : BaseType(rOther)
{
}

/***********************************************************************************/
/***********************************************************************************/

ContactMasterSlaveConstraint& ContactMasterSlaveConstraint::operator=(const ContactMasterSlaveConstraint& rOther)
{
    BaseType::operator=( rOther );
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

MasterSlaveConstraint::Pointer ContactMasterSlaveConstraint::Create(
    IndexType Id,
    DofPointerVectorType& rMasterDofsVector,
    DofPointerVectorType& rSlaveDofsVector,
    const MatrixType& rRelationMatrix,
    const VectorType& rConstantVector
    ) const
{
    return Kratos::make_shared<ContactMasterSlaveConstraint>(Id, rMasterDofsVector, rSlaveDofsVector, rRelationMatrix, rConstantVector);
}

/***********************************************************************************/
/***********************************************************************************/

MasterSlaveConstraint::Pointer ContactMasterSlaveConstraint::Create(
    IndexType Id,
    NodeType& rMasterNode,
    const VariableType& rMasterVariable,
    NodeType& rSlaveNode,
    const VariableType& rSlaveVariable,
    const double Weight,
    const double Constant
    ) const
{
    return Kratos::make_shared<ContactMasterSlaveConstraint>(Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
}

/***********************************************************************************/
/***********************************************************************************/

MasterSlaveConstraint::Pointer ContactMasterSlaveConstraint::Create(
    IndexType Id,
    NodeType& rMasterNode,
    const VariableComponentType& rMasterVariable,
    NodeType& rSlaveNode,
    const VariableComponentType& rSlaveVariable,
    const double Weight,
    const double Constant
    ) const
{
    return Kratos::make_shared<ContactMasterSlaveConstraint>(Id, rMasterNode, rMasterVariable, rSlaveNode, rSlaveVariable, Weight, Constant);
}

/***********************************************************************************/
/***********************************************************************************/

void ContactMasterSlaveConstraint::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{

}

/***********************************************************************************/
/***********************************************************************************/

std::string ContactMasterSlaveConstraint::GetInfo() const
{
    return "This is contact MPC !";
}

/***********************************************************************************/
/***********************************************************************************/

void ContactMasterSlaveConstraint::PrintInfo(std::ostream &rOStream) const
{
    rOStream << " ContactMasterSlaveConstraint Id  : " << this->Id() << std::endl;
    rOStream << " Number of Slaves          : " << BaseType::mSlaveDofsVector.size() << std::endl;
    rOStream << " Number of Masters         : " << BaseType::mMasterDofsVector.size() << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void ContactMasterSlaveConstraint::save(Serializer &rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType);
}

/***********************************************************************************/
/***********************************************************************************/

void ContactMasterSlaveConstraint::load(Serializer &rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType);
}

} // namespace Kratos
