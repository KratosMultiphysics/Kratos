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

// System includes


// External includes


// Project includes
#include "constraints/slip_constraint.h"

namespace Kratos
{

SlipConstraint::SlipConstraint(
    IndexType Id,
    Node<3>& rNode,
    const VariableComponentType& rVarX,
    const VariableComponentType& rVarY,
    const VariableComponentType& rVarZ,
    array_1d<double,3> NormalVector, //note that a copy is made
    const unsigned int dim
    ) : BaseType(Id)
{
    mRelationMatrix.resize(1,dim-1,false);
    mRelationMatrix.clear();

    mConstantVector.resize(1,false);
    mConstantVector.clear();

    DofPointerVectorType all_dofs;
    all_dofs.reserve(dim);
    all_dofs.push_back(rNode.pGetDof(rVarX));
    all_dofs.push_back(rNode.pGetDof(rVarY));
    if(dim == 3)
       all_dofs.push_back(rNode.pGetDof(rVarZ));

    const double n_norm = norm_2(NormalVector);
    KRATOS_ERROR_IF( n_norm < std::numeric_limits<double>::epsilon()) << "The norm of the normal vector is zero or almost zero" << std::endl;
    NormalVector /= n_norm;

    IndexType max_n_component_index = 0;
    double max_abs_n_component = std::abs(NormalVector[0]);
    for(IndexType i=1; i<dim; ++i) {
        if(std::abs(NormalVector[i]) > max_abs_n_component) {
            max_abs_n_component=std::abs(NormalVector[i]);
            max_n_component_index = i;
        }
    }

    mMasterDofsVector.clear();
    mSlaveDofsVector.clear();

    IndexType counter = 0;
    IndexType Tcounter = 0;
    for(auto& r_dof : all_dofs) {
        if(counter == max_n_component_index) {
            mSlaveDofsVector.push_back(r_dof);
        } else {
            mMasterDofsVector.push_back(r_dof);
            mRelationMatrix(0,Tcounter++) = -NormalVector[counter]/NormalVector[max_n_component_index];
        }

        counter++;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SlipConstraint::SetDofList(
    const DofPointerVectorType& rSlaveDofsVector,
    const DofPointerVectorType& rMasterDofsVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Slip constraint doesn't allow SetDofList" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SlipConstraint::SetSlaveDofsVector(const DofPointerVectorType& rSlaveDofsVector)
{
    KRATOS_ERROR << "Slip constraint doesn't allow SetSlaveDofsVector" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void SlipConstraint::SetMasterDofsVector(const DofPointerVectorType& rMasterDofsVector)
{
    KRATOS_ERROR << "Slip constraint doesn't allow SetMasterDofsVector" << std::endl;
}


/***********************************************************************************/
/***********************************************************************************/

void SlipConstraint::SetLocalSystem(
    const MatrixType& rRelationMatrix,
    const VectorType& rConstantVector,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR << "Slip constraint doesn't allow SetLocalSystem" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

std::string SlipConstraint::GetInfo() const
{
    return "SlipConstraint class !";
}

/***********************************************************************************/
/***********************************************************************************/

void SlipConstraint::PrintInfo(std::ostream &rOStream) const
{
    rOStream << " SlipConstraint Id  : " << this->Id() << std::endl;
    rOStream << " slave_dofs :" << std::endl;
    for(const auto& pDof : mSlaveDofsVector)
        rOStream << pDof->GetVariable().Name() << " of node : " << pDof->Id() << std::endl;
    rOStream << " master_dofs :" << std::endl;
    for(const auto& pDof : mMasterDofsVector)
        rOStream << pDof->GetVariable().Name() << " node : " << pDof->Id() << std::endl;
    rOStream << " relation matrix :" << std::endl;
    rOStream << mRelationMatrix << std::endl;
    // rOStream << " node Id " << mrNode.Id() << std::endl;
    // rOStream << " normal : " << mrNode.GetValue(*mpNormalVariable) << std::endl;
    // rOStream << " Number of Slaves          : " << mSlaveDofsVector.size() << std::endl;
    // rOStream << " Number of Masters         : " << mMasterDofsVector.size() << std::endl;
}

}  // namespace Kratos.

