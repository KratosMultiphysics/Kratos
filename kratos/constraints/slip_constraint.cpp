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
    Dof<double>* pDofX,
    Dof<double>* pDofY,
    array_1d<double,3> NormalVector //note that a copy is made
    ) : BaseType(Id)
{

    DofPointerVectorType all_dofs;
    all_dofs.reserve(2);
    all_dofs.push_back(pDofX);
    all_dofs.push_back(pDofY);

    ConstructorHelper(all_dofs, NormalVector);
}

SlipConstraint::SlipConstraint(
    IndexType Id,
    Dof<double>* pDofX,
    Dof<double>* pDofY,
    Dof<double>* pDofZ,
    array_1d<double,3> NormalVector //note that a copy is made
    ) : BaseType(Id)
{

    DofPointerVectorType all_dofs;
    all_dofs.reserve(3);
    all_dofs.push_back(pDofX);
    all_dofs.push_back(pDofY);
    all_dofs.push_back(pDofZ);

    ConstructorHelper(all_dofs, NormalVector);
}

void SlipConstraint::ConstructorHelper(
    DofPointerVectorType& rAllDofs,
    array_1d<double,3>& rNormalVector
    ) 
{
    const unsigned int dim = rAllDofs.size();
    mRelationMatrix.resize(1,dim-1,false);
    mRelationMatrix.clear();

    mConstantVector.resize(1,false);
    mConstantVector.clear();

    const double n_norm = norm_2(rNormalVector);
    KRATOS_ERROR_IF( n_norm < std::numeric_limits<double>::epsilon()) << "The norm of the normal vector is zero or almost zero" << std::endl;
    rNormalVector /= n_norm;

    IndexType max_n_component_index = 0;
    double max_abs_n_component = std::abs(rNormalVector[0]);
    for(IndexType i=1; i<dim; ++i) {
        if(std::abs(rNormalVector[i]) > max_abs_n_component) {
            max_abs_n_component=std::abs(rNormalVector[i]);
            max_n_component_index = i;
        }
    }

    mMasterDofsVector.clear();
    mSlaveDofsVector.clear();

    IndexType counter = 0;
    IndexType Tcounter = 0;
    for(auto& r_dof : rAllDofs) {
        if(counter == max_n_component_index) {
            mSlaveDofsVector.push_back(r_dof);
        } else {
            mMasterDofsVector.push_back(r_dof);
            mRelationMatrix(0,Tcounter++) = -rNormalVector[counter]/rNormalVector[max_n_component_index];
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

