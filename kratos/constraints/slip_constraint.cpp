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
    const Variable<array_1d<double,3> >& rNormalVar
    ) : BaseType(Id),
        mrNode(rNode),
        mpNormalVariable(&rNormalVar)
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
    KRATOS_ERROR_IF(norm_2(n) < std::numeric_limits<double>::epsilon()) << "The norm of the normal vector is zero or almost zero" << std::endl;
    n /= norm_2(n);

    array_1d<double,3> v;
    v[0] = all_dofs[0]->GetSolutionStepValue(rVarX);
    v[1] = all_dofs[1]->GetSolutionStepValue(rVarY);
    v[2] = all_dofs[2]->GetSolutionStepValue(rVarZ);

    IndexType max_n_component_index = 0;
    double max_abs_n_component = std::abs(n[0]);
    for(IndexType i=1; i<3; ++i) {
        if(std::abs(n[i]) > max_abs_n_component) {
            max_abs_n_component=std::abs(n[i]);
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
            mRelationMatrix(0,Tcounter++) = -n[counter]/n[max_n_component_index];
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
    rOStream << " node Id " << mrNode.Id() << std::endl;
    rOStream << " normal : " << mrNode.FastGetSolutionStepValue(*mpNormalVariable) << std::endl;
    rOStream << " Number of Slaves          : " << mSlaveDofsVector.size() << std::endl;
    rOStream << " Number of Masters         : " << mMasterDofsVector.size() << std::endl;
}

}  // namespace Kratos.

