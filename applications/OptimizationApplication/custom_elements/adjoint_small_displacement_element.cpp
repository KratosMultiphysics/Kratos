//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//

// System includes


// External includes


// Project includes
#include "custom_elements/adjoint_small_displacement_element.h"
#include "optimization_application_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
AdjointSmallDisplacementElement::AdjointSmallDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry, Element::Pointer pPrimal)
    : Element(NewId, pGeometry)
{
    mpPrimalElement = pPrimal;
}

//************************************************************************************
//************************************************************************************
// AdjointSmallDisplacementElement::AdjointSmallDisplacementElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
//     : Element(NewId, pGeometry, pProperties)
// {
//     const Element& rElem = KratosComponents<Element>::Get("SmallDisplacementElement3D4N");
//     mpPrimalElement = rElem.Create(NewId, pGeometry, pProperties);
// }

// Element::Pointer AdjointSmallDisplacementElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
// {
//     return Kratos::make_intrusive<AdjointSmallDisplacementElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
// }

// Element::Pointer AdjointSmallDisplacementElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
// {
//     return Kratos::make_intrusive<AdjointSmallDisplacementElement>(NewId, pGeom, pProperties);
// }

AdjointSmallDisplacementElement::~AdjointSmallDisplacementElement()
{
}

//************************************************************************************
//************************************************************************************
void AdjointSmallDisplacementElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                            rRightHandSideVector,
                                            rCurrentProcessInfo);


    // Some working variables
    const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
    const SizeType dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
    const SizeType num_dofs_per_node = dimension;
    const SizeType num_dofs = num_nodes * num_dofs_per_node;

    // Resizing as needed the RHS 
    if ( rRightHandSideVector.size() != num_dofs )
        rRightHandSideVector.resize( num_dofs, false );
    rRightHandSideVector = ZeroVector( num_dofs ); //resetting RHS


    Vector nodal_vals(num_dofs);
    for(unsigned int node_element = 0; node_element<num_nodes; node_element++)
    {
        const VectorType &source = GetGeometry()[node_element].FastGetSolutionStepValue(ADJOINT_RHS);
        auto node_weight = GetGeometry()[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        nodal_vals[3 * node_element + 0] = source[0]/node_weight;
        nodal_vals[3 * node_element + 1] = source[1]/node_weight;
        nodal_vals[3 * node_element + 2] = source[2]/node_weight;
    }

    noalias(rRightHandSideVector) += nodal_vals;

    //apply drichlet BC
    Vector temp;
    GetValuesVector(temp,0);    
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp); 

    KRATOS_CATCH("")
}

//******************************************************************************
//******************************************************************************
void AdjointSmallDisplacementElement::GetValuesVector(VectorType &rValues,
                                            int Step) const {
  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  if (dimension == 2) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"), Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"), Step);
    }
  } else if (dimension == 3) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"), Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"), Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z"), Step);
    }
  }
}

/***********************************************************************************/
/***********************************************************************************/

//************************************************************************************
//************************************************************************************
void AdjointSmallDisplacementElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void AdjointSmallDisplacementElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void AdjointSmallDisplacementElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"));

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 2;
            rResult[index] = GetGeometry()[i].GetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"),pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"),pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"),pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"),pos+1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z"),pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void AdjointSmallDisplacementElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{

    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X")));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y")));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X")));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y")));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z")));
        }
    }

    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
int AdjointSmallDisplacementElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = r_geometry[i];
        // KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VARS_SHAPE,rnode)

        KRATOS_CHECK_DOF_IN_NODE(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_X"), rnode)
        KRATOS_CHECK_DOF_IN_NODE(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Y"), rnode)
        KRATOS_CHECK_DOF_IN_NODE(KratosComponents<Variable<double>>::Get("ADJOINT_DISPLACEMENT_Z"), rnode)
    }

    return check;

    KRATOS_CATCH( "" );
}
/***********************************************************************************/
/***********************************************************************************/


} // Namespace Kratos
