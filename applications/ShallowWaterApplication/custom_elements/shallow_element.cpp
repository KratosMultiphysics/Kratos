//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "custom_elements/shallow_element.h"
#include "shallow_water_application_variables.h"

namespace Kratos
{

//----------------------------------------------------------------------

int ShallowElement::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    constexpr unsigned int nnodes = 3;

    // Check that all required variables have been registered
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(HEIGHT)
    KRATOS_CHECK_VARIABLE_KEY(PROJECTED_SCALAR1)
    KRATOS_CHECK_VARIABLE_KEY(PROJECTED_VECTOR1)
    KRATOS_CHECK_VARIABLE_KEY(BATHYMETRY)
    KRATOS_CHECK_VARIABLE_KEY(RAIN)
    KRATOS_CHECK_VARIABLE_KEY(MANNING)
    KRATOS_CHECK_VARIABLE_KEY(GRAVITY)
    KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME)
    KRATOS_CHECK_VARIABLE_KEY(WATER_HEIGHT_UNIT_CONVERTER)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( unsigned int i = 0; i < nnodes; i++ )
    {
        Node<3>& node = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HEIGHT, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_VECTOR1, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PROJECTED_SCALAR1, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(BATHYMETRY, node)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RAIN, node)

        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, node)
        KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, node)
        KRATOS_CHECK_DOF_IN_NODE(HEIGHT, node)
    }

    return ierr;

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int nnodes = 3;
    constexpr unsigned int elem_size = nnodes*3;

    if(rResult.size() != elem_size)
        rResult.resize(elem_size,false); // False says not to preserve existing storage!!

    GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (unsigned int i = 0; i < nnodes; i++)
    {
        rResult[counter++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
        rResult[counter++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
        rResult[counter++] = rGeom[i].GetDof(HEIGHT).EquationId();
    }

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowElement::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int nnodes = 3;
    constexpr unsigned int elem_size = nnodes*3;

    if(rElementalDofList.size() != elem_size)
        rElementalDofList.resize(elem_size);
    
    GeometryType& rGeom = GetGeometry();
    int counter=0;
    for (unsigned int i = 0; i < nnodes; i++)
    {
        rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
    }
    
    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int nnodes = 3;
    constexpr unsigned int elem_size = nnodes*3;

    // Resize of the Left and Right Hand side
    if(rLeftHandSideMatrix.size1() != elem_size)
        rLeftHandSideMatrix.resize(elem_size,elem_size,false); // False says not to preserve existing storage!!

    if(rRightHandSideVector.size() != elem_size)
        rRightHandSideVector.resize(elem_size,false);          // False says not to preserve existing storage!!

    KRATOS_CATCH("")
}

//----------------------------------------------------------------------

void ShallowElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
}

//----------------------------------------------------------------------

void ShallowElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
}

//----------------------------------------------------------------------

void ShallowElement::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == VEL_ART_VISC || rVariable == PR_ART_VISC || rVariable == RESIDUAL_NORM || rVariable == MIU)
    {
        for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++) 
            rValues[PointNumber] = double(this->GetValue(rVariable));
    }
}


} // namespace kratos