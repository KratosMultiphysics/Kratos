/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

/* ****************************************************************************
*  Projectname:         $KratosALEApplication
*  Last Modified by:    $Author: Andreas.Mini@tum.de $
*  Date:                $Date: November 2015 $
*  Revision:            $Revision: 1.4 $
* *****************************************************************************/


// System includes
#include <cmath>

// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/structural_meshmoving_element.h"
#include "ale_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//======================================================================
//Get displacement values
//======================================================================
template<>
void StructuralMeshMovingElement<2>::GetDisplacementValues(VectorType& rValues,
                                                           const int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;
    
    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);
    
    SizeType Index = 0;
    
    for (SizeType iNode = 0; iNode < NumNodes; ++iNode){

        rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_X,Step);
        rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_Y,Step);
    }
}

template<>
void StructuralMeshMovingElement<3>::GetDisplacementValues(VectorType& rValues,
                                                           const int Step)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*3;
    
    if (rValues.size() != LocalSize)
        rValues.resize(LocalSize,false);
    
    SizeType Index = 0;
    
    for (SizeType iNode = 0; iNode < NumNodes; ++iNode){

        rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_X,Step);
        rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_Y,Step);
        rValues[Index++] = rGeom[iNode].FastGetSolutionStepValue(DISPLACEMENT_Z,Step);
    }
}

//======================================================================
//Build up System Matrices
//======================================================================
template<>
void StructuralMeshMovingElement<2>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                          VectorType& rRightHandSideVector,
                                                          ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const SizeType LocalSize = NumNodes * 2;

    // Check sizes and initialize
    if(rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
    
    if(rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize,false);

    // Get the element's geometric parameters
    double Area;
    array_1d<double,3> N;
    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
    VectorType DetJ;
    this->GetGeometry().DeterminantOfJacobian(DetJ);

    // Plane strain ConstitutiveMatrix
    MatrixType ConstitutiveMatrix = ZeroMatrix(3,3);
    double YoungsModulus = 200000;
    double PoissonCoefficient = 0.3;

    // Stiffening of elements using Jacobian determinants and exponent between 0.0 and 2.0
    double J0 = 100;          // Factor influences how far the displacement spreads into the fluid mesh
    double Xi = 1.5;          // 1.5 Exponent influences stiffening of smaller elements; 0 = no stiffening
    double DetJMag = std::abs(DetJ[0]);
    double Quotient = J0 / DetJMag;
    YoungsModulus *= (DetJMag *  pow(Quotient, Xi));

    double Lambda = (YoungsModulus * PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient));
    double Mue = YoungsModulus/(2*(1-PoissonCoefficient)); // modified for artificial stiffening

    ConstitutiveMatrix ( 0 , 0 ) = Lambda + 2*Mue;
    ConstitutiveMatrix ( 1 , 1 ) = ConstitutiveMatrix ( 0 , 0 );
    ConstitutiveMatrix ( 2 , 2 ) = Mue;
    ConstitutiveMatrix ( 0 , 1 ) = Lambda;
    ConstitutiveMatrix ( 1 , 0 ) = Lambda;

    // B-matrix
    MatrixType B = ZeroMatrix(3,6);
    SizeType Index = 0;
    for(SizeType iNode=0; iNode<NumNodes; iNode++){

        B( 0, Index + 0 ) = DN_DX( iNode, 0 );
        B( 1, Index + 1 ) = DN_DX( iNode, 1 );
        B( 2, Index + 0 ) = DN_DX( iNode, 1 );
        B( 2, Index + 1 ) = DN_DX( iNode, 0 );
        Index += 2;
    }

    // Compute LHS
    MatrixType intermediateMatrix = prod(trans(B),ConstitutiveMatrix);
    noalias(rLeftHandSideMatrix) = prod(intermediateMatrix,B);

    // Compute RHS
    VectorType LastValues = ZeroVector(LocalSize);
    this->GetDisplacementValues(LastValues,0);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,LastValues);
    
    KRATOS_CATCH("");
}


template<>
void StructuralMeshMovingElement<3>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                          VectorType& rRightHandSideVector,
                                                          ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const SizeType LocalSize = NumNodes * 3;

    // Check sizes and initialize
    if(rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
    
    if(rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize,false);

    // Get the element's geometric parameters
    double Area;
    array_1d<double,4> N;
    boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
    Vector DetJ;
    this->GetGeometry().DeterminantOfJacobian(DetJ);

    // structural ConstitutiveMatrix
    MatrixType ConstitutiveMatrix = ZeroMatrix(6,6);
    double YoungsModulus = 200000;
    double PoissonCoefficient = 0.3;

    //Stiffening of elements using Jacobian determinants and exponent between 0.0 and 2.0
    double J0 = 100;          // Factor influences how far the displacement spreads into the fluid mesh
    double Xi = 1.5;          // 1.5 Exponent influences stiffening of smaller elements; 0 = no stiffening
    double DetJMag = std::abs(DetJ[0]);
    double Quotient = J0 / DetJMag;
    YoungsModulus *= (DetJMag *  pow(Quotient, Xi));

    double C00 = (YoungsModulus*(1-PoissonCoefficient)/((1+PoissonCoefficient)*(1-2*PoissonCoefficient)));
    double C33 = C00*(1-2*PoissonCoefficient)/(2*(1-PoissonCoefficient));
    double C01 = C00*PoissonCoefficient/(1-PoissonCoefficient);

    ConstitutiveMatrix ( 0 , 0 ) = C00;
    ConstitutiveMatrix ( 1 , 1 ) = C00;
    ConstitutiveMatrix ( 2 , 2 ) = C00;
    ConstitutiveMatrix ( 3 , 3 ) = C33;
    ConstitutiveMatrix ( 4 , 4 ) = C33;
    ConstitutiveMatrix ( 5 , 5 ) = C33;

    ConstitutiveMatrix ( 0 , 1 ) = C01;
    ConstitutiveMatrix ( 1 , 0 ) = C01;
    ConstitutiveMatrix ( 0 , 2 ) = C01;
    ConstitutiveMatrix ( 2 , 0 ) = C01;
    ConstitutiveMatrix ( 1 , 2 ) = C01;
    ConstitutiveMatrix ( 2 , 1 ) = C01;

    // B-matrix
    MatrixType B = ZeroMatrix(6,12);
    SizeType Index = 0;
    for(SizeType iNode=0; iNode<NumNodes; iNode++){

        B ( 0 , Index + 0 ) = DN_DX( iNode, 0 );
        B ( 1 , Index + 1 ) = DN_DX( iNode, 1 );
        B ( 2 , Index + 2 ) = DN_DX( iNode, 2 );
        B ( 3 , Index + 0 ) = DN_DX( iNode, 1 );
        B ( 3 , Index + 1 ) = DN_DX( iNode, 0 );
        B ( 4 , Index + 1 ) = DN_DX( iNode, 2 );
        B ( 4 , Index + 2 ) = DN_DX( iNode, 1 );
        B ( 5 , Index + 0 ) = DN_DX( iNode, 2 );
        B ( 5 , Index + 2 ) = DN_DX( iNode, 0 );

        Index += 3;
    }

    // Compute LHS
    MatrixType intermediateMatrix = prod(trans(B),ConstitutiveMatrix);
    noalias(rLeftHandSideMatrix) = prod(intermediateMatrix,B);

    //double prefactor = lambda + (2/dimension)*mue;
    //rLeftHandSideMatrix *= prefactor;

    // Compute RHS
    VectorType LastValues = ZeroVector(LocalSize);
    this->GetDisplacementValues(LastValues,0);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,LastValues);

    KRATOS_CATCH("");
}

//======================================================================
//Generate Equation Id Vector
//======================================================================
template<>
void StructuralMeshMovingElement<2>::EquationIdVector(EquationIdVectorType& rResult,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;
    
    if(rResult.size() != LocalSize)
        rResult.resize(LocalSize,false);
    
    SizeType Index = 0;
    
    for (SizeType iNode = 0; iNode < NumNodes; iNode++){

        rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_X).EquationId();
        rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_Y).EquationId();
    }
}

template<>
void StructuralMeshMovingElement<3>::EquationIdVector(EquationIdVectorType& rResult,
                                                      ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*3;
    
    if(rResult.size() != LocalSize)
        rResult.resize(LocalSize,false);
    
    SizeType Index = 0;
    
    for (SizeType iNode = 0; iNode < NumNodes; iNode++){

        rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_X).EquationId();
        rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[Index++] = rGeom[iNode].GetDof(DISPLACEMENT_Z).EquationId();
    }
}

//======================================================================
//Get dof List
//======================================================================
template<>
void StructuralMeshMovingElement<2>::GetDofList(DofsVectorType& rElementalDofList,
                                                ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;
    
    if(rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;
    
    for (SizeType iNode = 0; iNode < NumNodes; iNode++){

        rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_X);
        rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_Y);
    }
}


template<>
void StructuralMeshMovingElement<3>::GetDofList(DofsVectorType& rElementalDofList,
                                                ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*3;
    
    if(rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;
    
    for (SizeType iNode = 0; iNode < NumNodes; iNode++){

        rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_X);
        rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[LocalIndex++] = rGeom[iNode].pGetDof(DISPLACEMENT_Z);
    }
}

} // Namespace Kratos
