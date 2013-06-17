/*
==============================================================================
KratosStructuralApplication
A library based on:
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
//
//   Project Name:        Kratos
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2007-12-12 14:51:07 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/elastic_constraint.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
ElasticConstraint::ElasticConstraint(IndexType NewId, GeometryType::Pointer
                           pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
ElasticConstraint::ElasticConstraint(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

Condition::Pointer ElasticConstraint::Create(IndexType NewId, NodesArrayType
                                        const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new ElasticConstraint(NewId,
                              GetGeometry().Create(ThisNodes), pProperties));
}

ElasticConstraint::~ElasticConstraint()
{
}


//************************************************************************************
//************************************************************************************
void ElasticConstraint::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ElasticConstraint::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                  CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
        
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************


void ElasticConstraint::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag,
                                      bool CalculateResidualVectorFlag )
{
    KRATOS_TRY
    
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = 3;

    //resizing LHS and RHS where needed
    if(rLeftHandSideMatrix.size1() != number_of_nodes*dim)
        rLeftHandSideMatrix.resize(number_of_nodes*dim,number_of_nodes*dim,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes*dim,number_of_nodes*dim);
    if(rRightHandSideVector.size() != number_of_nodes*dim)
        rRightHandSideVector.resize(number_of_nodes*dim,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_nodes*dim);
    

    if( number_of_nodes == 1 )
    {
        array_1d<double,3>& bedding = GetGeometry()[0].GetSolutionStepValue(ELASTIC_BEDDING_STIFFNESS);
        array_1d<double,3>& displacement = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT);
        for( unsigned int i=0; i<3; i++ )
        {
            rLeftHandSideMatrix(i,i) = bedding[i];
            rRightHandSideVector[i] = -displacement[i]*bedding[i];
        }
    }
    else
    {
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( GetGeometry().GetDefaultIntegrationMethod() );
        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( GetGeometry().GetDefaultIntegrationMethod() );
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( GetGeometry().GetDefaultIntegrationMethod() );
        
        array_1d<double, 3 > tangent_1;
        array_1d<double, 3 > tangent_2;
        array_1d<double, 3 > normal_vector;
             
        for (unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++)
        {
            Vector bedding = ZeroVector(3);
            Vector displacement = ZeroVector(3);
            
            for( unsigned int i=0; i<3; i++ )
            {
                tangent_1[i] = 0.0;
                tangent_2[i] = 0.0;
                normal_vector[i] = 0.0;
            }
             
            for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
            {
                tangent_1 += GetGeometry()[node] * DN_De[PointNumber]( node, 0 );
                tangent_2 += GetGeometry()[node] * DN_De[PointNumber]( node, 1 );
            
                noalias( displacement ) += GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT ) * Ncontainer(PointNumber,node);
                noalias( bedding ) += GetGeometry()[node].GetSolutionStepValue( ELASTIC_BEDDING_STIFFNESS ) * Ncontainer(PointNumber,node);
            }
        
            normal_vector = MathUtils<double>::CrossProduct( tangent_1, tangent_2 );
            double IntegrationWeight = integration_points[PointNumber].Weight();
            double dA = MathUtils<double>::Norm3( normal_vector );
        
            for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
            {
                for( unsigned int i=0; i<dim; i++ )
                {
                    if( CalculateResidualVectorFlag )
                    {
                        rRightHandSideVector[node*dim+i] -= displacement[i]*bedding[i]*IntegrationWeight*dA;
                    }
                    if( CalculateStiffnessMatrixFlag )
                    {
                        rLeftHandSideMatrix(node*dim+i,node*dim+i) += bedding[i]*IntegrationWeight*dA;
                    }
                }
            }
        }
    }
    
    KRATOS_CATCH("")
    
}


//************************************************************************************
//************************************************************************************
void ElasticConstraint::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int index;
    unsigned int dim = 3;
    rResult.resize(number_of_nodes*dim);
    for (int i=0; i<number_of_nodes; i++)
    {
        index = i*dim;
        rResult[index] = (GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId());
        rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());
        rResult[index+2] = (GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId());
    }
}

//************************************************************************************
//************************************************************************************
void ElasticConstraint::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
    unsigned int dim = 3;
    ConditionalDofList.resize(GetGeometry().size()*dim);
    unsigned int index;
    for (unsigned int i=0; i<GetGeometry().size(); i++)
    {
        index = i*dim;
        ConditionalDofList[index] = (GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ConditionalDofList[index+1] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        ConditionalDofList[index+2] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }
}


} // Namespace Kratos



