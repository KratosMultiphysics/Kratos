// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/line_load_condition_2d.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
LineLoadCondition2D::LineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LineLoadCondition2D::LineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

Condition::Pointer LineLoadCondition2D::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return boost::make_shared<LineLoadCondition2D>(NewId, pGeom, pProperties);
}

Condition::Pointer LineLoadCondition2D::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return boost::make_shared<LineLoadCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

LineLoadCondition2D::~LineLoadCondition2D()
{
}


//************************************************************************************
//************************************************************************************
void LineLoadCondition2D::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************
void LineLoadCondition2D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//************************************************************************************
//************************************************************************************
void LineLoadCondition2D::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo,
                                        bool CalculateStiffnessMatrixFlag,
                                        bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = NumberOfNodes * dim;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
        {
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != MatSize )
        {
            rRightHandSideVector.resize( MatSize, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    //calculating actual jacobian
    GeometryType::JacobiansType J;

    J = GetGeometry().Jacobian( J );

    ////sizing work matrices
    Vector PressureOnNodes = ZeroVector( NumberOfNodes );

    // Pressure applied to the element itself
    double PressureOnCondition = 0.0;
    if( this->Has( PRESSURE ) )
    {
        PressureOnCondition += this->GetValue( PRESSURE );
    }
    if( this->Has( NEGATIVE_FACE_PRESSURE ) )
    {
        PressureOnCondition += this->GetValue( NEGATIVE_FACE_PRESSURE );
    }
    if( this->Has( POSITIVE_FACE_PRESSURE ) )
    {
        PressureOnCondition -= this->GetValue( POSITIVE_FACE_PRESSURE );
    }
    
    for ( unsigned int i = 0; i < PressureOnNodes.size(); i++ )
    {
        PressureOnNodes[i] = PressureOnCondition;
        if( GetGeometry()[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) )
        {
            PressureOnNodes[i] += GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE );
        }
        if( GetGeometry()[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) )
        {
            PressureOnNodes[i] -= GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
        }
    }

    // Vector with a loading applied to the elemnt
    array_1d<double, 3 > LineLoad = ZeroVector(3);
    if( this->Has( LINE_LOAD ) )
    {
        noalias(LineLoad) = this->GetValue( LINE_LOAD );
    }
    
    Vector v3 = ZeroVector( 2 ); //normal direction (not normalized)
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        double IntegrationWeight = integration_points[PointNumber].Weight();

        v3[0] = -J[PointNumber]( 1, 0 );
        v3[1] = J[PointNumber]( 0, 0 );

        // Calculating the pressure on the gauss point
        double GaussPressure = 0.0;
        for ( unsigned int ii = 0; ii < NumberOfNodes; ii++ )
        {
            GaussPressure += Ncontainer( PointNumber, ii ) * PressureOnNodes[ii];
        }

        if ( CalculateStiffnessMatrixFlag == true )
        {
            if ( GaussPressure != 0.0 )
            {
                CalculateAndSubKp( rLeftHandSideMatrix, DN_De[PointNumber], row( Ncontainer, PointNumber ), GaussPressure, IntegrationWeight );
            }
        }

        //adding contributions to the residual vector
        if ( CalculateResidualVectorFlag == true )
        {
            if ( GaussPressure != 0.0 )
            {
                CalculateAndAdd_PressureForce( rRightHandSideVector, row( Ncontainer, PointNumber ), v3, GaussPressure, IntegrationWeight );
            }
        }
        
        array_1d<double,3> GaussLoad = LineLoad;
        for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
        {
            if( GetGeometry()[ii].SolutionStepsDataHas( LINE_LOAD ) )
            {
                noalias(GaussLoad) += Ncontainer( PointNumber, ii )*GetGeometry()[ii].FastGetSolutionStepValue( LINE_LOAD );
            }
        }

        for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
        {
            unsigned int base = ii*3;
            for(unsigned int k=0; k<3; ++k)
            {
                rRightHandSideVector[base+k] += Ncontainer( PointNumber, ii )*GaussLoad[k];
            }
        }

    }

    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************

void LineLoadCondition2D::CalculateAndSubKp(
    Matrix& K,
    const Matrix& DN_De,
    const Vector& N,
    double pressure,
    double weight
)
{
    KRATOS_TRY

    Matrix Kij( 2, 2 );
    Matrix Cross_gn( 2, 2 );

    //double h0 = GetProperties()[THICKNESS];
    double h0 = 0.00;
    Cross_gn( 0, 0 ) = 0.0;
    Cross_gn( 0, 1 ) = -h0;
    Cross_gn( 1, 0 ) = -h0;
    Cross_gn( 1, 1 ) = 0.0;

    double coeff;

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        unsigned int RowIndex = i * 2;

        for ( unsigned int j = 0; j < GetGeometry().size(); j++ )
        {
            unsigned int ColIndex = j * 2;

            coeff = pressure * N[i] * DN_De( j, 0 ) * weight;
            Kij = -coeff * Cross_gn;

            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
            MathUtils<double>::SubtractMatrix( K, Kij, RowIndex, ColIndex );
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************
void LineLoadCondition2D::CalculateAndAdd_PressureForce(
    Vector& residualvector,
    const Vector& N,
    Vector& v3,
    double pressure,
    double weight )
{
    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int dim = 2;

    for ( unsigned int i = 0; i < NumberOfNodes; i++ )
    {
        int index = dim * i;
        double coeff = pressure * N[i] * weight;
        residualvector[index]   += coeff * v3[0];
        residualvector[index+1] += coeff * v3[1];
    }

}

} // Namespace Kratos


