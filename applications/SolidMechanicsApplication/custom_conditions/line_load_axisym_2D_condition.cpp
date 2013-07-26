//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/line_load_axisym_2D_condition.hpp"
#include "utilities/math_utils.h"
#include "solid_mechanics_application.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
LineLoadAxisym2DCondition::LineLoadAxisym2DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : LineLoad2DCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LineLoadAxisym2DCondition::LineLoadAxisym2DCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : LineLoad2DCondition( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
    mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
}


//************************************************************************************
//************************************************************************************
LineLoadAxisym2DCondition::LineLoadAxisym2DCondition(  LineLoadAxisym2DCondition const& rOther )
    : LineLoad2DCondition(rOther)
    , mThisIntegrationMethod(rOther.mThisIntegrationMethod)
{
}

//************************************************************************************
//************************************************************************************

Condition::Pointer LineLoadAxisym2DCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new LineLoadAxisym2DCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//************************************************************************************
//************************************************************************************

LineLoadAxisym2DCondition::~LineLoadAxisym2DCondition()
{
}


//************************************************************************************
//************************************************************************************
void LineLoadAxisym2DCondition::CalculateConditionalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = 2;

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    //calculating actual jacobian
    GeometryType::JacobiansType J;

    J = GetGeometry().Jacobian( J );

    //PRESSURE CONDITION:
    array_1d<double,3> PressureOnNodes (number_of_nodes,0.0);

    double PressureCondition = 0;
    //PressureCondition = GetValue( PRESSURE );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        PressureOnNodes[i] = PressureCondition;
        PressureOnNodes[i]+= GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE ) - GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );

    }

    //FORCE CONDITION:
    std::vector<array_1d<double,3> > ForceArray (number_of_nodes);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        ForceArray[i]=GetGeometry()[i].FastGetSolutionStepValue( FACE_LOAD );
    }


    //current gausspoint radius
    double CurrentRadius   = 0;
    double ReferenceRadius = 0;


    Vector NormalVector  = ZeroVector( 2 ); //normal direction (not normalized)
    Vector TangentVector = ZeroVector( 2 ); //normal direction (not normalized)

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //calculating tangent
        TangentVector[0] =  J[PointNumber]( 0, 0 ); // x_1,e
        TangentVector[1] =  J[PointNumber]( 1, 0 ); // x_2,e

        //Calculate IntegrationPoint radius
        Vector N=row(Ncontainer , PointNumber);
        CalculateRadius (CurrentRadius, ReferenceRadius, N);

        double IntegrationWeight = integration_points[PointNumber].Weight() * norm_2(TangentVector);

        if ( dimension == 2 ) IntegrationWeight *= 2 * 3.141592654 * CurrentRadius;

        //calculating normal
        NormalVector[0] = -J[PointNumber]( 1, 0 ); //-x_2,e
        NormalVector[1] =  J[PointNumber]( 0, 0 ); // x_1,e

        //calculating the pressure and force on the gauss point
        double gauss_pressure = 0.00;
        Vector ForceLoad = ZeroVector(dimension);

        for ( unsigned int ii = 0; ii < number_of_nodes; ii++ )
        {
            gauss_pressure += Ncontainer( PointNumber, ii ) * PressureOnNodes[ii];

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                ForceLoad[j] += Ncontainer( PointNumber, ii ) * ForceArray[ii][j];
            }
        }

        if ( CalculateStiffnessMatrixFlag == true )
        {
            if ( gauss_pressure != 0.00 )
            {
                CalculateAndSubKp( rLeftHandSideMatrix, DN_De[PointNumber], row( Ncontainer, PointNumber ), gauss_pressure, IntegrationWeight );
            }
        }

        //adding contributions to the residual vector
        if ( CalculateResidualVectorFlag == true )
        {
            if ( gauss_pressure != 0.00 )
                CalculateAndAddFacePressure( rRightHandSideVector, row( Ncontainer, PointNumber ), NormalVector, gauss_pressure, IntegrationWeight );

            if ( norm_2(ForceLoad) != 0.00)
                CalculateAndAddLineLoad( rRightHandSideVector, row( Ncontainer, PointNumber ), ForceLoad, IntegrationWeight );

        }

    }

    KRATOS_CATCH( "" )
}


//***********************************************************************
//***********************************************************************

void LineLoadAxisym2DCondition::CalculateRadius(double & rCurrentRadius,
        double & rReferenceRadius,
        const Vector& rN)


{

    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    rCurrentRadius=0;
    rReferenceRadius=0;

    if ( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            //Displacement from the reference to the current configuration
            array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);
            array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
            array_1d<double, 3 > & ReferencePosition    = GetGeometry()[i].Coordinates();
            array_1d<double, 3 > CurrentPosition        = ReferencePosition + DeltaDisplacement;

            rCurrentRadius   += CurrentPosition[0]*rN[i];
            rReferenceRadius += ReferencePosition[0]*rN[i];
            //std::cout<<" node "<<i<<" -> DeltaDisplacement : "<<DeltaDisplacement<<std::endl;
        }
    }


    if ( dimension == 3 )
    {
        std::cout<<" AXISYMMETRIC case and 3D is not possible "<<std::endl;
    }

    KRATOS_CATCH( "" )
}



} // Namespace Kratos


