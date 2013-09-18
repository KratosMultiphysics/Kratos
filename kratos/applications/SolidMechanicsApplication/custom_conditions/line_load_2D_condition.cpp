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
#include "custom_conditions/line_load_2D_condition.hpp"
#include "utilities/math_utils.h"
#include "solid_mechanics_application.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
LineLoad2DCondition::LineLoad2DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LineLoad2DCondition::LineLoad2DCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
LineLoad2DCondition::LineLoad2DCondition(  LineLoad2DCondition const& rOther )
    : Condition(rOther)
{

}

//************************************************************************************
//************************************************************************************

Condition::Pointer LineLoad2DCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new LineLoad2DCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//************************************************************************************
//************************************************************************************

LineLoad2DCondition::~LineLoad2DCondition()
{
}

//************************************************************************************
//************************************************************************************
void LineLoad2DCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateConditionalSystem( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    //KRATOS_WATCH(rRightHandSideVector);
}

//************************************************************************************
//************************************************************************************
void LineLoad2DCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateConditionalSystem( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

}

//************************************************************************************
//************************************************************************************
void LineLoad2DCondition::CalculateConditionalSystem( MatrixType& rLeftHandSideMatrix,
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


    Vector NormalVector  = ZeroVector( 2 ); //normal direction (not normalized)
    Vector TangentVector = ZeroVector( 2 ); //normal direction (not normalized)

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //calculating tangent
        TangentVector[0] =  J[PointNumber]( 0, 0 ); // x_1,e
        TangentVector[1] =  J[PointNumber]( 1, 0 ); // x_2,e

        double IntegrationWeight = integration_points[PointNumber].Weight() * norm_2(TangentVector);
	//double IntegrationWeight = integration_points[PointNumber].Weight() * 0.5 * GetGeometry().Length();

        if ( dimension == 2 ) IntegrationWeight *= GetProperties()[THICKNESS];

        //calculating normal
        NormalVector[0] = -J[PointNumber]( 1, 0 ); //-x_2,e
        NormalVector[1] =  J[PointNumber]( 0, 0 ); // x_1,e
	NormalVector/= norm_2(NormalVector);

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
void LineLoad2DCondition::CalculateAndSubKp( Matrix& rK,
        const Matrix& rDN_De,
        const Vector& rN,
        double rPressure,
        double rIntegrationWeight
                                           )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    Matrix Kij     ( 2, 2 );
    Matrix SkewSymmMatrix( 2, 2 );

    //Compute the K sub matrix
    SkewSymmMatrix( 0, 0 ) =  0.0;
    SkewSymmMatrix( 0, 1 ) = -1.0;
    SkewSymmMatrix( 1, 0 ) = -1.0;
    SkewSymmMatrix( 1, 1 ) =  0.0;

    double DiscretePressure=0;
    unsigned int RowIndex = 0;
    unsigned int ColIndex = 0;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        RowIndex = i * 2;

        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            ColIndex = j * 2;

            DiscretePressure = rPressure * rN[i] * rDN_De( j, 0 ) * rIntegrationWeight;
            Kij = - DiscretePressure * SkewSymmMatrix;

            //TAKE CARE: the load correction matrix should be SUBSTRACTED not added
            MathUtils<double>::SubtractMatrix( rK, Kij, RowIndex, ColIndex );
        }
    }


    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************
void LineLoad2DCondition::CalculateAndAddFacePressure(Vector& rF,
        const Vector& rN,
        Vector& rNormal,
        double rPressure,
        double rIntegrationWeight )
{

    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension  = 2;
    unsigned int index = 0;
    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);

        index = dimension * i;
        double DiscretePressure = rPressure * rN[i] * rIntegrationWeight;
        rF[index]   += DiscretePressure * rNormal[0];
        rF[index+1] += DiscretePressure * rNormal[1];

        ExternalForce[0] +=DiscretePressure * rNormal[0];
        ExternalForce[1] +=DiscretePressure * rNormal[1];

    }

    KRATOS_CATCH("")
}



//***********************************************************************
//***********************************************************************

void LineLoad2DCondition::CalculateAndAddLineLoad(Vector& rF,
        const Vector& rN,
        Vector& rForce,
        double  rIntegrationWeight )
{

    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = 2;
    unsigned int index = 0;


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        index = dimension * i;

        array_1d<double, 3 >& ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(FORCE_EXTERNAL);


        for ( unsigned int idim = 0; idim < dimension; idim++ )
        {
            rF[index+idim] += rN[i] * rForce[idim] * rIntegrationWeight;

            ExternalForce[idim] += rN[i] * rForce[idim] * rIntegrationWeight;
        }

	//std::cout<<" External Force stored in Node: "<<GetGeometry()[i].Id()<<" Force Load "<<rForce<<" norm "<<norm_2(rForce)<<std::endl;
	
    }

    

    KRATOS_CATCH("")

}


//************************************************************************************
//************************************************************************************
void LineLoad2DCondition::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    int number_of_nodes = GetGeometry().PointsNumber(); //=Geometry().size();
    const unsigned int dimension = 2;

    unsigned int conditiondimension = number_of_nodes*dimension;

    if ( rResult.size() != conditiondimension )
        rResult.resize( conditiondimension, false );

    unsigned int index = 0;
    for ( int i = 0; i < number_of_nodes; i++ )
    {
        index = i * dimension;
        rResult[index]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
void LineLoad2DCondition::GetDofList( DofsVectorType& rConditionalDofList, ProcessInfo& rCurrentProcessInfo )
{
    rConditionalDofList.resize( 0 );
    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        rConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        rConditionalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
    }

}

//************************************************************************************
//************************************************************************************
void  LineLoad2DCondition::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if ( rMassMatrix.size1() != 0 )
        rMassMatrix.resize( 0, 0, false );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void  LineLoad2DCondition::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if ( rDampMatrix.size1() != 0 )
        rDampMatrix.resize( 0, 0, false );

    KRATOS_CATCH( "" )
}

int LineLoad2DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

} // Namespace Kratos


