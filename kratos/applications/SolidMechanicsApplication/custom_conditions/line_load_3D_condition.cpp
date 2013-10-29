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
#include "custom_conditions/line_load_3D_condition.hpp"
#include "custom_utilities/solid_mechanics_math_utilities.hpp"

#include "solid_mechanics_application.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
LineLoad3DCondition::LineLoad3DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
}

//************************************************************************************
//************************************************************************************
LineLoad3DCondition::LineLoad3DCondition( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
}

//************************************************************************************
//************************************************************************************
LineLoad3DCondition::LineLoad3DCondition(  LineLoad3DCondition const& rOther )
    : Condition(rOther)
{

}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer LineLoad3DCondition::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new LineLoad3DCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}


//***********************************************************************************
//***********************************************************************************
// Destructor
LineLoad3DCondition::~LineLoad3DCondition()
{
}

//***********************************************************************************
//***********************************************************************************
void LineLoad3DCondition::EquationIdVector( EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int conditionaldimension = number_of_nodes * 3;

    if ( rResult.size() != conditionaldimension )
        rResult.resize( conditionaldimension );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * 3;
        rResult[index]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LineLoad3DCondition::GetDofList( DofsVectorType& ElementalDofList,
                                      ProcessInfo& rCurrentProcessInfo )
{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//***********************************************************************************
//***********************************************************************************
void LineLoad3DCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = 3;
    unsigned int MatSize = number_of_nodes * dimension;

    //resizing as needed the RHS
    if ( rRightHandSideVector.size() != MatSize )
        rRightHandSideVector.resize( MatSize, false );

    rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    unsigned int index = 0;
    //loop over integration points
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        //FORCE CONDITION:
        array_1d <double, 3> FaceLoad = ZeroVector( dimension );

        for ( unsigned int n = 0; n < GetGeometry().size(); n++ )
        {
            array_1d<double,3>&  LoadArray = ( GetGeometry()[n] ).GetSolutionStepValue( FACE_LOAD );

            for ( unsigned int idim = 0; idim < dimension; idim++ )
            {
                FaceLoad( idim ) += LoadArray( idim ) * Ncontainer( PointNumber, n );
            }
        }

        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

        IntegrationWeight *= GetGeometry().Length() * 0.5;


        // RIGHT HAND SIDE VECTOR
        for ( unsigned int m = 0; m < GetGeometry().size(); m++ )
        {
            array_1d<double, 3 > & ExternalForce = GetGeometry()[m].FastGetSolutionStepValue(FORCE_EXTERNAL);

            index = dimension * m;
	    GetGeometry()[m].SetLock();
            for ( unsigned int idim = 0; idim < dimension; idim++ )
            {
                rRightHandSideVector( index + idim ) += Ncontainer( PointNumber, m ) * FaceLoad( idim ) * IntegrationWeight;

                ExternalForce[idim] += Ncontainer( PointNumber, m ) * FaceLoad( idim ) * IntegrationWeight;
            }
	    GetGeometry()[m].UnSetLock();
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LineLoad3DCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo )
{

    rLeftHandSideMatrix = ZeroMatrix( 6, 6 );
    CalculateRightHandSide( rRightHandSideVector, rCurrentProcessInfo);
}

//***********************************************************************************
//***********************************************************************************
void LineLoad3DCondition::MassMatrix( MatrixType& rMassMatrix,
                                      ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rMassMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
void LineLoad3DCondition::DampMatrix( MatrixType& rDampMatrix,
                                      ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rDampMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************
//void LineLoad3DCondition::GetValuesVector( Vector& values, int Step )
//{
//    unsigned int number_of_nodes = GetGeometry().size();
//    unsigned int MatSize = number_of_nodes * 3;

//    if ( values.size() != MatSize )
//        values.resize( MatSize );

//    for ( unsigned int i = 0; i < number_of_nodes; i++ )
//    {
//        const array_1d<double, 3>& disp =
//            GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT, Step );
//        unsigned int index = i * 3;
//        values[index]   = disp[0];
//        values[index+1] = disp[1];
//        values[index+2] = disp[2];
//    }
//}

//***********************************************************************************
//***********************************************************************************
//void LineLoad3DCondition::GetFirstDerivativesVector( Vector& values, int Step )
//{
//    unsigned int number_of_nodes = GetGeometry().size();
//    unsigned int MatSize = number_of_nodes * 3;

//    if ( values.size() != MatSize )
//        values.resize( MatSize );

//    for ( unsigned int i = 0; i < number_of_nodes; i++ )
//    {
//        const array_1d<double, 3>& vel = GetGeometry()[i].FastGetSolutionStepValue( VELOCITY, Step );
//        unsigned int index = i * 3;
//        values[index]   = vel[0];
//        values[index+1] = vel[1];
//        values[index+2] = vel[2];
//    }
//}

//***********************************************************************************
//***********************************************************************************
//void LineLoad3DCondition::GetSecondDerivativesVector( Vector& values, int Step )
//{
//    unsigned int number_of_nodes = GetGeometry().size();
//    unsigned int MatSize = number_of_nodes * 3;

//    if ( values.size() != MatSize )
//        values.resize( MatSize );

//    for ( unsigned int i = 0; i < number_of_nodes; i++ )
//    {
//        const array_1d<double, 3>& acc =
//            GetGeometry()[i].FastGetSolutionStepValue( ACCELERATION, Step );
//        unsigned int index = i * 3;
//        values[index]   = acc[0];
//        values[index+1] = acc[1];
//        values[index+2] = acc[2];
//    }
//}

//***********************************************************************************
//***********************************************************************************
// --------- //
//  PRIVATE  //
// --------- //

//***********************************************************************************
//***********************************************************************************
//void LineLoad3DCondition::CalculateAndSubKp( Matrix& K, array_1d<double, 3>& ge,
//                                     array_1d<double, 3>& gn, const Matrix& DN_De,
//                                     const Vector& N, double pressure,
//                                     double weight )
//{
//    KRATOS_TRY
//    boost::numeric::ublas::bounded_matrix<double, 3, 3> Kij;
//    boost::numeric::ublas::bounded_matrix<double, 3, 3> Cross_ge;
//    boost::numeric::ublas::bounded_matrix<double, 3, 3> Cross_gn;
//    double coeff;
//    unsigned int number_of_nodes = GetGeometry().size();

//    MakeCrossMatrix( Cross_ge, ge );
//    MakeCrossMatrix( Cross_gn, gn );

//    for ( unsigned int i = 0; i < number_of_nodes; i++ )
//    {
//        unsigned int RowIndex = i * 3;

//        for ( unsigned int j = 0; j < number_of_nodes; j++ )
//        {
//            unsigned int ColIndex = j * 3;

//            coeff = pressure * N[i] * DN_De( j, 1 ) * weight;
//            noalias( Kij )  = coeff * Cross_ge;

//            coeff = pressure * N[i] * DN_De( j, 0 ) * weight;

//            noalias( Kij ) -= coeff * Cross_gn;
//            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
//            SubtractMatrix( K, Kij, RowIndex, ColIndex );
//        }
//    }

//    KRATOS_CATCH( "" )
//}

//***********************************************************************************
//***********************************************************************************
//void LineLoad3DCondition::MakeCrossMatrix( boost::numeric::ublas::bounded_matrix<double, 3, 3>& M,
//                                   array_1d<double, 3>& U )
//{
//    M( 0, 0 ) =  0.00;
//    M( 0, 1 ) = -U[2];
//    M( 0, 2 ) =  U[1];
//    M( 1, 0 ) =  U[2];
//    M( 1, 1 ) =  0.00;
//    M( 1, 2 ) = -U[0];
//    M( 2, 0 ) = -U[1];
//    M( 2, 1 ) =  U[0];
//    M( 2, 2 ) =  0.00;
//}

//***********************************************************************************
//***********************************************************************************
//void LineLoad3DCondition::CrossProduct( array_1d<double, 3>& cross,
//                                array_1d<double, 3>& a,
//                                array_1d<double, 3>& b )
//{
//    cross[0] = a[1] * b[2] - a[2] * b[1];
//    cross[1] = a[2] * b[0] - a[0] * b[2];
//    cross[2] = a[0] * b[1] - a[1] * b[0];
//}

//***********************************************************************************
//***********************************************************************************
//void  LineLoad3DCondition::ExpandReducedMatrix( Matrix& Destination, Matrix& ReducedMatrix )
//{
//    KRATOS_TRY
//    unsigned int size = ReducedMatrix.size2();

//    for ( unsigned int i = 0; i < size; i++ )
//    {
//        int rowindex = i * 3;

//        for ( unsigned int j = 0; j < size; j++ )
//        {
//            unsigned int colindex = j * 3;

//            for ( unsigned int ii = 0; ii < 3; ii++ )
//                Destination( rowindex + ii, colindex + ii ) += ReducedMatrix( i, j );
//        }
//    }

//    KRATOS_CATCH( "" )
//}

//***********************************************************************************
//***********************************************************************************
//void  LineLoad3DCondition::SubtractMatrix( MatrixType& Destination,
//                                   boost::numeric::ublas::bounded_matrix<double, 3, 3>& InputMatrix,
//                                   int InitialRow, int InitialCol )
//{
//    KRATOS_TRY

//    for ( unsigned int i = 0; i < 3; i++ )
//        for ( unsigned int j = 0; j < 3; j++ )
//            Destination( InitialRow + i, InitialCol + j ) -= InputMatrix( i, j );

//    KRATOS_CATCH( "" )
//}

//***********************************************************************************
//***********************************************************************************
//void LineLoad3DCondition::CalculateAndAdd_PressureForce( VectorType& residualvector,
//        const Vector& N,
//        const array_1d<double, 3>& v3,
//        double pressure,
//        double weight,
//        const ProcessInfo& rCurrentProcessInfo )
//{
//    KRATOS_TRY
//    unsigned int number_of_nodes = GetGeometry().size();

//    for ( unsigned int i = 0; i < number_of_nodes; i++ )
//    {
//        int index = 3 * i;
//        double coeff = pressure * N[i] * weight;
//        residualvector[index]   += coeff * v3[0];
//        residualvector[index+1] += coeff * v3[1];
//        residualvector[index+2] += coeff * v3[2];
//    }

//    KRATOS_CATCH( "" )
//}

//***********************************************************************************
//***********************************************************************************
//void LineLoad3DCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
//                                VectorType& rRightHandSideVector,
//                                const ProcessInfo& rCurrentProcessInfo,
//                                bool CalculateStiffnessMatrixFlag,
//                                bool CalculateResidualVectorFlag )
//{
//    KRATOS_TRY

//    const unsigned int number_of_nodes = GetGeometry().size();
//    unsigned int MatSize = number_of_nodes * 3;
//    const unsigned int dimension = 3;
//    //resizing as needed the LHS

//    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
//    {
//        if ( rLeftHandSideMatrix.size1() != MatSize )
//            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

//        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
//    }

//    //resizing as needed the RHS
//    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
//    {
//        if ( rRightHandSideVector.size() != MatSize )
//            rRightHandSideVector.resize( MatSize, false );

//        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
//    }

//    //reading integration points and local gradients
//    const GeometryType::IntegrationPointsArrayType& integration_points =
//        GetGeometry().IntegrationPoints();

//    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer =
//        GetGeometry().ShapeFunctionsLocalGradients();

//    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

//    //calculating actual jacobian
//    GeometryType::JacobiansType J;

//    J = GetGeometry().Jacobian( J );

//    //auxiliary terms
//    array_1d<double, 3> BodyForce;

//    //this should be used once the implementation of the Jacobian is correct in all geometries
////         array_1d<double,3> ge;
////         array_1d<double,3> gn;
////         array_1d<double,3> v3;
//    //loop over integration points
//    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
//    {
//        Vector Load( 3 );
//        noalias( Load ) = ZeroVector( 3 );
//        Vector temp( 3 );

//        for ( unsigned int n = 0; n < GetGeometry().size(); n++ )
//        {
//            noalias( temp ) = ( GetGeometry()[n] ).GetSolutionStepValue( FACE_LOAD );

//            for ( unsigned int i = 0; i < 3; i++ )
//            {
//                Load( i ) += temp( i ) * Ncontainer( PointNumber, n );
//            }
//        }

//        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

//        //to be replaced by the formulation in Face3D
//        Vector t1 = ZeroVector( 3 );//first tangential vector
//        Vector t2 = ZeroVector( 3 );//second tangential vector

//        for ( unsigned int n = 0; n < GetGeometry().size(); n++ )
//        {
//            t1[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 0 );
//            t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 0 );
//            t1[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 0 );
//            t2[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 1 );
//            t2[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 1 );
//            t2[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 1 );
//        }

//        //calculating normal
//        Vector v3 = ZeroVector( 3 );

//        v3[0] = t1[1] * t2[2] - t1[2] * t2[1];

//        v3[1] = t1[2] * t2[0] - t1[0] * t2[2];

//        v3[2] = t1[0] * t2[1] - t1[1] * t2[0];

//        double dA = sqrt( v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2] );

//        // RIGHT HAND SIDE VECTOR
//        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
//        {
//            for ( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
//                for ( unsigned int i = 0; i < 3; i++ )
//                    rRightHandSideVector( prim*dimension + i ) +=
//                        Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight * dA;
//        }
//    }

//    KRATOS_CATCH( "" )
//}

/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int LineLoad3DCondition::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

//***********************************************************************************
//***********************************************************************************
} // Namespace Kratos.
