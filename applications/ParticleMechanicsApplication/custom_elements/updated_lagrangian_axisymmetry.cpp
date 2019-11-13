//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//


// System includes
#include <omp.h>
#include <sstream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_axisymmetry.hpp"
#include "utilities/math_utils.h"
#include "includes/constitutive_law.h"
#include "custom_utilities/particle_mechanics_math_utilities.h"
#include "particle_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( )
    : UpdatedLagrangian( )
{
    //DO NOT CALL IT: only needed for Register and Serialization!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( IndexType NewId, GeometryType::Pointer pGeometry )
        : UpdatedLagrangian( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : UpdatedLagrangian( NewId, pGeometry, pProperties )
{
mFinalizedStep = true;
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::UpdatedLagrangianAxisymmetry( UpdatedLagrangianAxisymmetry const& rOther)
    :UpdatedLagrangian(rOther)
{
}

//*********************************OPERATIONS*****************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianAxisymmetry::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Element::Pointer( new UpdatedLagrangianAxisymmetry( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

//************************************CLONE*******************************************
//************************************************************************************

Element::Pointer UpdatedLagrangianAxisymmetry::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{

    UpdatedLagrangianAxisymmetry NewElement (NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    return Element::Pointer( new UpdatedLagrangianAxisymmetry(NewElement) );
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

UpdatedLagrangianAxisymmetry::~UpdatedLagrangianAxisymmetry()
{
}

//************************************************************************************
void UpdatedLagrangianAxisymmetry::Initialize()
{
    KRATOS_TRY
    UpdatedLagrangian::Initialize();

    const array_1d<double,3>& xg = this->GetValue(MP_COORD);
    const double mp_volume = this->GetValue(MP_VOLUME);
    const double pi = std::atan(1.0)*4.0;
    const double mp_mass = mp_volume * 2 * pi * xg[0] * GetProperties()[DENSITY];
    this->SetValue(MP_MASS, mp_mass);

    mDeterminantF0 = 1;
    mDeformationGradientF0 = IdentityMatrix(3);

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::InitializeGeneralVariables (GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    const unsigned int strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    rVariables.detF  = 1;

    rVariables.detF0 = 1;

    rVariables.detFT = 1;

    rVariables.B.resize( strain_size , number_of_nodes * dimension, false );

    rVariables.F.resize( 3, 3, false );
    rVariables.F = IdentityMatrix(3);

    rVariables.F0.resize( 3, 3, false );
    rVariables.F0 = IdentityMatrix(3);

    rVariables.FT.resize( 3, 3, false );
    rVariables.FT = IdentityMatrix(3);

    rVariables.ConstitutiveMatrix.resize( strain_size, strain_size, false );

    rVariables.StrainVector.resize( strain_size, false );

    rVariables.StressVector.resize( strain_size, false );

    rVariables.DN_DX.resize( number_of_nodes, dimension, false );
    rVariables.DN_De.resize( number_of_nodes, dimension, false );

    const array_1d<double,3>& xg = this->GetValue(MP_COORD);

    rVariables.N = this->MPMShapeFunctionPointValues(rVariables.N, xg);

    // Reading shape functions local gradients
    rVariables.DN_De = this->MPMShapeFunctionsLocalGradients( rVariables.DN_De);

    // CurrentDisp is the unknown variable. It represents the nodal delta displacement. When it is predicted is equal to zero.
    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateKinematics(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Define the stress measure
    rVariables.StressMeasure = ConstitutiveLaw::StressMeasure_Cauchy;

    // Calculating the reference jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n/d£]
    const array_1d<double,3>& xg = this->GetValue(MP_COORD);
    Matrix Jacobian;
    Jacobian = this->MPMJacobian( Jacobian, xg);

    // Calculating the inverse of the jacobian and the parameters needed [d£/dx_n]
    Matrix InvJ;
    double detJ;
    MathUtils<double>::InvertMatrix( Jacobian, InvJ, detJ);

    // Calculating the current jacobian from cartesian coordinates to parent coordinates for the MP element [dx_n+1/d£]
    Matrix jacobian;
    jacobian = this->MPMJacobianDelta( jacobian, xg, rVariables.CurrentDisp);

    //Calculating the inverse of the jacobian and the parameters needed [d£/(dx_n+1)]
    Matrix Invj;
    MathUtils<double>::InvertMatrix( jacobian, Invj, detJ); //overwrites detJ

    // Compute cartesian derivatives [dN/dx_n]
    rVariables.DN_DX = prod( rVariables.DN_De, InvJ);

    // Compute radius
    const double current_radius = ParticleMechanicsMathUtilities<double>::CalculateRadius(rVariables.N, GetGeometry());
    const double initial_radius = ParticleMechanicsMathUtilities<double>::CalculateRadius(rVariables.N, GetGeometry(), Initial);

    rVariables.CurrentDisp = CalculateCurrentDisp(rVariables.CurrentDisp, rCurrentProcessInfo);
    this->CalculateDeformationGradient (rVariables.DN_DX, rVariables.F, rVariables.CurrentDisp, current_radius, initial_radius);

    // Compute cartesian derivatives [dN/dx_n+1]
    rVariables.DN_DX = prod( rVariables.DN_De, Invj); //overwrites DX now is the current position dx

    rVariables.CurrentRadius   =  current_radius;
    rVariables.ReferenceRadius =  initial_radius;

    // Determinant of the previous Deformation Gradient F_n
    rVariables.detF0 = mDeterminantF0;
    rVariables.F0    = mDeformationGradientF0;

    // Compute the deformation matrix B
    this->CalculateDeformationMatrix(rVariables.B, rVariables.F, rVariables.DN_DX, rVariables.N);

    KRATOS_CATCH( "" )
}

//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateDeformationMatrix(Matrix& rB,
        Matrix& rF,
        Matrix& rDN_DX,
        Vector& rN)
{
    KRATOS_TRY

    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension       = r_geometry.WorkingSpaceDimension();
    const double radius = ParticleMechanicsMathUtilities<double>::CalculateRadius(rN, r_geometry);

    rB.clear(); // Set all components to zero

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        const unsigned int index = dimension * i;

        rB( 0, index + 0 ) = rDN_DX( i, 0 );
        rB( 1, index + 1 ) = rDN_DX( i, 1 );
        rB( 2, index + 0 ) = rN[i]/radius;
        rB( 3, index + 0 ) = rDN_DX( i, 1 );
        rB( 3, index + 1 ) = rDN_DX( i, 0 );
    }

    KRATOS_CATCH( "" )
}


//*************************COMPUTE DEFORMATION GRADIENT*******************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateDeformationGradient(const Matrix& rDN_DX,
        Matrix&  rF,
        Matrix&  rDeltaPosition,
        const double & rCurrentRadius,
        const double & rReferenceRadius)
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.PointsNumber();
    const unsigned int dimension       = r_geometry.WorkingSpaceDimension();

    rF = IdentityMatrix( 3 );

    if( dimension == 2 )
    {
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            rF ( 0 , 0 ) += rDeltaPosition(i,0)*rDN_DX ( i , 0 );
            rF ( 0 , 1 ) += rDeltaPosition(i,0)*rDN_DX ( i , 1 );
            rF ( 1 , 0 ) += rDeltaPosition(i,1)*rDN_DX ( i , 0 );
            rF ( 1 , 1 ) += rDeltaPosition(i,1)*rDN_DX ( i , 1 );
        }

        rF ( 2 , 2 ) = rCurrentRadius/rReferenceRadius;
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong!" << std::endl;
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
        GeneralVariables& rVariables,
        const double& rIntegrationWeight)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();

    // Axisymmetric geometric matrix
    double alpha_1 = 0;
    double alpha_2 = 0;
    double alpha_3 = 0;

    unsigned int index_i = 0;

    const double radius = ParticleMechanicsMathUtilities<double>::CalculateRadius(rVariables.N, GetGeometry());

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index_j =0;
        for ( unsigned int j = 0; j < number_of_nodes; j++ )
        {
            alpha_1 = rVariables.DN_DX(j,0) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[0] + rVariables.DN_DX(i,1) * rVariables.StressVector[3] );
            alpha_2 = rVariables.DN_DX(j,1) * ( rVariables.DN_DX(i,0) * rVariables.StressVector[3] + rVariables.DN_DX(i,1) * rVariables.StressVector[1] );
            alpha_3 = rVariables.N[i] * rVariables.N[j] * rVariables.StressVector[2] * (1.0/radius*radius);

            rLeftHandSideMatrix(index_i,index_j)     += (alpha_1 + alpha_2 + alpha_3) * rIntegrationWeight ;
            rLeftHandSideMatrix(index_i+1,index_j+1) += (alpha_1 + alpha_2) * rIntegrationWeight ;

            index_j+=2;
        }

        index_i+=2;

    }

    KRATOS_CATCH( "" )
}

//*************************COMPUTE ALMANSI STRAIN*************************************
//************************************************************************************
void UpdatedLagrangianAxisymmetry::CalculateAlmansiStrain(const Matrix& rF,
    Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Left Cauchy-Green Calculation
    Matrix LeftCauchyGreen = prod( rF, trans( rF ) );

    // Calculating the inverse of the jacobian
    Matrix InverseLeftCauchyGreen = ZeroMatrix( dimension, dimension );
    double det_b=0;
    MathUtils<double>::InvertMatrix( LeftCauchyGreen, InverseLeftCauchyGreen, det_b);

    if( dimension == 2 )
    {
        // Almansi Strain Calculation
        rStrainVector[0] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 0, 0 ) );
        rStrainVector[1] = 0.5 * (  1.00 - InverseLeftCauchyGreen( 1, 1 ) );
        rStrainVector[2] = - InverseLeftCauchyGreen( 0, 1 ); // xy
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong!" << std::endl;
    }


    KRATOS_CATCH( "" )
}

//*************************COMPUTE GREEN-LAGRANGE STRAIN*************************************
//************************************************************************************
// Green-Lagrange Strain: E = 0.5 * (U^2 - I) = 0.5 * (C - I)
void UpdatedLagrangianAxisymmetry::CalculateGreenLagrangeStrain(const Matrix& rF,
    Vector& rStrainVector )
{
    KRATOS_TRY

    const unsigned int dimension  = GetGeometry().WorkingSpaceDimension();

    // Right Cauchy-Green Calculation
    Matrix C ( dimension, dimension );
    noalias( C ) = prod( trans( rF ), rF );

    if( dimension == 2 )
    {
        //Green Lagrange Strain Calculation
        if ( rStrainVector.size() != 3 ) rStrainVector.resize( 3, false );
        rStrainVector[0] = 0.5 * ( C( 0, 0 ) - 1.00 );
        rStrainVector[1] = 0.5 * ( C( 1, 1 ) - 1.00 );
        rStrainVector[2] = C( 0, 1 ); // xy
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong!" << std::endl;
    }

    KRATOS_CATCH( "" )
}


//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    GeometryType& r_geometry = GetGeometry();
    const int number_of_nodes = r_geometry.size();
    const int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( rResult.size() != matrix_size )
        rResult.resize( matrix_size, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index] = r_geometry[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = r_geometry[i].GetDof( DISPLACEMENT_Y ).EquationId();
    }

}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::GetDofList( DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo )
{
    rElementalDofList.resize( 0 );

    GeometryType& r_geometry = GetGeometry();
    for ( unsigned int i = 0; i < r_geometry.size(); i++ )
    {
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( r_geometry[i].pGetDof( DISPLACEMENT_Y ) );
    }
}

//************************************************************************************
//************************************************************************************

// Function that return Jacobian matrix
Matrix& UpdatedLagrangianAxisymmetry::MPMJacobian( Matrix& rResult, const array_1d<double,3>& rPoint)
{

    KRATOS_TRY

    // Derivatives of shape functions
    Matrix shape_functions_gradients;
    shape_functions_gradients =this->MPMShapeFunctionsLocalGradients(
                                    shape_functions_gradients);

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int number_nodes = r_geometry.PointsNumber();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    if (dimension ==2)
    {
        rResult.resize( 2, 2, false );
        rResult = ZeroMatrix(2,2);

        for ( unsigned int i = 0; i < number_nodes; i++ )
        {
            rResult( 0, 0 ) += ( r_geometry.GetPoint( i ).X() *  shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( r_geometry.GetPoint( i ).X() *  shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( r_geometry.GetPoint( i ).Y() *  shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( r_geometry.GetPoint( i ).Y() *  shape_functions_gradients( i, 1 ) );
        }

    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong!" << std::endl;
    }

    return rResult;

    KRATOS_CATCH( "" )
}

/**
     * Jacobian in given point and given a delta position. This method calculate jacobian
     * matrix in given point and a given delta position.
     *
     * @param rPoint point which jacobians has to
* be calculated in it.
*
* @return Matrix of double which is jacobian matrix \f$ J \f$ in given point and a given delta position.
*
* @see DeterminantOfJacobian
* @see InverseOfJacobian
    */
Matrix& UpdatedLagrangianAxisymmetry::MPMJacobianDelta( Matrix& rResult, const array_1d<double,3>& rPoint, const Matrix & rDeltaPosition )
{
    KRATOS_TRY

    // Derivatives of shape functions
    Matrix shape_functions_gradients;
    shape_functions_gradients = this->MPMShapeFunctionsLocalGradients(
                                    shape_functions_gradients );

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    if (dimension ==2)
    {
        rResult.resize( 2, 2, false );
        rResult = ZeroMatrix(2,2);

        for ( unsigned int i = 0; i < r_geometry.size(); i++ )
        {
            rResult( 0, 0 ) += ( r_geometry.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 0, 1 ) += ( r_geometry.GetPoint( i ).X() + rDeltaPosition(i,0)) * ( shape_functions_gradients( i, 1 ) );
            rResult( 1, 0 ) += ( r_geometry.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 0 ) );
            rResult( 1, 1 ) += ( r_geometry.GetPoint( i ).Y() + rDeltaPosition(i,1)) * ( shape_functions_gradients( i, 1 ) );
        }
    }
    else
    {
        KRATOS_ERROR <<  "Dimension given is wrong!" << std::endl;
    }

    return rResult;

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

/**
     * Shape function values in given point. This method calculate the shape function
     * vector in given point.
     *
     * @param rPoint point which shape function values have to
* be calculated in it.
*
* @return Vector of double which is shape function vector \f$ N \f$ in given point.
*
    */
Vector& UpdatedLagrangianAxisymmetry::MPMShapeFunctionPointValues( Vector& rResult, const array_1d<double,3>& rPoint )
{
    KRATOS_TRY

    const GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    array_1d<double,3> rPointLocal = ZeroVector(3);
    rPointLocal = r_geometry.PointLocalCoordinates(rPointLocal, rPoint);

    // For Triangular 2D
    if (dimension == 2 && r_geometry.PointsNumber() == 3)
    {
        rResult.resize(3, false);
        rResult[0] = 1 - rPointLocal[0] - rPointLocal[1] ;
        rResult[1] = rPointLocal[0] ;
        rResult[2] = rPointLocal[1];
    }
    // For Quadrilateral 2D
    else if (dimension == 2 && r_geometry.PointsNumber() == 4)
    {
        rResult.resize(4, false);
        rResult[0] = 0.25 * (1 - rPointLocal[0]) * (1 - rPointLocal[1]) ;
        rResult[1] = 0.25 * (1 + rPointLocal[0]) * (1 - rPointLocal[1]) ;
        rResult[2] = 0.25 * (1 + rPointLocal[0]) * (1 + rPointLocal[1]) ;
        rResult[3] = 0.25 * (1 - rPointLocal[0]) * (1 + rPointLocal[1]) ;
    }

    return rResult;

    KRATOS_CATCH( "" )
}

// Function which return dN/de
Matrix& UpdatedLagrangianAxisymmetry::MPMShapeFunctionsLocalGradients( Matrix& rResult )
{
    const GeometryType& r_geometry = GetGeometry();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    array_1d<double,3> rPointLocal = ZeroVector(3);

    const array_1d<double,3>& xg = this->GetValue(MP_COORD);
    rPointLocal = r_geometry.PointLocalCoordinates(rPointLocal, xg);

    if (dimension == 2 && r_geometry.PointsNumber() == 3)
    {
        rResult = ZeroMatrix( 3, 2 );
        rResult( 0, 0 ) = -1.0;
        rResult( 0, 1 ) = -1.0;
        rResult( 1, 0 ) =  1.0;
        rResult( 1, 1 ) =  0.0;
        rResult( 2, 0 ) =  0.0;
        rResult( 2, 1 ) =  1.0;
    }
    else if (dimension == 2 && r_geometry.PointsNumber() == 4)
    {
        rResult = ZeroMatrix( 4, 2 );
        rResult( 0, 0 ) = -0.25 * (1 - rPointLocal[1]);
        rResult( 0, 1 ) = -0.25 * (1 - rPointLocal[0]);
        rResult( 1, 0 ) = 0.25 * (1 - rPointLocal[1]);
        rResult( 1, 1 ) = -0.25 * (1 + rPointLocal[0]);
        rResult( 2, 0 ) = 0.25 * (1 + rPointLocal[1]);
        rResult( 2, 1 ) = 0.25 * (1 + rPointLocal[0]);
        rResult( 3, 0 ) = -0.25 * (1 + rPointLocal[1]);
        rResult( 3, 1 ) = 0.25 * (1 - rPointLocal[0]);
    }

    return rResult;
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::GetValuesVector( Vector& values, int Step )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( DISPLACEMENT_Y, Step );
    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::GetFirstDerivativesVector( Vector& values, int Step )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( VELOCITY_Y, Step );

    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::GetSecondDerivativesVector( Vector& values, int Step )
{
    GeometryType& r_geometry = GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();
    unsigned int matrix_size = number_of_nodes * dimension;

    if ( values.size() != matrix_size ) values.resize( matrix_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        unsigned int index = i * dimension;
        values[index] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_X, Step );
        values[index + 1] = r_geometry[i].FastGetSolutionStepValue( ACCELERATION_Y, Step );

    }
}

//************************************************************************************
//************************************************************************************

void UpdatedLagrangianAxisymmetry::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    rSerializer.save("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.save("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.save("DeterminantF0",mDeterminantF0);
}

void UpdatedLagrangianAxisymmetry::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    rSerializer.load("ConstitutiveLawVector",mConstitutiveLawVector);
    rSerializer.load("DeformationGradientF0",mDeformationGradientF0);
    rSerializer.load("DeterminantF0",mDeterminantF0);
}

} // Namespace Kratos

