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
felix.nagel@rub.de
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
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:55:58 $
//   Revision:            $Revision: 1.11 $
//
//


// System includes
//#include <sys/sysinfo.h>

// External includes


// Project includes
#include "custom_elements/unsaturated_soils_element_3phase_small_strain.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/prism_3d_6.h"
#include "structural_application.h"

#include "boost/timer.hpp"

namespace Kratos
{
    UnsaturatedSoilsElement_3phase_SmallStrain::UnsaturatedSoilsElement_3phase_SmallStrain( IndexType NewId,
            GeometryType::Pointer pGeometry )
            : Element( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    UnsaturatedSoilsElement_3phase_SmallStrain::UnsaturatedSoilsElement_3phase_SmallStrain( IndexType NewId,
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
            : Element( NewId, pGeometry, pProperties )
    {
        //setting up the nodal degrees of freedom
        //with DISPLACEMENT_X, DISPLACEMENT_Y, DISPLACEMENT_Z
        //DOFs at the end of time step
        //All calculations are made on the general midpoint alpha
        // Variables DOF_ALPHA are updated in the scheme
        if ( GetGeometry().size() == 27 || GetGeometry().size() == 20 || GetGeometry().size() == 10 || GetGeometry().size() == 15 || GetGeometry().size() == 8 )
        {
            if ( GetGeometry().size() == 27 )
            {
                mNodesPressMin = 1;
                mNodesPressMax = 8;
                mNodesDispMin = 1;
                mNodesDispMax = 27;
                mpPressureGeometry = Geometry< Node<3> >::Pointer( new Hexahedra3D8 <Node<3> >(
                                         GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ),
                                         GetGeometry()( 4 ), GetGeometry()( 5 ), GetGeometry()( 6 ), GetGeometry()( 7 ) ) );
                mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            }

            if ( GetGeometry().size() == 20 )
            {
                mNodesPressMin = 1;
                mNodesPressMax = 8;
                mNodesDispMin = 1;
                mNodesDispMax = 20;
                mpPressureGeometry = Geometry< Node<3> >::Pointer( new Hexahedra3D8 <Node<3> >(
                                         GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ),
                                         GetGeometry()( 4 ), GetGeometry()( 5 ), GetGeometry()( 6 ), GetGeometry()( 7 ) ) );
                mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            }

            if ( GetGeometry().size() == 10 )
            {
                mNodesPressMin = 1;
                mNodesPressMax = 4;
                mNodesDispMin = 1;
                mNodesDispMax = 10;
                mpPressureGeometry = Geometry< Node<3> >::Pointer( new Tetrahedra3D4 <Node<3> >(
                                         GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ) ) );
                mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
            }

            if ( GetGeometry().size() == 15 )
            {
                mNodesPressMin = 1;
                mNodesPressMax = 6;
                mNodesDispMin = 1;
                mNodesDispMax = 15;
                mpPressureGeometry = Geometry< Node<3> >::Pointer( new Prism3D6 <Node<3> >(
                                         GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ), GetGeometry()( 4 ), GetGeometry()( 5 ) ) );
                mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            }

            // Attention this version does not fulfill the Babuska Brezzi Stability constraint and is therefore not an appropriate choice if consolidation under initially undrained conditions is analysed
            if ( GetGeometry().size() == 8 )
            {
                mNodesPressMin = 1;
                mNodesPressMax = 8;
                mNodesDispMin = 1;
                mNodesDispMax = 8;
                mpPressureGeometry = Geometry< Node<3> >::Pointer( new Hexahedra3D8 <Node<3> >(
                                         GetGeometry()( 0 ), GetGeometry()( 1 ), GetGeometry()( 2 ), GetGeometry()( 3 ),
                                         GetGeometry()( 4 ), GetGeometry()( 5 ), GetGeometry()( 6 ), GetGeometry()( 7 ) ) );
                mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
            }
        }
        else
            KRATOS_ERROR( std::logic_error, "This element matches only with a quadratic hexaeder (20 or 27), tetraeder (10) or prism (15) geometry" , *this );

    }

    Element::Pointer UnsaturatedSoilsElement_3phase_SmallStrain::Create( IndexType NewId,
            NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new UnsaturatedSoilsElement_3phase_SmallStrain( NewId, GetGeometry().Create( ThisNodes ),
                                 pProperties ) );
    }

    UnsaturatedSoilsElement_3phase_SmallStrain::~UnsaturatedSoilsElement_3phase_SmallStrain()
    {
    }


    //************************************************************************************
    //************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::Initialize()
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        mInvJ0.resize( integration_points.size() );

        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
            mInvJ0[i].resize( dim, dim, false );
            noalias( mInvJ0[i] ) = ZeroMatrix( dim, dim );
        }

        mDetJ0.resize( integration_points.size(), false );

        noalias( mDetJ0 ) = ZeroVector( integration_points.size() );

        GeometryType::JacobiansType J0( integration_points.size() );

        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );

        //calculating the inverse J0

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );
        }

        //Constitutive Law initialisation
        if ( mConstitutiveLawVector.size() != integration_points.size() )
        {
            mConstitutiveLawVector.resize( integration_points.size() );
            InitializeMaterial();
        }

        for ( unsigned int i = ( mNodesDispMin - 1 ) ; i < mNodesDispMax ; i++ )
        {
            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_NULL ) =
                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT );
            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_EINS ) =
                ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT );
            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_DT ) = ZeroVector( 3 );
            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_NULL_DT ) = ZeroVector( 3 );
            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_EINS_DT ) = ZeroVector( 3 );
            ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION ) = ZeroVector( 3 );
            ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION_NULL ) = ZeroVector( 3 );
            ( GetGeometry()[i] ).GetSolutionStepValue( ACCELERATION_EINS ) = ZeroVector( 3 );
            ( GetGeometry()[i] ).GetSolutionStepValue( DISPLACEMENT_OLD ) = ZeroVector( 3 );
        }

        for ( unsigned int i = ( mNodesPressMin - 1 ) ; i < mNodesPressMax ; i++ )
        {
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL ) =
                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE );
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS ) =
                ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE );
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_DT ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL_DT ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS_DT ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_ACCELERATION ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_NULL_ACCELERATION ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( WATER_PRESSURE_EINS_ACCELERATION ) = 0;

            ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE_NULL ) =
                ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE );
            ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE_EINS ) =
                ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE );
            ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE_DT ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE_NULL_DT ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE_EINS_DT ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE_ACCELERATION ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE_NULL_ACCELERATION ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( AIR_PRESSURE_EINS_ACCELERATION ) = 0;

            ( GetGeometry()[i] ).GetSolutionStepValue( REACTION_WATER_PRESSURE ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( REACTION_AIR_PRESSURE ) = 0;
            ( GetGeometry()[i] ).GetSolutionStepValue( PRESSURE ) = 0;
//             (GetGeometry()[i]).GetSolutionStepValue(DARCY_VELO_WATER)=ZeroVector(3);
        }

        KRATOS_CATCH( "" )
    }

    /**
    * THIS method is called from the scheme at the start of each solution step
    * @param rCurrentProcessInfo
    */
    void UnsaturatedSoilsElement_3phase_SmallStrain::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        for ( unsigned int Point = 0; Point < integration_points.size(); Point++ )
        {
            mConstitutiveLawVector[Point]->InitializeSolutionStep( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point ), CurrentProcessInfo );

//            mConstitutiveLawVector[Point]->SetValue(SUCTION, GetGeometry()[0].GetSolutionStepValue(TEMPERATURE),CurrentProcessInfo );
        }
    }

    //************************************************************************************
    //************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateAll( MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
            bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
    {
        KRATOS_TRY

        unsigned int number_of_nodes_disp = ( mNodesDispMax - mNodesDispMin + 1 );
        unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        if ( dim != 3 ) KRATOS_ERROR( std::logic_error, "UnsaturatedSoilsElement_3phase_SmallStrain cannot be used for other than 3D problems" , "" );

        Matrix B( 6, 6 );

        Matrix TanC_U( 6, 6 );

        Vector StrainVector( 6 );

        Vector StressVector( 6 );

        Matrix DN_DX_DISP( number_of_nodes_disp, dim );

        Matrix CurrentDisp( number_of_nodes_disp, dim );

        //resizing as needed the LHS
        unsigned int MatSize1 = ( number_of_nodes_disp * dim + number_of_nodes_press * 2 );

        unsigned int MatSizeU = number_of_nodes_disp * dim;

        unsigned int MatSizeP = number_of_nodes_press;

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != MatSize1 )
                rLeftHandSideMatrix.resize( MatSize1, MatSize1, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize1, MatSize1 ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            if ( rRightHandSideVector.size() != MatSize1 )
                rRightHandSideVector.resize( MatSize1, false );

            noalias( rRightHandSideVector ) = ZeroVector( MatSize1 ); //resetting RHS
        }

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =
            GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De_Pressure =
            mpPressureGeometry->ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer_Displacement = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        const Matrix& Ncontainer_Pressure = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

        double Weight;

        double capillaryPressure;

        double waterPressure;

        double airPressure;

        double capillaryPressure_Dt;

        double waterPressure_Dt;

        double airPressure_Dt;

        double porosity;

        double density;

        double DetJ = 0.0;

        Matrix Help_K_UU( MatSizeU, MatSizeU );

        Matrix Help_K_UW( MatSizeU, MatSizeP );

        Matrix Help_K_UA( MatSizeU, MatSizeP );

        Matrix Help_K_WU( MatSizeP, MatSizeU );

        Matrix Help_K_WW( MatSizeP, MatSizeP );

        Matrix Help_K_WA( MatSizeP, MatSizeP );

        Matrix Help_K_AU( MatSizeP, MatSizeU );

        Matrix Help_K_AW( MatSizeP, MatSizeP );

        Matrix Help_K_AA( MatSizeP, MatSizeP );

        Vector Help_R_U( MatSizeU );

        Vector Help_R_W( MatSizeP );

        Vector Help_R_A( MatSizeP );

//                 Matrix B_Operator(6,number_of_nodes_disp*3);
//                 noalias(B_Operator)= ZeroMatrix(6,number_of_nodes_disp*3);
//                 Matrix DN_DX_DISP(number_of_nodes_disp,3);
        Matrix DN_DX_PRESS( number_of_nodes_press, 3 );

        Vector N_DISP( number_of_nodes_disp );

        Vector N_PRESS( number_of_nodes_press );

//                 Matrix tanC_U(6,6);
        Matrix tanC_W( 3, 3 );

        Matrix tanC_A( 3, 3 );

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            noalias( Help_K_UU ) = ZeroMatrix( MatSizeU, MatSizeU );
            noalias( Help_K_UW ) = ZeroMatrix( MatSizeU, MatSizeP );
            noalias( Help_K_UA ) = ZeroMatrix( MatSizeU, MatSizeP );

            noalias( Help_K_WU ) = ZeroMatrix( MatSizeP, MatSizeU );
            noalias( Help_K_WW ) = ZeroMatrix( MatSizeP, MatSizeP );
            noalias( Help_K_WA ) = ZeroMatrix( MatSizeP, MatSizeP );

            noalias( Help_K_AU ) = ZeroMatrix( MatSizeP, MatSizeU );
            noalias( Help_K_AW ) = ZeroMatrix( MatSizeP, MatSizeP );
            noalias( Help_K_AA ) = ZeroMatrix( MatSizeP, MatSizeP );

        }

        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            noalias( Help_R_U ) = ZeroVector( MatSizeU );
            noalias( Help_R_W ) = ZeroVector( MatSizeP );
            noalias( Help_R_A ) = ZeroVector( MatSizeP );
        }

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
        {
            CurrentDisp( node, 0 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_X );
            CurrentDisp( node, 1 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_Y );
            CurrentDisp( node, 2 ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT_Z );
        }

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space sum_(beta=0)^(number of quadrature points)
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            noalias( DN_DX_PRESS ) = prod( DN_De_Pressure[PointNumber], mInvJ0[PointNumber] );

            noalias( DN_DX_DISP ) = prod( DN_De_Displacement[PointNumber], mInvJ0[PointNumber] );

            Weight = integration_points[PointNumber].Weight();

            DetJ = mDetJ0[PointNumber];

            // Shape Functions on current spatial quadrature point
            noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

            noalias( N_DISP ) = row( Ncontainer_Displacement, PointNumber );

            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, DN_DX_DISP );
            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );

            GetPressures( N_PRESS, capillaryPressure, waterPressure,
                          airPressure );

            GetDerivativeDPressuresDt( N_PRESS, capillaryPressure_Dt, waterPressure_Dt, airPressure_Dt );

            porosity = GetPorosity( DN_DX_DISP );

            density = GetAveragedDensity( capillaryPressure, airPressure, porosity );

            CalculateStressAndTangentialStiffnessUnsaturatedSoils( StressVector, TanC_U, tanC_W, tanC_A, StrainVector, waterPressure, airPressure, PointNumber, rCurrentProcessInfo );

            if ( CalculateStiffnessMatrixFlag == true )
            {
                //Calculation of spatial Stiffnes and Mass Matrix

                CalculateStiffnesMatrixUU( Help_K_UU, TanC_U, B, DN_DX_DISP, N_DISP, density, capillaryPressure, airPressure, Weight, DetJ );

                CalculateStiffnesMatrixUW( Help_K_UW, tanC_W, DN_DX_DISP, N_DISP,
                                           N_PRESS, capillaryPressure, airPressure, Weight, DetJ );

                CalculateStiffnesMatrixUA( Help_K_UA, tanC_A, DN_DX_DISP, N_DISP,
                                           N_PRESS, capillaryPressure, airPressure, Weight, DetJ );

                CalculateStiffnesMatrixWU( Help_K_WU, DN_DX_DISP, DN_DX_PRESS,
                                           N_PRESS, capillaryPressure,  Weight, DetJ );

                CalculateStiffnesMatrixWW( Help_K_WW, DN_DX_DISP,
                                           DN_DX_PRESS,
                                           N_PRESS, capillaryPressure, Weight, DetJ );

                CalculateStiffnesMatrixWA( Help_K_WA, DN_DX_DISP,
                                           DN_DX_PRESS,
                                           N_PRESS, capillaryPressure, Weight, DetJ );

                CalculateStiffnesMatrixAU( Help_K_AU, DN_DX_DISP, DN_DX_PRESS,
                                           N_PRESS, capillaryPressure, airPressure,
                                           airPressure_Dt, Weight, DetJ );

                CalculateStiffnesMatrixAW( Help_K_AW, DN_DX_DISP,
                                           DN_DX_PRESS,
                                           N_PRESS, capillaryPressure, airPressure,
                                           airPressure_Dt, Weight, DetJ );

                CalculateStiffnesMatrixAA( Help_K_AA, DN_DX_DISP,
                                           DN_DX_PRESS,
                                           N_PRESS, capillaryPressure, airPressure,
                                           airPressure_Dt, Weight, DetJ );
            }

            if ( CalculateResidualVectorFlag == true )
            {
                //Calculation of spatial Loadvector
                AddBodyForcesToRHSVectorU( Help_R_U, N_DISP, density,
                                           Weight, DetJ );
                AddInternalForcesToRHSU( Help_R_U, B, StressVector, Weight, DetJ );
                AddInternalForcesToRHSW( Help_R_W, DN_DX_DISP,
                                         DN_DX_PRESS, N_PRESS, capillaryPressure, Weight,
                                         DetJ );
                AddInternalForcesToRHSA( Help_R_A, DN_DX_DISP,
                                         DN_DX_PRESS, N_PRESS, capillaryPressure, airPressure, airPressure_Dt, Weight,
                                         DetJ );
            }

///////////////////////////////////////////////////////////////////////
// END Integration in space sum_(beta=0)^(number of quadrature points)
///////////////////////////////////////////////////////////////////////
        }

        if ( CalculateStiffnessMatrixFlag == true )
        {

            AssembleTimeSpaceStiffnessFromStiffSubMatrices( rLeftHandSideMatrix,
                    Help_K_UU, Help_K_UW, Help_K_UA,
                    Help_K_WU, Help_K_WW, Help_K_WA,
                    Help_K_AU, Help_K_AW, Help_K_AA );
        }

        //
        if ( CalculateResidualVectorFlag == true )
        {
            AssembleTimeSpaceRHSFromSubVectors( rRightHandSideVector, Help_R_U, Help_R_W, Help_R_A );
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,
                      CalculateResidualVectorFlag );
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    ////************************************************************************************
    ////************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int number_of_nodes_disp = ( mNodesDispMax - mNodesDispMin + 1 );
        unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        Matrix DN_DX_DISP( number_of_nodes_disp, dim );


        //resizing as needed the LHS
        unsigned int MatSize1 = ( number_of_nodes_disp * dim + number_of_nodes_press * 2 );
        unsigned int MatSizeU = number_of_nodes_disp * dim;
        unsigned int MatSizeP = number_of_nodes_press;


        if ( rDampMatrix.size1() != MatSize1 )
            rDampMatrix.resize( MatSize1, MatSize1, false );

        noalias( rDampMatrix ) = ZeroMatrix( MatSize1, MatSize1 ); //resetting LHS

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =
            GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer_Pressure = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

        double Weight;

        double capillaryPressure;

        double waterPressure;

        double airPressure;

        Matrix Help_D_UU( MatSizeU, MatSizeU );

        Matrix Help_D_UW( MatSizeU, MatSizeP );

        Matrix Help_D_UA( MatSizeU, MatSizeP );

        Matrix Help_D_WU( MatSizeP, MatSizeU );

        Matrix Help_D_WW( MatSizeP, MatSizeP );

        Matrix Help_D_WA( MatSizeP, MatSizeP );

        Matrix Help_D_AU( MatSizeP, MatSizeU );

        Matrix Help_D_AW( MatSizeP, MatSizeP );

        Matrix Help_D_AA( MatSizeP, MatSizeP );


//                 Matrix DN_DX_DISP(number_of_nodes_disp,3);

        Vector N_PRESS( number_of_nodes_press );

        noalias( Help_D_UU ) = ZeroMatrix( MatSizeU, MatSizeU );

        noalias( Help_D_UW ) = ZeroMatrix( MatSizeU, MatSizeP );

        noalias( Help_D_UA ) = ZeroMatrix( MatSizeU, MatSizeP );

        noalias( Help_D_WU ) = ZeroMatrix( MatSizeP, MatSizeU );

        noalias( Help_D_WW ) = ZeroMatrix( MatSizeP, MatSizeP );

        noalias( Help_D_WA ) = ZeroMatrix( MatSizeP, MatSizeP );

        noalias( Help_D_AU ) = ZeroMatrix( MatSizeP, MatSizeU );

        noalias( Help_D_AW ) = ZeroMatrix( MatSizeP, MatSizeP );

        noalias( Help_D_AA ) = ZeroMatrix( MatSizeP, MatSizeP );

        double DetJ = 0.0;

        /////////////////////////////////////////////////////////////////////////
//// Integration in space sum_(beta=0)^(number of quadrature points)
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            noalias( DN_DX_DISP ) = prod( DN_De_Displacement[PointNumber], mInvJ0[PointNumber] );
            Weight = integration_points[PointNumber].Weight();
            // Jacobian on current quadrature point

            DetJ = mDetJ0[PointNumber];
            // Shape Functions on current spatial quadrature point
            // Shape Functions on current spatial quadrature point
            noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );


            GetPressures( N_PRESS, capillaryPressure, waterPressure, airPressure );

//                         Calculation of spatial Stiffnes and Mass Matrix
            CalculateDampingMatrixWU( Help_D_WU, DN_DX_DISP, N_PRESS,
                                      capillaryPressure, Weight, DetJ );

            CalculateDampingMatrixWW( Help_D_WW, DN_DX_DISP, N_PRESS,
                                      capillaryPressure, Weight, DetJ );

            CalculateDampingMatrixWA( Help_D_WA, DN_DX_DISP, N_PRESS,
                                      capillaryPressure, Weight, DetJ );

            CalculateDampingMatrixAU( Help_D_AU, DN_DX_DISP, N_PRESS,
                                      capillaryPressure, Weight, DetJ );

            CalculateDampingMatrixAW( Help_D_AW, DN_DX_DISP, N_PRESS,
                                      capillaryPressure, Weight, DetJ );

            CalculateDampingMatrixAA( Help_D_AA, DN_DX_DISP, N_PRESS,
                                      capillaryPressure, airPressure, Weight, DetJ );
//                          CalculateMassMatrix(HelpMassMatrix, N, Weight,DetJ,density);

            ///////////////////////////////////////////////////////////////////////
// END Integration in space sum_(beta=0)^(number of quadrature points)
            ///////////////////////////////////////////////////////////////////////
        }


        AssembleTimeSpaceStiffnessFromDampSubMatrices( rDampMatrix,

                Help_D_UU, Help_D_UW, Help_D_UA,
                Help_D_WU, Help_D_WW, Help_D_WA,
                Help_D_AU, Help_D_AW, Help_D_AA );

        KRATOS_CATCH( "" )
    }

    ////************************************************************************************
    ////************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );
//             unsigned int number_of_nodes_disp = (mNodesDispMax-mNodesDispMin+1);
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        Matrix DN_DX_DISP( GetGeometry().size(), GetGeometry().WorkingSpaceDimension() );

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        Matrix nodesLocalCoords( number_of_nodes_press, dim );
        const GeometryType::ShapeFunctionsGradientsType& DN_DX_DISP_De =
            GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
        const GeometryType::ShapeFunctionsGradientsType& DN_DX_PRESS_De =
            mpPressureGeometry->ShapeFunctionsLocalGradients( mThisIntegrationMethod );

//             Matrix DN_DX_DISP(number_of_nodes_disp,3);
        Matrix DN_DX_PRESS( number_of_nodes_press, 3 );

        nodesLocalCoords( 0, 0 ) = -1.0;
        nodesLocalCoords( 0, 1 ) = -1.0;
        nodesLocalCoords( 0, 2 ) = -1.0;

        nodesLocalCoords( 1, 0 ) = 1.0;
        nodesLocalCoords( 1, 1 ) = -1.0;
        nodesLocalCoords( 1, 2 ) = -1.0;

        nodesLocalCoords( 2, 0 ) = 1.0;
        nodesLocalCoords( 2, 1 ) = 1.0;
        nodesLocalCoords( 2, 2 ) = -1.0;

        nodesLocalCoords( 3, 0 ) = -1.0;
        nodesLocalCoords( 3, 1 ) = 1.0;
        nodesLocalCoords( 3, 2 ) = -1.0;

        nodesLocalCoords( 4, 0 ) = -1.0;
        nodesLocalCoords( 4, 1 ) = -1.0;
        nodesLocalCoords( 4, 2 ) = 1.0;

        nodesLocalCoords( 5, 0 ) = 1.0;
        nodesLocalCoords( 5, 1 ) = -1.0;
        nodesLocalCoords( 5, 2 ) = 1.0;

        nodesLocalCoords( 6, 0 ) = 1.0;
        nodesLocalCoords( 6, 1 ) = 1.0;
        nodesLocalCoords( 6, 2 ) = 1.0;

        nodesLocalCoords( 7, 0 ) = -1.0;
        nodesLocalCoords( 7, 1 ) = 1.0;
        nodesLocalCoords( 7, 2 ) = 1.0;

        for ( unsigned int node = 0; node < number_of_nodes_press; node++ )
        {
            Vector local_coords( 3 );
            local_coords( 0 ) = nodesLocalCoords( node, 0 );
            local_coords( 1 ) = nodesLocalCoords( node, 1 );
            local_coords( 2 ) = nodesLocalCoords( node, 2 );

            noalias( DN_DX_PRESS ) = prod( DN_DX_PRESS_De[0], mInvJ0[0] );
            noalias( DN_DX_DISP ) = prod( DN_DX_DISP_De[0], mInvJ0[0] );

            ( GetGeometry()[node] ).GetSolutionStepValue( VISCOSITY ) = GetFlowWater( DN_DX_PRESS, DN_DX_DISP, ( GetGeometry()[node].GetSolutionStepValue( AIR_PRESSURE ) - GetGeometry()[node].GetSolutionStepValue( WATER_PRESSURE ) ) )( 0 );
        }


        for ( unsigned int Point = 0; Point < integration_points.size(); Point++ )
        {
            mConstitutiveLawVector[Point]->FinalizeSolutionStep( GetProperties(), GetGeometry(), ZeroVector( 0 ), CurrentProcessInfo );
        }

        Vector Dummy_Vector( 9 );

        noalias( Dummy_Vector ) = mConstitutiveLawVector[0]->GetValue( INTERNAL_VARIABLES, Dummy_Vector );

        for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
        {
            GetGeometry()[i].GetSolutionStepValue( MOMENTUM_X ) = Dummy_Vector( 0 );
            GetGeometry()[i].GetSolutionStepValue( MOMENTUM_Y ) = Dummy_Vector( 1 );
            GetGeometry()[i].GetSolutionStepValue( MOMENTUM_Z ) = Dummy_Vector( 2 );
            GetGeometry()[i].GetSolutionStepValue( PRESSURE ) = Dummy_Vector( 3 );
            GetGeometry()[i].GetSolutionStepValue( ERROR_RATIO ) = Dummy_Vector( 4 );
            GetGeometry()[i].GetSolutionStepValue( EXCESS_PORE_WATER_PRESSURE ) = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE ) - 9.81 * 1000.0 * ( 20.0 - GetGeometry()[i].Z() );
        }

//                 Vector Dummy_Vector(8);
//                 noalias(Dummy_Vector)= mConstitutiveLawVector[0]->GetValue(INTERNAL_VARIABLES);
//                 for(unsigned int i=0; i< GetGeometry().size(); i++)
//                 {
//                         GetGeometry()[i].GetSolutionStepValue(MOMENTUM_X)= Dummy_Vector(0);
//                         GetGeometry()[i].GetSolutionStepValue(MOMENTUM_Y)= Dummy_Vector(1);
//                         GetGeometry()[i].GetSolutionStepValue(MOMENTUM_Z)= Dummy_Vector(2);
//                 }
    }

    //************************************************************************************
    //************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateOnIntegrationPoints( const Variable<double >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        if ( Output.size() != integration_points.size() )
            Output.resize( integration_points.size() );

        const Matrix& Ncontainer_Pressure = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

        Vector N_PRESS( number_of_nodes_press );

        double capillaryPressure;

        double waterPressure;

        double airPressure;

        double saturation;

        /////////////////////////////////////////////////////////////////////////
//// Integration in space sum_(beta=0)^(number of quadrature points)
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Shape Functions on current spatial quadrature point
            if ( N_PRESS.size() != number_of_nodes_press )
                N_PRESS.resize( number_of_nodes_press );

            noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

            GetPressures( N_PRESS, capillaryPressure, waterPressure, airPressure );

            saturation = GetSaturation( capillaryPressure );

            if ( rVariable == SATURATION )
            {
                Output[PointNumber] = saturation;
            }

            if ( rVariable == WATER_PRESSURE )
            {
                Output[PointNumber] = waterPressure;
            }

            if ( rVariable == AIR_PRESSURE )
            {
                Output[PointNumber] = airPressure;
            }
        }

        KRATOS_CATCH( "" )
    }

    /**
    * Calculate Vector Variables at each integration point, used for postprocessing etc.
    * @param rVariable Global name of the variable to be calculated
    * @param output Vector to store the values on the qudrature points, output of the method
    * @param rCurrentProcessInfo
    */
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
            std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        GetValueOnIntegrationPoints( rVariable, Output, rCurrentProcessInfo );
    }

    //************************************************************************************
    //************************************************************************************

    inline void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateAndAddExtForceContribution( const Vector& N,
            const ProcessInfo& CurrentProcessInfo, Vector& BodyForce, VectorType& rRightHandSideVector,
            double weight )
    {
        KRATOS_TRY

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::EquationIdVector( EquationIdVectorType& rResult,
            ProcessInfo& CurrentProcessInfo )
    {
        unsigned int dim_press = 2;//two pressure dofs
        unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
        unsigned int MatSize =
            ( mNodesPressMax - mNodesPressMin + 1 ) * dim_press +
            ( mNodesDispMax - mNodesDispMin + 1 ) * dim_disp;

        if ( rResult.size() != MatSize )
            rResult.resize( MatSize, false );

        unsigned int adddisp = ( mNodesPressMax - mNodesPressMin + 1 ) * dim_press;

        for ( unsigned int i = ( mNodesPressMin - 1 );i < mNodesPressMax;i++ )
        {
            int index = ( i - mNodesPressMin + 1 ) * dim_press;
            rResult[index] = GetGeometry()[i].GetDof( WATER_PRESSURE ).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof( AIR_PRESSURE ).EquationId();
        }

        for ( unsigned int i = ( mNodesDispMin - 1 );i < mNodesDispMax;i++ )
        {
            unsigned int index = adddisp + ( i - mNodesDispMin + 1 ) * dim_disp;
            rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
            rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo&
            CurrentProcessInfo )
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = ( mNodesPressMin - 1 );i < mNodesPressMax;i++ )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( WATER_PRESSURE ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( AIR_PRESSURE ) );
        }

        for ( unsigned int i = ( mNodesDispMin - 1 );i < mNodesDispMax;i++ )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    //************************************************************************************
    //************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::GetValuesVector( Vector& values, int Step )
    {
        unsigned int dim_press = 2;//two pressure dofs two time nodes
        unsigned int dim_disp = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs two time nodes
        unsigned int MatSize =
            ( mNodesPressMax - mNodesPressMin + 1 ) * dim_press + ( mNodesDispMax -
                    mNodesDispMin + 1 ) * dim_disp;

        if ( values.size() != MatSize )
            values.resize( MatSize, false );

        unsigned int adddisp = ( mNodesPressMax - mNodesPressMin + 1 ) * dim_press;

        for ( unsigned int i = ( mNodesPressMin - 1 );i < mNodesPressMax;i++ )
        {
            int index = ( i - mNodesPressMin + 1 ) * dim_press;
            values( index ) =
                GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE, Step );
            values( index + 1 ) =
                GetGeometry()[i].GetSolutionStepValue( AIR_PRESSURE, Step );
        }

        for ( unsigned int i = ( mNodesDispMin - 1 );i < mNodesDispMax;i++ )
        {
            unsigned int index = adddisp + ( i - mNodesDispMin + 1 ) * dim_disp;

            values( index ) =
                GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
            values( index + 1 ) =
                GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );
            values( index + 2 ) =
                GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
        }

    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::AssembleTimeSpaceStiffnessFromDampSubMatrices
    ( MatrixType& rLeftHandSideMatrix,
      const Matrix& D_UU, const Matrix& D_UW, const Matrix& D_UA,
      const Matrix&  D_WU, const Matrix&  D_WW, const Matrix&  D_WA,
      const Matrix&  D_AU, const Matrix&  D_AW, const Matrix&  D_AA )
    {

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes_disp = mNodesDispMax - mNodesDispMin + 1;
        unsigned int number_of_nodes_press = mNodesPressMax - mNodesPressMin + 1;
        unsigned int dim_disp = dimension;
        unsigned int dim_press = 2;
        unsigned int index_time_prim;
        unsigned int index_time_sec;
        unsigned int index_space_prim;
        unsigned int index_space_sec;
        unsigned int addIndex_disp = number_of_nodes_press * dim_press;

        for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
        {
            for ( unsigned int i = 0; i < dim_disp; i++ )
            {
                index_space_prim = prim * dim_disp + i;

                index_time_prim = addIndex_disp + prim * dim_disp + i;

                for ( unsigned int sec = 0; sec < number_of_nodes_disp; sec++ )
                {
                    for ( unsigned int j = 0; j < dim_disp; j++ )
                    {
                        index_space_sec = sec * dim_disp + j;

                        index_time_sec = addIndex_disp + sec * dim_disp + j;

                        rLeftHandSideMatrix( index_time_prim,
                                             index_time_sec )
                        +=
                            ( -1 ) *
                            D_UU( index_space_prim, index_space_sec );
                    }
                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
        {
            for ( unsigned int i = 0; i < dim_disp; i++ )
            {
                index_space_prim = prim * dim_disp + i;

                index_time_prim = addIndex_disp + prim * dim_disp + i;

                for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
                {
                    index_space_sec = sec;

                    index_time_sec = sec * dim_press;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    +=
                        ( -1 ) *
                        D_UW( index_space_prim, index_space_sec );
                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
        {
            for ( unsigned int i = 0; i < dim_disp; i++ )
            {
                index_space_prim = prim * dim_disp + i;

                index_time_prim = addIndex_disp + prim * dim_disp + i;

                for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
                {
                    index_space_sec = sec;

                    index_time_sec = sec * dim_press + 1;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    +=
                        ( -1 ) *
                        D_UA( index_space_prim, index_space_sec );
                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {

            index_space_prim = prim;

            index_time_prim = prim * dim_press;

            for ( unsigned int sec = 0; sec < number_of_nodes_disp; sec++ )
            {
                for ( unsigned int j = 0; j < dim_disp; j++ )
                {
                    index_space_sec = sec * dim_disp + j;

                    index_time_sec = addIndex_disp + sec * dim_disp + j;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    +=
                        ( -1 ) *
                        D_WU( index_space_prim, index_space_sec );
                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_space_prim = prim;

            index_time_prim = prim * dim_press;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                +=
                    ( -1 ) *
                    D_WW( index_space_prim, index_space_sec );
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_space_prim = prim;

            index_time_prim = prim * dim_press;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press + 1;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                +=
                    ( -1 ) *
                    D_WA( index_space_prim, index_space_sec );
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {

            index_space_prim = prim;

            index_time_prim = prim * dim_press + 1;

            for ( unsigned int sec = 0; sec < number_of_nodes_disp; sec++ )
            {
                for ( unsigned int j = 0; j < dim_disp; j++ )
                {
                    index_space_sec = sec * dim_disp + j;

                    index_time_sec = addIndex_disp + sec * dim_disp + j;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    +=
                        ( -1 ) *
                        D_AU( index_space_prim, index_space_sec );
                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_space_prim = prim;

            index_time_prim = prim * dim_press + 1;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                +=
                    ( -1 ) *
                    D_AW( index_space_prim, index_space_sec );
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_space_prim = prim;

            index_time_prim = prim * dim_press + 1;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press + 1;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                +=
                    ( -1 ) *
                    D_AA( index_space_prim, index_space_sec );
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::AssembleTimeSpaceStiffnessFromStiffSubMatrices
    ( MatrixType& rLeftHandSideMatrix, const Matrix& K_UU, const Matrix& K_UW,
      const Matrix& K_UA, const Matrix& K_WU, const Matrix& K_WW,
      const Matrix& K_WA, const Matrix& K_AU, const Matrix& K_AW, const Matrix& K_AA
    )
    {
        KRATOS_TRY

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes_disp = mNodesDispMax - mNodesDispMin + 1;
        unsigned int number_of_nodes_press = mNodesPressMax - mNodesPressMin + 1;
        unsigned int dim_disp = dimension;
        unsigned int dim_press = 2;
        unsigned int index_time_prim;
        unsigned int index_time_sec;
        unsigned int index_space_prim;
        unsigned int index_space_sec;
        unsigned int addIndex_disp = number_of_nodes_press * dim_press;

        for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
        {
            for ( unsigned int i = 0; i < dim_disp; i++ )
            {
                index_space_prim = prim * dim_disp + i;

                index_time_prim = addIndex_disp + prim * dim_disp + i;

                for ( unsigned int sec = 0; sec < number_of_nodes_disp; sec++ )
                {
                    for ( unsigned int j = 0; j < dim_disp; j++ )
                    {
                        index_space_sec = sec * dim_disp + j;

                        index_time_sec = addIndex_disp + sec * dim_disp + j;

                        rLeftHandSideMatrix( index_time_prim, index_time_sec )
                        +=
                            ( -1 ) *
                            K_UU( index_space_prim, index_space_sec );
                    }
                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
        {
            for ( unsigned int i = 0; i < dim_disp; i++ )
            {
                index_space_prim = prim * dim_disp + i;

                index_time_prim = addIndex_disp + prim * dim_disp + i;

                for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
                {
                    index_space_sec = sec;

                    index_time_sec = sec * dim_press;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    +=
                        ( -1 ) *
                        K_UW( index_space_prim, index_space_sec );
                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
        {
            for ( unsigned int i = 0; i < dim_disp; i++ )
            {
                index_space_prim = prim * dim_disp + i;

                index_time_prim = addIndex_disp + prim * dim_disp + i;

                for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
                {
                    index_space_sec = sec;

                    index_time_sec = sec * dim_press + 1;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    +=
                        ( -1 ) *
                        K_UA( index_space_prim, index_space_sec );

                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {

            index_space_prim = prim;

            index_time_prim = prim * dim_press;

            for ( unsigned int sec = 0; sec < number_of_nodes_disp; sec++ )
            {
                for ( unsigned int j = 0; j < dim_disp; j++ )
                {
                    index_space_sec = sec * dim_disp + j;

                    index_time_sec = addIndex_disp + sec * dim_disp + j;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    +=
                        ( -1 ) *
                        K_WU( index_space_prim, index_space_sec );
                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_space_prim = prim;

            index_time_prim = prim * dim_press;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                +=
                    ( -1 ) *
                    K_WW( index_space_prim, index_space_sec );
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_space_prim = prim;

            index_time_prim = prim * dim_press;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press + 1;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                +=
                    ( -1 ) *
                    K_WA( index_space_prim, index_space_sec );

            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {

            index_space_prim = prim;

            index_time_prim = prim * dim_press + 1;

            for ( unsigned int sec = 0; sec < number_of_nodes_disp; sec++ )
            {
                for ( unsigned int j = 0; j < dim_disp; j++ )
                {
                    index_space_sec = sec * dim_disp + j;

                    index_time_sec = addIndex_disp + sec * dim_disp + j;

                    rLeftHandSideMatrix( index_time_prim, index_time_sec )
                    +=
                        ( -1 ) *
                        K_AU( index_space_prim, index_space_sec );

                }
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_space_prim = prim;

            index_time_prim = prim * dim_press + 1;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                +=
                    ( -1 ) *
                    K_AW( index_space_prim, index_space_sec );
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_space_prim = prim;

            index_time_prim = prim * dim_press + 1;

            for ( unsigned int sec = 0; sec < number_of_nodes_press; sec++ )
            {
                index_space_sec = sec;

                index_time_sec = sec * dim_press + 1;

                rLeftHandSideMatrix( index_time_prim, index_time_sec )
                +=
                    ( -1 ) *
                    K_AA( index_space_prim, index_space_sec );
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::AssembleTimeSpaceRHSFromSubVectors( VectorType&
            rRightHandSideVector, const Vector& R_U, const Vector& R_W, const Vector& R_A )
    {
        KRATOS_TRY

        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes_disp = ( mNodesDispMax - mNodesDispMin + 1 );
        unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );
        unsigned int dim_disp = dimension;
        unsigned int dim_press = 2;
        unsigned int index_time;
        unsigned int addIndex_disp = number_of_nodes_press * dim_press;

        for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
        {
            for ( unsigned int i = 0; i < dim_disp; i++ )
            {
                index_time = addIndex_disp + prim * dim_disp + i;

                rRightHandSideVector( index_time ) +=
                    R_U( prim * dim_disp + i );
            }
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_time = prim * dim_press;

            rRightHandSideVector( index_time ) +=
                R_W( prim );
        }

        for ( unsigned int prim = 0; prim < number_of_nodes_press; prim++ )
        {
            index_time = prim * dim_press + 1;

            rRightHandSideVector( index_time ) +=
                R_A( prim );
        }

        KRATOS_CATCH( "" )

    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE FORCEVECTORS DISPLACEMENT

    void UnsaturatedSoilsElement_3phase_SmallStrain::AddBodyForcesToRHSVectorU( Vector& R, Vector& N_DISP, double density, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        Vector gravity( dim );

        noalias( gravity ) = GetProperties()[GRAVITY];

        for ( unsigned int prim = 0; prim < ( mNodesDispMax - mNodesDispMin + 1 ); prim++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                R( prim*dim + i ) +=
                    N_DISP( prim ) * density * gravity( i ) *
                    detJ * Weight;
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::AddInternalForcesToRHSU
    ( Vector& R, const Matrix& B_Operator,
      Vector& StressVector, double Weight, double detJ )
    {
        KRATOS_TRY

        noalias( R ) -= detJ * Weight * prod( trans( B_Operator ), StressVector );

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //CALCULATE STIFFNESS MATRICES DISPLACEMENT

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixUU( Matrix& K,
            const Matrix& C, const Matrix& B_Operator, const Matrix& DN_DX_DISP,
            Vector& N_DISP, double density, double capillaryPressure,
            double airPressure, double Weight, double detJ )
    {
        KRATOS_TRY


        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int number_of_nodes_disp = mNodesDispMax - mNodesDispMin + 1;
        Vector Gravity( dim );
        noalias( Gravity ) = GetProperties()[GRAVITY];
        double density_soil = GetProperties()[DENSITY];
        double density_air = GetDensityAir( airPressure );
        double density_water = GetProperties()[DENSITY_WATER];
        double saturation = GetSaturation( capillaryPressure );
        double porosity_divu = GetDerivativeDPorosityDDivU( DN_DX_DISP );

        double DrhoDdivU =
            porosity_divu
            * ( -density_soil + ( 1 - saturation ) * density_air
                + saturation * density_water );

        for ( unsigned int prim = 0; prim < number_of_nodes_disp; prim++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                for ( unsigned int sec = 0; sec < number_of_nodes_disp; sec++ )
                {
                    for ( unsigned int j = 0; j < dim; j++ )
                    {
                        K( prim*dim + i, sec*dim + j ) += N_DISP( prim ) * DrhoDdivU * Gravity( i ) * DN_DX_DISP( sec, j ) * Weight * detJ;

//                                                 for(unsigned int alpha=0; alpha<6; alpha++)
//                                                         for(unsigned int beta=0; beta<6; beta++)
//                                                         K(prim*dim+i,sec*dim+j) += (-1)*(B_Operator[prim](alpha,i)
//                                                         *C(alpha,beta)*B_Operator[sec](beta,j))*detJ*Weight;
                    }
                }
            }
        }


        noalias( K ) -=

            prod( trans( B_Operator ), ( Weight * detJ ) * Matrix( prod( C, B_Operator ) ) );

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixUW( Matrix&
            Help_K_UW, Matrix& tanC_W, const Matrix& DN_DX_DISP, Vector& N_DISP, Vector& N_PRESS, double capillaryPressure, double airPressure,
            double Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

        Vector Gravity( dim );

        noalias( Gravity ) = GetProperties()[GRAVITY];

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double porosity = GetPorosity( DN_DX_DISP );

        double density_air = GetDensityAir( airPressure );

        double DrhoDp_w =
            porosity * ( density_air
                         - GetProperties()[DENSITY_WATER] ) * DSDpc;

        for ( unsigned int prim = 0; prim < displacement_size; prim++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                for ( unsigned int sec = 0; sec < pressure_size; sec++ )
                {
                    Help_K_UW( prim*dim + i, sec ) +=
                        N_DISP( prim ) * DrhoDp_w * Gravity( i ) * N_PRESS( sec )
                        * DetJ * Weight;

                    for ( unsigned int gamma = 0; gamma < dim; gamma++ )
                    {
                        Help_K_UW( prim*dim + i, sec ) +=
                            ( -1 ) * ( DN_DX_DISP( prim, gamma ) * tanC_W( i, gamma )
                                       * N_PRESS( sec ) )
                            * DetJ * Weight;
                    }

                }
            }
        }

    }

//************************************************************************************
//************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixUA( Matrix& Help_K_UA, Matrix& tanC_A, const Matrix& DN_DX_DISP, Vector& N_DISP, Vector& N_PRESS, double capillaryPressure, double airPressure, double Weight, double DetJ )
    {

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

        Vector Gravity( dim );

        noalias( Gravity ) = GetProperties()[GRAVITY];

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double porosity = GetPorosity( DN_DX_DISP );

        double density_air = GetDensityAir( airPressure );

        double saturation = GetSaturation( capillaryPressure );

        double DrhoDp_a =
            porosity * (( -density_air
                          + GetProperties()[DENSITY_WATER] ) * DSDpc + ( 1 - saturation ) * GetProperties()[BULK_AIR] );

        for ( unsigned int prim = 0; prim < displacement_size; prim++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                for ( unsigned int sec = 0; sec < pressure_size; sec++ )
                {
                    Help_K_UA( prim*dim + i, sec ) +=
                        N_DISP( prim ) * DrhoDp_a * Gravity( i ) * N_PRESS( sec )
                        * DetJ * Weight;

                    for ( unsigned int gamma = 0; gamma < dim; gamma++ )
                    {
                        Help_K_UA( prim*dim + i, sec ) +=
                            ( -1 ) * ( DN_DX_DISP( prim, gamma ) * tanC_A( i, gamma ) * N_PRESS( sec ) )
                            * DetJ * Weight;
                    }
                }
            }
        }

    }

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

    //CALCULATE FORCEVECTORS WATER

    void UnsaturatedSoilsElement_3phase_SmallStrain::AddInternalForcesToRHSW( Vector& Help_R_W, const
            Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS,
            double capillaryPressure, double Weight, double  DetJ )
    {
        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double porosity = GetPorosity( DN_DX_DISP );

        double DS_Dpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double saturation = GetSaturation( capillaryPressure );

        double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

        double div_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

        Vector flow_water( dim );
        noalias( flow_water ) = GetFlowWater( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            Help_R_W( prim ) +=
                N_PRESS( prim ) * porosity * DS_Dpc * Dpc_Dt
                * Weight * DetJ * GetProperties()[SCALE];
            Help_R_W( prim ) +=
                N_PRESS( prim ) * saturation * div_Dt
                * Weight * DetJ * GetProperties()[SCALE];

            for ( unsigned int gamma = 0; gamma < dim; gamma++ )
            {
                Help_R_W( prim ) +=
                    ( -1 ) * ( DN_DX_PRESS( prim, gamma ) * flow_water( gamma ) )
                    * Weight * DetJ * GetProperties()[SCALE];
            }
        }
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************


    //CALCULATE STIFFNESS MATRICES WATER

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixWU
    ( Matrix& Help_K_WU, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

        double DnDdivU = GetDerivativeDPorosityDDivU( DN_DX_DISP );

        Vector flow_water( dim );
        noalias( flow_water ) = GetFlowWater( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < displacement_size; sec++ )
            {
                for ( unsigned int j = 0; j < dim; j++ )
                {
                    Help_K_WU( prim, sec*dim + j ) +=
                        N_PRESS( prim )
                        * DnDdivU * DSDpc * Dpc_Dt *
                        DN_DX_DISP( sec, j ) * Weight * DetJ * GetProperties()[SCALE];
                }
            }
        }
    }

    //************************************************************************************
    //************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixWW
    ( Matrix& Help_K_WW, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS,
      Vector& N_PRESS, double capillaryPressure, double
      Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double porosity = GetPorosity( DN_DX_DISP );

        double D2S_Dpc2 = GetSecondDerivativeD2SaturationDpc2( capillaryPressure );

        double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

        Vector Dflow_waterDpw( dim );

        noalias( Dflow_waterDpw ) = GetDerivativeDWaterFlowDpw( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

        double Dflow_waterDgradpw = GetDerivativeDWaterFlowDGradpw( DN_DX_DISP, capillaryPressure );

        double Ddiv_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_K_WW( prim, sec ) +=
                    ( -1 ) * N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_K_WW( prim, sec ) +=
                    ( -1 ) * N_PRESS( prim ) * DSDpc * Ddiv_Dt * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                for ( unsigned int gamma = 0; gamma < dim; gamma++ )
                {
                    Help_K_WW( prim, sec ) +=
                        ( -1 ) * DN_DX_PRESS( prim, gamma ) * Dflow_waterDpw( gamma )
                        * N_PRESS( sec )
                        * Weight * DetJ * GetProperties()[SCALE];
                    Help_K_WW( prim, sec ) +=
                        ( -1 ) * DN_DX_PRESS( prim, gamma ) * Dflow_waterDgradpw
                        * DN_DX_PRESS( sec, gamma ) * Weight * DetJ * GetProperties()[SCALE];
                }

            }
        }

    }

//************************************************************************************
//************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixWA( Matrix& Help_K_WA, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double porosity = GetPorosity( DN_DX_DISP );

        double D2S_Dpc2 = GetSecondDerivativeD2SaturationDpc2( capillaryPressure );

        double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

        Vector Dflow_waterDpa( dim );

        noalias( Dflow_waterDpa ) = GetDerivativeDWaterFlowDpa( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

        double Ddiv_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_K_WA( prim, sec ) +=
                    N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_K_WA( prim, sec ) +=
                    N_PRESS( prim ) * DSDpc * Ddiv_Dt * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                for ( unsigned int gamma = 0; gamma < dim; gamma++ )
                {
                    Help_K_WA( prim, sec ) +=
                        ( -1 ) * DN_DX_PRESS( prim, gamma ) * Dflow_waterDpa( gamma )
                        * N_PRESS( sec )
                        * Weight * DetJ * GetProperties()[SCALE];
                }

            }
        }
    }

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

    //CALCULATE DAMPING MATRICES WATER
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateDampingMatrixWU
    ( Matrix& Help_D_WU, const Matrix&
      DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

        double saturation = GetSaturation( capillaryPressure );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < displacement_size; sec++ )
            {
                for ( unsigned int j = 0; j < dim; j++ )
                {
                    Help_D_WU( prim, sec*dim + j ) +=
                        N_PRESS( prim ) * saturation * DN_DX_DISP( sec, j )
                        * Weight * DetJ * GetProperties()[SCALE];
                }
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateDampingMatrixWW( Matrix& Help_D_WW, const Matrix&
            DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
    {
        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
        double porosity = GetPorosity( DN_DX_DISP );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_D_WW( prim, sec ) +=
                    ( -1 ) * N_PRESS( prim ) * porosity * DSDpc * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];
            }
        }
    }

//************************************************************************************
//************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateDampingMatrixWA( Matrix& Help_D_WA, const Matrix& DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
    {
        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
        double porosity = GetPorosity( DN_DX_DISP );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_D_WA( prim, sec ) +=
                    N_PRESS( prim ) * porosity * DSDpc * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];
            }
        }
    }

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::AddInternalForcesToRHSA( Vector& Help_R_A, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double airPressure, double airPressure_Dt, double Weight, double  DetJ )
    {
        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double porosity = GetPorosity( DN_DX_DISP );

        double DS_Dpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double saturation = GetSaturation( capillaryPressure );

        double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

        double density_air = GetDensityAir( airPressure );

        double div_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

        Vector grad_air( dim );
        noalias( grad_air ) = GetGradAirPressure( DN_DX_PRESS );

        Vector flow_air( dim );
        noalias( flow_air ) = GetFlowAir( DN_DX_PRESS, DN_DX_DISP, airPressure, capillaryPressure );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            Help_R_A( prim ) +=
                ( -1 ) * N_PRESS( prim ) * porosity * DS_Dpc * Dpc_Dt
                * Weight * DetJ * GetProperties()[SCALE];
            Help_R_A( prim ) +=
                N_PRESS( prim ) * ( 1 - saturation ) * div_Dt
                * Weight * DetJ * GetProperties()[SCALE];
            Help_R_A( prim ) +=
                N_PRESS( prim ) * porosity * ( 1 - saturation ) / density_air *
                GetProperties()[BULK_AIR] * airPressure_Dt
                * Weight * DetJ * GetProperties()[SCALE];

            for ( unsigned int gamma = 0; gamma < dim; gamma++ )
            {
                Help_R_A( prim ) +=
                    N_PRESS( prim ) * 1.0 / density_air *
                    GetProperties()[BULK_AIR] * grad_air( gamma ) * flow_air( gamma )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_R_A( prim ) +=
                    ( -1 ) * ( DN_DX_PRESS( prim, gamma ) * flow_air( gamma ) )
                    * Weight * DetJ * GetProperties()[SCALE];
            }
        }
    }

    //CALCULATE STIFFNESS MATRICES AIR

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixAU( Matrix& Help_K_AU, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS , double capillaryPressure, double airPressure, double airPressure_Dt, double Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double saturation = GetSaturation( capillaryPressure );

        double density_air = GetDensityAir( airPressure );

        double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

        double DnDdivU = GetDerivativeDPorosityDDivU( DN_DX_DISP );

        Vector grad_air( dim );
        noalias( grad_air ) = GetGradAirPressure( DN_DX_PRESS );

        Vector flow_air( dim );
        noalias( flow_air ) = GetFlowAir( DN_DX_PRESS, DN_DX_DISP, airPressure, capillaryPressure );

//            vector<unsigned int> help;
//            help.resize(2);

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < displacement_size; sec++ )
            {
                for ( unsigned int j = 0; j < dim; j++ )
                {
                    Help_K_AU( prim, sec*dim + j ) +=
                        N_PRESS( prim )
                        * DnDdivU * ( 1 - saturation ) / density_air * GetProperties()[BULK_AIR]
                        * airPressure_Dt * DN_DX_DISP( sec, j )
                        * Weight * DetJ * GetProperties()[SCALE];
                    Help_K_AU( prim, sec*dim + j ) +=
                        ( -1 ) * N_PRESS( prim )
                        * DnDdivU * DSDpc * Dpc_Dt *
                        DN_DX_DISP( sec, j )
                        * Weight * DetJ * GetProperties()[SCALE];
                }
            }
        }
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixAW( Matrix& Help_K_AW, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double airPressure, double airPressure_Dt,
            double Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double porosity = GetPorosity( DN_DX_DISP );

        double D2S_Dpc2 = GetSecondDerivativeD2SaturationDpc2( capillaryPressure );

        double density_air = GetDensityAir( airPressure );

        double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

        Vector Dflow_airDpw( dim );

        noalias( Dflow_airDpw ) = GetDerivativeDAirFlowDpw( DN_DX_PRESS, DN_DX_DISP, airPressure, capillaryPressure );

        Vector grad_air( dim );
        noalias( grad_air ) = GetGradAirPressure( DN_DX_PRESS );

        double Ddiv_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_K_AW( prim, sec ) +=
                    N_PRESS( prim ) * DSDpc * porosity / density_air * GetProperties()[BULK_AIR] * airPressure_Dt
                    * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_K_AW( prim, sec ) +=
                    N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_K_AW( prim, sec ) +=
                    N_PRESS( prim ) * DSDpc * Ddiv_Dt * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                for ( unsigned int gamma = 0; gamma < dim; gamma++ )
                {
                    Help_K_AW( prim, sec ) +=
                        ( -1 ) * DN_DX_PRESS( prim, gamma ) * Dflow_airDpw( gamma )
                        * N_PRESS( sec )
                        * Weight * DetJ * GetProperties()[SCALE];

                    Help_K_AW( prim, sec ) +=
                        N_PRESS( prim ) * 1 / density_air * GetProperties()[BULK_AIR]
                        * grad_air( gamma ) * Dflow_airDpw( gamma )
                        * N_PRESS( sec )
                        * Weight * DetJ * GetProperties()[SCALE];
                }

            }
        }
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStiffnesMatrixAA( Matrix& Help_K_AA, const Matrix& DN_DX_DISP, const Matrix& DN_DX_PRESS, Vector& N_PRESS, double capillaryPressure, double airPressure, double airPressure_Dt, double Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        double porosity = GetPorosity( DN_DX_DISP );

        double D2S_Dpc2 = GetSecondDerivativeD2SaturationDpc2( capillaryPressure );

        double Dpc_Dt = GetDerivativeDCapillaryPressureDt( N_PRESS );

        double saturation = GetSaturation( capillaryPressure );

        Vector Dflow_airDpa( dim );

        noalias( Dflow_airDpa ) = GetDerivativeDAirFlowDpa( DN_DX_PRESS, DN_DX_DISP, airPressure, capillaryPressure );

        Vector GradPa( dim );

        noalias( GradPa ) = GetGradAirPressure( DN_DX_PRESS );

        Vector flow_air( dim );

        noalias( flow_air ) = GetFlowAir( DN_DX_PRESS, DN_DX_DISP, airPressure, capillaryPressure );

        double Dflow_airDgradpa = GetDerivativeDAirFlowDGradpa( DN_DX_PRESS, airPressure, capillaryPressure );

        double Ddiv_Dt = GetDerivativeDDivUDt( DN_DX_DISP );

        double density_air = GetDensityAir( airPressure );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_K_AA( prim, sec ) +=
                    ( -1 ) * N_PRESS( prim ) * porosity * D2S_Dpc2 * Dpc_Dt * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_K_AA( prim, sec ) +=
                    ( -1 ) * N_PRESS( prim ) * DSDpc * Ddiv_Dt * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_K_AA( prim, sec ) +=
                    ( -1 ) * N_PRESS( prim ) *
                    DSDpc * porosity / density_air * GetProperties()[BULK_AIR] * airPressure_Dt
                    * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_K_AA( prim, sec ) +=
                    ( -1 ) * N_PRESS( prim ) *
                    porosity * ( 1 - saturation ) / ( density_air * density_air ) *
                    GetProperties()[BULK_AIR] * GetProperties()[BULK_AIR] * airPressure_Dt
                    * N_PRESS( sec ) * Weight * DetJ * GetProperties()[SCALE];

                for ( unsigned int gamma = 0; gamma < dim; gamma++ )
                {
//neglection of convective term
                    Help_K_AA( prim, sec ) +=
                        ( -1 ) * N_PRESS( prim ) *
                        1 / ( density_air * density_air ) * ( GetProperties()[BULK_AIR] * GetProperties()[BULK_AIR] ) * GradPa( gamma )
                        * flow_air( gamma )
                        * N_PRESS( sec )
                        * Weight * DetJ * GetProperties()[SCALE];

//neglection of convective term
                    Help_K_AA( prim, sec ) +=
                        N_PRESS( prim ) *
                        1 / density_air * GetProperties()[BULK_AIR] * flow_air( gamma )
                        * DN_DX_PRESS( sec, gamma )
                        * Weight * DetJ * GetProperties()[SCALE];

//neglection of convective term
                    Help_K_AA( prim, sec ) +=
                        N_PRESS( prim ) *
                        1 / density_air * GetProperties()[BULK_AIR] * GradPa( gamma ) * Dflow_airDpa( gamma )
                        * N_PRESS( sec )
                        * Weight * DetJ * GetProperties()[SCALE];

//neglection of convective term
                    Help_K_AA( prim, sec ) +=
                        N_PRESS( prim ) *
                        1 / density_air * GetProperties()[BULK_AIR] * GradPa( gamma ) * Dflow_airDgradpa
                        * DN_DX_PRESS( sec, gamma )
                        * Weight * DetJ * GetProperties()[SCALE];

                    Help_K_AA( prim, sec ) +=
                        ( -1 ) * DN_DX_PRESS( prim, gamma ) * Dflow_airDpa( gamma )
                        * N_PRESS( sec )
                        * Weight * DetJ * GetProperties()[SCALE];

                    Help_K_AA( prim, sec ) +=
                        ( -1 ) * DN_DX_PRESS( prim, gamma ) * Dflow_airDgradpa
                        * DN_DX_PRESS( sec, gamma )
                        * Weight * DetJ * GetProperties()[SCALE];
                }

            }
        }

    }

//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateDampingMatrixAU( Matrix& Help_D_AU, const Matrix& DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        unsigned int displacement_size = mNodesDispMax - mNodesDispMin + 1;

        double saturation = GetSaturation( capillaryPressure );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < displacement_size; sec++ )
            {
                for ( unsigned int j = 0; j < dim; j++ )
                {
                    Help_D_AU( prim, sec*dim + j ) +=
                        N_PRESS( prim ) * ( 1 - saturation ) * DN_DX_DISP( sec, j )
                        * Weight * DetJ * GetProperties()[SCALE];
                }
            }
        }
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateDampingMatrixAW( Matrix& Help_D_AW, const Matrix& DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double Weight, double DetJ )
    {
        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
        double porosity = GetPorosity( DN_DX_DISP );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_D_AW( prim, sec ) +=
                    N_PRESS( prim ) * porosity * DSDpc * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];
            }
        }
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateDampingMatrixAA( Matrix& Help_D_AA, const Matrix& DN_DX_DISP, Vector& N_PRESS, double capillaryPressure, double airPressure, double Weight, double DetJ )
    {
        unsigned int pressure_size = mNodesPressMax - mNodesPressMin + 1;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
        double porosity = GetPorosity( DN_DX_DISP );

        double saturation = GetSaturation( capillaryPressure );

        double density_air = GetDensityAir( airPressure );

        for ( unsigned int prim = 0; prim < pressure_size; prim++ )
        {
            for ( unsigned int sec = 0; sec < pressure_size; sec++ )
            {
                Help_D_AA( prim, sec ) +=
                    ( -1 ) * N_PRESS( prim ) * porosity * DSDpc * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];

                Help_D_AA( prim, sec ) +=
                    N_PRESS( prim ) * porosity * ( 1 - saturation ) / density_air
                    * GetProperties()[BULK_AIR] * N_PRESS( sec )
                    * Weight * DetJ * GetProperties()[SCALE];
            }
        }
    }

    //************************************************************************************
    //************************************************************************************

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetGradWaterPressure( const Matrix& DN_DX_PRESS )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        Vector result( dim );
        noalias( result ) = ZeroVector( dim );

        double presW_alpha;

        for ( unsigned int i = mNodesPressMin - 1 ; i < mNodesPressMax ;i++ )
        {

//                    //nodal Displacements
            presW_alpha = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

            //thus strain in vector format

            for ( unsigned int k = 0; k < 3; k++ )
            {
                result( k ) +=
                    presW_alpha
                    * DN_DX_PRESS(( i - mNodesPressMin + 1 ), k );
            }
        }

        return result;
    }

    //************************************************************************************
    //************************************************************************************

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetGradAirPressure( const Matrix& DN_DX_PRESS )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        Vector result( dim );
        noalias( result ) = ZeroVector( dim );

        double presA_alpha;

        for ( unsigned int i = mNodesPressMin - 1 ; i < mNodesPressMax ;i++ )
        {
//                    //nodal Displacements
            presA_alpha = GetGeometry()[i].GetSolutionStepValue( AIR_PRESSURE );

            //thus strain in vector format

            for ( unsigned int k = 0; k < 3; k++ )
            {
                result( k ) +=
                    presA_alpha
                    * DN_DX_PRESS(( i - mNodesPressMin + 1 ), k );
            }
        }

        return result;
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::GetPressures( const Vector& N_PRESS,
            double& capillaryPressure, double& waterPressure, double& airPressure )
    {
        capillaryPressure = 0.0;

        waterPressure = 0.0;

        airPressure = 0.0;

        double presW_alpha;

        double presA_alpha;


        //Calculating Strain on space points

        for ( unsigned int i = ( mNodesPressMin - 1 ) ; i < mNodesPressMax ;i++ )
        {

            presW_alpha = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE );

            presA_alpha = GetGeometry()[i].GetSolutionStepValue( AIR_PRESSURE );

            //thus strain in vector format
            capillaryPressure +=
                ( presA_alpha - presW_alpha )
                * N_PRESS( i - mNodesPressMin + 1 );

            waterPressure +=
                presW_alpha
                * N_PRESS( i - mNodesPressMin + 1 );

            airPressure +=
                presA_alpha
                * N_PRESS( i - mNodesPressMin + 1 );
        }
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDPressuresDt( const Vector& N_PRESS, double& capillaryPressure_Dt, double& waterPressure_Dt, double& airPressure_Dt )
    {
        capillaryPressure_Dt = 0.0;

        waterPressure_Dt = 0.0;

        airPressure_Dt = 0.0;

        double presW_alpha_Dt;

        double presA_alpha_Dt;
        //Calculating Strain on space points

        for ( unsigned int i = mNodesPressMin - 1 ; i < mNodesPressMax ;i++ )
        {
            presW_alpha_Dt =
                GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT );

            presA_alpha_Dt =
                GetGeometry()[i].GetSolutionStepValue( AIR_PRESSURE_DT );

            capillaryPressure_Dt +=
                ( presA_alpha_Dt - presW_alpha_Dt )
                * N_PRESS( i - mNodesPressMin + 1 );

            waterPressure_Dt +=
                presW_alpha_Dt
                * N_PRESS( i - mNodesPressMin + 1 );

            airPressure_Dt +=
                presA_alpha_Dt
                * N_PRESS( i - mNodesPressMin + 1 );
        }
    }

    //************************************************************************************
    //************************************************************************************

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDCapillaryPressureDt( const Vector& N_PRESS )
    {
        double capillaryPressure_Dt = 0.0;

        double presW_alpha_Dt;

        double presA_alpha_Dt;
        //Calculating Strain on space points

        for ( unsigned int i = mNodesPressMin - 1 ; i < mNodesPressMax ;i++ )
        {
            presW_alpha_Dt = GetGeometry()[i].GetSolutionStepValue( WATER_PRESSURE_DT );

            presA_alpha_Dt = GetGeometry()[i].GetSolutionStepValue( AIR_PRESSURE_DT );


            capillaryPressure_Dt +=
                ( presA_alpha_Dt - presW_alpha_Dt )
                * N_PRESS( i - mNodesPressMin + 1 );

        }

        return capillaryPressure_Dt;
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //POROSITY AND ITS DERIVATIVES

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetPorosity( const Matrix& DN_DX_DISP )
    {
        double initialPorosity = GetProperties()[POROSITY];

        double porosity = initialPorosity;
//             double div= GetDivU(DN_DX_DISP);
//
//             double porosity= 1-(1-initialPorosity)*exp(-div);
//
//             if(porosity < 0)
//             {
//     std::cout<<"porosity PROBLEM"<<std::endl;
//                 KRATOS_ERROR(std::logic_error, "Porosity is less than zero" , *this);
//             }
//             if(porosity > 1)
//             {
//     std::cout<<"porosity PROBLEM"<<std::endl;
//              KRATOS_ERROR(std::logic_error, "Porosity is bigger than one" , *this);
//             }

        return porosity;
    }

    //************************************************************************************
    //************************************************************************************

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDPorosityDDivU( const Matrix& DN_DX_DISP )
    {
//             double initialPorosity= GetProperties()[POROSITY];
//             double div= GetDivU(DN_DX_DISP);

        double porosity_divu = 0.0;
//             double porosity_divu= (1-initialPorosity)*exp(-div);

        return porosity_divu;
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //DENSITY AIR AND ITS DERIVATIVES

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetDensityAir( double airPressure )
    {

        double result = GetProperties()[DENSITY_AIR] + GetProperties()[BULK_AIR] * airPressure;
        return result;
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetGravity()
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        Vector gravity( dim );
        noalias( gravity ) = GetProperties()[GRAVITY];
        return gravity;
    }

    //************************************************************************************
    //************************************************************************************
    double UnsaturatedSoilsElement_3phase_SmallStrain::GetDivU( const Matrix& DN_DX_DISP )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double div = 0.0;

        Vector u_alpha( dim );

        for ( unsigned int i = mNodesDispMin - 1 ; i < mNodesDispMax ;i++ )
        {
//                     //nodal Displacements
            noalias( u_alpha ) =
                GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT );
            //thus strain in vector format

            for ( unsigned int k = 0; k < dim; k++ )
            {
                div += ( u_alpha( k ) )
                       * DN_DX_DISP( i - mNodesDispMin + 1, k );
            }
        }

        return div;
    }

    //************************************************************************************
    //************************************************************************************
    double UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDDivUDt( const Matrix& DN_DX_DISP )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double div = 0.0;

        Vector u_alpha_Dt( dim );

        for ( unsigned int i = mNodesDispMin - 1 ; i < mNodesDispMax ;i++ )
        {
//                                //nodal Displacements
            noalias( u_alpha_Dt ) =
                GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_DT );

            //thus strain in vector format

            for ( unsigned int k = 0; k < dim; k++ )
            {
                div += u_alpha_Dt( k ) * DN_DX_DISP( i - mNodesDispMin + 1, k );
            }
        }

        return div;
    }

    //AVERAGED DENSITY
    //************************************************************************************
    //************************************************************************************

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetAveragedDensity( double capillaryPressure,
            double airPressure, double porosity )
    {
        double result = 0.0;

        double density_soil = GetProperties()[DENSITY];
        double density_air = GetDensityAir( airPressure );
        double density_water = GetProperties()[DENSITY_WATER];
        double saturation = GetSaturation( capillaryPressure );

        result = ( 1.0 - porosity ) * density_soil +
                 porosity * ( saturation * density_water + ( 1.0 - saturation ) * density_air );

        return result;
    }


    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //SATURATION AND ITS DERIVATIVES

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetSaturation( double capillaryPressure )
    {
        double airEntryPressure = GetProperties()[AIR_ENTRY_VALUE];

        if ( airEntryPressure <= 0.0 )
            airEntryPressure = 1.0;

        double b = GetProperties()[FIRST_SATURATION_PARAM];

        double c = GetProperties()[SECOND_SATURATION_PARAM];

        double saturation = 0.0;

//
        if ( capillaryPressure < 0.0 )
            capillaryPressure = 0.0;

        saturation = pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c ) );

// For Liakopolous Benchmark
// saturation =  1.0-1.9722*1e-11*pow(capillaryPressure,2.4279);


        return saturation;
    }

//************************************************************************************
//************************************************************************************

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDSaturationDpc( double capillaryPressure )
    {
        double airEntryPressure = GetProperties()[AIR_ENTRY_VALUE];

        if ( airEntryPressure <= 0 )
            airEntryPressure = 1.0;

        double b = GetProperties()[FIRST_SATURATION_PARAM];

        double c = GetProperties()[SECOND_SATURATION_PARAM];

        double result = 0.0;

//
        if ( capillaryPressure < 0.0 )
            capillaryPressure = 0.0;

        result = ( -c ) * pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 1.0 ) ) * b *
                 pow(( capillaryPressure / airEntryPressure ), ( b - 1.0 ) ) * 1.0 / airEntryPressure;

// For Liakopolous Benchmark
// result =  -1.9722*2.4279*1e-11*pow(capillaryPressure,1.4279);

        return result;
    }

    //************************************************************************************
    //************************************************************************************

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetSecondDerivativeD2SaturationDpc2( double capillaryPressure )
    {
        double airEntryPressure = GetProperties()[AIR_ENTRY_VALUE];

        if ( airEntryPressure <= 0 )
            airEntryPressure = 1.0;

        double b = GetProperties()[FIRST_SATURATION_PARAM];

        double c = GetProperties()[SECOND_SATURATION_PARAM];

        double result = 0.0;

        if ( capillaryPressure < 0.0 )
            capillaryPressure = 0.0;

        result = ( -c ) * b / airEntryPressure * (
                     ( -c - 1.0 ) * pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 2.0 ) )
                     * b / airEntryPressure * pow(( capillaryPressure / airEntryPressure ), ( 2.0 * ( b - 1.0 ) ) )
                     + pow(( 1.0 + pow(( capillaryPressure / airEntryPressure ), b ) ), ( -c - 1.0 ) ) * ( b - 1.0 )
                     * pow(( capillaryPressure / airEntryPressure ), ( b - 2.0 ) ) * 1.0 / airEntryPressure );

// For Liakopolous Benschmark
// result =  -1.9722*2.4279*1.4279*1e-11*pow(capillaryPressure,0.4279);

        return result;
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //WATER FLOW AND ITS DERIVATIVES

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetFlowWater( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double capillaryPressure )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //Calculation of relative Permeability after Mualem
        double relPerm = GetSaturation( capillaryPressure );

        if ( relPerm <= 0.01 ) relPerm = 0.01;

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

        Vector gravity( dim );

        noalias( gravity ) = GetProperties()[GRAVITY];

        Vector result( dim );

        noalias( result ) = ZeroVector( dim );

        Vector grad_water( dim );

        noalias( grad_water ) = GetGradWaterPressure( DN_DX_PRESS );

        for ( unsigned int i = 0; i < dim; i++ )
        {
            result( i ) = -relPerm * GetProperties()[PERMEABILITY_WATER] /
                          ( GetProperties()[DENSITY_WATER] * 9.81 )
                          * ( grad_water( i ) - GetProperties()[DENSITY_WATER]
                              * gravity( i ) );
        }

        return result;
    }

//************************************************************************************
//************************************************************************************

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDWaterFlowDpw( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double capillaryPressure )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
        //Calculation of Derivative relative Permeability after Mualem
        double relPerm = GetSaturation( capillaryPressure );
        double relPerm_pw = ( -1.0 ) * DSDpc;

        if ( relPerm <= 0.01 )
        {
            relPerm = 0.01;
            relPerm_pw = 0.0;
        }

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm_pw= -2.207*1.0121*pow((1-saturation),0.0121)*(-DSDpc)*(-1);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

        Vector result( dim );

        noalias( result ) = ZeroVector( dim );

        Vector gravity( dim );

        noalias( gravity ) = GetProperties()[GRAVITY];

        Vector grad_water( dim );

        noalias( grad_water ) = GetGradWaterPressure( DN_DX_PRESS );

        for ( unsigned int i = 0; i < dim; i++ )
        {
            result( i ) = -relPerm_pw * GetProperties()[PERMEABILITY_WATER] /
                          ( GetProperties()[DENSITY_WATER] * 9.81 )
                          * ( grad_water( i ) - GetProperties()[DENSITY_WATER]
                              * gravity( i ) );
        }

        return result;
    }

    //************************************************************************************
    //************************************************************************************

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDWaterFlowDpa( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double capillaryPressure )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
        //Calculation of Derivative relative Permeability after Mualem
        double relPerm = GetSaturation( capillaryPressure );
        double relPerm_pa = DSDpc;

        if ( relPerm <= 0.01 )
        {
            relPerm = 0.01;
            relPerm_pa = 0.0;
        }

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm_pa= -2.207*1.0121*pow((1-saturation),0.0121)*(-DSDpc);
// // double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);


        Vector result( dim );

        noalias( result ) = ZeroVector( dim );

        Vector gravity( dim );

        noalias( gravity ) = GetProperties()[GRAVITY];

        Vector grad_water( dim );

        noalias( grad_water ) = GetGradWaterPressure( DN_DX_PRESS );

        for ( unsigned int i = 0; i < dim; i++ )
        {
            result( i ) = -relPerm_pa * GetProperties()[PERMEABILITY_WATER] /
                          ( GetProperties()[DENSITY_WATER] * 9.81 )
                          * ( grad_water( i ) - GetProperties()[DENSITY_WATER]
                              * gravity( i ) );
        }

        return result;
    }

    //************************************************************************************
    //************************************************************************************

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDWaterFlowDGradpw( const Matrix& DN_DX_DISP, double capillaryPressure )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //Calculation of Derivative relative Permeability after Mualem
        double relPerm = GetSaturation( capillaryPressure );

        if ( relPerm <= 0.01 )
            relPerm = 0.01;

// For Liakopolous Benschmark
// double saturation= GetSaturation(capillaryPressure);
// double relPerm= 1.0- 2.207*pow((1-saturation),1.0121);

        double result( dim );

        result = -relPerm * GetProperties()[PERMEABILITY_WATER] / ( GetProperties()[DENSITY_WATER] * 9.81 );

        return result;
    }

    //************************************************************************************
    //************************************************************************************
    //************************************************************************************
    //************************************************************************************

    //WATER FLOW AND ITS DERIVATIVES

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetFlowAir( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double airPressure, double capillaryPressure )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double saturation = GetSaturation( capillaryPressure );

// For Liakopolous Benschmark
// double relSat= (saturation-0.2)/0.8;
// double relPerm= 0.0001+pow((1.0-relSat),2)*(1-pow(relSat,5.0/3.0));

        double relPerm =  1.0 - saturation;

        if ( relPerm <= 0.01 ) relPerm = 0.01;

        Vector gravity( dim );

        noalias( gravity ) = GetProperties()[GRAVITY];

        Vector result( dim );

        noalias( result ) = ZeroVector( dim );

        double airDensity = GetDensityAir( airPressure );

        Vector grad_air( dim );

        noalias( grad_air ) = GetGradAirPressure( DN_DX_PRESS );

        for ( unsigned int i = 0; i < dim; i++ )
        {
            result( i ) = -relPerm * GetProperties()[PERMEABILITY_AIR] / ( airDensity * 9.81 ) *
                          ( grad_air( i ) - airDensity * gravity( i ) );
        }

        return result;
    }

//************************************************************************************
//************************************************************************************

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDAirFlowDpa( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double airPressure, double capillaryPressure )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double saturation = GetSaturation( capillaryPressure );

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
        //Calculation of Derivative relative Permeability after Mualem
        double relPerm_pa = -DSDpc;
        double relPerm = 1.0 - saturation;

        if ( relPerm <= 0.01 )
        {
            relPerm = 0.01;
            relPerm_pa = 0.0;
        }

// For Liakopolous Benschmark
// double relSat= (saturation-0.2)/0.8;
// double relPerm_pa= 2*(1.0-relSat)*(-1)/0.8*DSDpc*(1-pow(relSat,5.0/3.0))
//    +pow((1.0-relSat),2)*(-1)*5.0/3.0*pow(relSat,2.0/3.0)/0.8*DSDpc;
// double relPerm= 0.0001+ pow((1.0-relSat),2)*(1-pow(relSat,5.0/3.0));

        Vector gravity( dim );

        noalias( gravity ) = GetProperties()[GRAVITY];

        Vector result( dim );

        noalias( result ) = ZeroVector( dim );

        Vector grad_air( dim );

        noalias( grad_air ) = GetGradAirPressure( DN_DX_PRESS );

        double airDensity = GetDensityAir( airPressure );

        for ( unsigned int i = 0; i < dim; i++ )
        {
            result( i ) = -relPerm * GetProperties()[PERMEABILITY_AIR] / ( airDensity * 9.81 )
                          * ( -GetProperties()[BULK_AIR] * gravity( i ) )
                          - relPerm_pa * GetProperties()[PERMEABILITY_AIR] / ( airDensity * 9.81 )
                          * ( grad_air( i ) - airDensity * gravity( i ) )
                          + relPerm * GetProperties()[PERMEABILITY_AIR] /
                          ( airDensity * airDensity * 9.81 ) * GetProperties()[BULK_AIR]
                          * ( grad_air( i ) - airDensity * gravity( i ) );
        }

        return result;
    }

//************************************************************************************
//************************************************************************************

    Vector UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDAirFlowDpw( const Matrix& DN_DX_PRESS, const Matrix& DN_DX_DISP, double airPressure, double capillaryPressure )
    {
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double saturation = GetSaturation( capillaryPressure );
        double relPerm = 1.0 - saturation;

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );
// For Liakopolous Benschmark
// double relSat= (saturation-0.2)/0.8;
// double relPerm_pw= 2*(1.0-relSat)/0.8*DSDpc*(1-pow(relSat,5.0/3.0))
//    +pow((1.0-relSat),2)*5.0/3.0*pow(relSat,2.0/3.0)/0.8*DSDpc;
// double relPerm= 0.0001+ pow((1.0-relSat),2)*(1-pow(relSat,5.0/3.0));
//   double relPerm= 1.0-saturation;
        double relPerm_pw = DSDpc;

        if ( relPerm <= 0.01 )
        {
            relPerm = 0.0;
            relPerm_pw = 0.0;
        }

        Vector result( dim );

        noalias( result ) = ZeroVector( dim );

        Vector grad_air( dim );
        noalias( grad_air ) = GetGradAirPressure( DN_DX_PRESS );

        double airDensity = GetDensityAir( airPressure );

        Vector gravity( dim );
        noalias( gravity ) = GetProperties()[GRAVITY];

        for ( unsigned int i = 0; i < dim; i++ )
        {
            result( i ) = -relPerm_pw * GetProperties()[PERMEABILITY_AIR] / ( airDensity * 9.81 )
                          * ( grad_air( i ) - airDensity * gravity( i ) );
        }

        return result;
    }

//************************************************************************************
//************************************************************************************

    double UnsaturatedSoilsElement_3phase_SmallStrain::GetDerivativeDAirFlowDGradpa( const Matrix& DN_DX_DISP, double airPressure, double capillaryPressure )
    {

        double saturation = GetSaturation( capillaryPressure );

        //Calculation of Derivative relative Permeability after Mualem
        double relPerm = 1.0 - saturation;

        if ( relPerm <= 0.01 )
        {
            relPerm = 0.01;
        }

// For Liakopolous Benschmark
// double relSat= (saturation-0.2)/0.8;
// double relPerm= 0.0001+ pow((1.0-relSat),2)*(1-pow(relSat,5.0/3.0));

        double airDensity = GetDensityAir( airPressure );

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        double result( dim );

        result = -relPerm * GetProperties()[PERMEABILITY_AIR] / ( airDensity * 9.81 );

        return result;
    }

    //************************************************************************************
    //************************************************************************************
    //STRESSES, STRAINS AND CONSTITUTIVE MODELL (UNSATURATED CASE)
    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateEffectiveStress( Vector&
            StressVector, Matrix& tanC_W, Matrix& tanC_A, const double waterPressure, const double airPressure )
    {
        double capillaryPressure = airPressure - waterPressure;

        double saturation = GetSaturation( capillaryPressure );

        double DSDpc = GetDerivativeDSaturationDpc( capillaryPressure );

        noalias( tanC_W ) = ZeroMatrix( 3, 3 );

        noalias( tanC_A ) = ZeroMatrix( 3, 3 );

        for ( unsigned int i = 0; i < 3; i++ )
        {
            StressVector( i ) =
                StressVector( i ) -
                (( 1.0 - saturation ) * airPressure + saturation * waterPressure );
        }

        for ( unsigned int i = 0; i < 3; i++ )
        {
            tanC_W( i, i ) = ( DSDpc * ( waterPressure - airPressure )
                               - saturation );

            tanC_A( i, i ) =
                ( DSDpc * ( airPressure - waterPressure ) - ( 1.0 - saturation ) );
        }
    }

    //************************************************************************************
    //************************************************************************************

    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStressAndTangentialStiffnessUnsaturatedSoils( Vector& StressVector, Matrix& tanC_U, Matrix& tanC_W, Matrix& tanC_A, Vector& StrainVector, double waterPressure, double airPressure, int PointNumber, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        if ( tanC_W.size1() != dim || tanC_W.size2() != dim )
            tanC_W.resize( dim, dim, false );

        noalias( tanC_W ) = ZeroMatrix( dim, dim );

        if ( tanC_A.size1() != dim || tanC_A.size2() != dim )
            tanC_A.resize( dim, dim, false );

        noalias( tanC_A ) = ZeroMatrix( dim, dim );

        if ( tanC_U.size1() != 6 || tanC_U.size2() != 6 )
            tanC_U.resize( 6, 6 );

        noalias( tanC_U ) = ZeroMatrix( 6, 6 );

        if ( StressVector.size() != 6 )
            StressVector.resize( 6 );

        //Set suction in const. law
        mConstitutiveLawVector[PointNumber]->SetValue( SUCTION, ( airPressure - waterPressure ),  rCurrentProcessInfo );

        mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
            StrainVector,
            ZeroMatrix( 1 ),
            StressVector,
            tanC_U,
            rCurrentProcessInfo,
            GetProperties(),
            GetGeometry(),
            row( GetGeometry().ShapeFunctionsValues(), PointNumber ),
            true,
            1,
            true );

        CalculateEffectiveStress( StressVector, tanC_W, tanC_A, waterPressure, airPressure );

        KRATOS_CATCH( "" )
    }

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Computes the strain vector
     */
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector )
    {
        KRATOS_TRY
        noalias( StrainVector ) = ZeroVector( 6 );

        for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
        {
            for ( unsigned int item = 0; item < 6; item++ )
                for ( unsigned int dim = 0; dim < 3; dim++ )
                    StrainVector[item] += B( item, 3 * node + dim ) * Displacements( node, dim );
        }

        KRATOS_CATCH( "" )
    }


    //************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::CalculateBoperator
    ( Matrix& B_Operator, const Matrix& DN_DX )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();

//         if(B_Operator.size() != number_of_nodes)
//             B_Operator.resize(number_of_nodes);
        noalias( B_Operator ) = ZeroMatrix( 6, number_of_nodes * 3 );

        for ( unsigned int i = 0;i < number_of_nodes;i++ )
        {
//             if(B_Operator[i].size1() != 6 || B_Operator[i].size2() != 3)
//                 B_Operator[i].resize(6,3);

//             noalias(B_Operator[i])= ZeroMatrix(6,3);

            B_Operator( 0, i*3 ) = DN_DX( i, 0 );
            B_Operator( 1, i*3 + 1 ) = DN_DX( i, 1 );
            B_Operator( 2, i*3 + 2 ) = DN_DX( i, 2 );
            B_Operator( 3, i*3 ) = DN_DX( i, 1 );
            B_Operator( 3, i*3 + 1 ) = DN_DX( i, 0 );
            B_Operator( 4, i*3 + 1 ) = DN_DX( i, 2 );
            B_Operator( 4, i*3 + 2 ) = DN_DX( i, 1 );
            B_Operator( 5, i*3 ) = DN_DX( i, 2 );
            B_Operator( 5, i*3 + 2 ) = DN_DX( i, 0 );
        }

        return;

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    void UnsaturatedSoilsElement_3phase_SmallStrain::InitializeMaterial
    ()
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial(
                GetProperties(), GetGeometry(), ZeroVector( 0 ) );
        }

        KRATOS_CATCH( "" )
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        if ( rVariable == ELASTIC_LEFT_CAUCHY_GREEN_OLD )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size1() != 3 || rValues[i].size2() != 3 )
                    rValues[i].resize( 3, 3 );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( ELASTIC_LEFT_CAUCHY_GREEN_OLD, rValues[i] );
            }
        }
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
// std::cout<<"GetValue On Integration Points"<<std::endl;
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );
        Matrix DN_DX_DISP( GetGeometry().size(), GetGeometry().WorkingSpaceDimension() );

        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        if ( rVariable == INSITU_STRESS )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != 6 )
                    rValues[i].resize( 6 );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INSITU_STRESS, rValues[i] );
            }
        }

        //To Plot Internal variables
        if ( rVariable == INTERNAL_VARIABLES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != 9 )
                    rValues[i].resize( 9 );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INTERNAL_VARIABLES, rValues[i] );
            }
        }

        //To Plot Stresses
        if ( rVariable == STRESSES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != 6 )
                    rValues[i].resize( 6 );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );
            }
        }

        //To Plot Fluid Flows
        if ( rVariable == FLUID_FLOWS )
        {
            unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );

            Vector N_PRESS( number_of_nodes_press );

            Matrix DN_DX_PRESS( number_of_nodes_press, 3 );

            double capillaryPressure;

            double waterPressure;

            double airPressure;

            double saturation;

            Vector airFlow( 3 );

            Vector waterFlow( 3 );

            const GeometryType::ShapeFunctionsGradientsType& DN_De_Displacement =
                GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
            const GeometryType::ShapeFunctionsGradientsType& DN_De_Pressure =
                mpPressureGeometry->ShapeFunctionsLocalGradients( mThisIntegrationMethod );

            const Matrix& Ncontainer_Pressure = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

            for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
            {
                if ( rValues[PointNumber].size() != 9 )
                    rValues[PointNumber].resize( 9 );

                // Shape Functions on current spatial quadrature point
                if ( N_PRESS.size() != number_of_nodes_press )
                    N_PRESS.resize( number_of_nodes_press );

                noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

                GetPressures( N_PRESS, capillaryPressure, waterPressure, airPressure );

                saturation = GetSaturation( capillaryPressure );

                rValues[PointNumber][0] = waterPressure;

                rValues[PointNumber][1] = airPressure;

                rValues[PointNumber][2] = saturation;

                noalias( DN_DX_PRESS ) = prod( DN_De_Pressure[PointNumber], mInvJ0[PointNumber] );

                noalias( DN_DX_DISP ) = prod( DN_De_Displacement[PointNumber], mInvJ0[PointNumber] );

                noalias( waterFlow ) = GetFlowWater( DN_DX_PRESS, DN_DX_DISP, capillaryPressure );

                noalias( airFlow ) = GetFlowAir( DN_DX_PRESS, DN_DX_DISP, airPressure, capillaryPressure );

                rValues[PointNumber][3] = waterFlow( 0 );

                rValues[PointNumber][4] = waterFlow( 1 );

                rValues[PointNumber][5] = waterFlow( 2 );

                rValues[PointNumber][6] = airFlow( 0 );

                rValues[PointNumber][7] = airFlow( 1 );

                rValues[PointNumber][8] = airFlow( 2 );
            }
        }

        //To Plot Coordinates of Integration Points
        if ( rVariable == COORDINATES )
        {
            for ( unsigned int i = 0; i < integration_points.size(); i++ )
            {
                if ( rValues[i].size() != 3 )
                    rValues[i].resize( 3 );

                Geometry<Node<3> >::CoordinatesArrayType dummy;

                GetGeometry().GlobalCoordinates( dummy, integration_points[i] );

                noalias( rValues[i] ) = dummy;
            }
        }

// std::cout<<"END::GetValue On Integration Points"<<std::endl;
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int number_of_nodes_press = ( mNodesPressMax - mNodesPressMin + 1 );

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        if ( rValues.size() != integration_points.size() )
            rValues.resize( integration_points.size() );

        const Matrix& Ncontainer_Pressure = mpPressureGeometry->ShapeFunctionsValues( mThisIntegrationMethod );

        Vector N_PRESS( number_of_nodes_press );

        double capillaryPressure;

        double waterPressure;

        double airPressure;

        double saturation;

        /////////////////////////////////////////////////////////////////////////
//// Integration in space sum_(beta=0)^(number of quadrature points)
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            // Shape Functions on current spatial quadrature point
            if ( N_PRESS.size() != number_of_nodes_press )
                N_PRESS.resize( number_of_nodes_press );

            noalias( N_PRESS ) = row( Ncontainer_Pressure, PointNumber );

            GetPressures( N_PRESS, capillaryPressure, waterPressure, airPressure );

            saturation = GetSaturation( capillaryPressure );

            if ( rVariable == SATURATION )
            {
                rValues[PointNumber] = saturation;
            }

            if ( rVariable == WATER_PRESSURE )
            {
                rValues[PointNumber] = waterPressure;
            }

            if ( rVariable == AIR_PRESSURE )
            {
                rValues[PointNumber] = airPressure;
            }
        }

        KRATOS_CATCH( "" )
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            std::cout << "wrong size: " << rValues.size() << "!=" << mConstitutiveLawVector.size() << std::endl;
            return;
        }

        if ( rVariable == ELASTIC_LEFT_CAUCHY_GREEN_OLD )
        {

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                mConstitutiveLawVector[i]->SetValue( ELASTIC_LEFT_CAUCHY_GREEN_OLD, rValues[i], rCurrentProcessInfo );
            }
        }
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            std::cout << "wrong size: " << rValues.size() << "!=" << mConstitutiveLawVector.size() << std::endl;
            return;
        }

        if ( rVariable == INSITU_STRESS )
        {

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                mConstitutiveLawVector[i]->SetValue( INSITU_STRESS, rValues[i], rCurrentProcessInfo );
            }
        }
    }

    void UnsaturatedSoilsElement_3phase_SmallStrain::SetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == INSITU_STRESS_SCALE )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                mConstitutiveLawVector[i]->SetValue( INSITU_STRESS_SCALE, rValues[i], rCurrentProcessInfo );
            }
        }

        if ( rVariable == OVERCONSOLIDATION_RATIO )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                mConstitutiveLawVector[i]->SetValue( OVERCONSOLIDATION_RATIO, rValues[i], rCurrentProcessInfo );
            }
        }
    }

    UnsaturatedSoilsElement_3phase_SmallStrain::IntegrationMethod UnsaturatedSoilsElement_3phase_SmallStrain::GetIntegrationMethod()
    {
        return mThisIntegrationMethod;
    }
} // Namespace Kratos


