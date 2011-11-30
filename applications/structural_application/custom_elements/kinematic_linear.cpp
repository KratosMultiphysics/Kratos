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
/* *********************************************************
*
*   Last Modified by:    $Author: nagel $
*   Date:                $Date: 2009-03-20 08:55:48 $
*   Revision:            $Revision: 1.16 $
*
* ***********************************************************/

// System includes


// External includes

// Project includes
// #include "includes/define.h"
#include "custom_elements/kinematic_linear.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application.h"
#include <boost/timer.hpp>

namespace Kratos
{
    KinematicLinear::KinematicLinear( IndexType NewId,
                                      GeometryType::Pointer pGeometry )
            : Element( NewId, pGeometry )
    {
        mIsInitialized = false;
        //DO NOT ADD DOFS HERE!!!
        //THIS IS THE DEFAULT CONSTRUCTOR
    }

    /**
    * A simple kinematic linear 3D element for the solution
    * of the momentum balance in structural mechanics.
    * This element is used for students training at the Ruhr University Bochum.
    * Therefore it may includes comments that are obvious for the
    * experienced user.
    */
    KinematicLinear::KinematicLinear( IndexType NewId,
                                      GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
            : Element( NewId, pGeometry, pProperties )
    {
        mIsInitialized = false;
        //DEFINE HERE THE INTEGRATION METHOD (the number of integration points to be used)
        //see /kratos/source/integration_rules.cpp for further details
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();//default method
        //m... declares that this is a member variable initialized in the header file
        //select Integration Method for 8 node hexahedron

        if ( GetGeometry().size() == 8 )
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;//methods for hexahedra elements

        //select Integration Method for 27 and 20 nodes
        if ( GetGeometry().size() == 27 || GetGeometry().size() == 20 )
            mThisIntegrationMethod = GeometryData::GI_GAUSS_3;//methods for hexahedra elements

        //select Integration Method for 10 and 4 node tetrahedron
        if ( GetGeometry().size() == 4 || GetGeometry().size() == 10 )
            mThisIntegrationMethod = GeometryData::GI_GAUSS_2;//methods for tetrahedra elements
    }

    Element::Pointer KinematicLinear::Create( IndexType NewId,
            NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new KinematicLinear( NewId, GetGeometry().Create( ThisNodes ),
                                 pProperties ) );
    }

    KinematicLinear::~KinematicLinear()
    {
    }

    /**
    * Initialization of the element, called at the begin of each simulation.
    * Membervariables and the Material law are initialized here
    */
    void KinematicLinear::Initialize()
    {
        KRATOS_TRY//EXCEPTION HANDLING (see corresponing KRATOS_CATCH("") )

        //dimension of the problem
        unsigned int dim = GetGeometry().WorkingSpaceDimension();


        if ( mIsInitialized )
        {
            //Set Up Initial displacement for StressFreeActivation of Elements
            mInitialDisp.resize( GetGeometry().size(), dim, false );

            for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
                for ( unsigned int i = 0; i < 3; i++ )
                    mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];
            return;
        }

        //number of integration points used, mThisIntegrationMethod refers to the
        //integration method defined in the constructor
        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        //initializing the Jacobian, the inverse Jacobian and Jacobians determinant in the reference
        // configuration
        GeometryType::JacobiansType J0( integration_points.size() );

        mInvJ0.resize( integration_points.size() );

        mTotalDomainInitialSize = 0.00;

        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
            mInvJ0[i].resize( dim, dim, false );
            noalias( mInvJ0[i] ) = ZeroMatrix( dim, dim );
        }

        mDetJ0.resize( integration_points.size(), false );

        noalias( mDetJ0 ) = ZeroVector( integration_points.size() );

        //calculating the Jacobian
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );

        //calculating the inverse Jacobian

        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            //getting informations for integration
            double IntegrationWeight = integration_points[PointNumber].Weight();
            //calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );
            //calculating the total area
            mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
        }

        //Set Up Initial displacement for StressFreeActivation of Elements
        mInitialDisp.resize( GetGeometry().size(), dim, false );

        for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
            for ( unsigned int i = 0; i < 3; i++ )
                mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

//                         mInitialDisp(node, i)= 0.0;
        //Initialization of the constitutive law vector and
        // declaration, definition and initialization of the material
        // lwas at each integration point
        if ( mConstitutiveLawVector.size() != integration_points.size() )
        {
            mConstitutiveLawVector.resize( integration_points.size() );
        }

        InitializeMaterial();

        mIsInitialized = true;


        KRATOS_CATCH( "" )
    }

    /**
    * Calculate double Variables at each integration point, used for postprocessing etc.
    * @param rVariable Global name of the variable to be calculated
    * @param output Vector to store the values on the qudrature points, output of the method
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<double>& rVariable, Vector& Output,
            const ProcessInfo& rCurrentProcessInfo )
    {
        if ( Output.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
            Output.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
            Output[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, Output[ii] );
    }

    /**
    * Calculate Vector Variables at each integration point, used for postprocessing etc.
    * @param rVariable Global name of the variable to be calculated
    * @param output Vector to store the values on the qudrature points, output of the method
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<Vector>& rVariable,
            std::vector<Vector>& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        GetValueOnIntegrationPoints( rVariable, Output, rCurrentProcessInfo );
    }

    /**
    * Calculate Matrix Variables at each integration point, used for postprocessing etc.
    * @param rVariable Global name of the variable to be calculated
    * @param output Vector to store the values on the qudrature points, output of the method
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<Matrix>& rVariable,
            std::vector<Matrix>& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = 3;
        unsigned int StrainSize = 3 * ( dim - 1 );

        //Initialize local variables
        Matrix B( StrainSize, number_of_nodes*dim );
        Matrix TanC( StrainSize, StrainSize );
        Vector StrainVector( StrainSize );
        Vector StressVector( StrainSize );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix CurrentDisp( number_of_nodes, dim );

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( mThisIntegrationMethod );
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        if ( Output.size() != integration_points.size() )
            Output.resize( integration_points.size() );

        const GeometryType::ShapeFunctionsGradientsType& DN_De =
            GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        //Declaration of the integration weight
        double Weight;

        //loop over all integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, DN_DX );

            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );
            //assign the integration weight at the current integration point
            Weight = integration_points[PointNumber].Weight();

            //calculate material response
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                StrainVector,
                ZeroMatrix( 1 ),
                StressVector,
                TanC,
                rCurrentProcessInfo,
                GetProperties(),
                GetGeometry(),
                row( Ncontainer, PointNumber ),
                true,
                0,
                true );

            if ( Output[PointNumber].size2() != StrainVector.size() )
                Output[PointNumber].resize( 1, StrainVector.size(), false );

            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            {
                for ( unsigned int ii = 0; ii < StrainVector.size(); ii++ )
                    Output[PointNumber]( 0, ii ) = StrainVector[ii];
            }
            else if ( rVariable == PK2_STRESS_TENSOR )
            {
                for ( unsigned int ii = 0; ii < StrainVector.size(); ii++ )
                    Output[PointNumber]( 0, ii ) = StressVector[ii];
            }
            else if ( rVariable == INSITU_STRESS )
            {
                Vector dummy = row( Output[PointNumber], 0 );
                row( Output[PointNumber], 0 ) = mConstitutiveLawVector[PointNumber]->GetValue( INSITU_STRESS, dummy );
            }

//             std::cout << StressVector[2] << "\t";
        }

//         std::cout << std::endl;
        KRATOS_CATCH( "" )
    }

    /**
    * Initialization of the Material law at each integration point
    */
    void KinematicLinear::InitializeMaterial()
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }

        KRATOS_CATCH( "" )
    }

    void KinematicLinear::ResetConstitutiveLaw()
    {
        KRATOS_TRY

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
        }

        KRATOS_CATCH( "" )
    }

    /**
    * THIS is the main method here the integration in space (loop over the integration points) is done,
    * the algorithmic tangent and the (inner and outer) load vector is computed
    * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_nodes*dim)*(number_of_nodes*dim)
    * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
    * @param rCurrentProcessInfo
    * @param CalculateStiffnessMatrixFlag true: algorithmic tangent has to be computed
    * @param CalculateResidualVectorFlag true: load vector has to be computed
    */
    void KinematicLinear::CalculateAll( MatrixType& rLeftHandSideMatrix,
                                        VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo,
                                        bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = 3;
        unsigned int StrainSize = 3 * ( dim - 1 );

        //this is the size of the elements stiffness matrix/force vector
        unsigned int MatSize = ( GetGeometry().size() * GetGeometry().WorkingSpaceDimension() );

        //Initialize local variables
        Matrix B( StrainSize, MatSize );
        Matrix TanC( StrainSize, StrainSize );
        Vector StrainVector( StrainSize );
        Vector StressVector( StrainSize );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix CurrentDisp( number_of_nodes, dim );


        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            //resize the LHS=StiffnessMatrix if its size is not correct
            if ( rLeftHandSideMatrix.size1() != MatSize )
                rLeftHandSideMatrix.resize( MatSize, MatSize, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            //resize the RHS=force vector if its size is not correct
            if ( rRightHandSideVector.size() != MatSize )
                rRightHandSideVector.resize( MatSize, false );

            noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
        }

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
        {
            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, DN_DX );

            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );

            //calculate stress
//             mConstitutiveLawVector[PointNumber]->UpdateMaterial( StrainVector,
//                     GetProperties(), GetGeometry(), row(Ncontainer,PointNumber), rCurrentProcessInfo );
//             mConstitutiveLawVector[PointNumber]->CalculateStress(StrainVector, StressVector);
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                StrainVector,
                ZeroMatrix( 1 ),
                StressVector,
                TanC,
                rCurrentProcessInfo,
                GetProperties(),
                GetGeometry(),
                row( Ncontainer, PointNumber ),
                true,
                ( int )CalculateStiffnessMatrixFlag,
                true );

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                //calculate material tangent
//                 mConstitutiveLawVector[PointNumber]->CalculateConstitutiveMatrix(StrainVector, TanC);
                //calculate stiffness matrix
                noalias( rLeftHandSideMatrix ) +=
                    prod( trans( B ), ( integration_points[PointNumber].Weight() * mDetJ0[PointNumber] ) * Matrix( prod( TanC, B ) ) );
                //CalculateStiffnesMatrix(rLeftHandSideMatrix, TanC, msB, Weight, mDetJ0[PointNumber]);
            }

            if ( CalculateResidualVectorFlag == true )
            {
//                 //contribution to external forces
//                 BodyForce = GetProperties()[BODY_FORCE];
                //Calculation of Loadvector
                AddBodyForcesToRHS( rRightHandSideVector, row( Ncontainer, PointNumber ), integration_points[PointNumber].Weight(), mDetJ0[PointNumber] );
                AddInternalForcesToRHS( rRightHandSideVector, B, StressVector, integration_points[PointNumber].Weight(), mDetJ0[PointNumber] );
            }
        }//loop over integration points

        KRATOS_CATCH( "" )
    }

    /**
    * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
    * method with CalculateStiffnessMatrixFlag = false and CalculateResidualVectorFlag = true
    * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::CalculateRightHandSide( VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = false;
        bool CalculateResidualVectorFlag = true;
        MatrixType temp = Matrix();

        CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag );
    }

    /**
    * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
    * method with CalculateStiffnessMatrixFlag = true and CalculateResidualVectorFlag = true
    * @param rLeftHandSideMatrix algorithmic tangent, size (number_of_nodes*dim)*(number_of_nodes*dim)
    * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
            VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        //calculation flags
        bool CalculateStiffnessMatrixFlag = true;
        bool CalculateResidualVectorFlag = true;
        CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo,
                      CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
    }

    /**
    * THIS method is called from the scheme after each solution step, here the time step
    * start and end point variables can be transferred n --> n+1
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        for ( unsigned int Point = 0; Point < integration_points.size(); Point++ )
        {
            mConstitutiveLawVector[Point]->FinalizeSolutionStep( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point ), CurrentProcessInfo );
        }

//                 Matrix Dummy_Matrix(3,3);
//                 noalias(Dummy_Matrix)= mConstitutiveLawVector[0]->GetValue(PK2_STRESS_TENSOR);
//                 Vector Dummy_Vector(5);
//                 noalias(Dummy_Vector)= mConstitutiveLawVector[0]->GetValue(INTERNAL_VARIABLES);
//                 for(unsigned int i=0; i< GetGeometry().size(); i++)
//                 {
//                         GetGeometry()[i].GetSolutionStepValue(MOMENTUM_X)= Dummy_Vector(0);
//                         GetGeometry()[i].GetSolutionStepValue(MOMENTUM_Y)= Dummy_Vector(1);
//                         GetGeometry()[i].GetSolutionStepValue(MOMENTUM_Z)= Dummy_Vector(2);
//                         GetGeometry()[i].GetSolutionStepValue(PRESSURE)= Dummy_Vector(3);
//                         GetGeometry()[i].GetSolutionStepValue(ERROR_RATIO)= Dummy_Vector(4);
//                 }

    }

    /**
    * THIS method is called from the scheme at the start of each solution step
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        for ( unsigned int Point = 0; Point < integration_points.size(); Point++ )
        {
            mConstitutiveLawVector[Point]->InitializeSolutionStep( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point ), CurrentProcessInfo );

//            mConstitutiveLawVector[Point]->SetValue(SUCTION, GetGeometry()[0].GetSolutionStepValue(TEMPERATURE),CurrentProcessInfo );
        }
    }

    void KinematicLinear::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int MatSize = dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != MatSize )
            rMassMatrix.resize( MatSize, MatSize, false );

        rMassMatrix = ZeroMatrix( MatSize, MatSize );

        double TotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];;

        if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

        Vector LumpFact;

        LumpFact = GetGeometry().LumpingFactors( LumpFact );

        for ( unsigned int i = 0; i < NumberOfNodes; i++ )
        {
            double temp = LumpFact[i] * TotalMass;

            for ( unsigned int j = 0; j < dimension; j++ )
            {
                unsigned int index = i * dimension + j;
                rMassMatrix( index, index ) = temp;
            }
        }

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    void KinematicLinear::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int MatSize = number_of_nodes * dim;

        if ( rDampMatrix.size1() != MatSize )
            rDampMatrix.resize( MatSize, MatSize, false );

        noalias( rDampMatrix ) = ZeroMatrix( MatSize, MatSize );

        Matrix StiffnessMatrix = ZeroMatrix( MatSize, MatSize );

        Vector RHS_Vector = ZeroVector( MatSize );

        //rayleigh damping
        CalculateAll( StiffnessMatrix, RHS_Vector, rCurrentProcessInfo, true, false );

        double alpha = 0.001;

        double beta = 0.001;

        MassMatrix( rDampMatrix, rCurrentProcessInfo );

        noalias( rDampMatrix ) = alpha * StiffnessMatrix + beta * rDampMatrix;

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    void KinematicLinear::GetValuesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize )    values.resize( MatSize, false );

        for ( unsigned int i = 0;i < number_of_nodes;i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
        }
    }


    //************************************************************************************
    //************************************************************************************
    void KinematicLinear::GetFirstDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize )   values.resize( MatSize, false );

        for ( unsigned int i = 0;i < number_of_nodes;i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
        }
    }

    //************************************************************************************
    //************************************************************************************
    void KinematicLinear::GetSecondDerivativesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int MatSize = number_of_nodes * dim;

        if ( values.size() != MatSize ) values.resize( MatSize, false );

        for ( unsigned int i = 0;i < number_of_nodes;i++ )
        {
            unsigned int index = i * dim;
            values[index] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
            values[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

            if ( dim == 3 )
                values[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
        }
    }

    /**
    * returns the used integration method
    */
    KinematicLinear::IntegrationMethod KinematicLinear::GetIntegrationMethod()
    {
        return mThisIntegrationMethod;
    }

    /**
    * not used
    */
    void KinematicLinear::CalculateAndAddExtForceContribution( const Vector& N,
            const ProcessInfo& CurrentProcessInfo, Vector& BodyForce, VectorType& rRightHandSideVector,
            double weight )
    {
        KRATOS_TRY

        KRATOS_CATCH( "" )
    }

    /**
    * Informations for the builder and solver to assemble the global vectors and matrices.
    * Here a Vector containing the EquationIds of the differnt Dofs is created
    * @param rResult Vector of the EquationIds
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::EquationIdVector( EquationIdVectorType& rResult,
                                            ProcessInfo& CurrentProcessInfo )
    {
        unsigned int dim = ( GetGeometry().WorkingSpaceDimension() );//3 displacement dofs
        unsigned int MatSize = GetGeometry().size() * dim;

        if ( rResult.size() != MatSize )
            rResult.resize( MatSize, false );

        for ( unsigned int i = 0 ;i < GetGeometry().size() ;i++ )
        {
            int index = i * dim;
            rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
            rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
        }
    }

    /**
    * Informations for the builder and solver to assemble the global vectors and matrices.
    * Here a Container containing the pointers rto the DOFs of this element is created
    * @param ElementalDofList Container with of the DOFs associated with the nodes
    *                           of this element
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo&
                                      CurrentProcessInfo )
    {
        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ;i < GetGeometry().size() ;i++ )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
        }
    }

    /**
    * Adds the Body Forces to the load vector
    * @param R RHS Vector
    * @param N_DISP shape function values at the current integration points
    * @param Weight current integration weight
    * @param detJ current Determinant of the Jacobian
    */
    void KinematicLinear::AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        Vector gravity( dim );

        noalias( gravity ) = GetProperties()[GRAVITY];

        double density = GetProperties()[DENSITY];

        for ( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                R( prim*dim + i ) += N_DISP( prim ) * density * gravity( i ) * detJ * Weight;
            }
        }

        KRATOS_CATCH( "" )
    }

    /**
    * Adds the Internal Forces to the load vector
    * @param R RHS Vector
    * @param B_Operator B-Operator at the current integration point
    * @param StressVector current stress vector
    * @param Weight current integration weight
    * @param detJ current Determinant of the Jacobian
    */
    void KinematicLinear::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                for ( unsigned int gamma = 0; gamma < 6; gamma++ )
                {
                    R( prim*dim + i ) += ( -1 ) * ( B_Operator( gamma, prim * 3 + i ) * StressVector( gamma ) *
                                                    detJ * Weight );
                }
            }
        }

//         noalias(R) -= detJ*Weight* prod(trans(B_Operator),StressVector);

        KRATOS_CATCH( "" )
    }

    /**
    * Adds the Contribution of the current quadrature point to the load vector
    * @param K LHS Matrix
    * @param tan_C 6*6 algorithmic tangent of the materia law (derivation of stresses
    *               regarding strains
    * @param B_Operator B-Operator at the current integration point
    * @param Weight current integration weight
    * @param detJ current Determinant of the Jacobian
    */
    void KinematicLinear::CalculateStiffnesMatrix( Matrix& K, const
            Matrix& tan_C, const Matrix& B_Operator, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int prim = 0; prim < GetGeometry().size(); prim++ )
        {
            for ( unsigned int i = 0; i < dim; i++ )
            {
                for ( unsigned int sec = 0; sec < GetGeometry().size(); sec++ )
                {
                    for ( unsigned int j = 0; j < dim; j++ )
                    {
                        for ( unsigned int alpha = 0; alpha < 6; alpha++ )
                            for ( unsigned int beta = 0; beta < 6; beta++ )
                                K( prim*dim + i, sec*dim + j ) += B_Operator( alpha, 3 * prim + i )
                                                                  * tan_C( alpha, beta ) * B_Operator( beta, 3 * sec + j ) * detJ * Weight;
                    }
                }
            }
        }

//         noalias(K)-= prod(trans(B_Operator),(Weight*detJ)*Matrix(prod(C,B_Operator)) );

        KRATOS_CATCH( "" )
    }

    /**
    * Computes the stress vector and the algorithmic tangent at the current quadrature points
    * regarding the current strain vector
    * @param StressVector output of the method, Cauchy stress vector
    * @param tanC_U algorithmic tangen (derivative Cauchy stress vector/StrainVector) [6*6]
    * @param StrainVector output: current strain vector
    * @param B_Operator current B-operator
    * @param PointNumber number of the current integration point
    * @param CurrentProcessInfo
    */
    void KinematicLinear::CalculateStressAndTangentialStiffness( Vector& StressVector, Matrix& tanC_U,
            Vector& StrainVector, const Matrix& B_Operator, int PointNumber,
            const ProcessInfo& CurrentProcessInfo )
    {
        KRATOS_TRY

        //Calculate the current strain vector using the B-Operator
        //CalculateStrainVector(B_Operator, StrainVector, PointNumber);
        //Call the material law to get the current Stress Vector
        //mConstitutiveLawVector[PointNumber]->CalculateStress(StrainVector, StressVector);
        //Call the material law to get the current algorithmic tangent
        //mConstitutiveLawVector[PointNumber]->CalculateConstitutiveMatrix(StrainVector, tanC_U);

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the strain vector
     */
    void KinematicLinear::CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector )
    {
        KRATOS_TRY
        noalias( StrainVector ) = ZeroVector( 6 );

        for ( unsigned int node = 0; node < GetGeometry().size(); node++ )
        {
            for ( unsigned int item = 0; item < 6; item++ )
                for ( unsigned int dim = 0; dim < 3; dim++ )
                    StrainVector[item] += B( item, 3 * node + dim ) * ( Displacements( node, dim ) - mInitialDisp( node, dim ) );
        }

        KRATOS_CATCH( "" )
    }

    /**
    * Computes the B-Operator at the current quadrature point
    * @param B_Operator current B-operator
    * @param DN_DX shape function values at the current integration point
    */
    void KinematicLinear::CalculateBoperator
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

    /**
    * Calculate Matrix Variables at each integration point, used for postprocessing etc.
    * @param rVariable Global name of the variable to be calculated
    * @param rValues Vector to store the values on the qudrature points, output of the method
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == PK2_STRESS_TENSOR || rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

        return;
    }

    /**
    * Calculate Vector Variables at each integration point, used for postprocessing etc.
    * @param rVariable Global name of the variable to be calculated
    * @param rValues Vector to store the values on the qudrature points, output of the method
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
// std::cout<<"GetValue On Integration Points"<<std::endl;
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        if ( rVariable == MATERIAL_PARAMETERS )
        {
            if ( rValues.size() !=
                    GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
                rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ii++ )
                rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
        }

        if ( rVariable == INSITU_STRESS )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != 6 )
                    rValues[i].resize( 6 );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INSITU_STRESS, rValues[i] );
            }
        }

        if ( rVariable == INTERNAL_VARIABLES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
            {
                if ( rValues[i].size() != 8 )
                    rValues[i].resize( 8 );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INTERNAL_VARIABLES, rValues[i] );
            }
        }

// std::cout<<"END::GetValue On Integration Points"<<std::endl;
    }

    /**
    * Calculate double Variables at each integration point, used for postprocessing etc.
    * @param rVariable Global name of the variable to be calculated
    * @param rValues Vector to store the values on the qudrature points, output of the method
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( !( rVariable == DP_EPSILON ) )
            return;

        if ( rValues.size() != GetGeometry().IntegrationPoints().size() )
            rValues.resize( GetGeometry().IntegrationPoints().size() );

        //reading integration points and local gradients
        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
        {
            rValues[Point] = mConstitutiveLawVector[Point]->GetValue( rVariable, rValues[Point] );
        }
    }

    /**
    * Set a Matrix Variable from outside
    * @param rVariable Global name of the variable to be calculated
    * @param rValues Vector of the values on the quadrature points
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::SetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
        {
            std::cout << "In KinematicLinear:SetValueOnIntegrationPoints(...) Line 798, wrong size of the input vector, is = " << rValues.size() << " ; should be = " << mConstitutiveLawVector.size() << std::endl;
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

    /**
    * Set a Vector Variable from outside
    * @param rVariable Global name of the variable to be calculated
    * @param rValues Vector of the values on the quadrature points
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::SetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
//              std::cout << mConstitutiveLawVector[0] << std::endl;
        if ( rVariable == INSITU_STRESS )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( INSITU_STRESS, rValues[PointNumber], rCurrentProcessInfo );
            }
        }

        if ( rVariable == PRESTRESS )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( PRESTRESS, rValues[PointNumber], rCurrentProcessInfo );
            }
        }

        if ( rVariable == MATERIAL_PARAMETERS )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( MATERIAL_PARAMETERS,
                        rValues[PointNumber], rCurrentProcessInfo );
            }
        }
    }

    /**
    * Set a Double Variable from outside
    * @param rVariable Global name of the variable to be calculated
    * @param rValue value on the quadrature points
    * @param rCurrentProcessInfo
    */
    void KinematicLinear::SetValueOnIntegrationPoints( const Variable<double>& rVariable,
            std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
//              std::cout << mConstitutiveLawVector[0] << std::endl;
        if ( rVariable == SUCTION )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); PointNumber++ )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( SUCTION, rValues[PointNumber],
                        rCurrentProcessInfo );
            }
        }

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

    int KinematicLinear::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        if ( this->Id() < 1 )
        {
            KRATOS_ERROR( std::logic_error, "Element found with Id 0 or negative", "" )
        }

        if ( this->GetGeometry().Area() < 0 )
        {
            std::cout << "error on element -> " << this->Id() << std::endl;
            KRATOS_ERROR( std::logic_error, "Area can not be less than 0", "" )
        }

        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); i++ )
        {
            return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
        }

        return 0;

        KRATOS_CATCH( "" );

    }



} // Namespace Kratos

