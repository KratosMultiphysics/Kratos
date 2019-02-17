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
 *   Last Modified by:    $Author: hbui $
 *   Date:                $Date: 2013-02-22 16:16:48 $
 *   Revision:            $Revision: 1.2 $
 *
 * ***********************************************************/

// System includes


// External includes

// Project includes
// #include "includes/define.h"
#include "kinematic_linear.h"
#include "utilities/math_utils.h"
#include "custom_utilities/bathe_recover_stress_utility.h"
#include "custom_utilities/sd_math_utils.h"
#include "structural_application.h"

// #define ENABLE_DEBUG_CONSTITUTIVE_LAW

//TODO: there is a potential bug at the CalculateRightHandSide, which is used for calculating the reaction. In principal, it should not tell the material to update themself, however, CalculateRightHandSide indirectly call CalculateMaterialResponse. THis should be fixed, by introducing another abstract layer to update the material in the input parameters for CalculateAll

namespace Kratos
{
    KinematicLinear::KinematicLinear( IndexType NewId, GeometryType::Pointer pGeometry )
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
        mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();//default method
    }

    Element::Pointer KinematicLinear::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new KinematicLinear( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }
    
    Element::Pointer KinematicLinear::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new KinematicLinear( NewId, pGeom, pProperties ) );
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

            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
                for ( unsigned int i = 0; i < dim; ++i )
                    mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
//            std::cout << "Element " << Id() << " mInitialDisp is reinitialized to " << mInitialDisp << std::endl;
            #endif

            return;
        }

        if(GetProperties().Has( INTEGRATION_ORDER ) == true)
        {
            if(GetProperties()[INTEGRATION_ORDER] == 1)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_1;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 2)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_2;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 3)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_3;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 4)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_4;
            }
            else if(GetProperties()[INTEGRATION_ORDER] == 5)
            {
                mThisIntegrationMethod = GeometryData::GI_GAUSS_5;
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "Does not support for more integration points", *this)
        }
        else
            mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod(); // default method

        //number of integration points used, mThisIntegrationMethod refers to the
        //integration method defined in the constructor
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        //initializing the Jacobian, the inverse Jacobian and Jacobians determinant in the reference
        // configuration
        GeometryType::JacobiansType J0( integration_points.size() );

        mInvJ0.resize( integration_points.size() );

        mTotalDomainInitialSize = 0.00;

        for ( unsigned int i = 0; i < integration_points.size(); ++i )
        {
            mInvJ0[i].resize( dim, dim, false );
            noalias( mInvJ0[i] ) = ZeroMatrix( dim, dim );
        }

        mDetJ0.resize( integration_points.size(), false );

        noalias( mDetJ0 ) = ZeroVector( integration_points.size() );

        //calculating the Jacobian
        J0 = GetGeometry().Jacobian( J0, mThisIntegrationMethod );

        //calculating the inverse Jacobian
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            //getting informations for integration
            double IntegrationWeight = integration_points[PointNumber].Weight();
            //calculating and storing inverse of the jacobian and the parameters needed
            MathUtils<double>::InvertMatrix( J0[PointNumber], mInvJ0[PointNumber], mDetJ0[PointNumber] );
            //calculating the total domain size
            mTotalDomainInitialSize += mDetJ0[PointNumber] * IntegrationWeight;
        }
        this->SetValue(GEOMETRICAL_DOMAIN_SIZE, mTotalDomainInitialSize);

        //Set Up Initial displacement for StressFreeActivation of Elements
        mInitialDisp.resize( GetGeometry().size(), dim, false );

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            for ( unsigned int i = 0; i < dim; ++i )
                mInitialDisp( node, i ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT )[i];

        #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
//        std::cout << "Element " << Id() << " mInitialDisp is initialized to " << mInitialDisp << std::endl;
        #endif

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
    void KinematicLinear::CalculateOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& Output, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( Output.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
            Output.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(), false );

        for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
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

        // TODO: This needs to be reviewed (BUI)

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();;
        unsigned int strain_size = dim * ( dim + 1 ) / 2;

        //Initialize local variables
        Matrix B( strain_size, number_of_nodes*dim );
        Matrix TanC( strain_size, strain_size );
        Vector StrainVector( strain_size );
        Vector StressVector( strain_size );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix CurrentDisp( number_of_nodes, dim );

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points =
            GetGeometry().IntegrationPoints( mThisIntegrationMethod );
        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        if ( Output.size() != integration_points.size() )
            Output.resize( integration_points.size() );

        const GeometryType::ShapeFunctionsGradientsType& DN_De =
            GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );

        //Declaration of the integration weight
        //    double Weight;

        //loop over all integration points
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, DN_DX );

            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );
            //assign the integration weight at the current integration point
            //        Weight = integration_points[PointNumber].Weight();

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
                true
            );

            if ( Output[PointNumber].size2() != StrainVector.size() )
                Output[PointNumber].resize( 1, StrainVector.size(), false );

            if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
            {
                for ( unsigned int ii = 0; ii < StrainVector.size(); ++ii )
                    Output[PointNumber]( 0, ii ) = StrainVector[ii];
            }
            else if ( rVariable == PK2_STRESS_TENSOR )
            {
                for ( unsigned int ii = 0; ii < StrainVector.size(); ++ii )
                    Output[PointNumber]( 0, ii ) = StressVector[ii];
            }
            else if ( rVariable == INSITU_STRESS )
            {
                Vector dummy = row( Output[PointNumber], 0 );
                row( Output[PointNumber], 0 ) = mConstitutiveLawVector[PointNumber]->GetValue( INSITU_STRESS, dummy );
            }

            //             std::cout << StressVector[2] << "\t";
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        // clean the geometry
        GetGeometry().Clean();
        #endif

        KRATOS_CATCH( "" )
    }

    /**
     * Initialization of the Material law at each integration point
     */
    void KinematicLinear::InitializeMaterial()
    {
        KRATOS_TRY

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
                mConstitutiveLawVector[i]->SetValue( PARENT_ELEMENT_ID, this->Id(), *(ProcessInfo*)0);
                mConstitutiveLawVector[i]->SetValue( INTEGRATION_POINT_INDEX, i, *(ProcessInfo*)0);
                mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif

        KRATOS_CATCH( "" )
    }

    void KinematicLinear::ResetConstitutiveLaw()
    {
        KRATOS_TRY

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->ResetMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif

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
    void KinematicLinear::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                       VectorType& rRightHandSideVector,
                                       ProcessInfo& rCurrentProcessInfo,
                                       bool CalculateStiffnessMatrixFlag,
                                       bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY

        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();;
        unsigned int strain_size = dim * ( dim + 1 ) / 2;

        //this is the size of the elements stiffness matrix/force vector
        unsigned int mat_size = GetGeometry().size() * dim;

        //Initialize local variables
        Matrix B( strain_size, mat_size );
        Matrix TanC( strain_size, strain_size );
        Vector StrainVector( strain_size );
        Vector StressVector( strain_size );
        Matrix DN_DX( number_of_nodes, dim );
        Matrix CurrentDisp( number_of_nodes, dim );

        //resize the LHS=StiffnessMatrix if its size is not correct
        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            if ( rLeftHandSideMatrix.size1() != mat_size )
                rLeftHandSideMatrix.resize( mat_size, mat_size, false );

            noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
        }

        //resizing as needed the RHS
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            //resize the RHS=force vector if its size is not correct
            if ( rRightHandSideVector.size() != mat_size )
                rRightHandSideVector.resize( mat_size, false );

            noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //initialize the geometry
        GetGeometry().Initialize(mThisIntegrationMethod);
        #endif

        //reading integration points and local gradients
        const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints( mThisIntegrationMethod );

        const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );

        const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

        //Current displacements
        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            noalias( row( CurrentDisp, node ) ) = GetGeometry()[node].GetSolutionStepValue( DISPLACEMENT );
        
        //auxiliary terms
        const Vector& BodyForce = GetProperties()[BODY_FORCE];

        /////////////////////////////////////////////////////////////////////////
        //// Integration in space over quadrature points
        /////////////////////////////////////////////////////////////////////////
        for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); ++PointNumber )
        {
            noalias( DN_DX ) = prod( DN_De[PointNumber], mInvJ0[PointNumber] );
            //Initializing B_Operator at the current integration point
            CalculateBoperator( B, DN_DX );

            //calculate strain
            CalculateStrain( B, CurrentDisp, StrainVector );

            #ifdef ENABLE_DEBUG_CONSTITUTIVE_LAW
            mConstitutiveLawVector[PointNumber]->SetValue(PARENT_ELEMENT_ID, this->Id(), rCurrentProcessInfo);
            mConstitutiveLawVector[PointNumber]->SetValue(INTEGRATION_POINT_INDEX, PointNumber, rCurrentProcessInfo);
            std::cout << "At element " << Id() << " integration point " << PointNumber << ":" << std::endl;
            std::cout << "B: " << B << std::endl;
            std::cout << "CurrentDisp: " << CurrentDisp << std::endl;
            std::cout << "mInitialDisp: " << mInitialDisp << std::endl;
            std::cout << "StrainVector: " << StrainVector << std::endl;
            #endif
            
            //calculate stress and consistent tangent
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
                true
            );

            //calculating weights for integration on the reference configuration
            double IntToReferenceWeight = integration_points[PointNumber].Weight();

            //modify integration weight in case of 2D
            if ( dim == 2 ) IntToReferenceWeight *= GetProperties()[THICKNESS];

            if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
            {
                //calculate stiffness matrix
                noalias( rLeftHandSideMatrix ) +=
                    prod( trans( B ), ( IntToReferenceWeight * mDetJ0[PointNumber] ) * Matrix( prod( TanC, B ) ) );
            }

            if ( CalculateResidualVectorFlag == true )
            {
                //contribution of external forces
                CalculateAndAdd_ExtForceContribution( row( Ncontainer, PointNumber ), rCurrentProcessInfo, BodyForce, rRightHandSideVector, IntToReferenceWeight, mDetJ0[PointNumber]);

                //contribution of gravity (if there is)
                AddBodyForcesToRHS( rRightHandSideVector, row( Ncontainer, PointNumber ), IntToReferenceWeight, mDetJ0[PointNumber] );
                
                //contribution of internal forces
                AddInternalForcesToRHS( rRightHandSideVector, B, StressVector, IntToReferenceWeight, mDetJ0[PointNumber] );
            }
        }//loop over integration points

        // modify the right hand side to account for prescribed displacement
        // according to the book of Bazant & Jirasek, this scheme is more stable than the current scheme for prescribing displacement.
        // // However, I have to temporarily disable it to keep the consistency.
        if(CalculateStiffnessMatrixFlag && CalculateResidualVectorFlag)
        {
            for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
            {
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_X))
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_X);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim) * temp;
                }
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_Y))
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Y);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 1) * temp;
                }
                if(GetGeometry()[node].IsFixed(DISPLACEMENT_Z) && dim == 3)
                {
                    double temp = GetGeometry()[node].GetSolutionStepValue(PRESCRIBED_DELTA_DISPLACEMENT_Z);
                    for( unsigned int i = 0; i < mat_size; ++i )
                        rRightHandSideVector[i] -= rLeftHandSideMatrix(i, node * dim + 2) * temp;
                }
            }
        }

        #ifdef ENABLE_BEZIER_GEOMETRY
        //clean the internal data of the geometry
        GetGeometry().Clean();
        #endif

        KRATOS_CATCH( "" )
    }

    /**
     * THIS method is called from the scheme during the iteration loop, it calls the CalculateAll()
     * method with CalculateStiffnessMatrixFlag = false and CalculateResidualVectorFlag = true
     * @param rRightHandSideVector (inner and outer) load vector, size (number_of_nodes*dim)
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo
    )
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
     * THIS method is called from the scheme at the start of each solution step
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
        {
            std::vector<Vector> Values;
            this->GetValueOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
            }
        }

        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            Vector dummy;
//            dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point );
            mConstitutiveLawVector[Point]->InitializeSolutionStep( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
        }
    }

    void KinematicLinear::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
        //reset all resistant forces at node
        for ( unsigned int i = 0; i < GetGeometry().size(); ++i )
        {
            GetGeometry()[i].GetSolutionStepValue( REACTION_X ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Y ) = 0.0;
            GetGeometry()[i].GetSolutionStepValue( REACTION_Z ) = 0.0;
        }

        if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
        {
            std::vector<Vector> Values;
            this->GetValueOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
            }
        }

        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            Vector dummy;
//            dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point );
            mConstitutiveLawVector[Point]->InitializeNonLinearIteration( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
        }
    }
    
    void KinematicLinear::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
    {
        if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
        {
            std::vector<Vector> Values;
            this->GetValueOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
            }
        }

        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            Vector dummy;
//            dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point );
            mConstitutiveLawVector[Point]->FinalizeNonLinearIteration( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
        }
    }
    
    /**
     * THIS method is called from the scheme after each solution step, here the time step
     * start and end point variables can be transferred n --> n+1
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        if(mConstitutiveLawVector[0]->Has(CURRENT_STRAIN_VECTOR))
        {
            std::vector<Vector> Values;
            this->GetValueOnIntegrationPoints(CURRENT_STRAIN_VECTOR, Values, CurrentProcessInfo);
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
            {
                mConstitutiveLawVector[Point]->SetValue( CURRENT_STRAIN_VECTOR, Values[Point], CurrentProcessInfo );
            }
        }

        for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); ++Point )
        {
            Vector dummy;
//            dummy = row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), Point );
            mConstitutiveLawVector[Point]->FinalizeSolutionStep( GetProperties(), GetGeometry(), dummy, CurrentProcessInfo );
        }
    }

    void KinematicLinear::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Deprecated method", "")
    }

    void KinematicLinear::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int NumberOfNodes = GetGeometry().size();
        unsigned int mat_size = dimension * NumberOfNodes;

        if ( rMassMatrix.size1() != mat_size )
            rMassMatrix.resize( mat_size, mat_size, false );

        rMassMatrix = ZeroMatrix( mat_size, mat_size );

        double TotalMass = 0.0;
        if( GetValue(USE_DISTRIBUTED_PROPERTIES) )
        {
            TotalMass = mTotalDomainInitialSize * GetValue(DENSITY);
        }
        else
            TotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];

        if ( dimension == 2 ) TotalMass *= GetProperties()[THICKNESS];

        Vector LumpFact;

        LumpFact = GetGeometry().LumpingFactors( LumpFact );

        for ( unsigned int i = 0; i < NumberOfNodes; ++i )
        {
            double temp = LumpFact[i] * TotalMass;

            for ( unsigned int j = 0; j < dimension; ++j )
            {
                unsigned int index = i * dimension + j;
                rMassMatrix( index, index ) = temp;
            }
        }

        KRATOS_CATCH( "" )
    }

    void KinematicLinear::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_THROW_ERROR(std::logic_error, "Deprecated method", "")
    }

    void KinematicLinear::CalculateDampingMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        //resizing as needed the LHS
        unsigned int mat_size = number_of_nodes * dim;

        if ( rDampMatrix.size1() != mat_size )
            rDampMatrix.resize( mat_size, mat_size, false );

        noalias( rDampMatrix ) = ZeroMatrix( mat_size, mat_size );

        Matrix StiffnessMatrix = ZeroMatrix( mat_size, mat_size );

        Vector RHS_Vector = ZeroVector( mat_size );

        //rayleigh damping
        CalculateAll( StiffnessMatrix, RHS_Vector, rCurrentProcessInfo, true, false );

        double alpha = 0.001, beta = 0.001;

        if(GetProperties().Has(RAYLEIGH_DAMPING_ALPHA))
        {
            alpha = GetProperties()[RAYLEIGH_DAMPING_ALPHA];
        }

        if(GetProperties().Has(RAYLEIGH_DAMPING_BETA))
        {
            beta = GetProperties()[RAYLEIGH_DAMPING_BETA];
        }

        CalculateMassMatrix( rDampMatrix, rCurrentProcessInfo );

        noalias( rDampMatrix ) = alpha * rDampMatrix + beta * StiffnessMatrix;

        KRATOS_CATCH( "" )
    }

    //************************************************************************************
    //************************************************************************************
    void KinematicLinear::GetValuesVector( Vector& values, int Step )
    {
        const unsigned int number_of_nodes = GetGeometry().size();
        const unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )    values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
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
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size )   values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
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
        unsigned int mat_size = number_of_nodes * dim;

        if ( values.size() != mat_size ) values.resize( mat_size, false );

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
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
    KinematicLinear::IntegrationMethod KinematicLinear::GetIntegrationMethod() const
    {
        return mThisIntegrationMethod;
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
        unsigned int dim = ( GetGeometry().WorkingSpaceDimension() );
        unsigned int mat_size = GetGeometry().size() * dim;

        if ( rResult.size() != mat_size )
            rResult.resize( mat_size, false );

        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            int index = i * dim;
            rResult[index] = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
            rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
            if(dim == 3)
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
        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        ElementalDofList.resize( 0 );

        for ( unsigned int i = 0 ; i < GetGeometry().size() ; ++i )
        {
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
            ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
            if(dim == 3)
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
    inline void KinematicLinear::AddBodyForcesToRHS( Vector& R, const Vector& N_DISP, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();

        Vector gravity( dim );

        double density = 0.0;
        if( GetValue( USE_DISTRIBUTED_PROPERTIES ) )
        {
            noalias( gravity ) = GetValue(GRAVITY);
            density = GetValue(DENSITY);
        }
        else
        {
            noalias( gravity ) = GetProperties()[GRAVITY];
            density = GetProperties()[DENSITY];
        }

        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                R( prim*dim + i ) += N_DISP( prim ) * density * gravity( i ) * detJ * Weight;
            }
        }

        KRATOS_CATCH( "" )
    }

    inline void KinematicLinear::CalculateAndAdd_ExtForceContribution(
            const Vector& N,
            const ProcessInfo& CurrentProcessInfo,
            const Vector& BodyForce,
            VectorType& rRightHandSideVector,
            double weight,
            double detJ
            )
    {
        KRATOS_TRY
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        for ( unsigned int i = 0; i < number_of_nodes; ++i )
        {
            int index = dimension * i;

            for ( unsigned int j = 0; j < dimension; ++j )
                rRightHandSideVector[index + j] += weight * detJ * N[i] * BodyForce[j];
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
//    void KinematicLinear::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
//    {
//        KRATOS_TRY

//        unsigned int dim = GetGeometry().WorkingSpaceDimension();
//        unsigned int strain_size = dim * (dim + 1) / 2;
// 
//        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
//        {
//            for ( unsigned int i = 0; i < dim; ++i )
//            {
//                for ( unsigned int gamma = 0; gamma < strain_size; ++gamma )
//                {
//                    R( prim * dim + i ) += ( -1 ) * ( B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) * detJ * Weight );
//                }
//            }
//        }
//        //         noalias(R) -= detJ * Weight * prod(trans(B_Operator), StressVector);
// 
//        KRATOS_CATCH( "" )
//    }

    void KinematicLinear::AddInternalForcesToRHS( Vector& R, const Matrix& B_Operator, Vector& StressVector, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;
        Vector InternalForces(3);
        
        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                InternalForces(i) = 0.0;
                for ( unsigned int gamma = 0; gamma < strain_size; ++gamma )
                {
                    InternalForces(i) += B_Operator( gamma, prim * dim + i ) * StressVector( gamma ) * detJ * Weight;
                }
                
                R( prim * dim + i ) -= InternalForces(i);
            }
//            GetGeometry()[prim].GetSolutionStepValue( REACTION ) += InternalForces;
        }
        
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
    void KinematicLinear::CalculateStiffnesMatrix( Matrix& K, const Matrix& tan_C, const Matrix& B_Operator, double Weight, double detJ )
    {
        KRATOS_TRY

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2 ;

        for ( unsigned int prim = 0; prim < GetGeometry().size(); ++prim )
        {
            for ( unsigned int i = 0; i < dim; ++i )
            {
                for ( unsigned int sec = 0; sec < GetGeometry().size(); ++sec )
                {
                    for ( unsigned int j = 0; j < dim; ++j )
                    {
                        for ( unsigned int alpha = 0; alpha < strain_size; ++alpha )
                            for ( unsigned int beta = 0; beta < strain_size; ++beta )
                                K( prim*dim + i, sec*dim + j ) += B_Operator( alpha, dim * prim + i )
                                    * tan_C( alpha, beta ) * B_Operator( beta, dim * sec + j ) * detJ * Weight;
                    }
                }
            }
        }

//        noalias(K) -= prod( trans(B_Operator), (Weight * detJ) * Matrix(prod(C, B_Operator)) );

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the stress vector and the algorithmic tangent at the current quadrature points
     * regarding the current strain vector
     * @param StressVector output of the method, Cauchy stress vector
     * @param tanC_U algorithmic tangen (derivative Cauchy stress vector/StrainVector) [strain_size*strain_size]
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

            // TODO ?

            KRATOS_CATCH( "" )
    }

    /**
     * Computes the strain vector
     */
    void KinematicLinear::CalculateStrain( const Matrix& B, const Matrix& Displacements, Vector& StrainVector )
    {
        KRATOS_TRY
        unsigned int Dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = Dim * (Dim + 1) / 2;
        noalias( StrainVector ) = ZeroVector( strain_size );

        for ( unsigned int node = 0; node < GetGeometry().size(); ++node )
        {
            for ( unsigned int item = 0; item < strain_size; ++item )
                for ( unsigned int dim = 0; dim < Dim; ++dim )
                    StrainVector[item] += B( item, Dim * node + dim ) * ( Displacements( node, dim ) - mInitialDisp( node, dim ) );
        }

        KRATOS_CATCH( "" )
    }

    /**
     * Computes the B-Operator at the current quadrature point
     * @param B_Operator current B-operator
     * @param DN_DX shape function values at the current integration point
     */
    void KinematicLinear::CalculateBoperator( Matrix& B_Operator, const Matrix& DN_DX )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().PointsNumber();

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;

        noalias( B_Operator ) = ZeroMatrix( strain_size, number_of_nodes * dim );

        if(dim == 2)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
                B_Operator( 0, i*2 ) = DN_DX( i, 0 );
                B_Operator( 1, i*2 + 1 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 ) = DN_DX( i, 1 );
                B_Operator( 2, i*2 + 1 ) = DN_DX( i, 0 );
            }
        }
        else if(dim == 3)
        {
            for ( unsigned int i = 0; i < number_of_nodes; ++i )
            {
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
        }

        KRATOS_CATCH( "" )
    }

    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == PK2_STRESS_TENSOR || rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
        {
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }

        return;
    }

    /**
     * Calculate Matrix Variables at each integration point, used for postprocessing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the qudrature points, output of the method
     * @param rCurrentProcessInfo
     *
     * Calculate Vector Variables at each integration point, used for postprocessing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the qudrature points, output of the method
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != mConstitutiveLawVector.size() )
            rValues.resize( mConstitutiveLawVector.size() );

        unsigned int dim = GetGeometry().WorkingSpaceDimension();
        unsigned int strain_size = dim * (dim + 1) / 2;
        unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = dim * number_of_nodes;

        if ( rVariable == MATERIAL_PARAMETERS )
        {
            if ( rValues.size() != GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() )
                rValues.resize( GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size() );

            for ( unsigned int ii = 0; ii < mConstitutiveLawVector.size(); ++ii )
                rValues[ii] = mConstitutiveLawVector[ii]->GetValue( rVariable, rValues[ii] );
        }

        if ( rVariable == INSITU_STRESS || rVariable == PRESTRESS )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                if ( rValues[i].size() != strain_size )
                    rValues[i].resize( strain_size );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PRESTRESS, rValues[i] );
            }
        }

        if ( rVariable == PLASTIC_STRAIN_VECTOR )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                if ( rValues[i].size() != strain_size )
                    rValues[i].resize( strain_size );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( PLASTIC_STRAIN_VECTOR, rValues[i] );
            }
        }

        if ( rVariable == STRESSES )
        {
            Vector tmp(strain_size);
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                if ( rValues[i].size() != strain_size )
                    rValues[i].resize( strain_size );

                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( STRESSES, rValues[i] );
            }
        }

        if ( rVariable == RECOVERY_STRESSES )
        {
            /////////////////////////////////////////////////////////////////////////
            //// Calculate recover stresses
            /////////////////////////////////////////////////////////////////////////
            if( GetProperties().Has( STRESS_RECOVERY_TYPE ) == true )
            {
                int Type = GetProperties()[ STRESS_RECOVERY_TYPE ];
            
                if(Type == 0)
                {
                    // no recovery
                    GetValueOnIntegrationPoints( STRESSES, rValues, rCurrentProcessInfo );
                }
                else if(Type == 1)
                {
                    // new recovery method from Bathe
                    
                    int ExpansionLevel = GetProperties()[ NEIGHBOUR_EXPANSION_LEVEL ];
                    
                    BatheRecoverStressUtility StressUtils(ExpansionLevel);
                    StressUtils.CalculateImprovedStressOnIntegrationPoints( *this, rValues, rCurrentProcessInfo );
                }
                else
                    KRATOS_THROW_ERROR(std::logic_error, "The stress recovery type is not supported on element", Id());
                
            }
            else
                KRATOS_THROW_ERROR(std::logic_error, "The stress recovery method is not defined for element", Id());
            
        }
        
        if ( rVariable == INTERNAL_VARIABLES )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                noalias( rValues[i] ) = mConstitutiveLawVector[i]->GetValue( INTERNAL_VARIABLES, rValues[i] );
            }
        }

        if ( rVariable == STRAIN || rVariable == CURRENT_STRAIN_VECTOR )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            //initialize the geometry
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            // calculate shape function values and local gradients
            Matrix B(strain_size, mat_size);
            Vector StrainVector(strain_size);
            Matrix DN_DX(number_of_nodes, dim);
            Matrix CurrentDisp(number_of_nodes, dim);

            const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients( mThisIntegrationMethod );
            //const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod );

            // extract current displacements
            for (unsigned int node = 0; node < GetGeometry().size(); ++node)
                noalias(row(CurrentDisp, node)) =
                    GetGeometry()[node].GetSolutionStepValue(DISPLACEMENT);

            for (unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            {
                if (rValues[i].size() != 6)
                    rValues[i].resize(6);

                // compute B_Operator at the current integration point
                noalias(DN_DX) = prod(DN_De[i], mInvJ0[i]);
                CalculateBoperator(B, DN_DX);

                // compute the strain at integration point
                CalculateStrain(B, CurrentDisp, StrainVector);
                if(dim == 2)
                {
                    rValues[i](0) = StrainVector(0);
                    rValues[i](1) = StrainVector(1);
                    rValues[i](2) = 0.0; // note: it's only correct for plane strain, TODO: we must make this available for plane stress constitutive law
                    rValues[i](3) = StrainVector(2);
                    rValues[i](4) = 0.0;
                    rValues[i](5) = 0.0;
                }
                else if(dim == 3)
                    noalias(rValues[i]) = StrainVector;
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif
        }
    }

    /**
     * Calculate double Variables at each integration point, used for postprocessing etc.
     * @param rVariable Global name of the variable to be calculated
     * @param rValues Vector to store the values on the qudrature points, output of the method
     * @param rCurrentProcessInfo
     */
    void KinematicLinear::GetValueOnIntegrationPoints( const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if ( rValues.size() != GetGeometry().IntegrationPoints().size() )
            rValues.resize( GetGeometry().IntegrationPoints().size() );

        if( rVariable == K0 )
        {
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
            {
                rValues[Point] = GetValue( K0 );
            }
            return;
        }
        else if( rVariable == STRAIN_ENERGY )
        {
            std::vector<Vector> StrainList(rValues.size());
            GetValueOnIntegrationPoints(STRAIN, StrainList, rCurrentProcessInfo);

            std::vector<Vector> StressList(rValues.size());
            GetValueOnIntegrationPoints(STRESSES, StressList, rCurrentProcessInfo);

            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                // calculate strain energy as C = 0.5 * <epsilon, sigma>
                rValues[i] = 0.5 * inner_prod(StrainList[i], StressList[i]);
            }
        }
        else if( rVariable == JACOBIAN_0 )
        {
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = mDetJ0[i];
            }
        }
        else if( rVariable == MATERIAL_DENSITY )
        {
            for( unsigned int i = 0; i < rValues.size(); ++i )
            {
                rValues[i] = this->GetValue(MATERIAL_DENSITY);
            }
        }
        else
        {
            //reading integration points and local gradients
            for ( unsigned int Point = 0; Point < mConstitutiveLawVector.size(); Point++ )
            {
                rValues[Point] = mConstitutiveLawVector[Point]->GetValue( rVariable, rValues[Point] );
            }
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
            std::cout << "In KinematicLinear:SetValueOnIntegrationPoints(...) Line " << __LINE__ << ", wrong size of the input vector, is = " << rValues.size() << " ; should be = " << mConstitutiveLawVector.size() << std::endl;
            return;
        }

        if ( rVariable == ELASTIC_LEFT_CAUCHY_GREEN_OLD )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
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
        if ( rVariable == INSITU_STRESS )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); ++PointNumber )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( INSITU_STRESS, rValues[PointNumber], rCurrentProcessInfo );
            }
        }

        if ( rVariable == PRESTRESS )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); ++PointNumber )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( PRESTRESS, rValues[PointNumber], rCurrentProcessInfo );
            }
        }

        if ( rVariable == MATERIAL_PARAMETERS )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); ++PointNumber )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( MATERIAL_PARAMETERS, rValues[PointNumber], rCurrentProcessInfo );
            }
        }

        if ( rVariable == STRESSES )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); ++PointNumber )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( STRESSES, rValues[PointNumber], rCurrentProcessInfo );
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
        if ( rVariable == SUCTION )
        {
            for ( unsigned int PointNumber = 0; PointNumber < GetGeometry().IntegrationPoints( mThisIntegrationMethod ).size(); ++PointNumber )
            {
                mConstitutiveLawVector[PointNumber]->SetValue( SUCTION, rValues[PointNumber],
                        rCurrentProcessInfo );
            }
        }
        else if ( rVariable == K0 )
        {
            SetValue( K0, rValues[0] );
        }
        else if ( rVariable == PRESTRESS_FACTOR )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->SetValue( PRESTRESS_FACTOR, rValues[i], rCurrentProcessInfo );
            }
        }
        else if ( rVariable == OVERCONSOLIDATION_RATIO )
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->SetValue( OVERCONSOLIDATION_RATIO, rValues[i], rCurrentProcessInfo );
            }
        }
        else
        {
            for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
            {
                mConstitutiveLawVector[i]->SetValue( rVariable, rValues[i], rCurrentProcessInfo );
            }
        }
    }

    void KinematicLinear::SetValueOnIntegrationPoints( const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo )
    {
        if(rVariable == PARENT_ELEMENT_ID)
        {
            rValues.resize(mConstitutiveLawVector.size());
            std::fill(rValues.begin(), rValues.end(), Id());
        }
        
        if(rVariable == INTEGRATION_POINT_INDEX)
        {
            rValues.resize(mConstitutiveLawVector.size());
            for(unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i)
            {
                rValues[i] = i;
            }
        }
    }

    void KinematicLinear::SetValueOnIntegrationPoints( const Kratos::Variable< ConstitutiveLaw::Pointer >& rVariable, std::vector< ConstitutiveLaw::Pointer >& rValues, const Kratos::ProcessInfo& rCurrentProcessInfo )
    {
        if ( rVariable == CONSTITUTIVE_LAW )
        {
            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Initialize(mThisIntegrationMethod);
            #endif

            for ( unsigned int i = 0; i < rValues.size(); ++i )
            {
                mConstitutiveLawVector[i] = rValues[i];
                mConstitutiveLawVector[i]->InitializeMaterial( GetProperties(), GetGeometry(), row( GetGeometry().ShapeFunctionsValues( mThisIntegrationMethod ), i ) );
            }

            #ifdef ENABLE_BEZIER_GEOMETRY
            GetGeometry().Clean();
            #endif
        }
    }

    int KinematicLinear::Check( const Kratos::ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

//	    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

        if ( this->Id() < 1 )
        {
            KRATOS_THROW_ERROR( std::logic_error, "Element found with Id 0 or negative, Id() =", Id() );
        }

        if ( mTotalDomainInitialSize < 0.0 )
        {
            std::cout << "error on element -> " << this->Id() << std::endl;
            KRATOS_THROW_ERROR( std::logic_error, "Domain size can not be less than 0, mTotalDomainInitialSize =", mTotalDomainInitialSize );
        }

        //verify that the constitutive law exists
        if ( this->GetProperties().Has( CONSTITUTIVE_LAW ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "constitutive law not provided for property", this->GetProperties().Id() );
        }

        //Verify that the body force is defined
        if ( this->GetProperties().Has( BODY_FORCE ) == false )
        {
            KRATOS_THROW_ERROR( std::logic_error, "BODY_FORCE not provided for property", this->GetProperties().Id() )
        }

        //verify that the constitutive law has the correct dimension
//        if ( dimension == 2 )
//        {
//            if ( this->GetProperties().Has( THICKNESS ) == false )
//                KRATOS_THROW_ERROR( std::logic_error, "THICKNESS not provided for element ", this->Id() );

//            if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Getstrain_size() != 3 )
//                KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 2D element! expected strain size is 3 (el id = ) ", this->Id() );
//        }
//        else if(dimension == 3)
//        {
//            if ( this->GetProperties().GetValue( CONSTITUTIVE_LAW )->Getstrain_size() != 6 )
//                KRATOS_THROW_ERROR( std::logic_error, "wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) ", this->Id() );
//        }

        //check constitutive law
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            return mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
        }

        int ok = 0;
        for ( unsigned int i = 0; i < mConstitutiveLawVector.size(); ++i )
        {
            ok = mConstitutiveLawVector[i]->Check( GetProperties(), GetGeometry(), rCurrentProcessInfo );
            if( ok != 0 ) break;
//            if( mConstitutiveLawVector[i]->IsIncremental() )
//                KRATOS_THROW_ERROR( std::logic_error, "This element does not provide incremental strains!", "" );
//            if( mConstitutiveLawVector[i]->GetStrainMeasure() != ConstitutiveLaw::StrainMeasure_Linear )
//                KRATOS_THROW_ERROR( std::logic_error, "This element formulated in linear strain measure", "" );
//            if( mConstitutiveLawVector[i]->GetStressMeasure() != ConstitutiveLaw::StressMeasure_PK1 )
//                KRATOS_THROW_ERROR( std::logic_error, "This element is formulated in PK1 stresses", "" );
        }

        return ok;

        KRATOS_CATCH( "" );
    }

} // Namespace Kratos

#undef ENABLE_DEBUG_CONSTITUTIVE_LAW

