//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "tests/test_element.h"

namespace Kratos
{
    namespace Testing
    {
    //***********************DEFAULT CONSTRUCTOR******************************************//
    //************************************************************************************//

    TestElement::TestElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        const ResidualType TheResidualType
        )
        : Element( NewId, pGeometry )
        , mResidualType( TheResidualType )
    {
        //DO NOT ADD DOFS HERE!!!
    }


    //******************************CONSTRUCTOR*******************************************//
    //************************************************************************************//

    TestElement::TestElement(
        IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        const ResidualType TheResidualType
        )
        : Element( NewId, pGeometry, pProperties )
        , mResidualType( TheResidualType )
    {

    }

    //******************************COPY CONSTRUCTOR**************************************//
    //************************************************************************************//

    TestElement::TestElement( TestElement const& rOther)
        :Element(rOther)
        ,mResidualType(rOther.mResidualType)
    {

    }

    //*******************************ASSIGMENT OPERATOR***********************************//
    //************************************************************************************//

    TestElement&  TestElement::operator=(TestElement const& rOther)
    {
        //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

        Element::operator=(rOther);

        return *this;
    }

    //*********************************OPERATIONS*****************************************//
    //************************************************************************************//

    Element::Pointer TestElement::Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties
        ) const
    {
        //NEEDED TO CREATE AN ELEMENT
        return Kratos::make_shared<TestElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties, mResidualType );
    }


    //************************************CLONE*******************************************//
    //************************************************************************************//

    Element::Pointer TestElement::Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const
    {
        //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
        //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

        TestElement new_element(NewId, GetGeometry().Create( rThisNodes ), pGetProperties(), mResidualType );

        return Kratos::make_shared<TestElement>(new_element);
    }


    //*******************************DESTRUCTOR*******************************************//
    //************************************************************************************//

    TestElement::~TestElement()
    {
    }

    //************* GETTING METHODS
    //************************************************************************************//
    //************************************************************************************//

    void TestElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        //NEEDED TO DEFINE THE DOFS OF THE ELEMENT
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        rElementalDofList.resize( 0 );

        rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
        if( dimension == 3 )
            rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Z ) );
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
        ) const
    {
        //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( rResult.size() != dimension )
            rResult.resize( dimension, false );

        rResult[0] = GetGeometry()[0].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[1] = GetGeometry()[0].GetDof( DISPLACEMENT_Y ).EquationId();
        if( dimension == 3)
            rResult[2] = GetGeometry()[0].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    //*********************************DISPLACEMENT***************************************//
    //************************************************************************************//

    void TestElement::GetValuesVector( Vector& rValues, int Step ) const
    {
        //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( rValues.size() != dimension )
            rValues.resize( dimension, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 )
            rValues[2] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Z, Step );
    }


    //************************************VELOCITY****************************************//
    //************************************************************************************//

    void TestElement::GetFirstDerivativesVector( Vector& rValues, int Step ) const
    {
        //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( rValues.size() != dimension )
            rValues.resize( dimension, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 )
            rValues[2] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Z, Step );
    }

    //*********************************ACCELERATION***************************************//
    //************************************************************************************//

    void TestElement::GetSecondDerivativesVector( Vector& rValues, int Step ) const
    {
        //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        if ( rValues.size() != dimension )
            rValues.resize( dimension, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 )
            rValues[2] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Z, Step );
    }

    //************* COMPUTING  METHODS
    //************************************************************************************//
    //************************************************************************************//

    void TestElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {

        KRATOS_TRY;

        /* Calculate elemental system */

        // Compute RHS (RHS = rRightHandSideVector = Fext - Fint)
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

        // Compute LHS
        this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    //***********************************************************************************
    //***********************************************************************************

    void TestElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the RHS
        const unsigned int system_size = dimension;

        if ( rRightHandSideVector.size() != system_size )
            rRightHandSideVector.resize( system_size, false );

        rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

        const array_1d<double, 3 >& delta_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);

        switch ( mResidualType )
        {
            case ResidualType::LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rRightHandSideVector[j] -= delta_displacement[j] - 1.0;
                break;
            case ResidualType::NON_LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rRightHandSideVector[j] -= std::pow(delta_displacement[j], 2) - 1.0;
                break;
            default:
                KRATOS_ERROR << "NOT IMPLEMENTED" << std::endl;
        }
    }

    //***********************************************************************************
    //***********************************************************************************

    void TestElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int system_size = dimension;

        if ( rLeftHandSideMatrix.size1() != system_size )
            rLeftHandSideMatrix.resize( system_size, system_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

        const array_1d<double, 3 >& delta_displacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);

        switch ( mResidualType )
        {
            case ResidualType::LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rLeftHandSideMatrix(j, j) += 1.0;
                break;
            case ResidualType::NON_LINEAR:
                for ( unsigned int j = 0; j < dimension; ++j )
                    rLeftHandSideMatrix(j, j) += delta_displacement[j] * 2;
                break;
            default:
                KRATOS_ERROR << "NOT IMPLEMENTED" << std::endl;
        }
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestElement::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        //lumped
        unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        unsigned int system_size = dimension;

        if ( rMassMatrix.size1() != system_size )
            rMassMatrix.resize( system_size, system_size, false );

        rMassMatrix = ZeroMatrix( system_size, system_size );

        KRATOS_CATCH( "" );
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;

        //0.-Initialize the DampingMatrix:
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

        // Resizing as needed the LHS
        const unsigned int system_size = dimension;

        rDampingMatrix = ZeroMatrix( system_size, system_size );

        KRATOS_CATCH( "" );
    }

    //************************************************************************************//
    //************************************************************************************//

    int TestElement::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for ( std::size_t i = 0; i < this->GetGeometry().size(); ++i ) {
            Node<3>& rnode = this->GetGeometry()[i];

            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rnode)
            KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rnode)
        }

        return 0;

        KRATOS_CATCH( "Problem in the Check in the TestElement" )
    }


    //************************************************************************************//
    //************************************************************************************//

    void TestElement::save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void TestElement::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    } // Namespace Testing
} // Namespace Kratos


