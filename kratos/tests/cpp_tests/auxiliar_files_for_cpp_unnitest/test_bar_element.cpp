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
#include <math.h>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_bar_element.h"

namespace Kratos
{
    namespace Testing
    {
    //***********************DEFAULT CONSTRUCTOR******************************************//
    //************************************************************************************//

    TestBarElement::TestBarElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        )
        : Element( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }


    //******************************CONSTRUCTOR*******************************************//
    //************************************************************************************//

    TestBarElement::TestBarElement(
        IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        )
        : Element( NewId, pGeometry, pProperties )
    {

    }

    //******************************COPY CONSTRUCTOR**************************************//
    //************************************************************************************//

    TestBarElement::TestBarElement( TestBarElement const& rOther)
        :Element(rOther)
    {

    }

    //*******************************ASSIGMENT OPERATOR***********************************//
    //************************************************************************************//

    TestBarElement&  TestBarElement::operator=(TestBarElement const& rOther)
    {
        //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

        Element::operator=(rOther);

        return *this;
    }

    //*********************************OPERATIONS*****************************************//
    //************************************************************************************//

    Element::Pointer TestBarElement::Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties
        ) const
    {
        //NEEDED TO CREATE AN ELEMENT
        return Kratos::make_intrusive<TestBarElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties );
    }


    //************************************CLONE*******************************************//
    //************************************************************************************//

    Element::Pointer TestBarElement::Clone(
        IndexType NewId,
        NodesArrayType const& rThisNodes
        ) const
    {
        //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
        //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

        TestBarElement new_element(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

        return Kratos::make_intrusive<TestBarElement>(new_element);
    }


    //*******************************DESTRUCTOR*******************************************//
    //************************************************************************************//

    TestBarElement::~TestBarElement()
    {
    }

    //************* GETTING METHODS
    //************************************************************************************//
    //************************************************************************************//

    void TestBarElement::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        //NEEDED TO DEFINE THE DOFS OF THE ELEMENT

        rElementalDofList.resize( 0 );

        rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[0].pGetDof( DISPLACEMENT_Y ) );
        rElementalDofList.push_back( GetGeometry()[1].pGetDof( DISPLACEMENT_X ) );
        rElementalDofList.push_back( GetGeometry()[1].pGetDof( DISPLACEMENT_Y ) );
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestBarElement::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        //NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY

        if ( rResult.size() != 4 )
            rResult.resize( 4, false );

        rResult[0] = GetGeometry()[0].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[1] = GetGeometry()[0].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[2] = GetGeometry()[1].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[3] = GetGeometry()[1].GetDof( DISPLACEMENT_Y ).EquationId();
    }

    //*********************************DISPLACEMENT***************************************//
    //************************************************************************************//

    void TestBarElement::GetValuesVector( Vector& rValues, int Step )
    {
        //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
        if ( rValues.size() != 4 )
            rValues.resize( 4, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( DISPLACEMENT_Y, Step );
        rValues[2] = GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[3] = GetGeometry()[1].GetSolutionStepValue( DISPLACEMENT_Y, Step );
    }


    //************************************VELOCITY****************************************//
    //************************************************************************************//

    void TestBarElement::GetFirstDerivativesVector( Vector& rValues, int Step )
    {
        //GIVES THE VECTOR WITH THE TIME DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT VELOCITIES)
        if ( rValues.size() != 4 )
            rValues.resize( 4, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( VELOCITY_Y, Step );
        rValues[2] = GetGeometry()[1].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[3] = GetGeometry()[1].GetSolutionStepValue( VELOCITY_Y, Step );
    }

    //*********************************ACCELERATION***************************************//
    //************************************************************************************//

    void TestBarElement::GetSecondDerivativesVector( Vector& rValues, int Step )
    {
        //GIVES THE VECTOR WITH THE TIME SECOND DERIVATIVE OF THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT ACCELERATIONS)
        if ( rValues.size() != 4 )
            rValues.resize( 4, false );

        rValues[0] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[1] = GetGeometry()[0].GetSolutionStepValue( ACCELERATION_Y, Step );
        rValues[2] = GetGeometry()[1].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[3] = GetGeometry()[1].GetSolutionStepValue( ACCELERATION_Y, Step );
    }

    //************* COMPUTING  METHODS
    //************************************************************************************//
    //************************************************************************************//

    void TestBarElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
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

    void TestBarElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
    {
        // Resizing as needed the RHS
        const unsigned int system_size = 4;

        if ( rRightHandSideVector.size() != system_size )
            rRightHandSideVector.resize( system_size, false );

        rRightHandSideVector = ZeroVector( system_size ); //resetting RHS

        // We don't care about filling the RHS, we are just interested in the LHS for the test
    }

    //***********************************************************************************
    //***********************************************************************************

    void TestBarElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        // Resizing as needed the LHS
        const unsigned int system_size = 4;

        if ( rLeftHandSideMatrix.size1() != system_size )
            rLeftHandSideMatrix.resize( system_size, system_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( system_size, system_size ); //resetting LHS

        const double E = this->GetProperties()[YOUNG_MODULUS];
        const double A = this->GetProperties()[NODAL_AREA];
        const double L = this->GetGeometry().Length();
        const double k = E*A/L;

        const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
        const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
        const double theta = std::atan2(dy, dx);
        const double l = std::cos(theta);
        const double m = std::sin(theta);

        rLeftHandSideMatrix(0,0) = k * l * l;
        rLeftHandSideMatrix(0,1) = k * m * l;
        rLeftHandSideMatrix(0,2) = -k * l * l;
        rLeftHandSideMatrix(0,3) = -k * m * l;

        rLeftHandSideMatrix(1,0) = k * l * m;
        rLeftHandSideMatrix(1,1) = k * m * m;
        rLeftHandSideMatrix(1,2) = -k * m * l;
        rLeftHandSideMatrix(1,3) = -k * m * m;

        rLeftHandSideMatrix(2,0) = -k * l * l;
        rLeftHandSideMatrix(2,1) = -k * m * l;
        rLeftHandSideMatrix(2,2) = k * l * l;
        rLeftHandSideMatrix(2,3) = k * m * l;

        rLeftHandSideMatrix(3,0) = -k * l * m;
        rLeftHandSideMatrix(3,1) = -k * m * m;
        rLeftHandSideMatrix(3,2) = k * m * l;
        rLeftHandSideMatrix(3,3) = k * m * m;
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestBarElement::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        unsigned int system_size = 4;
        if ( rMassMatrix.size1() != system_size || rMassMatrix.size2() != system_size) {
            rMassMatrix.resize( system_size, system_size, false );
        }

        rMassMatrix = ZeroMatrix( system_size, system_size );
        const double rho = this->GetProperties()[DENSITY];
        const double A = this->GetProperties()[NODAL_AREA];
        const double L = this->GetGeometry().Length();
        const double aux_val = A * rho * L / 6.0;

        const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
        const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
        const double l = dx / L;
        const double m = dy / L;

        rMassMatrix(0,0) = 2.0 * aux_val * l * l;
        rMassMatrix(0,1) = 2.0 * aux_val * l * m;
        rMassMatrix(0,2) = aux_val * l * l;
        rMassMatrix(0,3) = aux_val * l * m;

        rMassMatrix(1,0) = 2.0 * aux_val * l * m;
        rMassMatrix(1,1) = 2.0 * aux_val * m * m;
        rMassMatrix(1,2) = aux_val * l * m;
        rMassMatrix(1,3) = aux_val * m * m;

        rMassMatrix(2,0) = aux_val * l * l;
        rMassMatrix(2,1) = aux_val * l * m;
        rMassMatrix(2,2) = 2.0 * aux_val * l * l;
        rMassMatrix(2,3) = 2.0 * aux_val * l * m;

        rMassMatrix(3,0) = aux_val * l * m;
        rMassMatrix(3,1) = aux_val * m * m;
        rMassMatrix(3,2) = 2.0 * aux_val * l * m;
        rMassMatrix(3,3) = 2.0 * aux_val * m * m;

        KRATOS_CATCH( "" );
    }

    //************************************************************************************//
    //************************************************************************************//

    void TestBarElement::CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        // Resizing as needed the LHS
        const unsigned int system_size = 4;

        rDampingMatrix = ZeroMatrix( system_size, system_size );

        KRATOS_CATCH( "" );
    }

    //************************************************************************************//
    //************************************************************************************//

    int TestBarElement::Check( const ProcessInfo& rCurrentProcessInfo )
    {
        KRATOS_TRY

        // Check that all required variables have been registered
        KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
        KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
        KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA)
        KRATOS_CHECK_VARIABLE_KEY(YOUNG_MODULUS)

        KRATOS_ERROR_IF_NOT(GetProperties().Has(YOUNG_MODULUS)) << "YOUNG_MODULUS not defined" << std::endl;
        KRATOS_ERROR_IF_NOT(GetProperties().Has(NODAL_AREA)) << "NODAL_AREA not defined" << std::endl;

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

        KRATOS_CATCH( "Problem in the Check in the TestBarElement" )
    }


    //************************************************************************************//
    //************************************************************************************//

    void TestBarElement::save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
    }

    void TestBarElement::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
    }

    } // Namespace Testing
} // Namespace Kratos
