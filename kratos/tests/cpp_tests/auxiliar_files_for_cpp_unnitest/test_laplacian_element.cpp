//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/variables.h"
#include "tests/cpp_tests/auxiliar_files_for_cpp_unnitest/test_laplacian_element.h"

namespace Kratos::Testing
{
//***********************DEFAULT CONSTRUCTOR******************************************//
//************************************************************************************//

TestLaplacianElement::TestLaplacianElement(
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

TestLaplacianElement::TestLaplacianElement(
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

TestLaplacianElement::TestLaplacianElement( TestLaplacianElement const& rOther)
    :Element(rOther)
    ,mResidualType(rOther.mResidualType)
{

}

//*******************************ASSIGMENT OPERATOR***********************************//
//************************************************************************************//

TestLaplacianElement&  TestLaplacianElement::operator=(TestLaplacianElement const& rOther)
{
    //ALL MEMBER VARIABLES THAT MUST BE KEPT IN AN "=" OPERATION NEEDS TO BE COPIED HERE

    Element::operator=(rOther);

    return *this;
}

//*********************************OPERATIONS*****************************************//
//************************************************************************************//

Element::Pointer TestLaplacianElement::Create(
    IndexType NewId,
    NodesArrayType const& rThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    //NEEDED TO CREATE AN ELEMENT
    return Kratos::make_intrusive<TestLaplacianElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties, mResidualType );
}


//************************************CLONE*******************************************//
//************************************************************************************//

Element::Pointer TestLaplacianElement::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    //YOU CREATE A NEW ELEMENT CLONING THEIR VARIABLES
    //ALL MEMBER VARIABLES THAT MUST BE CLONED HAVE TO BE DEFINED HERE

    TestLaplacianElement new_element(NewId, GetGeometry().Create( rThisNodes ), pGetProperties(), mResidualType );

    return Kratos::make_intrusive<TestLaplacianElement>(new_element);
}

//*******************************DESTRUCTOR*******************************************//
//************************************************************************************//

TestLaplacianElement::~TestLaplacianElement()
{
}

//******************************GETTING METHODS***************************************//
//************************************************************************************//

void TestLaplacianElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // NEEDED TO DEFINE THE DOFS OF THE ELEMENT
    const auto& r_geometry = this->GetGeometry();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(r_geometry.size());
    for (auto& r_node : r_geometry) {
        rElementalDofList.push_back(r_node.pGetDof(TEMPERATURE));
    }
}

//************************************************************************************//
//************************************************************************************//

void TestLaplacianElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // NEEDED TO DEFINE GLOBAL IDS FOR THE CORRECT ASSEMBLY
    const auto& r_geometry = this->GetGeometry();
    if (rResult.size() != r_geometry.size()) {
        rResult.resize( r_geometry.size());
    }

    unsigned int counter = 0;
    for (auto& r_node : r_geometry) {
        rResult[counter] = r_node.GetDof(TEMPERATURE).EquationId();
        ++counter;
    }
}

//*********************************TEMPERATURE****************************************//
//************************************************************************************//

void TestLaplacianElement::GetValuesVector( Vector& rValues, int Step ) const
{
    //GIVES THE VECTOR WITH THE DOFS VARIABLES OF THE ELEMENT (i.e. ELEMENT DISPLACEMENTS)
    const auto& r_geometry = this->GetGeometry();
    if (rValues.size() != r_geometry.size()) {
        rValues.resize( r_geometry.size());
    }

    unsigned int counter = 0;
    for (auto& r_node : r_geometry) {
        rValues[counter] = r_node.GetSolutionStepValue(TEMPERATURE, Step);
        ++counter;
    }
}

//************* COMPUTING  METHODS
//************************************************************************************//
//************************************************************************************//

void TestLaplacianElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{

    KRATOS_TRY;

    /* Calculate elemental system */

    // Compute RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // Compute LHS
    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

//***********************************************************************************
//***********************************************************************************

void TestLaplacianElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
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

void TestLaplacianElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
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

void TestLaplacianElement::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
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

void TestLaplacianElement::CalculateDampingMatrix(
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

int TestLaplacianElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( std::size_t i = 0; i < this->GetGeometry().size(); ++i ) {
        const Node& rnode = this->GetGeometry()[i];

        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,rnode)
    }

    return 0;

    KRATOS_CATCH( "Problem in the Check in the TestLaplacianElement" )
}


//************************************************************************************//
//************************************************************************************//

void TestLaplacianElement::CalculateOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
if (rVariable == CONSTITUTIVE_LAW) {
    const SizeType integration_points_number = mConstitutiveLawVector.size();
    if (rValues.size() != integration_points_number) {
        rValues.resize(integration_points_number);
    }
    for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
        rValues[point_number] = mConstitutiveLawVector[point_number];
    }
}
}

//************************************************************************************//
//************************************************************************************//

void TestLaplacianElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    const GeometryType& r_geometry = GetGeometry();
    const Properties& r_properties = GetProperties();
    const auto& N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLawVector[point_number]->InitializeMaterial( r_properties, r_geometry, row(N_values , point_number ));
    }
}

//************************************************************************************//
//************************************************************************************//

void TestLaplacianElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element )
}

//************************************************************************************//
//************************************************************************************//

void TestLaplacianElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element )
}

} // Namespace Kratos::Testing


