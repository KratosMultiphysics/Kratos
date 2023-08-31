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
#include "includes/expect.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "tests/test_utilities/test_laplacian_element.h"

namespace Kratos::Testing
{
//***********************DEFAULT CONSTRUCTOR******************************************//
//************************************************************************************//

TestLaplacianElement::TestLaplacianElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//******************************CONSTRUCTOR*******************************************//
//************************************************************************************//

TestLaplacianElement::TestLaplacianElement(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
    : Element( NewId, pGeometry, pProperties )
{
}

//******************************COPY CONSTRUCTOR**************************************//
//************************************************************************************//

TestLaplacianElement::TestLaplacianElement( TestLaplacianElement const& rOther)
    :Element(rOther)
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
    return Kratos::make_intrusive<TestLaplacianElement>( NewId, GetGeometry().Create( rThisNodes ), pProperties );
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

    TestLaplacianElement new_element(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

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
    if (rElementalDofList.size() != r_geometry.size()) {
        rElementalDofList.resize( r_geometry.size());
    }
    unsigned int counter = 0;
    for (auto& r_node : r_geometry) {
        rElementalDofList[counter] = r_node.pGetDof(TEMPERATURE);
        ++counter;
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

//*******************************COMPUTING  METHODS***********************************//
//************************************************************************************//

void TestLaplacianElement::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    // Geometry definition
    auto& r_geometry = this->GetGeometry();
    const unsigned int number_of_nodes = r_geometry.size();
    const unsigned int dimension = r_geometry.WorkingSpaceDimension();

    // Some definitions
    MatrixType DN_DX = ZeroMatrix(number_of_nodes, dimension);  // Gradients matrix
    MatrixType D = ZeroMatrix(dimension, dimension);            // Conductivity matrix
    VectorType N = ZeroVector(number_of_nodes);                 // Size = number of nodes . Position of the gauss point
    VectorType temp = ZeroVector(number_of_nodes);              // Dimension = number of nodes . . since we are using a residualbased approach

    if(rLeftHandSideMatrix.size1() != number_of_nodes)
        rLeftHandSideMatrix.resize(number_of_nodes,number_of_nodes,false); //resizing the system in case it does not have the right size
    rLeftHandSideMatrix.clear();

    if(rRightHandSideVector.size() != number_of_nodes)
        rRightHandSideVector.resize(number_of_nodes,false);

    // Reading properties and conditions
    const double conductivity = GetProperties()[CONDUCTIVITY];
    for (unsigned int i_dim = 0; i_dim < dimension; ++i_dim) {
        D(i_dim, i_dim) = conductivity;
    }

    // Gauss point loop
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
    double detJ = 0.0;
    MatrixType J, invJ, DN_De;
    for (unsigned int i_gp = 0; i_gp < r_integration_points.size(); ++i_gp) {
        auto& r_gp = r_integration_points[i_gp];
        GeometryUtils::JacobianOnInitialConfiguration(r_geometry, r_gp, J);
        MathUtils<double>::InvertMatrix(J, invJ, detJ);
        r_geometry.ShapeFunctionsLocalGradients(DN_De, r_gp);
        GeometryUtils::ShapeFunctionsGradients(DN_De, invJ, DN_DX);
        N = r_geometry.ShapeFunctionsValues(N, r_gp.Coordinates());
        noalias(rLeftHandSideMatrix) += r_gp.Weight() * detJ * prod(DN_DX, Matrix(prod(D, trans(DN_DX))));  // w detJ Bt D B
    }

    // Subtracting the dirichlet term
    // RHS -= LHS*DUMMY_UNKNOWNs
    for(unsigned int i_point = 0; i_point < number_of_nodes; ++i_point) {
        temp[i_point] = r_geometry[i_point].FastGetSolutionStepValue(TEMPERATURE);
    }
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, temp);

    KRATOS_CATCH( "" );
}

//***********************************************************************************
//***********************************************************************************

void TestLaplacianElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo )
{
    MatrixType lhs;
    CalculateLocalSystem(lhs, rRightHandSideVector, rCurrentProcessInfo);
}

//***********************************************************************************
//***********************************************************************************

void TestLaplacianElement::CalculateLeftHandSide( MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    VectorType rhs;
    CalculateLocalSystem(rLeftHandSideMatrix, rhs, rCurrentProcessInfo);
}

//************************************************************************************//
//************************************************************************************//

void TestLaplacianElement::CalculateMassMatrix( MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const unsigned int system_size = GetGeometry().size();

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

    const unsigned int system_size = GetGeometry().size();

    if ( rDampingMatrix.size1() != system_size )
        rDampingMatrix.resize( system_size, system_size, false );

    rDampingMatrix = ZeroMatrix( system_size, system_size );

    KRATOS_CATCH( "" );
}

//************************************************************************************//
//************************************************************************************//

int TestLaplacianElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (auto& r_node : this->GetGeometry()) {
        KRATOS_EXPECT_VARIABLE_IN_NODAL_DATA(TEMPERATURE,r_node)
        KRATOS_EXPECT_DOF_IN_NODE(TEMPERATURE,r_node)
    }

    return 0;

    KRATOS_CATCH( "Problem in the Check in the TestLaplacianElement" )
}

//************************************************************************************//
//************************************************************************************//

void TestLaplacianElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // ADD SOMETHING IF REQUIRED
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


