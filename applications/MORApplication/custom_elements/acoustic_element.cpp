// KRATOS
//
//  License:		BSD License
//					license: ../../license.txt
//
//  Main authors:   Elias Schwaiger
//                  Quirin Aumann
//
//

// System includes

// External includes


// Project includes
#include "utilities/geometry_utilities.h"
// Application includes
#include "custom_elements/acoustic_element.h"
#include "mor_application_variables.h"

namespace Kratos
{
AcousticElement::AcousticElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : Element( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

AcousticElement::AcousticElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Element( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<AcousticElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<AcousticElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

AcousticElement::~AcousticElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer AcousticElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    AcousticElement::Pointer p_new_elem = Kratos::make_intrusive<AcousticElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int AcousticElement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY
    const double numerical_limit = std::numeric_limits<double>::epsilon();

    const int ier = Element::Check(rCurrentProcessInfo);
    if( ier != 0 )
        return ier;

    // Check that all required variables have resonable values
    KRATOS_ERROR_IF(!GetProperties().Has(BULK_MODULUS) ||
                GetProperties()[BULK_MODULUS] <= numerical_limit)
        << "Please provide a reasonable value for \"BULK_MODULUS\" for element #"
        << Id() << std::endl;

    KRATOS_ERROR_IF(!GetProperties().Has(DENSITY) ||
                GetProperties()[DENSITY] <= numerical_limit)
        << "Please provide a reasonable value for \"DENSITY\" for element #"
        << Id() << std::endl;

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const SizeType number_of_points = GetGeometry().size();
    for( IndexType i = 0; i < number_of_points; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, rnode)
        KRATOS_CHECK_DOF_IN_NODE(PRESSURE, rnode)
    }

    return ier;

    KRATOS_CATCH( "" );
}

void AcousticElement::Initialize()
{
    KRATOS_TRY
    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::EquationIdVector(EquationIdVectorType& rResult,
                                        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_nodes = GetGeometry().PointsNumber();

    if( rResult.size() != num_nodes )
        rResult.resize(num_nodes,false);

    for( SizeType i_node = 0; i_node < num_nodes; ++i_node )
        rResult[i_node] = GetGeometry()[i_node].GetDof(PRESSURE).EquationId();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY;

    const SizeType num_nodes = GetGeometry().PointsNumber();

    if(rElementalDofList.size() != num_nodes)
        rElementalDofList.resize(num_nodes);

    for (SizeType i_node = 0; i_node < num_nodes; i_node++)
        rElementalDofList[i_node] = GetGeometry()[i_node].pGetDof(PRESSURE);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{

    const GeometryType& geom = GetGeometry();
    IntegrationMethod ThisIntegrationMethod = geom.GetDefaultIntegrationMethod();
    const SizeType number_of_nodes = geom.PointsNumber();
    const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(ThisIntegrationMethod);

    if( rLeftHandSideMatrix.size1() != number_of_nodes || rLeftHandSideMatrix.size2() != number_of_nodes ) {
        rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );

    ShapeFunctionDerivativesArrayType DN_DX;
    Vector DetJ;

    DN_DX = geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, ThisIntegrationMethod);
    const double density = GetProperties()[DENSITY];

    for( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        double int_weight = integration_points[point_number].Weight() * DetJ(point_number);
        noalias( rLeftHandSideMatrix ) += int_weight/density * prod( DN_DX[point_number], trans(DN_DX[point_number]) );
    }

}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& geom = GetGeometry();
    IntegrationMethod ThisIntegrationMethod = geom.GetDefaultIntegrationMethod();
    const SizeType number_of_nodes = geom.PointsNumber();
    const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(ThisIntegrationMethod);
    IndexType NumGauss = integration_points.size();
    const Matrix& NContainer = geom.ShapeFunctionsValues(ThisIntegrationMethod);

    if( rMassMatrix.size1() != number_of_nodes || rMassMatrix.size2() != number_of_nodes ) {
        rMassMatrix.resize(number_of_nodes, number_of_nodes, false);
    }
    noalias(rMassMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );

    const double bulk = GetProperties()[BULK_MODULUS];

    for( IndexType i = 0; i < number_of_nodes; ++i ) {
        for( IndexType j = 0; j < number_of_nodes; ++j ) {
            for( IndexType g = 0; g < NumGauss; ++g ) {
                const double DetJ = geom.DeterminantOfJacobian(g, ThisIntegrationMethod);
                const double GaussWeight = DetJ * integration_points[g].Weight();
                rMassMatrix(i,j) += NContainer(g, i) * NContainer(g, j) * GaussWeight * (1/bulk);
            }
        }
    }

}


/***********************************************************************************/
void AcousticElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& geom = GetGeometry();
    const SizeType number_of_nodes = geom.PointsNumber();
    if( rDampingMatrix.size1() != number_of_nodes || rDampingMatrix.size2() != number_of_nodes ) {
        rDampingMatrix.resize(number_of_nodes, number_of_nodes, false);
    }
    noalias(rDampingMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );

}

/***********************************************************************************/

void AcousticElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

    const GeometryType& geom = GetGeometry();
    const SizeType number_of_nodes = geom.PointsNumber();

    if ( rRightHandSideVector.size() != number_of_nodes ) {
        rRightHandSideVector.resize( number_of_nodes, false );
    }

    rRightHandSideVector = ZeroVector( number_of_nodes ); //resetting RHS

}


/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    this->CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH( "" )

}


/***********************************************************************************/


void AcousticElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
}

} // Namespace Kratos

