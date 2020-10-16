// KRATOS
//
//  License:		 BSD License
//					 license: ../../license.txt
//
//  Main authors:
//
//

// System includes

// External includes


// Project includes
#include "utilities/geometry_utilities.h"
// Application includes
#include "custom_elements/acoustic_element.h"

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

    // // The vector containing the constitutive laws
    // p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

int  AcousticElement::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    int ier = Element::Check(rCurrentProcessInfo);

    return ier;

    KRATOS_CATCH( "" );
}

void AcousticElement::Initialize()
{
    KRATOS_TRY
    std::cout << "i am initializing an acoustic element\n";
    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::EquationIdVector(EquationIdVectorType& rResult,
                                        ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    SizeType num_nodes = GetGeometry().PointsNumber();

    if(rResult.size() != num_nodes)
        rResult.resize(num_nodes,false);	

    for (SizeType i_node = 0; i_node < num_nodes; i_node++)
        rResult[i_node] = GetGeometry()[i_node].GetDof(PRESSURE).EquationId();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo)
{
    SizeType num_nodes = GetGeometry().PointsNumber();
  //  std::cout << "hello?\n";

    if(rElementalDofList.size() != num_nodes)
        rElementalDofList.resize(num_nodes);	

    for (SizeType i_node = 0; i_node < num_nodes; i_node++)
        rElementalDofList[i_node] = GetGeometry()[i_node].pGetDof(PRESSURE);

   // KRATOS_WATCH(rElementalDofList)
}

/***********************************************************************************/
/***********************************************************************************/

// double AcousticElement::CalculateDerivativesOnReferenceConfiguration(
//     Matrix& rJ0,
//     Matrix& rInvJ0,
//     Matrix& rDN_DX,
//     const IndexType PointNumber,
//     IntegrationMethod ThisIntegrationMethod
//     ) const
// {
//     const GeometryType& r_geom = GetGeometry();
//     GeometryUtils::JacobianOnInitialConfiguration(
//         r_geom,
//         r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ0);
//     double detJ0;
//     MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
//     const Matrix& rDN_De =
//         GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
//     GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
//     return detJ0;
// }

/***********************************************************************************/
/***********************************************************************************/

// void AcousticElement::CalculateB(
//     const Matrix& rDN_DX,
//     Matrix& rB
//     )
// {
//     const auto& r_geometry = GetGeometry();
//     const SizeType number_of_nodes = r_geometry.PointsNumber();
//     const SizeType dimension = r_geometry.WorkingSpaceDimension();

//     rB.clear();

//     if(dimension == 2) {
//         for ( IndexType i = 0; i < number_of_nodes; ++i ) {
//             rB(0, i    ) = rDN_DX(i, 0);
//             rB(1, i    ) = rDN_DX(i, 1);
//             rB(2, i    ) = rDN_DX(i, 2);
//         }
//     } else if(dimension == 3) {
//         for ( IndexType i = 0; i < number_of_nodes; ++i ) {
//             const IndexType initial_index = i*3;
//             rB(0, initial_index    ) = rDN_DX(i, 0);
//             rB(1, initial_index + 1) = rDN_DX(i, 1);
//             rB(2, initial_index + 2) = rDN_DX(i, 2);
//             rB(3, initial_index    ) = rDN_DX(i, 1);
//             rB(3, initial_index + 1) = rDN_DX(i, 0);
//             rB(4, initial_index + 1) = rDN_DX(i, 2);
//             rB(4, initial_index + 2) = rDN_DX(i, 1);
//             rB(5, initial_index    ) = rDN_DX(i, 2);
//             rB(5, initial_index + 2) = rDN_DX(i, 0);
//         }
//     }
// }

/***********************************************************************************/
/***********************************************************************************/

// void AcousticElement::CalculateKinematicVariables(
//     KinematicVariables& rThisKinematicVariables,
//     const IndexType PointNumber,
//     const GeometryType::IntegrationMethod& rIntegrationMethod
//     )
// {
//     const auto& r_geometry = GetGeometry();
//     const SizeType number_of_nodes = r_geometry.PointsNumber();
//     const SizeType dimension = r_geometry.WorkingSpaceDimension();

//     const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);
//     // Shape functions
//     rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(Vector &rResult, const CoordinatesArrayType& rCoordinates);
//     rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

//     KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

//     // Compute B

//     rThisKinematicVariables.B = rThisKinematicVariables.DN_DX

// }

/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{

        const GeometryType& geom = GetGeometry();
        IntegrationMethod ThisIntegrationMethod = geom.GetDefaultIntegrationMethod();
        const SizeType number_of_nodes = geom.PointsNumber();
        const GeometryType::IntegrationPointsArrayType& integration_points = geom.IntegrationPoints(ThisIntegrationMethod);
        SizeType dim1= geom.WorkingSpaceDimension();
        SizeType dim2= geom.LocalSpaceDimension();

        std::cout << "i am calculating an acoustic element\n";

        if( rLeftHandSideMatrix.size1() != number_of_nodes || rLeftHandSideMatrix.size2() != number_of_nodes )
        {
            rLeftHandSideMatrix.resize(number_of_nodes, number_of_nodes, false);
        }

        noalias(rLeftHandSideMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );

        ShapeFunctionDerivativesArrayType DN_DX;
        Vector DetJ;

        DN_DX = geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, ThisIntegrationMethod);


        // KRATOS_WATCH(integration_points.size())
        for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number )
        {
            double int_weight = integration_points[point_number].Weight() * DetJ(point_number);           
            noalias( rLeftHandSideMatrix ) += int_weight * prod( DN_DX[point_number], trans(DN_DX[point_number]));

        }







        Matrix J0;
        GeometryUtils::JacobianOnInitialConfiguration(geom, geom.IntegrationPoints(ThisIntegrationMethod)[2], J0);
        // double detJ0;
        // MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
        // const Matrix& rDN_De =
        // GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
        // GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
        // return detJ0;
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
    //const double freq2 =  std::pow(rCurrentProcessInfo[FREQUENCY], 2);
    if( rMassMatrix.size1() != number_of_nodes || rMassMatrix.size2() != number_of_nodes )
    {
        rMassMatrix.resize(number_of_nodes, number_of_nodes, false);
    }
    noalias(rMassMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );
 
    const double p = GetProperties()[DENSITY];
    const double G = GetProperties()[BULK_MODULUS];

    for (IndexType i = 0; i < number_of_nodes; i++)
    {
        for (IndexType j = 0; j < number_of_nodes; j++)
        { 
            for (IndexType g = 0; g < NumGauss; g++)
                {
                    double DetJ = geom.DeterminantOfJacobian(g, ThisIntegrationMethod);
                    double GaussWeight = DetJ * integration_points[g].Weight();
                    rMassMatrix(i,j) += NContainer(g, i) * NContainer(g, j) * GaussWeight * (p/G);
                }
        }        
    }

}


/***********************************************************************************/
void AcousticElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo)
{
    const GeometryType& geom = GetGeometry();
    const SizeType number_of_nodes = geom.PointsNumber();
    if( rDampingMatrix.size1() != number_of_nodes || rDampingMatrix.size2() != number_of_nodes )
    {
        rDampingMatrix.resize(number_of_nodes, number_of_nodes, false);
    }
    noalias(rDampingMatrix) = ZeroMatrix( number_of_nodes, number_of_nodes );

}

/***********************************************************************************/

void AcousticElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
 	
    const GeometryType& geom = GetGeometry();
    const SizeType number_of_nodes = geom.PointsNumber();
    // Resizing as needed the RHS
    // if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != number_of_nodes )
            rRightHandSideVector.resize( number_of_nodes, false );

        rRightHandSideVector = ZeroVector( number_of_nodes ); //resetting RHS
    // }
		
		
}


/***********************************************************************************/
/***********************************************************************************/

void AcousticElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

    // auto& r_geometry = this->GetGeometry();
    // const SizeType number_of_nodes = r_geometry.size();
    // const SizeType dimension = r_geometry.WorkingSpaceDimension();
    // const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();

    // KinematicVariables this_kinematic_variables(number_of_nodes, dimension);

    // const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());



    KRATOS_CATCH( "" )

}


/***********************************************************************************/


// void AcousticElement::save( Serializer& rSerializer ) const
// {
//     KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseSolidElement );
// }

/***********************************************************************************/
/***********************************************************************************/

// void AcousticElement::load( Serializer& rSerializer )
// {
//     KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseSolidElement );
// }

} // Namespace Kratos


