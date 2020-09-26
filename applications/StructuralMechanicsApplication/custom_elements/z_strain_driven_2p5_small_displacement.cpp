// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Miguel Angel Celigueta
//                   Ignasi de Pouplana
//

// Application includes
#include "custom_elements/z_strain_driven_2p5_small_displacement.h"

namespace Kratos
{
ZStrainDriven2p5DSmallDisplacement::ZStrainDriven2p5DSmallDisplacement( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseType( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

ZStrainDriven2p5DSmallDisplacement::ZStrainDriven2p5DSmallDisplacement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseType( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer ZStrainDriven2p5DSmallDisplacement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<ZStrainDriven2p5DSmallDisplacement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer ZStrainDriven2p5DSmallDisplacement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<ZStrainDriven2p5DSmallDisplacement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer ZStrainDriven2p5DSmallDisplacement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    ZStrainDriven2p5DSmallDisplacement::Pointer p_new_elem = Kratos::make_intrusive<ZStrainDriven2p5DSmallDisplacement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    // The vector containing the imposed Z strain
    p_new_elem->mImposedZStrainVector = mImposedZStrainVector;

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void ZStrainDriven2p5DSmallDisplacement::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Initialize(rCurrentProcessInfo);

    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    //Imposed Z strain vector initialisation
    if ( mImposedZStrainVector.size() != integration_points.size() )
        mImposedZStrainVector.resize( integration_points.size() );

    for ( IndexType point_number = 0; point_number < mImposedZStrainVector.size(); ++point_number ) {
        mImposedZStrainVector[point_number] = 0.0;
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void ZStrainDriven2p5DSmallDisplacement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    BaseType::SetConstitutiveVariables(rThisKinematicVariables,rThisConstitutiveVariables,
                                                rValues,PointNumber,IntegrationPoints);

    // StrainVector must have the shape of a 3D element
    const double eps_xy = rThisConstitutiveVariables.StrainVector[2];
    rThisConstitutiveVariables.StrainVector[2] = mImposedZStrainVector[PointNumber];
    rThisConstitutiveVariables.StrainVector[3] = eps_xy;

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // B matrix must have the shape of a 3D element
    for (IndexType i = 0; i < mat_size; ++i ) {
        rThisKinematicVariables.B(3,i) = rThisKinematicVariables.B(2,i);
        rThisKinematicVariables.B(2,i) = 0.0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ZStrainDriven2p5DSmallDisplacement::SetValuesOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == IMPOSED_Z_STRAIN_VALUE) {
        const SizeType integration_points_number = mImposedZStrainVector.size();
        for ( IndexType point_number = 0; point_number < integration_points_number; ++point_number ) {
            mImposedZStrainVector[point_number] = rValues[point_number];
        }
    } else {
        BaseType::SetValuesOnIntegrationPoints(rVariable,rValues,rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

int  ZStrainDriven2p5DSmallDisplacement::Check( const ProcessInfo& rCurrentProcessInfo ) const
{
    KRATOS_TRY

    int check = Element::Check(rCurrentProcessInfo);

    const auto& r_geometry = this->GetGeometry();
    const auto& r_properties = this->GetProperties();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(r_properties.Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << r_properties.Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strain_size = r_properties.GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( dimension == 2 ) {
        KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 2.5D element! expected strain size is 6 (el id = ) "<<  this->Id() << std::endl;
    } else {
        KRATOS_ERROR << "Wrong dimension. The 2.5D element must have a 2D geometry." << std::endl;
    }

    // Check constitutive law
    if ( BaseType::mConstitutiveLawVector.size() > 0 ) {
        check = BaseType::mConstitutiveLawVector[0]->Check( r_properties, r_geometry, rCurrentProcessInfo );
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void ZStrainDriven2p5DSmallDisplacement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    rSerializer.save("ImposedZStrainVector", mImposedZStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

void ZStrainDriven2p5DSmallDisplacement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    rSerializer.load("ImposedZStrainVector", mImposedZStrainVector);
}

} // Namespace Kratos


