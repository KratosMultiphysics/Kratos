// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "cs_dsg3_thick_shell_element_3D3N.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    // if (!rCurrentProcessInfo[IS_RESTARTED]) {
    //     if (this->UseGeometryIntegrationMethod()) {
    //         if (GetProperties().Has(INTEGRATION_ORDER) ) {
    //             mThisIntegrationMethod = static_cast<GeometryData::IntegrationMethod>(GetProperties()[INTEGRATION_ORDER] - 1);
    //         } else {
    //             mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_3;
    //         }
    //     }

    //     const auto& r_integration_points = this->IntegrationPoints(mThisIntegrationMethod);

    //     // Constitutive Law initialisation
    //     if (mConstitutiveLawVector.size() != r_integration_points.size())
    //         mConstitutiveLawVector.resize(r_integration_points.size());
    //     InitializeMaterial();
    // }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

// void CSDSG3ThickShellElement3D3N::InitializeMaterial()
// {
//     KRATOS_TRY

//     if (GetProperties()[CONSTITUTIVE_LAW] != nullptr) {
//         const auto& r_geometry   = GetGeometry();
//         const auto& r_properties = GetProperties();
//         auto N_values            = Vector();
//         for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
//             mConstitutiveLawVector[point_number] = r_properties[CONSTITUTIVE_LAW]->Clone();
//             mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
//         }
//     } else
//         KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

//     KRATOS_CATCH("");
// }

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer CSDSG3ThickShellElement3D3N::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    CSDSG3ThickShellElement3D3N::Pointer p_new_elem = Kratos::make_intrusive<CSDSG3ThickShellElement3D3N>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    // const auto& r_geometry = GetGeometry();
    // const SizeType number_of_nodes = r_geometry.size();
    // const SizeType dofs_per_node = GetDoFsPerNode(); // u, v, theta
    // const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // IndexType local_index = 0;

    // if (rResult.size() != dofs_per_node * number_of_nodes)
    //     rResult.resize(dofs_per_node * number_of_nodes, false);

    // const IndexType xpos    = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    // const IndexType rot_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    // for (IndexType i = 0; i < number_of_nodes; ++i) {
    //     rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos).EquationId();
    //     rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
    //     if (dimension == 3) {
    //         rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();
    //         rResult[local_index++] = r_geometry[i].GetDof(ROTATION_X, rot_pos).EquationId();
    //         rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Y, rot_pos + 1).EquationId();
    //     }
    //     rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z, rot_pos + 2).EquationId();
    // }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    // const auto& r_geom = GetGeometry();
    // const SizeType number_of_nodes = r_geom.size();
    // const SizeType dofs_per_node = GetDoFsPerNode(); // u, v, theta
    // const SizeType dimension = r_geom.WorkingSpaceDimension();
    // rElementalDofList.resize(dofs_per_node * number_of_nodes);
    // SizeType index = 0;

    // for (IndexType i = 0; i < number_of_nodes; ++i) {
    //     rElementalDofList[index++]   = r_geom[i].pGetDof(DISPLACEMENT_X);
    //     rElementalDofList[index++] = r_geom[i].pGetDof(DISPLACEMENT_Y);
    //     if (dimension == 3) {
    //         rElementalDofList[index++] = r_geom[i].pGetDof(DISPLACEMENT_Z);
    //         rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_X);
    //         rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_Y);
    //     }
    //     rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_Z);
    // }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::CalculateB(
    Matrix &rB,
    const GeometryType::Pointer pTriangleGeometry // the geometry of the sub-triangle
)
{
    const IndexType strain_size = GetStrainSize();
    const IndexType number_of_nodes = pTriangleGeometry->PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();


    if (rB.size1() != strain_size || rB.size2() != system_size)
        rB.resize(strain_size, system_size, false);
    rB.clear();

    const double area = pTriangleGeometry->Area();
    const auto& r_points = pTriangleGeometry->Points();
    const auto& coords_1 = r_points[0].Coordinates();
    const auto& coords_2 = r_points[1].Coordinates();
    const auto& coords_3 = r_points[2].Coordinates();

    const double a = coords_2[0] - coords_1[0];
    const double b = coords_2[1] - coords_1[1];
    const double c = coords_3[1] - coords_1[1];
    const double d = coords_3[0] - coords_1[0];

    const double aux_prod = 0.5 / area;
    const double N1_x = aux_prod * (coords_2[1] - coords_3[1]);
    const double N1_y = aux_prod * (coords_3[0] - coords_2[0]);
    const double N2_x = aux_prod * (coords_3[1] - coords_1[1]);
    const double N2_y = aux_prod * (coords_1[0] - coords_3[0]);
    const double N3_x = aux_prod * (coords_1[1] - coords_2[1]);
    const double N3_y = aux_prod * (coords_2[0] - coords_1[0]);
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY;

    KRATOS_CATCH("");
}


/***********************************************************************************/
/***********************************************************************************/

int CSDSG3ThickShellElement3D3N::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    return mConstitutiveLawVector[0]->Check(GetProperties(), GetGeometry(), rCurrentProcessInfo);

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

} // Namespace Kratos
