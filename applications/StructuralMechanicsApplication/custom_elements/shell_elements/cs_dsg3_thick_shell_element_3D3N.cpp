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
// #include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        if (UseGeometryIntegrationMethod()) {
            if (GetProperties().Has(INTEGRATION_ORDER) ) {
                mThisIntegrationMethod = static_cast<GeometryData::IntegrationMethod>(GetProperties()[INTEGRATION_ORDER] - 1);
            } else {
                mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
            }
        }

        const auto& r_integration_points = this->IntegrationPoints(mThisIntegrationMethod);

        // Constitutive Law initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size())
            mConstitutiveLawVector.resize(r_integration_points.size());
        InitializeMaterial();
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::InitializeMaterial()
{
    KRATOS_TRY

    // TODO: ensure retro-compatibility with older SHELLS using Section class

    const auto& r_properties = GetProperties();
    const auto& r_geometry = GetGeometry();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
        auto N_values = Vector();
        for (IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number) {
            mConstitutiveLawVector[point_number] = r_properties[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, N_values);
        }
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the CS-DSG3 shell element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH("")
}

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
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dofs_per_node = GetDoFsPerNode();

    IndexType local_index = 0;

    if (rResult.size() != dofs_per_node * number_of_nodes)
        rResult.resize(dofs_per_node * number_of_nodes, false);

    const IndexType xpos    = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType rot_pos = r_geometry[0].GetDofPosition(ROTATION_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_X, xpos    ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Y, xpos + 1).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(DISPLACEMENT_Z, xpos + 2).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_X, rot_pos     ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Y, rot_pos + 1 ).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(ROTATION_Z, rot_pos + 2 ).EquationId();
    }
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

    const auto& r_geom = GetGeometry();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType dofs_per_node = GetDoFsPerNode(); // u, v, w, theta_x, theta_y, theta_z
    rElementalDofList.resize(dofs_per_node * number_of_nodes);
    IndexType index = 0;

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rElementalDofList[index++] = r_geom[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[index++] = r_geom[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[index++] = r_geom[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_X);
        rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_Y);
        rElementalDofList[index++] = r_geom[i].pGetDof(ROTATION_Z);
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/
double CSDSG3ThickShellElement3D3N::CalculateArea(
    const array_3& r_coord_1, 
    const array_3& r_coord_2, 
    const array_3& r_coord_3 
) const
{
    const double x21 = r_coord_2[0] - r_coord_1[0];
    const double y21 = r_coord_2[1] - r_coord_1[1];
    const double x31 = r_coord_3[0] - r_coord_1[0];
    const double y31 = r_coord_3[1] - r_coord_1[1];
    return 0.5 * (x21 * y31 - y21 * x31);
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::CalculateRotationMatrixLocalToGlobal(
    bounded_3_matrix& rRotationMatrix
) const
{
    const auto& r_geometry = GetGeometry();
    array_3 v1, v2, v3; // basis vectors

    if (this->Has(LOCAL_AXIS_1)) {
        noalias(v1) = this->GetValue(LOCAL_AXIS_1); // We assume that the user has set a unit vector
        noalias(v2) = r_geometry[2] - r_geometry[0];
        v2 = v2 - inner_prod(v1, v2) * v1; // v2 orthogonal to v1
        const double norm_v2 = norm_2(v2);
        if (norm_v2 <= 1.0e-8) { // colineal
            noalias(v2) = r_geometry[1] - r_geometry[0];
            v2 = v2 - inner_prod(v1, v2) * v1; // v2 orthogonal to v1
        }
        v2 /= norm_2(v2);
    } else {
        noalias(v1) = r_geometry[1] - r_geometry[0];
        const double norm_v1 = norm_2(v1);
        KRATOS_DEBUG_ERROR_IF_NOT(norm_v1 > 0.0) << "Zero length local axis 1 for CSDSG3ThickShellElement3D3N " << this->Id() << std::endl;
        v1 /= norm_v1;
        noalias(v2) = r_geometry[2] - r_geometry[0];
        v2 = v2 - inner_prod(v1, v2) * v1; // v2 orthogonal to v1
        const double norm_v2 = norm_2(v2);
        KRATOS_DEBUG_ERROR_IF_NOT(norm_v2 > 0.0) << "Zero length local axis 2 for CSDSG3ThickShellElement3D3N " << this->Id() << std::endl;
        v2 /= norm_v2; 
    }
    noalias(v3) = MathUtils<double>::CrossProduct(v1, v2);

    // Assemble the basis vectors in the rotation matrix
    for (IndexType i = 0; i < 3; ++i) { // in rows, global to local
        rRotationMatrix(0, i) = v1[i];
        rRotationMatrix(1, i) = v2[i];
        rRotationMatrix(2, i) = v3[i];
    }

}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::RotateLHSToGlobal(
    MatrixType& rLHS,
    const bounded_3_matrix& rRotationMatrix // provided
) const
{
    // We will perform the rotation by blocks of 3x3 for efficiency
    bounded_3_matrix local_LHS_block;
    bounded_3_matrix temp_matrix;
    bounded_3_matrix temp_matrix2;

    // Loop over the six blocks in each direction
    for (IndexType i = 0; i < 6; ++i) {
        for (IndexType j = 0; j < 6; ++j) {
            // We extract the 3x3 block
            for (IndexType k = 0; k < 3; ++k) {
                for (IndexType l = 0; l < 3; ++l) {
                    local_LHS_block(k, l) = rLHS(i * 3 + k, j * 3 + l);
                }
            }

            // We perform the rotation to the block
            noalias(temp_matrix) = prod(trans(rRotationMatrix), local_LHS_block);
            noalias(temp_matrix2) = prod(temp_matrix, rRotationMatrix);

            // We store the values back in the LHS
            for (IndexType k = 0; k < 3; ++k) {
                for (IndexType l = 0; l < 3; ++l) {
                    rLHS(i * 3 + k, j * 3 + l) = temp_matrix2(k, l);
                }
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::RotateRHSToGlobal(
    VectorType& rRHS,
    const bounded_3_matrix& rRotationMatrix // provided
) const
{
    // We will perform the rotation by blocks of 3 for efficiency
    array_3 local_RHS_block;
    array_3 temp_vector;

    for (IndexType i = 0; i < 6; ++i) {
        // We extract the 3 components of the block
        for (IndexType k = 0; k < 3; ++k) {
            local_RHS_block[k] = rRHS[i * 3 + k];
        }
        // We perform the rotation to the block
        noalias(temp_vector) = prod(trans(rRotationMatrix), local_RHS_block);
        // We store the values back in the RHS
        for (IndexType k = 0; k < 3; ++k) {
            rRHS[i * 3 + k] = temp_vector[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::CalculateBTriangle(
    MatrixType& rB,
    const bounded_3_matrix& r_rotation_matrix,
    const array_3& r_coord_1, 
    const array_3& r_coord_2, 
    const array_3& r_coord_3,
    const double local_coord_1,
    const double local_coord_2,
    const double local_coord_3
)
{
    const IndexType strain_size = GetStrainSize();
    const IndexType number_of_nodes = GetGeometry().PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();

    if (rB.size1() != strain_size || rB.size2() != system_size)
        rB.resize(strain_size, system_size, false);
    rB.clear();

    const double alpha = 1.5;

    // beta parameters for the membrane part
    // const double b1 = 1.0;
    // const double b2 = 2.0;
    // const double b3 = 1.0;
    // const double b4 = 0.0;
    // const double b5 = 1.0;
    // const double b6 = -1.0;
    // const double b7 = -1.0;
    // const double b8 = -1.0;
    // const double b9 = -2.0;

    // Here we rotate to local coordinates (z is normal to the element)
    array_3 local_coords_1;
    array_3 local_coords_2;
    array_3 local_coords_3;
    noalias(local_coords_1) = prod(r_rotation_matrix, r_coord_1);
    noalias(local_coords_2) = prod(r_rotation_matrix, r_coord_2);
    noalias(local_coords_3) = prod(r_rotation_matrix, r_coord_3);

    const double x1 = local_coords_1[0];
    const double y1 = local_coords_1[1];
    const double x2 = local_coords_2[0];
    const double y2 = local_coords_2[1];
    const double x3 = local_coords_3[0];
    const double y3 = local_coords_3[1];

    const double area = CalculateArea(local_coords_1, local_coords_2, local_coords_3);

    const double x12 = x1 - x2;
    const double x23 = x2 - x3;
    const double x31 = x3 - x1;
    const double x32 = -x23;
    const double x13 = -x31;
    const double x21 = -x12;

    const double y12 = y1 - y2;
    const double y23 = y2 - y3;
    const double y31 = y3 - y1;
    const double y32 = -y23;
    const double y13 = -y31;
    const double y21 = -y12;

    const double aux_prod = 0.5 / area;

    // Membrane components with drilling rotations (Zhang et al 2011)
    // CST membrane part
    rB(0, 0) = y23;
    rB(0, 1) = x32;
    rB(2, 0) = x32;
    rB(2, 1) = y23;

    rB(0, 6) = y31;
    rB(0, 7) = x13;
    rB(2, 6) = x13;
    rB(2, 7) = y31;

    rB(0, 12) = y12;
    rB(0, 13) = x21;
    rB(2, 12) = x21;
    rB(2, 13) = y12;

    const double temp = alpha / 6.0;
    const double temp2 = 2.0 * temp;
    // Drilling rotation membrane components
    rB(0, 5) = temp * y23 * (y13 + y12);
    rB(1, 5) = temp * x32 * (x31 + x21);
    rB(2, 5) = temp2 * (x31 * y13 + x21 * y21);

    rB(0, 11) = temp * y31 * (y21 + y23);
    rB(1, 11) = temp * x13 * (x12 + x32);
    rB(2, 11) = temp2 * (x12 * y21 + x32 * y32);

    rB(0, 17) = temp * y12 * (y32 + y31);
    rB(1, 17) = temp * x21 * (x23 + x13);
    rB(2, 17) = temp2 * (x23 * y32 + x13 * y13);

    // SECOND ORDER MISSING!!!! 
    // TODO



    // ...

    // Bending components
    rB(3, 4) = y23;
    rB(4, 3) = x23;
    rB(5, 3) = y32;
    rB(5, 4) = x32;

    rB(3, 10) = y31;
    rB(4, 9) = x31;
    rB(5, 9) = y13;
    rB(5, 10) = x13;

    rB(3, 16) = y12;
    rB(4, 15) = x12;
    rB(5, 15) = y21;
    rB(5, 16) = x21;

    // Shear components
    // Bs1
    rB(6, 2) = y21 + y13;
    rB(6, 4) = area;
    rB(7, 2) = x31 + x12;
    rB(7, 3) = -area;
    // Bs2
    rB(6, 8) = y31;
    rB(6, 9) = y12 * y31 * 0.5;
    rB(6, 10) = x21 * y31 * 0.5;
    rB(7, 8) = x13;
    rB(7, 9) = y21 * x31 * 0.5;
    rB(7, 10) = x12 * x31 * 0.5;
    // Bs3
    rB(6, 14) = y12;
    rB(6, 15) = y21 * y31 * 0.5;
    rB(6, 16) = y21 * x31 * 0.5;
    rB(7, 14) = x21;
    rB(7, 15) = x12 * y31 * 0.5;
    rB(7, 16) = x21 * x31 * 0.5;

    rB *= aux_prod;
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::GetNodalValuesVector(VectorType& rNodalValues) const
{
    const auto& r_geometry = GetGeometry();

    IndexType index = 0;
    for (IndexType i = 0; i < r_geometry.PointsNumber(); ++i) {
        rNodalValues[index++] = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT_X);
        rNodalValues[index++] = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT_Y);
        rNodalValues[index++] = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT_Z);
        rNodalValues[index++] = r_geometry[i].FastGetSolutionStepValue(ROTATION_X);
        rNodalValues[index++] = r_geometry[i].FastGetSolutionStepValue(ROTATION_Y);
        rNodalValues[index++] = r_geometry[i].FastGetSolutionStepValue(ROTATION_Z);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void CSDSG3ThickShellElement3D3N::CalculateB(
    MatrixType &rB
)
{

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

    const IndexType strain_size = GetStrainSize();
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();
    const IndexType number_of_nodes = r_geometry.PointsNumber();
    const IndexType system_size = number_of_nodes * GetDoFsPerNode();
    
    bounded_3_matrix rotation_matrix;
    CalculateRotationMatrixLocalToGlobal(rotation_matrix);

    if (rLHS.size1() != system_size || rLHS.size2() != system_size)
        rLHS.resize(system_size, system_size, false);
    rLHS.clear();

    if (rRHS.size() != system_size)
        rRHS.resize(system_size, false);
    rRHS.clear();

    const double thickness = r_props[THICKNESS];
    const double area = r_geometry.Area();

    VectorType nodal_values(system_size);
    GetNodalValuesVector(nodal_values);
    // We rotate the nodal values to the local system of the shell
    RotateRHSToGlobal(nodal_values, rotation_matrix);

    ConstitutiveLaw::Parameters cl_values(r_geometry, r_props, rProcessInfo);
    auto &r_cl_options = cl_values.GetOptions();
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cl_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Let's initialize the constitutive law's values
    VectorType strain_vector(strain_size), stress_vector(strain_size);
    MatrixType constitutive_matrix(strain_size, strain_size);
    cl_values.SetStrainVector(strain_vector);
    cl_values.SetStressVector(stress_vector);
    cl_values.SetConstitutiveMatrix(constitutive_matrix);

    const auto& r_integration_points = CustomTriangleAreaCoordinatesQuadrature(area);
    double zeta1, zeta2, zeta3, weight;
    MatrixType B(strain_size, system_size);
    for (SizeType IP = 0; IP < r_integration_points.size(); ++IP) {
        zeta1 = r_integration_points[IP].X();
        zeta2 = r_integration_points[IP].Y();
        zeta3 = r_integration_points[IP].Z();
        weight = r_integration_points[IP].Weight();

        CalculateBTriangle(B, rotation_matrix, r_geometry[0].Coordinates(), r_geometry[1].Coordinates(), r_geometry[2].Coordinates(), zeta1, zeta2, zeta3);

        // We compute the strain at the integration point
        noalias(strain_vector) = prod(B, nodal_values);

        // We call the constitutive law to compute the stress
        cl_values.SetStrainVector(strain_vector);
        mConstitutiveLawVector[IP]->CalculateMaterialResponseCauchy(cl_values);
        noalias(stress_vector) = cl_values.GetStressVector();
        noalias(constitutive_matrix) = cl_values.GetConstitutiveMatrix();

        // We integrate the LHS and RHS
        noalias(rLHS) += weight * prod(trans(B), Matrix(prod(constitutive_matrix, B)));
        noalias(rRHS) -= weight * prod(trans(B), stress_vector);
    }

    RotateLHSToGlobal(rLHS, rotation_matrix);
    RotateRHSToGlobal(rRHS, rotation_matrix);

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

    const auto& r_properties = GetProperties();
    KRATOS_ERROR_IF_NOT(r_properties.Has(THICKNESS)) << "THICKNESS not provided for CSDSG3ThickShellElement3D3N " << this->Id() << std::endl;
    KRATOS_ERROR_IF_NOT(r_properties.GetValue(THICKNESS) > 0.0) << "Wrong value for THICKNESS in the  CSDSG3ThickShellElement3D3N " << this->Id() << std::endl;
    return mConstitutiveLawVector[0]->Check(r_properties, GetGeometry(), rCurrentProcessInfo);

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
