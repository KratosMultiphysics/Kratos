// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/checks.h"

// Application includes
#include "custom_elements/small_displacement_mixed_strain_displacement_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

Element::Pointer SmallDisplacementMixedStrainDisplacementElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedStrainDisplacementElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedStrainDisplacementElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<SmallDisplacementMixedStrainDisplacementElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer SmallDisplacementMixedStrainDisplacementElement::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    SmallDisplacementMixedStrainDisplacementElement::Pointer p_new_elem = Kratos::make_intrusive<SmallDisplacementMixedStrainDisplacementElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

void SmallDisplacementMixedStrainDisplacementElement::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry     = GetGeometry();
    const SizeType n_nodes     = r_geometry.PointsNumber();
    const SizeType dim         = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType dof_size    = n_nodes * (dim + strain_size);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size);
    }

    const IndexType displ_pos  = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType strain_pos = r_geometry[0].GetDofPosition(NODAL_STRAIN_VECTOR_XX);

    IndexType aux_index = 0;
    if (dim == 2) {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X,          displ_pos    ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y,          displ_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XX, strain_pos    ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_YY, strain_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XY, strain_pos + 2).EquationId();
        }
    } else {
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, displ_pos             ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, displ_pos + 1         ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, displ_pos + 2         ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XX, strain_pos    ).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_YY, strain_pos + 1).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_ZZ, strain_pos + 2).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XY, strain_pos + 3).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_YZ, strain_pos + 4).EquationId();
            rResult[aux_index++] = r_geometry[i_node].GetDof(NODAL_STRAIN_VECTOR_XZ, strain_pos + 5).EquationId();
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto& r_geometry     = GetGeometry();
    const SizeType n_nodes     = r_geometry.PointsNumber();
    const SizeType dim         = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType dof_size    = n_nodes * (dim + strain_size);

    if (rElementalDofList.size() != dof_size){
        rElementalDofList.resize(dof_size);
    }

    if (dim == 2) {
        for(IndexType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * dof_size    ] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * dof_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * dof_size + 2] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XX);
            rElementalDofList[i * dof_size + 3] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_YY);
            rElementalDofList[i * dof_size + 4] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XY);
        }
    } else if (dim == 3) {
        for(IndexType i = 0; i < n_nodes; ++i) {
            rElementalDofList[i * dof_size    ] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_X);
            rElementalDofList[i * dof_size + 1] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
            rElementalDofList[i * dof_size + 2] = this->GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
            rElementalDofList[i * dof_size + 3] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XX);
            rElementalDofList[i * dof_size + 4] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_YY);
            rElementalDofList[i * dof_size + 5] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_ZZ);
            rElementalDofList[i * dof_size + 6] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XY);
            rElementalDofList[i * dof_size + 7] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_YZ);
            rElementalDofList[i * dof_size + 8] = this->GetGeometry()[i].pGetDof(NODAL_STRAIN_VECTOR_XZ);
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::Initialize(const ProcessInfo &rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialization should not be done again in a restart!
    if (!rCurrentProcessInfo[IS_RESTARTED]) {
        // Integration method initialization
        mThisIntegrationMethod = GeometryData::IntegrationMethod::GI_GAUSS_2;
        const auto& r_integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

        // Constitutive Law Vector initialisation
        if (mConstitutiveLawVector.size() != r_integration_points.size()) {
            mConstitutiveLawVector.resize(r_integration_points.size());
        }

        // Initialize material
        InitializeMaterial();

    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // TODO

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // TODO

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    const SizeType block_size  = dim + strain_size;
    const SizeType matrix_size = block_size * n_nodes;

    // Check RHS size
    if (rRHS.size() != matrix_size) {
        rRHS.resize(matrix_size, false);
    }

    // Check LHS size
    if (rLHS.size1() != matrix_size || rLHS.size2() != matrix_size) {
        rLHS.resize(matrix_size, matrix_size, false);
    }
    noalias(rLHS) = ZeroMatrix(matrix_size);
    noalias(rRHS) = ZeroVector(matrix_size);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{

}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{

}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::InitializeMaterial()
{
    KRATOS_TRY

    const auto& r_properties = GetProperties();
    if (r_properties[CONSTITUTIVE_LAW] != nullptr) {
        const auto& r_geometry = GetGeometry();
        const auto& r_N_values = r_geometry.ShapeFunctionsValues(mThisIntegrationMethod);
        IndexType aux = 0;
        for (auto &it_gauss_pt : mConstitutiveLawVector) {
            it_gauss_pt = (r_properties[CONSTITUTIVE_LAW])->Clone();
            (it_gauss_pt)->InitializeMaterial(r_properties, r_geometry, row(r_N_values, aux));
            aux++;
        }
    } else {
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

bool SmallDisplacementMixedStrainDisplacementElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateN_EpsilonMatrix(
    Matrix &rN_Epsilon,
    const Vector &rN)
{
    const auto& r_geometry = GetGeometry();
    const SizeType dim     = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();
    rN_Epsilon.clear();

    if(dim == 2) {
        for (IndexType i = 0; i < n_nodes; ++i) {
            const IndexType initial_index = i * strain_size;
            rN_Epsilon(0, initial_index    ) = rN[i];
            rN_Epsilon(1, initial_index + 1) = rN[i];
            rN_Epsilon(2, initial_index + 2) = rN[i];
        }
    } else { // dim = 3
        for (IndexType i = 0; i < n_nodes; ++i) {
            const IndexType initial_index = i * strain_size;
            rN_Epsilon(0, initial_index    ) = rN[i];
            rN_Epsilon(1, initial_index + 1) = rN[i];
            rN_Epsilon(2, initial_index + 2) = rN[i];
            rN_Epsilon(3, initial_index + 3) = rN[i];
            rN_Epsilon(4, initial_index + 4) = rN[i];
            rN_Epsilon(5, initial_index + 5) = rN[i];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    ) const
{
    // Here we essentially set the input parameters
    rValues.SetStrainVector(rThisKinematicVariables.EquivalentStrain); // Equivalent stabilized total strain
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N);        // shape functions

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D);
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    ) const
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); // Here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

Vector SmallDisplacementMixedStrainDisplacementElement::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber) const
{
    const SizeType dim = GetGeometry().WorkingSpaceDimension();
    const auto aux_body_force = StructuralMechanicsElementUtilities::GetBodyForce(*this, rIntegrationPoints, PointNumber);
    Vector body_force(dim);
    for (IndexType d = 0; d < dim; ++d) {
        body_force(d) = aux_body_force(d);
    }
    return body_force;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod) const
{
    // const auto& r_geometry = GetGeometry();
    // const auto& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);

    // // Shape functions
    // rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    // // Calculate the inverse Jacobian
    // GeometryUtils::JacobianOnInitialConfiguration(
    //     r_geometry,
    //     r_integration_points[PointNumber],
    //     rThisKinematicVariables.J0);
    // MathUtils<double>::InvertMatrix(
    //     rThisKinematicVariables.J0,
    //     rThisKinematicVariables.InvJ0,
    //     rThisKinematicVariables.detJ0);
    // KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0)
    //     << "Element ID: " << this->Id() << " is inverted. det(J0) = " << rThisKinematicVariables.detJ0 << std::endl;

    // // Calculate the shape functions gradients
    // GeometryUtils::ShapeFunctionsGradients(
    //     r_geometry.ShapeFunctionsLocalGradients(rIntegrationMethod)[PointNumber],
    //     rThisKinematicVariables.InvJ0,
    //     rThisKinematicVariables.DN_DX);

    // // Compute B
    // CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.DN_DX);

    // // Calculate the equivalent total strain
    // CalculateEquivalentStrain(rThisKinematicVariables);

    // // Compute equivalent F
    // ComputeEquivalentF(rThisKinematicVariables.F, rThisKinematicVariables.EquivalentStrain);
    // rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX
    ) const
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateB(*this, rDN_DX, rB);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const
{
    // const auto& r_geom = GetGeometry();
    // const SizeType n_nodes = r_geom.PointsNumber();
    // const SizeType dim = r_geom.WorkingSpaceDimension();
    // const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // // Set the identity in Voigt notation
    // Vector voigt_identity = ZeroVector(strain_size);
    // for (IndexType d = 0; d < dim; ++d) {
    //     voigt_identity[d] = 1.0;
    // }

    // // Calculate the total strain
    // Vector eq_strain = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);

    // // Calculate the volumetric contribution
    // const Vector T_B_U = prod(mAnisotropyTensor, eq_strain);
    // double gauss_vol_strain = inner_prod(voigt_identity, T_B_U);
    // for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //     gauss_vol_strain -= rThisKinematicVariables.N[i_node] * rThisKinematicVariables.VolumetricNodalStrains[i_node];
    // }
    // gauss_vol_strain /= dim;

    // // Substract the volumetric contribution to the total strain
    // eq_strain -= gauss_vol_strain * prod(mInverseAnisotropyTensor, voigt_identity);

    // // Save the obtained value in the kinematics container
    // rThisKinematicVariables.EquivalentStrain = eq_strain;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::ComputeEquivalentF(
    Matrix& rF,
    const Vector& rStrainTensor
    ) const
{
    StructuralMechanicsElementUtilities::ComputeEquivalentF(*this, rStrainTensor, rF);
}

/***********************************************************************************/
/***********************************************************************************/

int  SmallDisplacementMixedStrainDisplacementElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // int check = SmallDisplacementMixedStrainDisplacementElement::BaseType::Check(rCurrentProcessInfo);

    // // Base check
    // check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    // // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    // const auto& r_geometry = this->GetGeometry();
    // for ( IndexType i = 0; i < r_geometry.size(); i++ ) {
    //     const NodeType& r_node = r_geometry[i];
    //     KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUMETRIC_STRAIN,r_node)

    //     KRATOS_CHECK_DOF_IN_NODE(VOLUMETRIC_STRAIN, r_node)
    // }

    // return check;

    return 1;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has(rVariable))
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    else
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // const auto& r_geometry = GetGeometry();
    // const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // const SizeType n_gauss = r_integration_points.size();
    // if (rOutput.size() != n_gauss) {
    //     rOutput.resize(n_gauss);
    // }

    // if (mConstitutiveLawVector[0]->Has( rVariable)) {
    //     GetValueOnConstitutiveLaw(rVariable, rOutput);
    // } else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
    //     // Create and initialize element variables:
    //     const SizeType n_nodes = r_geometry.PointsNumber();
    //     const SizeType dim = r_geometry.WorkingSpaceDimension();
    //     const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

    //     // Create the kinematics container and fill the nodal data
    //     KinematicVariables kinematic_variables(strain_size, dim, n_nodes);
    //     for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
    //         const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
    //         for (IndexType d = 0; d < dim; ++d) {
    //             kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
    //         }
    //         kinematic_variables.VolumetricNodalStrains[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    //     }

    //     // Create the constitutive variables and values containers
    //     ConstitutiveVariables constitutive_variables(strain_size);
    //     ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    //     auto& r_cons_law_options = cons_law_values.GetOptions();
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    //     r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);


    //     for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
    //         // Calculate kinematics
    //         CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

    //         // Call the constitutive law to update material variables
    //         if( rVariable == CAUCHY_STRESS_VECTOR) {
    //             // Compute material response
    //             CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
    //         } else {
    //             // Compute material response
    //             CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_PK2);
    //         }

    //         // Check sizes and save the output stress
    //         if (rOutput[i_gauss].size() != strain_size) {
    //             rOutput[i_gauss].resize(strain_size, false);
    //         }
    //         rOutput[i_gauss] = constitutive_variables.StressVector;
    //     }
    // } else {
    //     CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters SmallDisplacementMixedStrainDisplacementElement::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static"],
        "framework"                  : "lagrangian",
        "symmetric_lhs"              : true,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : ["CAUCHY_STRESS_VECTOR"],
            "nodal_historical"       : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","VOLUMETRIC_STRAIN"],
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3", "Quadrilateral2D4", "Tetrahedra3D4","Hexahedra3D8"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["PlaneStrain","PlaneStress","ThreeDimensional"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element implements a mixed displacement - volumetric strain formulation with Variational MultiScales (VMS) stabilization. This formulation is capable to deal with materials in the incompressible limit as well as with anisotropy."
    })");

    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    if (dimension == 2) {
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y","VOLUMETRIC_STRAIN"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","VOLUMETRIC_STRAIN"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SmallDisplacementMixedStrainDisplacementElement::BaseType);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementMixedStrainDisplacementElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SmallDisplacementMixedStrainDisplacementElement::BaseType);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

} // Namespace Kratos
