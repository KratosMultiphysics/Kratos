// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//                   Guglielmo Scovazzi
//

// System includes

// External includes

// Project includes
#include "utilities/element_size_calculator.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "includes/checks.h"

// Application includes
#include "custom_elements/total_lagrangian_mixed_detJ_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<>
double TotalLagrangianMixedDetJElement<2>::CalculateShearModulus(const Matrix &rC) const
{
    return 0.2*(rC(0,0) - 2.0*rC(0,1) + rC(1,1) + rC(2,2));
}

/***********************************************************************************/
/***********************************************************************************/

template<>
double TotalLagrangianMixedDetJElement<3>::CalculateShearModulus(const Matrix &rC) const
{
    return (4.0 / 33.0)*(rC(0,0) - rC(0,1) - rC(0,2) + rC(1,1) - rC(1,2) + rC(2,2) + (3.0/4.0)*(rC(3,3) + rC(4,4) + rC(5,5)));
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
double TotalLagrangianMixedDetJElement<TDim>::CalculateBulkModulus(const Matrix &rC) const
{
    double bulk_modulus = 0.0;
    for (SizeType i = 0; i < TDim; ++i) {
        for (SizeType j = 0; j < TDim; ++j) {
            bulk_modulus += rC(i,j);
        }
    }
    return bulk_modulus / std::pow(TDim, 2);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Element::Pointer TotalLagrangianMixedDetJElement<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<TotalLagrangianMixedDetJElement<TDim>>(NewId, GetGeometry().Create(ThisNodes), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Element::Pointer TotalLagrangianMixedDetJElement<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<TotalLagrangianMixedDetJElement<TDim>>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Element::Pointer TotalLagrangianMixedDetJElement<TDim>::Clone(
    IndexType NewId,
    NodesArrayType const& rThisNodes) const
{
    KRATOS_TRY

    TotalLagrangianMixedDetJElement::Pointer p_new_elem = Kratos::make_intrusive<TotalLagrangianMixedDetJElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
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

template<>
void TotalLagrangianMixedDetJElement<2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != LocalSize){
        rResult.resize(LocalSize);
    }

    const auto& r_geometry = GetGeometry();
    const IndexType disp_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType det_F_pos = r_geometry[0].GetDofPosition(VOLUMETRIC_STRAIN);

    IndexType aux_index = 0;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, det_F_pos).EquationId();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rResult.size() != LocalSize){
        rResult.resize(LocalSize);
    }

    const auto& r_geometry = GetGeometry();
    const IndexType disp_pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);
    const IndexType det_F_pos = r_geometry[0].GetDofPosition(VOLUMETRIC_STRAIN);

    IndexType aux_index = 0;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, disp_pos + 2).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(VOLUMETRIC_STRAIN, det_F_pos).EquationId();
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rElementalDofList.size() != LocalSize){
        rElementalDofList.resize(LocalSize);
    }

    const auto& r_geometry = GetGeometry();
    for(IndexType i = 0; i < NumNodes; ++i) {
        rElementalDofList[i * BlockSize] = r_geometry[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[i * BlockSize + 1] = r_geometry[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[i * BlockSize + 2] = r_geometry[i].pGetDof(VOLUMETRIC_STRAIN);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    if (rElementalDofList.size() != LocalSize){
        rElementalDofList.resize(LocalSize);
    }

    const auto& r_geometry = GetGeometry();
    for(IndexType i = 0; i < NumNodes; ++i){
        rElementalDofList[i * BlockSize] = r_geometry[i].pGetDof(DISPLACEMENT_X);
        rElementalDofList[i * BlockSize + 1] = r_geometry[i].pGetDof(DISPLACEMENT_Y);
        rElementalDofList[i * BlockSize + 2] = r_geometry[i].pGetDof(DISPLACEMENT_Z);
        rElementalDofList[i * BlockSize + 3] = r_geometry[i].pGetDof(VOLUMETRIC_STRAIN);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::Initialize(const ProcessInfo &rCurrentProcessInfo)
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

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < TDim; ++d) {
            kinematic_variables.Displacements(i_node, d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Set te constitutive law values
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(GetGeometry(), GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Recompute the kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Set the constitutive variables
        SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->InitializeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node, d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Set te constitutive law values
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Call the initialize material response
    for (IndexType i_gauss = 0; i_gauss < mConstitutiveLawVector.size(); ++i_gauss) {
        // Recompute the kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

        // Set the constitutive variables
        SetConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points);

        // Call the constitutive law to update material variables
        mConstitutiveLawVector[i_gauss]->FinalizeMaterialResponseCauchy(cons_law_values);
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();

    // Check RHS size
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < 2; ++d) {
            kinematic_variables.Displacements(i_node * 2 + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    // Calculate stabilization constant
    const double c_tau_u = 2.0;
    const double c_tau_th = 4.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau_u * std::pow(h,2) / 2.0;

    // Set the auxiliary references matching the automatic differentiation symbols
    array_1d<double,3> b_gauss;
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedVector<double,LocalSize> rhs;
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_PK2);

        // Calculate body force
        // Note that this already includes the density computed in the reference configuration
        b_gauss = StructuralMechanicsElementUtilities::GetBodyForce(*this, r_integration_points, i_gauss);

        // Calculate the stabilization constants
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
        double kappa = CalculateBulkModulus(constitutive_variables.ConstitutiveMatrix); // Equivalent bulk modulus
        const double tau_u = aux_tau / mu;
        const double tau_th = (c_tau_th * mu) / (mu + kappa);

        // Calculate and add the LHS Gauss point contributions
        const double clhs0 = DN(0,0)*u(0,1);
const double clhs1 = DN(1,0)*u(1,1);
const double clhs2 = DN(2,0)*u(2,1);
const double clhs3 = clhs0 + clhs1 + clhs2;
const double clhs4 = pow(clhs3, 2);
const double clhs5 = DN(0,0)*u(0,0);
const double clhs6 = DN(1,0)*u(1,0);
const double clhs7 = DN(2,0)*u(2,0);
const double clhs8 = clhs5 + clhs6 + clhs7;
const double clhs9 = clhs8 + 1;
const double clhs10 = pow(clhs9, 2);
const double clhs11 = clhs10 + clhs4;
const double clhs12 = C(0,0)*clhs11;
const double clhs13 = DN(0,1)*u(0,0);
const double clhs14 = DN(1,1)*u(1,0);
const double clhs15 = DN(2,1)*u(2,0);
const double clhs16 = clhs13 + clhs14 + clhs15;
const double clhs17 = pow(clhs16, 2);
const double clhs18 = DN(0,1)*u(0,1);
const double clhs19 = DN(1,1)*u(1,1);
const double clhs20 = DN(2,1)*u(2,1);
const double clhs21 = clhs18 + clhs19 + clhs20;
const double clhs22 = clhs21 + 1;
const double clhs23 = pow(clhs22, 2);
const double clhs24 = clhs17 + clhs23;
const double clhs25 = C(0,1)*clhs24;
const double clhs26 = clhs22*clhs3;
const double clhs27 = clhs16*clhs9;
const double clhs28 = clhs26 + clhs27;
const double clhs29 = 2*clhs28;
const double clhs30 = C(0,2)*clhs29 + clhs12 + clhs25;
const double clhs31 = DN(0,0)*clhs30;
const double clhs32 = clhs31*tau_th;
const double clhs33 = N[0]*th[0];
const double clhs34 = N[1]*th[1];
const double clhs35 = N[2]*th[2];
const double clhs36 = clhs33 + clhs34 + clhs35;
const double clhs37 = clhs36 + 1;
const double clhs38 = pow(clhs37, -0.33333333333333331);
const double clhs39 = 0.33333333333333331*clhs38;
const double clhs40 = clhs32*clhs39;
const double clhs41 = DN(0,0)*clhs19 + DN(0,0)*clhs20 + DN(0,0) - DN(0,1)*DN(1,0)*u(1,1) - DN(0,1)*DN(2,0)*u(2,1);
const double clhs42 = -clhs0*clhs14 - clhs0*clhs15 - clhs1*clhs13 - clhs1*clhs15 - clhs13*clhs2 - clhs14*clhs2 + clhs18*clhs6 + clhs18*clhs7 + clhs19*clhs5 + clhs19*clhs7 + clhs20*clhs5 + clhs20*clhs6 + clhs21;
const double clhs43 = clhs42 + clhs9;
const double clhs44 = pow(clhs43, -0.66666666666666663);
const double clhs45 = clhs41*clhs44;
const double clhs46 = C(0,2)*clhs11;
const double clhs47 = C(1,2)*clhs24;
const double clhs48 = C(2,2)*clhs29 + clhs46 + clhs47;
const double clhs49 = DN(0,1)*clhs48;
const double clhs50 = clhs39*tau_th;
const double clhs51 = clhs49*clhs50;
const double clhs52 = DN(0,1)*clhs16;
const double clhs53 = DN(0,0)*clhs9;
const double clhs54 = DN(0,0)*clhs16;
const double clhs55 = DN(0,1)*clhs9;
const double clhs56 = clhs54 + clhs55;
const double clhs57 = C(0,0)*clhs53 + C(0,1)*clhs52 + C(0,2)*clhs56;
const double clhs58 = 0.66666666666666663*DN(0,0);
const double clhs59 = clhs57*clhs58;
const double clhs60 = -clhs33 - clhs34 - clhs35 + clhs42 + clhs8;
const double clhs61 = clhs60*tau_th;
const double clhs62 = clhs38*clhs61;
const double clhs63 = clhs44*clhs62;
const double clhs64 = C(0,2)*clhs53 + C(1,2)*clhs52 + C(2,2)*clhs56;
const double clhs65 = 0.66666666666666663*DN(0,1);
const double clhs66 = clhs64*clhs65;
const double clhs67 = pow(clhs43, -1.6666666666666665);
const double clhs68 = 0.22222222222222221*clhs67;
const double clhs69 = clhs41*clhs68;
const double clhs70 = clhs32*clhs60;
const double clhs71 = clhs38*clhs70;
const double clhs72 = clhs49*clhs62;
const double clhs73 = 1.0*clhs44;
const double clhs74 = 0.33333333333333331*clhs67;
const double clhs75 = -clhs41*clhs74;
const double clhs76 = clhs17*clhs75 + clhs23*clhs75 + clhs52*clhs73;
const double clhs77 = clhs10*clhs75 + clhs4*clhs75 + clhs53*clhs73;
const double clhs78 = 0.5*clhs44;
const double clhs79 = 2*clhs26*clhs75 + 2*clhs27*clhs75 + 2*clhs54*clhs78 + 2*clhs55*clhs78;
const double clhs80 = C(0,0)*clhs77 + C(0,1)*clhs76 + C(0,2)*clhs79;
const double clhs81 = clhs36 + 1.0;
const double clhs82 = pow(clhs81, 0.66666666666666663);
const double clhs83 = clhs82*clhs9;
const double clhs84 = C(0,2)*clhs77 + C(1,2)*clhs76 + C(2,2)*clhs79;
const double clhs85 = clhs16*clhs82;
const double clhs86 = clhs44*clhs82;
const double clhs87 = clhs28*clhs86;
const double clhs88 = 0.5*clhs10*clhs86 + 0.5*clhs4*clhs86 - 0.5;
const double clhs89 = 0.5*clhs17*clhs86 + 0.5*clhs23*clhs86 - 0.5;
const double clhs90 = C(0,0)*clhs88 + C(0,1)*clhs89 + C(0,2)*clhs87;
const double clhs91 = C(0,2)*clhs88 + C(1,2)*clhs89 + C(2,2)*clhs87;
const double clhs92 = DN(0,0)*clhs90 + DN(0,1)*clhs91;
const double clhs93 = clhs80*clhs83 + clhs84*clhs85 + clhs92;
const double clhs94 = C(0,1)*clhs77 + C(1,1)*clhs76 + C(1,2)*clhs79;
const double clhs95 = C(0,1)*clhs88 + C(1,1)*clhs89 + C(1,2)*clhs87;
const double clhs96 = DN(0,0)*clhs91 + DN(0,1)*clhs95;
const double clhs97 = clhs83*clhs84 + clhs85*clhs94 + clhs96;
const double clhs98 = -DN(0,0)*DN(1,1)*u(1,0) - DN(0,0)*DN(2,1)*u(2,0) + DN(0,1)*clhs6 + DN(0,1)*clhs7 + DN(0,1);
const double clhs99 = clhs44*clhs98;
const double clhs100 = DN(0,0)*clhs3;
const double clhs101 = DN(0,1)*clhs22;
const double clhs102 = DN(0,1)*clhs3;
const double clhs103 = DN(0,0)*clhs22;
const double clhs104 = clhs102 + clhs103;
const double clhs105 = C(0,0)*clhs100 + C(0,1)*clhs101 + C(0,2)*clhs104;
const double clhs106 = clhs105*clhs58;
const double clhs107 = C(0,2)*clhs100 + C(1,2)*clhs101 + C(2,2)*clhs104;
const double clhs108 = clhs107*clhs65;
const double clhs109 = clhs68*clhs98;
const double clhs110 = -clhs74*clhs98;
const double clhs111 = clhs10*clhs110 + clhs100*clhs73 + clhs110*clhs4;
const double clhs112 = clhs101*clhs73 + clhs110*clhs17 + clhs110*clhs23;
const double clhs113 = 2*clhs102*clhs78 + 2*clhs103*clhs78 + 2*clhs110*clhs26 + 2*clhs110*clhs27;
const double clhs114 = C(0,2)*clhs111 + C(1,2)*clhs112 + C(2,2)*clhs113;
const double clhs115 = C(0,0)*clhs111 + C(0,1)*clhs112 + C(0,2)*clhs113;
const double clhs116 = clhs114*clhs16 + clhs115*clhs9;
const double clhs117 = DN(0,0)*clhs82;
const double clhs118 = C(0,1)*clhs111 + C(1,1)*clhs112 + C(1,2)*clhs113;
const double clhs119 = clhs114*clhs9 + clhs118*clhs16;
const double clhs120 = DN(0,1)*clhs82;
const double clhs121 = 0.1111111111111111*pow(clhs37, -1.3333333333333333);
const double clhs122 = clhs121*clhs61;
const double clhs123 = pow(clhs81, -0.33333333333333337);
const double clhs124 = 0.66666666666666663*clhs28;
const double clhs125 = C(2,2)*clhs124 + 0.33333333333333331*clhs46 + 0.33333333333333331*clhs47;
const double clhs126 = C(0,2)*clhs124 + 0.33333333333333331*clhs12 + 0.33333333333333331*clhs25;
const double clhs127 = clhs125*clhs16 + clhs126*clhs9;
const double clhs128 = C(0,1)*clhs11;
const double clhs129 = C(1,1)*clhs24;
const double clhs130 = C(1,2)*clhs124 + 0.33333333333333331*clhs128 + 0.33333333333333331*clhs129;
const double clhs131 = clhs125*clhs9 + clhs130*clhs16;
const double clhs132 = clhs44*(DN(0,0)*clhs123*clhs127 + DN(0,1)*clhs123*clhs131 - clhs121*clhs70 - clhs122*clhs49 - clhs40 - clhs51);
const double clhs133 = -DN(0,0)*DN(1,1)*u(0,1) + DN(1,0)*clhs18 + DN(1,0)*clhs20 + DN(1,0) - DN(1,1)*DN(2,0)*u(2,1);
const double clhs134 = clhs133*clhs44;
const double clhs135 = DN(1,1)*clhs16;
const double clhs136 = DN(1,0)*clhs9;
const double clhs137 = DN(1,0)*clhs16;
const double clhs138 = DN(1,1)*clhs9;
const double clhs139 = clhs137 + clhs138;
const double clhs140 = C(0,0)*clhs136 + C(0,1)*clhs135 + C(0,2)*clhs139;
const double clhs141 = clhs140*clhs58;
const double clhs142 = C(0,2)*clhs136 + C(1,2)*clhs135 + C(2,2)*clhs139;
const double clhs143 = clhs142*clhs65;
const double clhs144 = clhs133*clhs68;
const double clhs145 = -clhs133*clhs74;
const double clhs146 = clhs135*clhs73 + clhs145*clhs17 + clhs145*clhs23;
const double clhs147 = clhs10*clhs145 + clhs136*clhs73 + clhs145*clhs4;
const double clhs148 = 2*clhs137*clhs78 + 2*clhs138*clhs78 + 2*clhs145*clhs26 + 2*clhs145*clhs27;
const double clhs149 = C(0,0)*clhs147 + C(0,1)*clhs146 + C(0,2)*clhs148;
const double clhs150 = C(0,2)*clhs147 + C(1,2)*clhs146 + C(2,2)*clhs148;
const double clhs151 = DN(1,0)*clhs90 + DN(1,1)*clhs91;
const double clhs152 = clhs149*clhs83 + clhs150*clhs85 + clhs151;
const double clhs153 = C(0,1)*clhs147 + C(1,1)*clhs146 + C(1,2)*clhs148;
const double clhs154 = DN(1,0)*clhs91 + DN(1,1)*clhs95;
const double clhs155 = clhs150*clhs83 + clhs153*clhs85 + clhs154;
const double clhs156 = -DN(0,1)*DN(1,0)*u(0,0) - DN(1,0)*DN(2,1)*u(2,0) + DN(1,1)*clhs5 + DN(1,1)*clhs7 + DN(1,1);
const double clhs157 = clhs156*clhs44;
const double clhs158 = DN(1,0)*clhs3;
const double clhs159 = DN(1,1)*clhs22;
const double clhs160 = DN(1,1)*clhs3;
const double clhs161 = DN(1,0)*clhs22;
const double clhs162 = clhs160 + clhs161;
const double clhs163 = C(0,0)*clhs158 + C(0,1)*clhs159 + C(0,2)*clhs162;
const double clhs164 = clhs163*clhs58;
const double clhs165 = C(0,2)*clhs158 + C(1,2)*clhs159 + C(2,2)*clhs162;
const double clhs166 = clhs165*clhs65;
const double clhs167 = clhs156*clhs68;
const double clhs168 = -clhs156*clhs74;
const double clhs169 = clhs10*clhs168 + clhs158*clhs73 + clhs168*clhs4;
const double clhs170 = clhs159*clhs73 + clhs168*clhs17 + clhs168*clhs23;
const double clhs171 = 2*clhs160*clhs78 + 2*clhs161*clhs78 + 2*clhs168*clhs26 + 2*clhs168*clhs27;
const double clhs172 = C(0,2)*clhs169 + C(1,2)*clhs170 + C(2,2)*clhs171;
const double clhs173 = C(0,0)*clhs169 + C(0,1)*clhs170 + C(0,2)*clhs171;
const double clhs174 = clhs16*clhs172 + clhs173*clhs9;
const double clhs175 = C(0,1)*clhs169 + C(1,1)*clhs170 + C(1,2)*clhs171;
const double clhs176 = clhs16*clhs175 + clhs172*clhs9;
const double clhs177 = -DN(0,0)*DN(2,1)*u(0,1) - DN(1,0)*DN(2,1)*u(1,1) + DN(2,0)*clhs18 + DN(2,0)*clhs19 + DN(2,0);
const double clhs178 = clhs177*clhs44;
const double clhs179 = DN(2,1)*clhs16;
const double clhs180 = DN(2,0)*clhs9;
const double clhs181 = DN(2,0)*clhs16;
const double clhs182 = DN(2,1)*clhs9;
const double clhs183 = clhs181 + clhs182;
const double clhs184 = C(0,0)*clhs180 + C(0,1)*clhs179 + C(0,2)*clhs183;
const double clhs185 = clhs184*clhs58;
const double clhs186 = C(0,2)*clhs180 + C(1,2)*clhs179 + C(2,2)*clhs183;
const double clhs187 = clhs186*clhs65;
const double clhs188 = clhs177*clhs68;
const double clhs189 = -clhs177*clhs74;
const double clhs190 = clhs17*clhs189 + clhs179*clhs73 + clhs189*clhs23;
const double clhs191 = clhs10*clhs189 + clhs180*clhs73 + clhs189*clhs4;
const double clhs192 = 2*clhs181*clhs78 + 2*clhs182*clhs78 + 2*clhs189*clhs26 + 2*clhs189*clhs27;
const double clhs193 = C(0,0)*clhs191 + C(0,1)*clhs190 + C(0,2)*clhs192;
const double clhs194 = C(0,2)*clhs191 + C(1,2)*clhs190 + C(2,2)*clhs192;
const double clhs195 = DN(2,0)*clhs90 + DN(2,1)*clhs91;
const double clhs196 = clhs193*clhs83 + clhs194*clhs85 + clhs195;
const double clhs197 = C(0,1)*clhs191 + C(1,1)*clhs190 + C(1,2)*clhs192;
const double clhs198 = DN(2,0)*clhs91 + DN(2,1)*clhs95;
const double clhs199 = clhs194*clhs83 + clhs197*clhs85 + clhs198;
const double clhs200 = -DN(0,1)*DN(2,0)*u(0,0) - DN(1,1)*DN(2,0)*u(1,0) + DN(2,1)*clhs5 + DN(2,1)*clhs6 + DN(2,1);
const double clhs201 = clhs200*clhs44;
const double clhs202 = DN(2,0)*clhs3;
const double clhs203 = DN(2,1)*clhs22;
const double clhs204 = DN(2,1)*clhs3;
const double clhs205 = DN(2,0)*clhs22;
const double clhs206 = clhs204 + clhs205;
const double clhs207 = C(0,0)*clhs202 + C(0,1)*clhs203 + C(0,2)*clhs206;
const double clhs208 = clhs207*clhs58;
const double clhs209 = C(0,2)*clhs202 + C(1,2)*clhs203 + C(2,2)*clhs206;
const double clhs210 = clhs209*clhs65;
const double clhs211 = clhs200*clhs68;
const double clhs212 = -clhs200*clhs74;
const double clhs213 = clhs10*clhs212 + clhs202*clhs73 + clhs212*clhs4;
const double clhs214 = clhs17*clhs212 + clhs203*clhs73 + clhs212*clhs23;
const double clhs215 = 2*clhs204*clhs78 + 2*clhs205*clhs78 + 2*clhs212*clhs26 + 2*clhs212*clhs27;
const double clhs216 = C(0,2)*clhs213 + C(1,2)*clhs214 + C(2,2)*clhs215;
const double clhs217 = C(0,0)*clhs213 + C(0,1)*clhs214 + C(0,2)*clhs215;
const double clhs218 = clhs16*clhs216 + clhs217*clhs9;
const double clhs219 = C(0,1)*clhs213 + C(1,1)*clhs214 + C(1,2)*clhs215;
const double clhs220 = clhs16*clhs219 + clhs216*clhs9;
const double clhs221 = DN(0,0)*clhs48;
const double clhs222 = clhs221*clhs50;
const double clhs223 = C(1,2)*clhs29 + clhs128 + clhs129;
const double clhs224 = DN(0,1)*clhs223;
const double clhs225 = clhs224*clhs50;
const double clhs226 = clhs58*clhs64;
const double clhs227 = C(0,1)*clhs53 + C(1,1)*clhs52 + C(1,2)*clhs56;
const double clhs228 = clhs227*clhs65;
const double clhs229 = clhs62*clhs69;
const double clhs230 = clhs22*clhs84 + clhs3*clhs80;
const double clhs231 = clhs22*clhs94 + clhs3*clhs84;
const double clhs232 = clhs107*clhs58;
const double clhs233 = C(0,1)*clhs100 + C(1,1)*clhs101 + C(1,2)*clhs104;
const double clhs234 = clhs233*clhs65;
const double clhs235 = clhs109*clhs62;
const double clhs236 = clhs3*clhs82;
const double clhs237 = clhs22*clhs82;
const double clhs238 = clhs114*clhs237 + clhs115*clhs236 + clhs92;
const double clhs239 = clhs114*clhs236 + clhs118*clhs237 + clhs96;
const double clhs240 = clhs125*clhs22 + clhs126*clhs3;
const double clhs241 = clhs125*clhs3 + clhs130*clhs22;
const double clhs242 = clhs44*(DN(0,0)*clhs123*clhs240 + DN(0,1)*clhs123*clhs241 - clhs122*clhs221 - clhs122*clhs224 - clhs222 - clhs225);
const double clhs243 = clhs142*clhs58;
const double clhs244 = C(0,1)*clhs136 + C(1,1)*clhs135 + C(1,2)*clhs139;
const double clhs245 = clhs244*clhs65;
const double clhs246 = clhs144*clhs62;
const double clhs247 = clhs149*clhs3 + clhs150*clhs22;
const double clhs248 = clhs150*clhs3 + clhs153*clhs22;
const double clhs249 = clhs165*clhs58;
const double clhs250 = C(0,1)*clhs158 + C(1,1)*clhs159 + C(1,2)*clhs162;
const double clhs251 = clhs250*clhs65;
const double clhs252 = clhs167*clhs62;
const double clhs253 = clhs151 + clhs172*clhs237 + clhs173*clhs236;
const double clhs254 = clhs154 + clhs172*clhs236 + clhs175*clhs237;
const double clhs255 = clhs186*clhs58;
const double clhs256 = C(0,1)*clhs180 + C(1,1)*clhs179 + C(1,2)*clhs183;
const double clhs257 = clhs256*clhs65;
const double clhs258 = clhs188*clhs62;
const double clhs259 = clhs193*clhs3 + clhs194*clhs22;
const double clhs260 = clhs194*clhs3 + clhs197*clhs22;
const double clhs261 = clhs209*clhs58;
const double clhs262 = C(0,1)*clhs202 + C(1,1)*clhs203 + C(1,2)*clhs206;
const double clhs263 = clhs262*clhs65;
const double clhs264 = clhs211*clhs62;
const double clhs265 = clhs195 + clhs216*clhs237 + clhs217*clhs236;
const double clhs266 = clhs198 + clhs216*clhs236 + clhs219*clhs237;
const double clhs267 = 0.1111111111111111/clhs43;
const double clhs268 = clhs267*clhs41;
const double clhs269 = DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2];
const double clhs270 = 1.0/clhs81;
const double clhs271 = pow(clhs270*clhs43, 0.33333333333333331);
const double clhs272 = clhs271*tau_u;
const double clhs273 = clhs269*clhs272;
const double clhs274 = DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2];
const double clhs275 = clhs272*clhs274;
const double clhs276 = clhs267*clhs98;
const double clhs277 = clhs31 + clhs49;
const double clhs278 = 0.33333333333333331*clhs272;
const double clhs279 = clhs277*clhs278;
const double clhs280 = clhs221 + clhs224;
const double clhs281 = clhs278*clhs280;
const double clhs282 = tau_u*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs283 = b_gauss[1]*clhs282;
const double clhs284 = clhs133*clhs267;
const double clhs285 = b_gauss[0]*clhs282;
const double clhs286 = clhs156*clhs267;
const double clhs287 = N[0]*N[1];
const double clhs288 = tau_u*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs289 = b_gauss[1]*clhs288;
const double clhs290 = clhs177*clhs267;
const double clhs291 = b_gauss[0]*clhs288;
const double clhs292 = clhs200*clhs267;
const double clhs293 = N[0]*N[2];
const double clhs294 = DN(1,0)*clhs30;
const double clhs295 = clhs294*clhs50;
const double clhs296 = DN(1,1)*clhs48;
const double clhs297 = clhs296*clhs50;
const double clhs298 = 0.66666666666666663*DN(1,0);
const double clhs299 = clhs298*clhs57;
const double clhs300 = 0.66666666666666663*DN(1,1);
const double clhs301 = clhs300*clhs64;
const double clhs302 = clhs105*clhs298;
const double clhs303 = clhs107*clhs300;
const double clhs304 = DN(1,0)*clhs82;
const double clhs305 = DN(1,1)*clhs82;
const double clhs306 = clhs44*(DN(1,0)*clhs123*clhs127 + DN(1,1)*clhs123*clhs131 - clhs122*clhs294 - clhs122*clhs296 - clhs295 - clhs297);
const double clhs307 = clhs140*clhs298;
const double clhs308 = clhs142*clhs300;
const double clhs309 = clhs163*clhs298;
const double clhs310 = clhs165*clhs300;
const double clhs311 = clhs184*clhs298;
const double clhs312 = clhs186*clhs300;
const double clhs313 = clhs207*clhs298;
const double clhs314 = clhs209*clhs300;
const double clhs315 = DN(1,0)*clhs48;
const double clhs316 = clhs315*clhs50;
const double clhs317 = DN(1,1)*clhs223;
const double clhs318 = clhs317*clhs50;
const double clhs319 = clhs298*clhs64;
const double clhs320 = clhs227*clhs300;
const double clhs321 = clhs107*clhs298;
const double clhs322 = clhs233*clhs300;
const double clhs323 = clhs44*(DN(1,0)*clhs123*clhs240 + DN(1,1)*clhs123*clhs241 - clhs122*clhs315 - clhs122*clhs317 - clhs316 - clhs318);
const double clhs324 = clhs142*clhs298;
const double clhs325 = clhs244*clhs300;
const double clhs326 = clhs165*clhs298;
const double clhs327 = clhs250*clhs300;
const double clhs328 = clhs186*clhs298;
const double clhs329 = clhs256*clhs300;
const double clhs330 = clhs209*clhs298;
const double clhs331 = clhs262*clhs300;
const double clhs332 = clhs294 + clhs296;
const double clhs333 = clhs278*clhs332;
const double clhs334 = clhs315 + clhs317;
const double clhs335 = clhs278*clhs334;
const double clhs336 = tau_u*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs337 = b_gauss[1]*clhs336;
const double clhs338 = b_gauss[0]*clhs336;
const double clhs339 = N[1]*N[2];
const double clhs340 = DN(2,0)*clhs30;
const double clhs341 = clhs340*clhs50;
const double clhs342 = DN(2,1)*clhs48;
const double clhs343 = clhs342*clhs50;
const double clhs344 = 0.66666666666666663*DN(2,0);
const double clhs345 = clhs344*clhs57;
const double clhs346 = 0.66666666666666663*DN(2,1);
const double clhs347 = clhs346*clhs64;
const double clhs348 = clhs105*clhs344;
const double clhs349 = clhs107*clhs346;
const double clhs350 = DN(2,0)*clhs82;
const double clhs351 = DN(2,1)*clhs82;
const double clhs352 = clhs44*(DN(2,0)*clhs123*clhs127 + DN(2,1)*clhs123*clhs131 - clhs122*clhs340 - clhs122*clhs342 - clhs341 - clhs343);
const double clhs353 = clhs140*clhs344;
const double clhs354 = clhs142*clhs346;
const double clhs355 = clhs163*clhs344;
const double clhs356 = clhs165*clhs346;
const double clhs357 = clhs184*clhs344;
const double clhs358 = clhs186*clhs346;
const double clhs359 = clhs207*clhs344;
const double clhs360 = clhs209*clhs346;
const double clhs361 = DN(2,0)*clhs48;
const double clhs362 = clhs361*clhs50;
const double clhs363 = DN(2,1)*clhs223;
const double clhs364 = clhs363*clhs50;
const double clhs365 = clhs344*clhs64;
const double clhs366 = clhs227*clhs346;
const double clhs367 = clhs107*clhs344;
const double clhs368 = clhs233*clhs346;
const double clhs369 = clhs44*(DN(2,0)*clhs123*clhs240 + DN(2,1)*clhs123*clhs241 - clhs122*clhs361 - clhs122*clhs363 - clhs362 - clhs364);
const double clhs370 = clhs142*clhs344;
const double clhs371 = clhs244*clhs346;
const double clhs372 = clhs165*clhs344;
const double clhs373 = clhs250*clhs346;
const double clhs374 = clhs186*clhs344;
const double clhs375 = clhs256*clhs346;
const double clhs376 = clhs209*clhs344;
const double clhs377 = clhs262*clhs346;
const double clhs378 = clhs340 + clhs342;
const double clhs379 = clhs278*clhs378;
const double clhs380 = clhs361 + clhs363;
const double clhs381 = clhs278*clhs380;
lhs(0,0)=DN(0,0)*clhs93 + DN(0,1)*clhs97 + clhs40*clhs45 + clhs45*clhs51 + clhs59*clhs63 + clhs63*clhs66 - clhs69*clhs71 - clhs69*clhs72;
lhs(0,1)=clhs106*clhs63 + clhs108*clhs63 - clhs109*clhs71 - clhs109*clhs72 + clhs116*clhs117 + clhs119*clhs120 + clhs40*clhs99 + clhs51*clhs99;
lhs(0,2)=N[0]*clhs132;
lhs(0,3)=DN(0,0)*clhs152 + DN(0,1)*clhs155 + clhs134*clhs40 + clhs134*clhs51 + clhs141*clhs63 + clhs143*clhs63 - clhs144*clhs71 - clhs144*clhs72;
lhs(0,4)=clhs117*clhs174 + clhs120*clhs176 + clhs157*clhs40 + clhs157*clhs51 + clhs164*clhs63 + clhs166*clhs63 - clhs167*clhs71 - clhs167*clhs72;
lhs(0,5)=N[1]*clhs132;
lhs(0,6)=DN(0,0)*clhs196 + DN(0,1)*clhs199 + clhs178*clhs40 + clhs178*clhs51 + clhs185*clhs63 + clhs187*clhs63 - clhs188*clhs71 - clhs188*clhs72;
lhs(0,7)=clhs117*clhs218 + clhs120*clhs220 + clhs201*clhs40 + clhs201*clhs51 + clhs208*clhs63 + clhs210*clhs63 - clhs211*clhs71 - clhs211*clhs72;
lhs(0,8)=N[2]*clhs132;
lhs(1,0)=clhs117*clhs230 + clhs120*clhs231 - clhs221*clhs229 + clhs222*clhs45 - clhs224*clhs229 + clhs225*clhs45 + clhs226*clhs63 + clhs228*clhs63;
lhs(1,1)=DN(0,0)*clhs238 + DN(0,1)*clhs239 - clhs221*clhs235 + clhs222*clhs99 - clhs224*clhs235 + clhs225*clhs99 + clhs232*clhs63 + clhs234*clhs63;
lhs(1,2)=N[0]*clhs242;
lhs(1,3)=clhs117*clhs247 + clhs120*clhs248 + clhs134*clhs222 + clhs134*clhs225 - clhs221*clhs246 - clhs224*clhs246 + clhs243*clhs63 + clhs245*clhs63;
lhs(1,4)=DN(0,0)*clhs253 + DN(0,1)*clhs254 + clhs157*clhs222 + clhs157*clhs225 - clhs221*clhs252 - clhs224*clhs252 + clhs249*clhs63 + clhs251*clhs63;
lhs(1,5)=N[1]*clhs242;
lhs(1,6)=clhs117*clhs259 + clhs120*clhs260 + clhs178*clhs222 + clhs178*clhs225 - clhs221*clhs258 - clhs224*clhs258 + clhs255*clhs63 + clhs257*clhs63;
lhs(1,7)=DN(0,0)*clhs265 + DN(0,1)*clhs266 + clhs201*clhs222 + clhs201*clhs225 - clhs221*clhs264 - clhs224*clhs264 + clhs261*clhs63 + clhs263*clhs63;
lhs(1,8)=N[2]*clhs242;
lhs(2,0)=N[0]*clhs41 - clhs273*(clhs268*clhs31 + clhs268*clhs49 + clhs59 + clhs66) - clhs275*(clhs221*clhs268 + clhs224*clhs268 + clhs226 + clhs228);
lhs(2,1)=N[0]*clhs98 - clhs273*(clhs106 + clhs108 + clhs276*clhs31 + clhs276*clhs49) - clhs275*(clhs221*clhs276 + clhs224*clhs276 + clhs232 + clhs234);
lhs(2,2)=-DN(0,0)*clhs279 - DN(0,1)*clhs281 - pow(N[0], 2) + 0.1111111111111111*N[0]*clhs269*clhs270*clhs271*clhs277*tau_u + 0.1111111111111111*N[0]*clhs270*clhs271*clhs274*clhs280*tau_u;
lhs(2,3)=N[0]*clhs133 - clhs273*(clhs141 + clhs143 + clhs284*clhs31 + clhs284*clhs49) - clhs275*(clhs221*clhs284 + clhs224*clhs284 + clhs243 + clhs245) + clhs283;
lhs(2,4)=N[0]*clhs156 - clhs273*(clhs164 + clhs166 + clhs286*clhs31 + clhs286*clhs49) - clhs275*(clhs221*clhs286 + clhs224*clhs286 + clhs249 + clhs251) - clhs285;
lhs(2,5)=-DN(1,0)*clhs279 - DN(1,1)*clhs281 + 0.1111111111111111*N[1]*clhs269*clhs270*clhs271*clhs277*tau_u + 0.1111111111111111*N[1]*clhs270*clhs271*clhs274*clhs280*tau_u - clhs287;
lhs(2,6)=N[0]*clhs177 - clhs273*(clhs185 + clhs187 + clhs290*clhs31 + clhs290*clhs49) - clhs275*(clhs221*clhs290 + clhs224*clhs290 + clhs255 + clhs257) + clhs289;
lhs(2,7)=N[0]*clhs200 - clhs273*(clhs208 + clhs210 + clhs292*clhs31 + clhs292*clhs49) - clhs275*(clhs221*clhs292 + clhs224*clhs292 + clhs261 + clhs263) - clhs291;
lhs(2,8)=-DN(2,0)*clhs279 - DN(2,1)*clhs281 + 0.1111111111111111*N[2]*clhs269*clhs270*clhs271*clhs277*tau_u + 0.1111111111111111*N[2]*clhs270*clhs271*clhs274*clhs280*tau_u - clhs293;
lhs(3,0)=DN(1,0)*clhs93 + DN(1,1)*clhs97 - clhs229*clhs294 - clhs229*clhs296 + clhs295*clhs45 + clhs297*clhs45 + clhs299*clhs63 + clhs301*clhs63;
lhs(3,1)=clhs116*clhs304 + clhs119*clhs305 - clhs235*clhs294 - clhs235*clhs296 + clhs295*clhs99 + clhs297*clhs99 + clhs302*clhs63 + clhs303*clhs63;
lhs(3,2)=N[0]*clhs306;
lhs(3,3)=DN(1,0)*clhs152 + DN(1,1)*clhs155 + clhs134*clhs295 + clhs134*clhs297 - clhs246*clhs294 - clhs246*clhs296 + clhs307*clhs63 + clhs308*clhs63;
lhs(3,4)=clhs157*clhs295 + clhs157*clhs297 + clhs174*clhs304 + clhs176*clhs305 - clhs252*clhs294 - clhs252*clhs296 + clhs309*clhs63 + clhs310*clhs63;
lhs(3,5)=N[1]*clhs306;
lhs(3,6)=DN(1,0)*clhs196 + DN(1,1)*clhs199 + clhs178*clhs295 + clhs178*clhs297 - clhs258*clhs294 - clhs258*clhs296 + clhs311*clhs63 + clhs312*clhs63;
lhs(3,7)=clhs201*clhs295 + clhs201*clhs297 + clhs218*clhs304 + clhs220*clhs305 - clhs264*clhs294 - clhs264*clhs296 + clhs313*clhs63 + clhs314*clhs63;
lhs(3,8)=N[2]*clhs306;
lhs(4,0)=-clhs229*clhs315 - clhs229*clhs317 + clhs230*clhs304 + clhs231*clhs305 + clhs316*clhs45 + clhs318*clhs45 + clhs319*clhs63 + clhs320*clhs63;
lhs(4,1)=DN(1,0)*clhs238 + DN(1,1)*clhs239 - clhs235*clhs315 - clhs235*clhs317 + clhs316*clhs99 + clhs318*clhs99 + clhs321*clhs63 + clhs322*clhs63;
lhs(4,2)=N[0]*clhs323;
lhs(4,3)=clhs134*clhs316 + clhs134*clhs318 - clhs246*clhs315 - clhs246*clhs317 + clhs247*clhs304 + clhs248*clhs305 + clhs324*clhs63 + clhs325*clhs63;
lhs(4,4)=DN(1,0)*clhs253 + DN(1,1)*clhs254 + clhs157*clhs316 + clhs157*clhs318 - clhs252*clhs315 - clhs252*clhs317 + clhs326*clhs63 + clhs327*clhs63;
lhs(4,5)=N[1]*clhs323;
lhs(4,6)=clhs178*clhs316 + clhs178*clhs318 - clhs258*clhs315 - clhs258*clhs317 + clhs259*clhs304 + clhs260*clhs305 + clhs328*clhs63 + clhs329*clhs63;
lhs(4,7)=DN(1,0)*clhs265 + DN(1,1)*clhs266 + clhs201*clhs316 + clhs201*clhs318 - clhs264*clhs315 - clhs264*clhs317 + clhs330*clhs63 + clhs331*clhs63;
lhs(4,8)=N[2]*clhs323;
lhs(5,0)=N[1]*clhs41 - clhs273*(clhs268*clhs294 + clhs268*clhs296 + clhs299 + clhs301) - clhs275*(clhs268*clhs315 + clhs268*clhs317 + clhs319 + clhs320) - clhs283;
lhs(5,1)=N[1]*clhs98 - clhs273*(clhs276*clhs294 + clhs276*clhs296 + clhs302 + clhs303) - clhs275*(clhs276*clhs315 + clhs276*clhs317 + clhs321 + clhs322) + clhs285;
lhs(5,2)=-DN(0,0)*clhs333 - DN(0,1)*clhs335 + 0.1111111111111111*N[0]*clhs269*clhs270*clhs271*clhs332*tau_u + 0.1111111111111111*N[0]*clhs270*clhs271*clhs274*clhs334*tau_u - clhs287;
lhs(5,3)=N[1]*clhs133 - clhs273*(clhs284*clhs294 + clhs284*clhs296 + clhs307 + clhs308) - clhs275*(clhs284*clhs315 + clhs284*clhs317 + clhs324 + clhs325);
lhs(5,4)=N[1]*clhs156 - clhs273*(clhs286*clhs294 + clhs286*clhs296 + clhs309 + clhs310) - clhs275*(clhs286*clhs315 + clhs286*clhs317 + clhs326 + clhs327);
lhs(5,5)=-DN(1,0)*clhs333 - DN(1,1)*clhs335 - pow(N[1], 2) + 0.1111111111111111*N[1]*clhs269*clhs270*clhs271*clhs332*tau_u + 0.1111111111111111*N[1]*clhs270*clhs271*clhs274*clhs334*tau_u;
lhs(5,6)=N[1]*clhs177 - clhs273*(clhs290*clhs294 + clhs290*clhs296 + clhs311 + clhs312) - clhs275*(clhs290*clhs315 + clhs290*clhs317 + clhs328 + clhs329) + clhs337;
lhs(5,7)=N[1]*clhs200 - clhs273*(clhs292*clhs294 + clhs292*clhs296 + clhs313 + clhs314) - clhs275*(clhs292*clhs315 + clhs292*clhs317 + clhs330 + clhs331) - clhs338;
lhs(5,8)=-DN(2,0)*clhs333 - DN(2,1)*clhs335 + 0.1111111111111111*N[2]*clhs269*clhs270*clhs271*clhs332*tau_u + 0.1111111111111111*N[2]*clhs270*clhs271*clhs274*clhs334*tau_u - clhs339;
lhs(6,0)=DN(2,0)*clhs93 + DN(2,1)*clhs97 - clhs229*clhs340 - clhs229*clhs342 + clhs341*clhs45 + clhs343*clhs45 + clhs345*clhs63 + clhs347*clhs63;
lhs(6,1)=clhs116*clhs350 + clhs119*clhs351 - clhs235*clhs340 - clhs235*clhs342 + clhs341*clhs99 + clhs343*clhs99 + clhs348*clhs63 + clhs349*clhs63;
lhs(6,2)=N[0]*clhs352;
lhs(6,3)=DN(2,0)*clhs152 + DN(2,1)*clhs155 + clhs134*clhs341 + clhs134*clhs343 - clhs246*clhs340 - clhs246*clhs342 + clhs353*clhs63 + clhs354*clhs63;
lhs(6,4)=clhs157*clhs341 + clhs157*clhs343 + clhs174*clhs350 + clhs176*clhs351 - clhs252*clhs340 - clhs252*clhs342 + clhs355*clhs63 + clhs356*clhs63;
lhs(6,5)=N[1]*clhs352;
lhs(6,6)=DN(2,0)*clhs196 + DN(2,1)*clhs199 + clhs178*clhs341 + clhs178*clhs343 - clhs258*clhs340 - clhs258*clhs342 + clhs357*clhs63 + clhs358*clhs63;
lhs(6,7)=clhs201*clhs341 + clhs201*clhs343 + clhs218*clhs350 + clhs220*clhs351 - clhs264*clhs340 - clhs264*clhs342 + clhs359*clhs63 + clhs360*clhs63;
lhs(6,8)=N[2]*clhs352;
lhs(7,0)=-clhs229*clhs361 - clhs229*clhs363 + clhs230*clhs350 + clhs231*clhs351 + clhs362*clhs45 + clhs364*clhs45 + clhs365*clhs63 + clhs366*clhs63;
lhs(7,1)=DN(2,0)*clhs238 + DN(2,1)*clhs239 - clhs235*clhs361 - clhs235*clhs363 + clhs362*clhs99 + clhs364*clhs99 + clhs367*clhs63 + clhs368*clhs63;
lhs(7,2)=N[0]*clhs369;
lhs(7,3)=clhs134*clhs362 + clhs134*clhs364 - clhs246*clhs361 - clhs246*clhs363 + clhs247*clhs350 + clhs248*clhs351 + clhs370*clhs63 + clhs371*clhs63;
lhs(7,4)=DN(2,0)*clhs253 + DN(2,1)*clhs254 + clhs157*clhs362 + clhs157*clhs364 - clhs252*clhs361 - clhs252*clhs363 + clhs372*clhs63 + clhs373*clhs63;
lhs(7,5)=N[1]*clhs369;
lhs(7,6)=clhs178*clhs362 + clhs178*clhs364 - clhs258*clhs361 - clhs258*clhs363 + clhs259*clhs350 + clhs260*clhs351 + clhs374*clhs63 + clhs375*clhs63;
lhs(7,7)=DN(2,0)*clhs265 + DN(2,1)*clhs266 + clhs201*clhs362 + clhs201*clhs364 - clhs264*clhs361 - clhs264*clhs363 + clhs376*clhs63 + clhs377*clhs63;
lhs(7,8)=N[2]*clhs369;
lhs(8,0)=N[2]*clhs41 - clhs273*(clhs268*clhs340 + clhs268*clhs342 + clhs345 + clhs347) - clhs275*(clhs268*clhs361 + clhs268*clhs363 + clhs365 + clhs366) - clhs289;
lhs(8,1)=N[2]*clhs98 - clhs273*(clhs276*clhs340 + clhs276*clhs342 + clhs348 + clhs349) - clhs275*(clhs276*clhs361 + clhs276*clhs363 + clhs367 + clhs368) + clhs291;
lhs(8,2)=-DN(0,0)*clhs379 - DN(0,1)*clhs381 + 0.1111111111111111*N[0]*clhs269*clhs270*clhs271*clhs378*tau_u + 0.1111111111111111*N[0]*clhs270*clhs271*clhs274*clhs380*tau_u - clhs293;
lhs(8,3)=N[2]*clhs133 - clhs273*(clhs284*clhs340 + clhs284*clhs342 + clhs353 + clhs354) - clhs275*(clhs284*clhs361 + clhs284*clhs363 + clhs370 + clhs371) - clhs337;
lhs(8,4)=N[2]*clhs156 - clhs273*(clhs286*clhs340 + clhs286*clhs342 + clhs355 + clhs356) - clhs275*(clhs286*clhs361 + clhs286*clhs363 + clhs372 + clhs373) + clhs338;
lhs(8,5)=-DN(1,0)*clhs379 - DN(1,1)*clhs381 + 0.1111111111111111*N[1]*clhs269*clhs270*clhs271*clhs378*tau_u + 0.1111111111111111*N[1]*clhs270*clhs271*clhs274*clhs380*tau_u - clhs339;
lhs(8,6)=N[2]*clhs177 - clhs273*(clhs290*clhs340 + clhs290*clhs342 + clhs357 + clhs358) - clhs275*(clhs290*clhs361 + clhs290*clhs363 + clhs374 + clhs375);
lhs(8,7)=N[2]*clhs200 - clhs273*(clhs292*clhs340 + clhs292*clhs342 + clhs359 + clhs360) - clhs275*(clhs292*clhs361 + clhs292*clhs363 + clhs376 + clhs377);
lhs(8,8)=-DN(2,0)*clhs379 - DN(2,1)*clhs381 - pow(N[2], 2) + 0.1111111111111111*N[2]*clhs269*clhs270*clhs271*clhs378*tau_u + 0.1111111111111111*N[2]*clhs270*clhs271*clhs274*clhs380*tau_u;

        // Calculate and add the RHS Gauss point contribution
        const double crhs0 = DN(0,1)*u(0,0);
const double crhs1 = DN(1,1)*u(1,0);
const double crhs2 = DN(2,1)*u(2,0);
const double crhs3 = crhs0 + crhs1 + crhs2;
const double crhs4 = DN(0,0)*u(0,0);
const double crhs5 = DN(1,0)*u(1,0);
const double crhs6 = DN(2,0)*u(2,0);
const double crhs7 = crhs4 + crhs5 + crhs6;
const double crhs8 = crhs7 + 1;
const double crhs9 = S[0]*crhs8 + S[2]*crhs3;
const double crhs10 = S[1]*crhs3 + S[2]*crhs8;
const double crhs11 = DN(0,0)*u(0,1);
const double crhs12 = DN(1,0)*u(1,1);
const double crhs13 = DN(2,0)*u(2,1);
const double crhs14 = crhs11 + crhs12 + crhs13;
const double crhs15 = pow(crhs14, 2) + pow(crhs8, 2);
const double crhs16 = DN(0,1)*u(0,1);
const double crhs17 = DN(1,1)*u(1,1);
const double crhs18 = DN(2,1)*u(2,1);
const double crhs19 = crhs16 + crhs17 + crhs18;
const double crhs20 = crhs19 + 1;
const double crhs21 = pow(crhs20, 2) + pow(crhs3, 2);
const double crhs22 = 2*crhs14*crhs20 + 2*crhs3*crhs8;
const double crhs23 = C(0,0)*crhs15 + C(0,1)*crhs21 + C(0,2)*crhs22;
const double crhs24 = DN(0,0)*crhs23;
const double crhs25 = N[0]*th[0];
const double crhs26 = N[1]*th[1];
const double crhs27 = N[2]*th[2];
const double crhs28 = -crhs0*crhs12 - crhs0*crhs13 - crhs1*crhs11 - crhs1*crhs13 - crhs11*crhs2 - crhs12*crhs2 + crhs16*crhs5 + crhs16*crhs6 + crhs17*crhs4 + crhs17*crhs6 + crhs18*crhs4 + crhs18*crhs5 + crhs19;
const double crhs29 = -crhs25 - crhs26 - crhs27 + crhs28 + crhs7;
const double crhs30 = crhs25 + crhs26 + crhs27;
const double crhs31 = crhs28 + crhs8;
const double crhs32 = 0.33333333333333331*crhs29*pow(crhs31, -0.66666666666666663)*tau_th*pow(crhs30 + 1, -0.33333333333333331);
const double crhs33 = C(0,2)*crhs15 + C(1,2)*crhs21 + C(2,2)*crhs22;
const double crhs34 = DN(0,1)*crhs33;
const double crhs35 = S[0]*crhs14 + S[2]*crhs20;
const double crhs36 = S[1]*crhs20 + S[2]*crhs14;
const double crhs37 = DN(0,0)*crhs33;
const double crhs38 = C(0,1)*crhs15 + C(1,1)*crhs21 + C(1,2)*crhs22;
const double crhs39 = DN(0,1)*crhs38;
const double crhs40 = b_gauss[0]*tau_u;
const double crhs41 = b_gauss[1]*tau_u;
const double crhs42 = 0.33333333333333331*tau_u*pow(crhs31/(crhs30 + 1.0), 0.33333333333333331);
const double crhs43 = crhs42*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs44 = crhs42*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double crhs45 = DN(1,0)*crhs23;
const double crhs46 = DN(1,1)*crhs33;
const double crhs47 = DN(1,0)*crhs33;
const double crhs48 = DN(1,1)*crhs38;
const double crhs49 = DN(2,0)*crhs23;
const double crhs50 = DN(2,1)*crhs33;
const double crhs51 = DN(2,0)*crhs33;
const double crhs52 = DN(2,1)*crhs38;
rhs[0]=-DN(0,0)*crhs9 - DN(0,1)*crhs10 + N[0]*b_gauss[0] - crhs24*crhs32 - crhs32*crhs34;
rhs[1]=-DN(0,0)*crhs35 - DN(0,1)*crhs36 + N[0]*b_gauss[1] - crhs32*crhs37 - crhs32*crhs39;
rhs[2]=-N[0]*crhs29 + crhs40*(DN(0,0)*crhs20 - DN(0,1)*crhs14) - crhs41*(DN(0,0)*crhs3 - DN(0,1)*crhs8) + crhs43*(crhs24 + crhs34) + crhs44*(crhs37 + crhs39);
rhs[3]=-DN(1,0)*crhs9 - DN(1,1)*crhs10 + N[1]*b_gauss[0] - crhs32*crhs45 - crhs32*crhs46;
rhs[4]=-DN(1,0)*crhs35 - DN(1,1)*crhs36 + N[1]*b_gauss[1] - crhs32*crhs47 - crhs32*crhs48;
rhs[5]=-N[1]*crhs29 + crhs40*(DN(1,0)*crhs20 - DN(1,1)*crhs14) - crhs41*(DN(1,0)*crhs3 - DN(1,1)*crhs8) + crhs43*(crhs45 + crhs46) + crhs44*(crhs47 + crhs48);
rhs[6]=-DN(2,0)*crhs9 - DN(2,1)*crhs10 + N[2]*b_gauss[0] - crhs32*crhs49 - crhs32*crhs50;
rhs[7]=-DN(2,0)*crhs35 - DN(2,1)*crhs36 + N[2]*b_gauss[1] - crhs32*crhs51 - crhs32*crhs52;
rhs[8]=-N[2]*crhs29 + crhs40*(DN(2,0)*crhs20 - DN(2,1)*crhs14) - crhs41*(DN(2,0)*crhs3 - DN(2,1)*crhs8) + crhs43*(crhs49 + crhs50) + crhs44*(crhs51 + crhs52);

        //TODO: Amend this once the assembly is done in the input arrays
        rLeftHandSideMatrix += w_gauss * lhs;
        rRightHandSideVector += w_gauss * rhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();

    // Check RHS size
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < 3; ++d) {
            kinematic_variables.Displacements(i_node * 3 + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    // Calculate stabilization constant
    const double c_tau_u = 2.0;
    const double c_tau_th = 4.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau_u * std::pow(h,2) / 2.0;

    // Set the auxiliary references matching the automatic differentiation symbols
    array_1d<double,3> b_gauss;
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedVector<double,LocalSize> rhs;
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_PK2);

        // Calculate body force
        // Note that this already includes the density computed in the reference configuration
        b_gauss = StructuralMechanicsElementUtilities::GetBodyForce(*this, r_integration_points, i_gauss);

        // Calculate the stabilization constant
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
        double kappa = CalculateBulkModulus(constitutive_variables.ConstitutiveMatrix); // Equivalent bulk modulus
        const double tau_u = aux_tau / mu;
        const double tau_th = (c_tau_th * mu) / (mu + kappa);

        // Calculate and add the LHS Gauss point contributions
        //substitute_lhs_3D_4N
        // Calculate and add the RHS Gauss point contribution
        //substitute_rhs_3D_4N
        //TODO: Amend this once the assembly is done in the input arrays
        rLeftHandSideMatrix += w_gauss * lhs;
        rRightHandSideVector += w_gauss * rhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < 2; ++d) {
            kinematic_variables.Displacements(i_node * 2 + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();

    // Calculate stabilization constant
    const double c_tau_u = 2.0;
    const double c_tau_th = 4.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau_u * std::pow(h,2) / 2.0;

    // Set the auxiliary references matching the automatic differentiation symbols
    array_1d<double,3> b_gauss;
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_PK2);

        // Calculate body force
        // Note that this already includes the density computed in the reference configuration
        b_gauss = StructuralMechanicsElementUtilities::GetBodyForce(*this, r_integration_points, i_gauss);

        // Calculate the stabilization constant
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
        double kappa = CalculateBulkModulus(constitutive_variables.ConstitutiveMatrix); // Equivalent bulk modulus
        const double tau_u = aux_tau / mu;
        const double tau_th = (c_tau_th * mu) / (mu + kappa);

        // Calculate and add the LHS Gauss point contributions
        const double clhs0 = DN(0,0)*u(0,1);
const double clhs1 = DN(1,0)*u(1,1);
const double clhs2 = DN(2,0)*u(2,1);
const double clhs3 = clhs0 + clhs1 + clhs2;
const double clhs4 = pow(clhs3, 2);
const double clhs5 = DN(0,0)*u(0,0);
const double clhs6 = DN(1,0)*u(1,0);
const double clhs7 = DN(2,0)*u(2,0);
const double clhs8 = clhs5 + clhs6 + clhs7;
const double clhs9 = clhs8 + 1;
const double clhs10 = pow(clhs9, 2);
const double clhs11 = clhs10 + clhs4;
const double clhs12 = C(0,0)*clhs11;
const double clhs13 = DN(0,1)*u(0,0);
const double clhs14 = DN(1,1)*u(1,0);
const double clhs15 = DN(2,1)*u(2,0);
const double clhs16 = clhs13 + clhs14 + clhs15;
const double clhs17 = pow(clhs16, 2);
const double clhs18 = DN(0,1)*u(0,1);
const double clhs19 = DN(1,1)*u(1,1);
const double clhs20 = DN(2,1)*u(2,1);
const double clhs21 = clhs18 + clhs19 + clhs20;
const double clhs22 = clhs21 + 1;
const double clhs23 = pow(clhs22, 2);
const double clhs24 = clhs17 + clhs23;
const double clhs25 = C(0,1)*clhs24;
const double clhs26 = clhs22*clhs3;
const double clhs27 = clhs16*clhs9;
const double clhs28 = clhs26 + clhs27;
const double clhs29 = 2*clhs28;
const double clhs30 = C(0,2)*clhs29 + clhs12 + clhs25;
const double clhs31 = DN(0,0)*clhs30;
const double clhs32 = clhs31*tau_th;
const double clhs33 = N[0]*th[0];
const double clhs34 = N[1]*th[1];
const double clhs35 = N[2]*th[2];
const double clhs36 = clhs33 + clhs34 + clhs35;
const double clhs37 = clhs36 + 1;
const double clhs38 = pow(clhs37, -0.33333333333333331);
const double clhs39 = 0.33333333333333331*clhs38;
const double clhs40 = clhs32*clhs39;
const double clhs41 = DN(0,0)*clhs19 + DN(0,0)*clhs20 + DN(0,0) - DN(0,1)*DN(1,0)*u(1,1) - DN(0,1)*DN(2,0)*u(2,1);
const double clhs42 = -clhs0*clhs14 - clhs0*clhs15 - clhs1*clhs13 - clhs1*clhs15 - clhs13*clhs2 - clhs14*clhs2 + clhs18*clhs6 + clhs18*clhs7 + clhs19*clhs5 + clhs19*clhs7 + clhs20*clhs5 + clhs20*clhs6 + clhs21;
const double clhs43 = clhs42 + clhs9;
const double clhs44 = pow(clhs43, -0.66666666666666663);
const double clhs45 = clhs41*clhs44;
const double clhs46 = C(0,2)*clhs11;
const double clhs47 = C(1,2)*clhs24;
const double clhs48 = C(2,2)*clhs29 + clhs46 + clhs47;
const double clhs49 = DN(0,1)*clhs48;
const double clhs50 = clhs39*tau_th;
const double clhs51 = clhs49*clhs50;
const double clhs52 = DN(0,1)*clhs16;
const double clhs53 = DN(0,0)*clhs9;
const double clhs54 = DN(0,0)*clhs16;
const double clhs55 = DN(0,1)*clhs9;
const double clhs56 = clhs54 + clhs55;
const double clhs57 = C(0,0)*clhs53 + C(0,1)*clhs52 + C(0,2)*clhs56;
const double clhs58 = 0.66666666666666663*DN(0,0);
const double clhs59 = clhs57*clhs58;
const double clhs60 = -clhs33 - clhs34 - clhs35 + clhs42 + clhs8;
const double clhs61 = clhs60*tau_th;
const double clhs62 = clhs38*clhs61;
const double clhs63 = clhs44*clhs62;
const double clhs64 = C(0,2)*clhs53 + C(1,2)*clhs52 + C(2,2)*clhs56;
const double clhs65 = 0.66666666666666663*DN(0,1);
const double clhs66 = clhs64*clhs65;
const double clhs67 = pow(clhs43, -1.6666666666666665);
const double clhs68 = 0.22222222222222221*clhs67;
const double clhs69 = clhs41*clhs68;
const double clhs70 = clhs32*clhs60;
const double clhs71 = clhs38*clhs70;
const double clhs72 = clhs49*clhs62;
const double clhs73 = 1.0*clhs44;
const double clhs74 = 0.33333333333333331*clhs67;
const double clhs75 = -clhs41*clhs74;
const double clhs76 = clhs17*clhs75 + clhs23*clhs75 + clhs52*clhs73;
const double clhs77 = clhs10*clhs75 + clhs4*clhs75 + clhs53*clhs73;
const double clhs78 = 0.5*clhs44;
const double clhs79 = 2*clhs26*clhs75 + 2*clhs27*clhs75 + 2*clhs54*clhs78 + 2*clhs55*clhs78;
const double clhs80 = C(0,0)*clhs77 + C(0,1)*clhs76 + C(0,2)*clhs79;
const double clhs81 = clhs36 + 1.0;
const double clhs82 = pow(clhs81, 0.66666666666666663);
const double clhs83 = clhs82*clhs9;
const double clhs84 = C(0,2)*clhs77 + C(1,2)*clhs76 + C(2,2)*clhs79;
const double clhs85 = clhs16*clhs82;
const double clhs86 = clhs44*clhs82;
const double clhs87 = clhs28*clhs86;
const double clhs88 = 0.5*clhs10*clhs86 + 0.5*clhs4*clhs86 - 0.5;
const double clhs89 = 0.5*clhs17*clhs86 + 0.5*clhs23*clhs86 - 0.5;
const double clhs90 = C(0,0)*clhs88 + C(0,1)*clhs89 + C(0,2)*clhs87;
const double clhs91 = C(0,2)*clhs88 + C(1,2)*clhs89 + C(2,2)*clhs87;
const double clhs92 = DN(0,0)*clhs90 + DN(0,1)*clhs91;
const double clhs93 = clhs80*clhs83 + clhs84*clhs85 + clhs92;
const double clhs94 = C(0,1)*clhs77 + C(1,1)*clhs76 + C(1,2)*clhs79;
const double clhs95 = C(0,1)*clhs88 + C(1,1)*clhs89 + C(1,2)*clhs87;
const double clhs96 = DN(0,0)*clhs91 + DN(0,1)*clhs95;
const double clhs97 = clhs83*clhs84 + clhs85*clhs94 + clhs96;
const double clhs98 = -DN(0,0)*DN(1,1)*u(1,0) - DN(0,0)*DN(2,1)*u(2,0) + DN(0,1)*clhs6 + DN(0,1)*clhs7 + DN(0,1);
const double clhs99 = clhs44*clhs98;
const double clhs100 = DN(0,0)*clhs3;
const double clhs101 = DN(0,1)*clhs22;
const double clhs102 = DN(0,1)*clhs3;
const double clhs103 = DN(0,0)*clhs22;
const double clhs104 = clhs102 + clhs103;
const double clhs105 = C(0,0)*clhs100 + C(0,1)*clhs101 + C(0,2)*clhs104;
const double clhs106 = clhs105*clhs58;
const double clhs107 = C(0,2)*clhs100 + C(1,2)*clhs101 + C(2,2)*clhs104;
const double clhs108 = clhs107*clhs65;
const double clhs109 = clhs68*clhs98;
const double clhs110 = -clhs74*clhs98;
const double clhs111 = clhs10*clhs110 + clhs100*clhs73 + clhs110*clhs4;
const double clhs112 = clhs101*clhs73 + clhs110*clhs17 + clhs110*clhs23;
const double clhs113 = 2*clhs102*clhs78 + 2*clhs103*clhs78 + 2*clhs110*clhs26 + 2*clhs110*clhs27;
const double clhs114 = C(0,2)*clhs111 + C(1,2)*clhs112 + C(2,2)*clhs113;
const double clhs115 = C(0,0)*clhs111 + C(0,1)*clhs112 + C(0,2)*clhs113;
const double clhs116 = clhs114*clhs16 + clhs115*clhs9;
const double clhs117 = DN(0,0)*clhs82;
const double clhs118 = C(0,1)*clhs111 + C(1,1)*clhs112 + C(1,2)*clhs113;
const double clhs119 = clhs114*clhs9 + clhs118*clhs16;
const double clhs120 = DN(0,1)*clhs82;
const double clhs121 = 0.1111111111111111*pow(clhs37, -1.3333333333333333);
const double clhs122 = clhs121*clhs61;
const double clhs123 = pow(clhs81, -0.33333333333333337);
const double clhs124 = 0.66666666666666663*clhs28;
const double clhs125 = C(2,2)*clhs124 + 0.33333333333333331*clhs46 + 0.33333333333333331*clhs47;
const double clhs126 = C(0,2)*clhs124 + 0.33333333333333331*clhs12 + 0.33333333333333331*clhs25;
const double clhs127 = clhs125*clhs16 + clhs126*clhs9;
const double clhs128 = C(0,1)*clhs11;
const double clhs129 = C(1,1)*clhs24;
const double clhs130 = C(1,2)*clhs124 + 0.33333333333333331*clhs128 + 0.33333333333333331*clhs129;
const double clhs131 = clhs125*clhs9 + clhs130*clhs16;
const double clhs132 = clhs44*(DN(0,0)*clhs123*clhs127 + DN(0,1)*clhs123*clhs131 - clhs121*clhs70 - clhs122*clhs49 - clhs40 - clhs51);
const double clhs133 = -DN(0,0)*DN(1,1)*u(0,1) + DN(1,0)*clhs18 + DN(1,0)*clhs20 + DN(1,0) - DN(1,1)*DN(2,0)*u(2,1);
const double clhs134 = clhs133*clhs44;
const double clhs135 = DN(1,1)*clhs16;
const double clhs136 = DN(1,0)*clhs9;
const double clhs137 = DN(1,0)*clhs16;
const double clhs138 = DN(1,1)*clhs9;
const double clhs139 = clhs137 + clhs138;
const double clhs140 = C(0,0)*clhs136 + C(0,1)*clhs135 + C(0,2)*clhs139;
const double clhs141 = clhs140*clhs58;
const double clhs142 = C(0,2)*clhs136 + C(1,2)*clhs135 + C(2,2)*clhs139;
const double clhs143 = clhs142*clhs65;
const double clhs144 = clhs133*clhs68;
const double clhs145 = -clhs133*clhs74;
const double clhs146 = clhs135*clhs73 + clhs145*clhs17 + clhs145*clhs23;
const double clhs147 = clhs10*clhs145 + clhs136*clhs73 + clhs145*clhs4;
const double clhs148 = 2*clhs137*clhs78 + 2*clhs138*clhs78 + 2*clhs145*clhs26 + 2*clhs145*clhs27;
const double clhs149 = C(0,0)*clhs147 + C(0,1)*clhs146 + C(0,2)*clhs148;
const double clhs150 = C(0,2)*clhs147 + C(1,2)*clhs146 + C(2,2)*clhs148;
const double clhs151 = DN(1,0)*clhs90 + DN(1,1)*clhs91;
const double clhs152 = clhs149*clhs83 + clhs150*clhs85 + clhs151;
const double clhs153 = C(0,1)*clhs147 + C(1,1)*clhs146 + C(1,2)*clhs148;
const double clhs154 = DN(1,0)*clhs91 + DN(1,1)*clhs95;
const double clhs155 = clhs150*clhs83 + clhs153*clhs85 + clhs154;
const double clhs156 = -DN(0,1)*DN(1,0)*u(0,0) - DN(1,0)*DN(2,1)*u(2,0) + DN(1,1)*clhs5 + DN(1,1)*clhs7 + DN(1,1);
const double clhs157 = clhs156*clhs44;
const double clhs158 = DN(1,0)*clhs3;
const double clhs159 = DN(1,1)*clhs22;
const double clhs160 = DN(1,1)*clhs3;
const double clhs161 = DN(1,0)*clhs22;
const double clhs162 = clhs160 + clhs161;
const double clhs163 = C(0,0)*clhs158 + C(0,1)*clhs159 + C(0,2)*clhs162;
const double clhs164 = clhs163*clhs58;
const double clhs165 = C(0,2)*clhs158 + C(1,2)*clhs159 + C(2,2)*clhs162;
const double clhs166 = clhs165*clhs65;
const double clhs167 = clhs156*clhs68;
const double clhs168 = -clhs156*clhs74;
const double clhs169 = clhs10*clhs168 + clhs158*clhs73 + clhs168*clhs4;
const double clhs170 = clhs159*clhs73 + clhs168*clhs17 + clhs168*clhs23;
const double clhs171 = 2*clhs160*clhs78 + 2*clhs161*clhs78 + 2*clhs168*clhs26 + 2*clhs168*clhs27;
const double clhs172 = C(0,2)*clhs169 + C(1,2)*clhs170 + C(2,2)*clhs171;
const double clhs173 = C(0,0)*clhs169 + C(0,1)*clhs170 + C(0,2)*clhs171;
const double clhs174 = clhs16*clhs172 + clhs173*clhs9;
const double clhs175 = C(0,1)*clhs169 + C(1,1)*clhs170 + C(1,2)*clhs171;
const double clhs176 = clhs16*clhs175 + clhs172*clhs9;
const double clhs177 = -DN(0,0)*DN(2,1)*u(0,1) - DN(1,0)*DN(2,1)*u(1,1) + DN(2,0)*clhs18 + DN(2,0)*clhs19 + DN(2,0);
const double clhs178 = clhs177*clhs44;
const double clhs179 = DN(2,1)*clhs16;
const double clhs180 = DN(2,0)*clhs9;
const double clhs181 = DN(2,0)*clhs16;
const double clhs182 = DN(2,1)*clhs9;
const double clhs183 = clhs181 + clhs182;
const double clhs184 = C(0,0)*clhs180 + C(0,1)*clhs179 + C(0,2)*clhs183;
const double clhs185 = clhs184*clhs58;
const double clhs186 = C(0,2)*clhs180 + C(1,2)*clhs179 + C(2,2)*clhs183;
const double clhs187 = clhs186*clhs65;
const double clhs188 = clhs177*clhs68;
const double clhs189 = -clhs177*clhs74;
const double clhs190 = clhs17*clhs189 + clhs179*clhs73 + clhs189*clhs23;
const double clhs191 = clhs10*clhs189 + clhs180*clhs73 + clhs189*clhs4;
const double clhs192 = 2*clhs181*clhs78 + 2*clhs182*clhs78 + 2*clhs189*clhs26 + 2*clhs189*clhs27;
const double clhs193 = C(0,0)*clhs191 + C(0,1)*clhs190 + C(0,2)*clhs192;
const double clhs194 = C(0,2)*clhs191 + C(1,2)*clhs190 + C(2,2)*clhs192;
const double clhs195 = DN(2,0)*clhs90 + DN(2,1)*clhs91;
const double clhs196 = clhs193*clhs83 + clhs194*clhs85 + clhs195;
const double clhs197 = C(0,1)*clhs191 + C(1,1)*clhs190 + C(1,2)*clhs192;
const double clhs198 = DN(2,0)*clhs91 + DN(2,1)*clhs95;
const double clhs199 = clhs194*clhs83 + clhs197*clhs85 + clhs198;
const double clhs200 = -DN(0,1)*DN(2,0)*u(0,0) - DN(1,1)*DN(2,0)*u(1,0) + DN(2,1)*clhs5 + DN(2,1)*clhs6 + DN(2,1);
const double clhs201 = clhs200*clhs44;
const double clhs202 = DN(2,0)*clhs3;
const double clhs203 = DN(2,1)*clhs22;
const double clhs204 = DN(2,1)*clhs3;
const double clhs205 = DN(2,0)*clhs22;
const double clhs206 = clhs204 + clhs205;
const double clhs207 = C(0,0)*clhs202 + C(0,1)*clhs203 + C(0,2)*clhs206;
const double clhs208 = clhs207*clhs58;
const double clhs209 = C(0,2)*clhs202 + C(1,2)*clhs203 + C(2,2)*clhs206;
const double clhs210 = clhs209*clhs65;
const double clhs211 = clhs200*clhs68;
const double clhs212 = -clhs200*clhs74;
const double clhs213 = clhs10*clhs212 + clhs202*clhs73 + clhs212*clhs4;
const double clhs214 = clhs17*clhs212 + clhs203*clhs73 + clhs212*clhs23;
const double clhs215 = 2*clhs204*clhs78 + 2*clhs205*clhs78 + 2*clhs212*clhs26 + 2*clhs212*clhs27;
const double clhs216 = C(0,2)*clhs213 + C(1,2)*clhs214 + C(2,2)*clhs215;
const double clhs217 = C(0,0)*clhs213 + C(0,1)*clhs214 + C(0,2)*clhs215;
const double clhs218 = clhs16*clhs216 + clhs217*clhs9;
const double clhs219 = C(0,1)*clhs213 + C(1,1)*clhs214 + C(1,2)*clhs215;
const double clhs220 = clhs16*clhs219 + clhs216*clhs9;
const double clhs221 = DN(0,0)*clhs48;
const double clhs222 = clhs221*clhs50;
const double clhs223 = C(1,2)*clhs29 + clhs128 + clhs129;
const double clhs224 = DN(0,1)*clhs223;
const double clhs225 = clhs224*clhs50;
const double clhs226 = clhs58*clhs64;
const double clhs227 = C(0,1)*clhs53 + C(1,1)*clhs52 + C(1,2)*clhs56;
const double clhs228 = clhs227*clhs65;
const double clhs229 = clhs62*clhs69;
const double clhs230 = clhs22*clhs84 + clhs3*clhs80;
const double clhs231 = clhs22*clhs94 + clhs3*clhs84;
const double clhs232 = clhs107*clhs58;
const double clhs233 = C(0,1)*clhs100 + C(1,1)*clhs101 + C(1,2)*clhs104;
const double clhs234 = clhs233*clhs65;
const double clhs235 = clhs109*clhs62;
const double clhs236 = clhs3*clhs82;
const double clhs237 = clhs22*clhs82;
const double clhs238 = clhs114*clhs237 + clhs115*clhs236 + clhs92;
const double clhs239 = clhs114*clhs236 + clhs118*clhs237 + clhs96;
const double clhs240 = clhs125*clhs22 + clhs126*clhs3;
const double clhs241 = clhs125*clhs3 + clhs130*clhs22;
const double clhs242 = clhs44*(DN(0,0)*clhs123*clhs240 + DN(0,1)*clhs123*clhs241 - clhs122*clhs221 - clhs122*clhs224 - clhs222 - clhs225);
const double clhs243 = clhs142*clhs58;
const double clhs244 = C(0,1)*clhs136 + C(1,1)*clhs135 + C(1,2)*clhs139;
const double clhs245 = clhs244*clhs65;
const double clhs246 = clhs144*clhs62;
const double clhs247 = clhs149*clhs3 + clhs150*clhs22;
const double clhs248 = clhs150*clhs3 + clhs153*clhs22;
const double clhs249 = clhs165*clhs58;
const double clhs250 = C(0,1)*clhs158 + C(1,1)*clhs159 + C(1,2)*clhs162;
const double clhs251 = clhs250*clhs65;
const double clhs252 = clhs167*clhs62;
const double clhs253 = clhs151 + clhs172*clhs237 + clhs173*clhs236;
const double clhs254 = clhs154 + clhs172*clhs236 + clhs175*clhs237;
const double clhs255 = clhs186*clhs58;
const double clhs256 = C(0,1)*clhs180 + C(1,1)*clhs179 + C(1,2)*clhs183;
const double clhs257 = clhs256*clhs65;
const double clhs258 = clhs188*clhs62;
const double clhs259 = clhs193*clhs3 + clhs194*clhs22;
const double clhs260 = clhs194*clhs3 + clhs197*clhs22;
const double clhs261 = clhs209*clhs58;
const double clhs262 = C(0,1)*clhs202 + C(1,1)*clhs203 + C(1,2)*clhs206;
const double clhs263 = clhs262*clhs65;
const double clhs264 = clhs211*clhs62;
const double clhs265 = clhs195 + clhs216*clhs237 + clhs217*clhs236;
const double clhs266 = clhs198 + clhs216*clhs236 + clhs219*clhs237;
const double clhs267 = 0.1111111111111111/clhs43;
const double clhs268 = clhs267*clhs41;
const double clhs269 = DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2];
const double clhs270 = 1.0/clhs81;
const double clhs271 = pow(clhs270*clhs43, 0.33333333333333331);
const double clhs272 = clhs271*tau_u;
const double clhs273 = clhs269*clhs272;
const double clhs274 = DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2];
const double clhs275 = clhs272*clhs274;
const double clhs276 = clhs267*clhs98;
const double clhs277 = clhs31 + clhs49;
const double clhs278 = 0.33333333333333331*clhs272;
const double clhs279 = clhs277*clhs278;
const double clhs280 = clhs221 + clhs224;
const double clhs281 = clhs278*clhs280;
const double clhs282 = tau_u*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs283 = b_gauss[1]*clhs282;
const double clhs284 = clhs133*clhs267;
const double clhs285 = b_gauss[0]*clhs282;
const double clhs286 = clhs156*clhs267;
const double clhs287 = N[0]*N[1];
const double clhs288 = tau_u*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs289 = b_gauss[1]*clhs288;
const double clhs290 = clhs177*clhs267;
const double clhs291 = b_gauss[0]*clhs288;
const double clhs292 = clhs200*clhs267;
const double clhs293 = N[0]*N[2];
const double clhs294 = DN(1,0)*clhs30;
const double clhs295 = clhs294*clhs50;
const double clhs296 = DN(1,1)*clhs48;
const double clhs297 = clhs296*clhs50;
const double clhs298 = 0.66666666666666663*DN(1,0);
const double clhs299 = clhs298*clhs57;
const double clhs300 = 0.66666666666666663*DN(1,1);
const double clhs301 = clhs300*clhs64;
const double clhs302 = clhs105*clhs298;
const double clhs303 = clhs107*clhs300;
const double clhs304 = DN(1,0)*clhs82;
const double clhs305 = DN(1,1)*clhs82;
const double clhs306 = clhs44*(DN(1,0)*clhs123*clhs127 + DN(1,1)*clhs123*clhs131 - clhs122*clhs294 - clhs122*clhs296 - clhs295 - clhs297);
const double clhs307 = clhs140*clhs298;
const double clhs308 = clhs142*clhs300;
const double clhs309 = clhs163*clhs298;
const double clhs310 = clhs165*clhs300;
const double clhs311 = clhs184*clhs298;
const double clhs312 = clhs186*clhs300;
const double clhs313 = clhs207*clhs298;
const double clhs314 = clhs209*clhs300;
const double clhs315 = DN(1,0)*clhs48;
const double clhs316 = clhs315*clhs50;
const double clhs317 = DN(1,1)*clhs223;
const double clhs318 = clhs317*clhs50;
const double clhs319 = clhs298*clhs64;
const double clhs320 = clhs227*clhs300;
const double clhs321 = clhs107*clhs298;
const double clhs322 = clhs233*clhs300;
const double clhs323 = clhs44*(DN(1,0)*clhs123*clhs240 + DN(1,1)*clhs123*clhs241 - clhs122*clhs315 - clhs122*clhs317 - clhs316 - clhs318);
const double clhs324 = clhs142*clhs298;
const double clhs325 = clhs244*clhs300;
const double clhs326 = clhs165*clhs298;
const double clhs327 = clhs250*clhs300;
const double clhs328 = clhs186*clhs298;
const double clhs329 = clhs256*clhs300;
const double clhs330 = clhs209*clhs298;
const double clhs331 = clhs262*clhs300;
const double clhs332 = clhs294 + clhs296;
const double clhs333 = clhs278*clhs332;
const double clhs334 = clhs315 + clhs317;
const double clhs335 = clhs278*clhs334;
const double clhs336 = tau_u*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs337 = b_gauss[1]*clhs336;
const double clhs338 = b_gauss[0]*clhs336;
const double clhs339 = N[1]*N[2];
const double clhs340 = DN(2,0)*clhs30;
const double clhs341 = clhs340*clhs50;
const double clhs342 = DN(2,1)*clhs48;
const double clhs343 = clhs342*clhs50;
const double clhs344 = 0.66666666666666663*DN(2,0);
const double clhs345 = clhs344*clhs57;
const double clhs346 = 0.66666666666666663*DN(2,1);
const double clhs347 = clhs346*clhs64;
const double clhs348 = clhs105*clhs344;
const double clhs349 = clhs107*clhs346;
const double clhs350 = DN(2,0)*clhs82;
const double clhs351 = DN(2,1)*clhs82;
const double clhs352 = clhs44*(DN(2,0)*clhs123*clhs127 + DN(2,1)*clhs123*clhs131 - clhs122*clhs340 - clhs122*clhs342 - clhs341 - clhs343);
const double clhs353 = clhs140*clhs344;
const double clhs354 = clhs142*clhs346;
const double clhs355 = clhs163*clhs344;
const double clhs356 = clhs165*clhs346;
const double clhs357 = clhs184*clhs344;
const double clhs358 = clhs186*clhs346;
const double clhs359 = clhs207*clhs344;
const double clhs360 = clhs209*clhs346;
const double clhs361 = DN(2,0)*clhs48;
const double clhs362 = clhs361*clhs50;
const double clhs363 = DN(2,1)*clhs223;
const double clhs364 = clhs363*clhs50;
const double clhs365 = clhs344*clhs64;
const double clhs366 = clhs227*clhs346;
const double clhs367 = clhs107*clhs344;
const double clhs368 = clhs233*clhs346;
const double clhs369 = clhs44*(DN(2,0)*clhs123*clhs240 + DN(2,1)*clhs123*clhs241 - clhs122*clhs361 - clhs122*clhs363 - clhs362 - clhs364);
const double clhs370 = clhs142*clhs344;
const double clhs371 = clhs244*clhs346;
const double clhs372 = clhs165*clhs344;
const double clhs373 = clhs250*clhs346;
const double clhs374 = clhs186*clhs344;
const double clhs375 = clhs256*clhs346;
const double clhs376 = clhs209*clhs344;
const double clhs377 = clhs262*clhs346;
const double clhs378 = clhs340 + clhs342;
const double clhs379 = clhs278*clhs378;
const double clhs380 = clhs361 + clhs363;
const double clhs381 = clhs278*clhs380;
lhs(0,0)=DN(0,0)*clhs93 + DN(0,1)*clhs97 + clhs40*clhs45 + clhs45*clhs51 + clhs59*clhs63 + clhs63*clhs66 - clhs69*clhs71 - clhs69*clhs72;
lhs(0,1)=clhs106*clhs63 + clhs108*clhs63 - clhs109*clhs71 - clhs109*clhs72 + clhs116*clhs117 + clhs119*clhs120 + clhs40*clhs99 + clhs51*clhs99;
lhs(0,2)=N[0]*clhs132;
lhs(0,3)=DN(0,0)*clhs152 + DN(0,1)*clhs155 + clhs134*clhs40 + clhs134*clhs51 + clhs141*clhs63 + clhs143*clhs63 - clhs144*clhs71 - clhs144*clhs72;
lhs(0,4)=clhs117*clhs174 + clhs120*clhs176 + clhs157*clhs40 + clhs157*clhs51 + clhs164*clhs63 + clhs166*clhs63 - clhs167*clhs71 - clhs167*clhs72;
lhs(0,5)=N[1]*clhs132;
lhs(0,6)=DN(0,0)*clhs196 + DN(0,1)*clhs199 + clhs178*clhs40 + clhs178*clhs51 + clhs185*clhs63 + clhs187*clhs63 - clhs188*clhs71 - clhs188*clhs72;
lhs(0,7)=clhs117*clhs218 + clhs120*clhs220 + clhs201*clhs40 + clhs201*clhs51 + clhs208*clhs63 + clhs210*clhs63 - clhs211*clhs71 - clhs211*clhs72;
lhs(0,8)=N[2]*clhs132;
lhs(1,0)=clhs117*clhs230 + clhs120*clhs231 - clhs221*clhs229 + clhs222*clhs45 - clhs224*clhs229 + clhs225*clhs45 + clhs226*clhs63 + clhs228*clhs63;
lhs(1,1)=DN(0,0)*clhs238 + DN(0,1)*clhs239 - clhs221*clhs235 + clhs222*clhs99 - clhs224*clhs235 + clhs225*clhs99 + clhs232*clhs63 + clhs234*clhs63;
lhs(1,2)=N[0]*clhs242;
lhs(1,3)=clhs117*clhs247 + clhs120*clhs248 + clhs134*clhs222 + clhs134*clhs225 - clhs221*clhs246 - clhs224*clhs246 + clhs243*clhs63 + clhs245*clhs63;
lhs(1,4)=DN(0,0)*clhs253 + DN(0,1)*clhs254 + clhs157*clhs222 + clhs157*clhs225 - clhs221*clhs252 - clhs224*clhs252 + clhs249*clhs63 + clhs251*clhs63;
lhs(1,5)=N[1]*clhs242;
lhs(1,6)=clhs117*clhs259 + clhs120*clhs260 + clhs178*clhs222 + clhs178*clhs225 - clhs221*clhs258 - clhs224*clhs258 + clhs255*clhs63 + clhs257*clhs63;
lhs(1,7)=DN(0,0)*clhs265 + DN(0,1)*clhs266 + clhs201*clhs222 + clhs201*clhs225 - clhs221*clhs264 - clhs224*clhs264 + clhs261*clhs63 + clhs263*clhs63;
lhs(1,8)=N[2]*clhs242;
lhs(2,0)=N[0]*clhs41 - clhs273*(clhs268*clhs31 + clhs268*clhs49 + clhs59 + clhs66) - clhs275*(clhs221*clhs268 + clhs224*clhs268 + clhs226 + clhs228);
lhs(2,1)=N[0]*clhs98 - clhs273*(clhs106 + clhs108 + clhs276*clhs31 + clhs276*clhs49) - clhs275*(clhs221*clhs276 + clhs224*clhs276 + clhs232 + clhs234);
lhs(2,2)=-DN(0,0)*clhs279 - DN(0,1)*clhs281 - pow(N[0], 2) + 0.1111111111111111*N[0]*clhs269*clhs270*clhs271*clhs277*tau_u + 0.1111111111111111*N[0]*clhs270*clhs271*clhs274*clhs280*tau_u;
lhs(2,3)=N[0]*clhs133 - clhs273*(clhs141 + clhs143 + clhs284*clhs31 + clhs284*clhs49) - clhs275*(clhs221*clhs284 + clhs224*clhs284 + clhs243 + clhs245) + clhs283;
lhs(2,4)=N[0]*clhs156 - clhs273*(clhs164 + clhs166 + clhs286*clhs31 + clhs286*clhs49) - clhs275*(clhs221*clhs286 + clhs224*clhs286 + clhs249 + clhs251) - clhs285;
lhs(2,5)=-DN(1,0)*clhs279 - DN(1,1)*clhs281 + 0.1111111111111111*N[1]*clhs269*clhs270*clhs271*clhs277*tau_u + 0.1111111111111111*N[1]*clhs270*clhs271*clhs274*clhs280*tau_u - clhs287;
lhs(2,6)=N[0]*clhs177 - clhs273*(clhs185 + clhs187 + clhs290*clhs31 + clhs290*clhs49) - clhs275*(clhs221*clhs290 + clhs224*clhs290 + clhs255 + clhs257) + clhs289;
lhs(2,7)=N[0]*clhs200 - clhs273*(clhs208 + clhs210 + clhs292*clhs31 + clhs292*clhs49) - clhs275*(clhs221*clhs292 + clhs224*clhs292 + clhs261 + clhs263) - clhs291;
lhs(2,8)=-DN(2,0)*clhs279 - DN(2,1)*clhs281 + 0.1111111111111111*N[2]*clhs269*clhs270*clhs271*clhs277*tau_u + 0.1111111111111111*N[2]*clhs270*clhs271*clhs274*clhs280*tau_u - clhs293;
lhs(3,0)=DN(1,0)*clhs93 + DN(1,1)*clhs97 - clhs229*clhs294 - clhs229*clhs296 + clhs295*clhs45 + clhs297*clhs45 + clhs299*clhs63 + clhs301*clhs63;
lhs(3,1)=clhs116*clhs304 + clhs119*clhs305 - clhs235*clhs294 - clhs235*clhs296 + clhs295*clhs99 + clhs297*clhs99 + clhs302*clhs63 + clhs303*clhs63;
lhs(3,2)=N[0]*clhs306;
lhs(3,3)=DN(1,0)*clhs152 + DN(1,1)*clhs155 + clhs134*clhs295 + clhs134*clhs297 - clhs246*clhs294 - clhs246*clhs296 + clhs307*clhs63 + clhs308*clhs63;
lhs(3,4)=clhs157*clhs295 + clhs157*clhs297 + clhs174*clhs304 + clhs176*clhs305 - clhs252*clhs294 - clhs252*clhs296 + clhs309*clhs63 + clhs310*clhs63;
lhs(3,5)=N[1]*clhs306;
lhs(3,6)=DN(1,0)*clhs196 + DN(1,1)*clhs199 + clhs178*clhs295 + clhs178*clhs297 - clhs258*clhs294 - clhs258*clhs296 + clhs311*clhs63 + clhs312*clhs63;
lhs(3,7)=clhs201*clhs295 + clhs201*clhs297 + clhs218*clhs304 + clhs220*clhs305 - clhs264*clhs294 - clhs264*clhs296 + clhs313*clhs63 + clhs314*clhs63;
lhs(3,8)=N[2]*clhs306;
lhs(4,0)=-clhs229*clhs315 - clhs229*clhs317 + clhs230*clhs304 + clhs231*clhs305 + clhs316*clhs45 + clhs318*clhs45 + clhs319*clhs63 + clhs320*clhs63;
lhs(4,1)=DN(1,0)*clhs238 + DN(1,1)*clhs239 - clhs235*clhs315 - clhs235*clhs317 + clhs316*clhs99 + clhs318*clhs99 + clhs321*clhs63 + clhs322*clhs63;
lhs(4,2)=N[0]*clhs323;
lhs(4,3)=clhs134*clhs316 + clhs134*clhs318 - clhs246*clhs315 - clhs246*clhs317 + clhs247*clhs304 + clhs248*clhs305 + clhs324*clhs63 + clhs325*clhs63;
lhs(4,4)=DN(1,0)*clhs253 + DN(1,1)*clhs254 + clhs157*clhs316 + clhs157*clhs318 - clhs252*clhs315 - clhs252*clhs317 + clhs326*clhs63 + clhs327*clhs63;
lhs(4,5)=N[1]*clhs323;
lhs(4,6)=clhs178*clhs316 + clhs178*clhs318 - clhs258*clhs315 - clhs258*clhs317 + clhs259*clhs304 + clhs260*clhs305 + clhs328*clhs63 + clhs329*clhs63;
lhs(4,7)=DN(1,0)*clhs265 + DN(1,1)*clhs266 + clhs201*clhs316 + clhs201*clhs318 - clhs264*clhs315 - clhs264*clhs317 + clhs330*clhs63 + clhs331*clhs63;
lhs(4,8)=N[2]*clhs323;
lhs(5,0)=N[1]*clhs41 - clhs273*(clhs268*clhs294 + clhs268*clhs296 + clhs299 + clhs301) - clhs275*(clhs268*clhs315 + clhs268*clhs317 + clhs319 + clhs320) - clhs283;
lhs(5,1)=N[1]*clhs98 - clhs273*(clhs276*clhs294 + clhs276*clhs296 + clhs302 + clhs303) - clhs275*(clhs276*clhs315 + clhs276*clhs317 + clhs321 + clhs322) + clhs285;
lhs(5,2)=-DN(0,0)*clhs333 - DN(0,1)*clhs335 + 0.1111111111111111*N[0]*clhs269*clhs270*clhs271*clhs332*tau_u + 0.1111111111111111*N[0]*clhs270*clhs271*clhs274*clhs334*tau_u - clhs287;
lhs(5,3)=N[1]*clhs133 - clhs273*(clhs284*clhs294 + clhs284*clhs296 + clhs307 + clhs308) - clhs275*(clhs284*clhs315 + clhs284*clhs317 + clhs324 + clhs325);
lhs(5,4)=N[1]*clhs156 - clhs273*(clhs286*clhs294 + clhs286*clhs296 + clhs309 + clhs310) - clhs275*(clhs286*clhs315 + clhs286*clhs317 + clhs326 + clhs327);
lhs(5,5)=-DN(1,0)*clhs333 - DN(1,1)*clhs335 - pow(N[1], 2) + 0.1111111111111111*N[1]*clhs269*clhs270*clhs271*clhs332*tau_u + 0.1111111111111111*N[1]*clhs270*clhs271*clhs274*clhs334*tau_u;
lhs(5,6)=N[1]*clhs177 - clhs273*(clhs290*clhs294 + clhs290*clhs296 + clhs311 + clhs312) - clhs275*(clhs290*clhs315 + clhs290*clhs317 + clhs328 + clhs329) + clhs337;
lhs(5,7)=N[1]*clhs200 - clhs273*(clhs292*clhs294 + clhs292*clhs296 + clhs313 + clhs314) - clhs275*(clhs292*clhs315 + clhs292*clhs317 + clhs330 + clhs331) - clhs338;
lhs(5,8)=-DN(2,0)*clhs333 - DN(2,1)*clhs335 + 0.1111111111111111*N[2]*clhs269*clhs270*clhs271*clhs332*tau_u + 0.1111111111111111*N[2]*clhs270*clhs271*clhs274*clhs334*tau_u - clhs339;
lhs(6,0)=DN(2,0)*clhs93 + DN(2,1)*clhs97 - clhs229*clhs340 - clhs229*clhs342 + clhs341*clhs45 + clhs343*clhs45 + clhs345*clhs63 + clhs347*clhs63;
lhs(6,1)=clhs116*clhs350 + clhs119*clhs351 - clhs235*clhs340 - clhs235*clhs342 + clhs341*clhs99 + clhs343*clhs99 + clhs348*clhs63 + clhs349*clhs63;
lhs(6,2)=N[0]*clhs352;
lhs(6,3)=DN(2,0)*clhs152 + DN(2,1)*clhs155 + clhs134*clhs341 + clhs134*clhs343 - clhs246*clhs340 - clhs246*clhs342 + clhs353*clhs63 + clhs354*clhs63;
lhs(6,4)=clhs157*clhs341 + clhs157*clhs343 + clhs174*clhs350 + clhs176*clhs351 - clhs252*clhs340 - clhs252*clhs342 + clhs355*clhs63 + clhs356*clhs63;
lhs(6,5)=N[1]*clhs352;
lhs(6,6)=DN(2,0)*clhs196 + DN(2,1)*clhs199 + clhs178*clhs341 + clhs178*clhs343 - clhs258*clhs340 - clhs258*clhs342 + clhs357*clhs63 + clhs358*clhs63;
lhs(6,7)=clhs201*clhs341 + clhs201*clhs343 + clhs218*clhs350 + clhs220*clhs351 - clhs264*clhs340 - clhs264*clhs342 + clhs359*clhs63 + clhs360*clhs63;
lhs(6,8)=N[2]*clhs352;
lhs(7,0)=-clhs229*clhs361 - clhs229*clhs363 + clhs230*clhs350 + clhs231*clhs351 + clhs362*clhs45 + clhs364*clhs45 + clhs365*clhs63 + clhs366*clhs63;
lhs(7,1)=DN(2,0)*clhs238 + DN(2,1)*clhs239 - clhs235*clhs361 - clhs235*clhs363 + clhs362*clhs99 + clhs364*clhs99 + clhs367*clhs63 + clhs368*clhs63;
lhs(7,2)=N[0]*clhs369;
lhs(7,3)=clhs134*clhs362 + clhs134*clhs364 - clhs246*clhs361 - clhs246*clhs363 + clhs247*clhs350 + clhs248*clhs351 + clhs370*clhs63 + clhs371*clhs63;
lhs(7,4)=DN(2,0)*clhs253 + DN(2,1)*clhs254 + clhs157*clhs362 + clhs157*clhs364 - clhs252*clhs361 - clhs252*clhs363 + clhs372*clhs63 + clhs373*clhs63;
lhs(7,5)=N[1]*clhs369;
lhs(7,6)=clhs178*clhs362 + clhs178*clhs364 - clhs258*clhs361 - clhs258*clhs363 + clhs259*clhs350 + clhs260*clhs351 + clhs374*clhs63 + clhs375*clhs63;
lhs(7,7)=DN(2,0)*clhs265 + DN(2,1)*clhs266 + clhs201*clhs362 + clhs201*clhs364 - clhs264*clhs361 - clhs264*clhs363 + clhs376*clhs63 + clhs377*clhs63;
lhs(7,8)=N[2]*clhs369;
lhs(8,0)=N[2]*clhs41 - clhs273*(clhs268*clhs340 + clhs268*clhs342 + clhs345 + clhs347) - clhs275*(clhs268*clhs361 + clhs268*clhs363 + clhs365 + clhs366) - clhs289;
lhs(8,1)=N[2]*clhs98 - clhs273*(clhs276*clhs340 + clhs276*clhs342 + clhs348 + clhs349) - clhs275*(clhs276*clhs361 + clhs276*clhs363 + clhs367 + clhs368) + clhs291;
lhs(8,2)=-DN(0,0)*clhs379 - DN(0,1)*clhs381 + 0.1111111111111111*N[0]*clhs269*clhs270*clhs271*clhs378*tau_u + 0.1111111111111111*N[0]*clhs270*clhs271*clhs274*clhs380*tau_u - clhs293;
lhs(8,3)=N[2]*clhs133 - clhs273*(clhs284*clhs340 + clhs284*clhs342 + clhs353 + clhs354) - clhs275*(clhs284*clhs361 + clhs284*clhs363 + clhs370 + clhs371) - clhs337;
lhs(8,4)=N[2]*clhs156 - clhs273*(clhs286*clhs340 + clhs286*clhs342 + clhs355 + clhs356) - clhs275*(clhs286*clhs361 + clhs286*clhs363 + clhs372 + clhs373) + clhs338;
lhs(8,5)=-DN(1,0)*clhs379 - DN(1,1)*clhs381 + 0.1111111111111111*N[1]*clhs269*clhs270*clhs271*clhs378*tau_u + 0.1111111111111111*N[1]*clhs270*clhs271*clhs274*clhs380*tau_u - clhs339;
lhs(8,6)=N[2]*clhs177 - clhs273*(clhs290*clhs340 + clhs290*clhs342 + clhs357 + clhs358) - clhs275*(clhs290*clhs361 + clhs290*clhs363 + clhs374 + clhs375);
lhs(8,7)=N[2]*clhs200 - clhs273*(clhs292*clhs340 + clhs292*clhs342 + clhs359 + clhs360) - clhs275*(clhs292*clhs361 + clhs292*clhs363 + clhs376 + clhs377);
lhs(8,8)=-DN(2,0)*clhs379 - DN(2,1)*clhs381 - pow(N[2], 2) + 0.1111111111111111*N[2]*clhs269*clhs270*clhs271*clhs378*tau_u + 0.1111111111111111*N[2]*clhs270*clhs271*clhs274*clhs380*tau_u;

        //TODO: Amend this once the assembly is done in the input arrays
        rLeftHandSideMatrix += w_gauss * lhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < 3; ++d) {
            kinematic_variables.Displacements(i_node * 3 + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rLeftHandSideMatrix.clear();

    // Calculate stabilization constant
    const double c_tau_u = 2.0;
    const double c_tau_th = 4.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau_u * std::pow(h,2) / 2.0;

    // Set the auxiliary references matching the automatic differentiation symbols
    array_1d<double,3> b_gauss;
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_PK2);

        // Calculate body force
        // Note that this already includes the density computed in the reference configuration
        b_gauss = StructuralMechanicsElementUtilities::GetBodyForce(*this, r_integration_points, i_gauss);

        // Calculate the stabilization constant
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
        double kappa = CalculateBulkModulus(constitutive_variables.ConstitutiveMatrix); // Equivalent bulk modulus
        const double tau_u = aux_tau / mu;
        const double tau_th = (c_tau_th * mu) / (mu + kappa);

        // Calculate and add the LHS Gauss point contributions
        //substitute_lhs_3D_4N
        //TODO: Amend this once the assembly is done in the input arrays
        rLeftHandSideMatrix += w_gauss * lhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();

    // Check RHS size
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < 2; ++d) {
            kinematic_variables.Displacements(i_node * 2 + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rRightHandSideVector.clear();

    // Calculate stabilization constant
    const double c_tau_u = 2.0;
    const double c_tau_th = 4.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau_u * std::pow(h,2) / 2.0;

    // Set the auxiliary references matching the automatic differentiation symbols
    array_1d<double,3> b_gauss;
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedVector<double,LocalSize> rhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_PK2);

        // Calculate body force
        // Note that this already includes the density computed in the reference configuration
        b_gauss = StructuralMechanicsElementUtilities::GetBodyForce(*this, r_integration_points, i_gauss);

        // Calculate the stabilization constant
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
        double kappa = CalculateBulkModulus(constitutive_variables.ConstitutiveMatrix); // Equivalent bulk modulus
        const double tau_u = aux_tau / mu;
        const double tau_th = (c_tau_th * mu) / (mu + kappa);

        // Calculate and add the RHS Gauss point contribution
        const double crhs0 = DN(0,1)*u(0,0);
const double crhs1 = DN(1,1)*u(1,0);
const double crhs2 = DN(2,1)*u(2,0);
const double crhs3 = crhs0 + crhs1 + crhs2;
const double crhs4 = DN(0,0)*u(0,0);
const double crhs5 = DN(1,0)*u(1,0);
const double crhs6 = DN(2,0)*u(2,0);
const double crhs7 = crhs4 + crhs5 + crhs6;
const double crhs8 = crhs7 + 1;
const double crhs9 = S[0]*crhs8 + S[2]*crhs3;
const double crhs10 = S[1]*crhs3 + S[2]*crhs8;
const double crhs11 = DN(0,0)*u(0,1);
const double crhs12 = DN(1,0)*u(1,1);
const double crhs13 = DN(2,0)*u(2,1);
const double crhs14 = crhs11 + crhs12 + crhs13;
const double crhs15 = pow(crhs14, 2) + pow(crhs8, 2);
const double crhs16 = DN(0,1)*u(0,1);
const double crhs17 = DN(1,1)*u(1,1);
const double crhs18 = DN(2,1)*u(2,1);
const double crhs19 = crhs16 + crhs17 + crhs18;
const double crhs20 = crhs19 + 1;
const double crhs21 = pow(crhs20, 2) + pow(crhs3, 2);
const double crhs22 = 2*crhs14*crhs20 + 2*crhs3*crhs8;
const double crhs23 = C(0,0)*crhs15 + C(0,1)*crhs21 + C(0,2)*crhs22;
const double crhs24 = DN(0,0)*crhs23;
const double crhs25 = N[0]*th[0];
const double crhs26 = N[1]*th[1];
const double crhs27 = N[2]*th[2];
const double crhs28 = -crhs0*crhs12 - crhs0*crhs13 - crhs1*crhs11 - crhs1*crhs13 - crhs11*crhs2 - crhs12*crhs2 + crhs16*crhs5 + crhs16*crhs6 + crhs17*crhs4 + crhs17*crhs6 + crhs18*crhs4 + crhs18*crhs5 + crhs19;
const double crhs29 = -crhs25 - crhs26 - crhs27 + crhs28 + crhs7;
const double crhs30 = crhs25 + crhs26 + crhs27;
const double crhs31 = crhs28 + crhs8;
const double crhs32 = 0.33333333333333331*crhs29*pow(crhs31, -0.66666666666666663)*tau_th*pow(crhs30 + 1, -0.33333333333333331);
const double crhs33 = C(0,2)*crhs15 + C(1,2)*crhs21 + C(2,2)*crhs22;
const double crhs34 = DN(0,1)*crhs33;
const double crhs35 = S[0]*crhs14 + S[2]*crhs20;
const double crhs36 = S[1]*crhs20 + S[2]*crhs14;
const double crhs37 = DN(0,0)*crhs33;
const double crhs38 = C(0,1)*crhs15 + C(1,1)*crhs21 + C(1,2)*crhs22;
const double crhs39 = DN(0,1)*crhs38;
const double crhs40 = b_gauss[0]*tau_u;
const double crhs41 = b_gauss[1]*tau_u;
const double crhs42 = 0.33333333333333331*tau_u*pow(crhs31/(crhs30 + 1.0), 0.33333333333333331);
const double crhs43 = crhs42*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs44 = crhs42*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double crhs45 = DN(1,0)*crhs23;
const double crhs46 = DN(1,1)*crhs33;
const double crhs47 = DN(1,0)*crhs33;
const double crhs48 = DN(1,1)*crhs38;
const double crhs49 = DN(2,0)*crhs23;
const double crhs50 = DN(2,1)*crhs33;
const double crhs51 = DN(2,0)*crhs33;
const double crhs52 = DN(2,1)*crhs38;
rhs[0]=-DN(0,0)*crhs9 - DN(0,1)*crhs10 + N[0]*b_gauss[0] - crhs24*crhs32 - crhs32*crhs34;
rhs[1]=-DN(0,0)*crhs35 - DN(0,1)*crhs36 + N[0]*b_gauss[1] - crhs32*crhs37 - crhs32*crhs39;
rhs[2]=-N[0]*crhs29 + crhs40*(DN(0,0)*crhs20 - DN(0,1)*crhs14) - crhs41*(DN(0,0)*crhs3 - DN(0,1)*crhs8) + crhs43*(crhs24 + crhs34) + crhs44*(crhs37 + crhs39);
rhs[3]=-DN(1,0)*crhs9 - DN(1,1)*crhs10 + N[1]*b_gauss[0] - crhs32*crhs45 - crhs32*crhs46;
rhs[4]=-DN(1,0)*crhs35 - DN(1,1)*crhs36 + N[1]*b_gauss[1] - crhs32*crhs47 - crhs32*crhs48;
rhs[5]=-N[1]*crhs29 + crhs40*(DN(1,0)*crhs20 - DN(1,1)*crhs14) - crhs41*(DN(1,0)*crhs3 - DN(1,1)*crhs8) + crhs43*(crhs45 + crhs46) + crhs44*(crhs47 + crhs48);
rhs[6]=-DN(2,0)*crhs9 - DN(2,1)*crhs10 + N[2]*b_gauss[0] - crhs32*crhs49 - crhs32*crhs50;
rhs[7]=-DN(2,0)*crhs35 - DN(2,1)*crhs36 + N[2]*b_gauss[1] - crhs32*crhs51 - crhs32*crhs52;
rhs[8]=-N[2]*crhs29 + crhs40*(DN(2,0)*crhs20 - DN(2,1)*crhs14) - crhs41*(DN(2,0)*crhs3 - DN(2,1)*crhs8) + crhs43*(crhs49 + crhs50) + crhs44*(crhs51 + crhs52);

        //TODO: Amend this once the assembly is done in the input arrays
        rRightHandSideVector += w_gauss * rhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& r_geometry = GetGeometry();

    // Check RHS size
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < 3; ++d) {
            kinematic_variables.Displacements(i_node * 3 + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
    }

    // Create the constitutive variables and values containers
    ConstitutiveVariables constitutive_variables;
    ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
    auto& r_cons_law_options = cons_law_values.GetOptions();
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

    // Calculate the RHS and LHS contributions
    rRightHandSideVector.clear();

    // Calculate stabilization constant
    const double c_tau_u = 2.0;
    const double c_tau_th = 4.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau_u * std::pow(h,2) / 2.0;

    // Set the auxiliary references matching the automatic differentiation symbols
    array_1d<double,3> b_gauss;
    const auto& N = kinematic_variables.N;
    const auto& DN = kinematic_variables.DN_DX;
    const auto& u = kinematic_variables.Displacements;
    const auto& th = kinematic_variables.JacobianDeterminant;
    const auto& S = constitutive_variables.StressVector;
    const auto& C = constitutive_variables.ConstitutiveMatrix;

    // Aux RHS and LHS
    //TODO: To be removed
    BoundedVector<double,LocalSize> rhs;
    BoundedMatrix<double,LocalSize,LocalSize> lhs;

    const SizeType n_gauss = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());
    for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        // Calculate kinematics
        CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());
        const double w_gauss = kinematic_variables.detJ0 * r_integration_points[i_gauss].Weight();

        // Calculate the constitutive response
        CalculateConstitutiveVariables(
            kinematic_variables,
            constitutive_variables,
            cons_law_values,
            i_gauss,
            r_geometry.IntegrationPoints(this->GetIntegrationMethod()),
            ConstitutiveLaw::StressMeasure_PK2);

        // Calculate body force
        // Note that this already includes the density computed in the reference configuration
        b_gauss = StructuralMechanicsElementUtilities::GetBodyForce(*this, r_integration_points, i_gauss);

        // Calculate the stabilization constant
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
        double kappa = CalculateBulkModulus(constitutive_variables.ConstitutiveMatrix); // Equivalent bulk modulus
        const double tau_u = aux_tau / mu;
        const double tau_th = (c_tau_th * mu) / (mu + kappa);

        // Calculate and add the RHS Gauss point contribution
        //substitute_rhs_3D_4N
        //TODO: Amend this once the assembly is done in the input arrays
        rRightHandSideVector += w_gauss * rhs;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::InitializeMaterial()
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

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
bool TotalLagrangianMixedDetJElement<TDim>::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const
{
    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetStrainVector(rThisKinematicVariables.EquivalentStrain); // equivalent total strain
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // assuming that det(F) is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); // assuming that F is computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.ConstitutiveMatrix); //assuming the determinant is computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure) const
{
    // Set the constitutive variables
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod) const
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(rIntegrationMethod);

    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

    // Calculate the inverse Jacobian
    GeometryUtils::JacobianOnInitialConfiguration(
        r_geometry,
        r_integration_points[PointNumber],
        rThisKinematicVariables.J0);
    MathUtils<double>::InvertMatrix(
        rThisKinematicVariables.J0,
        rThisKinematicVariables.InvJ0,
        rThisKinematicVariables.detJ0);
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0)
        << "Element ID: " << this->Id() << " is inverted. det(J0) = " << rThisKinematicVariables.detJ0 << std::endl;

    // Calculate the shape functions gradients
    GeometryUtils::ShapeFunctionsGradients(
        r_geometry.ShapeFunctionsLocalGradients(rIntegrationMethod)[PointNumber],
        rThisKinematicVariables.InvJ0,
        rThisKinematicVariables.DN_DX);

    // Calculate the equivalent total strain
    CalculateEquivalentStrain(rThisKinematicVariables);

    // Compute equivalent F
    CalculateEquivalentF(rThisKinematicVariables);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const
{
    // Define references to the auxiliary symbols
    const auto& N = rThisKinematicVariables.N;
    const auto& DN = rThisKinematicVariables.DN_DX;
    const auto& u = rThisKinematicVariables.Displacements;
    const auto& th = rThisKinematicVariables.JacobianDeterminant;
    auto& r_eq_green_strain = rThisKinematicVariables.EquivalentStrain;

    // Fill the equivalent Green strain values
    const double cr_eq_green_strain0 = DN(0,0)*u(0,1);
const double cr_eq_green_strain1 = DN(1,0)*u(1,1);
const double cr_eq_green_strain2 = DN(2,0)*u(2,1);
const double cr_eq_green_strain3 = cr_eq_green_strain0 + cr_eq_green_strain1 + cr_eq_green_strain2;
const double cr_eq_green_strain4 = DN(0,0)*u(0,0);
const double cr_eq_green_strain5 = DN(1,1)*u(1,1);
const double cr_eq_green_strain6 = DN(2,1)*u(2,1);
const double cr_eq_green_strain7 = DN(0,1)*u(0,1);
const double cr_eq_green_strain8 = DN(1,0)*u(1,0);
const double cr_eq_green_strain9 = DN(2,0)*u(2,0);
const double cr_eq_green_strain10 = DN(1,1)*u(1,0);
const double cr_eq_green_strain11 = DN(2,1)*u(2,0);
const double cr_eq_green_strain12 = DN(0,1)*u(0,0);
const double cr_eq_green_strain13 = cr_eq_green_strain4 + cr_eq_green_strain8 + cr_eq_green_strain9 + 1;
const double cr_eq_green_strain14 = cr_eq_green_strain5 + cr_eq_green_strain6 + cr_eq_green_strain7;
const double cr_eq_green_strain15 = pow(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0, 0.66666666666666663)*pow(-cr_eq_green_strain0*cr_eq_green_strain10 - cr_eq_green_strain0*cr_eq_green_strain11 - cr_eq_green_strain1*cr_eq_green_strain11 - cr_eq_green_strain1*cr_eq_green_strain12 - cr_eq_green_strain10*cr_eq_green_strain2 - cr_eq_green_strain12*cr_eq_green_strain2 + cr_eq_green_strain13 + cr_eq_green_strain14 + cr_eq_green_strain4*cr_eq_green_strain5 + cr_eq_green_strain4*cr_eq_green_strain6 + cr_eq_green_strain5*cr_eq_green_strain9 + cr_eq_green_strain6*cr_eq_green_strain8 + cr_eq_green_strain7*cr_eq_green_strain8 + cr_eq_green_strain7*cr_eq_green_strain9, -0.66666666666666663);
const double cr_eq_green_strain16 = cr_eq_green_strain10 + cr_eq_green_strain11 + cr_eq_green_strain12;
const double cr_eq_green_strain17 = cr_eq_green_strain14 + 1;
r_eq_green_strain[0]=0.5*pow(cr_eq_green_strain13, 2)*cr_eq_green_strain15 + 0.5*cr_eq_green_strain15*pow(cr_eq_green_strain3, 2) - 0.5;
r_eq_green_strain[1]=0.5*cr_eq_green_strain15*pow(cr_eq_green_strain16, 2) + 0.5*cr_eq_green_strain15*pow(cr_eq_green_strain17, 2) - 0.5;
r_eq_green_strain[2]=1.0*cr_eq_green_strain15*(cr_eq_green_strain13*cr_eq_green_strain16 + cr_eq_green_strain17*cr_eq_green_strain3);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateEquivalentStrain(KinematicVariables& rThisKinematicVariables) const
{
    // Define references to the auxiliary symbols
    const auto& N = rThisKinematicVariables.N;
    const auto& DN = rThisKinematicVariables.DN_DX;
    const auto& u = rThisKinematicVariables.Displacements;
    const auto& th = rThisKinematicVariables.JacobianDeterminant;
    auto& r_eq_green_strain = rThisKinematicVariables.EquivalentStrain;

    // Fill the equivalent Green strain values
    //substitute_green_strain_3D_4N
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<2>::CalculateEquivalentF(KinematicVariables& rThisKinematicVariables) const
{
    // Define references to the auxiliary symbols
    const auto& N = rThisKinematicVariables.N;
    const auto& DN = rThisKinematicVariables.DN_DX;
    const auto& u = rThisKinematicVariables.Displacements;
    const auto& th = rThisKinematicVariables.JacobianDeterminant;
    auto& r_eq_def_gradient = rThisKinematicVariables.F;
    double& r_det_eq_def_gradient = rThisKinematicVariables.detF;

    // Fill the equivalent deformation gradient values
    const double cr_eq_def_gradient0 = DN(0,0)*u(0,0);
const double cr_eq_def_gradient1 = DN(1,0)*u(1,0);
const double cr_eq_def_gradient2 = DN(2,0)*u(2,0);
const double cr_eq_def_gradient3 = cr_eq_def_gradient0 + cr_eq_def_gradient1 + cr_eq_def_gradient2 + 1;
const double cr_eq_def_gradient4 = DN(1,1)*u(1,1);
const double cr_eq_def_gradient5 = DN(2,1)*u(2,1);
const double cr_eq_def_gradient6 = DN(0,1)*u(0,1);
const double cr_eq_def_gradient7 = DN(1,1)*u(1,0);
const double cr_eq_def_gradient8 = DN(0,0)*u(0,1);
const double cr_eq_def_gradient9 = DN(2,1)*u(2,0);
const double cr_eq_def_gradient10 = DN(0,1)*u(0,0);
const double cr_eq_def_gradient11 = DN(1,0)*u(1,1);
const double cr_eq_def_gradient12 = DN(2,0)*u(2,1);
const double cr_eq_def_gradient13 = cr_eq_def_gradient4 + cr_eq_def_gradient5 + cr_eq_def_gradient6;
const double cr_eq_def_gradient14 = pow(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0, 0.33333333333333331)*pow(cr_eq_def_gradient0*cr_eq_def_gradient4 + cr_eq_def_gradient0*cr_eq_def_gradient5 + cr_eq_def_gradient1*cr_eq_def_gradient5 + cr_eq_def_gradient1*cr_eq_def_gradient6 - cr_eq_def_gradient10*cr_eq_def_gradient11 - cr_eq_def_gradient10*cr_eq_def_gradient12 - cr_eq_def_gradient11*cr_eq_def_gradient9 - cr_eq_def_gradient12*cr_eq_def_gradient7 + cr_eq_def_gradient13 + cr_eq_def_gradient2*cr_eq_def_gradient4 + cr_eq_def_gradient2*cr_eq_def_gradient6 + cr_eq_def_gradient3 - cr_eq_def_gradient7*cr_eq_def_gradient8 - cr_eq_def_gradient8*cr_eq_def_gradient9, -0.33333333333333331);
r_eq_def_gradient(0,0)=cr_eq_def_gradient14*cr_eq_def_gradient3;
r_eq_def_gradient(0,1)=cr_eq_def_gradient14*(cr_eq_def_gradient10 + cr_eq_def_gradient7 + cr_eq_def_gradient9);
r_eq_def_gradient(1,0)=cr_eq_def_gradient14*(cr_eq_def_gradient11 + cr_eq_def_gradient12 + cr_eq_def_gradient8);
r_eq_def_gradient(1,1)=cr_eq_def_gradient14*(cr_eq_def_gradient13 + 1);

    const double cr_det_eq_def_gradient0 = DN(0,0)*u(0,0);
const double cr_det_eq_def_gradient1 = DN(0,1)*u(0,1);
const double cr_det_eq_def_gradient2 = DN(1,0)*u(1,0);
const double cr_det_eq_def_gradient3 = DN(1,1)*u(1,1);
const double cr_det_eq_def_gradient4 = DN(2,0)*u(2,0);
const double cr_det_eq_def_gradient5 = DN(2,1)*u(2,1);
const double cr_det_eq_def_gradient6 = DN(0,0)*u(0,1);
const double cr_det_eq_def_gradient7 = DN(1,1)*u(1,0);
const double cr_det_eq_def_gradient8 = DN(2,1)*u(2,0);
const double cr_det_eq_def_gradient9 = DN(0,1)*u(0,0);
const double cr_det_eq_def_gradient10 = DN(1,0)*u(1,1);
const double cr_det_eq_def_gradient11 = DN(2,0)*u(2,1);
r_det_eq_def_gradient=pow(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0, 0.66666666666666663)*pow(cr_det_eq_def_gradient0*cr_det_eq_def_gradient3 + cr_det_eq_def_gradient0*cr_det_eq_def_gradient5 + cr_det_eq_def_gradient0 + cr_det_eq_def_gradient1*cr_det_eq_def_gradient2 + cr_det_eq_def_gradient1*cr_det_eq_def_gradient4 + cr_det_eq_def_gradient1 - cr_det_eq_def_gradient10*cr_det_eq_def_gradient8 - cr_det_eq_def_gradient10*cr_det_eq_def_gradient9 - cr_det_eq_def_gradient11*cr_det_eq_def_gradient7 - cr_det_eq_def_gradient11*cr_det_eq_def_gradient9 + cr_det_eq_def_gradient2*cr_det_eq_def_gradient5 + cr_det_eq_def_gradient2 + cr_det_eq_def_gradient3*cr_det_eq_def_gradient4 + cr_det_eq_def_gradient3 + cr_det_eq_def_gradient4 + cr_det_eq_def_gradient5 - cr_det_eq_def_gradient6*cr_det_eq_def_gradient7 - cr_det_eq_def_gradient6*cr_det_eq_def_gradient8 + 1, 0.33333333333333337);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
void TotalLagrangianMixedDetJElement<3>::CalculateEquivalentF(KinematicVariables& rThisKinematicVariables) const
{
    // Define references to the auxiliary symbols
    const auto& N = rThisKinematicVariables.N;
    const auto& DN = rThisKinematicVariables.DN_DX;
    const auto& u = rThisKinematicVariables.Displacements;
    const auto& th = rThisKinematicVariables.JacobianDeterminant;
    auto& r_eq_def_gradient = rThisKinematicVariables.F;
    double& r_det_eq_def_gradient = rThisKinematicVariables.detF;

    // Fill the equivalent deformation gradient values
    //substitute_def_gradient_3D_4N
    //substitute_det_def_gradient_3D_4N
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
int  TotalLagrangianMixedDetJElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    int check = TotalLagrangianMixedDetJElement::BaseType::Check(rCurrentProcessInfo);

    // Base check
    check = StructuralMechanicsElementUtilities::SolidElementCheck(*this, rCurrentProcessInfo, mConstitutiveLawVector);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    const auto& r_geometry = this->GetGeometry();
    for ( IndexType i = 0; i < r_geometry.size(); i++ ) {
        const NodeType& r_node = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VOLUMETRIC_STRAIN,r_node)
        KRATOS_CHECK_DOF_IN_NODE(VOLUMETRIC_STRAIN, r_node)
    }

    return check;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::CalculateOnIntegrationPoints(
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

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = GetGeometry();
    const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

    const SizeType n_gauss = r_integration_points.size();
    if (rOutput.size() != n_gauss) {
        rOutput.resize(n_gauss);
    }

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else if (rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR) {
        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables;
        for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < TDim; ++d) {
                kinematic_variables.Displacements(i_node, d) = r_disp[d];
            }
            kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }

        // Create the constitutive variables and values containers
        ConstitutiveVariables constitutive_variables;
        ConstitutiveLaw::Parameters cons_law_values(r_geometry, GetProperties(), rCurrentProcessInfo);
        auto& r_cons_law_options = cons_law_values.GetOptions();
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_cons_law_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);


        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Call the constitutive law to update material variables
            if( rVariable == CAUCHY_STRESS_VECTOR) {
                // Compute material reponse
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
            } else {
                // Compute material reponse
                CalculateConstitutiveVariables(kinematic_variables, constitutive_variables, cons_law_values, i_gauss, r_integration_points, ConstitutiveLaw::StressMeasure_PK2);
            }

            // Check sizes and save the output stress
            if (rOutput[i_gauss].size() != StrainSize) {
                rOutput[i_gauss].resize(StrainSize, false);
            }
            rOutput[i_gauss] = constitutive_variables.StressVector;
        }
    } else if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables;
        for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < TDim; ++d) {
                kinematic_variables.Displacements(i_node, d) = r_disp[d];
            }
            kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(VOLUMETRIC_STRAIN);
        }

        for (IndexType i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Calculate kinematics
            CalculateKinematicVariables(kinematic_variables, i_gauss, GetIntegrationMethod());

            // Calculate Green-Lagrange strain
            TotalLagrangianMixedDetJElement<TDim>::CalculateEquivalentStrain(kinematic_variables);

            // Check sizes and save the output stress
            if (rOutput[i_gauss].size() != StrainSize) {
                rOutput[i_gauss].resize(StrainSize, false);
            }
            rOutput[i_gauss] = kinematic_variables.EquivalentStrain;
        }
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
const Parameters TotalLagrangianMixedDetJElement<TDim>::GetSpecifications() const
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
        "compatible_geometries"      : ["Triangle2D3", "Tetrahedra3D4"],
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
    if (TDim == 2) {
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

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TotalLagrangianMixedDetJElement::BaseType);
    int IntMethod = int(this->GetIntegrationMethod());
    rSerializer.save("IntegrationMethod",IntMethod);
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void TotalLagrangianMixedDetJElement<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TotalLagrangianMixedDetJElement::BaseType);
    int IntMethod;
    rSerializer.load("IntegrationMethod",IntMethod);
    mThisIntegrationMethod = IntegrationMethod(IntMethod);
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

// Explicit template instantiations
template class TotalLagrangianMixedDetJElement<2>;
template class TotalLagrangianMixedDetJElement<3>;

} // Namespace Kratos
