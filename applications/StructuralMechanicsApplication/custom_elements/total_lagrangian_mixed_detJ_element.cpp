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
const double clhs29 = C(0,2)*clhs28;
const double clhs30 = clhs12 + clhs25 + 2*clhs29;
const double clhs31 = DN(0,0)*clhs30;
const double clhs32 = 0.5*tau_th;
const double clhs33 = clhs31*clhs32;
const double clhs34 = -clhs0*clhs14 - clhs0*clhs15 - clhs1*clhs13 - clhs1*clhs15 - clhs13*clhs2 - clhs14*clhs2 + clhs18*clhs6 + clhs18*clhs7 + clhs19*clhs5 + clhs19*clhs7 + clhs20*clhs5 + clhs20*clhs6 + clhs21;
const double clhs35 = clhs34 + clhs9;
const double clhs36 = 1.0/clhs35;
const double clhs37 = DN(0,0)*clhs19 + DN(0,0)*clhs20 + DN(0,0) - DN(0,1)*DN(1,0)*u(1,1) - DN(0,1)*DN(2,0)*u(2,1);
const double clhs38 = clhs36*clhs37;
const double clhs39 = C(0,2)*clhs11;
const double clhs40 = C(1,2)*clhs24;
const double clhs41 = C(2,2)*clhs28;
const double clhs42 = clhs39 + clhs40 + 2*clhs41;
const double clhs43 = DN(0,1)*clhs42;
const double clhs44 = clhs32*clhs43;
const double clhs45 = DN(0,1)*clhs16;
const double clhs46 = DN(0,0)*clhs9;
const double clhs47 = DN(0,0)*clhs16 + DN(0,1)*clhs9;
const double clhs48 = C(0,0)*clhs46 + C(0,1)*clhs45 + C(0,2)*clhs47;
const double clhs49 = DN(0,0)*clhs48;
const double clhs50 = N[0]*th[0];
const double clhs51 = N[1]*th[1];
const double clhs52 = N[2]*th[2];
const double clhs53 = tau_th*(clhs34 - clhs50 - clhs51 - clhs52 + clhs8);
const double clhs54 = clhs36*clhs53;
const double clhs55 = C(0,2)*clhs46 + C(1,2)*clhs45 + C(2,2)*clhs47;
const double clhs56 = DN(0,1)*clhs55;
const double clhs57 = clhs31*clhs53;
const double clhs58 = 0.5/pow(clhs35, 2);
const double clhs59 = clhs37*clhs58;
const double clhs60 = clhs43*clhs53;
const double clhs61 = clhs50 + clhs51 + clhs52 + 1.0;
const double clhs62 = -clhs36*clhs37;
const double clhs63 = 0.5*clhs62;
const double clhs64 = clhs61*(clhs17*clhs63 + clhs23*clhs63 + 1.0*clhs45);
const double clhs65 = clhs61*(clhs10*clhs63 + clhs4*clhs63 + 1.0*clhs46);
const double clhs66 = clhs50 + clhs51 + clhs52 + 1.0;
const double clhs67 = clhs66*(clhs26*clhs62 + clhs27*clhs62 + clhs47);
const double clhs68 = C(0,2)*clhs65 + C(1,2)*clhs64 + C(2,2)*clhs67;
const double clhs69 = clhs16*clhs36;
const double clhs70 = clhs36*clhs9;
const double clhs71 = DN(0,0)*S[0] + DN(0,1)*S[2];
const double clhs72 = clhs68*clhs69 + clhs70*(C(0,0)*clhs65 + C(0,1)*clhs64 + C(0,2)*clhs67) + clhs71;
const double clhs73 = DN(0,0)*S[2] + DN(0,1)*S[1];
const double clhs74 = clhs68*clhs70 + clhs69*(C(0,1)*clhs65 + C(1,1)*clhs64 + C(1,2)*clhs67) + clhs73;
const double clhs75 = -DN(0,0)*clhs14 - DN(0,0)*clhs15 + DN(0,1)*clhs6 + DN(0,1)*clhs7 + DN(0,1);
const double clhs76 = DN(0,0)*clhs3;
const double clhs77 = DN(0,1)*clhs22;
const double clhs78 = DN(0,0)*clhs22 + DN(0,1)*clhs3;
const double clhs79 = C(0,0)*clhs76 + C(0,1)*clhs77 + C(0,2)*clhs78;
const double clhs80 = DN(0,0)*clhs79;
const double clhs81 = C(0,2)*clhs76 + C(1,2)*clhs77 + C(2,2)*clhs78;
const double clhs82 = DN(0,1)*clhs81;
const double clhs83 = clhs36*clhs75;
const double clhs84 = 0.5*clhs83;
const double clhs85 = clhs61*(1.0*DN(0,0)*clhs3 - clhs10*clhs84 - clhs4*clhs84);
const double clhs86 = clhs61*(1.0*DN(0,1)*clhs22 - clhs17*clhs84 - clhs23*clhs84);
const double clhs87 = clhs66*(-clhs26*clhs83 - clhs27*clhs83 + clhs78);
const double clhs88 = C(0,2)*clhs85 + C(1,2)*clhs86 + C(2,2)*clhs87;
const double clhs89 = clhs16*clhs88 + clhs9*(C(0,0)*clhs85 + C(0,1)*clhs86 + C(0,2)*clhs87);
const double clhs90 = clhs16*(C(0,1)*clhs85 + C(1,1)*clhs86 + C(1,2)*clhs87) + clhs88*clhs9;
const double clhs91 = 0.5*clhs39 + 0.5*clhs40 + clhs41;
const double clhs92 = 0.5*clhs12 + 0.5*clhs25 + clhs29;
const double clhs93 = clhs16*clhs91 + clhs9*clhs92;
const double clhs94 = C(0,1)*clhs11;
const double clhs95 = C(1,1)*clhs24;
const double clhs96 = C(1,2)*clhs28;
const double clhs97 = 0.5*clhs94 + 0.5*clhs95 + clhs96;
const double clhs98 = clhs16*clhs97 + clhs9*clhs91;
const double clhs99 = clhs36*(DN(0,0)*clhs93 + DN(0,1)*clhs98 - clhs33 - clhs44);
const double clhs100 = -DN(0,0)*DN(1,1)*u(0,1) + DN(1,0)*clhs18 + DN(1,0)*clhs20 + DN(1,0) - DN(1,1)*DN(2,0)*u(2,1);
const double clhs101 = clhs100*clhs36;
const double clhs102 = DN(1,1)*clhs16;
const double clhs103 = DN(1,0)*clhs9;
const double clhs104 = DN(1,0)*clhs16 + DN(1,1)*clhs9;
const double clhs105 = C(0,0)*clhs103 + C(0,1)*clhs102 + C(0,2)*clhs104;
const double clhs106 = DN(0,0)*clhs105;
const double clhs107 = C(0,2)*clhs103 + C(1,2)*clhs102 + C(2,2)*clhs104;
const double clhs108 = DN(0,1)*clhs107;
const double clhs109 = clhs100*clhs58;
const double clhs110 = -clhs100*clhs36;
const double clhs111 = 0.5*clhs110;
const double clhs112 = clhs61*(1.0*clhs102 + clhs111*clhs17 + clhs111*clhs23);
const double clhs113 = clhs61*(clhs10*clhs111 + 1.0*clhs103 + clhs111*clhs4);
const double clhs114 = clhs66*(clhs104 + clhs110*clhs26 + clhs110*clhs27);
const double clhs115 = C(0,2)*clhs113 + C(1,2)*clhs112 + C(2,2)*clhs114;
const double clhs116 = DN(1,0)*S[0] + DN(1,1)*S[2];
const double clhs117 = clhs115*clhs69 + clhs116 + clhs70*(C(0,0)*clhs113 + C(0,1)*clhs112 + C(0,2)*clhs114);
const double clhs118 = DN(1,0)*S[2] + DN(1,1)*S[1];
const double clhs119 = clhs115*clhs70 + clhs118 + clhs69*(C(0,1)*clhs113 + C(1,1)*clhs112 + C(1,2)*clhs114);
const double clhs120 = -DN(1,0)*clhs13 - DN(1,0)*clhs15 + DN(1,1)*clhs5 + DN(1,1)*clhs7 + DN(1,1);
const double clhs121 = DN(1,0)*clhs3;
const double clhs122 = DN(1,1)*clhs22;
const double clhs123 = DN(1,0)*clhs22 + DN(1,1)*clhs3;
const double clhs124 = C(0,0)*clhs121 + C(0,1)*clhs122 + C(0,2)*clhs123;
const double clhs125 = DN(0,0)*clhs124;
const double clhs126 = C(0,2)*clhs121 + C(1,2)*clhs122 + C(2,2)*clhs123;
const double clhs127 = DN(0,1)*clhs126;
const double clhs128 = clhs120*clhs36;
const double clhs129 = 0.5*clhs128;
const double clhs130 = clhs61*(1.0*DN(1,0)*clhs3 - clhs10*clhs129 - clhs129*clhs4);
const double clhs131 = clhs61*(1.0*DN(1,1)*clhs22 - clhs129*clhs17 - clhs129*clhs23);
const double clhs132 = clhs66*(clhs123 - clhs128*clhs26 - clhs128*clhs27);
const double clhs133 = C(0,2)*clhs130 + C(1,2)*clhs131 + C(2,2)*clhs132;
const double clhs134 = clhs133*clhs16 + clhs9*(C(0,0)*clhs130 + C(0,1)*clhs131 + C(0,2)*clhs132);
const double clhs135 = clhs133*clhs9 + clhs16*(C(0,1)*clhs130 + C(1,1)*clhs131 + C(1,2)*clhs132);
const double clhs136 = -DN(0,0)*DN(2,1)*u(0,1) - DN(1,0)*DN(2,1)*u(1,1) + DN(2,0)*clhs18 + DN(2,0)*clhs19 + DN(2,0);
const double clhs137 = clhs136*clhs36;
const double clhs138 = DN(2,1)*clhs16;
const double clhs139 = DN(2,0)*clhs9;
const double clhs140 = DN(2,0)*clhs16 + DN(2,1)*clhs9;
const double clhs141 = C(0,0)*clhs139 + C(0,1)*clhs138 + C(0,2)*clhs140;
const double clhs142 = DN(0,0)*clhs141;
const double clhs143 = C(0,2)*clhs139 + C(1,2)*clhs138 + C(2,2)*clhs140;
const double clhs144 = DN(0,1)*clhs143;
const double clhs145 = clhs136*clhs58;
const double clhs146 = -clhs136*clhs36;
const double clhs147 = 0.5*clhs146;
const double clhs148 = clhs61*(1.0*clhs138 + clhs147*clhs17 + clhs147*clhs23);
const double clhs149 = clhs61*(clhs10*clhs147 + 1.0*clhs139 + clhs147*clhs4);
const double clhs150 = clhs66*(clhs140 + clhs146*clhs26 + clhs146*clhs27);
const double clhs151 = C(0,2)*clhs149 + C(1,2)*clhs148 + C(2,2)*clhs150;
const double clhs152 = DN(2,0)*S[0] + DN(2,1)*S[2];
const double clhs153 = clhs151*clhs69 + clhs152 + clhs70*(C(0,0)*clhs149 + C(0,1)*clhs148 + C(0,2)*clhs150);
const double clhs154 = DN(2,0)*S[2] + DN(2,1)*S[1];
const double clhs155 = clhs151*clhs70 + clhs154 + clhs69*(C(0,1)*clhs149 + C(1,1)*clhs148 + C(1,2)*clhs150);
const double clhs156 = -DN(2,0)*clhs13 - DN(2,0)*clhs14 + DN(2,1)*clhs5 + DN(2,1)*clhs6 + DN(2,1);
const double clhs157 = DN(2,0)*clhs3;
const double clhs158 = DN(2,1)*clhs22;
const double clhs159 = DN(2,0)*clhs22 + DN(2,1)*clhs3;
const double clhs160 = C(0,0)*clhs157 + C(0,1)*clhs158 + C(0,2)*clhs159;
const double clhs161 = DN(0,0)*clhs160;
const double clhs162 = C(0,2)*clhs157 + C(1,2)*clhs158 + C(2,2)*clhs159;
const double clhs163 = DN(0,1)*clhs162;
const double clhs164 = clhs156*clhs36;
const double clhs165 = 0.5*clhs164;
const double clhs166 = clhs61*(1.0*DN(2,0)*clhs3 - clhs10*clhs165 - clhs165*clhs4);
const double clhs167 = clhs61*(1.0*DN(2,1)*clhs22 - clhs165*clhs17 - clhs165*clhs23);
const double clhs168 = clhs66*(clhs159 - clhs164*clhs26 - clhs164*clhs27);
const double clhs169 = C(0,2)*clhs166 + C(1,2)*clhs167 + C(2,2)*clhs168;
const double clhs170 = clhs16*clhs169 + clhs9*(C(0,0)*clhs166 + C(0,1)*clhs167 + C(0,2)*clhs168);
const double clhs171 = clhs16*(C(0,1)*clhs166 + C(1,1)*clhs167 + C(1,2)*clhs168) + clhs169*clhs9;
const double clhs172 = DN(0,0)*clhs42;
const double clhs173 = clhs172*clhs32;
const double clhs174 = clhs94 + clhs95 + 2*clhs96;
const double clhs175 = DN(0,1)*clhs174;
const double clhs176 = clhs175*clhs32;
const double clhs177 = DN(0,0)*clhs55;
const double clhs178 = C(0,1)*clhs46 + C(1,1)*clhs45 + C(1,2)*clhs47;
const double clhs179 = DN(0,1)*clhs178;
const double clhs180 = 0.5*clhs38;
const double clhs181 = clhs180*clhs53;
const double clhs182 = clhs61*(1.0*DN(0,1)*clhs16 - clhs17*clhs180 - clhs180*clhs23);
const double clhs183 = clhs61*(1.0*DN(0,0)*clhs9 - clhs10*clhs180 - clhs180*clhs4);
const double clhs184 = clhs66*(-clhs26*clhs38 - clhs27*clhs38 + clhs47);
const double clhs185 = C(0,2)*clhs183 + C(1,2)*clhs182 + C(2,2)*clhs184;
const double clhs186 = clhs185*clhs22 + clhs3*(C(0,0)*clhs183 + C(0,1)*clhs182 + C(0,2)*clhs184);
const double clhs187 = clhs185*clhs3 + clhs22*(C(0,1)*clhs183 + C(1,1)*clhs182 + C(1,2)*clhs184);
const double clhs188 = clhs84*tau_th;
const double clhs189 = DN(0,0)*clhs81;
const double clhs190 = C(0,1)*clhs76 + C(1,1)*clhs77 + C(1,2)*clhs78;
const double clhs191 = DN(0,1)*clhs190;
const double clhs192 = clhs172*clhs53;
const double clhs193 = clhs58*clhs75;
const double clhs194 = clhs175*clhs53;
const double clhs195 = -clhs36*clhs75;
const double clhs196 = 0.5*clhs195;
const double clhs197 = clhs61*(clhs10*clhs196 + clhs196*clhs4 + 1.0*clhs76);
const double clhs198 = clhs61*(clhs17*clhs196 + clhs196*clhs23 + 1.0*clhs77);
const double clhs199 = clhs66*(clhs195*clhs26 + clhs195*clhs27 + clhs78);
const double clhs200 = clhs3*clhs36;
const double clhs201 = C(0,2)*clhs197 + C(1,2)*clhs198 + C(2,2)*clhs199;
const double clhs202 = clhs22*clhs36;
const double clhs203 = clhs200*(C(0,0)*clhs197 + C(0,1)*clhs198 + C(0,2)*clhs199) + clhs201*clhs202 + clhs71;
const double clhs204 = clhs200*clhs201 + clhs202*(C(0,1)*clhs197 + C(1,1)*clhs198 + C(1,2)*clhs199) + clhs73;
const double clhs205 = clhs22*clhs91 + clhs3*clhs92;
const double clhs206 = clhs22*clhs97 + clhs3*clhs91;
const double clhs207 = clhs36*(DN(0,0)*clhs205 + DN(0,1)*clhs206 - clhs173 - clhs176);
const double clhs208 = DN(0,0)*clhs107;
const double clhs209 = C(0,1)*clhs103 + C(1,1)*clhs102 + C(1,2)*clhs104;
const double clhs210 = DN(0,1)*clhs209;
const double clhs211 = 0.5*clhs101;
const double clhs212 = clhs61*(1.0*DN(1,1)*clhs16 - clhs17*clhs211 - clhs211*clhs23);
const double clhs213 = clhs61*(1.0*DN(1,0)*clhs9 - clhs10*clhs211 - clhs211*clhs4);
const double clhs214 = clhs66*(-clhs101*clhs26 - clhs101*clhs27 + clhs104);
const double clhs215 = C(0,2)*clhs213 + C(1,2)*clhs212 + C(2,2)*clhs214;
const double clhs216 = clhs215*clhs22 + clhs3*(C(0,0)*clhs213 + C(0,1)*clhs212 + C(0,2)*clhs214);
const double clhs217 = clhs215*clhs3 + clhs22*(C(0,1)*clhs213 + C(1,1)*clhs212 + C(1,2)*clhs214);
const double clhs218 = DN(0,0)*clhs126;
const double clhs219 = C(0,1)*clhs121 + C(1,1)*clhs122 + C(1,2)*clhs123;
const double clhs220 = DN(0,1)*clhs219;
const double clhs221 = clhs120*clhs58;
const double clhs222 = -clhs120*clhs36;
const double clhs223 = 0.5*clhs222;
const double clhs224 = clhs61*(clhs10*clhs223 + 1.0*clhs121 + clhs223*clhs4);
const double clhs225 = clhs61*(1.0*clhs122 + clhs17*clhs223 + clhs223*clhs23);
const double clhs226 = clhs66*(clhs123 + clhs222*clhs26 + clhs222*clhs27);
const double clhs227 = C(0,2)*clhs224 + C(1,2)*clhs225 + C(2,2)*clhs226;
const double clhs228 = clhs116 + clhs200*(C(0,0)*clhs224 + C(0,1)*clhs225 + C(0,2)*clhs226) + clhs202*clhs227;
const double clhs229 = clhs118 + clhs200*clhs227 + clhs202*(C(0,1)*clhs224 + C(1,1)*clhs225 + C(1,2)*clhs226);
const double clhs230 = DN(0,0)*clhs143;
const double clhs231 = C(0,1)*clhs139 + C(1,1)*clhs138 + C(1,2)*clhs140;
const double clhs232 = DN(0,1)*clhs231;
const double clhs233 = 0.5*clhs137;
const double clhs234 = clhs61*(1.0*DN(2,1)*clhs16 - clhs17*clhs233 - clhs23*clhs233);
const double clhs235 = clhs61*(1.0*DN(2,0)*clhs9 - clhs10*clhs233 - clhs233*clhs4);
const double clhs236 = clhs66*(-clhs137*clhs26 - clhs137*clhs27 + clhs140);
const double clhs237 = C(0,2)*clhs235 + C(1,2)*clhs234 + C(2,2)*clhs236;
const double clhs238 = clhs22*clhs237 + clhs3*(C(0,0)*clhs235 + C(0,1)*clhs234 + C(0,2)*clhs236);
const double clhs239 = clhs22*(C(0,1)*clhs235 + C(1,1)*clhs234 + C(1,2)*clhs236) + clhs237*clhs3;
const double clhs240 = DN(0,0)*clhs162;
const double clhs241 = C(0,1)*clhs157 + C(1,1)*clhs158 + C(1,2)*clhs159;
const double clhs242 = DN(0,1)*clhs241;
const double clhs243 = clhs156*clhs58;
const double clhs244 = -clhs156*clhs36;
const double clhs245 = 0.5*clhs244;
const double clhs246 = clhs61*(clhs10*clhs245 + 1.0*clhs157 + clhs245*clhs4);
const double clhs247 = clhs61*(1.0*clhs158 + clhs17*clhs245 + clhs23*clhs245);
const double clhs248 = clhs66*(clhs159 + clhs244*clhs26 + clhs244*clhs27);
const double clhs249 = C(0,2)*clhs246 + C(1,2)*clhs247 + C(2,2)*clhs248;
const double clhs250 = clhs152 + clhs200*(C(0,0)*clhs246 + C(0,1)*clhs247 + C(0,2)*clhs248) + clhs202*clhs249;
const double clhs251 = clhs154 + clhs200*clhs249 + clhs202*(C(0,1)*clhs246 + C(1,1)*clhs247 + C(1,2)*clhs248);
const double clhs252 = tau_u*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double clhs253 = tau_u*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double clhs254 = (1.0/2.0)*tau_u;
const double clhs255 = clhs254*(clhs31 + clhs43);
const double clhs256 = clhs254*(clhs172 + clhs175);
const double clhs257 = tau_u*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs258 = b_gauss[1]*clhs257;
const double clhs259 = b_gauss[0]*clhs257;
const double clhs260 = N[0]*N[1];
const double clhs261 = tau_u*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs262 = b_gauss[1]*clhs261;
const double clhs263 = b_gauss[0]*clhs261;
const double clhs264 = N[0]*N[2];
const double clhs265 = DN(1,0)*clhs30;
const double clhs266 = clhs265*clhs32;
const double clhs267 = DN(1,1)*clhs42;
const double clhs268 = clhs267*clhs32;
const double clhs269 = DN(1,0)*clhs48;
const double clhs270 = DN(1,1)*clhs55;
const double clhs271 = clhs53*clhs59;
const double clhs272 = DN(1,0)*clhs79;
const double clhs273 = DN(1,1)*clhs81;
const double clhs274 = clhs53*clhs84;
const double clhs275 = clhs36*(DN(1,0)*clhs93 + DN(1,1)*clhs98 - clhs266 - clhs268);
const double clhs276 = DN(1,0)*clhs105;
const double clhs277 = DN(1,1)*clhs107;
const double clhs278 = clhs109*clhs53;
const double clhs279 = DN(1,0)*clhs124;
const double clhs280 = DN(1,1)*clhs126;
const double clhs281 = clhs129*clhs53;
const double clhs282 = DN(1,0)*clhs141;
const double clhs283 = DN(1,1)*clhs143;
const double clhs284 = clhs145*clhs53;
const double clhs285 = DN(1,0)*clhs160;
const double clhs286 = DN(1,1)*clhs162;
const double clhs287 = clhs165*clhs53;
const double clhs288 = DN(1,0)*clhs42;
const double clhs289 = clhs288*clhs32;
const double clhs290 = DN(1,1)*clhs174;
const double clhs291 = clhs290*clhs32;
const double clhs292 = DN(1,0)*clhs55;
const double clhs293 = DN(1,1)*clhs178;
const double clhs294 = DN(1,0)*clhs81;
const double clhs295 = DN(1,1)*clhs190;
const double clhs296 = clhs193*clhs53;
const double clhs297 = clhs36*(DN(1,0)*clhs205 + DN(1,1)*clhs206 - clhs289 - clhs291);
const double clhs298 = DN(1,0)*clhs107;
const double clhs299 = DN(1,1)*clhs209;
const double clhs300 = clhs211*clhs53;
const double clhs301 = DN(1,0)*clhs126;
const double clhs302 = DN(1,1)*clhs219;
const double clhs303 = clhs221*clhs53;
const double clhs304 = DN(1,0)*clhs143;
const double clhs305 = DN(1,1)*clhs231;
const double clhs306 = clhs233*clhs53;
const double clhs307 = DN(1,0)*clhs162;
const double clhs308 = DN(1,1)*clhs241;
const double clhs309 = clhs243*clhs53;
const double clhs310 = clhs254*(clhs265 + clhs267);
const double clhs311 = clhs254*(clhs288 + clhs290);
const double clhs312 = tau_u*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs313 = b_gauss[1]*clhs312;
const double clhs314 = b_gauss[0]*clhs312;
const double clhs315 = N[1]*N[2];
const double clhs316 = DN(2,0)*clhs30;
const double clhs317 = clhs316*clhs32;
const double clhs318 = DN(2,1)*clhs42;
const double clhs319 = clhs318*clhs32;
const double clhs320 = DN(2,0)*clhs48;
const double clhs321 = DN(2,1)*clhs55;
const double clhs322 = DN(2,0)*clhs79;
const double clhs323 = DN(2,1)*clhs81;
const double clhs324 = clhs36*(DN(2,0)*clhs93 + DN(2,1)*clhs98 - clhs317 - clhs319);
const double clhs325 = DN(2,0)*clhs105;
const double clhs326 = DN(2,1)*clhs107;
const double clhs327 = DN(2,0)*clhs124;
const double clhs328 = DN(2,1)*clhs126;
const double clhs329 = DN(2,0)*clhs141;
const double clhs330 = DN(2,1)*clhs143;
const double clhs331 = DN(2,0)*clhs160;
const double clhs332 = DN(2,1)*clhs162;
const double clhs333 = DN(2,0)*clhs42;
const double clhs334 = clhs32*clhs333;
const double clhs335 = DN(2,1)*clhs174;
const double clhs336 = clhs32*clhs335;
const double clhs337 = DN(2,0)*clhs55;
const double clhs338 = DN(2,1)*clhs178;
const double clhs339 = DN(2,0)*clhs81;
const double clhs340 = DN(2,1)*clhs190;
const double clhs341 = clhs36*(DN(2,0)*clhs205 + DN(2,1)*clhs206 - clhs334 - clhs336);
const double clhs342 = DN(2,0)*clhs107;
const double clhs343 = DN(2,1)*clhs209;
const double clhs344 = DN(2,0)*clhs126;
const double clhs345 = DN(2,1)*clhs219;
const double clhs346 = DN(2,0)*clhs143;
const double clhs347 = DN(2,1)*clhs231;
const double clhs348 = DN(2,0)*clhs162;
const double clhs349 = DN(2,1)*clhs241;
const double clhs350 = clhs254*(clhs316 + clhs318);
const double clhs351 = clhs254*(clhs333 + clhs335);
lhs(0,0)=DN(0,0)*clhs72 + DN(0,1)*clhs74 + clhs33*clhs38 + clhs38*clhs44 + clhs49*clhs54 + clhs54*clhs56 - clhs57*clhs59 - clhs59*clhs60;
lhs(0,1)=clhs36*(DN(0,0)*clhs89 + DN(0,1)*clhs90 + clhs33*clhs75 + clhs44*clhs75 + clhs53*clhs80 + clhs53*clhs82 - clhs57*clhs84 - clhs60*clhs84);
lhs(0,2)=N[0]*clhs99;
lhs(0,3)=DN(0,0)*clhs117 + DN(0,1)*clhs119 + clhs101*clhs33 + clhs101*clhs44 + clhs106*clhs54 + clhs108*clhs54 - clhs109*clhs57 - clhs109*clhs60;
lhs(0,4)=clhs36*(DN(0,0)*clhs134 + DN(0,1)*clhs135 + clhs120*clhs33 + clhs120*clhs44 + clhs125*clhs53 + clhs127*clhs53 - clhs129*clhs57 - clhs129*clhs60);
lhs(0,5)=N[1]*clhs99;
lhs(0,6)=DN(0,0)*clhs153 + DN(0,1)*clhs155 + clhs137*clhs33 + clhs137*clhs44 + clhs142*clhs54 + clhs144*clhs54 - clhs145*clhs57 - clhs145*clhs60;
lhs(0,7)=clhs36*(DN(0,0)*clhs170 + DN(0,1)*clhs171 + clhs156*clhs33 + clhs156*clhs44 + clhs161*clhs53 + clhs163*clhs53 - clhs165*clhs57 - clhs165*clhs60);
lhs(0,8)=N[2]*clhs99;
lhs(1,0)=clhs36*(DN(0,0)*clhs186 + DN(0,1)*clhs187 - clhs172*clhs181 + clhs173*clhs37 - clhs175*clhs181 + clhs176*clhs37 + clhs177*clhs53 + clhs179*clhs53);
lhs(1,1)=DN(0,0)*clhs203 + DN(0,1)*clhs204 + clhs172*clhs188 + clhs175*clhs188 + clhs189*clhs54 + clhs191*clhs54 - clhs192*clhs193 - clhs193*clhs194;
lhs(1,2)=N[0]*clhs207;
lhs(1,3)=clhs36*(DN(0,0)*clhs216 + DN(0,1)*clhs217 + clhs100*clhs173 + clhs100*clhs176 - clhs192*clhs211 - clhs194*clhs211 + clhs208*clhs53 + clhs210*clhs53);
lhs(1,4)=DN(0,0)*clhs228 + DN(0,1)*clhs229 + clhs128*clhs173 + clhs128*clhs176 - clhs192*clhs221 - clhs194*clhs221 + clhs218*clhs54 + clhs220*clhs54;
lhs(1,5)=N[1]*clhs207;
lhs(1,6)=clhs36*(DN(0,0)*clhs238 + DN(0,1)*clhs239 + clhs136*clhs173 + clhs136*clhs176 - clhs192*clhs233 - clhs194*clhs233 + clhs230*clhs53 + clhs232*clhs53);
lhs(1,7)=DN(0,0)*clhs250 + DN(0,1)*clhs251 + clhs164*clhs173 + clhs164*clhs176 - clhs192*clhs243 - clhs194*clhs243 + clhs240*clhs54 + clhs242*clhs54;
lhs(1,8)=N[2]*clhs207;
lhs(2,0)=N[0]*clhs37 - clhs252*(clhs49 + clhs56) - clhs253*(clhs177 + clhs179);
lhs(2,1)=N[0]*clhs75 - clhs252*(clhs80 + clhs82) - clhs253*(clhs189 + clhs191);
lhs(2,2)=-DN(0,0)*clhs255 - DN(0,1)*clhs256 - pow(N[0], 2);
lhs(2,3)=N[0]*clhs100 - clhs252*(clhs106 + clhs108) - clhs253*(clhs208 + clhs210) + clhs258;
lhs(2,4)=N[0]*clhs120 - clhs252*(clhs125 + clhs127) - clhs253*(clhs218 + clhs220) - clhs259;
lhs(2,5)=-DN(1,0)*clhs255 - DN(1,1)*clhs256 - clhs260;
lhs(2,6)=N[0]*clhs136 - clhs252*(clhs142 + clhs144) - clhs253*(clhs230 + clhs232) + clhs262;
lhs(2,7)=N[0]*clhs156 - clhs252*(clhs161 + clhs163) - clhs253*(clhs240 + clhs242) - clhs263;
lhs(2,8)=-DN(2,0)*clhs255 - DN(2,1)*clhs256 - clhs264;
lhs(3,0)=DN(1,0)*clhs72 + DN(1,1)*clhs74 - clhs265*clhs271 + clhs266*clhs38 - clhs267*clhs271 + clhs268*clhs38 + clhs269*clhs54 + clhs270*clhs54;
lhs(3,1)=clhs36*(DN(1,0)*clhs89 + DN(1,1)*clhs90 - clhs265*clhs274 + clhs266*clhs75 - clhs267*clhs274 + clhs268*clhs75 + clhs272*clhs53 + clhs273*clhs53);
lhs(3,2)=N[0]*clhs275;
lhs(3,3)=DN(1,0)*clhs117 + DN(1,1)*clhs119 + clhs101*clhs266 + clhs101*clhs268 - clhs265*clhs278 - clhs267*clhs278 + clhs276*clhs54 + clhs277*clhs54;
lhs(3,4)=clhs36*(DN(1,0)*clhs134 + DN(1,1)*clhs135 + clhs120*clhs266 + clhs120*clhs268 - clhs265*clhs281 - clhs267*clhs281 + clhs279*clhs53 + clhs280*clhs53);
lhs(3,5)=N[1]*clhs275;
lhs(3,6)=DN(1,0)*clhs153 + DN(1,1)*clhs155 + clhs137*clhs266 + clhs137*clhs268 - clhs265*clhs284 - clhs267*clhs284 + clhs282*clhs54 + clhs283*clhs54;
lhs(3,7)=clhs36*(DN(1,0)*clhs170 + DN(1,1)*clhs171 + clhs156*clhs266 + clhs156*clhs268 - clhs265*clhs287 - clhs267*clhs287 + clhs285*clhs53 + clhs286*clhs53);
lhs(3,8)=N[2]*clhs275;
lhs(4,0)=clhs36*(DN(1,0)*clhs186 + DN(1,1)*clhs187 - clhs181*clhs288 - clhs181*clhs290 + clhs289*clhs37 + clhs291*clhs37 + clhs292*clhs53 + clhs293*clhs53);
lhs(4,1)=DN(1,0)*clhs203 + DN(1,1)*clhs204 + clhs188*clhs288 + clhs188*clhs290 - clhs288*clhs296 - clhs290*clhs296 + clhs294*clhs54 + clhs295*clhs54;
lhs(4,2)=N[0]*clhs297;
lhs(4,3)=clhs36*(DN(1,0)*clhs216 + DN(1,1)*clhs217 + clhs100*clhs289 + clhs100*clhs291 - clhs288*clhs300 - clhs290*clhs300 + clhs298*clhs53 + clhs299*clhs53);
lhs(4,4)=DN(1,0)*clhs228 + DN(1,1)*clhs229 + clhs128*clhs289 + clhs128*clhs291 - clhs288*clhs303 - clhs290*clhs303 + clhs301*clhs54 + clhs302*clhs54;
lhs(4,5)=N[1]*clhs297;
lhs(4,6)=clhs36*(DN(1,0)*clhs238 + DN(1,1)*clhs239 + clhs136*clhs289 + clhs136*clhs291 - clhs288*clhs306 - clhs290*clhs306 + clhs304*clhs53 + clhs305*clhs53);
lhs(4,7)=DN(1,0)*clhs250 + DN(1,1)*clhs251 + clhs164*clhs289 + clhs164*clhs291 - clhs288*clhs309 - clhs290*clhs309 + clhs307*clhs54 + clhs308*clhs54;
lhs(4,8)=N[2]*clhs297;
lhs(5,0)=N[1]*clhs37 - clhs252*(clhs269 + clhs270) - clhs253*(clhs292 + clhs293) - clhs258;
lhs(5,1)=N[1]*clhs75 - clhs252*(clhs272 + clhs273) - clhs253*(clhs294 + clhs295) + clhs259;
lhs(5,2)=-DN(0,0)*clhs310 - DN(0,1)*clhs311 - clhs260;
lhs(5,3)=N[1]*clhs100 - clhs252*(clhs276 + clhs277) - clhs253*(clhs298 + clhs299);
lhs(5,4)=N[1]*clhs120 - clhs252*(clhs279 + clhs280) - clhs253*(clhs301 + clhs302);
lhs(5,5)=-DN(1,0)*clhs310 - DN(1,1)*clhs311 - pow(N[1], 2);
lhs(5,6)=N[1]*clhs136 - clhs252*(clhs282 + clhs283) - clhs253*(clhs304 + clhs305) + clhs313;
lhs(5,7)=N[1]*clhs156 - clhs252*(clhs285 + clhs286) - clhs253*(clhs307 + clhs308) - clhs314;
lhs(5,8)=-DN(2,0)*clhs310 - DN(2,1)*clhs311 - clhs315;
lhs(6,0)=DN(2,0)*clhs72 + DN(2,1)*clhs74 - clhs271*clhs316 - clhs271*clhs318 + clhs317*clhs38 + clhs319*clhs38 + clhs320*clhs54 + clhs321*clhs54;
lhs(6,1)=clhs36*(DN(2,0)*clhs89 + DN(2,1)*clhs90 - clhs274*clhs316 - clhs274*clhs318 + clhs317*clhs75 + clhs319*clhs75 + clhs322*clhs53 + clhs323*clhs53);
lhs(6,2)=N[0]*clhs324;
lhs(6,3)=DN(2,0)*clhs117 + DN(2,1)*clhs119 + clhs101*clhs317 + clhs101*clhs319 - clhs278*clhs316 - clhs278*clhs318 + clhs325*clhs54 + clhs326*clhs54;
lhs(6,4)=clhs36*(DN(2,0)*clhs134 + DN(2,1)*clhs135 + clhs120*clhs317 + clhs120*clhs319 - clhs281*clhs316 - clhs281*clhs318 + clhs327*clhs53 + clhs328*clhs53);
lhs(6,5)=N[1]*clhs324;
lhs(6,6)=DN(2,0)*clhs153 + DN(2,1)*clhs155 + clhs137*clhs317 + clhs137*clhs319 - clhs284*clhs316 - clhs284*clhs318 + clhs329*clhs54 + clhs330*clhs54;
lhs(6,7)=clhs36*(DN(2,0)*clhs170 + DN(2,1)*clhs171 + clhs156*clhs317 + clhs156*clhs319 - clhs287*clhs316 - clhs287*clhs318 + clhs331*clhs53 + clhs332*clhs53);
lhs(6,8)=N[2]*clhs324;
lhs(7,0)=clhs36*(DN(2,0)*clhs186 + DN(2,1)*clhs187 - clhs181*clhs333 - clhs181*clhs335 + clhs334*clhs37 + clhs336*clhs37 + clhs337*clhs53 + clhs338*clhs53);
lhs(7,1)=DN(2,0)*clhs203 + DN(2,1)*clhs204 + clhs188*clhs333 + clhs188*clhs335 - clhs296*clhs333 - clhs296*clhs335 + clhs339*clhs54 + clhs340*clhs54;
lhs(7,2)=N[0]*clhs341;
lhs(7,3)=clhs36*(DN(2,0)*clhs216 + DN(2,1)*clhs217 + clhs100*clhs334 + clhs100*clhs336 - clhs300*clhs333 - clhs300*clhs335 + clhs342*clhs53 + clhs343*clhs53);
lhs(7,4)=DN(2,0)*clhs228 + DN(2,1)*clhs229 + clhs128*clhs334 + clhs128*clhs336 - clhs303*clhs333 - clhs303*clhs335 + clhs344*clhs54 + clhs345*clhs54;
lhs(7,5)=N[1]*clhs341;
lhs(7,6)=clhs36*(DN(2,0)*clhs238 + DN(2,1)*clhs239 + clhs136*clhs334 + clhs136*clhs336 - clhs306*clhs333 - clhs306*clhs335 + clhs346*clhs53 + clhs347*clhs53);
lhs(7,7)=DN(2,0)*clhs250 + DN(2,1)*clhs251 + clhs164*clhs334 + clhs164*clhs336 - clhs309*clhs333 - clhs309*clhs335 + clhs348*clhs54 + clhs349*clhs54;
lhs(7,8)=N[2]*clhs341;
lhs(8,0)=N[2]*clhs37 - clhs252*(clhs320 + clhs321) - clhs253*(clhs337 + clhs338) - clhs262;
lhs(8,1)=N[2]*clhs75 - clhs252*(clhs322 + clhs323) - clhs253*(clhs339 + clhs340) + clhs263;
lhs(8,2)=-DN(0,0)*clhs350 - DN(0,1)*clhs351 - clhs264;
lhs(8,3)=N[2]*clhs100 - clhs252*(clhs325 + clhs326) - clhs253*(clhs342 + clhs343) - clhs313;
lhs(8,4)=N[2]*clhs120 - clhs252*(clhs327 + clhs328) - clhs253*(clhs344 + clhs345) + clhs314;
lhs(8,5)=-DN(1,0)*clhs350 - DN(1,1)*clhs351 - clhs315;
lhs(8,6)=N[2]*clhs136 - clhs252*(clhs329 + clhs330) - clhs253*(clhs346 + clhs347);
lhs(8,7)=N[2]*clhs156 - clhs252*(clhs331 + clhs332) - clhs253*(clhs348 + clhs349);
lhs(8,8)=-DN(2,0)*clhs350 - DN(2,1)*clhs351 - pow(N[2], 2);

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
const double crhs25 = -crhs0*crhs12 - crhs0*crhs13 - crhs1*crhs11 - crhs1*crhs13 - crhs11*crhs2 - crhs12*crhs2 + crhs16*crhs5 + crhs16*crhs6 + crhs17*crhs4 + crhs17*crhs6 + crhs18*crhs4 + crhs18*crhs5 + crhs19;
const double crhs26 = -N[0]*th[0] - N[1]*th[1] - N[2]*th[2] + crhs25 + crhs7;
const double crhs27 = 0.5*crhs26*tau_th/(crhs25 + crhs8);
const double crhs28 = C(0,2)*crhs15 + C(1,2)*crhs21 + C(2,2)*crhs22;
const double crhs29 = DN(0,1)*crhs28;
const double crhs30 = S[0]*crhs14 + S[2]*crhs20;
const double crhs31 = S[1]*crhs20 + S[2]*crhs14;
const double crhs32 = DN(0,0)*crhs28;
const double crhs33 = C(0,1)*crhs15 + C(1,1)*crhs21 + C(1,2)*crhs22;
const double crhs34 = DN(0,1)*crhs33;
const double crhs35 = b_gauss[0]*tau_u;
const double crhs36 = b_gauss[1]*tau_u;
const double crhs37 = (1.0/2.0)*tau_u;
const double crhs38 = crhs37*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs39 = crhs37*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double crhs40 = DN(1,0)*crhs23;
const double crhs41 = DN(1,1)*crhs28;
const double crhs42 = DN(1,0)*crhs28;
const double crhs43 = DN(1,1)*crhs33;
const double crhs44 = DN(2,0)*crhs23;
const double crhs45 = DN(2,1)*crhs28;
const double crhs46 = DN(2,0)*crhs28;
const double crhs47 = DN(2,1)*crhs33;
rhs[0]=-DN(0,0)*crhs9 - DN(0,1)*crhs10 + N[0]*b_gauss[0] - crhs24*crhs27 - crhs27*crhs29;
rhs[1]=-DN(0,0)*crhs30 - DN(0,1)*crhs31 + N[0]*b_gauss[1] - crhs27*crhs32 - crhs27*crhs34;
rhs[2]=-N[0]*crhs26 + crhs35*(DN(0,0)*crhs20 - DN(0,1)*crhs14) - crhs36*(DN(0,0)*crhs3 - DN(0,1)*crhs8) + crhs38*(crhs24 + crhs29) + crhs39*(crhs32 + crhs34);
rhs[3]=-DN(1,0)*crhs9 - DN(1,1)*crhs10 + N[1]*b_gauss[0] - crhs27*crhs40 - crhs27*crhs41;
rhs[4]=-DN(1,0)*crhs30 - DN(1,1)*crhs31 + N[1]*b_gauss[1] - crhs27*crhs42 - crhs27*crhs43;
rhs[5]=-N[1]*crhs26 + crhs35*(DN(1,0)*crhs20 - DN(1,1)*crhs14) - crhs36*(DN(1,0)*crhs3 - DN(1,1)*crhs8) + crhs38*(crhs40 + crhs41) + crhs39*(crhs42 + crhs43);
rhs[6]=-DN(2,0)*crhs9 - DN(2,1)*crhs10 + N[2]*b_gauss[0] - crhs27*crhs44 - crhs27*crhs45;
rhs[7]=-DN(2,0)*crhs30 - DN(2,1)*crhs31 + N[2]*b_gauss[1] - crhs27*crhs46 - crhs27*crhs47;
rhs[8]=-N[2]*crhs26 + crhs35*(DN(2,0)*crhs20 - DN(2,1)*crhs14) - crhs36*(DN(2,0)*crhs3 - DN(2,1)*crhs8) + crhs38*(crhs44 + crhs45) + crhs39*(crhs46 + crhs47);

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
const double clhs29 = C(0,2)*clhs28;
const double clhs30 = clhs12 + clhs25 + 2*clhs29;
const double clhs31 = DN(0,0)*clhs30;
const double clhs32 = 0.5*tau_th;
const double clhs33 = clhs31*clhs32;
const double clhs34 = -clhs0*clhs14 - clhs0*clhs15 - clhs1*clhs13 - clhs1*clhs15 - clhs13*clhs2 - clhs14*clhs2 + clhs18*clhs6 + clhs18*clhs7 + clhs19*clhs5 + clhs19*clhs7 + clhs20*clhs5 + clhs20*clhs6 + clhs21;
const double clhs35 = clhs34 + clhs9;
const double clhs36 = 1.0/clhs35;
const double clhs37 = DN(0,0)*clhs19 + DN(0,0)*clhs20 + DN(0,0) - DN(0,1)*DN(1,0)*u(1,1) - DN(0,1)*DN(2,0)*u(2,1);
const double clhs38 = clhs36*clhs37;
const double clhs39 = C(0,2)*clhs11;
const double clhs40 = C(1,2)*clhs24;
const double clhs41 = C(2,2)*clhs28;
const double clhs42 = clhs39 + clhs40 + 2*clhs41;
const double clhs43 = DN(0,1)*clhs42;
const double clhs44 = clhs32*clhs43;
const double clhs45 = DN(0,1)*clhs16;
const double clhs46 = DN(0,0)*clhs9;
const double clhs47 = DN(0,0)*clhs16 + DN(0,1)*clhs9;
const double clhs48 = C(0,0)*clhs46 + C(0,1)*clhs45 + C(0,2)*clhs47;
const double clhs49 = DN(0,0)*clhs48;
const double clhs50 = N[0]*th[0];
const double clhs51 = N[1]*th[1];
const double clhs52 = N[2]*th[2];
const double clhs53 = tau_th*(clhs34 - clhs50 - clhs51 - clhs52 + clhs8);
const double clhs54 = clhs36*clhs53;
const double clhs55 = C(0,2)*clhs46 + C(1,2)*clhs45 + C(2,2)*clhs47;
const double clhs56 = DN(0,1)*clhs55;
const double clhs57 = clhs31*clhs53;
const double clhs58 = 0.5/pow(clhs35, 2);
const double clhs59 = clhs37*clhs58;
const double clhs60 = clhs43*clhs53;
const double clhs61 = clhs50 + clhs51 + clhs52 + 1.0;
const double clhs62 = -clhs36*clhs37;
const double clhs63 = 0.5*clhs62;
const double clhs64 = clhs61*(clhs17*clhs63 + clhs23*clhs63 + 1.0*clhs45);
const double clhs65 = clhs61*(clhs10*clhs63 + clhs4*clhs63 + 1.0*clhs46);
const double clhs66 = clhs50 + clhs51 + clhs52 + 1.0;
const double clhs67 = clhs66*(clhs26*clhs62 + clhs27*clhs62 + clhs47);
const double clhs68 = C(0,2)*clhs65 + C(1,2)*clhs64 + C(2,2)*clhs67;
const double clhs69 = clhs16*clhs36;
const double clhs70 = clhs36*clhs9;
const double clhs71 = DN(0,0)*S[0] + DN(0,1)*S[2];
const double clhs72 = clhs68*clhs69 + clhs70*(C(0,0)*clhs65 + C(0,1)*clhs64 + C(0,2)*clhs67) + clhs71;
const double clhs73 = DN(0,0)*S[2] + DN(0,1)*S[1];
const double clhs74 = clhs68*clhs70 + clhs69*(C(0,1)*clhs65 + C(1,1)*clhs64 + C(1,2)*clhs67) + clhs73;
const double clhs75 = -DN(0,0)*clhs14 - DN(0,0)*clhs15 + DN(0,1)*clhs6 + DN(0,1)*clhs7 + DN(0,1);
const double clhs76 = DN(0,0)*clhs3;
const double clhs77 = DN(0,1)*clhs22;
const double clhs78 = DN(0,0)*clhs22 + DN(0,1)*clhs3;
const double clhs79 = C(0,0)*clhs76 + C(0,1)*clhs77 + C(0,2)*clhs78;
const double clhs80 = DN(0,0)*clhs79;
const double clhs81 = C(0,2)*clhs76 + C(1,2)*clhs77 + C(2,2)*clhs78;
const double clhs82 = DN(0,1)*clhs81;
const double clhs83 = clhs36*clhs75;
const double clhs84 = 0.5*clhs83;
const double clhs85 = clhs61*(1.0*DN(0,0)*clhs3 - clhs10*clhs84 - clhs4*clhs84);
const double clhs86 = clhs61*(1.0*DN(0,1)*clhs22 - clhs17*clhs84 - clhs23*clhs84);
const double clhs87 = clhs66*(-clhs26*clhs83 - clhs27*clhs83 + clhs78);
const double clhs88 = C(0,2)*clhs85 + C(1,2)*clhs86 + C(2,2)*clhs87;
const double clhs89 = clhs16*clhs88 + clhs9*(C(0,0)*clhs85 + C(0,1)*clhs86 + C(0,2)*clhs87);
const double clhs90 = clhs16*(C(0,1)*clhs85 + C(1,1)*clhs86 + C(1,2)*clhs87) + clhs88*clhs9;
const double clhs91 = 0.5*clhs39 + 0.5*clhs40 + clhs41;
const double clhs92 = 0.5*clhs12 + 0.5*clhs25 + clhs29;
const double clhs93 = clhs16*clhs91 + clhs9*clhs92;
const double clhs94 = C(0,1)*clhs11;
const double clhs95 = C(1,1)*clhs24;
const double clhs96 = C(1,2)*clhs28;
const double clhs97 = 0.5*clhs94 + 0.5*clhs95 + clhs96;
const double clhs98 = clhs16*clhs97 + clhs9*clhs91;
const double clhs99 = clhs36*(DN(0,0)*clhs93 + DN(0,1)*clhs98 - clhs33 - clhs44);
const double clhs100 = -DN(0,0)*DN(1,1)*u(0,1) + DN(1,0)*clhs18 + DN(1,0)*clhs20 + DN(1,0) - DN(1,1)*DN(2,0)*u(2,1);
const double clhs101 = clhs100*clhs36;
const double clhs102 = DN(1,1)*clhs16;
const double clhs103 = DN(1,0)*clhs9;
const double clhs104 = DN(1,0)*clhs16 + DN(1,1)*clhs9;
const double clhs105 = C(0,0)*clhs103 + C(0,1)*clhs102 + C(0,2)*clhs104;
const double clhs106 = DN(0,0)*clhs105;
const double clhs107 = C(0,2)*clhs103 + C(1,2)*clhs102 + C(2,2)*clhs104;
const double clhs108 = DN(0,1)*clhs107;
const double clhs109 = clhs100*clhs58;
const double clhs110 = -clhs100*clhs36;
const double clhs111 = 0.5*clhs110;
const double clhs112 = clhs61*(1.0*clhs102 + clhs111*clhs17 + clhs111*clhs23);
const double clhs113 = clhs61*(clhs10*clhs111 + 1.0*clhs103 + clhs111*clhs4);
const double clhs114 = clhs66*(clhs104 + clhs110*clhs26 + clhs110*clhs27);
const double clhs115 = C(0,2)*clhs113 + C(1,2)*clhs112 + C(2,2)*clhs114;
const double clhs116 = DN(1,0)*S[0] + DN(1,1)*S[2];
const double clhs117 = clhs115*clhs69 + clhs116 + clhs70*(C(0,0)*clhs113 + C(0,1)*clhs112 + C(0,2)*clhs114);
const double clhs118 = DN(1,0)*S[2] + DN(1,1)*S[1];
const double clhs119 = clhs115*clhs70 + clhs118 + clhs69*(C(0,1)*clhs113 + C(1,1)*clhs112 + C(1,2)*clhs114);
const double clhs120 = -DN(1,0)*clhs13 - DN(1,0)*clhs15 + DN(1,1)*clhs5 + DN(1,1)*clhs7 + DN(1,1);
const double clhs121 = DN(1,0)*clhs3;
const double clhs122 = DN(1,1)*clhs22;
const double clhs123 = DN(1,0)*clhs22 + DN(1,1)*clhs3;
const double clhs124 = C(0,0)*clhs121 + C(0,1)*clhs122 + C(0,2)*clhs123;
const double clhs125 = DN(0,0)*clhs124;
const double clhs126 = C(0,2)*clhs121 + C(1,2)*clhs122 + C(2,2)*clhs123;
const double clhs127 = DN(0,1)*clhs126;
const double clhs128 = clhs120*clhs36;
const double clhs129 = 0.5*clhs128;
const double clhs130 = clhs61*(1.0*DN(1,0)*clhs3 - clhs10*clhs129 - clhs129*clhs4);
const double clhs131 = clhs61*(1.0*DN(1,1)*clhs22 - clhs129*clhs17 - clhs129*clhs23);
const double clhs132 = clhs66*(clhs123 - clhs128*clhs26 - clhs128*clhs27);
const double clhs133 = C(0,2)*clhs130 + C(1,2)*clhs131 + C(2,2)*clhs132;
const double clhs134 = clhs133*clhs16 + clhs9*(C(0,0)*clhs130 + C(0,1)*clhs131 + C(0,2)*clhs132);
const double clhs135 = clhs133*clhs9 + clhs16*(C(0,1)*clhs130 + C(1,1)*clhs131 + C(1,2)*clhs132);
const double clhs136 = -DN(0,0)*DN(2,1)*u(0,1) - DN(1,0)*DN(2,1)*u(1,1) + DN(2,0)*clhs18 + DN(2,0)*clhs19 + DN(2,0);
const double clhs137 = clhs136*clhs36;
const double clhs138 = DN(2,1)*clhs16;
const double clhs139 = DN(2,0)*clhs9;
const double clhs140 = DN(2,0)*clhs16 + DN(2,1)*clhs9;
const double clhs141 = C(0,0)*clhs139 + C(0,1)*clhs138 + C(0,2)*clhs140;
const double clhs142 = DN(0,0)*clhs141;
const double clhs143 = C(0,2)*clhs139 + C(1,2)*clhs138 + C(2,2)*clhs140;
const double clhs144 = DN(0,1)*clhs143;
const double clhs145 = clhs136*clhs58;
const double clhs146 = -clhs136*clhs36;
const double clhs147 = 0.5*clhs146;
const double clhs148 = clhs61*(1.0*clhs138 + clhs147*clhs17 + clhs147*clhs23);
const double clhs149 = clhs61*(clhs10*clhs147 + 1.0*clhs139 + clhs147*clhs4);
const double clhs150 = clhs66*(clhs140 + clhs146*clhs26 + clhs146*clhs27);
const double clhs151 = C(0,2)*clhs149 + C(1,2)*clhs148 + C(2,2)*clhs150;
const double clhs152 = DN(2,0)*S[0] + DN(2,1)*S[2];
const double clhs153 = clhs151*clhs69 + clhs152 + clhs70*(C(0,0)*clhs149 + C(0,1)*clhs148 + C(0,2)*clhs150);
const double clhs154 = DN(2,0)*S[2] + DN(2,1)*S[1];
const double clhs155 = clhs151*clhs70 + clhs154 + clhs69*(C(0,1)*clhs149 + C(1,1)*clhs148 + C(1,2)*clhs150);
const double clhs156 = -DN(2,0)*clhs13 - DN(2,0)*clhs14 + DN(2,1)*clhs5 + DN(2,1)*clhs6 + DN(2,1);
const double clhs157 = DN(2,0)*clhs3;
const double clhs158 = DN(2,1)*clhs22;
const double clhs159 = DN(2,0)*clhs22 + DN(2,1)*clhs3;
const double clhs160 = C(0,0)*clhs157 + C(0,1)*clhs158 + C(0,2)*clhs159;
const double clhs161 = DN(0,0)*clhs160;
const double clhs162 = C(0,2)*clhs157 + C(1,2)*clhs158 + C(2,2)*clhs159;
const double clhs163 = DN(0,1)*clhs162;
const double clhs164 = clhs156*clhs36;
const double clhs165 = 0.5*clhs164;
const double clhs166 = clhs61*(1.0*DN(2,0)*clhs3 - clhs10*clhs165 - clhs165*clhs4);
const double clhs167 = clhs61*(1.0*DN(2,1)*clhs22 - clhs165*clhs17 - clhs165*clhs23);
const double clhs168 = clhs66*(clhs159 - clhs164*clhs26 - clhs164*clhs27);
const double clhs169 = C(0,2)*clhs166 + C(1,2)*clhs167 + C(2,2)*clhs168;
const double clhs170 = clhs16*clhs169 + clhs9*(C(0,0)*clhs166 + C(0,1)*clhs167 + C(0,2)*clhs168);
const double clhs171 = clhs16*(C(0,1)*clhs166 + C(1,1)*clhs167 + C(1,2)*clhs168) + clhs169*clhs9;
const double clhs172 = DN(0,0)*clhs42;
const double clhs173 = clhs172*clhs32;
const double clhs174 = clhs94 + clhs95 + 2*clhs96;
const double clhs175 = DN(0,1)*clhs174;
const double clhs176 = clhs175*clhs32;
const double clhs177 = DN(0,0)*clhs55;
const double clhs178 = C(0,1)*clhs46 + C(1,1)*clhs45 + C(1,2)*clhs47;
const double clhs179 = DN(0,1)*clhs178;
const double clhs180 = 0.5*clhs38;
const double clhs181 = clhs180*clhs53;
const double clhs182 = clhs61*(1.0*DN(0,1)*clhs16 - clhs17*clhs180 - clhs180*clhs23);
const double clhs183 = clhs61*(1.0*DN(0,0)*clhs9 - clhs10*clhs180 - clhs180*clhs4);
const double clhs184 = clhs66*(-clhs26*clhs38 - clhs27*clhs38 + clhs47);
const double clhs185 = C(0,2)*clhs183 + C(1,2)*clhs182 + C(2,2)*clhs184;
const double clhs186 = clhs185*clhs22 + clhs3*(C(0,0)*clhs183 + C(0,1)*clhs182 + C(0,2)*clhs184);
const double clhs187 = clhs185*clhs3 + clhs22*(C(0,1)*clhs183 + C(1,1)*clhs182 + C(1,2)*clhs184);
const double clhs188 = clhs84*tau_th;
const double clhs189 = DN(0,0)*clhs81;
const double clhs190 = C(0,1)*clhs76 + C(1,1)*clhs77 + C(1,2)*clhs78;
const double clhs191 = DN(0,1)*clhs190;
const double clhs192 = clhs172*clhs53;
const double clhs193 = clhs58*clhs75;
const double clhs194 = clhs175*clhs53;
const double clhs195 = -clhs36*clhs75;
const double clhs196 = 0.5*clhs195;
const double clhs197 = clhs61*(clhs10*clhs196 + clhs196*clhs4 + 1.0*clhs76);
const double clhs198 = clhs61*(clhs17*clhs196 + clhs196*clhs23 + 1.0*clhs77);
const double clhs199 = clhs66*(clhs195*clhs26 + clhs195*clhs27 + clhs78);
const double clhs200 = clhs3*clhs36;
const double clhs201 = C(0,2)*clhs197 + C(1,2)*clhs198 + C(2,2)*clhs199;
const double clhs202 = clhs22*clhs36;
const double clhs203 = clhs200*(C(0,0)*clhs197 + C(0,1)*clhs198 + C(0,2)*clhs199) + clhs201*clhs202 + clhs71;
const double clhs204 = clhs200*clhs201 + clhs202*(C(0,1)*clhs197 + C(1,1)*clhs198 + C(1,2)*clhs199) + clhs73;
const double clhs205 = clhs22*clhs91 + clhs3*clhs92;
const double clhs206 = clhs22*clhs97 + clhs3*clhs91;
const double clhs207 = clhs36*(DN(0,0)*clhs205 + DN(0,1)*clhs206 - clhs173 - clhs176);
const double clhs208 = DN(0,0)*clhs107;
const double clhs209 = C(0,1)*clhs103 + C(1,1)*clhs102 + C(1,2)*clhs104;
const double clhs210 = DN(0,1)*clhs209;
const double clhs211 = 0.5*clhs101;
const double clhs212 = clhs61*(1.0*DN(1,1)*clhs16 - clhs17*clhs211 - clhs211*clhs23);
const double clhs213 = clhs61*(1.0*DN(1,0)*clhs9 - clhs10*clhs211 - clhs211*clhs4);
const double clhs214 = clhs66*(-clhs101*clhs26 - clhs101*clhs27 + clhs104);
const double clhs215 = C(0,2)*clhs213 + C(1,2)*clhs212 + C(2,2)*clhs214;
const double clhs216 = clhs215*clhs22 + clhs3*(C(0,0)*clhs213 + C(0,1)*clhs212 + C(0,2)*clhs214);
const double clhs217 = clhs215*clhs3 + clhs22*(C(0,1)*clhs213 + C(1,1)*clhs212 + C(1,2)*clhs214);
const double clhs218 = DN(0,0)*clhs126;
const double clhs219 = C(0,1)*clhs121 + C(1,1)*clhs122 + C(1,2)*clhs123;
const double clhs220 = DN(0,1)*clhs219;
const double clhs221 = clhs120*clhs58;
const double clhs222 = -clhs120*clhs36;
const double clhs223 = 0.5*clhs222;
const double clhs224 = clhs61*(clhs10*clhs223 + 1.0*clhs121 + clhs223*clhs4);
const double clhs225 = clhs61*(1.0*clhs122 + clhs17*clhs223 + clhs223*clhs23);
const double clhs226 = clhs66*(clhs123 + clhs222*clhs26 + clhs222*clhs27);
const double clhs227 = C(0,2)*clhs224 + C(1,2)*clhs225 + C(2,2)*clhs226;
const double clhs228 = clhs116 + clhs200*(C(0,0)*clhs224 + C(0,1)*clhs225 + C(0,2)*clhs226) + clhs202*clhs227;
const double clhs229 = clhs118 + clhs200*clhs227 + clhs202*(C(0,1)*clhs224 + C(1,1)*clhs225 + C(1,2)*clhs226);
const double clhs230 = DN(0,0)*clhs143;
const double clhs231 = C(0,1)*clhs139 + C(1,1)*clhs138 + C(1,2)*clhs140;
const double clhs232 = DN(0,1)*clhs231;
const double clhs233 = 0.5*clhs137;
const double clhs234 = clhs61*(1.0*DN(2,1)*clhs16 - clhs17*clhs233 - clhs23*clhs233);
const double clhs235 = clhs61*(1.0*DN(2,0)*clhs9 - clhs10*clhs233 - clhs233*clhs4);
const double clhs236 = clhs66*(-clhs137*clhs26 - clhs137*clhs27 + clhs140);
const double clhs237 = C(0,2)*clhs235 + C(1,2)*clhs234 + C(2,2)*clhs236;
const double clhs238 = clhs22*clhs237 + clhs3*(C(0,0)*clhs235 + C(0,1)*clhs234 + C(0,2)*clhs236);
const double clhs239 = clhs22*(C(0,1)*clhs235 + C(1,1)*clhs234 + C(1,2)*clhs236) + clhs237*clhs3;
const double clhs240 = DN(0,0)*clhs162;
const double clhs241 = C(0,1)*clhs157 + C(1,1)*clhs158 + C(1,2)*clhs159;
const double clhs242 = DN(0,1)*clhs241;
const double clhs243 = clhs156*clhs58;
const double clhs244 = -clhs156*clhs36;
const double clhs245 = 0.5*clhs244;
const double clhs246 = clhs61*(clhs10*clhs245 + 1.0*clhs157 + clhs245*clhs4);
const double clhs247 = clhs61*(1.0*clhs158 + clhs17*clhs245 + clhs23*clhs245);
const double clhs248 = clhs66*(clhs159 + clhs244*clhs26 + clhs244*clhs27);
const double clhs249 = C(0,2)*clhs246 + C(1,2)*clhs247 + C(2,2)*clhs248;
const double clhs250 = clhs152 + clhs200*(C(0,0)*clhs246 + C(0,1)*clhs247 + C(0,2)*clhs248) + clhs202*clhs249;
const double clhs251 = clhs154 + clhs200*clhs249 + clhs202*(C(0,1)*clhs246 + C(1,1)*clhs247 + C(1,2)*clhs248);
const double clhs252 = tau_u*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double clhs253 = tau_u*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double clhs254 = (1.0/2.0)*tau_u;
const double clhs255 = clhs254*(clhs31 + clhs43);
const double clhs256 = clhs254*(clhs172 + clhs175);
const double clhs257 = tau_u*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs258 = b_gauss[1]*clhs257;
const double clhs259 = b_gauss[0]*clhs257;
const double clhs260 = N[0]*N[1];
const double clhs261 = tau_u*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs262 = b_gauss[1]*clhs261;
const double clhs263 = b_gauss[0]*clhs261;
const double clhs264 = N[0]*N[2];
const double clhs265 = DN(1,0)*clhs30;
const double clhs266 = clhs265*clhs32;
const double clhs267 = DN(1,1)*clhs42;
const double clhs268 = clhs267*clhs32;
const double clhs269 = DN(1,0)*clhs48;
const double clhs270 = DN(1,1)*clhs55;
const double clhs271 = clhs53*clhs59;
const double clhs272 = DN(1,0)*clhs79;
const double clhs273 = DN(1,1)*clhs81;
const double clhs274 = clhs53*clhs84;
const double clhs275 = clhs36*(DN(1,0)*clhs93 + DN(1,1)*clhs98 - clhs266 - clhs268);
const double clhs276 = DN(1,0)*clhs105;
const double clhs277 = DN(1,1)*clhs107;
const double clhs278 = clhs109*clhs53;
const double clhs279 = DN(1,0)*clhs124;
const double clhs280 = DN(1,1)*clhs126;
const double clhs281 = clhs129*clhs53;
const double clhs282 = DN(1,0)*clhs141;
const double clhs283 = DN(1,1)*clhs143;
const double clhs284 = clhs145*clhs53;
const double clhs285 = DN(1,0)*clhs160;
const double clhs286 = DN(1,1)*clhs162;
const double clhs287 = clhs165*clhs53;
const double clhs288 = DN(1,0)*clhs42;
const double clhs289 = clhs288*clhs32;
const double clhs290 = DN(1,1)*clhs174;
const double clhs291 = clhs290*clhs32;
const double clhs292 = DN(1,0)*clhs55;
const double clhs293 = DN(1,1)*clhs178;
const double clhs294 = DN(1,0)*clhs81;
const double clhs295 = DN(1,1)*clhs190;
const double clhs296 = clhs193*clhs53;
const double clhs297 = clhs36*(DN(1,0)*clhs205 + DN(1,1)*clhs206 - clhs289 - clhs291);
const double clhs298 = DN(1,0)*clhs107;
const double clhs299 = DN(1,1)*clhs209;
const double clhs300 = clhs211*clhs53;
const double clhs301 = DN(1,0)*clhs126;
const double clhs302 = DN(1,1)*clhs219;
const double clhs303 = clhs221*clhs53;
const double clhs304 = DN(1,0)*clhs143;
const double clhs305 = DN(1,1)*clhs231;
const double clhs306 = clhs233*clhs53;
const double clhs307 = DN(1,0)*clhs162;
const double clhs308 = DN(1,1)*clhs241;
const double clhs309 = clhs243*clhs53;
const double clhs310 = clhs254*(clhs265 + clhs267);
const double clhs311 = clhs254*(clhs288 + clhs290);
const double clhs312 = tau_u*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs313 = b_gauss[1]*clhs312;
const double clhs314 = b_gauss[0]*clhs312;
const double clhs315 = N[1]*N[2];
const double clhs316 = DN(2,0)*clhs30;
const double clhs317 = clhs316*clhs32;
const double clhs318 = DN(2,1)*clhs42;
const double clhs319 = clhs318*clhs32;
const double clhs320 = DN(2,0)*clhs48;
const double clhs321 = DN(2,1)*clhs55;
const double clhs322 = DN(2,0)*clhs79;
const double clhs323 = DN(2,1)*clhs81;
const double clhs324 = clhs36*(DN(2,0)*clhs93 + DN(2,1)*clhs98 - clhs317 - clhs319);
const double clhs325 = DN(2,0)*clhs105;
const double clhs326 = DN(2,1)*clhs107;
const double clhs327 = DN(2,0)*clhs124;
const double clhs328 = DN(2,1)*clhs126;
const double clhs329 = DN(2,0)*clhs141;
const double clhs330 = DN(2,1)*clhs143;
const double clhs331 = DN(2,0)*clhs160;
const double clhs332 = DN(2,1)*clhs162;
const double clhs333 = DN(2,0)*clhs42;
const double clhs334 = clhs32*clhs333;
const double clhs335 = DN(2,1)*clhs174;
const double clhs336 = clhs32*clhs335;
const double clhs337 = DN(2,0)*clhs55;
const double clhs338 = DN(2,1)*clhs178;
const double clhs339 = DN(2,0)*clhs81;
const double clhs340 = DN(2,1)*clhs190;
const double clhs341 = clhs36*(DN(2,0)*clhs205 + DN(2,1)*clhs206 - clhs334 - clhs336);
const double clhs342 = DN(2,0)*clhs107;
const double clhs343 = DN(2,1)*clhs209;
const double clhs344 = DN(2,0)*clhs126;
const double clhs345 = DN(2,1)*clhs219;
const double clhs346 = DN(2,0)*clhs143;
const double clhs347 = DN(2,1)*clhs231;
const double clhs348 = DN(2,0)*clhs162;
const double clhs349 = DN(2,1)*clhs241;
const double clhs350 = clhs254*(clhs316 + clhs318);
const double clhs351 = clhs254*(clhs333 + clhs335);
lhs(0,0)=DN(0,0)*clhs72 + DN(0,1)*clhs74 + clhs33*clhs38 + clhs38*clhs44 + clhs49*clhs54 + clhs54*clhs56 - clhs57*clhs59 - clhs59*clhs60;
lhs(0,1)=clhs36*(DN(0,0)*clhs89 + DN(0,1)*clhs90 + clhs33*clhs75 + clhs44*clhs75 + clhs53*clhs80 + clhs53*clhs82 - clhs57*clhs84 - clhs60*clhs84);
lhs(0,2)=N[0]*clhs99;
lhs(0,3)=DN(0,0)*clhs117 + DN(0,1)*clhs119 + clhs101*clhs33 + clhs101*clhs44 + clhs106*clhs54 + clhs108*clhs54 - clhs109*clhs57 - clhs109*clhs60;
lhs(0,4)=clhs36*(DN(0,0)*clhs134 + DN(0,1)*clhs135 + clhs120*clhs33 + clhs120*clhs44 + clhs125*clhs53 + clhs127*clhs53 - clhs129*clhs57 - clhs129*clhs60);
lhs(0,5)=N[1]*clhs99;
lhs(0,6)=DN(0,0)*clhs153 + DN(0,1)*clhs155 + clhs137*clhs33 + clhs137*clhs44 + clhs142*clhs54 + clhs144*clhs54 - clhs145*clhs57 - clhs145*clhs60;
lhs(0,7)=clhs36*(DN(0,0)*clhs170 + DN(0,1)*clhs171 + clhs156*clhs33 + clhs156*clhs44 + clhs161*clhs53 + clhs163*clhs53 - clhs165*clhs57 - clhs165*clhs60);
lhs(0,8)=N[2]*clhs99;
lhs(1,0)=clhs36*(DN(0,0)*clhs186 + DN(0,1)*clhs187 - clhs172*clhs181 + clhs173*clhs37 - clhs175*clhs181 + clhs176*clhs37 + clhs177*clhs53 + clhs179*clhs53);
lhs(1,1)=DN(0,0)*clhs203 + DN(0,1)*clhs204 + clhs172*clhs188 + clhs175*clhs188 + clhs189*clhs54 + clhs191*clhs54 - clhs192*clhs193 - clhs193*clhs194;
lhs(1,2)=N[0]*clhs207;
lhs(1,3)=clhs36*(DN(0,0)*clhs216 + DN(0,1)*clhs217 + clhs100*clhs173 + clhs100*clhs176 - clhs192*clhs211 - clhs194*clhs211 + clhs208*clhs53 + clhs210*clhs53);
lhs(1,4)=DN(0,0)*clhs228 + DN(0,1)*clhs229 + clhs128*clhs173 + clhs128*clhs176 - clhs192*clhs221 - clhs194*clhs221 + clhs218*clhs54 + clhs220*clhs54;
lhs(1,5)=N[1]*clhs207;
lhs(1,6)=clhs36*(DN(0,0)*clhs238 + DN(0,1)*clhs239 + clhs136*clhs173 + clhs136*clhs176 - clhs192*clhs233 - clhs194*clhs233 + clhs230*clhs53 + clhs232*clhs53);
lhs(1,7)=DN(0,0)*clhs250 + DN(0,1)*clhs251 + clhs164*clhs173 + clhs164*clhs176 - clhs192*clhs243 - clhs194*clhs243 + clhs240*clhs54 + clhs242*clhs54;
lhs(1,8)=N[2]*clhs207;
lhs(2,0)=N[0]*clhs37 - clhs252*(clhs49 + clhs56) - clhs253*(clhs177 + clhs179);
lhs(2,1)=N[0]*clhs75 - clhs252*(clhs80 + clhs82) - clhs253*(clhs189 + clhs191);
lhs(2,2)=-DN(0,0)*clhs255 - DN(0,1)*clhs256 - pow(N[0], 2);
lhs(2,3)=N[0]*clhs100 - clhs252*(clhs106 + clhs108) - clhs253*(clhs208 + clhs210) + clhs258;
lhs(2,4)=N[0]*clhs120 - clhs252*(clhs125 + clhs127) - clhs253*(clhs218 + clhs220) - clhs259;
lhs(2,5)=-DN(1,0)*clhs255 - DN(1,1)*clhs256 - clhs260;
lhs(2,6)=N[0]*clhs136 - clhs252*(clhs142 + clhs144) - clhs253*(clhs230 + clhs232) + clhs262;
lhs(2,7)=N[0]*clhs156 - clhs252*(clhs161 + clhs163) - clhs253*(clhs240 + clhs242) - clhs263;
lhs(2,8)=-DN(2,0)*clhs255 - DN(2,1)*clhs256 - clhs264;
lhs(3,0)=DN(1,0)*clhs72 + DN(1,1)*clhs74 - clhs265*clhs271 + clhs266*clhs38 - clhs267*clhs271 + clhs268*clhs38 + clhs269*clhs54 + clhs270*clhs54;
lhs(3,1)=clhs36*(DN(1,0)*clhs89 + DN(1,1)*clhs90 - clhs265*clhs274 + clhs266*clhs75 - clhs267*clhs274 + clhs268*clhs75 + clhs272*clhs53 + clhs273*clhs53);
lhs(3,2)=N[0]*clhs275;
lhs(3,3)=DN(1,0)*clhs117 + DN(1,1)*clhs119 + clhs101*clhs266 + clhs101*clhs268 - clhs265*clhs278 - clhs267*clhs278 + clhs276*clhs54 + clhs277*clhs54;
lhs(3,4)=clhs36*(DN(1,0)*clhs134 + DN(1,1)*clhs135 + clhs120*clhs266 + clhs120*clhs268 - clhs265*clhs281 - clhs267*clhs281 + clhs279*clhs53 + clhs280*clhs53);
lhs(3,5)=N[1]*clhs275;
lhs(3,6)=DN(1,0)*clhs153 + DN(1,1)*clhs155 + clhs137*clhs266 + clhs137*clhs268 - clhs265*clhs284 - clhs267*clhs284 + clhs282*clhs54 + clhs283*clhs54;
lhs(3,7)=clhs36*(DN(1,0)*clhs170 + DN(1,1)*clhs171 + clhs156*clhs266 + clhs156*clhs268 - clhs265*clhs287 - clhs267*clhs287 + clhs285*clhs53 + clhs286*clhs53);
lhs(3,8)=N[2]*clhs275;
lhs(4,0)=clhs36*(DN(1,0)*clhs186 + DN(1,1)*clhs187 - clhs181*clhs288 - clhs181*clhs290 + clhs289*clhs37 + clhs291*clhs37 + clhs292*clhs53 + clhs293*clhs53);
lhs(4,1)=DN(1,0)*clhs203 + DN(1,1)*clhs204 + clhs188*clhs288 + clhs188*clhs290 - clhs288*clhs296 - clhs290*clhs296 + clhs294*clhs54 + clhs295*clhs54;
lhs(4,2)=N[0]*clhs297;
lhs(4,3)=clhs36*(DN(1,0)*clhs216 + DN(1,1)*clhs217 + clhs100*clhs289 + clhs100*clhs291 - clhs288*clhs300 - clhs290*clhs300 + clhs298*clhs53 + clhs299*clhs53);
lhs(4,4)=DN(1,0)*clhs228 + DN(1,1)*clhs229 + clhs128*clhs289 + clhs128*clhs291 - clhs288*clhs303 - clhs290*clhs303 + clhs301*clhs54 + clhs302*clhs54;
lhs(4,5)=N[1]*clhs297;
lhs(4,6)=clhs36*(DN(1,0)*clhs238 + DN(1,1)*clhs239 + clhs136*clhs289 + clhs136*clhs291 - clhs288*clhs306 - clhs290*clhs306 + clhs304*clhs53 + clhs305*clhs53);
lhs(4,7)=DN(1,0)*clhs250 + DN(1,1)*clhs251 + clhs164*clhs289 + clhs164*clhs291 - clhs288*clhs309 - clhs290*clhs309 + clhs307*clhs54 + clhs308*clhs54;
lhs(4,8)=N[2]*clhs297;
lhs(5,0)=N[1]*clhs37 - clhs252*(clhs269 + clhs270) - clhs253*(clhs292 + clhs293) - clhs258;
lhs(5,1)=N[1]*clhs75 - clhs252*(clhs272 + clhs273) - clhs253*(clhs294 + clhs295) + clhs259;
lhs(5,2)=-DN(0,0)*clhs310 - DN(0,1)*clhs311 - clhs260;
lhs(5,3)=N[1]*clhs100 - clhs252*(clhs276 + clhs277) - clhs253*(clhs298 + clhs299);
lhs(5,4)=N[1]*clhs120 - clhs252*(clhs279 + clhs280) - clhs253*(clhs301 + clhs302);
lhs(5,5)=-DN(1,0)*clhs310 - DN(1,1)*clhs311 - pow(N[1], 2);
lhs(5,6)=N[1]*clhs136 - clhs252*(clhs282 + clhs283) - clhs253*(clhs304 + clhs305) + clhs313;
lhs(5,7)=N[1]*clhs156 - clhs252*(clhs285 + clhs286) - clhs253*(clhs307 + clhs308) - clhs314;
lhs(5,8)=-DN(2,0)*clhs310 - DN(2,1)*clhs311 - clhs315;
lhs(6,0)=DN(2,0)*clhs72 + DN(2,1)*clhs74 - clhs271*clhs316 - clhs271*clhs318 + clhs317*clhs38 + clhs319*clhs38 + clhs320*clhs54 + clhs321*clhs54;
lhs(6,1)=clhs36*(DN(2,0)*clhs89 + DN(2,1)*clhs90 - clhs274*clhs316 - clhs274*clhs318 + clhs317*clhs75 + clhs319*clhs75 + clhs322*clhs53 + clhs323*clhs53);
lhs(6,2)=N[0]*clhs324;
lhs(6,3)=DN(2,0)*clhs117 + DN(2,1)*clhs119 + clhs101*clhs317 + clhs101*clhs319 - clhs278*clhs316 - clhs278*clhs318 + clhs325*clhs54 + clhs326*clhs54;
lhs(6,4)=clhs36*(DN(2,0)*clhs134 + DN(2,1)*clhs135 + clhs120*clhs317 + clhs120*clhs319 - clhs281*clhs316 - clhs281*clhs318 + clhs327*clhs53 + clhs328*clhs53);
lhs(6,5)=N[1]*clhs324;
lhs(6,6)=DN(2,0)*clhs153 + DN(2,1)*clhs155 + clhs137*clhs317 + clhs137*clhs319 - clhs284*clhs316 - clhs284*clhs318 + clhs329*clhs54 + clhs330*clhs54;
lhs(6,7)=clhs36*(DN(2,0)*clhs170 + DN(2,1)*clhs171 + clhs156*clhs317 + clhs156*clhs319 - clhs287*clhs316 - clhs287*clhs318 + clhs331*clhs53 + clhs332*clhs53);
lhs(6,8)=N[2]*clhs324;
lhs(7,0)=clhs36*(DN(2,0)*clhs186 + DN(2,1)*clhs187 - clhs181*clhs333 - clhs181*clhs335 + clhs334*clhs37 + clhs336*clhs37 + clhs337*clhs53 + clhs338*clhs53);
lhs(7,1)=DN(2,0)*clhs203 + DN(2,1)*clhs204 + clhs188*clhs333 + clhs188*clhs335 - clhs296*clhs333 - clhs296*clhs335 + clhs339*clhs54 + clhs340*clhs54;
lhs(7,2)=N[0]*clhs341;
lhs(7,3)=clhs36*(DN(2,0)*clhs216 + DN(2,1)*clhs217 + clhs100*clhs334 + clhs100*clhs336 - clhs300*clhs333 - clhs300*clhs335 + clhs342*clhs53 + clhs343*clhs53);
lhs(7,4)=DN(2,0)*clhs228 + DN(2,1)*clhs229 + clhs128*clhs334 + clhs128*clhs336 - clhs303*clhs333 - clhs303*clhs335 + clhs344*clhs54 + clhs345*clhs54;
lhs(7,5)=N[1]*clhs341;
lhs(7,6)=clhs36*(DN(2,0)*clhs238 + DN(2,1)*clhs239 + clhs136*clhs334 + clhs136*clhs336 - clhs306*clhs333 - clhs306*clhs335 + clhs346*clhs53 + clhs347*clhs53);
lhs(7,7)=DN(2,0)*clhs250 + DN(2,1)*clhs251 + clhs164*clhs334 + clhs164*clhs336 - clhs309*clhs333 - clhs309*clhs335 + clhs348*clhs54 + clhs349*clhs54;
lhs(7,8)=N[2]*clhs341;
lhs(8,0)=N[2]*clhs37 - clhs252*(clhs320 + clhs321) - clhs253*(clhs337 + clhs338) - clhs262;
lhs(8,1)=N[2]*clhs75 - clhs252*(clhs322 + clhs323) - clhs253*(clhs339 + clhs340) + clhs263;
lhs(8,2)=-DN(0,0)*clhs350 - DN(0,1)*clhs351 - clhs264;
lhs(8,3)=N[2]*clhs100 - clhs252*(clhs325 + clhs326) - clhs253*(clhs342 + clhs343) - clhs313;
lhs(8,4)=N[2]*clhs120 - clhs252*(clhs327 + clhs328) - clhs253*(clhs344 + clhs345) + clhs314;
lhs(8,5)=-DN(1,0)*clhs350 - DN(1,1)*clhs351 - clhs315;
lhs(8,6)=N[2]*clhs136 - clhs252*(clhs329 + clhs330) - clhs253*(clhs346 + clhs347);
lhs(8,7)=N[2]*clhs156 - clhs252*(clhs331 + clhs332) - clhs253*(clhs348 + clhs349);
lhs(8,8)=-DN(2,0)*clhs350 - DN(2,1)*clhs351 - pow(N[2], 2);

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
const double crhs25 = -crhs0*crhs12 - crhs0*crhs13 - crhs1*crhs11 - crhs1*crhs13 - crhs11*crhs2 - crhs12*crhs2 + crhs16*crhs5 + crhs16*crhs6 + crhs17*crhs4 + crhs17*crhs6 + crhs18*crhs4 + crhs18*crhs5 + crhs19;
const double crhs26 = -N[0]*th[0] - N[1]*th[1] - N[2]*th[2] + crhs25 + crhs7;
const double crhs27 = 0.5*crhs26*tau_th/(crhs25 + crhs8);
const double crhs28 = C(0,2)*crhs15 + C(1,2)*crhs21 + C(2,2)*crhs22;
const double crhs29 = DN(0,1)*crhs28;
const double crhs30 = S[0]*crhs14 + S[2]*crhs20;
const double crhs31 = S[1]*crhs20 + S[2]*crhs14;
const double crhs32 = DN(0,0)*crhs28;
const double crhs33 = C(0,1)*crhs15 + C(1,1)*crhs21 + C(1,2)*crhs22;
const double crhs34 = DN(0,1)*crhs33;
const double crhs35 = b_gauss[0]*tau_u;
const double crhs36 = b_gauss[1]*tau_u;
const double crhs37 = (1.0/2.0)*tau_u;
const double crhs38 = crhs37*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs39 = crhs37*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double crhs40 = DN(1,0)*crhs23;
const double crhs41 = DN(1,1)*crhs28;
const double crhs42 = DN(1,0)*crhs28;
const double crhs43 = DN(1,1)*crhs33;
const double crhs44 = DN(2,0)*crhs23;
const double crhs45 = DN(2,1)*crhs28;
const double crhs46 = DN(2,0)*crhs28;
const double crhs47 = DN(2,1)*crhs33;
rhs[0]=-DN(0,0)*crhs9 - DN(0,1)*crhs10 + N[0]*b_gauss[0] - crhs24*crhs27 - crhs27*crhs29;
rhs[1]=-DN(0,0)*crhs30 - DN(0,1)*crhs31 + N[0]*b_gauss[1] - crhs27*crhs32 - crhs27*crhs34;
rhs[2]=-N[0]*crhs26 + crhs35*(DN(0,0)*crhs20 - DN(0,1)*crhs14) - crhs36*(DN(0,0)*crhs3 - DN(0,1)*crhs8) + crhs38*(crhs24 + crhs29) + crhs39*(crhs32 + crhs34);
rhs[3]=-DN(1,0)*crhs9 - DN(1,1)*crhs10 + N[1]*b_gauss[0] - crhs27*crhs40 - crhs27*crhs41;
rhs[4]=-DN(1,0)*crhs30 - DN(1,1)*crhs31 + N[1]*b_gauss[1] - crhs27*crhs42 - crhs27*crhs43;
rhs[5]=-N[1]*crhs26 + crhs35*(DN(1,0)*crhs20 - DN(1,1)*crhs14) - crhs36*(DN(1,0)*crhs3 - DN(1,1)*crhs8) + crhs38*(crhs40 + crhs41) + crhs39*(crhs42 + crhs43);
rhs[6]=-DN(2,0)*crhs9 - DN(2,1)*crhs10 + N[2]*b_gauss[0] - crhs27*crhs44 - crhs27*crhs45;
rhs[7]=-DN(2,0)*crhs30 - DN(2,1)*crhs31 + N[2]*b_gauss[1] - crhs27*crhs46 - crhs27*crhs47;
rhs[8]=-N[2]*crhs26 + crhs35*(DN(2,0)*crhs20 - DN(2,1)*crhs14) - crhs36*(DN(2,0)*crhs3 - DN(2,1)*crhs8) + crhs38*(crhs44 + crhs45) + crhs39*(crhs46 + crhs47);

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
const double cr_eq_green_strain15 = 1.0/(-cr_eq_green_strain0*cr_eq_green_strain10 - cr_eq_green_strain0*cr_eq_green_strain11 - cr_eq_green_strain1*cr_eq_green_strain11 - cr_eq_green_strain1*cr_eq_green_strain12 - cr_eq_green_strain10*cr_eq_green_strain2 - cr_eq_green_strain12*cr_eq_green_strain2 + cr_eq_green_strain13 + cr_eq_green_strain14 + cr_eq_green_strain4*cr_eq_green_strain5 + cr_eq_green_strain4*cr_eq_green_strain6 + cr_eq_green_strain5*cr_eq_green_strain9 + cr_eq_green_strain6*cr_eq_green_strain8 + cr_eq_green_strain7*cr_eq_green_strain8 + cr_eq_green_strain7*cr_eq_green_strain9);
const double cr_eq_green_strain16 = N[0]*th[0];
const double cr_eq_green_strain17 = N[1]*th[1];
const double cr_eq_green_strain18 = N[2]*th[2];
const double cr_eq_green_strain19 = cr_eq_green_strain15*(cr_eq_green_strain16 + cr_eq_green_strain17 + cr_eq_green_strain18 + 1.0);
const double cr_eq_green_strain20 = cr_eq_green_strain10 + cr_eq_green_strain11 + cr_eq_green_strain12;
const double cr_eq_green_strain21 = cr_eq_green_strain14 + 1;
r_eq_green_strain[0]=0.5*pow(cr_eq_green_strain13, 2)*cr_eq_green_strain19 + 0.5*cr_eq_green_strain19*pow(cr_eq_green_strain3, 2) - 0.5;
r_eq_green_strain[1]=0.5*cr_eq_green_strain19*pow(cr_eq_green_strain20, 2) + 0.5*cr_eq_green_strain19*pow(cr_eq_green_strain21, 2) - 0.5;
r_eq_green_strain[2]=cr_eq_green_strain15*(cr_eq_green_strain13*cr_eq_green_strain20 + cr_eq_green_strain21*cr_eq_green_strain3)*(1.0*cr_eq_green_strain16 + 1.0*cr_eq_green_strain17 + 1.0*cr_eq_green_strain18 + 1.0);

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
const double cr_eq_def_gradient14 = sqrt(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0)/sqrt(cr_eq_def_gradient0*cr_eq_def_gradient4 + cr_eq_def_gradient0*cr_eq_def_gradient5 + cr_eq_def_gradient1*cr_eq_def_gradient5 + cr_eq_def_gradient1*cr_eq_def_gradient6 - cr_eq_def_gradient10*cr_eq_def_gradient11 - cr_eq_def_gradient10*cr_eq_def_gradient12 - cr_eq_def_gradient11*cr_eq_def_gradient9 - cr_eq_def_gradient12*cr_eq_def_gradient7 + cr_eq_def_gradient13 + cr_eq_def_gradient2*cr_eq_def_gradient4 + cr_eq_def_gradient2*cr_eq_def_gradient6 + cr_eq_def_gradient3 - cr_eq_def_gradient7*cr_eq_def_gradient8 - cr_eq_def_gradient8*cr_eq_def_gradient9);
r_eq_def_gradient(0,0)=cr_eq_def_gradient14*cr_eq_def_gradient3;
r_eq_def_gradient(0,1)=cr_eq_def_gradient14*(cr_eq_def_gradient10 + cr_eq_def_gradient7 + cr_eq_def_gradient9);
r_eq_def_gradient(1,0)=cr_eq_def_gradient14*(cr_eq_def_gradient11 + cr_eq_def_gradient12 + cr_eq_def_gradient8);
r_eq_def_gradient(1,1)=cr_eq_def_gradient14*(cr_eq_def_gradient13 + 1);

    r_det_eq_def_gradient=1.0*(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1);

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
