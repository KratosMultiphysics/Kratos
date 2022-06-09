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
    const IndexType det_F_pos = r_geometry[0].GetDofPosition(DETERMINANT_F);

    IndexType aux_index = 0;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DETERMINANT_F, det_F_pos).EquationId();
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
    const IndexType det_F_pos = r_geometry[0].GetDofPosition(DETERMINANT_F);

    IndexType aux_index = 0;
    for (IndexType i_node = 0; i_node < NumNodes; ++i_node) {
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_X, disp_pos).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Y, disp_pos + 1).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DISPLACEMENT_Z, disp_pos + 2).EquationId();
        rResult[aux_index++] = r_geometry[i_node].GetDof(DETERMINANT_F, det_F_pos).EquationId();
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
        rElementalDofList[i * BlockSize + 2] = r_geometry[i].pGetDof(DETERMINANT_F);
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
        rElementalDofList[i * BlockSize + 3] = r_geometry[i].pGetDof(DETERMINANT_F);
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
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
    const SizeType dim = r_geometry.WorkingSpaceDimension();

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
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau * std::pow(h,2) / 2.0;

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
        // TODO: THIS MUST BE COMPUTED BY THE CONSTITUTIVE LAW, ASSUMING LINEAR ELASTIC MATERIAL SO FAR
        const auto& r_prop = GetProperties();
        double mu = r_prop[YOUNG_MODULUS]/(2*(1.0+r_prop[POISSON_RATIO])); // 2nd Lame constant (shear modulus)
        const double tau = aux_tau / mu;

        // Calculate and add the LHS Gauss point contributions
        const double clhs0 = DN(0,1)*u(0,0);
const double clhs1 = DN(1,1)*u(1,0);
const double clhs2 = DN(2,1)*u(2,0);
const double clhs3 = clhs0 + clhs1 + clhs2;
const double clhs4 = DN(0,1)*clhs3;
const double clhs5 = DN(0,0)*u(0,0);
const double clhs6 = DN(1,1)*u(1,1);
const double clhs7 = DN(2,1)*u(2,1);
const double clhs8 = DN(0,1)*u(0,1);
const double clhs9 = DN(1,0)*u(1,0);
const double clhs10 = DN(2,0)*u(2,0);
const double clhs11 = DN(0,0)*u(0,1);
const double clhs12 = DN(1,0)*u(1,1);
const double clhs13 = DN(2,0)*u(2,1);
const double clhs14 = clhs10 + clhs5 + clhs9 + 1;
const double clhs15 = clhs6 + clhs7 + clhs8;
const double clhs16 = -clhs0*clhs12 - clhs0*clhs13 - clhs1*clhs11 - clhs1*clhs13 + clhs10*clhs6 + clhs10*clhs8 - clhs11*clhs2 - clhs12*clhs2 + clhs14 + clhs15 + clhs5*clhs6 + clhs5*clhs7 + clhs7*clhs9 + clhs8*clhs9;
const double clhs17 = 1.0/clhs16;
const double clhs18 = 1.0*clhs17;
const double clhs19 = pow(clhs3, 2);
const double clhs20 = pow(clhs16, -2.0);
const double clhs21 = DN(0,1)*clhs12;
const double clhs22 = DN(0,1)*clhs13;
const double clhs23 = DN(0,0)*clhs6;
const double clhs24 = DN(0,0)*clhs7;
const double clhs25 = clhs20*(-DN(0,0) + clhs21 + clhs22 - clhs23 - clhs24);
const double clhs26 = 0.5*clhs25;
const double clhs27 = clhs15 + 1;
const double clhs28 = pow(clhs27, 2);
const double clhs29 = clhs18*clhs4 + clhs19*clhs26 + clhs26*clhs28;
const double clhs30 = DN(0,0)*clhs14;
const double clhs31 = clhs11 + clhs12 + clhs13;
const double clhs32 = pow(clhs31, 2);
const double clhs33 = pow(clhs14, 2);
const double clhs34 = clhs18*clhs30 + clhs26*clhs32 + clhs26*clhs33;
const double clhs35 = DN(0,0)*clhs3;
const double clhs36 = DN(0,1)*clhs14;
const double clhs37 = clhs27*clhs31;
const double clhs38 = clhs14*clhs3;
const double clhs39 = clhs17*clhs35 + clhs17*clhs36 + clhs25*clhs37 + clhs25*clhs38;
const double clhs40 = C(0,0)*clhs34 + C(0,1)*clhs29 + C(0,2)*clhs39;
const double clhs41 = N[0]*th[0] + N[1]*th[1] + N[2]*th[2];
const double clhs42 = pow(clhs41, 2);
const double clhs43 = clhs14*clhs42;
const double clhs44 = C(0,2)*clhs34 + C(1,2)*clhs29 + C(2,2)*clhs39;
const double clhs45 = clhs3*clhs42;
const double clhs46 = clhs37 + clhs38;
const double clhs47 = clhs17*clhs42;
const double clhs48 = clhs46*clhs47;
const double clhs49 = 0.5*clhs32*clhs47 + 0.5*clhs33*clhs47 - 0.5;
const double clhs50 = 0.5*clhs19*clhs47 + 0.5*clhs28*clhs47 - 0.5;
const double clhs51 = C(0,0)*clhs49 + C(0,1)*clhs50 + C(0,2)*clhs48;
const double clhs52 = C(0,2)*clhs49 + C(1,2)*clhs50 + C(2,2)*clhs48;
const double clhs53 = DN(0,0)*clhs51 + DN(0,1)*clhs52;
const double clhs54 = clhs40*clhs43 + clhs44*clhs45 + clhs53;
const double clhs55 = C(0,1)*clhs34 + C(1,1)*clhs29 + C(1,2)*clhs39;
const double clhs56 = C(0,1)*clhs49 + C(1,1)*clhs50 + C(1,2)*clhs48;
const double clhs57 = DN(0,0)*clhs52 + DN(0,1)*clhs56;
const double clhs58 = clhs43*clhs44 + clhs45*clhs55 + clhs57;
const double clhs59 = DN(0,0)*clhs31;
const double clhs60 = DN(0,0)*clhs1;
const double clhs61 = DN(0,0)*clhs2;
const double clhs62 = DN(0,1)*clhs9;
const double clhs63 = DN(0,1)*clhs10;
const double clhs64 = clhs20*(-DN(0,1) + clhs60 + clhs61 - clhs62 - clhs63);
const double clhs65 = 0.5*clhs64;
const double clhs66 = clhs18*clhs59 + clhs32*clhs65 + clhs33*clhs65;
const double clhs67 = DN(0,1)*clhs27;
const double clhs68 = clhs18*clhs67 + clhs19*clhs65 + clhs28*clhs65;
const double clhs69 = DN(0,1)*clhs31;
const double clhs70 = DN(0,0)*clhs27;
const double clhs71 = clhs17*clhs69 + clhs17*clhs70 + clhs37*clhs64 + clhs38*clhs64;
const double clhs72 = C(0,2)*clhs66 + C(1,2)*clhs68 + C(2,2)*clhs71;
const double clhs73 = C(0,0)*clhs66 + C(0,1)*clhs68 + C(0,2)*clhs71;
const double clhs74 = clhs14*clhs73 + clhs3*clhs72;
const double clhs75 = C(0,1)*clhs66 + C(1,1)*clhs68 + C(1,2)*clhs71;
const double clhs76 = clhs14*clhs72 + clhs3*clhs75;
const double clhs77 = clhs32 + clhs33;
const double clhs78 = C(0,2)*clhs77;
const double clhs79 = clhs19 + clhs28;
const double clhs80 = C(1,2)*clhs79;
const double clhs81 = 2.0*clhs46;
const double clhs82 = C(2,2)*clhs81 + 1.0*clhs78 + 1.0*clhs80;
const double clhs83 = C(0,0)*clhs77;
const double clhs84 = C(0,1)*clhs79;
const double clhs85 = C(0,2)*clhs81 + 1.0*clhs83 + 1.0*clhs84;
const double clhs86 = clhs14*clhs85 + clhs3*clhs82;
const double clhs87 = C(0,1)*clhs77;
const double clhs88 = C(1,1)*clhs79;
const double clhs89 = C(1,2)*clhs81 + 1.0*clhs87 + 1.0*clhs88;
const double clhs90 = clhs14*clhs82 + clhs3*clhs89;
const double clhs91 = clhs17*clhs41;
const double clhs92 = clhs91*(DN(0,0)*clhs86 + DN(0,1)*clhs90);
const double clhs93 = DN(1,1)*clhs3;
const double clhs94 = DN(1,1)*clhs11;
const double clhs95 = DN(1,1)*clhs13;
const double clhs96 = DN(1,0)*clhs8;
const double clhs97 = DN(1,0)*clhs7;
const double clhs98 = clhs20*(-DN(1,0) + clhs94 + clhs95 - clhs96 - clhs97);
const double clhs99 = 0.5*clhs98;
const double clhs100 = clhs18*clhs93 + clhs19*clhs99 + clhs28*clhs99;
const double clhs101 = DN(1,0)*clhs14;
const double clhs102 = clhs101*clhs18 + clhs32*clhs99 + clhs33*clhs99;
const double clhs103 = DN(1,0)*clhs3;
const double clhs104 = DN(1,1)*clhs14;
const double clhs105 = clhs103*clhs17 + clhs104*clhs17 + clhs37*clhs98 + clhs38*clhs98;
const double clhs106 = C(0,0)*clhs102 + C(0,1)*clhs100 + C(0,2)*clhs105;
const double clhs107 = C(0,2)*clhs102 + C(1,2)*clhs100 + C(2,2)*clhs105;
const double clhs108 = DN(1,0)*clhs51 + DN(1,1)*clhs52;
const double clhs109 = clhs106*clhs43 + clhs107*clhs45 + clhs108;
const double clhs110 = C(0,1)*clhs102 + C(1,1)*clhs100 + C(1,2)*clhs105;
const double clhs111 = DN(1,0)*clhs52 + DN(1,1)*clhs56;
const double clhs112 = clhs107*clhs43 + clhs110*clhs45 + clhs111;
const double clhs113 = DN(1,0)*clhs31;
const double clhs114 = DN(1,0)*clhs0;
const double clhs115 = DN(1,0)*clhs2;
const double clhs116 = DN(1,1)*clhs5;
const double clhs117 = DN(1,1)*clhs10;
const double clhs118 = clhs20*(-DN(1,1) + clhs114 + clhs115 - clhs116 - clhs117);
const double clhs119 = 0.5*clhs118;
const double clhs120 = clhs113*clhs18 + clhs119*clhs32 + clhs119*clhs33;
const double clhs121 = DN(1,1)*clhs27;
const double clhs122 = clhs119*clhs19 + clhs119*clhs28 + clhs121*clhs18;
const double clhs123 = DN(1,1)*clhs31;
const double clhs124 = DN(1,0)*clhs27;
const double clhs125 = clhs118*clhs37 + clhs118*clhs38 + clhs123*clhs17 + clhs124*clhs17;
const double clhs126 = C(0,2)*clhs120 + C(1,2)*clhs122 + C(2,2)*clhs125;
const double clhs127 = C(0,0)*clhs120 + C(0,1)*clhs122 + C(0,2)*clhs125;
const double clhs128 = clhs126*clhs3 + clhs127*clhs14;
const double clhs129 = C(0,1)*clhs120 + C(1,1)*clhs122 + C(1,2)*clhs125;
const double clhs130 = clhs126*clhs14 + clhs129*clhs3;
const double clhs131 = DN(2,1)*clhs3;
const double clhs132 = DN(2,1)*clhs11;
const double clhs133 = DN(2,1)*clhs12;
const double clhs134 = DN(2,0)*clhs8;
const double clhs135 = DN(2,0)*clhs6;
const double clhs136 = clhs20*(-DN(2,0) + clhs132 + clhs133 - clhs134 - clhs135);
const double clhs137 = 0.5*clhs136;
const double clhs138 = clhs131*clhs18 + clhs137*clhs19 + clhs137*clhs28;
const double clhs139 = DN(2,0)*clhs14;
const double clhs140 = clhs137*clhs32 + clhs137*clhs33 + clhs139*clhs18;
const double clhs141 = DN(2,0)*clhs3;
const double clhs142 = DN(2,1)*clhs14;
const double clhs143 = clhs136*clhs37 + clhs136*clhs38 + clhs141*clhs17 + clhs142*clhs17;
const double clhs144 = C(0,0)*clhs140 + C(0,1)*clhs138 + C(0,2)*clhs143;
const double clhs145 = C(0,2)*clhs140 + C(1,2)*clhs138 + C(2,2)*clhs143;
const double clhs146 = DN(2,0)*clhs51 + DN(2,1)*clhs52;
const double clhs147 = clhs144*clhs43 + clhs145*clhs45 + clhs146;
const double clhs148 = C(0,1)*clhs140 + C(1,1)*clhs138 + C(1,2)*clhs143;
const double clhs149 = DN(2,0)*clhs52 + DN(2,1)*clhs56;
const double clhs150 = clhs145*clhs43 + clhs148*clhs45 + clhs149;
const double clhs151 = DN(2,0)*clhs31;
const double clhs152 = DN(2,0)*clhs0;
const double clhs153 = DN(2,0)*clhs1;
const double clhs154 = DN(2,1)*clhs5;
const double clhs155 = DN(2,1)*clhs9;
const double clhs156 = clhs20*(-DN(2,1) + clhs152 + clhs153 - clhs154 - clhs155);
const double clhs157 = 0.5*clhs156;
const double clhs158 = clhs151*clhs18 + clhs157*clhs32 + clhs157*clhs33;
const double clhs159 = DN(2,1)*clhs27;
const double clhs160 = clhs157*clhs19 + clhs157*clhs28 + clhs159*clhs18;
const double clhs161 = DN(2,1)*clhs31;
const double clhs162 = DN(2,0)*clhs27;
const double clhs163 = clhs156*clhs37 + clhs156*clhs38 + clhs161*clhs17 + clhs162*clhs17;
const double clhs164 = C(0,2)*clhs158 + C(1,2)*clhs160 + C(2,2)*clhs163;
const double clhs165 = C(0,0)*clhs158 + C(0,1)*clhs160 + C(0,2)*clhs163;
const double clhs166 = clhs14*clhs165 + clhs164*clhs3;
const double clhs167 = C(0,1)*clhs158 + C(1,1)*clhs160 + C(1,2)*clhs163;
const double clhs168 = clhs14*clhs164 + clhs167*clhs3;
const double clhs169 = clhs27*clhs44 + clhs31*clhs40;
const double clhs170 = clhs27*clhs55 + clhs31*clhs44;
const double clhs171 = clhs31*clhs42;
const double clhs172 = clhs27*clhs42;
const double clhs173 = clhs171*clhs73 + clhs172*clhs72 + clhs53;
const double clhs174 = clhs171*clhs72 + clhs172*clhs75 + clhs57;
const double clhs175 = clhs27*clhs82 + clhs31*clhs85;
const double clhs176 = clhs27*clhs89 + clhs31*clhs82;
const double clhs177 = clhs91*(DN(0,0)*clhs175 + DN(0,1)*clhs176);
const double clhs178 = clhs106*clhs31 + clhs107*clhs27;
const double clhs179 = clhs107*clhs31 + clhs110*clhs27;
const double clhs180 = clhs108 + clhs126*clhs172 + clhs127*clhs171;
const double clhs181 = clhs111 + clhs126*clhs171 + clhs129*clhs172;
const double clhs182 = clhs144*clhs31 + clhs145*clhs27;
const double clhs183 = clhs145*clhs31 + clhs148*clhs27;
const double clhs184 = clhs146 + clhs164*clhs172 + clhs165*clhs171;
const double clhs185 = clhs149 + clhs164*clhs171 + clhs167*clhs172;
const double clhs186 = DN(0,0) - clhs21 - clhs22 + clhs23 + clhs24;
const double clhs187 = clhs35 + clhs36;
const double clhs188 = C(0,0)*clhs30 + C(0,1)*clhs4 + C(0,2)*clhs187;
const double clhs189 = C(0,2)*clhs30 + C(1,2)*clhs4 + C(2,2)*clhs187;
const double clhs190 = 2*clhs46;
const double clhs191 = C(0,2)*clhs190 + clhs83 + clhs84;
const double clhs192 = DN(0,0)*clhs191;
const double clhs193 = 0.25/clhs16;
const double clhs194 = clhs186*clhs193;
const double clhs195 = C(2,2)*clhs190 + clhs78 + clhs80;
const double clhs196 = DN(0,1)*clhs195;
const double clhs197 = 1.0/clhs41;
const double clhs198 = tau*sqrt(clhs16*clhs197);
const double clhs199 = clhs198*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double clhs200 = C(0,1)*clhs30 + C(1,1)*clhs4 + C(1,2)*clhs187;
const double clhs201 = DN(0,0)*clhs195;
const double clhs202 = C(1,2)*clhs190 + clhs87 + clhs88;
const double clhs203 = DN(0,1)*clhs202;
const double clhs204 = clhs198*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double clhs205 = DN(0,1) - clhs60 - clhs61 + clhs62 + clhs63;
const double clhs206 = clhs69 + clhs70;
const double clhs207 = C(0,0)*clhs59 + C(0,1)*clhs67 + C(0,2)*clhs206;
const double clhs208 = C(0,2)*clhs59 + C(1,2)*clhs67 + C(2,2)*clhs206;
const double clhs209 = clhs193*clhs205;
const double clhs210 = C(0,1)*clhs59 + C(1,1)*clhs67 + C(1,2)*clhs206;
const double clhs211 = clhs192 + clhs196;
const double clhs212 = (1.0/2.0)*clhs198;
const double clhs213 = clhs211*clhs212;
const double clhs214 = clhs201 + clhs203;
const double clhs215 = clhs212*clhs214;
const double clhs216 = 0.25*clhs197;
const double clhs217 = N[0]*clhs216;
const double clhs218 = clhs199*clhs211;
const double clhs219 = clhs204*clhs214;
const double clhs220 = tau*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs221 = b_gauss[1]*clhs220;
const double clhs222 = DN(1,0) - clhs94 - clhs95 + clhs96 + clhs97;
const double clhs223 = clhs103 + clhs104;
const double clhs224 = C(0,0)*clhs101 + C(0,1)*clhs93 + C(0,2)*clhs223;
const double clhs225 = C(0,2)*clhs101 + C(1,2)*clhs93 + C(2,2)*clhs223;
const double clhs226 = clhs193*clhs222;
const double clhs227 = C(0,1)*clhs101 + C(1,1)*clhs93 + C(1,2)*clhs223;
const double clhs228 = b_gauss[0]*clhs220;
const double clhs229 = DN(1,1) - clhs114 - clhs115 + clhs116 + clhs117;
const double clhs230 = clhs123 + clhs124;
const double clhs231 = C(0,0)*clhs113 + C(0,1)*clhs121 + C(0,2)*clhs230;
const double clhs232 = C(0,2)*clhs113 + C(1,2)*clhs121 + C(2,2)*clhs230;
const double clhs233 = clhs193*clhs229;
const double clhs234 = C(0,1)*clhs113 + C(1,1)*clhs121 + C(1,2)*clhs230;
const double clhs235 = -N[0]*N[1];
const double clhs236 = N[1]*clhs216;
const double clhs237 = tau*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs238 = b_gauss[1]*clhs237;
const double clhs239 = DN(2,0) - clhs132 - clhs133 + clhs134 + clhs135;
const double clhs240 = clhs141 + clhs142;
const double clhs241 = C(0,0)*clhs139 + C(0,1)*clhs131 + C(0,2)*clhs240;
const double clhs242 = C(0,2)*clhs139 + C(1,2)*clhs131 + C(2,2)*clhs240;
const double clhs243 = clhs193*clhs239;
const double clhs244 = C(0,1)*clhs139 + C(1,1)*clhs131 + C(1,2)*clhs240;
const double clhs245 = b_gauss[0]*clhs237;
const double clhs246 = DN(2,1) - clhs152 - clhs153 + clhs154 + clhs155;
const double clhs247 = clhs161 + clhs162;
const double clhs248 = C(0,0)*clhs151 + C(0,1)*clhs159 + C(0,2)*clhs247;
const double clhs249 = C(0,2)*clhs151 + C(1,2)*clhs159 + C(2,2)*clhs247;
const double clhs250 = clhs193*clhs246;
const double clhs251 = C(0,1)*clhs151 + C(1,1)*clhs159 + C(1,2)*clhs247;
const double clhs252 = -N[0]*N[2];
const double clhs253 = N[2]*clhs216;
const double clhs254 = clhs91*(DN(1,0)*clhs86 + DN(1,1)*clhs90);
const double clhs255 = clhs91*(DN(1,0)*clhs175 + DN(1,1)*clhs176);
const double clhs256 = DN(1,0)*clhs191;
const double clhs257 = DN(1,1)*clhs195;
const double clhs258 = DN(1,0)*clhs195;
const double clhs259 = DN(1,1)*clhs202;
const double clhs260 = clhs256 + clhs257;
const double clhs261 = clhs212*clhs260;
const double clhs262 = clhs258 + clhs259;
const double clhs263 = clhs212*clhs262;
const double clhs264 = clhs199*clhs260;
const double clhs265 = clhs204*clhs262;
const double clhs266 = tau*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs267 = b_gauss[1]*clhs266;
const double clhs268 = b_gauss[0]*clhs266;
const double clhs269 = -N[1]*N[2];
const double clhs270 = clhs91*(DN(2,0)*clhs86 + DN(2,1)*clhs90);
const double clhs271 = clhs91*(DN(2,0)*clhs175 + DN(2,1)*clhs176);
const double clhs272 = DN(2,0)*clhs191;
const double clhs273 = DN(2,1)*clhs195;
const double clhs274 = DN(2,0)*clhs195;
const double clhs275 = DN(2,1)*clhs202;
const double clhs276 = clhs272 + clhs273;
const double clhs277 = clhs212*clhs276;
const double clhs278 = clhs274 + clhs275;
const double clhs279 = clhs212*clhs278;
const double clhs280 = clhs199*clhs276;
const double clhs281 = clhs204*clhs278;
lhs(0,0)=DN(0,0)*clhs54 + DN(0,1)*clhs58;
lhs(0,1)=clhs42*(DN(0,0)*clhs74 + DN(0,1)*clhs76);
lhs(0,2)=N[0]*clhs92;
lhs(0,3)=DN(0,0)*clhs109 + DN(0,1)*clhs112;
lhs(0,4)=clhs42*(DN(0,0)*clhs128 + DN(0,1)*clhs130);
lhs(0,5)=N[1]*clhs92;
lhs(0,6)=DN(0,0)*clhs147 + DN(0,1)*clhs150;
lhs(0,7)=clhs42*(DN(0,0)*clhs166 + DN(0,1)*clhs168);
lhs(0,8)=N[2]*clhs92;
lhs(1,0)=clhs42*(DN(0,0)*clhs169 + DN(0,1)*clhs170);
lhs(1,1)=DN(0,0)*clhs173 + DN(0,1)*clhs174;
lhs(1,2)=N[0]*clhs177;
lhs(1,3)=clhs42*(DN(0,0)*clhs178 + DN(0,1)*clhs179);
lhs(1,4)=DN(0,0)*clhs180 + DN(0,1)*clhs181;
lhs(1,5)=N[1]*clhs177;
lhs(1,6)=clhs42*(DN(0,0)*clhs182 + DN(0,1)*clhs183);
lhs(1,7)=DN(0,0)*clhs184 + DN(0,1)*clhs185;
lhs(1,8)=N[2]*clhs177;
lhs(2,0)=N[0]*clhs186 - clhs199*(DN(0,0)*clhs188 + DN(0,1)*clhs189 + clhs192*clhs194 + clhs194*clhs196) - clhs204*(DN(0,0)*clhs189 + DN(0,1)*clhs200 + clhs194*clhs201 + clhs194*clhs203);
lhs(2,1)=N[0]*clhs205 - clhs199*(DN(0,0)*clhs207 + DN(0,1)*clhs208 + clhs192*clhs209 + clhs196*clhs209) - clhs204*(DN(0,0)*clhs208 + DN(0,1)*clhs210 + clhs201*clhs209 + clhs203*clhs209);
lhs(2,2)=-DN(0,0)*clhs213 - DN(0,1)*clhs215 - pow(N[0], 2) + clhs217*clhs218 + clhs217*clhs219;
lhs(2,3)=N[0]*clhs222 - clhs199*(DN(0,0)*clhs224 + DN(0,1)*clhs225 + clhs192*clhs226 + clhs196*clhs226) - clhs204*(DN(0,0)*clhs225 + DN(0,1)*clhs227 + clhs201*clhs226 + clhs203*clhs226) + clhs221;
lhs(2,4)=N[0]*clhs229 - clhs199*(DN(0,0)*clhs231 + DN(0,1)*clhs232 + clhs192*clhs233 + clhs196*clhs233) - clhs204*(DN(0,0)*clhs232 + DN(0,1)*clhs234 + clhs201*clhs233 + clhs203*clhs233) - clhs228;
lhs(2,5)=-DN(1,0)*clhs213 - DN(1,1)*clhs215 + clhs218*clhs236 + clhs219*clhs236 + clhs235;
lhs(2,6)=N[0]*clhs239 - clhs199*(DN(0,0)*clhs241 + DN(0,1)*clhs242 + clhs192*clhs243 + clhs196*clhs243) - clhs204*(DN(0,0)*clhs242 + DN(0,1)*clhs244 + clhs201*clhs243 + clhs203*clhs243) + clhs238;
lhs(2,7)=N[0]*clhs246 - clhs199*(DN(0,0)*clhs248 + DN(0,1)*clhs249 + clhs192*clhs250 + clhs196*clhs250) - clhs204*(DN(0,0)*clhs249 + DN(0,1)*clhs251 + clhs201*clhs250 + clhs203*clhs250) - clhs245;
lhs(2,8)=-DN(2,0)*clhs213 - DN(2,1)*clhs215 + clhs218*clhs253 + clhs219*clhs253 + clhs252;
lhs(3,0)=DN(1,0)*clhs54 + DN(1,1)*clhs58;
lhs(3,1)=clhs42*(DN(1,0)*clhs74 + DN(1,1)*clhs76);
lhs(3,2)=N[0]*clhs254;
lhs(3,3)=DN(1,0)*clhs109 + DN(1,1)*clhs112;
lhs(3,4)=clhs42*(DN(1,0)*clhs128 + DN(1,1)*clhs130);
lhs(3,5)=N[1]*clhs254;
lhs(3,6)=DN(1,0)*clhs147 + DN(1,1)*clhs150;
lhs(3,7)=clhs42*(DN(1,0)*clhs166 + DN(1,1)*clhs168);
lhs(3,8)=N[2]*clhs254;
lhs(4,0)=clhs42*(DN(1,0)*clhs169 + DN(1,1)*clhs170);
lhs(4,1)=DN(1,0)*clhs173 + DN(1,1)*clhs174;
lhs(4,2)=N[0]*clhs255;
lhs(4,3)=clhs42*(DN(1,0)*clhs178 + DN(1,1)*clhs179);
lhs(4,4)=DN(1,0)*clhs180 + DN(1,1)*clhs181;
lhs(4,5)=N[1]*clhs255;
lhs(4,6)=clhs42*(DN(1,0)*clhs182 + DN(1,1)*clhs183);
lhs(4,7)=DN(1,0)*clhs184 + DN(1,1)*clhs185;
lhs(4,8)=N[2]*clhs255;
lhs(5,0)=N[1]*clhs186 - clhs199*(DN(1,0)*clhs188 + DN(1,1)*clhs189 + clhs194*clhs256 + clhs194*clhs257) - clhs204*(DN(1,0)*clhs189 + DN(1,1)*clhs200 + clhs194*clhs258 + clhs194*clhs259) - clhs221;
lhs(5,1)=N[1]*clhs205 - clhs199*(DN(1,0)*clhs207 + DN(1,1)*clhs208 + clhs209*clhs256 + clhs209*clhs257) - clhs204*(DN(1,0)*clhs208 + DN(1,1)*clhs210 + clhs209*clhs258 + clhs209*clhs259) + clhs228;
lhs(5,2)=-DN(0,0)*clhs261 - DN(0,1)*clhs263 + clhs217*clhs264 + clhs217*clhs265 + clhs235;
lhs(5,3)=N[1]*clhs222 - clhs199*(DN(1,0)*clhs224 + DN(1,1)*clhs225 + clhs226*clhs256 + clhs226*clhs257) - clhs204*(DN(1,0)*clhs225 + DN(1,1)*clhs227 + clhs226*clhs258 + clhs226*clhs259);
lhs(5,4)=N[1]*clhs229 - clhs199*(DN(1,0)*clhs231 + DN(1,1)*clhs232 + clhs233*clhs256 + clhs233*clhs257) - clhs204*(DN(1,0)*clhs232 + DN(1,1)*clhs234 + clhs233*clhs258 + clhs233*clhs259);
lhs(5,5)=-DN(1,0)*clhs261 - DN(1,1)*clhs263 - pow(N[1], 2) + clhs236*clhs264 + clhs236*clhs265;
lhs(5,6)=N[1]*clhs239 - clhs199*(DN(1,0)*clhs241 + DN(1,1)*clhs242 + clhs243*clhs256 + clhs243*clhs257) - clhs204*(DN(1,0)*clhs242 + DN(1,1)*clhs244 + clhs243*clhs258 + clhs243*clhs259) + clhs267;
lhs(5,7)=N[1]*clhs246 - clhs199*(DN(1,0)*clhs248 + DN(1,1)*clhs249 + clhs250*clhs256 + clhs250*clhs257) - clhs204*(DN(1,0)*clhs249 + DN(1,1)*clhs251 + clhs250*clhs258 + clhs250*clhs259) - clhs268;
lhs(5,8)=-DN(2,0)*clhs261 - DN(2,1)*clhs263 + clhs253*clhs264 + clhs253*clhs265 + clhs269;
lhs(6,0)=DN(2,0)*clhs54 + DN(2,1)*clhs58;
lhs(6,1)=clhs42*(DN(2,0)*clhs74 + DN(2,1)*clhs76);
lhs(6,2)=N[0]*clhs270;
lhs(6,3)=DN(2,0)*clhs109 + DN(2,1)*clhs112;
lhs(6,4)=clhs42*(DN(2,0)*clhs128 + DN(2,1)*clhs130);
lhs(6,5)=N[1]*clhs270;
lhs(6,6)=DN(2,0)*clhs147 + DN(2,1)*clhs150;
lhs(6,7)=clhs42*(DN(2,0)*clhs166 + DN(2,1)*clhs168);
lhs(6,8)=N[2]*clhs270;
lhs(7,0)=clhs42*(DN(2,0)*clhs169 + DN(2,1)*clhs170);
lhs(7,1)=DN(2,0)*clhs173 + DN(2,1)*clhs174;
lhs(7,2)=N[0]*clhs271;
lhs(7,3)=clhs42*(DN(2,0)*clhs178 + DN(2,1)*clhs179);
lhs(7,4)=DN(2,0)*clhs180 + DN(2,1)*clhs181;
lhs(7,5)=N[1]*clhs271;
lhs(7,6)=clhs42*(DN(2,0)*clhs182 + DN(2,1)*clhs183);
lhs(7,7)=DN(2,0)*clhs184 + DN(2,1)*clhs185;
lhs(7,8)=N[2]*clhs271;
lhs(8,0)=N[2]*clhs186 - clhs199*(DN(2,0)*clhs188 + DN(2,1)*clhs189 + clhs194*clhs272 + clhs194*clhs273) - clhs204*(DN(2,0)*clhs189 + DN(2,1)*clhs200 + clhs194*clhs274 + clhs194*clhs275) - clhs238;
lhs(8,1)=N[2]*clhs205 - clhs199*(DN(2,0)*clhs207 + DN(2,1)*clhs208 + clhs209*clhs272 + clhs209*clhs273) - clhs204*(DN(2,0)*clhs208 + DN(2,1)*clhs210 + clhs209*clhs274 + clhs209*clhs275) + clhs245;
lhs(8,2)=-DN(0,0)*clhs277 - DN(0,1)*clhs279 + clhs217*clhs280 + clhs217*clhs281 + clhs252;
lhs(8,3)=N[2]*clhs222 - clhs199*(DN(2,0)*clhs224 + DN(2,1)*clhs225 + clhs226*clhs272 + clhs226*clhs273) - clhs204*(DN(2,0)*clhs225 + DN(2,1)*clhs227 + clhs226*clhs274 + clhs226*clhs275) - clhs267;
lhs(8,4)=N[2]*clhs229 - clhs199*(DN(2,0)*clhs231 + DN(2,1)*clhs232 + clhs233*clhs272 + clhs233*clhs273) - clhs204*(DN(2,0)*clhs232 + DN(2,1)*clhs234 + clhs233*clhs274 + clhs233*clhs275) + clhs268;
lhs(8,5)=-DN(1,0)*clhs277 - DN(1,1)*clhs279 + clhs236*clhs280 + clhs236*clhs281 + clhs269;
lhs(8,6)=N[2]*clhs239 - clhs199*(DN(2,0)*clhs241 + DN(2,1)*clhs242 + clhs243*clhs272 + clhs243*clhs273) - clhs204*(DN(2,0)*clhs242 + DN(2,1)*clhs244 + clhs243*clhs274 + clhs243*clhs275);
lhs(8,7)=N[2]*clhs246 - clhs199*(DN(2,0)*clhs248 + DN(2,1)*clhs249 + clhs250*clhs272 + clhs250*clhs273) - clhs204*(DN(2,0)*clhs249 + DN(2,1)*clhs251 + clhs250*clhs274 + clhs250*clhs275);
lhs(8,8)=-DN(2,0)*clhs277 - DN(2,1)*clhs279 - pow(N[2], 2) + clhs253*clhs280 + clhs253*clhs281;

        // Calculate and add the RHS Gauss point contribution
        const double crhs0 = DN(0,1)*u(0,0);
const double crhs1 = DN(1,1)*u(1,0);
const double crhs2 = DN(2,1)*u(2,0);
const double crhs3 = crhs0 + crhs1 + crhs2;
const double crhs4 = DN(0,0)*u(0,0);
const double crhs5 = DN(1,0)*u(1,0);
const double crhs6 = DN(2,0)*u(2,0);
const double crhs7 = crhs4 + crhs5 + crhs6 + 1;
const double crhs8 = S[0]*crhs7 + S[2]*crhs3;
const double crhs9 = S[1]*crhs3 + S[2]*crhs7;
const double crhs10 = DN(0,0)*u(0,1);
const double crhs11 = DN(1,0)*u(1,1);
const double crhs12 = DN(2,0)*u(2,1);
const double crhs13 = crhs10 + crhs11 + crhs12;
const double crhs14 = DN(0,1)*u(0,1);
const double crhs15 = DN(1,1)*u(1,1);
const double crhs16 = DN(2,1)*u(2,1);
const double crhs17 = crhs14 + crhs15 + crhs16;
const double crhs18 = crhs17 + 1;
const double crhs19 = S[0]*crhs13 + S[2]*crhs18;
const double crhs20 = S[1]*crhs18 + S[2]*crhs13;
const double crhs21 = b_gauss[0]*tau;
const double crhs22 = b_gauss[1]*tau;
const double crhs23 = N[0]*th[0];
const double crhs24 = N[1]*th[1];
const double crhs25 = N[2]*th[2];
const double crhs26 = -crhs0*crhs11 - crhs0*crhs12 - crhs1*crhs10 - crhs1*crhs12 - crhs10*crhs2 - crhs11*crhs2 + crhs14*crhs5 + crhs14*crhs6 + crhs15*crhs4 + crhs15*crhs6 + crhs16*crhs4 + crhs16*crhs5 + crhs17 + crhs7;
const double crhs27 = -crhs23 - crhs24 - crhs25 + crhs26;
const double crhs28 = pow(crhs13, 2) + pow(crhs7, 2);
const double crhs29 = pow(crhs18, 2) + pow(crhs3, 2);
const double crhs30 = 2*crhs13*crhs18 + 2*crhs3*crhs7;
const double crhs31 = C(0,0)*crhs28 + C(0,1)*crhs29 + C(0,2)*crhs30;
const double crhs32 = C(0,2)*crhs28 + C(1,2)*crhs29 + C(2,2)*crhs30;
const double crhs33 = (1.0/2.0)*tau*sqrt(crhs26/(crhs23 + crhs24 + crhs25));
const double crhs34 = crhs33*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs35 = C(0,1)*crhs28 + C(1,1)*crhs29 + C(1,2)*crhs30;
const double crhs36 = crhs33*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
rhs[0]=-DN(0,0)*crhs8 - DN(0,1)*crhs9 + N[0]*b_gauss[0];
rhs[1]=-DN(0,0)*crhs19 - DN(0,1)*crhs20 + N[0]*b_gauss[1];
rhs[2]=-N[0]*crhs27 + crhs21*(DN(0,0)*crhs18 - DN(0,1)*crhs13) - crhs22*(DN(0,0)*crhs3 - DN(0,1)*crhs7) + crhs34*(DN(0,0)*crhs31 + DN(0,1)*crhs32) + crhs36*(DN(0,0)*crhs32 + DN(0,1)*crhs35);
rhs[3]=-DN(1,0)*crhs8 - DN(1,1)*crhs9 + N[1]*b_gauss[0];
rhs[4]=-DN(1,0)*crhs19 - DN(1,1)*crhs20 + N[1]*b_gauss[1];
rhs[5]=-N[1]*crhs27 + crhs21*(DN(1,0)*crhs18 - DN(1,1)*crhs13) - crhs22*(DN(1,0)*crhs3 - DN(1,1)*crhs7) + crhs34*(DN(1,0)*crhs31 + DN(1,1)*crhs32) + crhs36*(DN(1,0)*crhs32 + DN(1,1)*crhs35);
rhs[6]=-DN(2,0)*crhs8 - DN(2,1)*crhs9 + N[2]*b_gauss[0];
rhs[7]=-DN(2,0)*crhs19 - DN(2,1)*crhs20 + N[2]*b_gauss[1];
rhs[8]=-N[2]*crhs27 + crhs21*(DN(2,0)*crhs18 - DN(2,1)*crhs13) - crhs22*(DN(2,0)*crhs3 - DN(2,1)*crhs7) + crhs34*(DN(2,0)*crhs31 + DN(2,1)*crhs32) + crhs36*(DN(2,0)*crhs32 + DN(2,1)*crhs35);

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
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau * std::pow(h,2) / 2.0;

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
        // TODO: THIS MUST BE COMPUTED BY THE CONSTITUTIVE LAW, ASSUMING LINEAR ELASTIC MATERIAL SO FAR
        const auto& r_prop = GetProperties();
        double mu = r_prop[YOUNG_MODULUS]/(2*(1.0+r_prop[POISSON_RATIO])); // 2nd Lame constant (shear modulus)
        const double tau = aux_tau / mu;

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
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau * std::pow(h,2) / 2.0;

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
        // TODO: THIS MUST BE COMPUTED BY THE CONSTITUTIVE LAW, ASSUMING LINEAR ELASTIC MATERIAL SO FAR
        const auto& r_prop = GetProperties();
        double mu = r_prop[YOUNG_MODULUS]/(2*(1.0+r_prop[POISSON_RATIO])); // 2nd Lame constant (shear modulus)
        const double tau = aux_tau / mu;

        // Calculate and add the LHS Gauss point contributions
        const double clhs0 = DN(0,1)*u(0,0);
const double clhs1 = DN(1,1)*u(1,0);
const double clhs2 = DN(2,1)*u(2,0);
const double clhs3 = clhs0 + clhs1 + clhs2;
const double clhs4 = DN(0,1)*clhs3;
const double clhs5 = DN(0,0)*u(0,0);
const double clhs6 = DN(1,1)*u(1,1);
const double clhs7 = DN(2,1)*u(2,1);
const double clhs8 = DN(0,1)*u(0,1);
const double clhs9 = DN(1,0)*u(1,0);
const double clhs10 = DN(2,0)*u(2,0);
const double clhs11 = DN(0,0)*u(0,1);
const double clhs12 = DN(1,0)*u(1,1);
const double clhs13 = DN(2,0)*u(2,1);
const double clhs14 = clhs10 + clhs5 + clhs9 + 1;
const double clhs15 = clhs6 + clhs7 + clhs8;
const double clhs16 = -clhs0*clhs12 - clhs0*clhs13 - clhs1*clhs11 - clhs1*clhs13 + clhs10*clhs6 + clhs10*clhs8 - clhs11*clhs2 - clhs12*clhs2 + clhs14 + clhs15 + clhs5*clhs6 + clhs5*clhs7 + clhs7*clhs9 + clhs8*clhs9;
const double clhs17 = 1.0/clhs16;
const double clhs18 = 1.0*clhs17;
const double clhs19 = pow(clhs3, 2);
const double clhs20 = pow(clhs16, -2.0);
const double clhs21 = DN(0,1)*clhs12;
const double clhs22 = DN(0,1)*clhs13;
const double clhs23 = DN(0,0)*clhs6;
const double clhs24 = DN(0,0)*clhs7;
const double clhs25 = clhs20*(-DN(0,0) + clhs21 + clhs22 - clhs23 - clhs24);
const double clhs26 = 0.5*clhs25;
const double clhs27 = clhs15 + 1;
const double clhs28 = pow(clhs27, 2);
const double clhs29 = clhs18*clhs4 + clhs19*clhs26 + clhs26*clhs28;
const double clhs30 = DN(0,0)*clhs14;
const double clhs31 = clhs11 + clhs12 + clhs13;
const double clhs32 = pow(clhs31, 2);
const double clhs33 = pow(clhs14, 2);
const double clhs34 = clhs18*clhs30 + clhs26*clhs32 + clhs26*clhs33;
const double clhs35 = DN(0,0)*clhs3;
const double clhs36 = DN(0,1)*clhs14;
const double clhs37 = clhs27*clhs31;
const double clhs38 = clhs14*clhs3;
const double clhs39 = clhs17*clhs35 + clhs17*clhs36 + clhs25*clhs37 + clhs25*clhs38;
const double clhs40 = C(0,0)*clhs34 + C(0,1)*clhs29 + C(0,2)*clhs39;
const double clhs41 = N[0]*th[0] + N[1]*th[1] + N[2]*th[2];
const double clhs42 = pow(clhs41, 2);
const double clhs43 = clhs14*clhs42;
const double clhs44 = C(0,2)*clhs34 + C(1,2)*clhs29 + C(2,2)*clhs39;
const double clhs45 = clhs3*clhs42;
const double clhs46 = clhs37 + clhs38;
const double clhs47 = clhs17*clhs42;
const double clhs48 = clhs46*clhs47;
const double clhs49 = 0.5*clhs32*clhs47 + 0.5*clhs33*clhs47 - 0.5;
const double clhs50 = 0.5*clhs19*clhs47 + 0.5*clhs28*clhs47 - 0.5;
const double clhs51 = C(0,0)*clhs49 + C(0,1)*clhs50 + C(0,2)*clhs48;
const double clhs52 = C(0,2)*clhs49 + C(1,2)*clhs50 + C(2,2)*clhs48;
const double clhs53 = DN(0,0)*clhs51 + DN(0,1)*clhs52;
const double clhs54 = clhs40*clhs43 + clhs44*clhs45 + clhs53;
const double clhs55 = C(0,1)*clhs34 + C(1,1)*clhs29 + C(1,2)*clhs39;
const double clhs56 = C(0,1)*clhs49 + C(1,1)*clhs50 + C(1,2)*clhs48;
const double clhs57 = DN(0,0)*clhs52 + DN(0,1)*clhs56;
const double clhs58 = clhs43*clhs44 + clhs45*clhs55 + clhs57;
const double clhs59 = DN(0,0)*clhs31;
const double clhs60 = DN(0,0)*clhs1;
const double clhs61 = DN(0,0)*clhs2;
const double clhs62 = DN(0,1)*clhs9;
const double clhs63 = DN(0,1)*clhs10;
const double clhs64 = clhs20*(-DN(0,1) + clhs60 + clhs61 - clhs62 - clhs63);
const double clhs65 = 0.5*clhs64;
const double clhs66 = clhs18*clhs59 + clhs32*clhs65 + clhs33*clhs65;
const double clhs67 = DN(0,1)*clhs27;
const double clhs68 = clhs18*clhs67 + clhs19*clhs65 + clhs28*clhs65;
const double clhs69 = DN(0,1)*clhs31;
const double clhs70 = DN(0,0)*clhs27;
const double clhs71 = clhs17*clhs69 + clhs17*clhs70 + clhs37*clhs64 + clhs38*clhs64;
const double clhs72 = C(0,2)*clhs66 + C(1,2)*clhs68 + C(2,2)*clhs71;
const double clhs73 = C(0,0)*clhs66 + C(0,1)*clhs68 + C(0,2)*clhs71;
const double clhs74 = clhs14*clhs73 + clhs3*clhs72;
const double clhs75 = C(0,1)*clhs66 + C(1,1)*clhs68 + C(1,2)*clhs71;
const double clhs76 = clhs14*clhs72 + clhs3*clhs75;
const double clhs77 = clhs32 + clhs33;
const double clhs78 = C(0,2)*clhs77;
const double clhs79 = clhs19 + clhs28;
const double clhs80 = C(1,2)*clhs79;
const double clhs81 = 2.0*clhs46;
const double clhs82 = C(2,2)*clhs81 + 1.0*clhs78 + 1.0*clhs80;
const double clhs83 = C(0,0)*clhs77;
const double clhs84 = C(0,1)*clhs79;
const double clhs85 = C(0,2)*clhs81 + 1.0*clhs83 + 1.0*clhs84;
const double clhs86 = clhs14*clhs85 + clhs3*clhs82;
const double clhs87 = C(0,1)*clhs77;
const double clhs88 = C(1,1)*clhs79;
const double clhs89 = C(1,2)*clhs81 + 1.0*clhs87 + 1.0*clhs88;
const double clhs90 = clhs14*clhs82 + clhs3*clhs89;
const double clhs91 = clhs17*clhs41;
const double clhs92 = clhs91*(DN(0,0)*clhs86 + DN(0,1)*clhs90);
const double clhs93 = DN(1,1)*clhs3;
const double clhs94 = DN(1,1)*clhs11;
const double clhs95 = DN(1,1)*clhs13;
const double clhs96 = DN(1,0)*clhs8;
const double clhs97 = DN(1,0)*clhs7;
const double clhs98 = clhs20*(-DN(1,0) + clhs94 + clhs95 - clhs96 - clhs97);
const double clhs99 = 0.5*clhs98;
const double clhs100 = clhs18*clhs93 + clhs19*clhs99 + clhs28*clhs99;
const double clhs101 = DN(1,0)*clhs14;
const double clhs102 = clhs101*clhs18 + clhs32*clhs99 + clhs33*clhs99;
const double clhs103 = DN(1,0)*clhs3;
const double clhs104 = DN(1,1)*clhs14;
const double clhs105 = clhs103*clhs17 + clhs104*clhs17 + clhs37*clhs98 + clhs38*clhs98;
const double clhs106 = C(0,0)*clhs102 + C(0,1)*clhs100 + C(0,2)*clhs105;
const double clhs107 = C(0,2)*clhs102 + C(1,2)*clhs100 + C(2,2)*clhs105;
const double clhs108 = DN(1,0)*clhs51 + DN(1,1)*clhs52;
const double clhs109 = clhs106*clhs43 + clhs107*clhs45 + clhs108;
const double clhs110 = C(0,1)*clhs102 + C(1,1)*clhs100 + C(1,2)*clhs105;
const double clhs111 = DN(1,0)*clhs52 + DN(1,1)*clhs56;
const double clhs112 = clhs107*clhs43 + clhs110*clhs45 + clhs111;
const double clhs113 = DN(1,0)*clhs31;
const double clhs114 = DN(1,0)*clhs0;
const double clhs115 = DN(1,0)*clhs2;
const double clhs116 = DN(1,1)*clhs5;
const double clhs117 = DN(1,1)*clhs10;
const double clhs118 = clhs20*(-DN(1,1) + clhs114 + clhs115 - clhs116 - clhs117);
const double clhs119 = 0.5*clhs118;
const double clhs120 = clhs113*clhs18 + clhs119*clhs32 + clhs119*clhs33;
const double clhs121 = DN(1,1)*clhs27;
const double clhs122 = clhs119*clhs19 + clhs119*clhs28 + clhs121*clhs18;
const double clhs123 = DN(1,1)*clhs31;
const double clhs124 = DN(1,0)*clhs27;
const double clhs125 = clhs118*clhs37 + clhs118*clhs38 + clhs123*clhs17 + clhs124*clhs17;
const double clhs126 = C(0,2)*clhs120 + C(1,2)*clhs122 + C(2,2)*clhs125;
const double clhs127 = C(0,0)*clhs120 + C(0,1)*clhs122 + C(0,2)*clhs125;
const double clhs128 = clhs126*clhs3 + clhs127*clhs14;
const double clhs129 = C(0,1)*clhs120 + C(1,1)*clhs122 + C(1,2)*clhs125;
const double clhs130 = clhs126*clhs14 + clhs129*clhs3;
const double clhs131 = DN(2,1)*clhs3;
const double clhs132 = DN(2,1)*clhs11;
const double clhs133 = DN(2,1)*clhs12;
const double clhs134 = DN(2,0)*clhs8;
const double clhs135 = DN(2,0)*clhs6;
const double clhs136 = clhs20*(-DN(2,0) + clhs132 + clhs133 - clhs134 - clhs135);
const double clhs137 = 0.5*clhs136;
const double clhs138 = clhs131*clhs18 + clhs137*clhs19 + clhs137*clhs28;
const double clhs139 = DN(2,0)*clhs14;
const double clhs140 = clhs137*clhs32 + clhs137*clhs33 + clhs139*clhs18;
const double clhs141 = DN(2,0)*clhs3;
const double clhs142 = DN(2,1)*clhs14;
const double clhs143 = clhs136*clhs37 + clhs136*clhs38 + clhs141*clhs17 + clhs142*clhs17;
const double clhs144 = C(0,0)*clhs140 + C(0,1)*clhs138 + C(0,2)*clhs143;
const double clhs145 = C(0,2)*clhs140 + C(1,2)*clhs138 + C(2,2)*clhs143;
const double clhs146 = DN(2,0)*clhs51 + DN(2,1)*clhs52;
const double clhs147 = clhs144*clhs43 + clhs145*clhs45 + clhs146;
const double clhs148 = C(0,1)*clhs140 + C(1,1)*clhs138 + C(1,2)*clhs143;
const double clhs149 = DN(2,0)*clhs52 + DN(2,1)*clhs56;
const double clhs150 = clhs145*clhs43 + clhs148*clhs45 + clhs149;
const double clhs151 = DN(2,0)*clhs31;
const double clhs152 = DN(2,0)*clhs0;
const double clhs153 = DN(2,0)*clhs1;
const double clhs154 = DN(2,1)*clhs5;
const double clhs155 = DN(2,1)*clhs9;
const double clhs156 = clhs20*(-DN(2,1) + clhs152 + clhs153 - clhs154 - clhs155);
const double clhs157 = 0.5*clhs156;
const double clhs158 = clhs151*clhs18 + clhs157*clhs32 + clhs157*clhs33;
const double clhs159 = DN(2,1)*clhs27;
const double clhs160 = clhs157*clhs19 + clhs157*clhs28 + clhs159*clhs18;
const double clhs161 = DN(2,1)*clhs31;
const double clhs162 = DN(2,0)*clhs27;
const double clhs163 = clhs156*clhs37 + clhs156*clhs38 + clhs161*clhs17 + clhs162*clhs17;
const double clhs164 = C(0,2)*clhs158 + C(1,2)*clhs160 + C(2,2)*clhs163;
const double clhs165 = C(0,0)*clhs158 + C(0,1)*clhs160 + C(0,2)*clhs163;
const double clhs166 = clhs14*clhs165 + clhs164*clhs3;
const double clhs167 = C(0,1)*clhs158 + C(1,1)*clhs160 + C(1,2)*clhs163;
const double clhs168 = clhs14*clhs164 + clhs167*clhs3;
const double clhs169 = clhs27*clhs44 + clhs31*clhs40;
const double clhs170 = clhs27*clhs55 + clhs31*clhs44;
const double clhs171 = clhs31*clhs42;
const double clhs172 = clhs27*clhs42;
const double clhs173 = clhs171*clhs73 + clhs172*clhs72 + clhs53;
const double clhs174 = clhs171*clhs72 + clhs172*clhs75 + clhs57;
const double clhs175 = clhs27*clhs82 + clhs31*clhs85;
const double clhs176 = clhs27*clhs89 + clhs31*clhs82;
const double clhs177 = clhs91*(DN(0,0)*clhs175 + DN(0,1)*clhs176);
const double clhs178 = clhs106*clhs31 + clhs107*clhs27;
const double clhs179 = clhs107*clhs31 + clhs110*clhs27;
const double clhs180 = clhs108 + clhs126*clhs172 + clhs127*clhs171;
const double clhs181 = clhs111 + clhs126*clhs171 + clhs129*clhs172;
const double clhs182 = clhs144*clhs31 + clhs145*clhs27;
const double clhs183 = clhs145*clhs31 + clhs148*clhs27;
const double clhs184 = clhs146 + clhs164*clhs172 + clhs165*clhs171;
const double clhs185 = clhs149 + clhs164*clhs171 + clhs167*clhs172;
const double clhs186 = DN(0,0) - clhs21 - clhs22 + clhs23 + clhs24;
const double clhs187 = clhs35 + clhs36;
const double clhs188 = C(0,0)*clhs30 + C(0,1)*clhs4 + C(0,2)*clhs187;
const double clhs189 = C(0,2)*clhs30 + C(1,2)*clhs4 + C(2,2)*clhs187;
const double clhs190 = 2*clhs46;
const double clhs191 = C(0,2)*clhs190 + clhs83 + clhs84;
const double clhs192 = DN(0,0)*clhs191;
const double clhs193 = 0.25/clhs16;
const double clhs194 = clhs186*clhs193;
const double clhs195 = C(2,2)*clhs190 + clhs78 + clhs80;
const double clhs196 = DN(0,1)*clhs195;
const double clhs197 = 1.0/clhs41;
const double clhs198 = tau*sqrt(clhs16*clhs197);
const double clhs199 = clhs198*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double clhs200 = C(0,1)*clhs30 + C(1,1)*clhs4 + C(1,2)*clhs187;
const double clhs201 = DN(0,0)*clhs195;
const double clhs202 = C(1,2)*clhs190 + clhs87 + clhs88;
const double clhs203 = DN(0,1)*clhs202;
const double clhs204 = clhs198*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double clhs205 = DN(0,1) - clhs60 - clhs61 + clhs62 + clhs63;
const double clhs206 = clhs69 + clhs70;
const double clhs207 = C(0,0)*clhs59 + C(0,1)*clhs67 + C(0,2)*clhs206;
const double clhs208 = C(0,2)*clhs59 + C(1,2)*clhs67 + C(2,2)*clhs206;
const double clhs209 = clhs193*clhs205;
const double clhs210 = C(0,1)*clhs59 + C(1,1)*clhs67 + C(1,2)*clhs206;
const double clhs211 = clhs192 + clhs196;
const double clhs212 = (1.0/2.0)*clhs198;
const double clhs213 = clhs211*clhs212;
const double clhs214 = clhs201 + clhs203;
const double clhs215 = clhs212*clhs214;
const double clhs216 = 0.25*clhs197;
const double clhs217 = N[0]*clhs216;
const double clhs218 = clhs199*clhs211;
const double clhs219 = clhs204*clhs214;
const double clhs220 = tau*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs221 = b_gauss[1]*clhs220;
const double clhs222 = DN(1,0) - clhs94 - clhs95 + clhs96 + clhs97;
const double clhs223 = clhs103 + clhs104;
const double clhs224 = C(0,0)*clhs101 + C(0,1)*clhs93 + C(0,2)*clhs223;
const double clhs225 = C(0,2)*clhs101 + C(1,2)*clhs93 + C(2,2)*clhs223;
const double clhs226 = clhs193*clhs222;
const double clhs227 = C(0,1)*clhs101 + C(1,1)*clhs93 + C(1,2)*clhs223;
const double clhs228 = b_gauss[0]*clhs220;
const double clhs229 = DN(1,1) - clhs114 - clhs115 + clhs116 + clhs117;
const double clhs230 = clhs123 + clhs124;
const double clhs231 = C(0,0)*clhs113 + C(0,1)*clhs121 + C(0,2)*clhs230;
const double clhs232 = C(0,2)*clhs113 + C(1,2)*clhs121 + C(2,2)*clhs230;
const double clhs233 = clhs193*clhs229;
const double clhs234 = C(0,1)*clhs113 + C(1,1)*clhs121 + C(1,2)*clhs230;
const double clhs235 = -N[0]*N[1];
const double clhs236 = N[1]*clhs216;
const double clhs237 = tau*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs238 = b_gauss[1]*clhs237;
const double clhs239 = DN(2,0) - clhs132 - clhs133 + clhs134 + clhs135;
const double clhs240 = clhs141 + clhs142;
const double clhs241 = C(0,0)*clhs139 + C(0,1)*clhs131 + C(0,2)*clhs240;
const double clhs242 = C(0,2)*clhs139 + C(1,2)*clhs131 + C(2,2)*clhs240;
const double clhs243 = clhs193*clhs239;
const double clhs244 = C(0,1)*clhs139 + C(1,1)*clhs131 + C(1,2)*clhs240;
const double clhs245 = b_gauss[0]*clhs237;
const double clhs246 = DN(2,1) - clhs152 - clhs153 + clhs154 + clhs155;
const double clhs247 = clhs161 + clhs162;
const double clhs248 = C(0,0)*clhs151 + C(0,1)*clhs159 + C(0,2)*clhs247;
const double clhs249 = C(0,2)*clhs151 + C(1,2)*clhs159 + C(2,2)*clhs247;
const double clhs250 = clhs193*clhs246;
const double clhs251 = C(0,1)*clhs151 + C(1,1)*clhs159 + C(1,2)*clhs247;
const double clhs252 = -N[0]*N[2];
const double clhs253 = N[2]*clhs216;
const double clhs254 = clhs91*(DN(1,0)*clhs86 + DN(1,1)*clhs90);
const double clhs255 = clhs91*(DN(1,0)*clhs175 + DN(1,1)*clhs176);
const double clhs256 = DN(1,0)*clhs191;
const double clhs257 = DN(1,1)*clhs195;
const double clhs258 = DN(1,0)*clhs195;
const double clhs259 = DN(1,1)*clhs202;
const double clhs260 = clhs256 + clhs257;
const double clhs261 = clhs212*clhs260;
const double clhs262 = clhs258 + clhs259;
const double clhs263 = clhs212*clhs262;
const double clhs264 = clhs199*clhs260;
const double clhs265 = clhs204*clhs262;
const double clhs266 = tau*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs267 = b_gauss[1]*clhs266;
const double clhs268 = b_gauss[0]*clhs266;
const double clhs269 = -N[1]*N[2];
const double clhs270 = clhs91*(DN(2,0)*clhs86 + DN(2,1)*clhs90);
const double clhs271 = clhs91*(DN(2,0)*clhs175 + DN(2,1)*clhs176);
const double clhs272 = DN(2,0)*clhs191;
const double clhs273 = DN(2,1)*clhs195;
const double clhs274 = DN(2,0)*clhs195;
const double clhs275 = DN(2,1)*clhs202;
const double clhs276 = clhs272 + clhs273;
const double clhs277 = clhs212*clhs276;
const double clhs278 = clhs274 + clhs275;
const double clhs279 = clhs212*clhs278;
const double clhs280 = clhs199*clhs276;
const double clhs281 = clhs204*clhs278;
lhs(0,0)=DN(0,0)*clhs54 + DN(0,1)*clhs58;
lhs(0,1)=clhs42*(DN(0,0)*clhs74 + DN(0,1)*clhs76);
lhs(0,2)=N[0]*clhs92;
lhs(0,3)=DN(0,0)*clhs109 + DN(0,1)*clhs112;
lhs(0,4)=clhs42*(DN(0,0)*clhs128 + DN(0,1)*clhs130);
lhs(0,5)=N[1]*clhs92;
lhs(0,6)=DN(0,0)*clhs147 + DN(0,1)*clhs150;
lhs(0,7)=clhs42*(DN(0,0)*clhs166 + DN(0,1)*clhs168);
lhs(0,8)=N[2]*clhs92;
lhs(1,0)=clhs42*(DN(0,0)*clhs169 + DN(0,1)*clhs170);
lhs(1,1)=DN(0,0)*clhs173 + DN(0,1)*clhs174;
lhs(1,2)=N[0]*clhs177;
lhs(1,3)=clhs42*(DN(0,0)*clhs178 + DN(0,1)*clhs179);
lhs(1,4)=DN(0,0)*clhs180 + DN(0,1)*clhs181;
lhs(1,5)=N[1]*clhs177;
lhs(1,6)=clhs42*(DN(0,0)*clhs182 + DN(0,1)*clhs183);
lhs(1,7)=DN(0,0)*clhs184 + DN(0,1)*clhs185;
lhs(1,8)=N[2]*clhs177;
lhs(2,0)=N[0]*clhs186 - clhs199*(DN(0,0)*clhs188 + DN(0,1)*clhs189 + clhs192*clhs194 + clhs194*clhs196) - clhs204*(DN(0,0)*clhs189 + DN(0,1)*clhs200 + clhs194*clhs201 + clhs194*clhs203);
lhs(2,1)=N[0]*clhs205 - clhs199*(DN(0,0)*clhs207 + DN(0,1)*clhs208 + clhs192*clhs209 + clhs196*clhs209) - clhs204*(DN(0,0)*clhs208 + DN(0,1)*clhs210 + clhs201*clhs209 + clhs203*clhs209);
lhs(2,2)=-DN(0,0)*clhs213 - DN(0,1)*clhs215 - pow(N[0], 2) + clhs217*clhs218 + clhs217*clhs219;
lhs(2,3)=N[0]*clhs222 - clhs199*(DN(0,0)*clhs224 + DN(0,1)*clhs225 + clhs192*clhs226 + clhs196*clhs226) - clhs204*(DN(0,0)*clhs225 + DN(0,1)*clhs227 + clhs201*clhs226 + clhs203*clhs226) + clhs221;
lhs(2,4)=N[0]*clhs229 - clhs199*(DN(0,0)*clhs231 + DN(0,1)*clhs232 + clhs192*clhs233 + clhs196*clhs233) - clhs204*(DN(0,0)*clhs232 + DN(0,1)*clhs234 + clhs201*clhs233 + clhs203*clhs233) - clhs228;
lhs(2,5)=-DN(1,0)*clhs213 - DN(1,1)*clhs215 + clhs218*clhs236 + clhs219*clhs236 + clhs235;
lhs(2,6)=N[0]*clhs239 - clhs199*(DN(0,0)*clhs241 + DN(0,1)*clhs242 + clhs192*clhs243 + clhs196*clhs243) - clhs204*(DN(0,0)*clhs242 + DN(0,1)*clhs244 + clhs201*clhs243 + clhs203*clhs243) + clhs238;
lhs(2,7)=N[0]*clhs246 - clhs199*(DN(0,0)*clhs248 + DN(0,1)*clhs249 + clhs192*clhs250 + clhs196*clhs250) - clhs204*(DN(0,0)*clhs249 + DN(0,1)*clhs251 + clhs201*clhs250 + clhs203*clhs250) - clhs245;
lhs(2,8)=-DN(2,0)*clhs213 - DN(2,1)*clhs215 + clhs218*clhs253 + clhs219*clhs253 + clhs252;
lhs(3,0)=DN(1,0)*clhs54 + DN(1,1)*clhs58;
lhs(3,1)=clhs42*(DN(1,0)*clhs74 + DN(1,1)*clhs76);
lhs(3,2)=N[0]*clhs254;
lhs(3,3)=DN(1,0)*clhs109 + DN(1,1)*clhs112;
lhs(3,4)=clhs42*(DN(1,0)*clhs128 + DN(1,1)*clhs130);
lhs(3,5)=N[1]*clhs254;
lhs(3,6)=DN(1,0)*clhs147 + DN(1,1)*clhs150;
lhs(3,7)=clhs42*(DN(1,0)*clhs166 + DN(1,1)*clhs168);
lhs(3,8)=N[2]*clhs254;
lhs(4,0)=clhs42*(DN(1,0)*clhs169 + DN(1,1)*clhs170);
lhs(4,1)=DN(1,0)*clhs173 + DN(1,1)*clhs174;
lhs(4,2)=N[0]*clhs255;
lhs(4,3)=clhs42*(DN(1,0)*clhs178 + DN(1,1)*clhs179);
lhs(4,4)=DN(1,0)*clhs180 + DN(1,1)*clhs181;
lhs(4,5)=N[1]*clhs255;
lhs(4,6)=clhs42*(DN(1,0)*clhs182 + DN(1,1)*clhs183);
lhs(4,7)=DN(1,0)*clhs184 + DN(1,1)*clhs185;
lhs(4,8)=N[2]*clhs255;
lhs(5,0)=N[1]*clhs186 - clhs199*(DN(1,0)*clhs188 + DN(1,1)*clhs189 + clhs194*clhs256 + clhs194*clhs257) - clhs204*(DN(1,0)*clhs189 + DN(1,1)*clhs200 + clhs194*clhs258 + clhs194*clhs259) - clhs221;
lhs(5,1)=N[1]*clhs205 - clhs199*(DN(1,0)*clhs207 + DN(1,1)*clhs208 + clhs209*clhs256 + clhs209*clhs257) - clhs204*(DN(1,0)*clhs208 + DN(1,1)*clhs210 + clhs209*clhs258 + clhs209*clhs259) + clhs228;
lhs(5,2)=-DN(0,0)*clhs261 - DN(0,1)*clhs263 + clhs217*clhs264 + clhs217*clhs265 + clhs235;
lhs(5,3)=N[1]*clhs222 - clhs199*(DN(1,0)*clhs224 + DN(1,1)*clhs225 + clhs226*clhs256 + clhs226*clhs257) - clhs204*(DN(1,0)*clhs225 + DN(1,1)*clhs227 + clhs226*clhs258 + clhs226*clhs259);
lhs(5,4)=N[1]*clhs229 - clhs199*(DN(1,0)*clhs231 + DN(1,1)*clhs232 + clhs233*clhs256 + clhs233*clhs257) - clhs204*(DN(1,0)*clhs232 + DN(1,1)*clhs234 + clhs233*clhs258 + clhs233*clhs259);
lhs(5,5)=-DN(1,0)*clhs261 - DN(1,1)*clhs263 - pow(N[1], 2) + clhs236*clhs264 + clhs236*clhs265;
lhs(5,6)=N[1]*clhs239 - clhs199*(DN(1,0)*clhs241 + DN(1,1)*clhs242 + clhs243*clhs256 + clhs243*clhs257) - clhs204*(DN(1,0)*clhs242 + DN(1,1)*clhs244 + clhs243*clhs258 + clhs243*clhs259) + clhs267;
lhs(5,7)=N[1]*clhs246 - clhs199*(DN(1,0)*clhs248 + DN(1,1)*clhs249 + clhs250*clhs256 + clhs250*clhs257) - clhs204*(DN(1,0)*clhs249 + DN(1,1)*clhs251 + clhs250*clhs258 + clhs250*clhs259) - clhs268;
lhs(5,8)=-DN(2,0)*clhs261 - DN(2,1)*clhs263 + clhs253*clhs264 + clhs253*clhs265 + clhs269;
lhs(6,0)=DN(2,0)*clhs54 + DN(2,1)*clhs58;
lhs(6,1)=clhs42*(DN(2,0)*clhs74 + DN(2,1)*clhs76);
lhs(6,2)=N[0]*clhs270;
lhs(6,3)=DN(2,0)*clhs109 + DN(2,1)*clhs112;
lhs(6,4)=clhs42*(DN(2,0)*clhs128 + DN(2,1)*clhs130);
lhs(6,5)=N[1]*clhs270;
lhs(6,6)=DN(2,0)*clhs147 + DN(2,1)*clhs150;
lhs(6,7)=clhs42*(DN(2,0)*clhs166 + DN(2,1)*clhs168);
lhs(6,8)=N[2]*clhs270;
lhs(7,0)=clhs42*(DN(2,0)*clhs169 + DN(2,1)*clhs170);
lhs(7,1)=DN(2,0)*clhs173 + DN(2,1)*clhs174;
lhs(7,2)=N[0]*clhs271;
lhs(7,3)=clhs42*(DN(2,0)*clhs178 + DN(2,1)*clhs179);
lhs(7,4)=DN(2,0)*clhs180 + DN(2,1)*clhs181;
lhs(7,5)=N[1]*clhs271;
lhs(7,6)=clhs42*(DN(2,0)*clhs182 + DN(2,1)*clhs183);
lhs(7,7)=DN(2,0)*clhs184 + DN(2,1)*clhs185;
lhs(7,8)=N[2]*clhs271;
lhs(8,0)=N[2]*clhs186 - clhs199*(DN(2,0)*clhs188 + DN(2,1)*clhs189 + clhs194*clhs272 + clhs194*clhs273) - clhs204*(DN(2,0)*clhs189 + DN(2,1)*clhs200 + clhs194*clhs274 + clhs194*clhs275) - clhs238;
lhs(8,1)=N[2]*clhs205 - clhs199*(DN(2,0)*clhs207 + DN(2,1)*clhs208 + clhs209*clhs272 + clhs209*clhs273) - clhs204*(DN(2,0)*clhs208 + DN(2,1)*clhs210 + clhs209*clhs274 + clhs209*clhs275) + clhs245;
lhs(8,2)=-DN(0,0)*clhs277 - DN(0,1)*clhs279 + clhs217*clhs280 + clhs217*clhs281 + clhs252;
lhs(8,3)=N[2]*clhs222 - clhs199*(DN(2,0)*clhs224 + DN(2,1)*clhs225 + clhs226*clhs272 + clhs226*clhs273) - clhs204*(DN(2,0)*clhs225 + DN(2,1)*clhs227 + clhs226*clhs274 + clhs226*clhs275) - clhs267;
lhs(8,4)=N[2]*clhs229 - clhs199*(DN(2,0)*clhs231 + DN(2,1)*clhs232 + clhs233*clhs272 + clhs233*clhs273) - clhs204*(DN(2,0)*clhs232 + DN(2,1)*clhs234 + clhs233*clhs274 + clhs233*clhs275) + clhs268;
lhs(8,5)=-DN(1,0)*clhs277 - DN(1,1)*clhs279 + clhs236*clhs280 + clhs236*clhs281 + clhs269;
lhs(8,6)=N[2]*clhs239 - clhs199*(DN(2,0)*clhs241 + DN(2,1)*clhs242 + clhs243*clhs272 + clhs243*clhs273) - clhs204*(DN(2,0)*clhs242 + DN(2,1)*clhs244 + clhs243*clhs274 + clhs243*clhs275);
lhs(8,7)=N[2]*clhs246 - clhs199*(DN(2,0)*clhs248 + DN(2,1)*clhs249 + clhs250*clhs272 + clhs250*clhs273) - clhs204*(DN(2,0)*clhs249 + DN(2,1)*clhs251 + clhs250*clhs274 + clhs250*clhs275);
lhs(8,8)=-DN(2,0)*clhs277 - DN(2,1)*clhs279 - pow(N[2], 2) + clhs253*clhs280 + clhs253*clhs281;

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
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check LHS size
    if (rLeftHandSideMatrix.size1() != matrix_size || rLeftHandSideMatrix.size2() != matrix_size) {
        rLeftHandSideMatrix.resize(matrix_size, matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau * std::pow(h,2) / 2.0;

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
        // TODO: THIS MUST BE COMPUTED BY THE CONSTITUTIVE LAW, ASSUMING LINEAR ELASTIC MATERIAL SO FAR
        const auto& r_prop = GetProperties();
        double mu = r_prop[YOUNG_MODULUS]/(2*(1.0+r_prop[POISSON_RATIO])); // 2nd Lame constant (shear modulus)
        const double tau = aux_tau / mu;

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
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau * std::pow(h,2) / 2.0;

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
        // TODO: THIS MUST BE COMPUTED BY THE CONSTITUTIVE LAW, ASSUMING LINEAR ELASTIC MATERIAL SO FAR
        const auto& r_prop = GetProperties();
        double mu = r_prop[YOUNG_MODULUS]/(2*(1.0+r_prop[POISSON_RATIO])); // 2nd Lame constant (shear modulus)
        const double tau = aux_tau / mu;

        // Calculate and add the RHS Gauss point contribution
        const double crhs0 = DN(0,1)*u(0,0);
const double crhs1 = DN(1,1)*u(1,0);
const double crhs2 = DN(2,1)*u(2,0);
const double crhs3 = crhs0 + crhs1 + crhs2;
const double crhs4 = DN(0,0)*u(0,0);
const double crhs5 = DN(1,0)*u(1,0);
const double crhs6 = DN(2,0)*u(2,0);
const double crhs7 = crhs4 + crhs5 + crhs6 + 1;
const double crhs8 = S[0]*crhs7 + S[2]*crhs3;
const double crhs9 = S[1]*crhs3 + S[2]*crhs7;
const double crhs10 = DN(0,0)*u(0,1);
const double crhs11 = DN(1,0)*u(1,1);
const double crhs12 = DN(2,0)*u(2,1);
const double crhs13 = crhs10 + crhs11 + crhs12;
const double crhs14 = DN(0,1)*u(0,1);
const double crhs15 = DN(1,1)*u(1,1);
const double crhs16 = DN(2,1)*u(2,1);
const double crhs17 = crhs14 + crhs15 + crhs16;
const double crhs18 = crhs17 + 1;
const double crhs19 = S[0]*crhs13 + S[2]*crhs18;
const double crhs20 = S[1]*crhs18 + S[2]*crhs13;
const double crhs21 = b_gauss[0]*tau;
const double crhs22 = b_gauss[1]*tau;
const double crhs23 = N[0]*th[0];
const double crhs24 = N[1]*th[1];
const double crhs25 = N[2]*th[2];
const double crhs26 = -crhs0*crhs11 - crhs0*crhs12 - crhs1*crhs10 - crhs1*crhs12 - crhs10*crhs2 - crhs11*crhs2 + crhs14*crhs5 + crhs14*crhs6 + crhs15*crhs4 + crhs15*crhs6 + crhs16*crhs4 + crhs16*crhs5 + crhs17 + crhs7;
const double crhs27 = -crhs23 - crhs24 - crhs25 + crhs26;
const double crhs28 = pow(crhs13, 2) + pow(crhs7, 2);
const double crhs29 = pow(crhs18, 2) + pow(crhs3, 2);
const double crhs30 = 2*crhs13*crhs18 + 2*crhs3*crhs7;
const double crhs31 = C(0,0)*crhs28 + C(0,1)*crhs29 + C(0,2)*crhs30;
const double crhs32 = C(0,2)*crhs28 + C(1,2)*crhs29 + C(2,2)*crhs30;
const double crhs33 = (1.0/2.0)*tau*sqrt(crhs26/(crhs23 + crhs24 + crhs25));
const double crhs34 = crhs33*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs35 = C(0,1)*crhs28 + C(1,1)*crhs29 + C(1,2)*crhs30;
const double crhs36 = crhs33*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
rhs[0]=-DN(0,0)*crhs8 - DN(0,1)*crhs9 + N[0]*b_gauss[0];
rhs[1]=-DN(0,0)*crhs19 - DN(0,1)*crhs20 + N[0]*b_gauss[1];
rhs[2]=-N[0]*crhs27 + crhs21*(DN(0,0)*crhs18 - DN(0,1)*crhs13) - crhs22*(DN(0,0)*crhs3 - DN(0,1)*crhs7) + crhs34*(DN(0,0)*crhs31 + DN(0,1)*crhs32) + crhs36*(DN(0,0)*crhs32 + DN(0,1)*crhs35);
rhs[3]=-DN(1,0)*crhs8 - DN(1,1)*crhs9 + N[1]*b_gauss[0];
rhs[4]=-DN(1,0)*crhs19 - DN(1,1)*crhs20 + N[1]*b_gauss[1];
rhs[5]=-N[1]*crhs27 + crhs21*(DN(1,0)*crhs18 - DN(1,1)*crhs13) - crhs22*(DN(1,0)*crhs3 - DN(1,1)*crhs7) + crhs34*(DN(1,0)*crhs31 + DN(1,1)*crhs32) + crhs36*(DN(1,0)*crhs32 + DN(1,1)*crhs35);
rhs[6]=-DN(2,0)*crhs8 - DN(2,1)*crhs9 + N[2]*b_gauss[0];
rhs[7]=-DN(2,0)*crhs19 - DN(2,1)*crhs20 + N[2]*b_gauss[1];
rhs[8]=-N[2]*crhs27 + crhs21*(DN(2,0)*crhs18 - DN(2,1)*crhs13) - crhs22*(DN(2,0)*crhs3 - DN(2,1)*crhs7) + crhs34*(DN(2,0)*crhs31 + DN(2,1)*crhs32) + crhs36*(DN(2,0)*crhs32 + DN(2,1)*crhs35);

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
    const SizeType dim = r_geometry.WorkingSpaceDimension();
    const SizeType n_nodes = r_geometry.PointsNumber();
    const SizeType block_size = dim + 1;
    const SizeType matrix_size = block_size * n_nodes;
    const SizeType strain_size = GetProperties().GetValue(CONSTITUTIVE_LAW)->GetStrainSize();

    // Check RHS size
    if (rRightHandSideVector.size() != matrix_size) {
        rRightHandSideVector.resize(matrix_size, false);
    }

    // Create the kinematics container and fill the nodal data
    KinematicVariables kinematic_variables;
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
        for (IndexType d = 0; d < dim; ++d) {
            kinematic_variables.Displacements(i_node * dim + d) = r_disp[d];
        }
        kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
    const double c_tau = 2.0;
    const double h = ElementSizeCalculator<2,NumNodes>::AverageElementSize(r_geometry);
    const double aux_tau = c_tau * std::pow(h,2) / 2.0;

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
        // TODO: THIS MUST BE COMPUTED BY THE CONSTITUTIVE LAW, ASSUMING LINEAR ELASTIC MATERIAL SO FAR
        const auto& r_prop = GetProperties();
        double mu = r_prop[YOUNG_MODULUS]/(2*(1.0+r_prop[POISSON_RATIO])); // 2nd Lame constant (shear modulus)
        const double tau = aux_tau / mu;

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
const double cr_eq_green_strain15 = pow(N[0]*th[0] + N[1]*th[1] + N[2]*th[2], 2)*1.0/(-cr_eq_green_strain0*cr_eq_green_strain10 - cr_eq_green_strain0*cr_eq_green_strain11 - cr_eq_green_strain1*cr_eq_green_strain11 - cr_eq_green_strain1*cr_eq_green_strain12 - cr_eq_green_strain10*cr_eq_green_strain2 - cr_eq_green_strain12*cr_eq_green_strain2 + cr_eq_green_strain13 + cr_eq_green_strain14 + cr_eq_green_strain4*cr_eq_green_strain5 + cr_eq_green_strain4*cr_eq_green_strain6 + cr_eq_green_strain5*cr_eq_green_strain9 + cr_eq_green_strain6*cr_eq_green_strain8 + cr_eq_green_strain7*cr_eq_green_strain8 + cr_eq_green_strain7*cr_eq_green_strain9);
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
const double cr_eq_def_gradient14 = (N[0]*th[0] + N[1]*th[1] + N[2]*th[2])*pow(cr_eq_def_gradient0*cr_eq_def_gradient4 + cr_eq_def_gradient0*cr_eq_def_gradient5 + cr_eq_def_gradient1*cr_eq_def_gradient5 + cr_eq_def_gradient1*cr_eq_def_gradient6 - cr_eq_def_gradient10*cr_eq_def_gradient11 - cr_eq_def_gradient10*cr_eq_def_gradient12 - cr_eq_def_gradient11*cr_eq_def_gradient9 - cr_eq_def_gradient12*cr_eq_def_gradient7 + cr_eq_def_gradient13 + cr_eq_def_gradient2*cr_eq_def_gradient4 + cr_eq_def_gradient2*cr_eq_def_gradient6 + cr_eq_def_gradient3 - cr_eq_def_gradient7*cr_eq_def_gradient8 - cr_eq_def_gradient8*cr_eq_def_gradient9, -0.5);
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
const double cr_det_eq_def_gradient6 = cr_det_eq_def_gradient0*cr_det_eq_def_gradient3;
const double cr_det_eq_def_gradient7 = cr_det_eq_def_gradient0*cr_det_eq_def_gradient5;
const double cr_det_eq_def_gradient8 = cr_det_eq_def_gradient1*cr_det_eq_def_gradient2;
const double cr_det_eq_def_gradient9 = cr_det_eq_def_gradient1*cr_det_eq_def_gradient4;
const double cr_det_eq_def_gradient10 = cr_det_eq_def_gradient2*cr_det_eq_def_gradient5;
const double cr_det_eq_def_gradient11 = cr_det_eq_def_gradient3*cr_det_eq_def_gradient4;
const double cr_det_eq_def_gradient12 = DN(0,0)*u(0,1);
const double cr_det_eq_def_gradient13 = DN(1,1)*u(1,0);
const double cr_det_eq_def_gradient14 = cr_det_eq_def_gradient12*cr_det_eq_def_gradient13;
const double cr_det_eq_def_gradient15 = DN(2,1)*u(2,0);
const double cr_det_eq_def_gradient16 = cr_det_eq_def_gradient12*cr_det_eq_def_gradient15;
const double cr_det_eq_def_gradient17 = DN(0,1)*u(0,0);
const double cr_det_eq_def_gradient18 = DN(1,0)*u(1,1);
const double cr_det_eq_def_gradient19 = cr_det_eq_def_gradient17*cr_det_eq_def_gradient18;
const double cr_det_eq_def_gradient20 = DN(2,0)*u(2,1);
const double cr_det_eq_def_gradient21 = cr_det_eq_def_gradient17*cr_det_eq_def_gradient20;
const double cr_det_eq_def_gradient22 = cr_det_eq_def_gradient15*cr_det_eq_def_gradient18;
const double cr_det_eq_def_gradient23 = cr_det_eq_def_gradient13*cr_det_eq_def_gradient20;
const double cr_det_eq_def_gradient24 = 2*N[0]*th[0];
const double cr_det_eq_def_gradient25 = N[1]*th[1];
const double cr_det_eq_def_gradient26 = cr_det_eq_def_gradient24*cr_det_eq_def_gradient25;
const double cr_det_eq_def_gradient27 = N[2]*th[2];
const double cr_det_eq_def_gradient28 = cr_det_eq_def_gradient24*cr_det_eq_def_gradient27;
const double cr_det_eq_def_gradient29 = 2*cr_det_eq_def_gradient25*cr_det_eq_def_gradient27;
const double cr_det_eq_def_gradient30 = pow(N[0], 2)*pow(th[0], 2);
const double cr_det_eq_def_gradient31 = pow(N[1], 2)*pow(th[1], 2);
const double cr_det_eq_def_gradient32 = pow(N[2], 2)*pow(th[2], 2);
r_det_eq_def_gradient=1.0/(cr_det_eq_def_gradient0 + cr_det_eq_def_gradient1 + cr_det_eq_def_gradient10 + cr_det_eq_def_gradient11 - cr_det_eq_def_gradient14 - cr_det_eq_def_gradient16 - cr_det_eq_def_gradient19 + cr_det_eq_def_gradient2 - cr_det_eq_def_gradient21 - cr_det_eq_def_gradient22 - cr_det_eq_def_gradient23 + cr_det_eq_def_gradient3 + cr_det_eq_def_gradient4 + cr_det_eq_def_gradient5 + cr_det_eq_def_gradient6 + cr_det_eq_def_gradient7 + cr_det_eq_def_gradient8 + cr_det_eq_def_gradient9 + 1)*(cr_det_eq_def_gradient0*cr_det_eq_def_gradient26 + cr_det_eq_def_gradient0*cr_det_eq_def_gradient28 + cr_det_eq_def_gradient0*cr_det_eq_def_gradient29 + cr_det_eq_def_gradient0*cr_det_eq_def_gradient30 + cr_det_eq_def_gradient0*cr_det_eq_def_gradient31 + cr_det_eq_def_gradient0*cr_det_eq_def_gradient32 + cr_det_eq_def_gradient1*cr_det_eq_def_gradient26 + cr_det_eq_def_gradient1*cr_det_eq_def_gradient28 + cr_det_eq_def_gradient1*cr_det_eq_def_gradient29 + cr_det_eq_def_gradient1*cr_det_eq_def_gradient30 + cr_det_eq_def_gradient1*cr_det_eq_def_gradient31 + cr_det_eq_def_gradient1*cr_det_eq_def_gradient32 + cr_det_eq_def_gradient10*cr_det_eq_def_gradient26 + cr_det_eq_def_gradient10*cr_det_eq_def_gradient28 + cr_det_eq_def_gradient10*cr_det_eq_def_gradient29 + cr_det_eq_def_gradient10*cr_det_eq_def_gradient30 + cr_det_eq_def_gradient10*cr_det_eq_def_gradient31 + cr_det_eq_def_gradient10*cr_det_eq_def_gradient32 + cr_det_eq_def_gradient11*cr_det_eq_def_gradient26 + cr_det_eq_def_gradient11*cr_det_eq_def_gradient28 + cr_det_eq_def_gradient11*cr_det_eq_def_gradient29 + cr_det_eq_def_gradient11*cr_det_eq_def_gradient30 + cr_det_eq_def_gradient11*cr_det_eq_def_gradient31 + cr_det_eq_def_gradient11*cr_det_eq_def_gradient32 - cr_det_eq_def_gradient14*cr_det_eq_def_gradient26 - cr_det_eq_def_gradient14*cr_det_eq_def_gradient28 - cr_det_eq_def_gradient14*cr_det_eq_def_gradient29 - cr_det_eq_def_gradient14*cr_det_eq_def_gradient30 - cr_det_eq_def_gradient14*cr_det_eq_def_gradient31 - cr_det_eq_def_gradient14*cr_det_eq_def_gradient32 - cr_det_eq_def_gradient16*cr_det_eq_def_gradient26 - cr_det_eq_def_gradient16*cr_det_eq_def_gradient28 - cr_det_eq_def_gradient16*cr_det_eq_def_gradient29 - cr_det_eq_def_gradient16*cr_det_eq_def_gradient30 - cr_det_eq_def_gradient16*cr_det_eq_def_gradient31 - cr_det_eq_def_gradient16*cr_det_eq_def_gradient32 - cr_det_eq_def_gradient19*cr_det_eq_def_gradient26 - cr_det_eq_def_gradient19*cr_det_eq_def_gradient28 - cr_det_eq_def_gradient19*cr_det_eq_def_gradient29 - cr_det_eq_def_gradient19*cr_det_eq_def_gradient30 - cr_det_eq_def_gradient19*cr_det_eq_def_gradient31 - cr_det_eq_def_gradient19*cr_det_eq_def_gradient32 + cr_det_eq_def_gradient2*cr_det_eq_def_gradient26 + cr_det_eq_def_gradient2*cr_det_eq_def_gradient28 + cr_det_eq_def_gradient2*cr_det_eq_def_gradient29 + cr_det_eq_def_gradient2*cr_det_eq_def_gradient30 + cr_det_eq_def_gradient2*cr_det_eq_def_gradient31 + cr_det_eq_def_gradient2*cr_det_eq_def_gradient32 - cr_det_eq_def_gradient21*cr_det_eq_def_gradient26 - cr_det_eq_def_gradient21*cr_det_eq_def_gradient28 - cr_det_eq_def_gradient21*cr_det_eq_def_gradient29 - cr_det_eq_def_gradient21*cr_det_eq_def_gradient30 - cr_det_eq_def_gradient21*cr_det_eq_def_gradient31 - cr_det_eq_def_gradient21*cr_det_eq_def_gradient32 - cr_det_eq_def_gradient22*cr_det_eq_def_gradient26 - cr_det_eq_def_gradient22*cr_det_eq_def_gradient28 - cr_det_eq_def_gradient22*cr_det_eq_def_gradient29 - cr_det_eq_def_gradient22*cr_det_eq_def_gradient30 - cr_det_eq_def_gradient22*cr_det_eq_def_gradient31 - cr_det_eq_def_gradient22*cr_det_eq_def_gradient32 - cr_det_eq_def_gradient23*cr_det_eq_def_gradient26 - cr_det_eq_def_gradient23*cr_det_eq_def_gradient28 - cr_det_eq_def_gradient23*cr_det_eq_def_gradient29 - cr_det_eq_def_gradient23*cr_det_eq_def_gradient30 - cr_det_eq_def_gradient23*cr_det_eq_def_gradient31 - cr_det_eq_def_gradient23*cr_det_eq_def_gradient32 + cr_det_eq_def_gradient26*cr_det_eq_def_gradient3 + cr_det_eq_def_gradient26*cr_det_eq_def_gradient4 + cr_det_eq_def_gradient26*cr_det_eq_def_gradient5 + cr_det_eq_def_gradient26*cr_det_eq_def_gradient6 + cr_det_eq_def_gradient26*cr_det_eq_def_gradient7 + cr_det_eq_def_gradient26*cr_det_eq_def_gradient8 + cr_det_eq_def_gradient26*cr_det_eq_def_gradient9 + cr_det_eq_def_gradient26 + cr_det_eq_def_gradient28*cr_det_eq_def_gradient3 + cr_det_eq_def_gradient28*cr_det_eq_def_gradient4 + cr_det_eq_def_gradient28*cr_det_eq_def_gradient5 + cr_det_eq_def_gradient28*cr_det_eq_def_gradient6 + cr_det_eq_def_gradient28*cr_det_eq_def_gradient7 + cr_det_eq_def_gradient28*cr_det_eq_def_gradient8 + cr_det_eq_def_gradient28*cr_det_eq_def_gradient9 + cr_det_eq_def_gradient28 + cr_det_eq_def_gradient29*cr_det_eq_def_gradient3 + cr_det_eq_def_gradient29*cr_det_eq_def_gradient4 + cr_det_eq_def_gradient29*cr_det_eq_def_gradient5 + cr_det_eq_def_gradient29*cr_det_eq_def_gradient6 + cr_det_eq_def_gradient29*cr_det_eq_def_gradient7 + cr_det_eq_def_gradient29*cr_det_eq_def_gradient8 + cr_det_eq_def_gradient29*cr_det_eq_def_gradient9 + cr_det_eq_def_gradient29 + cr_det_eq_def_gradient3*cr_det_eq_def_gradient30 + cr_det_eq_def_gradient3*cr_det_eq_def_gradient31 + cr_det_eq_def_gradient3*cr_det_eq_def_gradient32 + cr_det_eq_def_gradient30*cr_det_eq_def_gradient4 + cr_det_eq_def_gradient30*cr_det_eq_def_gradient5 + cr_det_eq_def_gradient30*cr_det_eq_def_gradient6 + cr_det_eq_def_gradient30*cr_det_eq_def_gradient7 + cr_det_eq_def_gradient30*cr_det_eq_def_gradient8 + cr_det_eq_def_gradient30*cr_det_eq_def_gradient9 + cr_det_eq_def_gradient30 + cr_det_eq_def_gradient31*cr_det_eq_def_gradient4 + cr_det_eq_def_gradient31*cr_det_eq_def_gradient5 + cr_det_eq_def_gradient31*cr_det_eq_def_gradient6 + cr_det_eq_def_gradient31*cr_det_eq_def_gradient7 + cr_det_eq_def_gradient31*cr_det_eq_def_gradient8 + cr_det_eq_def_gradient31*cr_det_eq_def_gradient9 + cr_det_eq_def_gradient31 + cr_det_eq_def_gradient32*cr_det_eq_def_gradient4 + cr_det_eq_def_gradient32*cr_det_eq_def_gradient5 + cr_det_eq_def_gradient32*cr_det_eq_def_gradient6 + cr_det_eq_def_gradient32*cr_det_eq_def_gradient7 + cr_det_eq_def_gradient32*cr_det_eq_def_gradient8 + cr_det_eq_def_gradient32*cr_det_eq_def_gradient9 + cr_det_eq_def_gradient32);

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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DETERMINANT_F,r_node)
        KRATOS_CHECK_DOF_IN_NODE(DETERMINANT_F, r_node)
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
        // Create and initialize element variables:
        const SizeType n_nodes = r_geometry.PointsNumber();
        const SizeType dim = r_geometry.WorkingSpaceDimension();
        const SizeType strain_size = mConstitutiveLawVector[0]->GetStrainSize();

        // Create the kinematics container and fill the nodal data
        KinematicVariables kinematic_variables;
        for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_disp = r_geometry[i_node].FastGetSolutionStepValue(DISPLACEMENT);
            for (IndexType d = 0; d < dim; ++d) {
                kinematic_variables.Displacements(i_node, d) = r_disp[d];
            }
            kinematic_variables.JacobianDeterminant[i_node] = r_geometry[i_node].FastGetSolutionStepValue(DETERMINANT_F);
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
            if (rOutput[i_gauss].size() != strain_size) {
                rOutput[i_gauss].resize(strain_size, false);
            }
            rOutput[i_gauss] = constitutive_variables.StressVector;
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
            "nodal_historical"       : ["DISPLACEMENT","DETERMINANT_F"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["DISPLACEMENT","DETERMINANT_F"],
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
        std::vector<std::string> dofs_2d({"DISPLACEMENT_X","DISPLACEMENT_Y","DETERMINANT_F"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"DISPLACEMENT_X","DISPLACEMENT_Y","DISPLACEMENT_Z","DETERMINANT_F"});
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
