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
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
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
const double clhs21 = DN(0,0)*clhs6 + DN(0,0)*clhs7 + DN(0,0) - DN(0,1)*DN(1,0)*u(1,1) - DN(0,1)*DN(2,0)*u(2,1);
const double clhs22 = -clhs20*clhs21;
const double clhs23 = 0.5*clhs22;
const double clhs24 = clhs15 + 1;
const double clhs25 = pow(clhs24, 2);
const double clhs26 = clhs18*clhs4 + clhs19*clhs23 + clhs23*clhs25;
const double clhs27 = DN(0,0)*clhs14;
const double clhs28 = clhs11 + clhs12 + clhs13;
const double clhs29 = pow(clhs28, 2);
const double clhs30 = pow(clhs14, 2);
const double clhs31 = clhs18*clhs27 + clhs23*clhs29 + clhs23*clhs30;
const double clhs32 = DN(0,0)*clhs3;
const double clhs33 = DN(0,1)*clhs14;
const double clhs34 = clhs24*clhs28;
const double clhs35 = clhs14*clhs3;
const double clhs36 = clhs17*clhs32 + clhs17*clhs33 + clhs22*clhs34 + clhs22*clhs35;
const double clhs37 = C(0,0)*clhs31 + C(0,1)*clhs26 + C(0,2)*clhs36;
const double clhs38 = pow(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0, 1.0);
const double clhs39 = clhs14*clhs38;
const double clhs40 = C(0,2)*clhs31 + C(1,2)*clhs26 + C(2,2)*clhs36;
const double clhs41 = clhs3*clhs38;
const double clhs42 = clhs34 + clhs35;
const double clhs43 = C(0,2)*clhs42;
const double clhs44 = clhs17*clhs38;
const double clhs45 = 0.5*clhs29*clhs44 + 0.5*clhs30*clhs44 - 0.5;
const double clhs46 = 0.5*clhs19*clhs44 + 0.5*clhs25*clhs44 - 0.5;
const double clhs47 = C(0,0)*clhs45 + C(0,1)*clhs46 + clhs43*clhs44;
const double clhs48 = C(2,2)*clhs42;
const double clhs49 = C(0,2)*clhs45 + C(1,2)*clhs46 + clhs44*clhs48;
const double clhs50 = DN(0,0)*clhs47 + DN(0,1)*clhs49;
const double clhs51 = clhs37*clhs39 + clhs40*clhs41 + clhs50;
const double clhs52 = C(0,1)*clhs31 + C(1,1)*clhs26 + C(1,2)*clhs36;
const double clhs53 = C(1,2)*clhs42;
const double clhs54 = C(0,1)*clhs45 + C(1,1)*clhs46 + clhs44*clhs53;
const double clhs55 = DN(0,0)*clhs49 + DN(0,1)*clhs54;
const double clhs56 = clhs39*clhs40 + clhs41*clhs52 + clhs55;
const double clhs57 = DN(0,0)*clhs28;
const double clhs58 = -DN(0,0)*DN(1,1)*u(1,0) - DN(0,0)*DN(2,1)*u(2,0) + DN(0,1)*clhs10 + DN(0,1)*clhs9 + DN(0,1);
const double clhs59 = -clhs20*clhs58;
const double clhs60 = 0.5*clhs59;
const double clhs61 = clhs18*clhs57 + clhs29*clhs60 + clhs30*clhs60;
const double clhs62 = DN(0,1)*clhs24;
const double clhs63 = clhs18*clhs62 + clhs19*clhs60 + clhs25*clhs60;
const double clhs64 = DN(0,1)*clhs28;
const double clhs65 = DN(0,0)*clhs24;
const double clhs66 = clhs17*clhs64 + clhs17*clhs65 + clhs34*clhs59 + clhs35*clhs59;
const double clhs67 = C(0,2)*clhs61 + C(1,2)*clhs63 + C(2,2)*clhs66;
const double clhs68 = C(0,0)*clhs61 + C(0,1)*clhs63 + C(0,2)*clhs66;
const double clhs69 = clhs14*clhs68 + clhs3*clhs67;
const double clhs70 = C(0,1)*clhs61 + C(1,1)*clhs63 + C(1,2)*clhs66;
const double clhs71 = clhs14*clhs67 + clhs3*clhs70;
const double clhs72 = clhs29 + clhs30;
const double clhs73 = C(0,2)*clhs72;
const double clhs74 = clhs19 + clhs25;
const double clhs75 = C(1,2)*clhs74;
const double clhs76 = clhs48 + 0.5*clhs73 + 0.5*clhs75;
const double clhs77 = C(0,0)*clhs72;
const double clhs78 = C(0,1)*clhs74;
const double clhs79 = clhs43 + 0.5*clhs77 + 0.5*clhs78;
const double clhs80 = clhs14*clhs79 + clhs3*clhs76;
const double clhs81 = C(0,1)*clhs72;
const double clhs82 = C(1,1)*clhs74;
const double clhs83 = clhs53 + 0.5*clhs81 + 0.5*clhs82;
const double clhs84 = clhs14*clhs76 + clhs3*clhs83;
const double clhs85 = clhs17*(DN(0,0)*clhs80 + DN(0,1)*clhs84);
const double clhs86 = DN(1,1)*clhs3;
const double clhs87 = -DN(0,0)*DN(1,1)*u(0,1) + DN(1,0)*clhs7 + DN(1,0)*clhs8 + DN(1,0) - DN(1,1)*DN(2,0)*u(2,1);
const double clhs88 = -clhs20*clhs87;
const double clhs89 = 0.5*clhs88;
const double clhs90 = clhs18*clhs86 + clhs19*clhs89 + clhs25*clhs89;
const double clhs91 = DN(1,0)*clhs14;
const double clhs92 = clhs18*clhs91 + clhs29*clhs89 + clhs30*clhs89;
const double clhs93 = DN(1,0)*clhs3;
const double clhs94 = DN(1,1)*clhs14;
const double clhs95 = clhs17*clhs93 + clhs17*clhs94 + clhs34*clhs88 + clhs35*clhs88;
const double clhs96 = C(0,0)*clhs92 + C(0,1)*clhs90 + C(0,2)*clhs95;
const double clhs97 = C(0,2)*clhs92 + C(1,2)*clhs90 + C(2,2)*clhs95;
const double clhs98 = DN(1,0)*clhs47 + DN(1,1)*clhs49;
const double clhs99 = clhs39*clhs96 + clhs41*clhs97 + clhs98;
const double clhs100 = C(0,1)*clhs92 + C(1,1)*clhs90 + C(1,2)*clhs95;
const double clhs101 = DN(1,0)*clhs49 + DN(1,1)*clhs54;
const double clhs102 = clhs100*clhs41 + clhs101 + clhs39*clhs97;
const double clhs103 = DN(1,0)*clhs28;
const double clhs104 = -DN(0,1)*DN(1,0)*u(0,0) - DN(1,0)*DN(2,1)*u(2,0) + DN(1,1)*clhs10 + DN(1,1)*clhs5 + DN(1,1);
const double clhs105 = -clhs104*clhs20;
const double clhs106 = 0.5*clhs105;
const double clhs107 = clhs103*clhs18 + clhs106*clhs29 + clhs106*clhs30;
const double clhs108 = DN(1,1)*clhs24;
const double clhs109 = clhs106*clhs19 + clhs106*clhs25 + clhs108*clhs18;
const double clhs110 = DN(1,1)*clhs28;
const double clhs111 = DN(1,0)*clhs24;
const double clhs112 = clhs105*clhs34 + clhs105*clhs35 + clhs110*clhs17 + clhs111*clhs17;
const double clhs113 = C(0,2)*clhs107 + C(1,2)*clhs109 + C(2,2)*clhs112;
const double clhs114 = C(0,0)*clhs107 + C(0,1)*clhs109 + C(0,2)*clhs112;
const double clhs115 = clhs113*clhs3 + clhs114*clhs14;
const double clhs116 = C(0,1)*clhs107 + C(1,1)*clhs109 + C(1,2)*clhs112;
const double clhs117 = clhs113*clhs14 + clhs116*clhs3;
const double clhs118 = DN(2,1)*clhs3;
const double clhs119 = -DN(0,0)*DN(2,1)*u(0,1) - DN(1,0)*DN(2,1)*u(1,1) + DN(2,0)*clhs6 + DN(2,0)*clhs8 + DN(2,0);
const double clhs120 = -clhs119*clhs20;
const double clhs121 = 0.5*clhs120;
const double clhs122 = clhs118*clhs18 + clhs121*clhs19 + clhs121*clhs25;
const double clhs123 = DN(2,0)*clhs14;
const double clhs124 = clhs121*clhs29 + clhs121*clhs30 + clhs123*clhs18;
const double clhs125 = DN(2,0)*clhs3;
const double clhs126 = DN(2,1)*clhs14;
const double clhs127 = clhs120*clhs34 + clhs120*clhs35 + clhs125*clhs17 + clhs126*clhs17;
const double clhs128 = C(0,0)*clhs124 + C(0,1)*clhs122 + C(0,2)*clhs127;
const double clhs129 = C(0,2)*clhs124 + C(1,2)*clhs122 + C(2,2)*clhs127;
const double clhs130 = DN(2,0)*clhs47 + DN(2,1)*clhs49;
const double clhs131 = clhs128*clhs39 + clhs129*clhs41 + clhs130;
const double clhs132 = C(0,1)*clhs124 + C(1,1)*clhs122 + C(1,2)*clhs127;
const double clhs133 = DN(2,0)*clhs49 + DN(2,1)*clhs54;
const double clhs134 = clhs129*clhs39 + clhs132*clhs41 + clhs133;
const double clhs135 = DN(2,0)*clhs28;
const double clhs136 = -DN(0,1)*DN(2,0)*u(0,0) - DN(1,1)*DN(2,0)*u(1,0) + DN(2,1)*clhs5 + DN(2,1)*clhs9 + DN(2,1);
const double clhs137 = -clhs136*clhs20;
const double clhs138 = 0.5*clhs137;
const double clhs139 = clhs135*clhs18 + clhs138*clhs29 + clhs138*clhs30;
const double clhs140 = DN(2,1)*clhs24;
const double clhs141 = clhs138*clhs19 + clhs138*clhs25 + clhs140*clhs18;
const double clhs142 = DN(2,1)*clhs28;
const double clhs143 = DN(2,0)*clhs24;
const double clhs144 = clhs137*clhs34 + clhs137*clhs35 + clhs142*clhs17 + clhs143*clhs17;
const double clhs145 = C(0,2)*clhs139 + C(1,2)*clhs141 + C(2,2)*clhs144;
const double clhs146 = C(0,0)*clhs139 + C(0,1)*clhs141 + C(0,2)*clhs144;
const double clhs147 = clhs14*clhs146 + clhs145*clhs3;
const double clhs148 = C(0,1)*clhs139 + C(1,1)*clhs141 + C(1,2)*clhs144;
const double clhs149 = clhs14*clhs145 + clhs148*clhs3;
const double clhs150 = clhs24*clhs40 + clhs28*clhs37;
const double clhs151 = clhs24*clhs52 + clhs28*clhs40;
const double clhs152 = clhs28*clhs38;
const double clhs153 = clhs24*clhs38;
const double clhs154 = clhs152*clhs68 + clhs153*clhs67 + clhs50;
const double clhs155 = clhs152*clhs67 + clhs153*clhs70 + clhs55;
const double clhs156 = clhs24*clhs76 + clhs28*clhs79;
const double clhs157 = clhs24*clhs83 + clhs28*clhs76;
const double clhs158 = clhs17*(DN(0,0)*clhs156 + DN(0,1)*clhs157);
const double clhs159 = clhs24*clhs97 + clhs28*clhs96;
const double clhs160 = clhs100*clhs24 + clhs28*clhs97;
const double clhs161 = clhs113*clhs153 + clhs114*clhs152 + clhs98;
const double clhs162 = clhs101 + clhs113*clhs152 + clhs116*clhs153;
const double clhs163 = clhs128*clhs28 + clhs129*clhs24;
const double clhs164 = clhs129*clhs28 + clhs132*clhs24;
const double clhs165 = clhs130 + clhs145*clhs153 + clhs146*clhs152;
const double clhs166 = clhs133 + clhs145*clhs152 + clhs148*clhs153;
const double clhs167 = clhs32 + clhs33;
const double clhs168 = C(0,0)*clhs27 + C(0,1)*clhs4 + C(0,2)*clhs167;
const double clhs169 = C(0,2)*clhs27 + C(1,2)*clhs4 + C(2,2)*clhs167;
const double clhs170 = tau*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double clhs171 = C(0,1)*clhs27 + C(1,1)*clhs4 + C(1,2)*clhs167;
const double clhs172 = tau*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double clhs173 = clhs64 + clhs65;
const double clhs174 = C(0,0)*clhs57 + C(0,1)*clhs62 + C(0,2)*clhs173;
const double clhs175 = C(0,2)*clhs57 + C(1,2)*clhs62 + C(2,2)*clhs173;
const double clhs176 = C(0,1)*clhs57 + C(1,1)*clhs62 + C(1,2)*clhs173;
const double clhs177 = 2*clhs43 + clhs77 + clhs78;
const double clhs178 = 2*clhs48 + clhs73 + clhs75;
const double clhs179 = (1.0/2.0)*tau;
const double clhs180 = clhs179*(DN(0,0)*clhs177 + DN(0,1)*clhs178);
const double clhs181 = 2*clhs53 + clhs81 + clhs82;
const double clhs182 = clhs179*(DN(0,0)*clhs178 + DN(0,1)*clhs181);
const double clhs183 = tau*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs184 = b_gauss[1]*clhs183;
const double clhs185 = clhs93 + clhs94;
const double clhs186 = C(0,0)*clhs91 + C(0,1)*clhs86 + C(0,2)*clhs185;
const double clhs187 = C(0,2)*clhs91 + C(1,2)*clhs86 + C(2,2)*clhs185;
const double clhs188 = C(0,1)*clhs91 + C(1,1)*clhs86 + C(1,2)*clhs185;
const double clhs189 = b_gauss[0]*clhs183;
const double clhs190 = clhs110 + clhs111;
const double clhs191 = C(0,0)*clhs103 + C(0,1)*clhs108 + C(0,2)*clhs190;
const double clhs192 = C(0,2)*clhs103 + C(1,2)*clhs108 + C(2,2)*clhs190;
const double clhs193 = C(0,1)*clhs103 + C(1,1)*clhs108 + C(1,2)*clhs190;
const double clhs194 = N[0]*N[1];
const double clhs195 = tau*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs196 = b_gauss[1]*clhs195;
const double clhs197 = clhs125 + clhs126;
const double clhs198 = C(0,0)*clhs123 + C(0,1)*clhs118 + C(0,2)*clhs197;
const double clhs199 = C(0,2)*clhs123 + C(1,2)*clhs118 + C(2,2)*clhs197;
const double clhs200 = C(0,1)*clhs123 + C(1,1)*clhs118 + C(1,2)*clhs197;
const double clhs201 = b_gauss[0]*clhs195;
const double clhs202 = clhs142 + clhs143;
const double clhs203 = C(0,0)*clhs135 + C(0,1)*clhs140 + C(0,2)*clhs202;
const double clhs204 = C(0,2)*clhs135 + C(1,2)*clhs140 + C(2,2)*clhs202;
const double clhs205 = C(0,1)*clhs135 + C(1,1)*clhs140 + C(1,2)*clhs202;
const double clhs206 = N[0]*N[2];
const double clhs207 = clhs17*(DN(1,0)*clhs80 + DN(1,1)*clhs84);
const double clhs208 = clhs17*(DN(1,0)*clhs156 + DN(1,1)*clhs157);
const double clhs209 = clhs179*(DN(1,0)*clhs177 + DN(1,1)*clhs178);
const double clhs210 = clhs179*(DN(1,0)*clhs178 + DN(1,1)*clhs181);
const double clhs211 = tau*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs212 = b_gauss[1]*clhs211;
const double clhs213 = b_gauss[0]*clhs211;
const double clhs214 = N[1]*N[2];
const double clhs215 = clhs17*(DN(2,0)*clhs80 + DN(2,1)*clhs84);
const double clhs216 = clhs17*(DN(2,0)*clhs156 + DN(2,1)*clhs157);
const double clhs217 = clhs179*(DN(2,0)*clhs177 + DN(2,1)*clhs178);
const double clhs218 = clhs179*(DN(2,0)*clhs178 + DN(2,1)*clhs181);
lhs(0,0)=DN(0,0)*clhs51 + DN(0,1)*clhs56;
lhs(0,1)=clhs38*(DN(0,0)*clhs69 + DN(0,1)*clhs71);
lhs(0,2)=N[0]*clhs85;
lhs(0,3)=DN(0,0)*clhs99 + DN(0,1)*clhs102;
lhs(0,4)=clhs38*(DN(0,0)*clhs115 + DN(0,1)*clhs117);
lhs(0,5)=N[1]*clhs85;
lhs(0,6)=DN(0,0)*clhs131 + DN(0,1)*clhs134;
lhs(0,7)=clhs38*(DN(0,0)*clhs147 + DN(0,1)*clhs149);
lhs(0,8)=N[2]*clhs85;
lhs(1,0)=clhs38*(DN(0,0)*clhs150 + DN(0,1)*clhs151);
lhs(1,1)=DN(0,0)*clhs154 + DN(0,1)*clhs155;
lhs(1,2)=N[0]*clhs158;
lhs(1,3)=clhs38*(DN(0,0)*clhs159 + DN(0,1)*clhs160);
lhs(1,4)=DN(0,0)*clhs161 + DN(0,1)*clhs162;
lhs(1,5)=N[1]*clhs158;
lhs(1,6)=clhs38*(DN(0,0)*clhs163 + DN(0,1)*clhs164);
lhs(1,7)=DN(0,0)*clhs165 + DN(0,1)*clhs166;
lhs(1,8)=N[2]*clhs158;
lhs(2,0)=N[0]*clhs21 - clhs170*(DN(0,0)*clhs168 + DN(0,1)*clhs169) - clhs172*(DN(0,0)*clhs169 + DN(0,1)*clhs171);
lhs(2,1)=N[0]*clhs58 - clhs170*(DN(0,0)*clhs174 + DN(0,1)*clhs175) - clhs172*(DN(0,0)*clhs175 + DN(0,1)*clhs176);
lhs(2,2)=-DN(0,0)*clhs180 - DN(0,1)*clhs182 - pow(N[0], 2);
lhs(2,3)=N[0]*clhs87 - clhs170*(DN(0,0)*clhs186 + DN(0,1)*clhs187) - clhs172*(DN(0,0)*clhs187 + DN(0,1)*clhs188) + clhs184;
lhs(2,4)=N[0]*clhs104 - clhs170*(DN(0,0)*clhs191 + DN(0,1)*clhs192) - clhs172*(DN(0,0)*clhs192 + DN(0,1)*clhs193) - clhs189;
lhs(2,5)=-DN(1,0)*clhs180 - DN(1,1)*clhs182 - clhs194;
lhs(2,6)=N[0]*clhs119 - clhs170*(DN(0,0)*clhs198 + DN(0,1)*clhs199) - clhs172*(DN(0,0)*clhs199 + DN(0,1)*clhs200) + clhs196;
lhs(2,7)=N[0]*clhs136 - clhs170*(DN(0,0)*clhs203 + DN(0,1)*clhs204) - clhs172*(DN(0,0)*clhs204 + DN(0,1)*clhs205) - clhs201;
lhs(2,8)=-DN(2,0)*clhs180 - DN(2,1)*clhs182 - clhs206;
lhs(3,0)=DN(1,0)*clhs51 + DN(1,1)*clhs56;
lhs(3,1)=clhs38*(DN(1,0)*clhs69 + DN(1,1)*clhs71);
lhs(3,2)=N[0]*clhs207;
lhs(3,3)=DN(1,0)*clhs99 + DN(1,1)*clhs102;
lhs(3,4)=clhs38*(DN(1,0)*clhs115 + DN(1,1)*clhs117);
lhs(3,5)=N[1]*clhs207;
lhs(3,6)=DN(1,0)*clhs131 + DN(1,1)*clhs134;
lhs(3,7)=clhs38*(DN(1,0)*clhs147 + DN(1,1)*clhs149);
lhs(3,8)=N[2]*clhs207;
lhs(4,0)=clhs38*(DN(1,0)*clhs150 + DN(1,1)*clhs151);
lhs(4,1)=DN(1,0)*clhs154 + DN(1,1)*clhs155;
lhs(4,2)=N[0]*clhs208;
lhs(4,3)=clhs38*(DN(1,0)*clhs159 + DN(1,1)*clhs160);
lhs(4,4)=DN(1,0)*clhs161 + DN(1,1)*clhs162;
lhs(4,5)=N[1]*clhs208;
lhs(4,6)=clhs38*(DN(1,0)*clhs163 + DN(1,1)*clhs164);
lhs(4,7)=DN(1,0)*clhs165 + DN(1,1)*clhs166;
lhs(4,8)=N[2]*clhs208;
lhs(5,0)=N[1]*clhs21 - clhs170*(DN(1,0)*clhs168 + DN(1,1)*clhs169) - clhs172*(DN(1,0)*clhs169 + DN(1,1)*clhs171) - clhs184;
lhs(5,1)=N[1]*clhs58 - clhs170*(DN(1,0)*clhs174 + DN(1,1)*clhs175) - clhs172*(DN(1,0)*clhs175 + DN(1,1)*clhs176) + clhs189;
lhs(5,2)=-DN(0,0)*clhs209 - DN(0,1)*clhs210 - clhs194;
lhs(5,3)=N[1]*clhs87 - clhs170*(DN(1,0)*clhs186 + DN(1,1)*clhs187) - clhs172*(DN(1,0)*clhs187 + DN(1,1)*clhs188);
lhs(5,4)=N[1]*clhs104 - clhs170*(DN(1,0)*clhs191 + DN(1,1)*clhs192) - clhs172*(DN(1,0)*clhs192 + DN(1,1)*clhs193);
lhs(5,5)=-DN(1,0)*clhs209 - DN(1,1)*clhs210 - pow(N[1], 2);
lhs(5,6)=N[1]*clhs119 - clhs170*(DN(1,0)*clhs198 + DN(1,1)*clhs199) - clhs172*(DN(1,0)*clhs199 + DN(1,1)*clhs200) + clhs212;
lhs(5,7)=N[1]*clhs136 - clhs170*(DN(1,0)*clhs203 + DN(1,1)*clhs204) - clhs172*(DN(1,0)*clhs204 + DN(1,1)*clhs205) - clhs213;
lhs(5,8)=-DN(2,0)*clhs209 - DN(2,1)*clhs210 - clhs214;
lhs(6,0)=DN(2,0)*clhs51 + DN(2,1)*clhs56;
lhs(6,1)=clhs38*(DN(2,0)*clhs69 + DN(2,1)*clhs71);
lhs(6,2)=N[0]*clhs215;
lhs(6,3)=DN(2,0)*clhs99 + DN(2,1)*clhs102;
lhs(6,4)=clhs38*(DN(2,0)*clhs115 + DN(2,1)*clhs117);
lhs(6,5)=N[1]*clhs215;
lhs(6,6)=DN(2,0)*clhs131 + DN(2,1)*clhs134;
lhs(6,7)=clhs38*(DN(2,0)*clhs147 + DN(2,1)*clhs149);
lhs(6,8)=N[2]*clhs215;
lhs(7,0)=clhs38*(DN(2,0)*clhs150 + DN(2,1)*clhs151);
lhs(7,1)=DN(2,0)*clhs154 + DN(2,1)*clhs155;
lhs(7,2)=N[0]*clhs216;
lhs(7,3)=clhs38*(DN(2,0)*clhs159 + DN(2,1)*clhs160);
lhs(7,4)=DN(2,0)*clhs161 + DN(2,1)*clhs162;
lhs(7,5)=N[1]*clhs216;
lhs(7,6)=clhs38*(DN(2,0)*clhs163 + DN(2,1)*clhs164);
lhs(7,7)=DN(2,0)*clhs165 + DN(2,1)*clhs166;
lhs(7,8)=N[2]*clhs216;
lhs(8,0)=N[2]*clhs21 - clhs170*(DN(2,0)*clhs168 + DN(2,1)*clhs169) - clhs172*(DN(2,0)*clhs169 + DN(2,1)*clhs171) - clhs196;
lhs(8,1)=N[2]*clhs58 - clhs170*(DN(2,0)*clhs174 + DN(2,1)*clhs175) - clhs172*(DN(2,0)*clhs175 + DN(2,1)*clhs176) + clhs201;
lhs(8,2)=-DN(0,0)*clhs217 - DN(0,1)*clhs218 - clhs206;
lhs(8,3)=N[2]*clhs87 - clhs170*(DN(2,0)*clhs186 + DN(2,1)*clhs187) - clhs172*(DN(2,0)*clhs187 + DN(2,1)*clhs188) - clhs212;
lhs(8,4)=N[2]*clhs104 - clhs170*(DN(2,0)*clhs191 + DN(2,1)*clhs192) - clhs172*(DN(2,0)*clhs192 + DN(2,1)*clhs193) + clhs213;
lhs(8,5)=-DN(1,0)*clhs217 - DN(1,1)*clhs218 - clhs214;
lhs(8,6)=N[2]*clhs119 - clhs170*(DN(2,0)*clhs198 + DN(2,1)*clhs199) - clhs172*(DN(2,0)*clhs199 + DN(2,1)*clhs200);
lhs(8,7)=N[2]*clhs136 - clhs170*(DN(2,0)*clhs203 + DN(2,1)*clhs204) - clhs172*(DN(2,0)*clhs204 + DN(2,1)*clhs205);
lhs(8,8)=-DN(2,0)*clhs217 - DN(2,1)*clhs218 - pow(N[2], 2);

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
const double crhs15 = DN(0,1)*u(0,1);
const double crhs16 = DN(1,1)*u(1,1);
const double crhs17 = DN(2,1)*u(2,1);
const double crhs18 = crhs15 + crhs16 + crhs17;
const double crhs19 = crhs18 + 1;
const double crhs20 = S[0]*crhs14 + S[2]*crhs19;
const double crhs21 = S[1]*crhs19 + S[2]*crhs14;
const double crhs22 = b_gauss[0]*tau;
const double crhs23 = b_gauss[1]*tau;
const double crhs24 = -N[0]*th[0] - N[1]*th[1] - N[2]*th[2] - crhs0*crhs12 - crhs0*crhs13 - crhs1*crhs11 - crhs1*crhs13 - crhs11*crhs2 - crhs12*crhs2 + crhs15*crhs5 + crhs15*crhs6 + crhs16*crhs4 + crhs16*crhs6 + crhs17*crhs4 + crhs17*crhs5 + crhs18 + crhs7;
const double crhs25 = pow(crhs14, 2) + pow(crhs8, 2);
const double crhs26 = pow(crhs19, 2) + pow(crhs3, 2);
const double crhs27 = 2*crhs14*crhs19 + 2*crhs3*crhs8;
const double crhs28 = C(0,0)*crhs25 + C(0,1)*crhs26 + C(0,2)*crhs27;
const double crhs29 = C(0,2)*crhs25 + C(1,2)*crhs26 + C(2,2)*crhs27;
const double crhs30 = (1.0/2.0)*tau;
const double crhs31 = crhs30*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs32 = C(0,1)*crhs25 + C(1,1)*crhs26 + C(1,2)*crhs27;
const double crhs33 = crhs30*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
rhs[0]=-DN(0,0)*crhs9 - DN(0,1)*crhs10 + N[0]*b_gauss[0];
rhs[1]=-DN(0,0)*crhs20 - DN(0,1)*crhs21 + N[0]*b_gauss[1];
rhs[2]=-N[0]*crhs24 + crhs22*(DN(0,0)*crhs19 - DN(0,1)*crhs14) - crhs23*(DN(0,0)*crhs3 - DN(0,1)*crhs8) + crhs31*(DN(0,0)*crhs28 + DN(0,1)*crhs29) + crhs33*(DN(0,0)*crhs29 + DN(0,1)*crhs32);
rhs[3]=-DN(1,0)*crhs9 - DN(1,1)*crhs10 + N[1]*b_gauss[0];
rhs[4]=-DN(1,0)*crhs20 - DN(1,1)*crhs21 + N[1]*b_gauss[1];
rhs[5]=-N[1]*crhs24 + crhs22*(DN(1,0)*crhs19 - DN(1,1)*crhs14) - crhs23*(DN(1,0)*crhs3 - DN(1,1)*crhs8) + crhs31*(DN(1,0)*crhs28 + DN(1,1)*crhs29) + crhs33*(DN(1,0)*crhs29 + DN(1,1)*crhs32);
rhs[6]=-DN(2,0)*crhs9 - DN(2,1)*crhs10 + N[2]*b_gauss[0];
rhs[7]=-DN(2,0)*crhs20 - DN(2,1)*crhs21 + N[2]*b_gauss[1];
rhs[8]=-N[2]*crhs24 + crhs22*(DN(2,0)*crhs19 - DN(2,1)*crhs14) - crhs23*(DN(2,0)*crhs3 - DN(2,1)*crhs8) + crhs31*(DN(2,0)*crhs28 + DN(2,1)*crhs29) + crhs33*(DN(2,0)*crhs29 + DN(2,1)*crhs32);

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
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
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
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
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
const double clhs21 = DN(0,0)*clhs6 + DN(0,0)*clhs7 + DN(0,0) - DN(0,1)*DN(1,0)*u(1,1) - DN(0,1)*DN(2,0)*u(2,1);
const double clhs22 = -clhs20*clhs21;
const double clhs23 = 0.5*clhs22;
const double clhs24 = clhs15 + 1;
const double clhs25 = pow(clhs24, 2);
const double clhs26 = clhs18*clhs4 + clhs19*clhs23 + clhs23*clhs25;
const double clhs27 = DN(0,0)*clhs14;
const double clhs28 = clhs11 + clhs12 + clhs13;
const double clhs29 = pow(clhs28, 2);
const double clhs30 = pow(clhs14, 2);
const double clhs31 = clhs18*clhs27 + clhs23*clhs29 + clhs23*clhs30;
const double clhs32 = DN(0,0)*clhs3;
const double clhs33 = DN(0,1)*clhs14;
const double clhs34 = clhs24*clhs28;
const double clhs35 = clhs14*clhs3;
const double clhs36 = clhs17*clhs32 + clhs17*clhs33 + clhs22*clhs34 + clhs22*clhs35;
const double clhs37 = C(0,0)*clhs31 + C(0,1)*clhs26 + C(0,2)*clhs36;
const double clhs38 = pow(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0, 1.0);
const double clhs39 = clhs14*clhs38;
const double clhs40 = C(0,2)*clhs31 + C(1,2)*clhs26 + C(2,2)*clhs36;
const double clhs41 = clhs3*clhs38;
const double clhs42 = clhs34 + clhs35;
const double clhs43 = C(0,2)*clhs42;
const double clhs44 = clhs17*clhs38;
const double clhs45 = 0.5*clhs29*clhs44 + 0.5*clhs30*clhs44 - 0.5;
const double clhs46 = 0.5*clhs19*clhs44 + 0.5*clhs25*clhs44 - 0.5;
const double clhs47 = C(0,0)*clhs45 + C(0,1)*clhs46 + clhs43*clhs44;
const double clhs48 = C(2,2)*clhs42;
const double clhs49 = C(0,2)*clhs45 + C(1,2)*clhs46 + clhs44*clhs48;
const double clhs50 = DN(0,0)*clhs47 + DN(0,1)*clhs49;
const double clhs51 = clhs37*clhs39 + clhs40*clhs41 + clhs50;
const double clhs52 = C(0,1)*clhs31 + C(1,1)*clhs26 + C(1,2)*clhs36;
const double clhs53 = C(1,2)*clhs42;
const double clhs54 = C(0,1)*clhs45 + C(1,1)*clhs46 + clhs44*clhs53;
const double clhs55 = DN(0,0)*clhs49 + DN(0,1)*clhs54;
const double clhs56 = clhs39*clhs40 + clhs41*clhs52 + clhs55;
const double clhs57 = DN(0,0)*clhs28;
const double clhs58 = -DN(0,0)*DN(1,1)*u(1,0) - DN(0,0)*DN(2,1)*u(2,0) + DN(0,1)*clhs10 + DN(0,1)*clhs9 + DN(0,1);
const double clhs59 = -clhs20*clhs58;
const double clhs60 = 0.5*clhs59;
const double clhs61 = clhs18*clhs57 + clhs29*clhs60 + clhs30*clhs60;
const double clhs62 = DN(0,1)*clhs24;
const double clhs63 = clhs18*clhs62 + clhs19*clhs60 + clhs25*clhs60;
const double clhs64 = DN(0,1)*clhs28;
const double clhs65 = DN(0,0)*clhs24;
const double clhs66 = clhs17*clhs64 + clhs17*clhs65 + clhs34*clhs59 + clhs35*clhs59;
const double clhs67 = C(0,2)*clhs61 + C(1,2)*clhs63 + C(2,2)*clhs66;
const double clhs68 = C(0,0)*clhs61 + C(0,1)*clhs63 + C(0,2)*clhs66;
const double clhs69 = clhs14*clhs68 + clhs3*clhs67;
const double clhs70 = C(0,1)*clhs61 + C(1,1)*clhs63 + C(1,2)*clhs66;
const double clhs71 = clhs14*clhs67 + clhs3*clhs70;
const double clhs72 = clhs29 + clhs30;
const double clhs73 = C(0,2)*clhs72;
const double clhs74 = clhs19 + clhs25;
const double clhs75 = C(1,2)*clhs74;
const double clhs76 = clhs48 + 0.5*clhs73 + 0.5*clhs75;
const double clhs77 = C(0,0)*clhs72;
const double clhs78 = C(0,1)*clhs74;
const double clhs79 = clhs43 + 0.5*clhs77 + 0.5*clhs78;
const double clhs80 = clhs14*clhs79 + clhs3*clhs76;
const double clhs81 = C(0,1)*clhs72;
const double clhs82 = C(1,1)*clhs74;
const double clhs83 = clhs53 + 0.5*clhs81 + 0.5*clhs82;
const double clhs84 = clhs14*clhs76 + clhs3*clhs83;
const double clhs85 = clhs17*(DN(0,0)*clhs80 + DN(0,1)*clhs84);
const double clhs86 = DN(1,1)*clhs3;
const double clhs87 = -DN(0,0)*DN(1,1)*u(0,1) + DN(1,0)*clhs7 + DN(1,0)*clhs8 + DN(1,0) - DN(1,1)*DN(2,0)*u(2,1);
const double clhs88 = -clhs20*clhs87;
const double clhs89 = 0.5*clhs88;
const double clhs90 = clhs18*clhs86 + clhs19*clhs89 + clhs25*clhs89;
const double clhs91 = DN(1,0)*clhs14;
const double clhs92 = clhs18*clhs91 + clhs29*clhs89 + clhs30*clhs89;
const double clhs93 = DN(1,0)*clhs3;
const double clhs94 = DN(1,1)*clhs14;
const double clhs95 = clhs17*clhs93 + clhs17*clhs94 + clhs34*clhs88 + clhs35*clhs88;
const double clhs96 = C(0,0)*clhs92 + C(0,1)*clhs90 + C(0,2)*clhs95;
const double clhs97 = C(0,2)*clhs92 + C(1,2)*clhs90 + C(2,2)*clhs95;
const double clhs98 = DN(1,0)*clhs47 + DN(1,1)*clhs49;
const double clhs99 = clhs39*clhs96 + clhs41*clhs97 + clhs98;
const double clhs100 = C(0,1)*clhs92 + C(1,1)*clhs90 + C(1,2)*clhs95;
const double clhs101 = DN(1,0)*clhs49 + DN(1,1)*clhs54;
const double clhs102 = clhs100*clhs41 + clhs101 + clhs39*clhs97;
const double clhs103 = DN(1,0)*clhs28;
const double clhs104 = -DN(0,1)*DN(1,0)*u(0,0) - DN(1,0)*DN(2,1)*u(2,0) + DN(1,1)*clhs10 + DN(1,1)*clhs5 + DN(1,1);
const double clhs105 = -clhs104*clhs20;
const double clhs106 = 0.5*clhs105;
const double clhs107 = clhs103*clhs18 + clhs106*clhs29 + clhs106*clhs30;
const double clhs108 = DN(1,1)*clhs24;
const double clhs109 = clhs106*clhs19 + clhs106*clhs25 + clhs108*clhs18;
const double clhs110 = DN(1,1)*clhs28;
const double clhs111 = DN(1,0)*clhs24;
const double clhs112 = clhs105*clhs34 + clhs105*clhs35 + clhs110*clhs17 + clhs111*clhs17;
const double clhs113 = C(0,2)*clhs107 + C(1,2)*clhs109 + C(2,2)*clhs112;
const double clhs114 = C(0,0)*clhs107 + C(0,1)*clhs109 + C(0,2)*clhs112;
const double clhs115 = clhs113*clhs3 + clhs114*clhs14;
const double clhs116 = C(0,1)*clhs107 + C(1,1)*clhs109 + C(1,2)*clhs112;
const double clhs117 = clhs113*clhs14 + clhs116*clhs3;
const double clhs118 = DN(2,1)*clhs3;
const double clhs119 = -DN(0,0)*DN(2,1)*u(0,1) - DN(1,0)*DN(2,1)*u(1,1) + DN(2,0)*clhs6 + DN(2,0)*clhs8 + DN(2,0);
const double clhs120 = -clhs119*clhs20;
const double clhs121 = 0.5*clhs120;
const double clhs122 = clhs118*clhs18 + clhs121*clhs19 + clhs121*clhs25;
const double clhs123 = DN(2,0)*clhs14;
const double clhs124 = clhs121*clhs29 + clhs121*clhs30 + clhs123*clhs18;
const double clhs125 = DN(2,0)*clhs3;
const double clhs126 = DN(2,1)*clhs14;
const double clhs127 = clhs120*clhs34 + clhs120*clhs35 + clhs125*clhs17 + clhs126*clhs17;
const double clhs128 = C(0,0)*clhs124 + C(0,1)*clhs122 + C(0,2)*clhs127;
const double clhs129 = C(0,2)*clhs124 + C(1,2)*clhs122 + C(2,2)*clhs127;
const double clhs130 = DN(2,0)*clhs47 + DN(2,1)*clhs49;
const double clhs131 = clhs128*clhs39 + clhs129*clhs41 + clhs130;
const double clhs132 = C(0,1)*clhs124 + C(1,1)*clhs122 + C(1,2)*clhs127;
const double clhs133 = DN(2,0)*clhs49 + DN(2,1)*clhs54;
const double clhs134 = clhs129*clhs39 + clhs132*clhs41 + clhs133;
const double clhs135 = DN(2,0)*clhs28;
const double clhs136 = -DN(0,1)*DN(2,0)*u(0,0) - DN(1,1)*DN(2,0)*u(1,0) + DN(2,1)*clhs5 + DN(2,1)*clhs9 + DN(2,1);
const double clhs137 = -clhs136*clhs20;
const double clhs138 = 0.5*clhs137;
const double clhs139 = clhs135*clhs18 + clhs138*clhs29 + clhs138*clhs30;
const double clhs140 = DN(2,1)*clhs24;
const double clhs141 = clhs138*clhs19 + clhs138*clhs25 + clhs140*clhs18;
const double clhs142 = DN(2,1)*clhs28;
const double clhs143 = DN(2,0)*clhs24;
const double clhs144 = clhs137*clhs34 + clhs137*clhs35 + clhs142*clhs17 + clhs143*clhs17;
const double clhs145 = C(0,2)*clhs139 + C(1,2)*clhs141 + C(2,2)*clhs144;
const double clhs146 = C(0,0)*clhs139 + C(0,1)*clhs141 + C(0,2)*clhs144;
const double clhs147 = clhs14*clhs146 + clhs145*clhs3;
const double clhs148 = C(0,1)*clhs139 + C(1,1)*clhs141 + C(1,2)*clhs144;
const double clhs149 = clhs14*clhs145 + clhs148*clhs3;
const double clhs150 = clhs24*clhs40 + clhs28*clhs37;
const double clhs151 = clhs24*clhs52 + clhs28*clhs40;
const double clhs152 = clhs28*clhs38;
const double clhs153 = clhs24*clhs38;
const double clhs154 = clhs152*clhs68 + clhs153*clhs67 + clhs50;
const double clhs155 = clhs152*clhs67 + clhs153*clhs70 + clhs55;
const double clhs156 = clhs24*clhs76 + clhs28*clhs79;
const double clhs157 = clhs24*clhs83 + clhs28*clhs76;
const double clhs158 = clhs17*(DN(0,0)*clhs156 + DN(0,1)*clhs157);
const double clhs159 = clhs24*clhs97 + clhs28*clhs96;
const double clhs160 = clhs100*clhs24 + clhs28*clhs97;
const double clhs161 = clhs113*clhs153 + clhs114*clhs152 + clhs98;
const double clhs162 = clhs101 + clhs113*clhs152 + clhs116*clhs153;
const double clhs163 = clhs128*clhs28 + clhs129*clhs24;
const double clhs164 = clhs129*clhs28 + clhs132*clhs24;
const double clhs165 = clhs130 + clhs145*clhs153 + clhs146*clhs152;
const double clhs166 = clhs133 + clhs145*clhs152 + clhs148*clhs153;
const double clhs167 = clhs32 + clhs33;
const double clhs168 = C(0,0)*clhs27 + C(0,1)*clhs4 + C(0,2)*clhs167;
const double clhs169 = C(0,2)*clhs27 + C(1,2)*clhs4 + C(2,2)*clhs167;
const double clhs170 = tau*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double clhs171 = C(0,1)*clhs27 + C(1,1)*clhs4 + C(1,2)*clhs167;
const double clhs172 = tau*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double clhs173 = clhs64 + clhs65;
const double clhs174 = C(0,0)*clhs57 + C(0,1)*clhs62 + C(0,2)*clhs173;
const double clhs175 = C(0,2)*clhs57 + C(1,2)*clhs62 + C(2,2)*clhs173;
const double clhs176 = C(0,1)*clhs57 + C(1,1)*clhs62 + C(1,2)*clhs173;
const double clhs177 = 2*clhs43 + clhs77 + clhs78;
const double clhs178 = 2*clhs48 + clhs73 + clhs75;
const double clhs179 = (1.0/2.0)*tau;
const double clhs180 = clhs179*(DN(0,0)*clhs177 + DN(0,1)*clhs178);
const double clhs181 = 2*clhs53 + clhs81 + clhs82;
const double clhs182 = clhs179*(DN(0,0)*clhs178 + DN(0,1)*clhs181);
const double clhs183 = tau*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs184 = b_gauss[1]*clhs183;
const double clhs185 = clhs93 + clhs94;
const double clhs186 = C(0,0)*clhs91 + C(0,1)*clhs86 + C(0,2)*clhs185;
const double clhs187 = C(0,2)*clhs91 + C(1,2)*clhs86 + C(2,2)*clhs185;
const double clhs188 = C(0,1)*clhs91 + C(1,1)*clhs86 + C(1,2)*clhs185;
const double clhs189 = b_gauss[0]*clhs183;
const double clhs190 = clhs110 + clhs111;
const double clhs191 = C(0,0)*clhs103 + C(0,1)*clhs108 + C(0,2)*clhs190;
const double clhs192 = C(0,2)*clhs103 + C(1,2)*clhs108 + C(2,2)*clhs190;
const double clhs193 = C(0,1)*clhs103 + C(1,1)*clhs108 + C(1,2)*clhs190;
const double clhs194 = N[0]*N[1];
const double clhs195 = tau*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs196 = b_gauss[1]*clhs195;
const double clhs197 = clhs125 + clhs126;
const double clhs198 = C(0,0)*clhs123 + C(0,1)*clhs118 + C(0,2)*clhs197;
const double clhs199 = C(0,2)*clhs123 + C(1,2)*clhs118 + C(2,2)*clhs197;
const double clhs200 = C(0,1)*clhs123 + C(1,1)*clhs118 + C(1,2)*clhs197;
const double clhs201 = b_gauss[0]*clhs195;
const double clhs202 = clhs142 + clhs143;
const double clhs203 = C(0,0)*clhs135 + C(0,1)*clhs140 + C(0,2)*clhs202;
const double clhs204 = C(0,2)*clhs135 + C(1,2)*clhs140 + C(2,2)*clhs202;
const double clhs205 = C(0,1)*clhs135 + C(1,1)*clhs140 + C(1,2)*clhs202;
const double clhs206 = N[0]*N[2];
const double clhs207 = clhs17*(DN(1,0)*clhs80 + DN(1,1)*clhs84);
const double clhs208 = clhs17*(DN(1,0)*clhs156 + DN(1,1)*clhs157);
const double clhs209 = clhs179*(DN(1,0)*clhs177 + DN(1,1)*clhs178);
const double clhs210 = clhs179*(DN(1,0)*clhs178 + DN(1,1)*clhs181);
const double clhs211 = tau*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs212 = b_gauss[1]*clhs211;
const double clhs213 = b_gauss[0]*clhs211;
const double clhs214 = N[1]*N[2];
const double clhs215 = clhs17*(DN(2,0)*clhs80 + DN(2,1)*clhs84);
const double clhs216 = clhs17*(DN(2,0)*clhs156 + DN(2,1)*clhs157);
const double clhs217 = clhs179*(DN(2,0)*clhs177 + DN(2,1)*clhs178);
const double clhs218 = clhs179*(DN(2,0)*clhs178 + DN(2,1)*clhs181);
lhs(0,0)=DN(0,0)*clhs51 + DN(0,1)*clhs56;
lhs(0,1)=clhs38*(DN(0,0)*clhs69 + DN(0,1)*clhs71);
lhs(0,2)=N[0]*clhs85;
lhs(0,3)=DN(0,0)*clhs99 + DN(0,1)*clhs102;
lhs(0,4)=clhs38*(DN(0,0)*clhs115 + DN(0,1)*clhs117);
lhs(0,5)=N[1]*clhs85;
lhs(0,6)=DN(0,0)*clhs131 + DN(0,1)*clhs134;
lhs(0,7)=clhs38*(DN(0,0)*clhs147 + DN(0,1)*clhs149);
lhs(0,8)=N[2]*clhs85;
lhs(1,0)=clhs38*(DN(0,0)*clhs150 + DN(0,1)*clhs151);
lhs(1,1)=DN(0,0)*clhs154 + DN(0,1)*clhs155;
lhs(1,2)=N[0]*clhs158;
lhs(1,3)=clhs38*(DN(0,0)*clhs159 + DN(0,1)*clhs160);
lhs(1,4)=DN(0,0)*clhs161 + DN(0,1)*clhs162;
lhs(1,5)=N[1]*clhs158;
lhs(1,6)=clhs38*(DN(0,0)*clhs163 + DN(0,1)*clhs164);
lhs(1,7)=DN(0,0)*clhs165 + DN(0,1)*clhs166;
lhs(1,8)=N[2]*clhs158;
lhs(2,0)=N[0]*clhs21 - clhs170*(DN(0,0)*clhs168 + DN(0,1)*clhs169) - clhs172*(DN(0,0)*clhs169 + DN(0,1)*clhs171);
lhs(2,1)=N[0]*clhs58 - clhs170*(DN(0,0)*clhs174 + DN(0,1)*clhs175) - clhs172*(DN(0,0)*clhs175 + DN(0,1)*clhs176);
lhs(2,2)=-DN(0,0)*clhs180 - DN(0,1)*clhs182 - pow(N[0], 2);
lhs(2,3)=N[0]*clhs87 - clhs170*(DN(0,0)*clhs186 + DN(0,1)*clhs187) - clhs172*(DN(0,0)*clhs187 + DN(0,1)*clhs188) + clhs184;
lhs(2,4)=N[0]*clhs104 - clhs170*(DN(0,0)*clhs191 + DN(0,1)*clhs192) - clhs172*(DN(0,0)*clhs192 + DN(0,1)*clhs193) - clhs189;
lhs(2,5)=-DN(1,0)*clhs180 - DN(1,1)*clhs182 - clhs194;
lhs(2,6)=N[0]*clhs119 - clhs170*(DN(0,0)*clhs198 + DN(0,1)*clhs199) - clhs172*(DN(0,0)*clhs199 + DN(0,1)*clhs200) + clhs196;
lhs(2,7)=N[0]*clhs136 - clhs170*(DN(0,0)*clhs203 + DN(0,1)*clhs204) - clhs172*(DN(0,0)*clhs204 + DN(0,1)*clhs205) - clhs201;
lhs(2,8)=-DN(2,0)*clhs180 - DN(2,1)*clhs182 - clhs206;
lhs(3,0)=DN(1,0)*clhs51 + DN(1,1)*clhs56;
lhs(3,1)=clhs38*(DN(1,0)*clhs69 + DN(1,1)*clhs71);
lhs(3,2)=N[0]*clhs207;
lhs(3,3)=DN(1,0)*clhs99 + DN(1,1)*clhs102;
lhs(3,4)=clhs38*(DN(1,0)*clhs115 + DN(1,1)*clhs117);
lhs(3,5)=N[1]*clhs207;
lhs(3,6)=DN(1,0)*clhs131 + DN(1,1)*clhs134;
lhs(3,7)=clhs38*(DN(1,0)*clhs147 + DN(1,1)*clhs149);
lhs(3,8)=N[2]*clhs207;
lhs(4,0)=clhs38*(DN(1,0)*clhs150 + DN(1,1)*clhs151);
lhs(4,1)=DN(1,0)*clhs154 + DN(1,1)*clhs155;
lhs(4,2)=N[0]*clhs208;
lhs(4,3)=clhs38*(DN(1,0)*clhs159 + DN(1,1)*clhs160);
lhs(4,4)=DN(1,0)*clhs161 + DN(1,1)*clhs162;
lhs(4,5)=N[1]*clhs208;
lhs(4,6)=clhs38*(DN(1,0)*clhs163 + DN(1,1)*clhs164);
lhs(4,7)=DN(1,0)*clhs165 + DN(1,1)*clhs166;
lhs(4,8)=N[2]*clhs208;
lhs(5,0)=N[1]*clhs21 - clhs170*(DN(1,0)*clhs168 + DN(1,1)*clhs169) - clhs172*(DN(1,0)*clhs169 + DN(1,1)*clhs171) - clhs184;
lhs(5,1)=N[1]*clhs58 - clhs170*(DN(1,0)*clhs174 + DN(1,1)*clhs175) - clhs172*(DN(1,0)*clhs175 + DN(1,1)*clhs176) + clhs189;
lhs(5,2)=-DN(0,0)*clhs209 - DN(0,1)*clhs210 - clhs194;
lhs(5,3)=N[1]*clhs87 - clhs170*(DN(1,0)*clhs186 + DN(1,1)*clhs187) - clhs172*(DN(1,0)*clhs187 + DN(1,1)*clhs188);
lhs(5,4)=N[1]*clhs104 - clhs170*(DN(1,0)*clhs191 + DN(1,1)*clhs192) - clhs172*(DN(1,0)*clhs192 + DN(1,1)*clhs193);
lhs(5,5)=-DN(1,0)*clhs209 - DN(1,1)*clhs210 - pow(N[1], 2);
lhs(5,6)=N[1]*clhs119 - clhs170*(DN(1,0)*clhs198 + DN(1,1)*clhs199) - clhs172*(DN(1,0)*clhs199 + DN(1,1)*clhs200) + clhs212;
lhs(5,7)=N[1]*clhs136 - clhs170*(DN(1,0)*clhs203 + DN(1,1)*clhs204) - clhs172*(DN(1,0)*clhs204 + DN(1,1)*clhs205) - clhs213;
lhs(5,8)=-DN(2,0)*clhs209 - DN(2,1)*clhs210 - clhs214;
lhs(6,0)=DN(2,0)*clhs51 + DN(2,1)*clhs56;
lhs(6,1)=clhs38*(DN(2,0)*clhs69 + DN(2,1)*clhs71);
lhs(6,2)=N[0]*clhs215;
lhs(6,3)=DN(2,0)*clhs99 + DN(2,1)*clhs102;
lhs(6,4)=clhs38*(DN(2,0)*clhs115 + DN(2,1)*clhs117);
lhs(6,5)=N[1]*clhs215;
lhs(6,6)=DN(2,0)*clhs131 + DN(2,1)*clhs134;
lhs(6,7)=clhs38*(DN(2,0)*clhs147 + DN(2,1)*clhs149);
lhs(6,8)=N[2]*clhs215;
lhs(7,0)=clhs38*(DN(2,0)*clhs150 + DN(2,1)*clhs151);
lhs(7,1)=DN(2,0)*clhs154 + DN(2,1)*clhs155;
lhs(7,2)=N[0]*clhs216;
lhs(7,3)=clhs38*(DN(2,0)*clhs159 + DN(2,1)*clhs160);
lhs(7,4)=DN(2,0)*clhs161 + DN(2,1)*clhs162;
lhs(7,5)=N[1]*clhs216;
lhs(7,6)=clhs38*(DN(2,0)*clhs163 + DN(2,1)*clhs164);
lhs(7,7)=DN(2,0)*clhs165 + DN(2,1)*clhs166;
lhs(7,8)=N[2]*clhs216;
lhs(8,0)=N[2]*clhs21 - clhs170*(DN(2,0)*clhs168 + DN(2,1)*clhs169) - clhs172*(DN(2,0)*clhs169 + DN(2,1)*clhs171) - clhs196;
lhs(8,1)=N[2]*clhs58 - clhs170*(DN(2,0)*clhs174 + DN(2,1)*clhs175) - clhs172*(DN(2,0)*clhs175 + DN(2,1)*clhs176) + clhs201;
lhs(8,2)=-DN(0,0)*clhs217 - DN(0,1)*clhs218 - clhs206;
lhs(8,3)=N[2]*clhs87 - clhs170*(DN(2,0)*clhs186 + DN(2,1)*clhs187) - clhs172*(DN(2,0)*clhs187 + DN(2,1)*clhs188) - clhs212;
lhs(8,4)=N[2]*clhs104 - clhs170*(DN(2,0)*clhs191 + DN(2,1)*clhs192) - clhs172*(DN(2,0)*clhs192 + DN(2,1)*clhs193) + clhs213;
lhs(8,5)=-DN(1,0)*clhs217 - DN(1,1)*clhs218 - clhs214;
lhs(8,6)=N[2]*clhs119 - clhs170*(DN(2,0)*clhs198 + DN(2,1)*clhs199) - clhs172*(DN(2,0)*clhs199 + DN(2,1)*clhs200);
lhs(8,7)=N[2]*clhs136 - clhs170*(DN(2,0)*clhs203 + DN(2,1)*clhs204) - clhs172*(DN(2,0)*clhs204 + DN(2,1)*clhs205);
lhs(8,8)=-DN(2,0)*clhs217 - DN(2,1)*clhs218 - pow(N[2], 2);

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
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
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
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
        const double tau = aux_tau / mu;

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
const double crhs15 = DN(0,1)*u(0,1);
const double crhs16 = DN(1,1)*u(1,1);
const double crhs17 = DN(2,1)*u(2,1);
const double crhs18 = crhs15 + crhs16 + crhs17;
const double crhs19 = crhs18 + 1;
const double crhs20 = S[0]*crhs14 + S[2]*crhs19;
const double crhs21 = S[1]*crhs19 + S[2]*crhs14;
const double crhs22 = b_gauss[0]*tau;
const double crhs23 = b_gauss[1]*tau;
const double crhs24 = -N[0]*th[0] - N[1]*th[1] - N[2]*th[2] - crhs0*crhs12 - crhs0*crhs13 - crhs1*crhs11 - crhs1*crhs13 - crhs11*crhs2 - crhs12*crhs2 + crhs15*crhs5 + crhs15*crhs6 + crhs16*crhs4 + crhs16*crhs6 + crhs17*crhs4 + crhs17*crhs5 + crhs18 + crhs7;
const double crhs25 = pow(crhs14, 2) + pow(crhs8, 2);
const double crhs26 = pow(crhs19, 2) + pow(crhs3, 2);
const double crhs27 = 2*crhs14*crhs19 + 2*crhs3*crhs8;
const double crhs28 = C(0,0)*crhs25 + C(0,1)*crhs26 + C(0,2)*crhs27;
const double crhs29 = C(0,2)*crhs25 + C(1,2)*crhs26 + C(2,2)*crhs27;
const double crhs30 = (1.0/2.0)*tau;
const double crhs31 = crhs30*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double crhs32 = C(0,1)*crhs25 + C(1,1)*crhs26 + C(1,2)*crhs27;
const double crhs33 = crhs30*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
rhs[0]=-DN(0,0)*crhs9 - DN(0,1)*crhs10 + N[0]*b_gauss[0];
rhs[1]=-DN(0,0)*crhs20 - DN(0,1)*crhs21 + N[0]*b_gauss[1];
rhs[2]=-N[0]*crhs24 + crhs22*(DN(0,0)*crhs19 - DN(0,1)*crhs14) - crhs23*(DN(0,0)*crhs3 - DN(0,1)*crhs8) + crhs31*(DN(0,0)*crhs28 + DN(0,1)*crhs29) + crhs33*(DN(0,0)*crhs29 + DN(0,1)*crhs32);
rhs[3]=-DN(1,0)*crhs9 - DN(1,1)*crhs10 + N[1]*b_gauss[0];
rhs[4]=-DN(1,0)*crhs20 - DN(1,1)*crhs21 + N[1]*b_gauss[1];
rhs[5]=-N[1]*crhs24 + crhs22*(DN(1,0)*crhs19 - DN(1,1)*crhs14) - crhs23*(DN(1,0)*crhs3 - DN(1,1)*crhs8) + crhs31*(DN(1,0)*crhs28 + DN(1,1)*crhs29) + crhs33*(DN(1,0)*crhs29 + DN(1,1)*crhs32);
rhs[6]=-DN(2,0)*crhs9 - DN(2,1)*crhs10 + N[2]*b_gauss[0];
rhs[7]=-DN(2,0)*crhs20 - DN(2,1)*crhs21 + N[2]*b_gauss[1];
rhs[8]=-N[2]*crhs24 + crhs22*(DN(2,0)*crhs19 - DN(2,1)*crhs14) - crhs23*(DN(2,0)*crhs3 - DN(2,1)*crhs8) + crhs31*(DN(2,0)*crhs28 + DN(2,1)*crhs29) + crhs33*(DN(2,0)*crhs29 + DN(2,1)*crhs32);

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
        double mu = CalculateShearModulus(constitutive_variables.ConstitutiveMatrix); // 2nd Lame constant (shear modulus)
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
const double cr_eq_green_strain15 = pow(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0, 1.0)*1.0/(-cr_eq_green_strain0*cr_eq_green_strain10 - cr_eq_green_strain0*cr_eq_green_strain11 - cr_eq_green_strain1*cr_eq_green_strain11 - cr_eq_green_strain1*cr_eq_green_strain12 - cr_eq_green_strain10*cr_eq_green_strain2 - cr_eq_green_strain12*cr_eq_green_strain2 + cr_eq_green_strain13 + cr_eq_green_strain14 + cr_eq_green_strain4*cr_eq_green_strain5 + cr_eq_green_strain4*cr_eq_green_strain6 + cr_eq_green_strain5*cr_eq_green_strain9 + cr_eq_green_strain6*cr_eq_green_strain8 + cr_eq_green_strain7*cr_eq_green_strain8 + cr_eq_green_strain7*cr_eq_green_strain9);
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
const double cr_eq_def_gradient14 = sqrt(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0)*pow(cr_eq_def_gradient0*cr_eq_def_gradient4 + cr_eq_def_gradient0*cr_eq_def_gradient5 + cr_eq_def_gradient1*cr_eq_def_gradient5 + cr_eq_def_gradient1*cr_eq_def_gradient6 - cr_eq_def_gradient10*cr_eq_def_gradient11 - cr_eq_def_gradient10*cr_eq_def_gradient12 - cr_eq_def_gradient11*cr_eq_def_gradient9 - cr_eq_def_gradient12*cr_eq_def_gradient7 + cr_eq_def_gradient13 + cr_eq_def_gradient2*cr_eq_def_gradient4 + cr_eq_def_gradient2*cr_eq_def_gradient6 + cr_eq_def_gradient3 - cr_eq_def_gradient7*cr_eq_def_gradient8 - cr_eq_def_gradient8*cr_eq_def_gradient9, -0.5);
r_eq_def_gradient(0,0)=cr_eq_def_gradient14*cr_eq_def_gradient3;
r_eq_def_gradient(0,1)=cr_eq_def_gradient14*(cr_eq_def_gradient10 + cr_eq_def_gradient7 + cr_eq_def_gradient9);
r_eq_def_gradient(1,0)=cr_eq_def_gradient14*(cr_eq_def_gradient11 + cr_eq_def_gradient12 + cr_eq_def_gradient8);
r_eq_def_gradient(1,1)=cr_eq_def_gradient14*(cr_eq_def_gradient13 + 1);

    r_det_eq_def_gradient=pow(N[0]*th[0] + N[1]*th[1] + N[2]*th[2] + 1.0, 1.0);

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
