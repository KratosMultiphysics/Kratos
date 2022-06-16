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
const double clhs34 = DN(0,0)*clhs19 + DN(0,0)*clhs20 + DN(0,0) - DN(0,1)*DN(1,0)*u(1,1) - DN(0,1)*DN(2,0)*u(2,1);
const double clhs35 = -clhs0*clhs14 - clhs0*clhs15 - clhs1*clhs13 - clhs1*clhs15 - clhs13*clhs2 - clhs14*clhs2 + clhs18*clhs6 + clhs18*clhs7 + clhs19*clhs5 + clhs19*clhs7 + clhs20*clhs5 + clhs20*clhs6 + clhs21;
const double clhs36 = clhs35 + clhs9;
const double clhs37 = 1.0/clhs36;
const double clhs38 = clhs34*clhs37;
const double clhs39 = C(0,2)*clhs11;
const double clhs40 = C(1,2)*clhs24;
const double clhs41 = C(2,2)*clhs28;
const double clhs42 = clhs39 + clhs40 + 2*clhs41;
const double clhs43 = DN(0,1)*clhs42;
const double clhs44 = clhs32*clhs43;
const double clhs45 = DN(0,1)*clhs16;
const double clhs46 = DN(0,0)*clhs9;
const double clhs47 = DN(0,0)*clhs16;
const double clhs48 = DN(0,1)*clhs9;
const double clhs49 = clhs47 + clhs48;
const double clhs50 = C(0,0)*clhs46 + C(0,1)*clhs45 + C(0,2)*clhs49;
const double clhs51 = DN(0,0)*clhs50;
const double clhs52 = N[0]*th[0];
const double clhs53 = N[1]*th[1];
const double clhs54 = N[2]*th[2];
const double clhs55 = clhs35 - clhs52 - clhs53 - clhs54 + clhs8;
const double clhs56 = clhs37*clhs55*tau_th;
const double clhs57 = C(0,2)*clhs46 + C(1,2)*clhs45 + C(2,2)*clhs49;
const double clhs58 = DN(0,1)*clhs57;
const double clhs59 = pow(clhs36, -2.0);
const double clhs60 = clhs55*clhs59;
const double clhs61 = clhs34*clhs60;
const double clhs62 = 1.0*clhs37;
const double clhs63 = -clhs34*clhs59;
const double clhs64 = 0.5*clhs63;
const double clhs65 = clhs17*clhs64 + clhs23*clhs64 + clhs45*clhs62;
const double clhs66 = clhs10*clhs64 + clhs4*clhs64 + clhs46*clhs62;
const double clhs67 = clhs26*clhs63 + clhs27*clhs63 + clhs37*clhs47 + clhs37*clhs48;
const double clhs68 = C(0,0)*clhs66 + C(0,1)*clhs65 + C(0,2)*clhs67;
const double clhs69 = pow(clhs52 + clhs53 + clhs54 + 1.0, 1.0);
const double clhs70 = clhs69*clhs9;
const double clhs71 = C(0,2)*clhs66 + C(1,2)*clhs65 + C(2,2)*clhs67;
const double clhs72 = clhs16*clhs69;
const double clhs73 = clhs37*clhs69;
const double clhs74 = 0.5*clhs10*clhs73 + 0.5*clhs4*clhs73 - 0.5;
const double clhs75 = 0.5*clhs17*clhs73 + 0.5*clhs23*clhs73 - 0.5;
const double clhs76 = C(0,0)*clhs74 + C(0,1)*clhs75 + clhs29*clhs73;
const double clhs77 = C(0,2)*clhs74 + C(1,2)*clhs75 + clhs41*clhs73;
const double clhs78 = DN(0,0)*clhs76 + DN(0,1)*clhs77;
const double clhs79 = clhs68*clhs70 + clhs71*clhs72 + clhs78;
const double clhs80 = C(0,1)*clhs66 + C(1,1)*clhs65 + C(1,2)*clhs67;
const double clhs81 = C(1,2)*clhs28;
const double clhs82 = C(0,1)*clhs74 + C(1,1)*clhs75 + clhs73*clhs81;
const double clhs83 = DN(0,0)*clhs77 + DN(0,1)*clhs82;
const double clhs84 = clhs70*clhs71 + clhs72*clhs80 + clhs83;
const double clhs85 = -DN(0,0)*DN(1,1)*u(1,0) - DN(0,0)*DN(2,1)*u(2,0) + DN(0,1)*clhs6 + DN(0,1)*clhs7 + DN(0,1);
const double clhs86 = clhs37*clhs85;
const double clhs87 = DN(0,0)*clhs3;
const double clhs88 = DN(0,1)*clhs22;
const double clhs89 = DN(0,1)*clhs3;
const double clhs90 = DN(0,0)*clhs22;
const double clhs91 = clhs89 + clhs90;
const double clhs92 = C(0,0)*clhs87 + C(0,1)*clhs88 + C(0,2)*clhs91;
const double clhs93 = DN(0,0)*clhs92;
const double clhs94 = C(0,2)*clhs87 + C(1,2)*clhs88 + C(2,2)*clhs91;
const double clhs95 = DN(0,1)*clhs94;
const double clhs96 = clhs60*clhs85;
const double clhs97 = -clhs59*clhs85;
const double clhs98 = 0.5*clhs97;
const double clhs99 = clhs10*clhs98 + clhs4*clhs98 + clhs62*clhs87;
const double clhs100 = clhs17*clhs98 + clhs23*clhs98 + clhs62*clhs88;
const double clhs101 = clhs26*clhs97 + clhs27*clhs97 + clhs37*clhs89 + clhs37*clhs90;
const double clhs102 = C(0,2)*clhs99 + C(1,2)*clhs100 + C(2,2)*clhs101;
const double clhs103 = C(0,0)*clhs99 + C(0,1)*clhs100 + C(0,2)*clhs101;
const double clhs104 = clhs102*clhs16 + clhs103*clhs9;
const double clhs105 = DN(0,0)*clhs69;
const double clhs106 = C(0,1)*clhs99 + C(1,1)*clhs100 + C(1,2)*clhs101;
const double clhs107 = clhs102*clhs9 + clhs106*clhs16;
const double clhs108 = DN(0,1)*clhs69;
const double clhs109 = 0.5*clhs39 + 0.5*clhs40 + clhs41;
const double clhs110 = 0.5*clhs12 + 0.5*clhs25 + clhs29;
const double clhs111 = clhs109*clhs16 + clhs110*clhs9;
const double clhs112 = C(0,1)*clhs11;
const double clhs113 = C(1,1)*clhs24;
const double clhs114 = 0.5*clhs112 + 0.5*clhs113 + clhs81;
const double clhs115 = clhs109*clhs9 + clhs114*clhs16;
const double clhs116 = clhs37*(DN(0,0)*clhs111 + DN(0,1)*clhs115 - clhs33 - clhs44);
const double clhs117 = -DN(0,0)*DN(1,1)*u(0,1) + DN(1,0)*clhs18 + DN(1,0)*clhs20 + DN(1,0) - DN(1,1)*DN(2,0)*u(2,1);
const double clhs118 = clhs117*clhs37;
const double clhs119 = DN(1,1)*clhs16;
const double clhs120 = DN(1,0)*clhs9;
const double clhs121 = DN(1,0)*clhs16;
const double clhs122 = DN(1,1)*clhs9;
const double clhs123 = clhs121 + clhs122;
const double clhs124 = C(0,0)*clhs120 + C(0,1)*clhs119 + C(0,2)*clhs123;
const double clhs125 = DN(0,0)*clhs124;
const double clhs126 = C(0,2)*clhs120 + C(1,2)*clhs119 + C(2,2)*clhs123;
const double clhs127 = DN(0,1)*clhs126;
const double clhs128 = clhs117*clhs60;
const double clhs129 = -clhs117*clhs59;
const double clhs130 = 0.5*clhs129;
const double clhs131 = clhs119*clhs62 + clhs130*clhs17 + clhs130*clhs23;
const double clhs132 = clhs10*clhs130 + clhs120*clhs62 + clhs130*clhs4;
const double clhs133 = clhs121*clhs37 + clhs122*clhs37 + clhs129*clhs26 + clhs129*clhs27;
const double clhs134 = C(0,0)*clhs132 + C(0,1)*clhs131 + C(0,2)*clhs133;
const double clhs135 = C(0,2)*clhs132 + C(1,2)*clhs131 + C(2,2)*clhs133;
const double clhs136 = DN(1,0)*clhs76 + DN(1,1)*clhs77;
const double clhs137 = clhs134*clhs70 + clhs135*clhs72 + clhs136;
const double clhs138 = C(0,1)*clhs132 + C(1,1)*clhs131 + C(1,2)*clhs133;
const double clhs139 = DN(1,0)*clhs77 + DN(1,1)*clhs82;
const double clhs140 = clhs135*clhs70 + clhs138*clhs72 + clhs139;
const double clhs141 = -DN(0,1)*DN(1,0)*u(0,0) - DN(1,0)*DN(2,1)*u(2,0) + DN(1,1)*clhs5 + DN(1,1)*clhs7 + DN(1,1);
const double clhs142 = clhs141*clhs37;
const double clhs143 = DN(1,0)*clhs3;
const double clhs144 = DN(1,1)*clhs22;
const double clhs145 = DN(1,1)*clhs3;
const double clhs146 = DN(1,0)*clhs22;
const double clhs147 = clhs145 + clhs146;
const double clhs148 = C(0,0)*clhs143 + C(0,1)*clhs144 + C(0,2)*clhs147;
const double clhs149 = DN(0,0)*clhs148;
const double clhs150 = C(0,2)*clhs143 + C(1,2)*clhs144 + C(2,2)*clhs147;
const double clhs151 = DN(0,1)*clhs150;
const double clhs152 = clhs141*clhs60;
const double clhs153 = -clhs141*clhs59;
const double clhs154 = 0.5*clhs153;
const double clhs155 = clhs10*clhs154 + clhs143*clhs62 + clhs154*clhs4;
const double clhs156 = clhs144*clhs62 + clhs154*clhs17 + clhs154*clhs23;
const double clhs157 = clhs145*clhs37 + clhs146*clhs37 + clhs153*clhs26 + clhs153*clhs27;
const double clhs158 = C(0,2)*clhs155 + C(1,2)*clhs156 + C(2,2)*clhs157;
const double clhs159 = C(0,0)*clhs155 + C(0,1)*clhs156 + C(0,2)*clhs157;
const double clhs160 = clhs158*clhs16 + clhs159*clhs9;
const double clhs161 = C(0,1)*clhs155 + C(1,1)*clhs156 + C(1,2)*clhs157;
const double clhs162 = clhs158*clhs9 + clhs16*clhs161;
const double clhs163 = -DN(0,0)*DN(2,1)*u(0,1) - DN(1,0)*DN(2,1)*u(1,1) + DN(2,0)*clhs18 + DN(2,0)*clhs19 + DN(2,0);
const double clhs164 = clhs163*clhs37;
const double clhs165 = DN(2,1)*clhs16;
const double clhs166 = DN(2,0)*clhs9;
const double clhs167 = DN(2,0)*clhs16;
const double clhs168 = DN(2,1)*clhs9;
const double clhs169 = clhs167 + clhs168;
const double clhs170 = C(0,0)*clhs166 + C(0,1)*clhs165 + C(0,2)*clhs169;
const double clhs171 = DN(0,0)*clhs170;
const double clhs172 = C(0,2)*clhs166 + C(1,2)*clhs165 + C(2,2)*clhs169;
const double clhs173 = DN(0,1)*clhs172;
const double clhs174 = clhs163*clhs60;
const double clhs175 = -clhs163*clhs59;
const double clhs176 = 0.5*clhs175;
const double clhs177 = clhs165*clhs62 + clhs17*clhs176 + clhs176*clhs23;
const double clhs178 = clhs10*clhs176 + clhs166*clhs62 + clhs176*clhs4;
const double clhs179 = clhs167*clhs37 + clhs168*clhs37 + clhs175*clhs26 + clhs175*clhs27;
const double clhs180 = C(0,0)*clhs178 + C(0,1)*clhs177 + C(0,2)*clhs179;
const double clhs181 = C(0,2)*clhs178 + C(1,2)*clhs177 + C(2,2)*clhs179;
const double clhs182 = DN(2,0)*clhs76 + DN(2,1)*clhs77;
const double clhs183 = clhs180*clhs70 + clhs181*clhs72 + clhs182;
const double clhs184 = C(0,1)*clhs178 + C(1,1)*clhs177 + C(1,2)*clhs179;
const double clhs185 = DN(2,0)*clhs77 + DN(2,1)*clhs82;
const double clhs186 = clhs181*clhs70 + clhs184*clhs72 + clhs185;
const double clhs187 = -DN(0,1)*DN(2,0)*u(0,0) - DN(1,1)*DN(2,0)*u(1,0) + DN(2,1)*clhs5 + DN(2,1)*clhs6 + DN(2,1);
const double clhs188 = clhs187*clhs37;
const double clhs189 = DN(2,0)*clhs3;
const double clhs190 = DN(2,1)*clhs22;
const double clhs191 = DN(2,1)*clhs3;
const double clhs192 = DN(2,0)*clhs22;
const double clhs193 = clhs191 + clhs192;
const double clhs194 = C(0,0)*clhs189 + C(0,1)*clhs190 + C(0,2)*clhs193;
const double clhs195 = DN(0,0)*clhs194;
const double clhs196 = C(0,2)*clhs189 + C(1,2)*clhs190 + C(2,2)*clhs193;
const double clhs197 = DN(0,1)*clhs196;
const double clhs198 = clhs187*clhs60;
const double clhs199 = -clhs187*clhs59;
const double clhs200 = 0.5*clhs199;
const double clhs201 = clhs10*clhs200 + clhs189*clhs62 + clhs200*clhs4;
const double clhs202 = clhs17*clhs200 + clhs190*clhs62 + clhs200*clhs23;
const double clhs203 = clhs191*clhs37 + clhs192*clhs37 + clhs199*clhs26 + clhs199*clhs27;
const double clhs204 = C(0,2)*clhs201 + C(1,2)*clhs202 + C(2,2)*clhs203;
const double clhs205 = C(0,0)*clhs201 + C(0,1)*clhs202 + C(0,2)*clhs203;
const double clhs206 = clhs16*clhs204 + clhs205*clhs9;
const double clhs207 = C(0,1)*clhs201 + C(1,1)*clhs202 + C(1,2)*clhs203;
const double clhs208 = clhs16*clhs207 + clhs204*clhs9;
const double clhs209 = DN(0,0)*clhs42;
const double clhs210 = clhs209*clhs32;
const double clhs211 = clhs112 + clhs113 + 2*clhs81;
const double clhs212 = DN(0,1)*clhs211;
const double clhs213 = clhs212*clhs32;
const double clhs214 = DN(0,0)*clhs57;
const double clhs215 = C(0,1)*clhs46 + C(1,1)*clhs45 + C(1,2)*clhs49;
const double clhs216 = DN(0,1)*clhs215;
const double clhs217 = clhs22*clhs71 + clhs3*clhs68;
const double clhs218 = clhs22*clhs80 + clhs3*clhs71;
const double clhs219 = DN(0,0)*clhs94;
const double clhs220 = C(0,1)*clhs87 + C(1,1)*clhs88 + C(1,2)*clhs91;
const double clhs221 = DN(0,1)*clhs220;
const double clhs222 = clhs3*clhs69;
const double clhs223 = clhs22*clhs69;
const double clhs224 = clhs102*clhs223 + clhs103*clhs222 + clhs78;
const double clhs225 = clhs102*clhs222 + clhs106*clhs223 + clhs83;
const double clhs226 = clhs109*clhs22 + clhs110*clhs3;
const double clhs227 = clhs109*clhs3 + clhs114*clhs22;
const double clhs228 = clhs37*(DN(0,0)*clhs226 + DN(0,1)*clhs227 - clhs210 - clhs213);
const double clhs229 = DN(0,0)*clhs126;
const double clhs230 = C(0,1)*clhs120 + C(1,1)*clhs119 + C(1,2)*clhs123;
const double clhs231 = DN(0,1)*clhs230;
const double clhs232 = clhs134*clhs3 + clhs135*clhs22;
const double clhs233 = clhs135*clhs3 + clhs138*clhs22;
const double clhs234 = DN(0,0)*clhs150;
const double clhs235 = C(0,1)*clhs143 + C(1,1)*clhs144 + C(1,2)*clhs147;
const double clhs236 = DN(0,1)*clhs235;
const double clhs237 = clhs136 + clhs158*clhs223 + clhs159*clhs222;
const double clhs238 = clhs139 + clhs158*clhs222 + clhs161*clhs223;
const double clhs239 = DN(0,0)*clhs172;
const double clhs240 = C(0,1)*clhs166 + C(1,1)*clhs165 + C(1,2)*clhs169;
const double clhs241 = DN(0,1)*clhs240;
const double clhs242 = clhs180*clhs3 + clhs181*clhs22;
const double clhs243 = clhs181*clhs3 + clhs184*clhs22;
const double clhs244 = DN(0,0)*clhs196;
const double clhs245 = C(0,1)*clhs189 + C(1,1)*clhs190 + C(1,2)*clhs193;
const double clhs246 = DN(0,1)*clhs245;
const double clhs247 = clhs182 + clhs204*clhs223 + clhs205*clhs222;
const double clhs248 = clhs185 + clhs204*clhs222 + clhs207*clhs223;
const double clhs249 = tau_u*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double clhs250 = tau_u*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double clhs251 = (1.0/2.0)*tau_u;
const double clhs252 = clhs251*(clhs31 + clhs43);
const double clhs253 = clhs251*(clhs209 + clhs212);
const double clhs254 = tau_u*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs255 = b_gauss[1]*clhs254;
const double clhs256 = b_gauss[0]*clhs254;
const double clhs257 = N[0]*N[1];
const double clhs258 = tau_u*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs259 = b_gauss[1]*clhs258;
const double clhs260 = b_gauss[0]*clhs258;
const double clhs261 = N[0]*N[2];
const double clhs262 = DN(1,0)*clhs30;
const double clhs263 = clhs262*clhs32;
const double clhs264 = DN(1,1)*clhs42;
const double clhs265 = clhs264*clhs32;
const double clhs266 = DN(1,0)*clhs50;
const double clhs267 = DN(1,1)*clhs57;
const double clhs268 = DN(1,0)*clhs92;
const double clhs269 = DN(1,1)*clhs94;
const double clhs270 = DN(1,0)*clhs69;
const double clhs271 = DN(1,1)*clhs69;
const double clhs272 = clhs37*(DN(1,0)*clhs111 + DN(1,1)*clhs115 - clhs263 - clhs265);
const double clhs273 = DN(1,0)*clhs124;
const double clhs274 = DN(1,1)*clhs126;
const double clhs275 = DN(1,0)*clhs148;
const double clhs276 = DN(1,1)*clhs150;
const double clhs277 = DN(1,0)*clhs170;
const double clhs278 = DN(1,1)*clhs172;
const double clhs279 = DN(1,0)*clhs194;
const double clhs280 = DN(1,1)*clhs196;
const double clhs281 = DN(1,0)*clhs42;
const double clhs282 = clhs281*clhs32;
const double clhs283 = DN(1,1)*clhs211;
const double clhs284 = clhs283*clhs32;
const double clhs285 = DN(1,0)*clhs57;
const double clhs286 = DN(1,1)*clhs215;
const double clhs287 = DN(1,0)*clhs94;
const double clhs288 = DN(1,1)*clhs220;
const double clhs289 = clhs37*(DN(1,0)*clhs226 + DN(1,1)*clhs227 - clhs282 - clhs284);
const double clhs290 = DN(1,0)*clhs126;
const double clhs291 = DN(1,1)*clhs230;
const double clhs292 = DN(1,0)*clhs150;
const double clhs293 = DN(1,1)*clhs235;
const double clhs294 = DN(1,0)*clhs172;
const double clhs295 = DN(1,1)*clhs240;
const double clhs296 = DN(1,0)*clhs196;
const double clhs297 = DN(1,1)*clhs245;
const double clhs298 = clhs251*(clhs262 + clhs264);
const double clhs299 = clhs251*(clhs281 + clhs283);
const double clhs300 = tau_u*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs301 = b_gauss[1]*clhs300;
const double clhs302 = b_gauss[0]*clhs300;
const double clhs303 = N[1]*N[2];
const double clhs304 = DN(2,0)*clhs30;
const double clhs305 = clhs304*clhs32;
const double clhs306 = DN(2,1)*clhs42;
const double clhs307 = clhs306*clhs32;
const double clhs308 = DN(2,0)*clhs50;
const double clhs309 = DN(2,1)*clhs57;
const double clhs310 = DN(2,0)*clhs92;
const double clhs311 = DN(2,1)*clhs94;
const double clhs312 = DN(2,0)*clhs69;
const double clhs313 = DN(2,1)*clhs69;
const double clhs314 = clhs37*(DN(2,0)*clhs111 + DN(2,1)*clhs115 - clhs305 - clhs307);
const double clhs315 = DN(2,0)*clhs124;
const double clhs316 = DN(2,1)*clhs126;
const double clhs317 = DN(2,0)*clhs148;
const double clhs318 = DN(2,1)*clhs150;
const double clhs319 = DN(2,0)*clhs170;
const double clhs320 = DN(2,1)*clhs172;
const double clhs321 = DN(2,0)*clhs194;
const double clhs322 = DN(2,1)*clhs196;
const double clhs323 = DN(2,0)*clhs42;
const double clhs324 = clhs32*clhs323;
const double clhs325 = DN(2,1)*clhs211;
const double clhs326 = clhs32*clhs325;
const double clhs327 = DN(2,0)*clhs57;
const double clhs328 = DN(2,1)*clhs215;
const double clhs329 = DN(2,0)*clhs94;
const double clhs330 = DN(2,1)*clhs220;
const double clhs331 = clhs37*(DN(2,0)*clhs226 + DN(2,1)*clhs227 - clhs324 - clhs326);
const double clhs332 = DN(2,0)*clhs126;
const double clhs333 = DN(2,1)*clhs230;
const double clhs334 = DN(2,0)*clhs150;
const double clhs335 = DN(2,1)*clhs235;
const double clhs336 = DN(2,0)*clhs172;
const double clhs337 = DN(2,1)*clhs240;
const double clhs338 = DN(2,0)*clhs196;
const double clhs339 = DN(2,1)*clhs245;
const double clhs340 = clhs251*(clhs304 + clhs306);
const double clhs341 = clhs251*(clhs323 + clhs325);
lhs(0,0)=DN(0,0)*clhs79 + DN(0,1)*clhs84 + clhs33*clhs38 - clhs33*clhs61 + clhs38*clhs44 - clhs44*clhs61 + clhs51*clhs56 + clhs56*clhs58;
lhs(0,1)=clhs104*clhs105 + clhs107*clhs108 + clhs33*clhs86 - clhs33*clhs96 + clhs44*clhs86 - clhs44*clhs96 + clhs56*clhs93 + clhs56*clhs95;
lhs(0,2)=N[0]*clhs116;
lhs(0,3)=DN(0,0)*clhs137 + DN(0,1)*clhs140 + clhs118*clhs33 + clhs118*clhs44 + clhs125*clhs56 + clhs127*clhs56 - clhs128*clhs33 - clhs128*clhs44;
lhs(0,4)=clhs105*clhs160 + clhs108*clhs162 + clhs142*clhs33 + clhs142*clhs44 + clhs149*clhs56 + clhs151*clhs56 - clhs152*clhs33 - clhs152*clhs44;
lhs(0,5)=N[1]*clhs116;
lhs(0,6)=DN(0,0)*clhs183 + DN(0,1)*clhs186 + clhs164*clhs33 + clhs164*clhs44 + clhs171*clhs56 + clhs173*clhs56 - clhs174*clhs33 - clhs174*clhs44;
lhs(0,7)=clhs105*clhs206 + clhs108*clhs208 + clhs188*clhs33 + clhs188*clhs44 + clhs195*clhs56 + clhs197*clhs56 - clhs198*clhs33 - clhs198*clhs44;
lhs(0,8)=N[2]*clhs116;
lhs(1,0)=clhs105*clhs217 + clhs108*clhs218 + clhs210*clhs38 - clhs210*clhs61 + clhs213*clhs38 - clhs213*clhs61 + clhs214*clhs56 + clhs216*clhs56;
lhs(1,1)=DN(0,0)*clhs224 + DN(0,1)*clhs225 + clhs210*clhs86 - clhs210*clhs96 + clhs213*clhs86 - clhs213*clhs96 + clhs219*clhs56 + clhs221*clhs56;
lhs(1,2)=N[0]*clhs228;
lhs(1,3)=clhs105*clhs232 + clhs108*clhs233 + clhs118*clhs210 + clhs118*clhs213 - clhs128*clhs210 - clhs128*clhs213 + clhs229*clhs56 + clhs231*clhs56;
lhs(1,4)=DN(0,0)*clhs237 + DN(0,1)*clhs238 + clhs142*clhs210 + clhs142*clhs213 - clhs152*clhs210 - clhs152*clhs213 + clhs234*clhs56 + clhs236*clhs56;
lhs(1,5)=N[1]*clhs228;
lhs(1,6)=clhs105*clhs242 + clhs108*clhs243 + clhs164*clhs210 + clhs164*clhs213 - clhs174*clhs210 - clhs174*clhs213 + clhs239*clhs56 + clhs241*clhs56;
lhs(1,7)=DN(0,0)*clhs247 + DN(0,1)*clhs248 + clhs188*clhs210 + clhs188*clhs213 - clhs198*clhs210 - clhs198*clhs213 + clhs244*clhs56 + clhs246*clhs56;
lhs(1,8)=N[2]*clhs228;
lhs(2,0)=N[0]*clhs34 - clhs249*(clhs51 + clhs58) - clhs250*(clhs214 + clhs216);
lhs(2,1)=N[0]*clhs85 - clhs249*(clhs93 + clhs95) - clhs250*(clhs219 + clhs221);
lhs(2,2)=-DN(0,0)*clhs252 - DN(0,1)*clhs253 - pow(N[0], 2);
lhs(2,3)=N[0]*clhs117 - clhs249*(clhs125 + clhs127) - clhs250*(clhs229 + clhs231) + clhs255;
lhs(2,4)=N[0]*clhs141 - clhs249*(clhs149 + clhs151) - clhs250*(clhs234 + clhs236) - clhs256;
lhs(2,5)=-DN(1,0)*clhs252 - DN(1,1)*clhs253 - clhs257;
lhs(2,6)=N[0]*clhs163 - clhs249*(clhs171 + clhs173) - clhs250*(clhs239 + clhs241) + clhs259;
lhs(2,7)=N[0]*clhs187 - clhs249*(clhs195 + clhs197) - clhs250*(clhs244 + clhs246) - clhs260;
lhs(2,8)=-DN(2,0)*clhs252 - DN(2,1)*clhs253 - clhs261;
lhs(3,0)=DN(1,0)*clhs79 + DN(1,1)*clhs84 + clhs263*clhs38 - clhs263*clhs61 + clhs265*clhs38 - clhs265*clhs61 + clhs266*clhs56 + clhs267*clhs56;
lhs(3,1)=clhs104*clhs270 + clhs107*clhs271 + clhs263*clhs86 - clhs263*clhs96 + clhs265*clhs86 - clhs265*clhs96 + clhs268*clhs56 + clhs269*clhs56;
lhs(3,2)=N[0]*clhs272;
lhs(3,3)=DN(1,0)*clhs137 + DN(1,1)*clhs140 + clhs118*clhs263 + clhs118*clhs265 - clhs128*clhs263 - clhs128*clhs265 + clhs273*clhs56 + clhs274*clhs56;
lhs(3,4)=clhs142*clhs263 + clhs142*clhs265 - clhs152*clhs263 - clhs152*clhs265 + clhs160*clhs270 + clhs162*clhs271 + clhs275*clhs56 + clhs276*clhs56;
lhs(3,5)=N[1]*clhs272;
lhs(3,6)=DN(1,0)*clhs183 + DN(1,1)*clhs186 + clhs164*clhs263 + clhs164*clhs265 - clhs174*clhs263 - clhs174*clhs265 + clhs277*clhs56 + clhs278*clhs56;
lhs(3,7)=clhs188*clhs263 + clhs188*clhs265 - clhs198*clhs263 - clhs198*clhs265 + clhs206*clhs270 + clhs208*clhs271 + clhs279*clhs56 + clhs280*clhs56;
lhs(3,8)=N[2]*clhs272;
lhs(4,0)=clhs217*clhs270 + clhs218*clhs271 + clhs282*clhs38 - clhs282*clhs61 + clhs284*clhs38 - clhs284*clhs61 + clhs285*clhs56 + clhs286*clhs56;
lhs(4,1)=DN(1,0)*clhs224 + DN(1,1)*clhs225 + clhs282*clhs86 - clhs282*clhs96 + clhs284*clhs86 - clhs284*clhs96 + clhs287*clhs56 + clhs288*clhs56;
lhs(4,2)=N[0]*clhs289;
lhs(4,3)=clhs118*clhs282 + clhs118*clhs284 - clhs128*clhs282 - clhs128*clhs284 + clhs232*clhs270 + clhs233*clhs271 + clhs290*clhs56 + clhs291*clhs56;
lhs(4,4)=DN(1,0)*clhs237 + DN(1,1)*clhs238 + clhs142*clhs282 + clhs142*clhs284 - clhs152*clhs282 - clhs152*clhs284 + clhs292*clhs56 + clhs293*clhs56;
lhs(4,5)=N[1]*clhs289;
lhs(4,6)=clhs164*clhs282 + clhs164*clhs284 - clhs174*clhs282 - clhs174*clhs284 + clhs242*clhs270 + clhs243*clhs271 + clhs294*clhs56 + clhs295*clhs56;
lhs(4,7)=DN(1,0)*clhs247 + DN(1,1)*clhs248 + clhs188*clhs282 + clhs188*clhs284 - clhs198*clhs282 - clhs198*clhs284 + clhs296*clhs56 + clhs297*clhs56;
lhs(4,8)=N[2]*clhs289;
lhs(5,0)=N[1]*clhs34 - clhs249*(clhs266 + clhs267) - clhs250*(clhs285 + clhs286) - clhs255;
lhs(5,1)=N[1]*clhs85 - clhs249*(clhs268 + clhs269) - clhs250*(clhs287 + clhs288) + clhs256;
lhs(5,2)=-DN(0,0)*clhs298 - DN(0,1)*clhs299 - clhs257;
lhs(5,3)=N[1]*clhs117 - clhs249*(clhs273 + clhs274) - clhs250*(clhs290 + clhs291);
lhs(5,4)=N[1]*clhs141 - clhs249*(clhs275 + clhs276) - clhs250*(clhs292 + clhs293);
lhs(5,5)=-DN(1,0)*clhs298 - DN(1,1)*clhs299 - pow(N[1], 2);
lhs(5,6)=N[1]*clhs163 - clhs249*(clhs277 + clhs278) - clhs250*(clhs294 + clhs295) + clhs301;
lhs(5,7)=N[1]*clhs187 - clhs249*(clhs279 + clhs280) - clhs250*(clhs296 + clhs297) - clhs302;
lhs(5,8)=-DN(2,0)*clhs298 - DN(2,1)*clhs299 - clhs303;
lhs(6,0)=DN(2,0)*clhs79 + DN(2,1)*clhs84 + clhs305*clhs38 - clhs305*clhs61 + clhs307*clhs38 - clhs307*clhs61 + clhs308*clhs56 + clhs309*clhs56;
lhs(6,1)=clhs104*clhs312 + clhs107*clhs313 + clhs305*clhs86 - clhs305*clhs96 + clhs307*clhs86 - clhs307*clhs96 + clhs310*clhs56 + clhs311*clhs56;
lhs(6,2)=N[0]*clhs314;
lhs(6,3)=DN(2,0)*clhs137 + DN(2,1)*clhs140 + clhs118*clhs305 + clhs118*clhs307 - clhs128*clhs305 - clhs128*clhs307 + clhs315*clhs56 + clhs316*clhs56;
lhs(6,4)=clhs142*clhs305 + clhs142*clhs307 - clhs152*clhs305 - clhs152*clhs307 + clhs160*clhs312 + clhs162*clhs313 + clhs317*clhs56 + clhs318*clhs56;
lhs(6,5)=N[1]*clhs314;
lhs(6,6)=DN(2,0)*clhs183 + DN(2,1)*clhs186 + clhs164*clhs305 + clhs164*clhs307 - clhs174*clhs305 - clhs174*clhs307 + clhs319*clhs56 + clhs320*clhs56;
lhs(6,7)=clhs188*clhs305 + clhs188*clhs307 - clhs198*clhs305 - clhs198*clhs307 + clhs206*clhs312 + clhs208*clhs313 + clhs321*clhs56 + clhs322*clhs56;
lhs(6,8)=N[2]*clhs314;
lhs(7,0)=clhs217*clhs312 + clhs218*clhs313 + clhs324*clhs38 - clhs324*clhs61 + clhs326*clhs38 - clhs326*clhs61 + clhs327*clhs56 + clhs328*clhs56;
lhs(7,1)=DN(2,0)*clhs224 + DN(2,1)*clhs225 + clhs324*clhs86 - clhs324*clhs96 + clhs326*clhs86 - clhs326*clhs96 + clhs329*clhs56 + clhs330*clhs56;
lhs(7,2)=N[0]*clhs331;
lhs(7,3)=clhs118*clhs324 + clhs118*clhs326 - clhs128*clhs324 - clhs128*clhs326 + clhs232*clhs312 + clhs233*clhs313 + clhs332*clhs56 + clhs333*clhs56;
lhs(7,4)=DN(2,0)*clhs237 + DN(2,1)*clhs238 + clhs142*clhs324 + clhs142*clhs326 - clhs152*clhs324 - clhs152*clhs326 + clhs334*clhs56 + clhs335*clhs56;
lhs(7,5)=N[1]*clhs331;
lhs(7,6)=clhs164*clhs324 + clhs164*clhs326 - clhs174*clhs324 - clhs174*clhs326 + clhs242*clhs312 + clhs243*clhs313 + clhs336*clhs56 + clhs337*clhs56;
lhs(7,7)=DN(2,0)*clhs247 + DN(2,1)*clhs248 + clhs188*clhs324 + clhs188*clhs326 - clhs198*clhs324 - clhs198*clhs326 + clhs338*clhs56 + clhs339*clhs56;
lhs(7,8)=N[2]*clhs331;
lhs(8,0)=N[2]*clhs34 - clhs249*(clhs308 + clhs309) - clhs250*(clhs327 + clhs328) - clhs259;
lhs(8,1)=N[2]*clhs85 - clhs249*(clhs310 + clhs311) - clhs250*(clhs329 + clhs330) + clhs260;
lhs(8,2)=-DN(0,0)*clhs340 - DN(0,1)*clhs341 - clhs261;
lhs(8,3)=N[2]*clhs117 - clhs249*(clhs315 + clhs316) - clhs250*(clhs332 + clhs333) - clhs301;
lhs(8,4)=N[2]*clhs141 - clhs249*(clhs317 + clhs318) - clhs250*(clhs334 + clhs335) + clhs302;
lhs(8,5)=-DN(1,0)*clhs340 - DN(1,1)*clhs341 - clhs303;
lhs(8,6)=N[2]*clhs163 - clhs249*(clhs319 + clhs320) - clhs250*(clhs336 + clhs337);
lhs(8,7)=N[2]*clhs187 - clhs249*(clhs321 + clhs322) - clhs250*(clhs338 + clhs339);
lhs(8,8)=-DN(2,0)*clhs340 - DN(2,1)*clhs341 - pow(N[2], 2);

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
const double crhs27 = 0.5*crhs26*tau_th*1.0/(crhs25 + crhs8);
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
const double clhs34 = DN(0,0)*clhs19 + DN(0,0)*clhs20 + DN(0,0) - DN(0,1)*DN(1,0)*u(1,1) - DN(0,1)*DN(2,0)*u(2,1);
const double clhs35 = -clhs0*clhs14 - clhs0*clhs15 - clhs1*clhs13 - clhs1*clhs15 - clhs13*clhs2 - clhs14*clhs2 + clhs18*clhs6 + clhs18*clhs7 + clhs19*clhs5 + clhs19*clhs7 + clhs20*clhs5 + clhs20*clhs6 + clhs21;
const double clhs36 = clhs35 + clhs9;
const double clhs37 = 1.0/clhs36;
const double clhs38 = clhs34*clhs37;
const double clhs39 = C(0,2)*clhs11;
const double clhs40 = C(1,2)*clhs24;
const double clhs41 = C(2,2)*clhs28;
const double clhs42 = clhs39 + clhs40 + 2*clhs41;
const double clhs43 = DN(0,1)*clhs42;
const double clhs44 = clhs32*clhs43;
const double clhs45 = DN(0,1)*clhs16;
const double clhs46 = DN(0,0)*clhs9;
const double clhs47 = DN(0,0)*clhs16;
const double clhs48 = DN(0,1)*clhs9;
const double clhs49 = clhs47 + clhs48;
const double clhs50 = C(0,0)*clhs46 + C(0,1)*clhs45 + C(0,2)*clhs49;
const double clhs51 = DN(0,0)*clhs50;
const double clhs52 = N[0]*th[0];
const double clhs53 = N[1]*th[1];
const double clhs54 = N[2]*th[2];
const double clhs55 = clhs35 - clhs52 - clhs53 - clhs54 + clhs8;
const double clhs56 = clhs37*clhs55*tau_th;
const double clhs57 = C(0,2)*clhs46 + C(1,2)*clhs45 + C(2,2)*clhs49;
const double clhs58 = DN(0,1)*clhs57;
const double clhs59 = pow(clhs36, -2.0);
const double clhs60 = clhs55*clhs59;
const double clhs61 = clhs34*clhs60;
const double clhs62 = 1.0*clhs37;
const double clhs63 = -clhs34*clhs59;
const double clhs64 = 0.5*clhs63;
const double clhs65 = clhs17*clhs64 + clhs23*clhs64 + clhs45*clhs62;
const double clhs66 = clhs10*clhs64 + clhs4*clhs64 + clhs46*clhs62;
const double clhs67 = clhs26*clhs63 + clhs27*clhs63 + clhs37*clhs47 + clhs37*clhs48;
const double clhs68 = C(0,0)*clhs66 + C(0,1)*clhs65 + C(0,2)*clhs67;
const double clhs69 = pow(clhs52 + clhs53 + clhs54 + 1.0, 1.0);
const double clhs70 = clhs69*clhs9;
const double clhs71 = C(0,2)*clhs66 + C(1,2)*clhs65 + C(2,2)*clhs67;
const double clhs72 = clhs16*clhs69;
const double clhs73 = clhs37*clhs69;
const double clhs74 = 0.5*clhs10*clhs73 + 0.5*clhs4*clhs73 - 0.5;
const double clhs75 = 0.5*clhs17*clhs73 + 0.5*clhs23*clhs73 - 0.5;
const double clhs76 = C(0,0)*clhs74 + C(0,1)*clhs75 + clhs29*clhs73;
const double clhs77 = C(0,2)*clhs74 + C(1,2)*clhs75 + clhs41*clhs73;
const double clhs78 = DN(0,0)*clhs76 + DN(0,1)*clhs77;
const double clhs79 = clhs68*clhs70 + clhs71*clhs72 + clhs78;
const double clhs80 = C(0,1)*clhs66 + C(1,1)*clhs65 + C(1,2)*clhs67;
const double clhs81 = C(1,2)*clhs28;
const double clhs82 = C(0,1)*clhs74 + C(1,1)*clhs75 + clhs73*clhs81;
const double clhs83 = DN(0,0)*clhs77 + DN(0,1)*clhs82;
const double clhs84 = clhs70*clhs71 + clhs72*clhs80 + clhs83;
const double clhs85 = -DN(0,0)*DN(1,1)*u(1,0) - DN(0,0)*DN(2,1)*u(2,0) + DN(0,1)*clhs6 + DN(0,1)*clhs7 + DN(0,1);
const double clhs86 = clhs37*clhs85;
const double clhs87 = DN(0,0)*clhs3;
const double clhs88 = DN(0,1)*clhs22;
const double clhs89 = DN(0,1)*clhs3;
const double clhs90 = DN(0,0)*clhs22;
const double clhs91 = clhs89 + clhs90;
const double clhs92 = C(0,0)*clhs87 + C(0,1)*clhs88 + C(0,2)*clhs91;
const double clhs93 = DN(0,0)*clhs92;
const double clhs94 = C(0,2)*clhs87 + C(1,2)*clhs88 + C(2,2)*clhs91;
const double clhs95 = DN(0,1)*clhs94;
const double clhs96 = clhs60*clhs85;
const double clhs97 = -clhs59*clhs85;
const double clhs98 = 0.5*clhs97;
const double clhs99 = clhs10*clhs98 + clhs4*clhs98 + clhs62*clhs87;
const double clhs100 = clhs17*clhs98 + clhs23*clhs98 + clhs62*clhs88;
const double clhs101 = clhs26*clhs97 + clhs27*clhs97 + clhs37*clhs89 + clhs37*clhs90;
const double clhs102 = C(0,2)*clhs99 + C(1,2)*clhs100 + C(2,2)*clhs101;
const double clhs103 = C(0,0)*clhs99 + C(0,1)*clhs100 + C(0,2)*clhs101;
const double clhs104 = clhs102*clhs16 + clhs103*clhs9;
const double clhs105 = DN(0,0)*clhs69;
const double clhs106 = C(0,1)*clhs99 + C(1,1)*clhs100 + C(1,2)*clhs101;
const double clhs107 = clhs102*clhs9 + clhs106*clhs16;
const double clhs108 = DN(0,1)*clhs69;
const double clhs109 = 0.5*clhs39 + 0.5*clhs40 + clhs41;
const double clhs110 = 0.5*clhs12 + 0.5*clhs25 + clhs29;
const double clhs111 = clhs109*clhs16 + clhs110*clhs9;
const double clhs112 = C(0,1)*clhs11;
const double clhs113 = C(1,1)*clhs24;
const double clhs114 = 0.5*clhs112 + 0.5*clhs113 + clhs81;
const double clhs115 = clhs109*clhs9 + clhs114*clhs16;
const double clhs116 = clhs37*(DN(0,0)*clhs111 + DN(0,1)*clhs115 - clhs33 - clhs44);
const double clhs117 = -DN(0,0)*DN(1,1)*u(0,1) + DN(1,0)*clhs18 + DN(1,0)*clhs20 + DN(1,0) - DN(1,1)*DN(2,0)*u(2,1);
const double clhs118 = clhs117*clhs37;
const double clhs119 = DN(1,1)*clhs16;
const double clhs120 = DN(1,0)*clhs9;
const double clhs121 = DN(1,0)*clhs16;
const double clhs122 = DN(1,1)*clhs9;
const double clhs123 = clhs121 + clhs122;
const double clhs124 = C(0,0)*clhs120 + C(0,1)*clhs119 + C(0,2)*clhs123;
const double clhs125 = DN(0,0)*clhs124;
const double clhs126 = C(0,2)*clhs120 + C(1,2)*clhs119 + C(2,2)*clhs123;
const double clhs127 = DN(0,1)*clhs126;
const double clhs128 = clhs117*clhs60;
const double clhs129 = -clhs117*clhs59;
const double clhs130 = 0.5*clhs129;
const double clhs131 = clhs119*clhs62 + clhs130*clhs17 + clhs130*clhs23;
const double clhs132 = clhs10*clhs130 + clhs120*clhs62 + clhs130*clhs4;
const double clhs133 = clhs121*clhs37 + clhs122*clhs37 + clhs129*clhs26 + clhs129*clhs27;
const double clhs134 = C(0,0)*clhs132 + C(0,1)*clhs131 + C(0,2)*clhs133;
const double clhs135 = C(0,2)*clhs132 + C(1,2)*clhs131 + C(2,2)*clhs133;
const double clhs136 = DN(1,0)*clhs76 + DN(1,1)*clhs77;
const double clhs137 = clhs134*clhs70 + clhs135*clhs72 + clhs136;
const double clhs138 = C(0,1)*clhs132 + C(1,1)*clhs131 + C(1,2)*clhs133;
const double clhs139 = DN(1,0)*clhs77 + DN(1,1)*clhs82;
const double clhs140 = clhs135*clhs70 + clhs138*clhs72 + clhs139;
const double clhs141 = -DN(0,1)*DN(1,0)*u(0,0) - DN(1,0)*DN(2,1)*u(2,0) + DN(1,1)*clhs5 + DN(1,1)*clhs7 + DN(1,1);
const double clhs142 = clhs141*clhs37;
const double clhs143 = DN(1,0)*clhs3;
const double clhs144 = DN(1,1)*clhs22;
const double clhs145 = DN(1,1)*clhs3;
const double clhs146 = DN(1,0)*clhs22;
const double clhs147 = clhs145 + clhs146;
const double clhs148 = C(0,0)*clhs143 + C(0,1)*clhs144 + C(0,2)*clhs147;
const double clhs149 = DN(0,0)*clhs148;
const double clhs150 = C(0,2)*clhs143 + C(1,2)*clhs144 + C(2,2)*clhs147;
const double clhs151 = DN(0,1)*clhs150;
const double clhs152 = clhs141*clhs60;
const double clhs153 = -clhs141*clhs59;
const double clhs154 = 0.5*clhs153;
const double clhs155 = clhs10*clhs154 + clhs143*clhs62 + clhs154*clhs4;
const double clhs156 = clhs144*clhs62 + clhs154*clhs17 + clhs154*clhs23;
const double clhs157 = clhs145*clhs37 + clhs146*clhs37 + clhs153*clhs26 + clhs153*clhs27;
const double clhs158 = C(0,2)*clhs155 + C(1,2)*clhs156 + C(2,2)*clhs157;
const double clhs159 = C(0,0)*clhs155 + C(0,1)*clhs156 + C(0,2)*clhs157;
const double clhs160 = clhs158*clhs16 + clhs159*clhs9;
const double clhs161 = C(0,1)*clhs155 + C(1,1)*clhs156 + C(1,2)*clhs157;
const double clhs162 = clhs158*clhs9 + clhs16*clhs161;
const double clhs163 = -DN(0,0)*DN(2,1)*u(0,1) - DN(1,0)*DN(2,1)*u(1,1) + DN(2,0)*clhs18 + DN(2,0)*clhs19 + DN(2,0);
const double clhs164 = clhs163*clhs37;
const double clhs165 = DN(2,1)*clhs16;
const double clhs166 = DN(2,0)*clhs9;
const double clhs167 = DN(2,0)*clhs16;
const double clhs168 = DN(2,1)*clhs9;
const double clhs169 = clhs167 + clhs168;
const double clhs170 = C(0,0)*clhs166 + C(0,1)*clhs165 + C(0,2)*clhs169;
const double clhs171 = DN(0,0)*clhs170;
const double clhs172 = C(0,2)*clhs166 + C(1,2)*clhs165 + C(2,2)*clhs169;
const double clhs173 = DN(0,1)*clhs172;
const double clhs174 = clhs163*clhs60;
const double clhs175 = -clhs163*clhs59;
const double clhs176 = 0.5*clhs175;
const double clhs177 = clhs165*clhs62 + clhs17*clhs176 + clhs176*clhs23;
const double clhs178 = clhs10*clhs176 + clhs166*clhs62 + clhs176*clhs4;
const double clhs179 = clhs167*clhs37 + clhs168*clhs37 + clhs175*clhs26 + clhs175*clhs27;
const double clhs180 = C(0,0)*clhs178 + C(0,1)*clhs177 + C(0,2)*clhs179;
const double clhs181 = C(0,2)*clhs178 + C(1,2)*clhs177 + C(2,2)*clhs179;
const double clhs182 = DN(2,0)*clhs76 + DN(2,1)*clhs77;
const double clhs183 = clhs180*clhs70 + clhs181*clhs72 + clhs182;
const double clhs184 = C(0,1)*clhs178 + C(1,1)*clhs177 + C(1,2)*clhs179;
const double clhs185 = DN(2,0)*clhs77 + DN(2,1)*clhs82;
const double clhs186 = clhs181*clhs70 + clhs184*clhs72 + clhs185;
const double clhs187 = -DN(0,1)*DN(2,0)*u(0,0) - DN(1,1)*DN(2,0)*u(1,0) + DN(2,1)*clhs5 + DN(2,1)*clhs6 + DN(2,1);
const double clhs188 = clhs187*clhs37;
const double clhs189 = DN(2,0)*clhs3;
const double clhs190 = DN(2,1)*clhs22;
const double clhs191 = DN(2,1)*clhs3;
const double clhs192 = DN(2,0)*clhs22;
const double clhs193 = clhs191 + clhs192;
const double clhs194 = C(0,0)*clhs189 + C(0,1)*clhs190 + C(0,2)*clhs193;
const double clhs195 = DN(0,0)*clhs194;
const double clhs196 = C(0,2)*clhs189 + C(1,2)*clhs190 + C(2,2)*clhs193;
const double clhs197 = DN(0,1)*clhs196;
const double clhs198 = clhs187*clhs60;
const double clhs199 = -clhs187*clhs59;
const double clhs200 = 0.5*clhs199;
const double clhs201 = clhs10*clhs200 + clhs189*clhs62 + clhs200*clhs4;
const double clhs202 = clhs17*clhs200 + clhs190*clhs62 + clhs200*clhs23;
const double clhs203 = clhs191*clhs37 + clhs192*clhs37 + clhs199*clhs26 + clhs199*clhs27;
const double clhs204 = C(0,2)*clhs201 + C(1,2)*clhs202 + C(2,2)*clhs203;
const double clhs205 = C(0,0)*clhs201 + C(0,1)*clhs202 + C(0,2)*clhs203;
const double clhs206 = clhs16*clhs204 + clhs205*clhs9;
const double clhs207 = C(0,1)*clhs201 + C(1,1)*clhs202 + C(1,2)*clhs203;
const double clhs208 = clhs16*clhs207 + clhs204*clhs9;
const double clhs209 = DN(0,0)*clhs42;
const double clhs210 = clhs209*clhs32;
const double clhs211 = clhs112 + clhs113 + 2*clhs81;
const double clhs212 = DN(0,1)*clhs211;
const double clhs213 = clhs212*clhs32;
const double clhs214 = DN(0,0)*clhs57;
const double clhs215 = C(0,1)*clhs46 + C(1,1)*clhs45 + C(1,2)*clhs49;
const double clhs216 = DN(0,1)*clhs215;
const double clhs217 = clhs22*clhs71 + clhs3*clhs68;
const double clhs218 = clhs22*clhs80 + clhs3*clhs71;
const double clhs219 = DN(0,0)*clhs94;
const double clhs220 = C(0,1)*clhs87 + C(1,1)*clhs88 + C(1,2)*clhs91;
const double clhs221 = DN(0,1)*clhs220;
const double clhs222 = clhs3*clhs69;
const double clhs223 = clhs22*clhs69;
const double clhs224 = clhs102*clhs223 + clhs103*clhs222 + clhs78;
const double clhs225 = clhs102*clhs222 + clhs106*clhs223 + clhs83;
const double clhs226 = clhs109*clhs22 + clhs110*clhs3;
const double clhs227 = clhs109*clhs3 + clhs114*clhs22;
const double clhs228 = clhs37*(DN(0,0)*clhs226 + DN(0,1)*clhs227 - clhs210 - clhs213);
const double clhs229 = DN(0,0)*clhs126;
const double clhs230 = C(0,1)*clhs120 + C(1,1)*clhs119 + C(1,2)*clhs123;
const double clhs231 = DN(0,1)*clhs230;
const double clhs232 = clhs134*clhs3 + clhs135*clhs22;
const double clhs233 = clhs135*clhs3 + clhs138*clhs22;
const double clhs234 = DN(0,0)*clhs150;
const double clhs235 = C(0,1)*clhs143 + C(1,1)*clhs144 + C(1,2)*clhs147;
const double clhs236 = DN(0,1)*clhs235;
const double clhs237 = clhs136 + clhs158*clhs223 + clhs159*clhs222;
const double clhs238 = clhs139 + clhs158*clhs222 + clhs161*clhs223;
const double clhs239 = DN(0,0)*clhs172;
const double clhs240 = C(0,1)*clhs166 + C(1,1)*clhs165 + C(1,2)*clhs169;
const double clhs241 = DN(0,1)*clhs240;
const double clhs242 = clhs180*clhs3 + clhs181*clhs22;
const double clhs243 = clhs181*clhs3 + clhs184*clhs22;
const double clhs244 = DN(0,0)*clhs196;
const double clhs245 = C(0,1)*clhs189 + C(1,1)*clhs190 + C(1,2)*clhs193;
const double clhs246 = DN(0,1)*clhs245;
const double clhs247 = clhs182 + clhs204*clhs223 + clhs205*clhs222;
const double clhs248 = clhs185 + clhs204*clhs222 + clhs207*clhs223;
const double clhs249 = tau_u*(DN(0,0)*th[0] + DN(1,0)*th[1] + DN(2,0)*th[2]);
const double clhs250 = tau_u*(DN(0,1)*th[0] + DN(1,1)*th[1] + DN(2,1)*th[2]);
const double clhs251 = (1.0/2.0)*tau_u;
const double clhs252 = clhs251*(clhs31 + clhs43);
const double clhs253 = clhs251*(clhs209 + clhs212);
const double clhs254 = tau_u*(DN(0,0)*DN(1,1) - DN(0,1)*DN(1,0));
const double clhs255 = b_gauss[1]*clhs254;
const double clhs256 = b_gauss[0]*clhs254;
const double clhs257 = N[0]*N[1];
const double clhs258 = tau_u*(DN(0,0)*DN(2,1) - DN(0,1)*DN(2,0));
const double clhs259 = b_gauss[1]*clhs258;
const double clhs260 = b_gauss[0]*clhs258;
const double clhs261 = N[0]*N[2];
const double clhs262 = DN(1,0)*clhs30;
const double clhs263 = clhs262*clhs32;
const double clhs264 = DN(1,1)*clhs42;
const double clhs265 = clhs264*clhs32;
const double clhs266 = DN(1,0)*clhs50;
const double clhs267 = DN(1,1)*clhs57;
const double clhs268 = DN(1,0)*clhs92;
const double clhs269 = DN(1,1)*clhs94;
const double clhs270 = DN(1,0)*clhs69;
const double clhs271 = DN(1,1)*clhs69;
const double clhs272 = clhs37*(DN(1,0)*clhs111 + DN(1,1)*clhs115 - clhs263 - clhs265);
const double clhs273 = DN(1,0)*clhs124;
const double clhs274 = DN(1,1)*clhs126;
const double clhs275 = DN(1,0)*clhs148;
const double clhs276 = DN(1,1)*clhs150;
const double clhs277 = DN(1,0)*clhs170;
const double clhs278 = DN(1,1)*clhs172;
const double clhs279 = DN(1,0)*clhs194;
const double clhs280 = DN(1,1)*clhs196;
const double clhs281 = DN(1,0)*clhs42;
const double clhs282 = clhs281*clhs32;
const double clhs283 = DN(1,1)*clhs211;
const double clhs284 = clhs283*clhs32;
const double clhs285 = DN(1,0)*clhs57;
const double clhs286 = DN(1,1)*clhs215;
const double clhs287 = DN(1,0)*clhs94;
const double clhs288 = DN(1,1)*clhs220;
const double clhs289 = clhs37*(DN(1,0)*clhs226 + DN(1,1)*clhs227 - clhs282 - clhs284);
const double clhs290 = DN(1,0)*clhs126;
const double clhs291 = DN(1,1)*clhs230;
const double clhs292 = DN(1,0)*clhs150;
const double clhs293 = DN(1,1)*clhs235;
const double clhs294 = DN(1,0)*clhs172;
const double clhs295 = DN(1,1)*clhs240;
const double clhs296 = DN(1,0)*clhs196;
const double clhs297 = DN(1,1)*clhs245;
const double clhs298 = clhs251*(clhs262 + clhs264);
const double clhs299 = clhs251*(clhs281 + clhs283);
const double clhs300 = tau_u*(DN(1,0)*DN(2,1) - DN(1,1)*DN(2,0));
const double clhs301 = b_gauss[1]*clhs300;
const double clhs302 = b_gauss[0]*clhs300;
const double clhs303 = N[1]*N[2];
const double clhs304 = DN(2,0)*clhs30;
const double clhs305 = clhs304*clhs32;
const double clhs306 = DN(2,1)*clhs42;
const double clhs307 = clhs306*clhs32;
const double clhs308 = DN(2,0)*clhs50;
const double clhs309 = DN(2,1)*clhs57;
const double clhs310 = DN(2,0)*clhs92;
const double clhs311 = DN(2,1)*clhs94;
const double clhs312 = DN(2,0)*clhs69;
const double clhs313 = DN(2,1)*clhs69;
const double clhs314 = clhs37*(DN(2,0)*clhs111 + DN(2,1)*clhs115 - clhs305 - clhs307);
const double clhs315 = DN(2,0)*clhs124;
const double clhs316 = DN(2,1)*clhs126;
const double clhs317 = DN(2,0)*clhs148;
const double clhs318 = DN(2,1)*clhs150;
const double clhs319 = DN(2,0)*clhs170;
const double clhs320 = DN(2,1)*clhs172;
const double clhs321 = DN(2,0)*clhs194;
const double clhs322 = DN(2,1)*clhs196;
const double clhs323 = DN(2,0)*clhs42;
const double clhs324 = clhs32*clhs323;
const double clhs325 = DN(2,1)*clhs211;
const double clhs326 = clhs32*clhs325;
const double clhs327 = DN(2,0)*clhs57;
const double clhs328 = DN(2,1)*clhs215;
const double clhs329 = DN(2,0)*clhs94;
const double clhs330 = DN(2,1)*clhs220;
const double clhs331 = clhs37*(DN(2,0)*clhs226 + DN(2,1)*clhs227 - clhs324 - clhs326);
const double clhs332 = DN(2,0)*clhs126;
const double clhs333 = DN(2,1)*clhs230;
const double clhs334 = DN(2,0)*clhs150;
const double clhs335 = DN(2,1)*clhs235;
const double clhs336 = DN(2,0)*clhs172;
const double clhs337 = DN(2,1)*clhs240;
const double clhs338 = DN(2,0)*clhs196;
const double clhs339 = DN(2,1)*clhs245;
const double clhs340 = clhs251*(clhs304 + clhs306);
const double clhs341 = clhs251*(clhs323 + clhs325);
lhs(0,0)=DN(0,0)*clhs79 + DN(0,1)*clhs84 + clhs33*clhs38 - clhs33*clhs61 + clhs38*clhs44 - clhs44*clhs61 + clhs51*clhs56 + clhs56*clhs58;
lhs(0,1)=clhs104*clhs105 + clhs107*clhs108 + clhs33*clhs86 - clhs33*clhs96 + clhs44*clhs86 - clhs44*clhs96 + clhs56*clhs93 + clhs56*clhs95;
lhs(0,2)=N[0]*clhs116;
lhs(0,3)=DN(0,0)*clhs137 + DN(0,1)*clhs140 + clhs118*clhs33 + clhs118*clhs44 + clhs125*clhs56 + clhs127*clhs56 - clhs128*clhs33 - clhs128*clhs44;
lhs(0,4)=clhs105*clhs160 + clhs108*clhs162 + clhs142*clhs33 + clhs142*clhs44 + clhs149*clhs56 + clhs151*clhs56 - clhs152*clhs33 - clhs152*clhs44;
lhs(0,5)=N[1]*clhs116;
lhs(0,6)=DN(0,0)*clhs183 + DN(0,1)*clhs186 + clhs164*clhs33 + clhs164*clhs44 + clhs171*clhs56 + clhs173*clhs56 - clhs174*clhs33 - clhs174*clhs44;
lhs(0,7)=clhs105*clhs206 + clhs108*clhs208 + clhs188*clhs33 + clhs188*clhs44 + clhs195*clhs56 + clhs197*clhs56 - clhs198*clhs33 - clhs198*clhs44;
lhs(0,8)=N[2]*clhs116;
lhs(1,0)=clhs105*clhs217 + clhs108*clhs218 + clhs210*clhs38 - clhs210*clhs61 + clhs213*clhs38 - clhs213*clhs61 + clhs214*clhs56 + clhs216*clhs56;
lhs(1,1)=DN(0,0)*clhs224 + DN(0,1)*clhs225 + clhs210*clhs86 - clhs210*clhs96 + clhs213*clhs86 - clhs213*clhs96 + clhs219*clhs56 + clhs221*clhs56;
lhs(1,2)=N[0]*clhs228;
lhs(1,3)=clhs105*clhs232 + clhs108*clhs233 + clhs118*clhs210 + clhs118*clhs213 - clhs128*clhs210 - clhs128*clhs213 + clhs229*clhs56 + clhs231*clhs56;
lhs(1,4)=DN(0,0)*clhs237 + DN(0,1)*clhs238 + clhs142*clhs210 + clhs142*clhs213 - clhs152*clhs210 - clhs152*clhs213 + clhs234*clhs56 + clhs236*clhs56;
lhs(1,5)=N[1]*clhs228;
lhs(1,6)=clhs105*clhs242 + clhs108*clhs243 + clhs164*clhs210 + clhs164*clhs213 - clhs174*clhs210 - clhs174*clhs213 + clhs239*clhs56 + clhs241*clhs56;
lhs(1,7)=DN(0,0)*clhs247 + DN(0,1)*clhs248 + clhs188*clhs210 + clhs188*clhs213 - clhs198*clhs210 - clhs198*clhs213 + clhs244*clhs56 + clhs246*clhs56;
lhs(1,8)=N[2]*clhs228;
lhs(2,0)=N[0]*clhs34 - clhs249*(clhs51 + clhs58) - clhs250*(clhs214 + clhs216);
lhs(2,1)=N[0]*clhs85 - clhs249*(clhs93 + clhs95) - clhs250*(clhs219 + clhs221);
lhs(2,2)=-DN(0,0)*clhs252 - DN(0,1)*clhs253 - pow(N[0], 2);
lhs(2,3)=N[0]*clhs117 - clhs249*(clhs125 + clhs127) - clhs250*(clhs229 + clhs231) + clhs255;
lhs(2,4)=N[0]*clhs141 - clhs249*(clhs149 + clhs151) - clhs250*(clhs234 + clhs236) - clhs256;
lhs(2,5)=-DN(1,0)*clhs252 - DN(1,1)*clhs253 - clhs257;
lhs(2,6)=N[0]*clhs163 - clhs249*(clhs171 + clhs173) - clhs250*(clhs239 + clhs241) + clhs259;
lhs(2,7)=N[0]*clhs187 - clhs249*(clhs195 + clhs197) - clhs250*(clhs244 + clhs246) - clhs260;
lhs(2,8)=-DN(2,0)*clhs252 - DN(2,1)*clhs253 - clhs261;
lhs(3,0)=DN(1,0)*clhs79 + DN(1,1)*clhs84 + clhs263*clhs38 - clhs263*clhs61 + clhs265*clhs38 - clhs265*clhs61 + clhs266*clhs56 + clhs267*clhs56;
lhs(3,1)=clhs104*clhs270 + clhs107*clhs271 + clhs263*clhs86 - clhs263*clhs96 + clhs265*clhs86 - clhs265*clhs96 + clhs268*clhs56 + clhs269*clhs56;
lhs(3,2)=N[0]*clhs272;
lhs(3,3)=DN(1,0)*clhs137 + DN(1,1)*clhs140 + clhs118*clhs263 + clhs118*clhs265 - clhs128*clhs263 - clhs128*clhs265 + clhs273*clhs56 + clhs274*clhs56;
lhs(3,4)=clhs142*clhs263 + clhs142*clhs265 - clhs152*clhs263 - clhs152*clhs265 + clhs160*clhs270 + clhs162*clhs271 + clhs275*clhs56 + clhs276*clhs56;
lhs(3,5)=N[1]*clhs272;
lhs(3,6)=DN(1,0)*clhs183 + DN(1,1)*clhs186 + clhs164*clhs263 + clhs164*clhs265 - clhs174*clhs263 - clhs174*clhs265 + clhs277*clhs56 + clhs278*clhs56;
lhs(3,7)=clhs188*clhs263 + clhs188*clhs265 - clhs198*clhs263 - clhs198*clhs265 + clhs206*clhs270 + clhs208*clhs271 + clhs279*clhs56 + clhs280*clhs56;
lhs(3,8)=N[2]*clhs272;
lhs(4,0)=clhs217*clhs270 + clhs218*clhs271 + clhs282*clhs38 - clhs282*clhs61 + clhs284*clhs38 - clhs284*clhs61 + clhs285*clhs56 + clhs286*clhs56;
lhs(4,1)=DN(1,0)*clhs224 + DN(1,1)*clhs225 + clhs282*clhs86 - clhs282*clhs96 + clhs284*clhs86 - clhs284*clhs96 + clhs287*clhs56 + clhs288*clhs56;
lhs(4,2)=N[0]*clhs289;
lhs(4,3)=clhs118*clhs282 + clhs118*clhs284 - clhs128*clhs282 - clhs128*clhs284 + clhs232*clhs270 + clhs233*clhs271 + clhs290*clhs56 + clhs291*clhs56;
lhs(4,4)=DN(1,0)*clhs237 + DN(1,1)*clhs238 + clhs142*clhs282 + clhs142*clhs284 - clhs152*clhs282 - clhs152*clhs284 + clhs292*clhs56 + clhs293*clhs56;
lhs(4,5)=N[1]*clhs289;
lhs(4,6)=clhs164*clhs282 + clhs164*clhs284 - clhs174*clhs282 - clhs174*clhs284 + clhs242*clhs270 + clhs243*clhs271 + clhs294*clhs56 + clhs295*clhs56;
lhs(4,7)=DN(1,0)*clhs247 + DN(1,1)*clhs248 + clhs188*clhs282 + clhs188*clhs284 - clhs198*clhs282 - clhs198*clhs284 + clhs296*clhs56 + clhs297*clhs56;
lhs(4,8)=N[2]*clhs289;
lhs(5,0)=N[1]*clhs34 - clhs249*(clhs266 + clhs267) - clhs250*(clhs285 + clhs286) - clhs255;
lhs(5,1)=N[1]*clhs85 - clhs249*(clhs268 + clhs269) - clhs250*(clhs287 + clhs288) + clhs256;
lhs(5,2)=-DN(0,0)*clhs298 - DN(0,1)*clhs299 - clhs257;
lhs(5,3)=N[1]*clhs117 - clhs249*(clhs273 + clhs274) - clhs250*(clhs290 + clhs291);
lhs(5,4)=N[1]*clhs141 - clhs249*(clhs275 + clhs276) - clhs250*(clhs292 + clhs293);
lhs(5,5)=-DN(1,0)*clhs298 - DN(1,1)*clhs299 - pow(N[1], 2);
lhs(5,6)=N[1]*clhs163 - clhs249*(clhs277 + clhs278) - clhs250*(clhs294 + clhs295) + clhs301;
lhs(5,7)=N[1]*clhs187 - clhs249*(clhs279 + clhs280) - clhs250*(clhs296 + clhs297) - clhs302;
lhs(5,8)=-DN(2,0)*clhs298 - DN(2,1)*clhs299 - clhs303;
lhs(6,0)=DN(2,0)*clhs79 + DN(2,1)*clhs84 + clhs305*clhs38 - clhs305*clhs61 + clhs307*clhs38 - clhs307*clhs61 + clhs308*clhs56 + clhs309*clhs56;
lhs(6,1)=clhs104*clhs312 + clhs107*clhs313 + clhs305*clhs86 - clhs305*clhs96 + clhs307*clhs86 - clhs307*clhs96 + clhs310*clhs56 + clhs311*clhs56;
lhs(6,2)=N[0]*clhs314;
lhs(6,3)=DN(2,0)*clhs137 + DN(2,1)*clhs140 + clhs118*clhs305 + clhs118*clhs307 - clhs128*clhs305 - clhs128*clhs307 + clhs315*clhs56 + clhs316*clhs56;
lhs(6,4)=clhs142*clhs305 + clhs142*clhs307 - clhs152*clhs305 - clhs152*clhs307 + clhs160*clhs312 + clhs162*clhs313 + clhs317*clhs56 + clhs318*clhs56;
lhs(6,5)=N[1]*clhs314;
lhs(6,6)=DN(2,0)*clhs183 + DN(2,1)*clhs186 + clhs164*clhs305 + clhs164*clhs307 - clhs174*clhs305 - clhs174*clhs307 + clhs319*clhs56 + clhs320*clhs56;
lhs(6,7)=clhs188*clhs305 + clhs188*clhs307 - clhs198*clhs305 - clhs198*clhs307 + clhs206*clhs312 + clhs208*clhs313 + clhs321*clhs56 + clhs322*clhs56;
lhs(6,8)=N[2]*clhs314;
lhs(7,0)=clhs217*clhs312 + clhs218*clhs313 + clhs324*clhs38 - clhs324*clhs61 + clhs326*clhs38 - clhs326*clhs61 + clhs327*clhs56 + clhs328*clhs56;
lhs(7,1)=DN(2,0)*clhs224 + DN(2,1)*clhs225 + clhs324*clhs86 - clhs324*clhs96 + clhs326*clhs86 - clhs326*clhs96 + clhs329*clhs56 + clhs330*clhs56;
lhs(7,2)=N[0]*clhs331;
lhs(7,3)=clhs118*clhs324 + clhs118*clhs326 - clhs128*clhs324 - clhs128*clhs326 + clhs232*clhs312 + clhs233*clhs313 + clhs332*clhs56 + clhs333*clhs56;
lhs(7,4)=DN(2,0)*clhs237 + DN(2,1)*clhs238 + clhs142*clhs324 + clhs142*clhs326 - clhs152*clhs324 - clhs152*clhs326 + clhs334*clhs56 + clhs335*clhs56;
lhs(7,5)=N[1]*clhs331;
lhs(7,6)=clhs164*clhs324 + clhs164*clhs326 - clhs174*clhs324 - clhs174*clhs326 + clhs242*clhs312 + clhs243*clhs313 + clhs336*clhs56 + clhs337*clhs56;
lhs(7,7)=DN(2,0)*clhs247 + DN(2,1)*clhs248 + clhs188*clhs324 + clhs188*clhs326 - clhs198*clhs324 - clhs198*clhs326 + clhs338*clhs56 + clhs339*clhs56;
lhs(7,8)=N[2]*clhs331;
lhs(8,0)=N[2]*clhs34 - clhs249*(clhs308 + clhs309) - clhs250*(clhs327 + clhs328) - clhs259;
lhs(8,1)=N[2]*clhs85 - clhs249*(clhs310 + clhs311) - clhs250*(clhs329 + clhs330) + clhs260;
lhs(8,2)=-DN(0,0)*clhs340 - DN(0,1)*clhs341 - clhs261;
lhs(8,3)=N[2]*clhs117 - clhs249*(clhs315 + clhs316) - clhs250*(clhs332 + clhs333) - clhs301;
lhs(8,4)=N[2]*clhs141 - clhs249*(clhs317 + clhs318) - clhs250*(clhs334 + clhs335) + clhs302;
lhs(8,5)=-DN(1,0)*clhs340 - DN(1,1)*clhs341 - clhs303;
lhs(8,6)=N[2]*clhs163 - clhs249*(clhs319 + clhs320) - clhs250*(clhs336 + clhs337);
lhs(8,7)=N[2]*clhs187 - clhs249*(clhs321 + clhs322) - clhs250*(clhs338 + clhs339);
lhs(8,8)=-DN(2,0)*clhs340 - DN(2,1)*clhs341 - pow(N[2], 2);

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
const double crhs27 = 0.5*crhs26*tau_th*1.0/(crhs25 + crhs8);
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
