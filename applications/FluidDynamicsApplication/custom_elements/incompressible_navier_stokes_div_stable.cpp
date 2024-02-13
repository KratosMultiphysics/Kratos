//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "utilities/element_size_calculator.h"

// Application includes
#include "incompressible_navier_stokes_div_stable.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
IncompressibleNavierStokesDivStable<TDim>::IncompressibleNavierStokesDivStable(IndexType NewId)
    : Element(NewId)
{}

template< unsigned int TDim >
IncompressibleNavierStokesDivStable<TDim>::IncompressibleNavierStokesDivStable(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{}

template< unsigned int TDim >
IncompressibleNavierStokesDivStable<TDim>::IncompressibleNavierStokesDivStable(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{}

template< unsigned int TDim >
IncompressibleNavierStokesDivStable<TDim>::IncompressibleNavierStokesDivStable(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{}

template< unsigned int TDim >
IncompressibleNavierStokesDivStable<TDim>::~IncompressibleNavierStokesDivStable()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer IncompressibleNavierStokesDivStable<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<IncompressibleNavierStokesDivStable<TDim>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< unsigned int TDim >
Element::Pointer IncompressibleNavierStokesDivStable<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<IncompressibleNavierStokesDivStable<TDim>>(NewId, pGeom, pProperties);
}

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // If we are restarting, the constitutive law will be already defined
    if (mpConstitutiveLaw == nullptr) {
        const auto& r_properties = this->GetProperties();
        KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
            << "In initialization of Element " << this->Info()
            << ": No CONSTITUTIVE_LAW defined for property "
            << r_properties.Id() << "." << std::endl;

        mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();

        const auto& r_geometry = this->GetGeometry();
        const auto& r_shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);
        mpConstitutiveLaw->InitializeMaterial(r_properties, r_geometry, row(r_shape_functions,0));
    }

    KRATOS_CATCH("");
}

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rResult.size() != LocalSize) {
        rResult.resize(LocalSize, false);
    }

    IndexType local_index = 0;
    const auto& r_geometry = this->GetGeometry();
    const IndexType x_pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_X, x_pos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Y, x_pos+1).EquationId();
        if constexpr (TDim == 3) {
            rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Z, x_pos+2).EquationId();
        }
    }

    const IndexType p_pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(PRESSURE, p_pos).EquationId();
    }
}

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::GetDofList(
    DofsVectorType &rElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != LocalSize) {
        rElementalDofList.resize(LocalSize);

    }

    IndexType local_index = 0;
    const auto& r_geometry = this->GetGeometry();
    const IndexType x_pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_X, x_pos);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Y, x_pos+1);
        if constexpr (TDim == 3) {
            rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Z, x_pos+2);
        }
    }

    const IndexType p_pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(PRESSURE, p_pos);
    }
}

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Initialize element data
    ElementDataContainer aux_data;
    SetElementData(rCurrentProcessInfo, aux_data);

    // Initialize constitutive law parameters
    const auto& r_geom = this->GetGeometry();
    const auto p_prop = this->GetProperties();
    ConstitutiveLaw::Parameters cons_law_params(r_geom, p_prop, rCurrentProcessInfo);

    cons_law_params.SetStrainVector(aux_data.StrainRate);
    cons_law_params.SetStressVector(aux_data.ShearStress);
    cons_law_params.SetConstitutiveMatrix(aux_data.ConstitutiveMatrix);

    auto& cons_law_options = cons_law_params.GetOptions();
    cons_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS);
    cons_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Calculate kinematics
    Vector weights;
    Matrix velocity_N;
    Matrix pressure_N;
    GeometryType::ShapeFunctionsGradientsType velocity_DN;
    GeometryType::ShapeFunctionsGradientsType pressure_DN;
    Vector velocity_bubble;
    std::vector<BoundedMatrix<double, 1, TDim>> velocity_bubble_grad;
    DenseVector<GeometryType::ShapeFunctionsSecondDerivativesType> velocity_DDN;
    CalculateKinematics(weights, velocity_N, pressure_N, velocity_DN, pressure_DN, velocity_bubble, velocity_bubble_grad, velocity_DDN);

    // Initialize enrichment arrays
    array_1d<double, TDim> RHSee = ZeroVector(TDim);
    BoundedMatrix<double, TDim, TDim> Kee = ZeroMatrix(TDim, TDim);
    BoundedMatrix<double, TDim, LocalSize> Keu = ZeroMatrix(TDim, LocalSize);
    BoundedMatrix<double, LocalSize, TDim> Kue = ZeroMatrix(LocalSize, TDim);

    // Loop Gauss points
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(IntegrationMethod);
    for (IndexType g = 0; g < n_gauss; ++g) {
        // set current Gauss point kinematics
        noalias(aux_data.N_v) = row(velocity_N, g);
        noalias(aux_data.N_p) = row(pressure_N, g);
        noalias(aux_data.DN_v) = velocity_DN[g];
        noalias(aux_data.DN_p) = pressure_DN[g];
        aux_data.N_e = velocity_bubble[g];
        noalias(aux_data.DN_e) = velocity_bubble_grad[g];
        aux_data.DDN_v = velocity_DDN[g];
        aux_data.Weight = weights[g];

        // Calculate current Gauss point material response
        CalculateStrainRate(aux_data);
        mpConstitutiveLaw->CalculateMaterialResponseCauchy(cons_law_params);
        mpConstitutiveLaw->CalculateValue(cons_law_params, EFFECTIVE_VISCOSITY, aux_data.EffectiveViscosity);

        // Assemble standard Galerkin contribution
        ComputeGaussPointLHSContribution(aux_data, rLeftHandSideMatrix);
        ComputeGaussPointRHSContribution(aux_data, rRightHandSideVector);

        // Assemble bubble function contributions (to be condensed)
        ComputeGaussPointEnrichmentContribution(aux_data, RHSee, Kue, Keu, Kee);
    }

    // Condense bubble function contribution
    double det;
    BoundedMatrix<double, TDim, TDim> invKee;
    MathUtils<double>::InvertMatrix(Kee, invKee, det);
    const BoundedMatrix<double, LocalSize, TDim> KueinvKee = prod(Kue, invKee);
    noalias(rLeftHandSideMatrix) -= prod(KueinvKee, Keu);
    noalias(rRightHandSideVector) -= prod(KueinvKee, RHSee);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template< unsigned int TDim >
int IncompressibleNavierStokesDivStable<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int out = Element::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template< unsigned int TDim >
const Parameters IncompressibleNavierStokesDivStable<TDim>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [""],
            "nodal_historical"       : ["VELOCITY","PRESSURE"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","ACCELERATION","MESH_VELOCITY","PRESSURE","IS_STRUCTURE","DISPLACEMENT","BODY_FORCE","NODAL_AREA","NODAL_H","ADVPROJ","DIVPROJ","REACTION","REACTION_WATER_PRESSURE","EXTERNAL_PRESSURE","NORMAL","Y_WALL","Q_VALUE"]
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D6","Tetrahedra3D10"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["Newtonian2DLaw","Newtonian3DLaw","NewtonianTemperatureDependent2DLaw","NewtonianTemperatureDependent3DLaw","Euler2DLaw","Euler3DLaw"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 2,
        "documentation"   :
            "This implements a div-stable incompressible Navier-Stokes element with bubble function enrichment. No viscous behavior is hardcoded, meaning that any fluid constitutive model can be used through a constitutive law."
    })");

    if (TDim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template< unsigned int TDim >
std::string IncompressibleNavierStokesDivStable<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressibleNavierStokesDivStable" << TDim << "D" << VelocityNumNodes << "N #" << this->Id();
    return buffer.str();
}

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;

    if (mpConstitutiveLaw != nullptr) {
        rOStream << "with constitutive law " << std::endl;
        mpConstitutiveLaw->PrintInfo(rOStream);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private operations

template <unsigned int TDim>
void IncompressibleNavierStokesDivStable<TDim>::SetElementData(
    const ProcessInfo& rProcessInfo,
    ElementDataContainer &rData)
{
    // Set nodal data
    const auto& r_geom = this->GetGeometry();
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        const auto& r_v = r_geom[i].FastGetSolutionStepValue(VELOCITY);
        const auto& r_v_n = r_geom[i].FastGetSolutionStepValue(VELOCITY, 1);
        const auto& r_v_nn = r_geom[i].FastGetSolutionStepValue(VELOCITY, 2);
        const auto& r_v_mesh = r_geom[i].FastGetSolutionStepValue(MESH_VELOCITY);
        const auto& r_body_force = r_geom[i].FastGetSolutionStepValue(BODY_FORCE);

        for (IndexType d = 0; d < TDim; ++d) {
            rData.Velocity(i, d) = r_v[d];
            rData.VelocityOld1(i, d) = r_v_n[d];
            rData.VelocityOld2(i, d) = r_v_nn[d];
            rData.MeshVelocity(i, d) = r_v_mesh[d];
            rData.BodyForce(i, d) = r_body_force[d];
        }
    }

    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rData.Pressure[i] = r_geom[i].FastGetSolutionStepValue(PRESSURE);
    }

    // Set material values
    rData.Density = this->GetProperties().GetValue(DENSITY);

    // Set stabilization values
    rData.StabC1 = 2.0;
    rData.StabC2 = 1.0;
    rData.DynamicTau = rProcessInfo[DYNAMIC_TAU];
    rData.ElementSize = ElementSizeCalculator<TDim, VelocityNumNodes>::AverageElementSize(r_geom);

    // Set ProcessInfo data
    rData.DeltaTime = rProcessInfo[DELTA_TIME];
    const auto& r_bdf_coefs = rProcessInfo[BDF_COEFFICIENTS];
    rData.BDF0 = r_bdf_coefs[0];
    rData.BDF1 = r_bdf_coefs[1];
    rData.BDF2 = r_bdf_coefs[2];
}

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::CalculateKinematics(
    Vector& rGaussWeights,
    Matrix& rVelocityN,
    Matrix& rPressureN,
    GeometryType::ShapeFunctionsGradientsType& rVelocityDNDX,
    GeometryType::ShapeFunctionsGradientsType& rPressureDNDX,
    Vector& rVelocityBubble,
    std::vector<BoundedMatrix<double, 1, TDim>>& rVelocityBubbleGrad,
    DenseVector<GeometryType::ShapeFunctionsSecondDerivativesType>& rVelocityDDNDDX)
{
    // Get element geometry
    const auto& r_geom = this->GetGeometry();

    // Integration rule data
    // Note that we use the same for both velocity and pressure interpolations
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(IntegrationMethod);
    const auto integration_points = r_geom.IntegrationPoints(IntegrationMethod);

    // Calculate Jacobians at integration points
    Matrix J;
    Matrix inv_J;
    double det_J;
    Vector det_J_vect(n_gauss);
    std::vector<BoundedMatrix<double, TDim, TDim>> inv_J_vect(n_gauss);
    for (IndexType g = 0; g < n_gauss; ++g) {
        r_geom.Jacobian(J, g, IntegrationMethod);
        MathUtils<double>::InvertMatrix(J, inv_J, det_J);
        det_J_vect[g] = det_J;
        noalias(inv_J_vect[g]) = inv_J;
    }

    // Calculate velocity kinematics from the geometry (P2 interpolation)
    rVelocityN = r_geom.ShapeFunctionsValues(IntegrationMethod);
    const auto& r_DN_De_v = r_geom.ShapeFunctionsLocalGradients(IntegrationMethod);
    if (rVelocityDNDX.size() != n_gauss) {
        rVelocityDNDX.resize(n_gauss, false);
    }
    for (IndexType g = 0; g < n_gauss; ++g) {
        rVelocityDNDX[g] = prod(r_DN_De_v[g], inv_J_vect[g]);
    }
    GeometryUtils::ShapeFunctionsSecondDerivativesTransformOnAllIntegrationPoints(rVelocityDDNDDX, r_geom, IntegrationMethod);

    // Calculate pressure kinematics from an auxiliary geometry (P1 interpolation)
    GeometryType::UniquePointer p_aux_geom = nullptr;
    if constexpr (TDim == 2) {
        p_aux_geom = Kratos::make_unique<Triangle2D3<NodeType>>(r_geom(0), r_geom(1), r_geom(2));
    } else {
        p_aux_geom = Kratos::make_unique<Tetrahedra3D4<NodeType>>(r_geom(0), r_geom(1), r_geom(2), r_geom(3));
    }
    rPressureN = p_aux_geom->ShapeFunctionsValues(IntegrationMethod);
    if (rPressureDNDX.size() != n_gauss) {
        rPressureDNDX.resize(n_gauss, false);
    }
    const auto& r_DN_De_p = p_aux_geom->ShapeFunctionsLocalGradients(IntegrationMethod);
    for (IndexType g = 0; g < n_gauss; ++g) {
        rPressureDNDX[g] = prod(r_DN_De_p[g], inv_J_vect[g]);
    }

    // Calculate enrichment bubble kinematics using the previous auxiliary element barycentric coordinates
    if (rVelocityBubble.size() != n_gauss) {
        rVelocityBubble.resize(n_gauss, false);
    }
    if (rVelocityBubbleGrad.size() != n_gauss) {
        rVelocityBubbleGrad.resize(n_gauss);
    }
    for (IndexType g = 0; g < n_gauss; ++g) {
        const auto& r_bar_coords = row(rPressureN, g);
        const auto& r_bar_coords_grads = rPressureDNDX[g];
        auto& r_vel_bub_grad_g = rVelocityBubbleGrad[g];
        if constexpr (TDim == 2) {
            rVelocityBubble[g] = 27.0*(r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2]);
            r_vel_bub_grad_g(0,0) = 27.0*(r_bar_coords_grads(0,0) * r_bar_coords[1] * r_bar_coords[2] + r_bar_coords[0] * r_bar_coords_grads(1,0) * r_bar_coords[2] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,0));
            r_vel_bub_grad_g(0,1) = 27.0*(r_bar_coords_grads(0,1) * r_bar_coords[1] * r_bar_coords[2] + r_bar_coords[0] * r_bar_coords_grads(1,1) * r_bar_coords[2] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,1));
        } else {
            rVelocityBubble[g] = 256.0*(r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords[3]);
            r_vel_bub_grad_g(0,0) = 256.0*(r_bar_coords_grads(0,0) * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords_grads(1,0) * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,0) * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords_grads(3,0));
            r_vel_bub_grad_g(0,1) = 256.0*(r_bar_coords_grads(0,1) * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords_grads(1,1) * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,1) * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords_grads(3,1));
            r_vel_bub_grad_g(0,2) = 256.0*(r_bar_coords_grads(0,2) * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords_grads(1,2) * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,2) * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords_grads(3,2));
        }
    }

    // Calculate integration points weight
    if (rGaussWeights.size() != n_gauss) {
        rGaussWeights.resize(n_gauss, false);
    }
    for (IndexType g = 0; g < n_gauss; ++g) {
        rGaussWeights[g] = det_J_vect[g] * integration_points[g].Weight();
    }
}

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::CalculateStrainRate(ElementDataContainer& rData)
{
    //TODO: Add here the bubble contribution to the enrichment!!!

    if (rData.StrainRate.size() != StrainSize) {
        rData.StrainRate.resize(StrainSize, false);
    }
    noalias(rData.StrainRate) = ZeroVector(StrainSize);

    if constexpr (TDim == 2) {
        for (IndexType i = 0; i < VelocityNumNodes; ++i) {
            rData.StrainRate[0] += rData.DN_v(i,0)*rData.Velocity(i,0);
            rData.StrainRate[1] += rData.DN_v(i,1)*rData.Velocity(i,1);
            rData.StrainRate[2] += rData.DN_v(i,0)*rData.Velocity(i,1) + rData.DN_v(i,1)*rData.Velocity(i,0);
        }
    } else {
        for (IndexType i = 0; i < VelocityNumNodes; ++i) {
            rData.StrainRate[0] += rData.DN_v(i,0)*rData.Velocity(i,0);
            rData.StrainRate[1] += rData.DN_v(i,1)*rData.Velocity(i,1);
            rData.StrainRate[2] += rData.DN_v(i,2)*rData.Velocity(i,2);
            rData.StrainRate[3] += rData.DN_v(i,0)*rData.Velocity(i,1) + rData.DN_v(i,1)*rData.Velocity(i,0);
            rData.StrainRate[4] += rData.DN_v(i,1)*rData.Velocity(i,2) + rData.DN_v(i,2)*rData.Velocity(i,1);
            rData.StrainRate[5] += rData.DN_v(i,0)*rData.Velocity(i,2) + rData.DN_v(i,2)*rData.Velocity(i,0);
        }
    }

}

template <>
void IncompressibleNavierStokesDivStable<2>::ComputeGaussPointLHSContribution(
    const ElementDataContainer& rData,
    MatrixType& rLHS)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double dt = rData.DeltaTime;

    // Get stabilization data
    const double h = rData.ElementSize;
    const double stab_c1 = rData.StabC1;
    const double stab_c2 = rData.StabC2;
    const double dyn_tau = rData.DynamicTau;

    // Calculate convective velocity
    const BoundedMatrix<double,2,6> vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const auto& C = rData.ConstitutiveMatrix;

    // Get shape function values
    const auto& N_p = rData.N_p;
    const auto& N_v = rData.N_v;
    const auto& DN_p = rData.DN_p;
    const auto& DN_v = rData.DN_v;
    const auto& DDN_v = rData.DDN_v;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;
    const double crLHS0 = C(0,0)*DN_v(0,0) + C(0,2)*DN_v(0,1);
const double crLHS1 = C(0,2)*DN_v(0,0);
const double crLHS2 = C(2,2)*DN_v(0,1) + crLHS1;
const double crLHS3 = 1.0*C(0,0);
const double crLHS4 = C(0,2)*DDN_v[0](0,0);
const double crLHS5 = 1.0*DDN_v[0](0,1);
const double crLHS6 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crLHS7 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crLHS8 = rho*(DN_v(0,0)*crLHS6 + DN_v(0,1)*crLHS7);
const double crLHS9 = rData.BDF0*rho;
const double crLHS10 = N_v[0]*crLHS9;
const double crLHS11 = -crLHS10 - crLHS8;
const double crLHS12 = C(0,2)*crLHS5 + C(2,2)*crLHS5 + DDN_v[0](0,0)*crLHS3 + crLHS11 + 1.0*crLHS4;
const double crLHS13 = 1.0/(dyn_tau*rho/rData.DeltaTime + rho*stab_c2*sqrt(pow(crLHS6, 2) + pow(crLHS7, 2))/h + mu*stab_c1/pow(h, 2));
const double crLHS14 = crLHS13*crLHS8;
const double crLHS15 = rho*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1));
const double crLHS16 = N_v[0]*crLHS15;
const double crLHS17 = crLHS13*crLHS16;
const double crLHS18 = pow(N_v[0], 2)*crLHS9 + N_v[0]*crLHS8;
const double crLHS19 = C(0,1)*DN_v(0,1) + crLHS1;
const double crLHS20 = C(1,2)*DN_v(0,1);
const double crLHS21 = C(2,2)*DN_v(0,0) + crLHS20;
const double crLHS22 = C(0,1)*DDN_v[0](0,1) + C(1,2)*DDN_v[0](0,1) + C(2,2)*DDN_v[0](0,0) + crLHS4;
const double crLHS23 = C(0,0)*DN_v(1,0) + C(0,2)*DN_v(1,1);
const double crLHS24 = C(0,2)*DN_v(1,0);
const double crLHS25 = C(2,2)*DN_v(1,1) + crLHS24;
const double crLHS26 = C(0,2)*DDN_v[1](0,0);
const double crLHS27 = 1.0*DDN_v[1](0,1);
const double crLHS28 = rho*(DN_v(1,0)*crLHS6 + DN_v(1,1)*crLHS7);
const double crLHS29 = N_v[1]*crLHS9;
const double crLHS30 = -crLHS28 - crLHS29;
const double crLHS31 = C(0,2)*crLHS27 + C(2,2)*crLHS27 + DDN_v[1](0,0)*crLHS3 + 1.0*crLHS26 + crLHS30;
const double crLHS32 = N_v[1]*crLHS10;
const double crLHS33 = N_v[0]*crLHS28 + crLHS32;
const double crLHS34 = C(0,1)*DN_v(1,1) + crLHS24;
const double crLHS35 = C(1,2)*DN_v(1,1);
const double crLHS36 = C(2,2)*DN_v(1,0) + crLHS35;
const double crLHS37 = C(0,1)*DDN_v[1](0,1) + C(1,2)*DDN_v[1](0,1) + C(2,2)*DDN_v[1](0,0) + crLHS26;
const double crLHS38 = C(0,0)*DN_v(2,0) + C(0,2)*DN_v(2,1);
const double crLHS39 = C(0,2)*DN_v(2,0);
const double crLHS40 = C(2,2)*DN_v(2,1) + crLHS39;
const double crLHS41 = C(0,2)*DDN_v[2](0,0);
const double crLHS42 = 1.0*DDN_v[2](0,1);
const double crLHS43 = rho*(DN_v(2,0)*crLHS6 + DN_v(2,1)*crLHS7);
const double crLHS44 = N_v[2]*crLHS9;
const double crLHS45 = -crLHS43 - crLHS44;
const double crLHS46 = C(0,2)*crLHS42 + C(2,2)*crLHS42 + DDN_v[2](0,0)*crLHS3 + 1.0*crLHS41 + crLHS45;
const double crLHS47 = N_v[2]*crLHS10;
const double crLHS48 = N_v[0]*crLHS43 + crLHS47;
const double crLHS49 = C(0,1)*DN_v(2,1) + crLHS39;
const double crLHS50 = C(1,2)*DN_v(2,1);
const double crLHS51 = C(2,2)*DN_v(2,0) + crLHS50;
const double crLHS52 = C(0,1)*DDN_v[2](0,1) + C(1,2)*DDN_v[2](0,1) + C(2,2)*DDN_v[2](0,0) + crLHS41;
const double crLHS53 = C(0,0)*DN_v(3,0) + C(0,2)*DN_v(3,1);
const double crLHS54 = C(0,2)*DN_v(3,0);
const double crLHS55 = C(2,2)*DN_v(3,1) + crLHS54;
const double crLHS56 = C(0,2)*DDN_v[3](0,0);
const double crLHS57 = 1.0*DDN_v[3](0,1);
const double crLHS58 = rho*(DN_v(3,0)*crLHS6 + DN_v(3,1)*crLHS7);
const double crLHS59 = N_v[3]*crLHS9;
const double crLHS60 = -crLHS58 - crLHS59;
const double crLHS61 = C(0,2)*crLHS57 + C(2,2)*crLHS57 + DDN_v[3](0,0)*crLHS3 + 1.0*crLHS56 + crLHS60;
const double crLHS62 = N_v[3]*crLHS10;
const double crLHS63 = N_v[0]*crLHS58 + crLHS62;
const double crLHS64 = C(0,1)*DN_v(3,1) + crLHS54;
const double crLHS65 = C(1,2)*DN_v(3,1);
const double crLHS66 = C(2,2)*DN_v(3,0) + crLHS65;
const double crLHS67 = C(0,1)*DDN_v[3](0,1) + C(1,2)*DDN_v[3](0,1) + C(2,2)*DDN_v[3](0,0) + crLHS56;
const double crLHS68 = C(0,0)*DN_v(4,0) + C(0,2)*DN_v(4,1);
const double crLHS69 = C(0,2)*DN_v(4,0);
const double crLHS70 = C(2,2)*DN_v(4,1) + crLHS69;
const double crLHS71 = C(0,2)*DDN_v[4](0,0);
const double crLHS72 = 1.0*DDN_v[4](0,1);
const double crLHS73 = rho*(DN_v(4,0)*crLHS6 + DN_v(4,1)*crLHS7);
const double crLHS74 = N_v[4]*crLHS9;
const double crLHS75 = -crLHS73 - crLHS74;
const double crLHS76 = C(0,2)*crLHS72 + C(2,2)*crLHS72 + DDN_v[4](0,0)*crLHS3 + 1.0*crLHS71 + crLHS75;
const double crLHS77 = N_v[4]*crLHS10;
const double crLHS78 = N_v[0]*crLHS73 + crLHS77;
const double crLHS79 = C(0,1)*DN_v(4,1) + crLHS69;
const double crLHS80 = C(1,2)*DN_v(4,1);
const double crLHS81 = C(2,2)*DN_v(4,0) + crLHS80;
const double crLHS82 = C(0,1)*DDN_v[4](0,1) + C(1,2)*DDN_v[4](0,1) + C(2,2)*DDN_v[4](0,0) + crLHS71;
const double crLHS83 = C(0,0)*DN_v(5,0) + C(0,2)*DN_v(5,1);
const double crLHS84 = C(0,2)*DN_v(5,0);
const double crLHS85 = C(2,2)*DN_v(5,1) + crLHS84;
const double crLHS86 = C(0,2)*DDN_v[5](0,0);
const double crLHS87 = 1.0*DDN_v[5](0,1);
const double crLHS88 = rho*(DN_v(5,0)*crLHS6 + DN_v(5,1)*crLHS7);
const double crLHS89 = -N_v[5]*crLHS9 - crLHS88;
const double crLHS90 = C(0,2)*crLHS87 + C(2,2)*crLHS87 + DDN_v[5](0,0)*crLHS3 + 1.0*crLHS86 + crLHS89;
const double crLHS91 = N_v[5]*crLHS10;
const double crLHS92 = N_v[0]*crLHS88 + crLHS91;
const double crLHS93 = C(0,1)*DN_v(5,1) + crLHS84;
const double crLHS94 = C(1,2)*DN_v(5,1);
const double crLHS95 = C(2,2)*DN_v(5,0) + crLHS94;
const double crLHS96 = C(0,1)*DDN_v[5](0,1) + C(1,2)*DDN_v[5](0,1) + C(2,2)*DDN_v[5](0,0) + crLHS86;
const double crLHS97 = -DN_v(0,0)*N_p[0];
const double crLHS98 = DN_p(0,0)*crLHS13;
const double crLHS99 = -DN_v(0,0)*N_p[1];
const double crLHS100 = DN_p(1,0)*crLHS13;
const double crLHS101 = -DN_v(0,0)*N_p[2];
const double crLHS102 = DN_p(2,0)*crLHS13;
const double crLHS103 = C(0,1)*DN_v(0,0) + crLHS20;
const double crLHS104 = C(1,2)*DDN_v[0](1,1);
const double crLHS105 = C(0,1)*DDN_v[0](1,0) + C(0,2)*DDN_v[0](1,0) + C(2,2)*DDN_v[0](1,1) + crLHS104;
const double crLHS106 = C(1,1)*DN_v(0,1) + C(1,2)*DN_v(0,0);
const double crLHS107 = 1.0*C(1,1);
const double crLHS108 = 1.0*DDN_v[0](1,0);
const double crLHS109 = C(1,2)*crLHS108 + C(2,2)*crLHS108 + DDN_v[0](1,1)*crLHS107 + 1.0*crLHS104 + crLHS11;
const double crLHS110 = C(0,1)*DN_v(1,0) + crLHS35;
const double crLHS111 = C(1,2)*DDN_v[1](1,1);
const double crLHS112 = C(0,1)*DDN_v[1](1,0) + C(0,2)*DDN_v[1](1,0) + C(2,2)*DDN_v[1](1,1) + crLHS111;
const double crLHS113 = C(1,1)*DN_v(1,1) + C(1,2)*DN_v(1,0);
const double crLHS114 = 1.0*DDN_v[1](1,0);
const double crLHS115 = C(1,2)*crLHS114 + C(2,2)*crLHS114 + DDN_v[1](1,1)*crLHS107 + 1.0*crLHS111 + crLHS30;
const double crLHS116 = C(0,1)*DN_v(2,0) + crLHS50;
const double crLHS117 = C(1,2)*DDN_v[2](1,1);
const double crLHS118 = C(0,1)*DDN_v[2](1,0) + C(0,2)*DDN_v[2](1,0) + C(2,2)*DDN_v[2](1,1) + crLHS117;
const double crLHS119 = C(1,1)*DN_v(2,1) + C(1,2)*DN_v(2,0);
const double crLHS120 = 1.0*DDN_v[2](1,0);
const double crLHS121 = C(1,2)*crLHS120 + C(2,2)*crLHS120 + DDN_v[2](1,1)*crLHS107 + 1.0*crLHS117 + crLHS45;
const double crLHS122 = C(0,1)*DN_v(3,0) + crLHS65;
const double crLHS123 = C(1,2)*DDN_v[3](1,1);
const double crLHS124 = C(0,1)*DDN_v[3](1,0) + C(0,2)*DDN_v[3](1,0) + C(2,2)*DDN_v[3](1,1) + crLHS123;
const double crLHS125 = C(1,1)*DN_v(3,1) + C(1,2)*DN_v(3,0);
const double crLHS126 = 1.0*DDN_v[3](1,0);
const double crLHS127 = C(1,2)*crLHS126 + C(2,2)*crLHS126 + DDN_v[3](1,1)*crLHS107 + 1.0*crLHS123 + crLHS60;
const double crLHS128 = C(0,1)*DN_v(4,0) + crLHS80;
const double crLHS129 = C(1,2)*DDN_v[4](1,1);
const double crLHS130 = C(0,1)*DDN_v[4](1,0) + C(0,2)*DDN_v[4](1,0) + C(2,2)*DDN_v[4](1,1) + crLHS129;
const double crLHS131 = C(1,1)*DN_v(4,1) + C(1,2)*DN_v(4,0);
const double crLHS132 = 1.0*DDN_v[4](1,0);
const double crLHS133 = C(1,2)*crLHS132 + C(2,2)*crLHS132 + DDN_v[4](1,1)*crLHS107 + 1.0*crLHS129 + crLHS75;
const double crLHS134 = C(0,1)*DN_v(5,0) + crLHS94;
const double crLHS135 = C(1,2)*DDN_v[5](1,1);
const double crLHS136 = C(0,1)*DDN_v[5](1,0) + C(0,2)*DDN_v[5](1,0) + C(2,2)*DDN_v[5](1,1) + crLHS135;
const double crLHS137 = C(1,1)*DN_v(5,1) + C(1,2)*DN_v(5,0);
const double crLHS138 = 1.0*DDN_v[5](1,0);
const double crLHS139 = C(1,2)*crLHS138 + C(2,2)*crLHS138 + DDN_v[5](1,1)*crLHS107 + 1.0*crLHS135 + crLHS89;
const double crLHS140 = -DN_v(0,1)*N_p[0];
const double crLHS141 = DN_p(0,1)*crLHS13;
const double crLHS142 = -DN_v(0,1)*N_p[1];
const double crLHS143 = DN_p(1,1)*crLHS13;
const double crLHS144 = -DN_v(0,1)*N_p[2];
const double crLHS145 = DN_p(2,1)*crLHS13;
const double crLHS146 = crLHS13*crLHS28;
const double crLHS147 = N_v[1]*crLHS15;
const double crLHS148 = crLHS13*crLHS147;
const double crLHS149 = N_v[1]*crLHS8 + crLHS32;
const double crLHS150 = pow(N_v[1], 2)*crLHS9 + N_v[1]*crLHS28;
const double crLHS151 = N_v[2]*crLHS29;
const double crLHS152 = N_v[1]*crLHS43 + crLHS151;
const double crLHS153 = N_v[3]*crLHS29;
const double crLHS154 = N_v[1]*crLHS58 + crLHS153;
const double crLHS155 = N_v[4]*crLHS29;
const double crLHS156 = N_v[1]*crLHS73 + crLHS155;
const double crLHS157 = N_v[5]*crLHS29;
const double crLHS158 = N_v[1]*crLHS88 + crLHS157;
const double crLHS159 = -DN_v(1,0)*N_p[0];
const double crLHS160 = -DN_v(1,0)*N_p[1];
const double crLHS161 = -DN_v(1,0)*N_p[2];
const double crLHS162 = -DN_v(1,1)*N_p[0];
const double crLHS163 = -DN_v(1,1)*N_p[1];
const double crLHS164 = -DN_v(1,1)*N_p[2];
const double crLHS165 = crLHS13*crLHS43;
const double crLHS166 = N_v[2]*crLHS15;
const double crLHS167 = crLHS13*crLHS166;
const double crLHS168 = N_v[2]*crLHS8 + crLHS47;
const double crLHS169 = N_v[2]*crLHS28 + crLHS151;
const double crLHS170 = pow(N_v[2], 2)*crLHS9 + N_v[2]*crLHS43;
const double crLHS171 = N_v[3]*crLHS44;
const double crLHS172 = N_v[2]*crLHS58 + crLHS171;
const double crLHS173 = N_v[4]*crLHS44;
const double crLHS174 = N_v[2]*crLHS73 + crLHS173;
const double crLHS175 = N_v[5]*crLHS44;
const double crLHS176 = N_v[2]*crLHS88 + crLHS175;
const double crLHS177 = -DN_v(2,0)*N_p[0];
const double crLHS178 = -DN_v(2,0)*N_p[1];
const double crLHS179 = -DN_v(2,0)*N_p[2];
const double crLHS180 = -DN_v(2,1)*N_p[0];
const double crLHS181 = -DN_v(2,1)*N_p[1];
const double crLHS182 = -DN_v(2,1)*N_p[2];
const double crLHS183 = crLHS13*crLHS58;
const double crLHS184 = N_v[3]*crLHS15;
const double crLHS185 = crLHS13*crLHS184;
const double crLHS186 = N_v[3]*crLHS8 + crLHS62;
const double crLHS187 = N_v[3]*crLHS28 + crLHS153;
const double crLHS188 = N_v[3]*crLHS43 + crLHS171;
const double crLHS189 = pow(N_v[3], 2)*crLHS9 + N_v[3]*crLHS58;
const double crLHS190 = N_v[4]*crLHS59;
const double crLHS191 = N_v[3]*crLHS73 + crLHS190;
const double crLHS192 = N_v[5]*crLHS59;
const double crLHS193 = N_v[3]*crLHS88 + crLHS192;
const double crLHS194 = -DN_v(3,0)*N_p[0];
const double crLHS195 = -DN_v(3,0)*N_p[1];
const double crLHS196 = -DN_v(3,0)*N_p[2];
const double crLHS197 = -DN_v(3,1)*N_p[0];
const double crLHS198 = -DN_v(3,1)*N_p[1];
const double crLHS199 = -DN_v(3,1)*N_p[2];
const double crLHS200 = crLHS13*crLHS73;
const double crLHS201 = N_v[4]*crLHS15;
const double crLHS202 = crLHS13*crLHS201;
const double crLHS203 = N_v[4]*crLHS8 + crLHS77;
const double crLHS204 = N_v[4]*crLHS28 + crLHS155;
const double crLHS205 = N_v[4]*crLHS43 + crLHS173;
const double crLHS206 = N_v[4]*crLHS58 + crLHS190;
const double crLHS207 = pow(N_v[4], 2)*crLHS9 + N_v[4]*crLHS73;
const double crLHS208 = N_v[5]*crLHS74;
const double crLHS209 = N_v[4]*crLHS88 + crLHS208;
const double crLHS210 = -DN_v(4,0)*N_p[0];
const double crLHS211 = -DN_v(4,0)*N_p[1];
const double crLHS212 = -DN_v(4,0)*N_p[2];
const double crLHS213 = -DN_v(4,1)*N_p[0];
const double crLHS214 = -DN_v(4,1)*N_p[1];
const double crLHS215 = -DN_v(4,1)*N_p[2];
const double crLHS216 = crLHS13*crLHS88;
const double crLHS217 = N_v[5]*crLHS15;
const double crLHS218 = crLHS13*crLHS217;
const double crLHS219 = N_v[5]*crLHS8 + crLHS91;
const double crLHS220 = N_v[5]*crLHS28 + crLHS157;
const double crLHS221 = N_v[5]*crLHS43 + crLHS175;
const double crLHS222 = N_v[5]*crLHS58 + crLHS192;
const double crLHS223 = N_v[5]*crLHS73 + crLHS208;
const double crLHS224 = pow(N_v[5], 2)*crLHS9 + N_v[5]*crLHS88;
const double crLHS225 = -DN_v(5,0)*N_p[0];
const double crLHS226 = -DN_v(5,0)*N_p[1];
const double crLHS227 = -DN_v(5,0)*N_p[2];
const double crLHS228 = -DN_v(5,1)*N_p[0];
const double crLHS229 = -DN_v(5,1)*N_p[1];
const double crLHS230 = -DN_v(5,1)*N_p[2];
const double crLHS231 = crLHS13*gauss_weight;
const double crLHS232 = crLHS231*(DN_p(0,0)*DN_p(1,0) + DN_p(0,1)*DN_p(1,1));
const double crLHS233 = crLHS231*(DN_p(0,0)*DN_p(2,0) + DN_p(0,1)*DN_p(2,1));
const double crLHS234 = crLHS231*(DN_p(1,0)*DN_p(2,0) + DN_p(1,1)*DN_p(2,1));
rLHS(0,0)+=gauss_weight*(DN_v(0,0)*crLHS0 + DN_v(0,1)*crLHS2 - crLHS12*crLHS14 - crLHS12*crLHS17 + crLHS18);
rLHS(0,1)+=gauss_weight*(DN_v(0,0)*crLHS19 + DN_v(0,1)*crLHS21 - crLHS14*crLHS22 - crLHS17*crLHS22);
rLHS(0,2)+=gauss_weight*(DN_v(0,0)*crLHS23 + DN_v(0,1)*crLHS25 - crLHS14*crLHS31 - crLHS17*crLHS31 + crLHS33);
rLHS(0,3)+=gauss_weight*(DN_v(0,0)*crLHS34 + DN_v(0,1)*crLHS36 - crLHS14*crLHS37 - crLHS17*crLHS37);
rLHS(0,4)+=gauss_weight*(DN_v(0,0)*crLHS38 + DN_v(0,1)*crLHS40 - crLHS14*crLHS46 - crLHS17*crLHS46 + crLHS48);
rLHS(0,5)+=gauss_weight*(DN_v(0,0)*crLHS49 + DN_v(0,1)*crLHS51 - crLHS14*crLHS52 - crLHS17*crLHS52);
rLHS(0,6)+=gauss_weight*(DN_v(0,0)*crLHS53 + DN_v(0,1)*crLHS55 - crLHS14*crLHS61 - crLHS17*crLHS61 + crLHS63);
rLHS(0,7)+=gauss_weight*(DN_v(0,0)*crLHS64 + DN_v(0,1)*crLHS66 - crLHS14*crLHS67 - crLHS17*crLHS67);
rLHS(0,8)+=gauss_weight*(DN_v(0,0)*crLHS68 + DN_v(0,1)*crLHS70 - crLHS14*crLHS76 - crLHS17*crLHS76 + crLHS78);
rLHS(0,9)+=gauss_weight*(DN_v(0,0)*crLHS79 + DN_v(0,1)*crLHS81 - crLHS14*crLHS82 - crLHS17*crLHS82);
rLHS(0,10)+=gauss_weight*(DN_v(0,0)*crLHS83 + DN_v(0,1)*crLHS85 - crLHS14*crLHS90 - crLHS17*crLHS90 + crLHS92);
rLHS(0,11)+=gauss_weight*(DN_v(0,0)*crLHS93 + DN_v(0,1)*crLHS95 - crLHS14*crLHS96 - crLHS17*crLHS96);
rLHS(0,12)+=gauss_weight*(crLHS16*crLHS98 + crLHS8*crLHS98 + crLHS97);
rLHS(0,13)+=gauss_weight*(crLHS100*crLHS16 + crLHS100*crLHS8 + crLHS99);
rLHS(0,14)+=gauss_weight*(crLHS101 + crLHS102*crLHS16 + crLHS102*crLHS8);
rLHS(1,0)+=gauss_weight*(DN_v(0,0)*crLHS2 + DN_v(0,1)*crLHS103 - crLHS105*crLHS14 - crLHS105*crLHS17);
rLHS(1,1)+=gauss_weight*(DN_v(0,0)*crLHS21 + DN_v(0,1)*crLHS106 - crLHS109*crLHS14 - crLHS109*crLHS17 + crLHS18);
rLHS(1,2)+=gauss_weight*(DN_v(0,0)*crLHS25 + DN_v(0,1)*crLHS110 - crLHS112*crLHS14 - crLHS112*crLHS17);
rLHS(1,3)+=gauss_weight*(DN_v(0,0)*crLHS36 + DN_v(0,1)*crLHS113 - crLHS115*crLHS14 - crLHS115*crLHS17 + crLHS33);
rLHS(1,4)+=gauss_weight*(DN_v(0,0)*crLHS40 + DN_v(0,1)*crLHS116 - crLHS118*crLHS14 - crLHS118*crLHS17);
rLHS(1,5)+=gauss_weight*(DN_v(0,0)*crLHS51 + DN_v(0,1)*crLHS119 - crLHS121*crLHS14 - crLHS121*crLHS17 + crLHS48);
rLHS(1,6)+=gauss_weight*(DN_v(0,0)*crLHS55 + DN_v(0,1)*crLHS122 - crLHS124*crLHS14 - crLHS124*crLHS17);
rLHS(1,7)+=gauss_weight*(DN_v(0,0)*crLHS66 + DN_v(0,1)*crLHS125 - crLHS127*crLHS14 - crLHS127*crLHS17 + crLHS63);
rLHS(1,8)+=gauss_weight*(DN_v(0,0)*crLHS70 + DN_v(0,1)*crLHS128 - crLHS130*crLHS14 - crLHS130*crLHS17);
rLHS(1,9)+=gauss_weight*(DN_v(0,0)*crLHS81 + DN_v(0,1)*crLHS131 - crLHS133*crLHS14 - crLHS133*crLHS17 + crLHS78);
rLHS(1,10)+=gauss_weight*(DN_v(0,0)*crLHS85 + DN_v(0,1)*crLHS134 - crLHS136*crLHS14 - crLHS136*crLHS17);
rLHS(1,11)+=gauss_weight*(DN_v(0,0)*crLHS95 + DN_v(0,1)*crLHS137 - crLHS139*crLHS14 - crLHS139*crLHS17 + crLHS92);
rLHS(1,12)+=gauss_weight*(crLHS140 + crLHS141*crLHS16 + crLHS141*crLHS8);
rLHS(1,13)+=gauss_weight*(crLHS142 + crLHS143*crLHS16 + crLHS143*crLHS8);
rLHS(1,14)+=gauss_weight*(crLHS144 + crLHS145*crLHS16 + crLHS145*crLHS8);
rLHS(2,0)+=gauss_weight*(DN_v(1,0)*crLHS0 + DN_v(1,1)*crLHS2 - crLHS12*crLHS146 - crLHS12*crLHS148 + crLHS149);
rLHS(2,1)+=gauss_weight*(DN_v(1,0)*crLHS19 + DN_v(1,1)*crLHS21 - crLHS146*crLHS22 - crLHS148*crLHS22);
rLHS(2,2)+=gauss_weight*(DN_v(1,0)*crLHS23 + DN_v(1,1)*crLHS25 - crLHS146*crLHS31 - crLHS148*crLHS31 + crLHS150);
rLHS(2,3)+=gauss_weight*(DN_v(1,0)*crLHS34 + DN_v(1,1)*crLHS36 - crLHS146*crLHS37 - crLHS148*crLHS37);
rLHS(2,4)+=gauss_weight*(DN_v(1,0)*crLHS38 + DN_v(1,1)*crLHS40 - crLHS146*crLHS46 - crLHS148*crLHS46 + crLHS152);
rLHS(2,5)+=gauss_weight*(DN_v(1,0)*crLHS49 + DN_v(1,1)*crLHS51 - crLHS146*crLHS52 - crLHS148*crLHS52);
rLHS(2,6)+=gauss_weight*(DN_v(1,0)*crLHS53 + DN_v(1,1)*crLHS55 - crLHS146*crLHS61 - crLHS148*crLHS61 + crLHS154);
rLHS(2,7)+=gauss_weight*(DN_v(1,0)*crLHS64 + DN_v(1,1)*crLHS66 - crLHS146*crLHS67 - crLHS148*crLHS67);
rLHS(2,8)+=gauss_weight*(DN_v(1,0)*crLHS68 + DN_v(1,1)*crLHS70 - crLHS146*crLHS76 - crLHS148*crLHS76 + crLHS156);
rLHS(2,9)+=gauss_weight*(DN_v(1,0)*crLHS79 + DN_v(1,1)*crLHS81 - crLHS146*crLHS82 - crLHS148*crLHS82);
rLHS(2,10)+=gauss_weight*(DN_v(1,0)*crLHS83 + DN_v(1,1)*crLHS85 - crLHS146*crLHS90 - crLHS148*crLHS90 + crLHS158);
rLHS(2,11)+=gauss_weight*(DN_v(1,0)*crLHS93 + DN_v(1,1)*crLHS95 - crLHS146*crLHS96 - crLHS148*crLHS96);
rLHS(2,12)+=gauss_weight*(crLHS147*crLHS98 + crLHS159 + crLHS28*crLHS98);
rLHS(2,13)+=gauss_weight*(crLHS100*crLHS147 + crLHS100*crLHS28 + crLHS160);
rLHS(2,14)+=gauss_weight*(crLHS102*crLHS147 + crLHS102*crLHS28 + crLHS161);
rLHS(3,0)+=gauss_weight*(DN_v(1,0)*crLHS2 + DN_v(1,1)*crLHS103 - crLHS105*crLHS146 - crLHS105*crLHS148);
rLHS(3,1)+=gauss_weight*(DN_v(1,0)*crLHS21 + DN_v(1,1)*crLHS106 - crLHS109*crLHS146 - crLHS109*crLHS148 + crLHS149);
rLHS(3,2)+=gauss_weight*(DN_v(1,0)*crLHS25 + DN_v(1,1)*crLHS110 - crLHS112*crLHS146 - crLHS112*crLHS148);
rLHS(3,3)+=gauss_weight*(DN_v(1,0)*crLHS36 + DN_v(1,1)*crLHS113 - crLHS115*crLHS146 - crLHS115*crLHS148 + crLHS150);
rLHS(3,4)+=gauss_weight*(DN_v(1,0)*crLHS40 + DN_v(1,1)*crLHS116 - crLHS118*crLHS146 - crLHS118*crLHS148);
rLHS(3,5)+=gauss_weight*(DN_v(1,0)*crLHS51 + DN_v(1,1)*crLHS119 - crLHS121*crLHS146 - crLHS121*crLHS148 + crLHS152);
rLHS(3,6)+=gauss_weight*(DN_v(1,0)*crLHS55 + DN_v(1,1)*crLHS122 - crLHS124*crLHS146 - crLHS124*crLHS148);
rLHS(3,7)+=gauss_weight*(DN_v(1,0)*crLHS66 + DN_v(1,1)*crLHS125 - crLHS127*crLHS146 - crLHS127*crLHS148 + crLHS154);
rLHS(3,8)+=gauss_weight*(DN_v(1,0)*crLHS70 + DN_v(1,1)*crLHS128 - crLHS130*crLHS146 - crLHS130*crLHS148);
rLHS(3,9)+=gauss_weight*(DN_v(1,0)*crLHS81 + DN_v(1,1)*crLHS131 - crLHS133*crLHS146 - crLHS133*crLHS148 + crLHS156);
rLHS(3,10)+=gauss_weight*(DN_v(1,0)*crLHS85 + DN_v(1,1)*crLHS134 - crLHS136*crLHS146 - crLHS136*crLHS148);
rLHS(3,11)+=gauss_weight*(DN_v(1,0)*crLHS95 + DN_v(1,1)*crLHS137 - crLHS139*crLHS146 - crLHS139*crLHS148 + crLHS158);
rLHS(3,12)+=gauss_weight*(crLHS141*crLHS147 + crLHS141*crLHS28 + crLHS162);
rLHS(3,13)+=gauss_weight*(crLHS143*crLHS147 + crLHS143*crLHS28 + crLHS163);
rLHS(3,14)+=gauss_weight*(crLHS145*crLHS147 + crLHS145*crLHS28 + crLHS164);
rLHS(4,0)+=gauss_weight*(DN_v(2,0)*crLHS0 + DN_v(2,1)*crLHS2 - crLHS12*crLHS165 - crLHS12*crLHS167 + crLHS168);
rLHS(4,1)+=gauss_weight*(DN_v(2,0)*crLHS19 + DN_v(2,1)*crLHS21 - crLHS165*crLHS22 - crLHS167*crLHS22);
rLHS(4,2)+=gauss_weight*(DN_v(2,0)*crLHS23 + DN_v(2,1)*crLHS25 - crLHS165*crLHS31 - crLHS167*crLHS31 + crLHS169);
rLHS(4,3)+=gauss_weight*(DN_v(2,0)*crLHS34 + DN_v(2,1)*crLHS36 - crLHS165*crLHS37 - crLHS167*crLHS37);
rLHS(4,4)+=gauss_weight*(DN_v(2,0)*crLHS38 + DN_v(2,1)*crLHS40 - crLHS165*crLHS46 - crLHS167*crLHS46 + crLHS170);
rLHS(4,5)+=gauss_weight*(DN_v(2,0)*crLHS49 + DN_v(2,1)*crLHS51 - crLHS165*crLHS52 - crLHS167*crLHS52);
rLHS(4,6)+=gauss_weight*(DN_v(2,0)*crLHS53 + DN_v(2,1)*crLHS55 - crLHS165*crLHS61 - crLHS167*crLHS61 + crLHS172);
rLHS(4,7)+=gauss_weight*(DN_v(2,0)*crLHS64 + DN_v(2,1)*crLHS66 - crLHS165*crLHS67 - crLHS167*crLHS67);
rLHS(4,8)+=gauss_weight*(DN_v(2,0)*crLHS68 + DN_v(2,1)*crLHS70 - crLHS165*crLHS76 - crLHS167*crLHS76 + crLHS174);
rLHS(4,9)+=gauss_weight*(DN_v(2,0)*crLHS79 + DN_v(2,1)*crLHS81 - crLHS165*crLHS82 - crLHS167*crLHS82);
rLHS(4,10)+=gauss_weight*(DN_v(2,0)*crLHS83 + DN_v(2,1)*crLHS85 - crLHS165*crLHS90 - crLHS167*crLHS90 + crLHS176);
rLHS(4,11)+=gauss_weight*(DN_v(2,0)*crLHS93 + DN_v(2,1)*crLHS95 - crLHS165*crLHS96 - crLHS167*crLHS96);
rLHS(4,12)+=gauss_weight*(crLHS166*crLHS98 + crLHS177 + crLHS43*crLHS98);
rLHS(4,13)+=gauss_weight*(crLHS100*crLHS166 + crLHS100*crLHS43 + crLHS178);
rLHS(4,14)+=gauss_weight*(crLHS102*crLHS166 + crLHS102*crLHS43 + crLHS179);
rLHS(5,0)+=gauss_weight*(DN_v(2,0)*crLHS2 + DN_v(2,1)*crLHS103 - crLHS105*crLHS165 - crLHS105*crLHS167);
rLHS(5,1)+=gauss_weight*(DN_v(2,0)*crLHS21 + DN_v(2,1)*crLHS106 - crLHS109*crLHS165 - crLHS109*crLHS167 + crLHS168);
rLHS(5,2)+=gauss_weight*(DN_v(2,0)*crLHS25 + DN_v(2,1)*crLHS110 - crLHS112*crLHS165 - crLHS112*crLHS167);
rLHS(5,3)+=gauss_weight*(DN_v(2,0)*crLHS36 + DN_v(2,1)*crLHS113 - crLHS115*crLHS165 - crLHS115*crLHS167 + crLHS169);
rLHS(5,4)+=gauss_weight*(DN_v(2,0)*crLHS40 + DN_v(2,1)*crLHS116 - crLHS118*crLHS165 - crLHS118*crLHS167);
rLHS(5,5)+=gauss_weight*(DN_v(2,0)*crLHS51 + DN_v(2,1)*crLHS119 - crLHS121*crLHS165 - crLHS121*crLHS167 + crLHS170);
rLHS(5,6)+=gauss_weight*(DN_v(2,0)*crLHS55 + DN_v(2,1)*crLHS122 - crLHS124*crLHS165 - crLHS124*crLHS167);
rLHS(5,7)+=gauss_weight*(DN_v(2,0)*crLHS66 + DN_v(2,1)*crLHS125 - crLHS127*crLHS165 - crLHS127*crLHS167 + crLHS172);
rLHS(5,8)+=gauss_weight*(DN_v(2,0)*crLHS70 + DN_v(2,1)*crLHS128 - crLHS130*crLHS165 - crLHS130*crLHS167);
rLHS(5,9)+=gauss_weight*(DN_v(2,0)*crLHS81 + DN_v(2,1)*crLHS131 - crLHS133*crLHS165 - crLHS133*crLHS167 + crLHS174);
rLHS(5,10)+=gauss_weight*(DN_v(2,0)*crLHS85 + DN_v(2,1)*crLHS134 - crLHS136*crLHS165 - crLHS136*crLHS167);
rLHS(5,11)+=gauss_weight*(DN_v(2,0)*crLHS95 + DN_v(2,1)*crLHS137 - crLHS139*crLHS165 - crLHS139*crLHS167 + crLHS176);
rLHS(5,12)+=gauss_weight*(crLHS141*crLHS166 + crLHS141*crLHS43 + crLHS180);
rLHS(5,13)+=gauss_weight*(crLHS143*crLHS166 + crLHS143*crLHS43 + crLHS181);
rLHS(5,14)+=gauss_weight*(crLHS145*crLHS166 + crLHS145*crLHS43 + crLHS182);
rLHS(6,0)+=gauss_weight*(DN_v(3,0)*crLHS0 + DN_v(3,1)*crLHS2 - crLHS12*crLHS183 - crLHS12*crLHS185 + crLHS186);
rLHS(6,1)+=gauss_weight*(DN_v(3,0)*crLHS19 + DN_v(3,1)*crLHS21 - crLHS183*crLHS22 - crLHS185*crLHS22);
rLHS(6,2)+=gauss_weight*(DN_v(3,0)*crLHS23 + DN_v(3,1)*crLHS25 - crLHS183*crLHS31 - crLHS185*crLHS31 + crLHS187);
rLHS(6,3)+=gauss_weight*(DN_v(3,0)*crLHS34 + DN_v(3,1)*crLHS36 - crLHS183*crLHS37 - crLHS185*crLHS37);
rLHS(6,4)+=gauss_weight*(DN_v(3,0)*crLHS38 + DN_v(3,1)*crLHS40 - crLHS183*crLHS46 - crLHS185*crLHS46 + crLHS188);
rLHS(6,5)+=gauss_weight*(DN_v(3,0)*crLHS49 + DN_v(3,1)*crLHS51 - crLHS183*crLHS52 - crLHS185*crLHS52);
rLHS(6,6)+=gauss_weight*(DN_v(3,0)*crLHS53 + DN_v(3,1)*crLHS55 - crLHS183*crLHS61 - crLHS185*crLHS61 + crLHS189);
rLHS(6,7)+=gauss_weight*(DN_v(3,0)*crLHS64 + DN_v(3,1)*crLHS66 - crLHS183*crLHS67 - crLHS185*crLHS67);
rLHS(6,8)+=gauss_weight*(DN_v(3,0)*crLHS68 + DN_v(3,1)*crLHS70 - crLHS183*crLHS76 - crLHS185*crLHS76 + crLHS191);
rLHS(6,9)+=gauss_weight*(DN_v(3,0)*crLHS79 + DN_v(3,1)*crLHS81 - crLHS183*crLHS82 - crLHS185*crLHS82);
rLHS(6,10)+=gauss_weight*(DN_v(3,0)*crLHS83 + DN_v(3,1)*crLHS85 - crLHS183*crLHS90 - crLHS185*crLHS90 + crLHS193);
rLHS(6,11)+=gauss_weight*(DN_v(3,0)*crLHS93 + DN_v(3,1)*crLHS95 - crLHS183*crLHS96 - crLHS185*crLHS96);
rLHS(6,12)+=gauss_weight*(crLHS184*crLHS98 + crLHS194 + crLHS58*crLHS98);
rLHS(6,13)+=gauss_weight*(crLHS100*crLHS184 + crLHS100*crLHS58 + crLHS195);
rLHS(6,14)+=gauss_weight*(crLHS102*crLHS184 + crLHS102*crLHS58 + crLHS196);
rLHS(7,0)+=gauss_weight*(DN_v(3,0)*crLHS2 + DN_v(3,1)*crLHS103 - crLHS105*crLHS183 - crLHS105*crLHS185);
rLHS(7,1)+=gauss_weight*(DN_v(3,0)*crLHS21 + DN_v(3,1)*crLHS106 - crLHS109*crLHS183 - crLHS109*crLHS185 + crLHS186);
rLHS(7,2)+=gauss_weight*(DN_v(3,0)*crLHS25 + DN_v(3,1)*crLHS110 - crLHS112*crLHS183 - crLHS112*crLHS185);
rLHS(7,3)+=gauss_weight*(DN_v(3,0)*crLHS36 + DN_v(3,1)*crLHS113 - crLHS115*crLHS183 - crLHS115*crLHS185 + crLHS187);
rLHS(7,4)+=gauss_weight*(DN_v(3,0)*crLHS40 + DN_v(3,1)*crLHS116 - crLHS118*crLHS183 - crLHS118*crLHS185);
rLHS(7,5)+=gauss_weight*(DN_v(3,0)*crLHS51 + DN_v(3,1)*crLHS119 - crLHS121*crLHS183 - crLHS121*crLHS185 + crLHS188);
rLHS(7,6)+=gauss_weight*(DN_v(3,0)*crLHS55 + DN_v(3,1)*crLHS122 - crLHS124*crLHS183 - crLHS124*crLHS185);
rLHS(7,7)+=gauss_weight*(DN_v(3,0)*crLHS66 + DN_v(3,1)*crLHS125 - crLHS127*crLHS183 - crLHS127*crLHS185 + crLHS189);
rLHS(7,8)+=gauss_weight*(DN_v(3,0)*crLHS70 + DN_v(3,1)*crLHS128 - crLHS130*crLHS183 - crLHS130*crLHS185);
rLHS(7,9)+=gauss_weight*(DN_v(3,0)*crLHS81 + DN_v(3,1)*crLHS131 - crLHS133*crLHS183 - crLHS133*crLHS185 + crLHS191);
rLHS(7,10)+=gauss_weight*(DN_v(3,0)*crLHS85 + DN_v(3,1)*crLHS134 - crLHS136*crLHS183 - crLHS136*crLHS185);
rLHS(7,11)+=gauss_weight*(DN_v(3,0)*crLHS95 + DN_v(3,1)*crLHS137 - crLHS139*crLHS183 - crLHS139*crLHS185 + crLHS193);
rLHS(7,12)+=gauss_weight*(crLHS141*crLHS184 + crLHS141*crLHS58 + crLHS197);
rLHS(7,13)+=gauss_weight*(crLHS143*crLHS184 + crLHS143*crLHS58 + crLHS198);
rLHS(7,14)+=gauss_weight*(crLHS145*crLHS184 + crLHS145*crLHS58 + crLHS199);
rLHS(8,0)+=gauss_weight*(DN_v(4,0)*crLHS0 + DN_v(4,1)*crLHS2 - crLHS12*crLHS200 - crLHS12*crLHS202 + crLHS203);
rLHS(8,1)+=gauss_weight*(DN_v(4,0)*crLHS19 + DN_v(4,1)*crLHS21 - crLHS200*crLHS22 - crLHS202*crLHS22);
rLHS(8,2)+=gauss_weight*(DN_v(4,0)*crLHS23 + DN_v(4,1)*crLHS25 - crLHS200*crLHS31 - crLHS202*crLHS31 + crLHS204);
rLHS(8,3)+=gauss_weight*(DN_v(4,0)*crLHS34 + DN_v(4,1)*crLHS36 - crLHS200*crLHS37 - crLHS202*crLHS37);
rLHS(8,4)+=gauss_weight*(DN_v(4,0)*crLHS38 + DN_v(4,1)*crLHS40 - crLHS200*crLHS46 - crLHS202*crLHS46 + crLHS205);
rLHS(8,5)+=gauss_weight*(DN_v(4,0)*crLHS49 + DN_v(4,1)*crLHS51 - crLHS200*crLHS52 - crLHS202*crLHS52);
rLHS(8,6)+=gauss_weight*(DN_v(4,0)*crLHS53 + DN_v(4,1)*crLHS55 - crLHS200*crLHS61 - crLHS202*crLHS61 + crLHS206);
rLHS(8,7)+=gauss_weight*(DN_v(4,0)*crLHS64 + DN_v(4,1)*crLHS66 - crLHS200*crLHS67 - crLHS202*crLHS67);
rLHS(8,8)+=gauss_weight*(DN_v(4,0)*crLHS68 + DN_v(4,1)*crLHS70 - crLHS200*crLHS76 - crLHS202*crLHS76 + crLHS207);
rLHS(8,9)+=gauss_weight*(DN_v(4,0)*crLHS79 + DN_v(4,1)*crLHS81 - crLHS200*crLHS82 - crLHS202*crLHS82);
rLHS(8,10)+=gauss_weight*(DN_v(4,0)*crLHS83 + DN_v(4,1)*crLHS85 - crLHS200*crLHS90 - crLHS202*crLHS90 + crLHS209);
rLHS(8,11)+=gauss_weight*(DN_v(4,0)*crLHS93 + DN_v(4,1)*crLHS95 - crLHS200*crLHS96 - crLHS202*crLHS96);
rLHS(8,12)+=gauss_weight*(crLHS201*crLHS98 + crLHS210 + crLHS73*crLHS98);
rLHS(8,13)+=gauss_weight*(crLHS100*crLHS201 + crLHS100*crLHS73 + crLHS211);
rLHS(8,14)+=gauss_weight*(crLHS102*crLHS201 + crLHS102*crLHS73 + crLHS212);
rLHS(9,0)+=gauss_weight*(DN_v(4,0)*crLHS2 + DN_v(4,1)*crLHS103 - crLHS105*crLHS200 - crLHS105*crLHS202);
rLHS(9,1)+=gauss_weight*(DN_v(4,0)*crLHS21 + DN_v(4,1)*crLHS106 - crLHS109*crLHS200 - crLHS109*crLHS202 + crLHS203);
rLHS(9,2)+=gauss_weight*(DN_v(4,0)*crLHS25 + DN_v(4,1)*crLHS110 - crLHS112*crLHS200 - crLHS112*crLHS202);
rLHS(9,3)+=gauss_weight*(DN_v(4,0)*crLHS36 + DN_v(4,1)*crLHS113 - crLHS115*crLHS200 - crLHS115*crLHS202 + crLHS204);
rLHS(9,4)+=gauss_weight*(DN_v(4,0)*crLHS40 + DN_v(4,1)*crLHS116 - crLHS118*crLHS200 - crLHS118*crLHS202);
rLHS(9,5)+=gauss_weight*(DN_v(4,0)*crLHS51 + DN_v(4,1)*crLHS119 - crLHS121*crLHS200 - crLHS121*crLHS202 + crLHS205);
rLHS(9,6)+=gauss_weight*(DN_v(4,0)*crLHS55 + DN_v(4,1)*crLHS122 - crLHS124*crLHS200 - crLHS124*crLHS202);
rLHS(9,7)+=gauss_weight*(DN_v(4,0)*crLHS66 + DN_v(4,1)*crLHS125 - crLHS127*crLHS200 - crLHS127*crLHS202 + crLHS206);
rLHS(9,8)+=gauss_weight*(DN_v(4,0)*crLHS70 + DN_v(4,1)*crLHS128 - crLHS130*crLHS200 - crLHS130*crLHS202);
rLHS(9,9)+=gauss_weight*(DN_v(4,0)*crLHS81 + DN_v(4,1)*crLHS131 - crLHS133*crLHS200 - crLHS133*crLHS202 + crLHS207);
rLHS(9,10)+=gauss_weight*(DN_v(4,0)*crLHS85 + DN_v(4,1)*crLHS134 - crLHS136*crLHS200 - crLHS136*crLHS202);
rLHS(9,11)+=gauss_weight*(DN_v(4,0)*crLHS95 + DN_v(4,1)*crLHS137 - crLHS139*crLHS200 - crLHS139*crLHS202 + crLHS209);
rLHS(9,12)+=gauss_weight*(crLHS141*crLHS201 + crLHS141*crLHS73 + crLHS213);
rLHS(9,13)+=gauss_weight*(crLHS143*crLHS201 + crLHS143*crLHS73 + crLHS214);
rLHS(9,14)+=gauss_weight*(crLHS145*crLHS201 + crLHS145*crLHS73 + crLHS215);
rLHS(10,0)+=gauss_weight*(DN_v(5,0)*crLHS0 + DN_v(5,1)*crLHS2 - crLHS12*crLHS216 - crLHS12*crLHS218 + crLHS219);
rLHS(10,1)+=gauss_weight*(DN_v(5,0)*crLHS19 + DN_v(5,1)*crLHS21 - crLHS216*crLHS22 - crLHS218*crLHS22);
rLHS(10,2)+=gauss_weight*(DN_v(5,0)*crLHS23 + DN_v(5,1)*crLHS25 - crLHS216*crLHS31 - crLHS218*crLHS31 + crLHS220);
rLHS(10,3)+=gauss_weight*(DN_v(5,0)*crLHS34 + DN_v(5,1)*crLHS36 - crLHS216*crLHS37 - crLHS218*crLHS37);
rLHS(10,4)+=gauss_weight*(DN_v(5,0)*crLHS38 + DN_v(5,1)*crLHS40 - crLHS216*crLHS46 - crLHS218*crLHS46 + crLHS221);
rLHS(10,5)+=gauss_weight*(DN_v(5,0)*crLHS49 + DN_v(5,1)*crLHS51 - crLHS216*crLHS52 - crLHS218*crLHS52);
rLHS(10,6)+=gauss_weight*(DN_v(5,0)*crLHS53 + DN_v(5,1)*crLHS55 - crLHS216*crLHS61 - crLHS218*crLHS61 + crLHS222);
rLHS(10,7)+=gauss_weight*(DN_v(5,0)*crLHS64 + DN_v(5,1)*crLHS66 - crLHS216*crLHS67 - crLHS218*crLHS67);
rLHS(10,8)+=gauss_weight*(DN_v(5,0)*crLHS68 + DN_v(5,1)*crLHS70 - crLHS216*crLHS76 - crLHS218*crLHS76 + crLHS223);
rLHS(10,9)+=gauss_weight*(DN_v(5,0)*crLHS79 + DN_v(5,1)*crLHS81 - crLHS216*crLHS82 - crLHS218*crLHS82);
rLHS(10,10)+=gauss_weight*(DN_v(5,0)*crLHS83 + DN_v(5,1)*crLHS85 - crLHS216*crLHS90 - crLHS218*crLHS90 + crLHS224);
rLHS(10,11)+=gauss_weight*(DN_v(5,0)*crLHS93 + DN_v(5,1)*crLHS95 - crLHS216*crLHS96 - crLHS218*crLHS96);
rLHS(10,12)+=gauss_weight*(crLHS217*crLHS98 + crLHS225 + crLHS88*crLHS98);
rLHS(10,13)+=gauss_weight*(crLHS100*crLHS217 + crLHS100*crLHS88 + crLHS226);
rLHS(10,14)+=gauss_weight*(crLHS102*crLHS217 + crLHS102*crLHS88 + crLHS227);
rLHS(11,0)+=gauss_weight*(DN_v(5,0)*crLHS2 + DN_v(5,1)*crLHS103 - crLHS105*crLHS216 - crLHS105*crLHS218);
rLHS(11,1)+=gauss_weight*(DN_v(5,0)*crLHS21 + DN_v(5,1)*crLHS106 - crLHS109*crLHS216 - crLHS109*crLHS218 + crLHS219);
rLHS(11,2)+=gauss_weight*(DN_v(5,0)*crLHS25 + DN_v(5,1)*crLHS110 - crLHS112*crLHS216 - crLHS112*crLHS218);
rLHS(11,3)+=gauss_weight*(DN_v(5,0)*crLHS36 + DN_v(5,1)*crLHS113 - crLHS115*crLHS216 - crLHS115*crLHS218 + crLHS220);
rLHS(11,4)+=gauss_weight*(DN_v(5,0)*crLHS40 + DN_v(5,1)*crLHS116 - crLHS118*crLHS216 - crLHS118*crLHS218);
rLHS(11,5)+=gauss_weight*(DN_v(5,0)*crLHS51 + DN_v(5,1)*crLHS119 - crLHS121*crLHS216 - crLHS121*crLHS218 + crLHS221);
rLHS(11,6)+=gauss_weight*(DN_v(5,0)*crLHS55 + DN_v(5,1)*crLHS122 - crLHS124*crLHS216 - crLHS124*crLHS218);
rLHS(11,7)+=gauss_weight*(DN_v(5,0)*crLHS66 + DN_v(5,1)*crLHS125 - crLHS127*crLHS216 - crLHS127*crLHS218 + crLHS222);
rLHS(11,8)+=gauss_weight*(DN_v(5,0)*crLHS70 + DN_v(5,1)*crLHS128 - crLHS130*crLHS216 - crLHS130*crLHS218);
rLHS(11,9)+=gauss_weight*(DN_v(5,0)*crLHS81 + DN_v(5,1)*crLHS131 - crLHS133*crLHS216 - crLHS133*crLHS218 + crLHS223);
rLHS(11,10)+=gauss_weight*(DN_v(5,0)*crLHS85 + DN_v(5,1)*crLHS134 - crLHS136*crLHS216 - crLHS136*crLHS218);
rLHS(11,11)+=gauss_weight*(DN_v(5,0)*crLHS95 + DN_v(5,1)*crLHS137 - crLHS139*crLHS216 - crLHS139*crLHS218 + crLHS224);
rLHS(11,12)+=gauss_weight*(crLHS141*crLHS217 + crLHS141*crLHS88 + crLHS228);
rLHS(11,13)+=gauss_weight*(crLHS143*crLHS217 + crLHS143*crLHS88 + crLHS229);
rLHS(11,14)+=gauss_weight*(crLHS145*crLHS217 + crLHS145*crLHS88 + crLHS230);
rLHS(12,0)+=-gauss_weight*(crLHS105*crLHS141 + crLHS12*crLHS98 + crLHS97);
rLHS(12,1)+=-gauss_weight*(crLHS109*crLHS141 + crLHS140 + crLHS22*crLHS98);
rLHS(12,2)+=-gauss_weight*(crLHS112*crLHS141 + crLHS159 + crLHS31*crLHS98);
rLHS(12,3)+=-gauss_weight*(crLHS115*crLHS141 + crLHS162 + crLHS37*crLHS98);
rLHS(12,4)+=-gauss_weight*(crLHS118*crLHS141 + crLHS177 + crLHS46*crLHS98);
rLHS(12,5)+=-gauss_weight*(crLHS121*crLHS141 + crLHS180 + crLHS52*crLHS98);
rLHS(12,6)+=-gauss_weight*(crLHS124*crLHS141 + crLHS194 + crLHS61*crLHS98);
rLHS(12,7)+=-gauss_weight*(crLHS127*crLHS141 + crLHS197 + crLHS67*crLHS98);
rLHS(12,8)+=-gauss_weight*(crLHS130*crLHS141 + crLHS210 + crLHS76*crLHS98);
rLHS(12,9)+=-gauss_weight*(crLHS133*crLHS141 + crLHS213 + crLHS82*crLHS98);
rLHS(12,10)+=-gauss_weight*(crLHS136*crLHS141 + crLHS225 + crLHS90*crLHS98);
rLHS(12,11)+=-gauss_weight*(crLHS139*crLHS141 + crLHS228 + crLHS96*crLHS98);
rLHS(12,12)+=crLHS231*(pow(DN_p(0,0), 2) + pow(DN_p(0,1), 2));
rLHS(12,13)+=crLHS232;
rLHS(12,14)+=crLHS233;
rLHS(13,0)+=-gauss_weight*(crLHS100*crLHS12 + crLHS105*crLHS143 + crLHS99);
rLHS(13,1)+=-gauss_weight*(crLHS100*crLHS22 + crLHS109*crLHS143 + crLHS142);
rLHS(13,2)+=-gauss_weight*(crLHS100*crLHS31 + crLHS112*crLHS143 + crLHS160);
rLHS(13,3)+=-gauss_weight*(crLHS100*crLHS37 + crLHS115*crLHS143 + crLHS163);
rLHS(13,4)+=-gauss_weight*(crLHS100*crLHS46 + crLHS118*crLHS143 + crLHS178);
rLHS(13,5)+=-gauss_weight*(crLHS100*crLHS52 + crLHS121*crLHS143 + crLHS181);
rLHS(13,6)+=-gauss_weight*(crLHS100*crLHS61 + crLHS124*crLHS143 + crLHS195);
rLHS(13,7)+=-gauss_weight*(crLHS100*crLHS67 + crLHS127*crLHS143 + crLHS198);
rLHS(13,8)+=-gauss_weight*(crLHS100*crLHS76 + crLHS130*crLHS143 + crLHS211);
rLHS(13,9)+=-gauss_weight*(crLHS100*crLHS82 + crLHS133*crLHS143 + crLHS214);
rLHS(13,10)+=-gauss_weight*(crLHS100*crLHS90 + crLHS136*crLHS143 + crLHS226);
rLHS(13,11)+=-gauss_weight*(crLHS100*crLHS96 + crLHS139*crLHS143 + crLHS229);
rLHS(13,12)+=crLHS232;
rLHS(13,13)+=crLHS231*(pow(DN_p(1,0), 2) + pow(DN_p(1,1), 2));
rLHS(13,14)+=crLHS234;
rLHS(14,0)+=-gauss_weight*(crLHS101 + crLHS102*crLHS12 + crLHS105*crLHS145);
rLHS(14,1)+=-gauss_weight*(crLHS102*crLHS22 + crLHS109*crLHS145 + crLHS144);
rLHS(14,2)+=-gauss_weight*(crLHS102*crLHS31 + crLHS112*crLHS145 + crLHS161);
rLHS(14,3)+=-gauss_weight*(crLHS102*crLHS37 + crLHS115*crLHS145 + crLHS164);
rLHS(14,4)+=-gauss_weight*(crLHS102*crLHS46 + crLHS118*crLHS145 + crLHS179);
rLHS(14,5)+=-gauss_weight*(crLHS102*crLHS52 + crLHS121*crLHS145 + crLHS182);
rLHS(14,6)+=-gauss_weight*(crLHS102*crLHS61 + crLHS124*crLHS145 + crLHS196);
rLHS(14,7)+=-gauss_weight*(crLHS102*crLHS67 + crLHS127*crLHS145 + crLHS199);
rLHS(14,8)+=-gauss_weight*(crLHS102*crLHS76 + crLHS130*crLHS145 + crLHS212);
rLHS(14,9)+=-gauss_weight*(crLHS102*crLHS82 + crLHS133*crLHS145 + crLHS215);
rLHS(14,10)+=-gauss_weight*(crLHS102*crLHS90 + crLHS136*crLHS145 + crLHS227);
rLHS(14,11)+=-gauss_weight*(crLHS102*crLHS96 + crLHS139*crLHS145 + crLHS230);
rLHS(14,12)+=crLHS233;
rLHS(14,13)+=crLHS234;
rLHS(14,14)+=crLHS231*(pow(DN_p(2,0), 2) + pow(DN_p(2,1), 2));

}

template <>
void IncompressibleNavierStokesDivStable<3>::ComputeGaussPointLHSContribution(
    const ElementDataContainer& rData,
    MatrixType& rLHS)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double dt = rData.DeltaTime;

    // Get stabilization data
    const double h = rData.ElementSize;
    const double stab_c1 = rData.StabC1;
    const double stab_c2 = rData.StabC2;
    const double dyn_tau = rData.DynamicTau;

    // Calculate convective velocity
    const BoundedMatrix<double,3,10> vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const auto& C = rData.ConstitutiveMatrix;

    // Get shape function values
    const auto& N_p = rData.N_p;
    const auto& N_v = rData.N_v;
    const auto& DN_p = rData.DN_p;
    const auto& DN_v = rData.DN_v;
    const auto& DDN_v = rData.DDN_v;

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;
    //substitute_lhs_3D
}

template <>
void IncompressibleNavierStokesDivStable<2>::ComputeGaussPointRHSContribution(
    const ElementDataContainer& rData,
    VectorType& rRHS)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double bdf1 = rData.BDF1;
    const double bdf2 = rData.BDF2;
    const double dt = rData.DeltaTime;

    // Get stabilization data
    const double h = rData.ElementSize;
    const double stab_c1 = rData.StabC1;
    const double stab_c2 = rData.StabC2;
    const double dyn_tau = rData.DynamicTau;

    // Get nodal data
    const auto& r_v = rData.Velocity;
    const auto& r_vn = rData.VelocityOld1;
    const auto& r_vnn = rData.VelocityOld2;
    const auto& r_vmesh = rData.MeshVelocity;
    const auto& r_f = rData.BodyForce;
    const auto& r_p = rData.Pressure;

    // Calculate convective velocity
    const BoundedMatrix<double,2,6> vconv = r_v - r_vmesh;

    // Get stress from material response
    const auto& r_stress = rData.ShearStress;
    const auto& C = rData.ConstitutiveMatrix;

    // Get shape function values
    const auto& N_p = rData.N_p;
    const auto& N_v = rData.N_v;
    const auto& DN_p = rData.DN_p;
    const auto& DN_v = rData.DN_v;
    const auto& DDN_v = rData.DDN_v;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;
    const double crRHS0 = N_p[0]*r_p[0] + N_p[1]*r_p[1] + N_p[2]*r_p[2];
const double crRHS1 = rho*(N_v[0]*r_f(0,0) + N_v[1]*r_f(1,0) + N_v[2]*r_f(2,0) + N_v[3]*r_f(3,0) + N_v[4]*r_f(4,0) + N_v[5]*r_f(5,0));
const double crRHS2 = rho*(N_v[0]*(rData.BDF0*r_v(0,0) + rData.BDF1*r_vn(0,0) + rData.BDF2*r_vnn(0,0)) + N_v[1]*(rData.BDF0*r_v(1,0) + rData.BDF1*r_vn(1,0) + rData.BDF2*r_vnn(1,0)) + N_v[2]*(rData.BDF0*r_v(2,0) + rData.BDF1*r_vn(2,0) + rData.BDF2*r_vnn(2,0)) + N_v[3]*(rData.BDF0*r_v(3,0) + rData.BDF1*r_vn(3,0) + rData.BDF2*r_vnn(3,0)) + N_v[4]*(rData.BDF0*r_v(4,0) + rData.BDF1*r_vn(4,0) + rData.BDF2*r_vnn(4,0)) + N_v[5]*(rData.BDF0*r_v(5,0) + rData.BDF1*r_vn(5,0) + rData.BDF2*r_vnn(5,0)));
const double crRHS3 = DN_v(0,0)*r_v(0,0) + DN_v(1,0)*r_v(1,0) + DN_v(2,0)*r_v(2,0) + DN_v(3,0)*r_v(3,0) + DN_v(4,0)*r_v(4,0) + DN_v(5,0)*r_v(5,0);
const double crRHS4 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crRHS5 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crRHS6 = rho*(crRHS3*crRHS4 + crRHS5*(DN_v(0,1)*r_v(0,0) + DN_v(1,1)*r_v(1,0) + DN_v(2,1)*r_v(2,0) + DN_v(3,1)*r_v(3,0) + DN_v(4,1)*r_v(4,0) + DN_v(5,1)*r_v(5,0)));
const double crRHS7 = 1.0*C(0,0);
const double crRHS8 = DDN_v[0](0,0)*r_v(0,0);
const double crRHS9 = DDN_v[1](0,0)*r_v(1,0);
const double crRHS10 = DDN_v[2](0,0)*r_v(2,0);
const double crRHS11 = DDN_v[3](0,0)*r_v(3,0);
const double crRHS12 = DDN_v[4](0,0)*r_v(4,0);
const double crRHS13 = DDN_v[5](0,0)*r_v(5,0);
const double crRHS14 = 1.0*C(0,1);
const double crRHS15 = DDN_v[0](0,1)*r_v(0,1);
const double crRHS16 = DDN_v[1](0,1)*r_v(1,1);
const double crRHS17 = DDN_v[2](0,1)*r_v(2,1);
const double crRHS18 = DDN_v[3](0,1)*r_v(3,1);
const double crRHS19 = DDN_v[4](0,1)*r_v(4,1);
const double crRHS20 = DDN_v[5](0,1)*r_v(5,1);
const double crRHS21 = 1.0*C(0,2);
const double crRHS22 = 1.0*C(1,2);
const double crRHS23 = DDN_v[0](0,0)*r_v(0,1) + DDN_v[0](0,1)*r_v(0,0);
const double crRHS24 = DDN_v[1](0,0)*r_v(1,1) + DDN_v[1](0,1)*r_v(1,0);
const double crRHS25 = DDN_v[2](0,0)*r_v(2,1) + DDN_v[2](0,1)*r_v(2,0);
const double crRHS26 = DDN_v[3](0,0)*r_v(3,1) + DDN_v[3](0,1)*r_v(3,0);
const double crRHS27 = DDN_v[4](0,0)*r_v(4,1) + DDN_v[4](0,1)*r_v(4,0);
const double crRHS28 = DDN_v[5](0,0)*r_v(5,1) + DDN_v[5](0,1)*r_v(5,0);
const double crRHS29 = 1.0*C(2,2);
const double crRHS30 = 1.0/(dyn_tau*rho/rData.DeltaTime + rho*stab_c2*sqrt(pow(crRHS4, 2) + pow(crRHS5, 2))/h + mu*stab_c1/pow(h, 2));
const double crRHS31 = crRHS30*(-DN_p(0,0)*r_p[0] - DN_p(1,0)*r_p[1] - DN_p(2,0)*r_p[2] + crRHS1 + crRHS10*crRHS21 + crRHS10*crRHS7 + crRHS11*crRHS21 + crRHS11*crRHS7 + crRHS12*crRHS21 + crRHS12*crRHS7 + crRHS13*crRHS21 + crRHS13*crRHS7 + crRHS14*crRHS15 + crRHS14*crRHS16 + crRHS14*crRHS17 + crRHS14*crRHS18 + crRHS14*crRHS19 + crRHS14*crRHS20 + crRHS15*crRHS22 + crRHS16*crRHS22 + crRHS17*crRHS22 + crRHS18*crRHS22 + crRHS19*crRHS22 - crRHS2 + crRHS20*crRHS22 + crRHS21*crRHS23 + crRHS21*crRHS24 + crRHS21*crRHS25 + crRHS21*crRHS26 + crRHS21*crRHS27 + crRHS21*crRHS28 + crRHS21*crRHS8 + crRHS21*crRHS9 + crRHS23*crRHS29 + crRHS24*crRHS29 + crRHS25*crRHS29 + crRHS26*crRHS29 + crRHS27*crRHS29 + crRHS28*crRHS29 - crRHS6 + crRHS7*crRHS8 + crRHS7*crRHS9);
const double crRHS32 = rho*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1));
const double crRHS33 = N_v[0]*crRHS32;
const double crRHS34 = rho*(DN_v(0,0)*crRHS4 + DN_v(0,1)*crRHS5);
const double crRHS35 = rho*(N_v[0]*r_f(0,1) + N_v[1]*r_f(1,1) + N_v[2]*r_f(2,1) + N_v[3]*r_f(3,1) + N_v[4]*r_f(4,1) + N_v[5]*r_f(5,1));
const double crRHS36 = rho*(N_v[0]*(rData.BDF0*r_v(0,1) + rData.BDF1*r_vn(0,1) + rData.BDF2*r_vnn(0,1)) + N_v[1]*(rData.BDF0*r_v(1,1) + rData.BDF1*r_vn(1,1) + rData.BDF2*r_vnn(1,1)) + N_v[2]*(rData.BDF0*r_v(2,1) + rData.BDF1*r_vn(2,1) + rData.BDF2*r_vnn(2,1)) + N_v[3]*(rData.BDF0*r_v(3,1) + rData.BDF1*r_vn(3,1) + rData.BDF2*r_vnn(3,1)) + N_v[4]*(rData.BDF0*r_v(4,1) + rData.BDF1*r_vn(4,1) + rData.BDF2*r_vnn(4,1)) + N_v[5]*(rData.BDF0*r_v(5,1) + rData.BDF1*r_vn(5,1) + rData.BDF2*r_vnn(5,1)));
const double crRHS37 = DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1);
const double crRHS38 = rho*(crRHS37*crRHS5 + crRHS4*(DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1)));
const double crRHS39 = DDN_v[0](1,0)*r_v(0,0);
const double crRHS40 = DDN_v[1](1,0)*r_v(1,0);
const double crRHS41 = DDN_v[2](1,0)*r_v(2,0);
const double crRHS42 = DDN_v[3](1,0)*r_v(3,0);
const double crRHS43 = DDN_v[4](1,0)*r_v(4,0);
const double crRHS44 = DDN_v[5](1,0)*r_v(5,0);
const double crRHS45 = 1.0*C(1,1);
const double crRHS46 = DDN_v[0](1,1)*r_v(0,1);
const double crRHS47 = DDN_v[1](1,1)*r_v(1,1);
const double crRHS48 = DDN_v[2](1,1)*r_v(2,1);
const double crRHS49 = DDN_v[3](1,1)*r_v(3,1);
const double crRHS50 = DDN_v[4](1,1)*r_v(4,1);
const double crRHS51 = DDN_v[5](1,1)*r_v(5,1);
const double crRHS52 = DDN_v[0](1,0)*r_v(0,1) + DDN_v[0](1,1)*r_v(0,0);
const double crRHS53 = DDN_v[1](1,0)*r_v(1,1) + DDN_v[1](1,1)*r_v(1,0);
const double crRHS54 = DDN_v[2](1,0)*r_v(2,1) + DDN_v[2](1,1)*r_v(2,0);
const double crRHS55 = DDN_v[3](1,0)*r_v(3,1) + DDN_v[3](1,1)*r_v(3,0);
const double crRHS56 = DDN_v[4](1,0)*r_v(4,1) + DDN_v[4](1,1)*r_v(4,0);
const double crRHS57 = DDN_v[5](1,0)*r_v(5,1) + DDN_v[5](1,1)*r_v(5,0);
const double crRHS58 = crRHS30*(-DN_p(0,1)*r_p[0] - DN_p(1,1)*r_p[1] - DN_p(2,1)*r_p[2] + crRHS14*crRHS39 + crRHS14*crRHS40 + crRHS14*crRHS41 + crRHS14*crRHS42 + crRHS14*crRHS43 + crRHS14*crRHS44 + crRHS21*crRHS39 + crRHS21*crRHS40 + crRHS21*crRHS41 + crRHS21*crRHS42 + crRHS21*crRHS43 + crRHS21*crRHS44 + crRHS22*crRHS46 + crRHS22*crRHS47 + crRHS22*crRHS48 + crRHS22*crRHS49 + crRHS22*crRHS50 + crRHS22*crRHS51 + crRHS22*crRHS52 + crRHS22*crRHS53 + crRHS22*crRHS54 + crRHS22*crRHS55 + crRHS22*crRHS56 + crRHS22*crRHS57 + crRHS29*crRHS52 + crRHS29*crRHS53 + crRHS29*crRHS54 + crRHS29*crRHS55 + crRHS29*crRHS56 + crRHS29*crRHS57 + crRHS35 - crRHS36 - crRHS38 + crRHS45*crRHS46 + crRHS45*crRHS47 + crRHS45*crRHS48 + crRHS45*crRHS49 + crRHS45*crRHS50 + crRHS45*crRHS51);
const double crRHS59 = N_v[1]*crRHS32;
const double crRHS60 = rho*(DN_v(1,0)*crRHS4 + DN_v(1,1)*crRHS5);
const double crRHS61 = N_v[2]*crRHS32;
const double crRHS62 = rho*(DN_v(2,0)*crRHS4 + DN_v(2,1)*crRHS5);
const double crRHS63 = N_v[3]*crRHS32;
const double crRHS64 = rho*(DN_v(3,0)*crRHS4 + DN_v(3,1)*crRHS5);
const double crRHS65 = N_v[4]*crRHS32;
const double crRHS66 = rho*(DN_v(4,0)*crRHS4 + DN_v(4,1)*crRHS5);
const double crRHS67 = N_v[5]*crRHS32;
const double crRHS68 = rho*(DN_v(5,0)*crRHS4 + DN_v(5,1)*crRHS5);
const double crRHS69 = crRHS3 + crRHS37;
rRHS[0]+=-gauss_weight*(-DN_v(0,0)*crRHS0 + DN_v(0,0)*r_stress[0] + DN_v(0,1)*r_stress[2] - N_v[0]*crRHS1 + N_v[0]*crRHS2 + N_v[0]*crRHS6 - crRHS31*crRHS33 - crRHS31*crRHS34);
rRHS[1]+=-gauss_weight*(DN_v(0,0)*r_stress[2] - DN_v(0,1)*crRHS0 + DN_v(0,1)*r_stress[1] - N_v[0]*crRHS35 + N_v[0]*crRHS36 + N_v[0]*crRHS38 - crRHS33*crRHS58 - crRHS34*crRHS58);
rRHS[2]+=-gauss_weight*(-DN_v(1,0)*crRHS0 + DN_v(1,0)*r_stress[0] + DN_v(1,1)*r_stress[2] - N_v[1]*crRHS1 + N_v[1]*crRHS2 + N_v[1]*crRHS6 - crRHS31*crRHS59 - crRHS31*crRHS60);
rRHS[3]+=-gauss_weight*(DN_v(1,0)*r_stress[2] - DN_v(1,1)*crRHS0 + DN_v(1,1)*r_stress[1] - N_v[1]*crRHS35 + N_v[1]*crRHS36 + N_v[1]*crRHS38 - crRHS58*crRHS59 - crRHS58*crRHS60);
rRHS[4]+=-gauss_weight*(-DN_v(2,0)*crRHS0 + DN_v(2,0)*r_stress[0] + DN_v(2,1)*r_stress[2] - N_v[2]*crRHS1 + N_v[2]*crRHS2 + N_v[2]*crRHS6 - crRHS31*crRHS61 - crRHS31*crRHS62);
rRHS[5]+=-gauss_weight*(DN_v(2,0)*r_stress[2] - DN_v(2,1)*crRHS0 + DN_v(2,1)*r_stress[1] - N_v[2]*crRHS35 + N_v[2]*crRHS36 + N_v[2]*crRHS38 - crRHS58*crRHS61 - crRHS58*crRHS62);
rRHS[6]+=-gauss_weight*(-DN_v(3,0)*crRHS0 + DN_v(3,0)*r_stress[0] + DN_v(3,1)*r_stress[2] - N_v[3]*crRHS1 + N_v[3]*crRHS2 + N_v[3]*crRHS6 - crRHS31*crRHS63 - crRHS31*crRHS64);
rRHS[7]+=-gauss_weight*(DN_v(3,0)*r_stress[2] - DN_v(3,1)*crRHS0 + DN_v(3,1)*r_stress[1] - N_v[3]*crRHS35 + N_v[3]*crRHS36 + N_v[3]*crRHS38 - crRHS58*crRHS63 - crRHS58*crRHS64);
rRHS[8]+=-gauss_weight*(-DN_v(4,0)*crRHS0 + DN_v(4,0)*r_stress[0] + DN_v(4,1)*r_stress[2] - N_v[4]*crRHS1 + N_v[4]*crRHS2 + N_v[4]*crRHS6 - crRHS31*crRHS65 - crRHS31*crRHS66);
rRHS[9]+=-gauss_weight*(DN_v(4,0)*r_stress[2] - DN_v(4,1)*crRHS0 + DN_v(4,1)*r_stress[1] - N_v[4]*crRHS35 + N_v[4]*crRHS36 + N_v[4]*crRHS38 - crRHS58*crRHS65 - crRHS58*crRHS66);
rRHS[10]+=-gauss_weight*(-DN_v(5,0)*crRHS0 + DN_v(5,0)*r_stress[0] + DN_v(5,1)*r_stress[2] - N_v[5]*crRHS1 + N_v[5]*crRHS2 + N_v[5]*crRHS6 - crRHS31*crRHS67 - crRHS31*crRHS68);
rRHS[11]+=-gauss_weight*(DN_v(5,0)*r_stress[2] - DN_v(5,1)*crRHS0 + DN_v(5,1)*r_stress[1] - N_v[5]*crRHS35 + N_v[5]*crRHS36 + N_v[5]*crRHS38 - crRHS58*crRHS67 - crRHS58*crRHS68);
rRHS[12]+=gauss_weight*(DN_p(0,0)*crRHS31 + DN_p(0,1)*crRHS58 - N_p[0]*crRHS69);
rRHS[13]+=gauss_weight*(DN_p(1,0)*crRHS31 + DN_p(1,1)*crRHS58 - N_p[1]*crRHS69);
rRHS[14]+=gauss_weight*(DN_p(2,0)*crRHS31 + DN_p(2,1)*crRHS58 - N_p[2]*crRHS69);

}

template <>
void IncompressibleNavierStokesDivStable<3>::ComputeGaussPointRHSContribution(
    const ElementDataContainer& rData,
    VectorType& rRHS)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double bdf1 = rData.BDF1;
    const double bdf2 = rData.BDF2;
    const double dt = rData.DeltaTime;

    // Get stabilization data
    const double h = rData.ElementSize;
    const double stab_c1 = rData.StabC1;
    const double stab_c2 = rData.StabC2;
    const double dyn_tau = rData.DynamicTau;

    // Get nodal data
    const auto& r_v = rData.Velocity;
    const auto& r_vn = rData.VelocityOld1;
    const auto& r_vnn = rData.VelocityOld2;
    const auto& r_vmesh = rData.MeshVelocity;
    const auto& r_f = rData.BodyForce;
    const auto& r_p = rData.Pressure;

    // Calculate convective velocity
    const BoundedMatrix<double,3,10> vconv = r_v - r_vmesh;

    // Get stress from material response
    const auto& r_stress = rData.ShearStress;
    const auto& C = rData.ConstitutiveMatrix;

    // Get shape function values
    const auto& N_p = rData.N_p;
    const auto& N_v = rData.N_v;
    const auto& DN_p = rData.DN_p;
    const auto& DN_v = rData.DN_v;
    const auto& DDN_v = rData.DDN_v;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;
    //substitute_rhs_3D
}

template <>
void IncompressibleNavierStokesDivStable<2>::ComputeGaussPointEnrichmentContribution(
    const ElementDataContainer& rData,
    array_1d<double, 2>& rRHSee,
    BoundedMatrix<double, LocalSize, 2>& rKue,
    BoundedMatrix<double, 2, LocalSize>& rKeu,
    BoundedMatrix<double, 2, 2>& rKee)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double bdf1 = rData.BDF1;
    const double bdf2 = rData.BDF2;
    const double dt = rData.DeltaTime;

    // Get stabilization data
    const double h = rData.ElementSize;
    const double stab_c1 = rData.StabC1;
    const double stab_c2 = rData.StabC2;
    const double dyn_tau = rData.DynamicTau;

    // Get nodal data
    const auto& r_v = rData.Velocity;
    const auto& r_vn = rData.VelocityOld1;
    const auto& r_vnn = rData.VelocityOld2;
    const auto& r_vmesh = rData.MeshVelocity;
    const auto& r_f = rData.BodyForce;
    const auto& r_p = rData.Pressure;

    // Calculate convective velocity
    const BoundedMatrix<double,3,10> vconv = r_v - r_vmesh;

    // Get stress from material response
    const auto& r_stress = rData.ShearStress;
    const auto& C = rData.ConstitutiveMatrix;

    // Get shape function values
    const auto& N_p = rData.N_p;
    const auto& N_v = rData.N_v;
    const auto& DN_p = rData.DN_p;
    const auto& DN_v = rData.DN_v;
    const double N_e = rData.N_e;
    const auto& DN_e = rData.DN_e;
    const auto& DDN_v = rData.DDN_v;

    // Assemble enrichment contribution
    const double gauss_weight = rData.Weight;
    const BoundedMatrix<double,1,2> v_e = ZeroMatrix(1,2);

    const double crRHSee0 = N_p[0]*r_p[0] + N_p[1]*r_p[1] + N_p[2]*r_p[2];
const double crRHSee1 = DN_e(0,0)*v_e(0,0);
const double crRHSee2 = DN_e(0,1)*v_e(0,1);
const double crRHSee3 = DN_e(0,0)*v_e(0,1);
const double crRHSee4 = DN_e(0,1)*v_e(0,0);
const double crRHSee5 = crRHSee3 + crRHSee4;
const double crRHSee6 = C(0,2)*crRHSee1 + C(1,2)*crRHSee2 + C(2,2)*crRHSee5 + r_stress[2];
const double crRHSee7 = rho*(N_v[0]*r_f(0,0) + N_v[1]*r_f(1,0) + N_v[2]*r_f(2,0) + N_v[3]*r_f(3,0) + N_v[4]*r_f(4,0) + N_v[5]*r_f(5,0));
const double crRHSee8 = rho*(N_v[0]*(rData.BDF0*r_v(0,0) + rData.BDF1*r_vn(0,0) + rData.BDF2*r_vnn(0,0)) + N_v[1]*(rData.BDF0*r_v(1,0) + rData.BDF1*r_vn(1,0) + rData.BDF2*r_vnn(1,0)) + N_v[2]*(rData.BDF0*r_v(2,0) + rData.BDF1*r_vn(2,0) + rData.BDF2*r_vnn(2,0)) + N_v[3]*(rData.BDF0*r_v(3,0) + rData.BDF1*r_vn(3,0) + rData.BDF2*r_vnn(3,0)) + N_v[4]*(rData.BDF0*r_v(4,0) + rData.BDF1*r_vn(4,0) + rData.BDF2*r_vnn(4,0)) + N_v[5]*(rData.BDF0*r_v(5,0) + rData.BDF1*r_vn(5,0) + rData.BDF2*r_vnn(5,0)));
const double crRHSee9 = N_e*v_e(0,1);
const double crRHSee10 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crRHSee11 = DN_e(0,0)*crRHSee10;
const double crRHSee12 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crRHSee13 = DN_v(0,0)*r_v(0,0) + DN_v(1,0)*r_v(1,0) + DN_v(2,0)*r_v(2,0) + DN_v(3,0)*r_v(3,0) + DN_v(4,0)*r_v(4,0) + DN_v(5,0)*r_v(5,0);
const double crRHSee14 = N_e*v_e(0,0);
const double crRHSee15 = DN_v(0,1)*r_v(0,0) + DN_v(1,1)*r_v(1,0) + DN_v(2,1)*r_v(2,0) + DN_v(3,1)*r_v(3,0) + DN_v(4,1)*r_v(4,0) + DN_v(5,1)*r_v(5,0);
const double crRHSee16 = rho*(DN_e(0,0)*N_e*pow(v_e(0,0), 2) + crRHSee10*crRHSee13 + crRHSee11*v_e(0,0) + crRHSee12*crRHSee15 + crRHSee12*crRHSee4 + crRHSee13*crRHSee14 + crRHSee15*crRHSee9 + crRHSee4*crRHSee9);
const double crRHSee17 = 1.0*C(0,0);
const double crRHSee18 = DDN_v[0](0,0)*r_v(0,0);
const double crRHSee19 = DDN_v[1](0,0)*r_v(1,0);
const double crRHSee20 = DDN_v[2](0,0)*r_v(2,0);
const double crRHSee21 = DDN_v[3](0,0)*r_v(3,0);
const double crRHSee22 = DDN_v[4](0,0)*r_v(4,0);
const double crRHSee23 = DDN_v[5](0,0)*r_v(5,0);
const double crRHSee24 = 1.0*C(0,1);
const double crRHSee25 = DDN_v[0](0,1)*r_v(0,1);
const double crRHSee26 = DDN_v[1](0,1)*r_v(1,1);
const double crRHSee27 = DDN_v[2](0,1)*r_v(2,1);
const double crRHSee28 = DDN_v[3](0,1)*r_v(3,1);
const double crRHSee29 = DDN_v[4](0,1)*r_v(4,1);
const double crRHSee30 = DDN_v[5](0,1)*r_v(5,1);
const double crRHSee31 = 1.0*C(0,2);
const double crRHSee32 = 1.0*C(1,2);
const double crRHSee33 = DDN_v[0](0,0)*r_v(0,1) + DDN_v[0](0,1)*r_v(0,0);
const double crRHSee34 = DDN_v[1](0,0)*r_v(1,1) + DDN_v[1](0,1)*r_v(1,0);
const double crRHSee35 = DDN_v[2](0,0)*r_v(2,1) + DDN_v[2](0,1)*r_v(2,0);
const double crRHSee36 = DDN_v[3](0,0)*r_v(3,1) + DDN_v[3](0,1)*r_v(3,0);
const double crRHSee37 = DDN_v[4](0,0)*r_v(4,1) + DDN_v[4](0,1)*r_v(4,0);
const double crRHSee38 = DDN_v[5](0,0)*r_v(5,1) + DDN_v[5](0,1)*r_v(5,0);
const double crRHSee39 = 1.0*C(2,2);
const double crRHSee40 = rho/(dyn_tau*rho/rData.DeltaTime + rho*stab_c2*sqrt(pow(crRHSee10, 2) + pow(crRHSee12, 2))/h + mu*stab_c1/pow(h, 2));
const double crRHSee41 = crRHSee40*(-DN_p(0,0)*r_p[0] - DN_p(1,0)*r_p[1] - DN_p(2,0)*r_p[2] - crRHSee16 + crRHSee17*crRHSee18 + crRHSee17*crRHSee19 + crRHSee17*crRHSee20 + crRHSee17*crRHSee21 + crRHSee17*crRHSee22 + crRHSee17*crRHSee23 + crRHSee18*crRHSee31 + crRHSee19*crRHSee31 + crRHSee20*crRHSee31 + crRHSee21*crRHSee31 + crRHSee22*crRHSee31 + crRHSee23*crRHSee31 + crRHSee24*crRHSee25 + crRHSee24*crRHSee26 + crRHSee24*crRHSee27 + crRHSee24*crRHSee28 + crRHSee24*crRHSee29 + crRHSee24*crRHSee30 + crRHSee25*crRHSee32 + crRHSee26*crRHSee32 + crRHSee27*crRHSee32 + crRHSee28*crRHSee32 + crRHSee29*crRHSee32 + crRHSee30*crRHSee32 + crRHSee31*crRHSee33 + crRHSee31*crRHSee34 + crRHSee31*crRHSee35 + crRHSee31*crRHSee36 + crRHSee31*crRHSee37 + crRHSee31*crRHSee38 + crRHSee33*crRHSee39 + crRHSee34*crRHSee39 + crRHSee35*crRHSee39 + crRHSee36*crRHSee39 + crRHSee37*crRHSee39 + crRHSee38*crRHSee39 + crRHSee7 - crRHSee8);
const double crRHSee42 = N_e*crRHSee41;
const double crRHSee43 = 2.0*crRHSee1 + 2.0*crRHSee2;
const double crRHSee44 = 1.0*DN_v(0,0)*vconv(0,0) + 1.0*DN_v(0,1)*vconv(0,1) + 1.0*DN_v(1,0)*vconv(1,0) + 1.0*DN_v(1,1)*vconv(1,1) + 1.0*DN_v(2,0)*vconv(2,0) + 1.0*DN_v(2,1)*vconv(2,1) + 1.0*DN_v(3,0)*vconv(3,0) + 1.0*DN_v(3,1)*vconv(3,1) + 1.0*DN_v(4,0)*vconv(4,0) + 1.0*DN_v(4,1)*vconv(4,1) + 1.0*DN_v(5,0)*vconv(5,0) + 1.0*DN_v(5,1)*vconv(5,1);
const double crRHSee45 = DN_e(0,1)*crRHSee12;
const double crRHSee46 = 1.0*crRHSee11 + 1.0*crRHSee45;
const double crRHSee47 = rho*(N_v[0]*r_f(0,1) + N_v[1]*r_f(1,1) + N_v[2]*r_f(2,1) + N_v[3]*r_f(3,1) + N_v[4]*r_f(4,1) + N_v[5]*r_f(5,1));
const double crRHSee48 = rho*(N_v[0]*(rData.BDF0*r_v(0,1) + rData.BDF1*r_vn(0,1) + rData.BDF2*r_vnn(0,1)) + N_v[1]*(rData.BDF0*r_v(1,1) + rData.BDF1*r_vn(1,1) + rData.BDF2*r_vnn(1,1)) + N_v[2]*(rData.BDF0*r_v(2,1) + rData.BDF1*r_vn(2,1) + rData.BDF2*r_vnn(2,1)) + N_v[3]*(rData.BDF0*r_v(3,1) + rData.BDF1*r_vn(3,1) + rData.BDF2*r_vnn(3,1)) + N_v[4]*(rData.BDF0*r_v(4,1) + rData.BDF1*r_vn(4,1) + rData.BDF2*r_vnn(4,1)) + N_v[5]*(rData.BDF0*r_v(5,1) + rData.BDF1*r_vn(5,1) + rData.BDF2*r_vnn(5,1)));
const double crRHSee49 = DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1);
const double crRHSee50 = DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1);
const double crRHSee51 = rho*(DN_e(0,1)*N_e*pow(v_e(0,1), 2) + crRHSee10*crRHSee3 + crRHSee10*crRHSee49 + crRHSee12*crRHSee50 + crRHSee14*crRHSee3 + crRHSee14*crRHSee49 + crRHSee45*v_e(0,1) + crRHSee50*crRHSee9);
const double crRHSee52 = DDN_v[0](1,0)*r_v(0,0);
const double crRHSee53 = DDN_v[1](1,0)*r_v(1,0);
const double crRHSee54 = DDN_v[2](1,0)*r_v(2,0);
const double crRHSee55 = DDN_v[3](1,0)*r_v(3,0);
const double crRHSee56 = DDN_v[4](1,0)*r_v(4,0);
const double crRHSee57 = DDN_v[5](1,0)*r_v(5,0);
const double crRHSee58 = 1.0*C(1,1);
const double crRHSee59 = DDN_v[0](1,1)*r_v(0,1);
const double crRHSee60 = DDN_v[1](1,1)*r_v(1,1);
const double crRHSee61 = DDN_v[2](1,1)*r_v(2,1);
const double crRHSee62 = DDN_v[3](1,1)*r_v(3,1);
const double crRHSee63 = DDN_v[4](1,1)*r_v(4,1);
const double crRHSee64 = DDN_v[5](1,1)*r_v(5,1);
const double crRHSee65 = DDN_v[0](1,0)*r_v(0,1) + DDN_v[0](1,1)*r_v(0,0);
const double crRHSee66 = DDN_v[1](1,0)*r_v(1,1) + DDN_v[1](1,1)*r_v(1,0);
const double crRHSee67 = DDN_v[2](1,0)*r_v(2,1) + DDN_v[2](1,1)*r_v(2,0);
const double crRHSee68 = DDN_v[3](1,0)*r_v(3,1) + DDN_v[3](1,1)*r_v(3,0);
const double crRHSee69 = DDN_v[4](1,0)*r_v(4,1) + DDN_v[4](1,1)*r_v(4,0);
const double crRHSee70 = DDN_v[5](1,0)*r_v(5,1) + DDN_v[5](1,1)*r_v(5,0);
const double crRHSee71 = crRHSee40*(-DN_p(0,1)*r_p[0] - DN_p(1,1)*r_p[1] - DN_p(2,1)*r_p[2] + crRHSee24*crRHSee52 + crRHSee24*crRHSee53 + crRHSee24*crRHSee54 + crRHSee24*crRHSee55 + crRHSee24*crRHSee56 + crRHSee24*crRHSee57 + crRHSee31*crRHSee52 + crRHSee31*crRHSee53 + crRHSee31*crRHSee54 + crRHSee31*crRHSee55 + crRHSee31*crRHSee56 + crRHSee31*crRHSee57 + crRHSee32*crRHSee59 + crRHSee32*crRHSee60 + crRHSee32*crRHSee61 + crRHSee32*crRHSee62 + crRHSee32*crRHSee63 + crRHSee32*crRHSee64 + crRHSee32*crRHSee65 + crRHSee32*crRHSee66 + crRHSee32*crRHSee67 + crRHSee32*crRHSee68 + crRHSee32*crRHSee69 + crRHSee32*crRHSee70 + crRHSee39*crRHSee65 + crRHSee39*crRHSee66 + crRHSee39*crRHSee67 + crRHSee39*crRHSee68 + crRHSee39*crRHSee69 + crRHSee39*crRHSee70 + crRHSee47 - crRHSee48 - crRHSee51 + crRHSee58*crRHSee59 + crRHSee58*crRHSee60 + crRHSee58*crRHSee61 + crRHSee58*crRHSee62 + crRHSee58*crRHSee63 + crRHSee58*crRHSee64);
const double crRHSee72 = N_e*crRHSee71;
rRHSee[0]+=gauss_weight*(DN_e(0,0)*crRHSee0 - DN_e(0,0)*(C(0,0)*crRHSee1 + C(0,1)*crRHSee2 + C(0,2)*crRHSee5 + r_stress[0]) - DN_e(0,1)*crRHSee6 - N_e*crRHSee16 + N_e*crRHSee7 - N_e*crRHSee8 + crRHSee41*crRHSee46 + crRHSee42*crRHSee43 + crRHSee42*crRHSee44);
rRHSee[1]+=gauss_weight*(-DN_e(0,0)*crRHSee6 + DN_e(0,1)*crRHSee0 - DN_e(0,1)*(C(0,1)*crRHSee1 + C(1,1)*crRHSee2 + C(1,2)*crRHSee5 + r_stress[1]) + N_e*crRHSee47 - N_e*crRHSee48 - N_e*crRHSee51 + crRHSee43*crRHSee72 + crRHSee44*crRHSee72 + crRHSee46*crRHSee71);


    const double crKue0 = C(0,2)*DN_v(0,0);
const double crKue1 = C(2,2)*DN_v(0,1) + crKue0;
const double crKue2 = DN_e(0,1)*N_e*v_e(0,1);
const double crKue3 = DN_e(0,0)*N_e*v_e(0,0);
const double crKue4 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crKue5 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crKue6 = DN_e(0,0)*crKue4 + DN_e(0,1)*crKue5;
const double crKue7 = N_e*(DN_v(0,0)*r_v(0,0) + DN_v(1,0)*r_v(1,0) + DN_v(2,0)*r_v(2,0) + DN_v(3,0)*r_v(3,0) + DN_v(4,0)*r_v(4,0) + DN_v(5,0)*r_v(5,0)) + crKue2 + 2*crKue3 + crKue6;
const double crKue8 = N_v[0]*rho;
const double crKue9 = 1.0/(dyn_tau*rho/rData.DeltaTime + rho*stab_c2*sqrt(pow(crKue4, 2) + pow(crKue5, 2))/h + mu*stab_c1/pow(h, 2));
const double crKue10 = crKue9*pow(rho, 2);
const double crKue11 = crKue10*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1));
const double crKue12 = N_v[0]*crKue11;
const double crKue13 = crKue10*(DN_v(0,0)*crKue4 + DN_v(0,1)*crKue5);
const double crKue14 = C(1,2)*DN_v(0,1);
const double crKue15 = DN_e(0,1)*v_e(0,0) + DN_v(0,1)*r_v(0,0) + DN_v(1,1)*r_v(1,0) + DN_v(2,1)*r_v(2,0) + DN_v(3,1)*r_v(3,0) + DN_v(4,1)*r_v(4,0) + DN_v(5,1)*r_v(5,0);
const double crKue16 = N_e*crKue8;
const double crKue17 = N_e*crKue12;
const double crKue18 = N_e*crKue13;
const double crKue19 = C(2,2)*DN_v(0,0) + crKue14;
const double crKue20 = DN_e(0,0)*v_e(0,1) + DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1);
const double crKue21 = N_e*(DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1)) + 2*crKue2 + crKue3 + crKue6;
const double crKue22 = C(0,2)*DN_v(1,0);
const double crKue23 = C(2,2)*DN_v(1,1) + crKue22;
const double crKue24 = N_v[1]*rho;
const double crKue25 = N_v[1]*crKue11;
const double crKue26 = crKue10*(DN_v(1,0)*crKue4 + DN_v(1,1)*crKue5);
const double crKue27 = C(1,2)*DN_v(1,1);
const double crKue28 = N_e*crKue24;
const double crKue29 = N_e*crKue25;
const double crKue30 = N_e*crKue26;
const double crKue31 = C(2,2)*DN_v(1,0) + crKue27;
const double crKue32 = C(0,2)*DN_v(2,0);
const double crKue33 = C(2,2)*DN_v(2,1) + crKue32;
const double crKue34 = N_v[2]*rho;
const double crKue35 = N_v[2]*crKue11;
const double crKue36 = crKue10*(DN_v(2,0)*crKue4 + DN_v(2,1)*crKue5);
const double crKue37 = C(1,2)*DN_v(2,1);
const double crKue38 = N_e*crKue34;
const double crKue39 = N_e*crKue35;
const double crKue40 = N_e*crKue36;
const double crKue41 = C(2,2)*DN_v(2,0) + crKue37;
const double crKue42 = C(0,2)*DN_v(3,0);
const double crKue43 = C(2,2)*DN_v(3,1) + crKue42;
const double crKue44 = N_v[3]*rho;
const double crKue45 = N_v[3]*crKue11;
const double crKue46 = crKue10*(DN_v(3,0)*crKue4 + DN_v(3,1)*crKue5);
const double crKue47 = C(1,2)*DN_v(3,1);
const double crKue48 = N_e*crKue44;
const double crKue49 = N_e*crKue45;
const double crKue50 = N_e*crKue46;
const double crKue51 = C(2,2)*DN_v(3,0) + crKue47;
const double crKue52 = C(0,2)*DN_v(4,0);
const double crKue53 = C(2,2)*DN_v(4,1) + crKue52;
const double crKue54 = N_v[4]*rho;
const double crKue55 = N_v[4]*crKue11;
const double crKue56 = crKue10*(DN_v(4,0)*crKue4 + DN_v(4,1)*crKue5);
const double crKue57 = C(1,2)*DN_v(4,1);
const double crKue58 = N_e*crKue54;
const double crKue59 = N_e*crKue55;
const double crKue60 = N_e*crKue56;
const double crKue61 = C(2,2)*DN_v(4,0) + crKue57;
const double crKue62 = C(0,2)*DN_v(5,0);
const double crKue63 = C(2,2)*DN_v(5,1) + crKue62;
const double crKue64 = N_v[5]*rho;
const double crKue65 = N_v[5]*crKue11;
const double crKue66 = crKue10*(DN_v(5,0)*crKue4 + DN_v(5,1)*crKue5);
const double crKue67 = C(1,2)*DN_v(5,1);
const double crKue68 = N_e*crKue64;
const double crKue69 = N_e*crKue65;
const double crKue70 = N_e*crKue66;
const double crKue71 = C(2,2)*DN_v(5,0) + crKue67;
const double crKue72 = crKue9*rho;
const double crKue73 = N_e*crKue72;
const double crKue74 = crKue20*crKue73;
const double crKue75 = crKue7*crKue72;
const double crKue76 = crKue15*crKue73;
const double crKue77 = crKue21*crKue72;
rKue(0,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(0,0) + C(0,2)*DN_v(0,1)) + DN_e(0,1)*crKue1 + crKue12*crKue7 + crKue13*crKue7 + crKue7*crKue8);
rKue(0,1)+=gauss_weight*(DN_e(0,0)*crKue1 + DN_e(0,1)*(C(0,1)*DN_v(0,0) + crKue14) + crKue15*crKue16 + crKue15*crKue17 + crKue15*crKue18);
rKue(1,0)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(0,1) + crKue0) + DN_e(0,1)*crKue19 + crKue16*crKue20 + crKue17*crKue20 + crKue18*crKue20);
rKue(1,1)+=gauss_weight*(DN_e(0,0)*crKue19 + DN_e(0,1)*(C(1,1)*DN_v(0,1) + C(1,2)*DN_v(0,0)) + crKue12*crKue21 + crKue13*crKue21 + crKue21*crKue8);
rKue(2,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(1,0) + C(0,2)*DN_v(1,1)) + DN_e(0,1)*crKue23 + crKue24*crKue7 + crKue25*crKue7 + crKue26*crKue7);
rKue(2,1)+=gauss_weight*(DN_e(0,0)*crKue23 + DN_e(0,1)*(C(0,1)*DN_v(1,0) + crKue27) + crKue15*crKue28 + crKue15*crKue29 + crKue15*crKue30);
rKue(3,0)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(1,1) + crKue22) + DN_e(0,1)*crKue31 + crKue20*crKue28 + crKue20*crKue29 + crKue20*crKue30);
rKue(3,1)+=gauss_weight*(DN_e(0,0)*crKue31 + DN_e(0,1)*(C(1,1)*DN_v(1,1) + C(1,2)*DN_v(1,0)) + crKue21*crKue24 + crKue21*crKue25 + crKue21*crKue26);
rKue(4,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(2,0) + C(0,2)*DN_v(2,1)) + DN_e(0,1)*crKue33 + crKue34*crKue7 + crKue35*crKue7 + crKue36*crKue7);
rKue(4,1)+=gauss_weight*(DN_e(0,0)*crKue33 + DN_e(0,1)*(C(0,1)*DN_v(2,0) + crKue37) + crKue15*crKue38 + crKue15*crKue39 + crKue15*crKue40);
rKue(5,0)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(2,1) + crKue32) + DN_e(0,1)*crKue41 + crKue20*crKue38 + crKue20*crKue39 + crKue20*crKue40);
rKue(5,1)+=gauss_weight*(DN_e(0,0)*crKue41 + DN_e(0,1)*(C(1,1)*DN_v(2,1) + C(1,2)*DN_v(2,0)) + crKue21*crKue34 + crKue21*crKue35 + crKue21*crKue36);
rKue(6,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(3,0) + C(0,2)*DN_v(3,1)) + DN_e(0,1)*crKue43 + crKue44*crKue7 + crKue45*crKue7 + crKue46*crKue7);
rKue(6,1)+=gauss_weight*(DN_e(0,0)*crKue43 + DN_e(0,1)*(C(0,1)*DN_v(3,0) + crKue47) + crKue15*crKue48 + crKue15*crKue49 + crKue15*crKue50);
rKue(7,0)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(3,1) + crKue42) + DN_e(0,1)*crKue51 + crKue20*crKue48 + crKue20*crKue49 + crKue20*crKue50);
rKue(7,1)+=gauss_weight*(DN_e(0,0)*crKue51 + DN_e(0,1)*(C(1,1)*DN_v(3,1) + C(1,2)*DN_v(3,0)) + crKue21*crKue44 + crKue21*crKue45 + crKue21*crKue46);
rKue(8,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(4,0) + C(0,2)*DN_v(4,1)) + DN_e(0,1)*crKue53 + crKue54*crKue7 + crKue55*crKue7 + crKue56*crKue7);
rKue(8,1)+=gauss_weight*(DN_e(0,0)*crKue53 + DN_e(0,1)*(C(0,1)*DN_v(4,0) + crKue57) + crKue15*crKue58 + crKue15*crKue59 + crKue15*crKue60);
rKue(9,0)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(4,1) + crKue52) + DN_e(0,1)*crKue61 + crKue20*crKue58 + crKue20*crKue59 + crKue20*crKue60);
rKue(9,1)+=gauss_weight*(DN_e(0,0)*crKue61 + DN_e(0,1)*(C(1,1)*DN_v(4,1) + C(1,2)*DN_v(4,0)) + crKue21*crKue54 + crKue21*crKue55 + crKue21*crKue56);
rKue(10,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(5,0) + C(0,2)*DN_v(5,1)) + DN_e(0,1)*crKue63 + crKue64*crKue7 + crKue65*crKue7 + crKue66*crKue7);
rKue(10,1)+=gauss_weight*(DN_e(0,0)*crKue63 + DN_e(0,1)*(C(0,1)*DN_v(5,0) + crKue67) + crKue15*crKue68 + crKue15*crKue69 + crKue15*crKue70);
rKue(11,0)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(5,1) + crKue62) + DN_e(0,1)*crKue71 + crKue20*crKue68 + crKue20*crKue69 + crKue20*crKue70);
rKue(11,1)+=gauss_weight*(DN_e(0,0)*crKue71 + DN_e(0,1)*(C(1,1)*DN_v(5,1) + C(1,2)*DN_v(5,0)) + crKue21*crKue64 + crKue21*crKue65 + crKue21*crKue66);
rKue(12,0)+=gauss_weight*(DN_e(0,0)*N_p[0] + DN_p(0,0)*crKue75 + DN_p(0,1)*crKue74);
rKue(12,1)+=gauss_weight*(DN_e(0,1)*N_p[0] + DN_p(0,0)*crKue76 + DN_p(0,1)*crKue77);
rKue(13,0)+=gauss_weight*(DN_e(0,0)*N_p[1] + DN_p(1,0)*crKue75 + DN_p(1,1)*crKue74);
rKue(13,1)+=gauss_weight*(DN_e(0,1)*N_p[1] + DN_p(1,0)*crKue76 + DN_p(1,1)*crKue77);
rKue(14,0)+=gauss_weight*(DN_e(0,0)*N_p[2] + DN_p(2,0)*crKue75 + DN_p(2,1)*crKue74);
rKue(14,1)+=gauss_weight*(DN_e(0,1)*N_p[2] + DN_p(2,0)*crKue76 + DN_p(2,1)*crKue77);


    const double crKeu0 = C(0,2)*DN_v(0,0);
const double crKeu1 = C(2,2)*DN_v(0,1) + crKeu0;
const double crKeu2 = 1.0*C(0,0);
const double crKeu3 = C(0,2)*DDN_v[0](0,0);
const double crKeu4 = 1.0*DDN_v[0](0,1);
const double crKeu5 = N_e*v_e(0,0);
const double crKeu6 = N_e*v_e(0,1);
const double crKeu7 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crKeu8 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crKeu9 = rho*(DN_v(0,0)*crKeu5 + DN_v(0,0)*crKeu7 + DN_v(0,1)*crKeu6 + DN_v(0,1)*crKeu8);
const double crKeu10 = rData.BDF0*rho;
const double crKeu11 = N_v[0]*crKeu10;
const double crKeu12 = -crKeu11 - crKeu9;
const double crKeu13 = rho/(dyn_tau*rho/rData.DeltaTime + rho*stab_c2*sqrt(pow(crKeu7, 2) + pow(crKeu8, 2))/h + mu*stab_c1/pow(h, 2));
const double crKeu14 = crKeu13*(C(0,2)*crKeu4 + C(2,2)*crKeu4 + DDN_v[0](0,0)*crKeu2 + crKeu12 + 1.0*crKeu3);
const double crKeu15 = 1.0*DN_e(0,0)*crKeu7 + 1.0*DN_e(0,1)*crKeu8;
const double crKeu16 = N_e*crKeu14;
const double crKeu17 = 1.0*DN_v(0,0)*vconv(0,0) + 1.0*DN_v(0,1)*vconv(0,1) + 1.0*DN_v(1,0)*vconv(1,0) + 1.0*DN_v(1,1)*vconv(1,1) + 1.0*DN_v(2,0)*vconv(2,0) + 1.0*DN_v(2,1)*vconv(2,1) + 1.0*DN_v(3,0)*vconv(3,0) + 1.0*DN_v(3,1)*vconv(3,1) + 1.0*DN_v(4,0)*vconv(4,0) + 1.0*DN_v(4,1)*vconv(4,1) + 1.0*DN_v(5,0)*vconv(5,0) + 1.0*DN_v(5,1)*vconv(5,1);
const double crKeu18 = 2.0*DN_e(0,0)*v_e(0,0) + 2.0*DN_e(0,1)*v_e(0,1);
const double crKeu19 = N_e*crKeu11 + N_e*crKeu9;
const double crKeu20 = C(1,2)*DN_v(0,1);
const double crKeu21 = C(2,2)*DN_v(0,0) + crKeu20;
const double crKeu22 = C(0,1)*DDN_v[0](0,1) + C(1,2)*DDN_v[0](0,1) + C(2,2)*DDN_v[0](0,0) + crKeu3;
const double crKeu23 = N_e*crKeu13;
const double crKeu24 = crKeu22*crKeu23;
const double crKeu25 = crKeu13*crKeu15;
const double crKeu26 = C(0,2)*DN_v(1,0);
const double crKeu27 = C(2,2)*DN_v(1,1) + crKeu26;
const double crKeu28 = C(0,2)*DDN_v[1](0,0);
const double crKeu29 = 1.0*DDN_v[1](0,1);
const double crKeu30 = rho*(DN_v(1,0)*crKeu5 + DN_v(1,0)*crKeu7 + DN_v(1,1)*crKeu6 + DN_v(1,1)*crKeu8);
const double crKeu31 = N_v[1]*crKeu10;
const double crKeu32 = -crKeu30 - crKeu31;
const double crKeu33 = C(0,2)*crKeu29 + C(2,2)*crKeu29 + DDN_v[1](0,0)*crKeu2 + 1.0*crKeu28 + crKeu32;
const double crKeu34 = crKeu23*crKeu33;
const double crKeu35 = N_e*crKeu30 + N_e*crKeu31;
const double crKeu36 = C(1,2)*DN_v(1,1);
const double crKeu37 = C(2,2)*DN_v(1,0) + crKeu36;
const double crKeu38 = C(0,1)*DDN_v[1](0,1) + C(1,2)*DDN_v[1](0,1) + C(2,2)*DDN_v[1](0,0) + crKeu28;
const double crKeu39 = crKeu23*crKeu38;
const double crKeu40 = C(0,2)*DN_v(2,0);
const double crKeu41 = C(2,2)*DN_v(2,1) + crKeu40;
const double crKeu42 = C(0,2)*DDN_v[2](0,0);
const double crKeu43 = 1.0*DDN_v[2](0,1);
const double crKeu44 = rho*(DN_v(2,0)*crKeu5 + DN_v(2,0)*crKeu7 + DN_v(2,1)*crKeu6 + DN_v(2,1)*crKeu8);
const double crKeu45 = N_v[2]*crKeu10;
const double crKeu46 = -crKeu44 - crKeu45;
const double crKeu47 = C(0,2)*crKeu43 + C(2,2)*crKeu43 + DDN_v[2](0,0)*crKeu2 + 1.0*crKeu42 + crKeu46;
const double crKeu48 = crKeu23*crKeu47;
const double crKeu49 = N_e*crKeu44 + N_e*crKeu45;
const double crKeu50 = C(1,2)*DN_v(2,1);
const double crKeu51 = C(2,2)*DN_v(2,0) + crKeu50;
const double crKeu52 = C(0,1)*DDN_v[2](0,1) + C(1,2)*DDN_v[2](0,1) + C(2,2)*DDN_v[2](0,0) + crKeu42;
const double crKeu53 = crKeu23*crKeu52;
const double crKeu54 = C(0,2)*DN_v(3,0);
const double crKeu55 = C(2,2)*DN_v(3,1) + crKeu54;
const double crKeu56 = C(0,2)*DDN_v[3](0,0);
const double crKeu57 = 1.0*DDN_v[3](0,1);
const double crKeu58 = rho*(DN_v(3,0)*crKeu5 + DN_v(3,0)*crKeu7 + DN_v(3,1)*crKeu6 + DN_v(3,1)*crKeu8);
const double crKeu59 = N_v[3]*crKeu10;
const double crKeu60 = -crKeu58 - crKeu59;
const double crKeu61 = C(0,2)*crKeu57 + C(2,2)*crKeu57 + DDN_v[3](0,0)*crKeu2 + 1.0*crKeu56 + crKeu60;
const double crKeu62 = crKeu23*crKeu61;
const double crKeu63 = N_e*crKeu58 + N_e*crKeu59;
const double crKeu64 = C(1,2)*DN_v(3,1);
const double crKeu65 = C(2,2)*DN_v(3,0) + crKeu64;
const double crKeu66 = C(0,1)*DDN_v[3](0,1) + C(1,2)*DDN_v[3](0,1) + C(2,2)*DDN_v[3](0,0) + crKeu56;
const double crKeu67 = crKeu23*crKeu66;
const double crKeu68 = C(0,2)*DN_v(4,0);
const double crKeu69 = C(2,2)*DN_v(4,1) + crKeu68;
const double crKeu70 = C(0,2)*DDN_v[4](0,0);
const double crKeu71 = 1.0*DDN_v[4](0,1);
const double crKeu72 = rho*(DN_v(4,0)*crKeu5 + DN_v(4,0)*crKeu7 + DN_v(4,1)*crKeu6 + DN_v(4,1)*crKeu8);
const double crKeu73 = N_v[4]*crKeu10;
const double crKeu74 = -crKeu72 - crKeu73;
const double crKeu75 = C(0,2)*crKeu71 + C(2,2)*crKeu71 + DDN_v[4](0,0)*crKeu2 + 1.0*crKeu70 + crKeu74;
const double crKeu76 = crKeu23*crKeu75;
const double crKeu77 = N_e*crKeu72 + N_e*crKeu73;
const double crKeu78 = C(1,2)*DN_v(4,1);
const double crKeu79 = C(2,2)*DN_v(4,0) + crKeu78;
const double crKeu80 = C(0,1)*DDN_v[4](0,1) + C(1,2)*DDN_v[4](0,1) + C(2,2)*DDN_v[4](0,0) + crKeu70;
const double crKeu81 = crKeu23*crKeu80;
const double crKeu82 = C(0,2)*DN_v(5,0);
const double crKeu83 = C(2,2)*DN_v(5,1) + crKeu82;
const double crKeu84 = C(0,2)*DDN_v[5](0,0);
const double crKeu85 = 1.0*DDN_v[5](0,1);
const double crKeu86 = rho*(DN_v(5,0)*crKeu5 + DN_v(5,0)*crKeu7 + DN_v(5,1)*crKeu6 + DN_v(5,1)*crKeu8);
const double crKeu87 = N_v[5]*crKeu10;
const double crKeu88 = -crKeu86 - crKeu87;
const double crKeu89 = C(0,2)*crKeu85 + C(2,2)*crKeu85 + DDN_v[5](0,0)*crKeu2 + 1.0*crKeu84 + crKeu88;
const double crKeu90 = crKeu23*crKeu89;
const double crKeu91 = N_e*crKeu86 + N_e*crKeu87;
const double crKeu92 = C(1,2)*DN_v(5,1);
const double crKeu93 = C(2,2)*DN_v(5,0) + crKeu92;
const double crKeu94 = C(0,1)*DDN_v[5](0,1) + C(1,2)*DDN_v[5](0,1) + C(2,2)*DDN_v[5](0,0) + crKeu84;
const double crKeu95 = crKeu23*crKeu94;
const double crKeu96 = DN_p(0,0)*crKeu23;
const double crKeu97 = DN_p(1,0)*crKeu23;
const double crKeu98 = DN_p(2,0)*crKeu23;
const double crKeu99 = C(1,2)*DDN_v[0](1,1);
const double crKeu100 = C(0,1)*DDN_v[0](1,0) + C(0,2)*DDN_v[0](1,0) + C(2,2)*DDN_v[0](1,1) + crKeu99;
const double crKeu101 = crKeu100*crKeu23;
const double crKeu102 = 1.0*C(1,1);
const double crKeu103 = 1.0*DDN_v[0](1,0);
const double crKeu104 = C(1,2)*crKeu103 + C(2,2)*crKeu103 + DDN_v[0](1,1)*crKeu102 + crKeu12 + 1.0*crKeu99;
const double crKeu105 = crKeu104*crKeu23;
const double crKeu106 = C(1,2)*DDN_v[1](1,1);
const double crKeu107 = C(0,1)*DDN_v[1](1,0) + C(0,2)*DDN_v[1](1,0) + C(2,2)*DDN_v[1](1,1) + crKeu106;
const double crKeu108 = crKeu107*crKeu23;
const double crKeu109 = 1.0*DDN_v[1](1,0);
const double crKeu110 = C(1,2)*crKeu109 + C(2,2)*crKeu109 + DDN_v[1](1,1)*crKeu102 + 1.0*crKeu106 + crKeu32;
const double crKeu111 = crKeu110*crKeu23;
const double crKeu112 = C(1,2)*DDN_v[2](1,1);
const double crKeu113 = C(0,1)*DDN_v[2](1,0) + C(0,2)*DDN_v[2](1,0) + C(2,2)*DDN_v[2](1,1) + crKeu112;
const double crKeu114 = crKeu113*crKeu23;
const double crKeu115 = 1.0*DDN_v[2](1,0);
const double crKeu116 = C(1,2)*crKeu115 + C(2,2)*crKeu115 + DDN_v[2](1,1)*crKeu102 + 1.0*crKeu112 + crKeu46;
const double crKeu117 = crKeu116*crKeu23;
const double crKeu118 = C(1,2)*DDN_v[3](1,1);
const double crKeu119 = C(0,1)*DDN_v[3](1,0) + C(0,2)*DDN_v[3](1,0) + C(2,2)*DDN_v[3](1,1) + crKeu118;
const double crKeu120 = crKeu119*crKeu23;
const double crKeu121 = 1.0*DDN_v[3](1,0);
const double crKeu122 = C(1,2)*crKeu121 + C(2,2)*crKeu121 + DDN_v[3](1,1)*crKeu102 + 1.0*crKeu118 + crKeu60;
const double crKeu123 = crKeu122*crKeu23;
const double crKeu124 = C(1,2)*DDN_v[4](1,1);
const double crKeu125 = C(0,1)*DDN_v[4](1,0) + C(0,2)*DDN_v[4](1,0) + C(2,2)*DDN_v[4](1,1) + crKeu124;
const double crKeu126 = crKeu125*crKeu23;
const double crKeu127 = 1.0*DDN_v[4](1,0);
const double crKeu128 = C(1,2)*crKeu127 + C(2,2)*crKeu127 + DDN_v[4](1,1)*crKeu102 + 1.0*crKeu124 + crKeu74;
const double crKeu129 = crKeu128*crKeu23;
const double crKeu130 = C(1,2)*DDN_v[5](1,1);
const double crKeu131 = C(0,1)*DDN_v[5](1,0) + C(0,2)*DDN_v[5](1,0) + C(2,2)*DDN_v[5](1,1) + crKeu130;
const double crKeu132 = crKeu131*crKeu23;
const double crKeu133 = 1.0*DDN_v[5](1,0);
const double crKeu134 = C(1,2)*crKeu133 + C(2,2)*crKeu133 + DDN_v[5](1,1)*crKeu102 + 1.0*crKeu130 + crKeu88;
const double crKeu135 = crKeu134*crKeu23;
const double crKeu136 = DN_p(0,1)*crKeu23;
const double crKeu137 = DN_p(1,1)*crKeu23;
const double crKeu138 = DN_p(2,1)*crKeu23;
rKeu(0,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(0,0) + C(0,2)*DN_v(0,1)) + DN_e(0,1)*crKeu1 - crKeu14*crKeu15 - crKeu16*crKeu17 - crKeu16*crKeu18 + crKeu19);
rKeu(0,1)+=-gauss_weight*(-DN_e(0,0)*(C(0,1)*DN_v(0,1) + crKeu0) - DN_e(0,1)*crKeu21 + crKeu17*crKeu24 + crKeu18*crKeu24 + crKeu22*crKeu25);
rKeu(0,2)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(1,0) + C(0,2)*DN_v(1,1)) + DN_e(0,1)*crKeu27 - crKeu17*crKeu34 - crKeu18*crKeu34 - crKeu25*crKeu33 + crKeu35);
rKeu(0,3)+=-gauss_weight*(-DN_e(0,0)*(C(0,1)*DN_v(1,1) + crKeu26) - DN_e(0,1)*crKeu37 + crKeu17*crKeu39 + crKeu18*crKeu39 + crKeu25*crKeu38);
rKeu(0,4)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(2,0) + C(0,2)*DN_v(2,1)) + DN_e(0,1)*crKeu41 - crKeu17*crKeu48 - crKeu18*crKeu48 - crKeu25*crKeu47 + crKeu49);
rKeu(0,5)+=-gauss_weight*(-DN_e(0,0)*(C(0,1)*DN_v(2,1) + crKeu40) - DN_e(0,1)*crKeu51 + crKeu17*crKeu53 + crKeu18*crKeu53 + crKeu25*crKeu52);
rKeu(0,6)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(3,0) + C(0,2)*DN_v(3,1)) + DN_e(0,1)*crKeu55 - crKeu17*crKeu62 - crKeu18*crKeu62 - crKeu25*crKeu61 + crKeu63);
rKeu(0,7)+=-gauss_weight*(-DN_e(0,0)*(C(0,1)*DN_v(3,1) + crKeu54) - DN_e(0,1)*crKeu65 + crKeu17*crKeu67 + crKeu18*crKeu67 + crKeu25*crKeu66);
rKeu(0,8)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(4,0) + C(0,2)*DN_v(4,1)) + DN_e(0,1)*crKeu69 - crKeu17*crKeu76 - crKeu18*crKeu76 - crKeu25*crKeu75 + crKeu77);
rKeu(0,9)+=-gauss_weight*(-DN_e(0,0)*(C(0,1)*DN_v(4,1) + crKeu68) - DN_e(0,1)*crKeu79 + crKeu17*crKeu81 + crKeu18*crKeu81 + crKeu25*crKeu80);
rKeu(0,10)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(5,0) + C(0,2)*DN_v(5,1)) + DN_e(0,1)*crKeu83 - crKeu17*crKeu90 - crKeu18*crKeu90 - crKeu25*crKeu89 + crKeu91);
rKeu(0,11)+=-gauss_weight*(-DN_e(0,0)*(C(0,1)*DN_v(5,1) + crKeu82) - DN_e(0,1)*crKeu93 + crKeu17*crKeu95 + crKeu18*crKeu95 + crKeu25*crKeu94);
rKeu(0,12)+=gauss_weight*(-DN_e(0,0)*N_p[0] + DN_p(0,0)*crKeu25 + crKeu17*crKeu96 + crKeu18*crKeu96);
rKeu(0,13)+=gauss_weight*(-DN_e(0,0)*N_p[1] + DN_p(1,0)*crKeu25 + crKeu17*crKeu97 + crKeu18*crKeu97);
rKeu(0,14)+=gauss_weight*(-DN_e(0,0)*N_p[2] + DN_p(2,0)*crKeu25 + crKeu17*crKeu98 + crKeu18*crKeu98);
rKeu(1,0)+=-gauss_weight*(-DN_e(0,0)*crKeu1 - DN_e(0,1)*(C(0,1)*DN_v(0,0) + crKeu20) + crKeu100*crKeu25 + crKeu101*crKeu17 + crKeu101*crKeu18);
rKeu(1,1)+=gauss_weight*(DN_e(0,0)*crKeu21 + DN_e(0,1)*(C(1,1)*DN_v(0,1) + C(1,2)*DN_v(0,0)) - crKeu104*crKeu25 - crKeu105*crKeu17 - crKeu105*crKeu18 + crKeu19);
rKeu(1,2)+=-gauss_weight*(-DN_e(0,0)*crKeu27 - DN_e(0,1)*(C(0,1)*DN_v(1,0) + crKeu36) + crKeu107*crKeu25 + crKeu108*crKeu17 + crKeu108*crKeu18);
rKeu(1,3)+=gauss_weight*(DN_e(0,0)*crKeu37 + DN_e(0,1)*(C(1,1)*DN_v(1,1) + C(1,2)*DN_v(1,0)) - crKeu110*crKeu25 - crKeu111*crKeu17 - crKeu111*crKeu18 + crKeu35);
rKeu(1,4)+=-gauss_weight*(-DN_e(0,0)*crKeu41 - DN_e(0,1)*(C(0,1)*DN_v(2,0) + crKeu50) + crKeu113*crKeu25 + crKeu114*crKeu17 + crKeu114*crKeu18);
rKeu(1,5)+=gauss_weight*(DN_e(0,0)*crKeu51 + DN_e(0,1)*(C(1,1)*DN_v(2,1) + C(1,2)*DN_v(2,0)) - crKeu116*crKeu25 - crKeu117*crKeu17 - crKeu117*crKeu18 + crKeu49);
rKeu(1,6)+=-gauss_weight*(-DN_e(0,0)*crKeu55 - DN_e(0,1)*(C(0,1)*DN_v(3,0) + crKeu64) + crKeu119*crKeu25 + crKeu120*crKeu17 + crKeu120*crKeu18);
rKeu(1,7)+=gauss_weight*(DN_e(0,0)*crKeu65 + DN_e(0,1)*(C(1,1)*DN_v(3,1) + C(1,2)*DN_v(3,0)) - crKeu122*crKeu25 - crKeu123*crKeu17 - crKeu123*crKeu18 + crKeu63);
rKeu(1,8)+=-gauss_weight*(-DN_e(0,0)*crKeu69 - DN_e(0,1)*(C(0,1)*DN_v(4,0) + crKeu78) + crKeu125*crKeu25 + crKeu126*crKeu17 + crKeu126*crKeu18);
rKeu(1,9)+=gauss_weight*(DN_e(0,0)*crKeu79 + DN_e(0,1)*(C(1,1)*DN_v(4,1) + C(1,2)*DN_v(4,0)) - crKeu128*crKeu25 - crKeu129*crKeu17 - crKeu129*crKeu18 + crKeu77);
rKeu(1,10)+=-gauss_weight*(-DN_e(0,0)*crKeu83 - DN_e(0,1)*(C(0,1)*DN_v(5,0) + crKeu92) + crKeu131*crKeu25 + crKeu132*crKeu17 + crKeu132*crKeu18);
rKeu(1,11)+=gauss_weight*(DN_e(0,0)*crKeu93 + DN_e(0,1)*(C(1,1)*DN_v(5,1) + C(1,2)*DN_v(5,0)) - crKeu134*crKeu25 - crKeu135*crKeu17 - crKeu135*crKeu18 + crKeu91);
rKeu(1,12)+=gauss_weight*(-DN_e(0,1)*N_p[0] + DN_p(0,1)*crKeu25 + crKeu136*crKeu17 + crKeu136*crKeu18);
rKeu(1,13)+=gauss_weight*(-DN_e(0,1)*N_p[1] + DN_p(1,1)*crKeu25 + crKeu137*crKeu17 + crKeu137*crKeu18);
rKeu(1,14)+=gauss_weight*(-DN_e(0,1)*N_p[2] + DN_p(2,1)*crKeu25 + crKeu138*crKeu17 + crKeu138*crKeu18);


    const double crKee0 = C(0,2)*DN_e(0,0);
const double crKee1 = C(2,2)*DN_e(0,1) + crKee0;
const double crKee2 = DN_v(0,0)*r_v(0,0) + DN_v(1,0)*r_v(1,0) + DN_v(2,0)*r_v(2,0) + DN_v(3,0)*r_v(3,0) + DN_v(4,0)*r_v(4,0) + DN_v(5,0)*r_v(5,0);
const double crKee3 = N_e*crKee2;
const double crKee4 = DN_e(0,1)*v_e(0,1);
const double crKee5 = N_e*crKee4;
const double crKee6 = DN_e(0,0)*v_e(0,0);
const double crKee7 = N_e*crKee6;
const double crKee8 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crKee9 = DN_e(0,0)*crKee8;
const double crKee10 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crKee11 = DN_e(0,1)*crKee10;
const double crKee12 = crKee11 + crKee9;
const double crKee13 = crKee12 + crKee3 + crKee5 + 2*crKee7;
const double crKee14 = N_e*rho;
const double crKee15 = 1.0/(dyn_tau*rho/rData.DeltaTime + rho*stab_c2*sqrt(pow(crKee10, 2) + pow(crKee8, 2))/h + mu*stab_c1/pow(h, 2));
const double crKee16 = crKee15*pow(rho, 2);
const double crKee17 = crKee13*crKee16;
const double crKee18 = N_e*crKee17;
const double crKee19 = 2.0*crKee4 + 2.0*crKee6;
const double crKee20 = 1.0*DN_v(0,0)*vconv(0,0) + 1.0*DN_v(0,1)*vconv(0,1) + 1.0*DN_v(1,0)*vconv(1,0) + 1.0*DN_v(1,1)*vconv(1,1) + 1.0*DN_v(2,0)*vconv(2,0) + 1.0*DN_v(2,1)*vconv(2,1) + 1.0*DN_v(3,0)*vconv(3,0) + 1.0*DN_v(3,1)*vconv(3,1) + 1.0*DN_v(4,0)*vconv(4,0) + 1.0*DN_v(4,1)*vconv(4,1) + 1.0*DN_v(5,0)*vconv(5,0) + 1.0*DN_v(5,1)*vconv(5,1);
const double crKee21 = 1.0*crKee12;
const double crKee22 = 1.0*C(0,0);
const double crKee23 = DDN_v[0](0,0)*r_v(0,0);
const double crKee24 = DDN_v[1](0,0)*r_v(1,0);
const double crKee25 = DDN_v[2](0,0)*r_v(2,0);
const double crKee26 = DDN_v[3](0,0)*r_v(3,0);
const double crKee27 = DDN_v[4](0,0)*r_v(4,0);
const double crKee28 = DDN_v[5](0,0)*r_v(5,0);
const double crKee29 = 1.0*C(0,1);
const double crKee30 = DDN_v[0](0,1)*r_v(0,1);
const double crKee31 = DDN_v[1](0,1)*r_v(1,1);
const double crKee32 = DDN_v[2](0,1)*r_v(2,1);
const double crKee33 = DDN_v[3](0,1)*r_v(3,1);
const double crKee34 = DDN_v[4](0,1)*r_v(4,1);
const double crKee35 = DDN_v[5](0,1)*r_v(5,1);
const double crKee36 = 1.0*C(0,2);
const double crKee37 = 1.0*C(1,2);
const double crKee38 = DDN_v[0](0,0)*r_v(0,1) + DDN_v[0](0,1)*r_v(0,0);
const double crKee39 = DDN_v[1](0,0)*r_v(1,1) + DDN_v[1](0,1)*r_v(1,0);
const double crKee40 = DDN_v[2](0,0)*r_v(2,1) + DDN_v[2](0,1)*r_v(2,0);
const double crKee41 = DDN_v[3](0,0)*r_v(3,1) + DDN_v[3](0,1)*r_v(3,0);
const double crKee42 = DDN_v[4](0,0)*r_v(4,1) + DDN_v[4](0,1)*r_v(4,0);
const double crKee43 = DDN_v[5](0,0)*r_v(5,1) + DDN_v[5](0,1)*r_v(5,0);
const double crKee44 = 1.0*C(2,2);
const double crKee45 = DN_v(0,1)*r_v(0,0) + DN_v(1,1)*r_v(1,0) + DN_v(2,1)*r_v(2,0) + DN_v(3,1)*r_v(3,0) + DN_v(4,1)*r_v(4,0) + DN_v(5,1)*r_v(5,0);
const double crKee46 = N_e*v_e(0,1);
const double crKee47 = 2.0*crKee14*crKee15;
const double crKee48 = crKee47*(-DN_p(0,0)*r_p[0] - DN_p(1,0)*r_p[1] - DN_p(2,0)*r_p[2] + crKee22*crKee23 + crKee22*crKee24 + crKee22*crKee25 + crKee22*crKee26 + crKee22*crKee27 + crKee22*crKee28 + crKee23*crKee36 + crKee24*crKee36 + crKee25*crKee36 + crKee26*crKee36 + crKee27*crKee36 + crKee28*crKee36 + crKee29*crKee30 + crKee29*crKee31 + crKee29*crKee32 + crKee29*crKee33 + crKee29*crKee34 + crKee29*crKee35 + crKee30*crKee37 + crKee31*crKee37 + crKee32*crKee37 + crKee33*crKee37 + crKee34*crKee37 + crKee35*crKee37 + crKee36*crKee38 + crKee36*crKee39 + crKee36*crKee40 + crKee36*crKee41 + crKee36*crKee42 + crKee36*crKee43 + crKee38*crKee44 + crKee39*crKee44 + crKee40*crKee44 + crKee41*crKee44 + crKee42*crKee44 + crKee43*crKee44 + rho*(N_v[0]*r_f(0,0) + N_v[1]*r_f(1,0) + N_v[2]*r_f(2,0) + N_v[3]*r_f(3,0) + N_v[4]*r_f(4,0) + N_v[5]*r_f(5,0)) - rho*(N_v[0]*(rData.BDF0*r_v(0,0) + rData.BDF1*r_vn(0,0) + rData.BDF2*r_vnn(0,0)) + N_v[1]*(rData.BDF0*r_v(1,0) + rData.BDF1*r_vn(1,0) + rData.BDF2*r_vnn(1,0)) + N_v[2]*(rData.BDF0*r_v(2,0) + rData.BDF1*r_vn(2,0) + rData.BDF2*r_vnn(2,0)) + N_v[3]*(rData.BDF0*r_v(3,0) + rData.BDF1*r_vn(3,0) + rData.BDF2*r_vnn(3,0)) + N_v[4]*(rData.BDF0*r_v(4,0) + rData.BDF1*r_vn(4,0) + rData.BDF2*r_vnn(4,0)) + N_v[5]*(rData.BDF0*r_v(5,0) + rData.BDF1*r_vn(5,0) + rData.BDF2*r_vnn(5,0))) - rho*(DN_e(0,0)*N_e*pow(v_e(0,0), 2) + crKee10*crKee45 + crKee11*v_e(0,0) + crKee2*crKee8 + crKee3*v_e(0,0) + crKee45*crKee46 + crKee5*v_e(0,0) + crKee9*v_e(0,0)));
const double crKee49 = C(1,2)*DN_e(0,1);
const double crKee50 = C(2,2)*DN_e(0,0) + crKee49;
const double crKee51 = DN_e(0,1)*v_e(0,0) + crKee45;
const double crKee52 = pow(N_e, 2);
const double crKee53 = crKee52*rho;
const double crKee54 = crKee16*crKee52;
const double crKee55 = crKee51*crKee54;
const double crKee56 = crKee16*crKee21;
const double crKee57 = N_e*crKee56;
const double crKee58 = DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1);
const double crKee59 = DN_e(0,0)*v_e(0,1) + crKee58;
const double crKee60 = crKee54*crKee59;
const double crKee61 = DDN_v[0](1,0)*r_v(0,0);
const double crKee62 = DDN_v[1](1,0)*r_v(1,0);
const double crKee63 = DDN_v[2](1,0)*r_v(2,0);
const double crKee64 = DDN_v[3](1,0)*r_v(3,0);
const double crKee65 = DDN_v[4](1,0)*r_v(4,0);
const double crKee66 = DDN_v[5](1,0)*r_v(5,0);
const double crKee67 = 1.0*C(1,1);
const double crKee68 = DDN_v[0](1,1)*r_v(0,1);
const double crKee69 = DDN_v[1](1,1)*r_v(1,1);
const double crKee70 = DDN_v[2](1,1)*r_v(2,1);
const double crKee71 = DDN_v[3](1,1)*r_v(3,1);
const double crKee72 = DDN_v[4](1,1)*r_v(4,1);
const double crKee73 = DDN_v[5](1,1)*r_v(5,1);
const double crKee74 = DDN_v[0](1,0)*r_v(0,1) + DDN_v[0](1,1)*r_v(0,0);
const double crKee75 = DDN_v[1](1,0)*r_v(1,1) + DDN_v[1](1,1)*r_v(1,0);
const double crKee76 = DDN_v[2](1,0)*r_v(2,1) + DDN_v[2](1,1)*r_v(2,0);
const double crKee77 = DDN_v[3](1,0)*r_v(3,1) + DDN_v[3](1,1)*r_v(3,0);
const double crKee78 = DDN_v[4](1,0)*r_v(4,1) + DDN_v[4](1,1)*r_v(4,0);
const double crKee79 = DDN_v[5](1,0)*r_v(5,1) + DDN_v[5](1,1)*r_v(5,0);
const double crKee80 = DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1);
const double crKee81 = N_e*crKee80;
const double crKee82 = crKee47*(-DN_p(0,1)*r_p[0] - DN_p(1,1)*r_p[1] - DN_p(2,1)*r_p[2] + crKee29*crKee61 + crKee29*crKee62 + crKee29*crKee63 + crKee29*crKee64 + crKee29*crKee65 + crKee29*crKee66 + crKee36*crKee61 + crKee36*crKee62 + crKee36*crKee63 + crKee36*crKee64 + crKee36*crKee65 + crKee36*crKee66 + crKee37*crKee68 + crKee37*crKee69 + crKee37*crKee70 + crKee37*crKee71 + crKee37*crKee72 + crKee37*crKee73 + crKee37*crKee74 + crKee37*crKee75 + crKee37*crKee76 + crKee37*crKee77 + crKee37*crKee78 + crKee37*crKee79 + crKee44*crKee74 + crKee44*crKee75 + crKee44*crKee76 + crKee44*crKee77 + crKee44*crKee78 + crKee44*crKee79 + crKee67*crKee68 + crKee67*crKee69 + crKee67*crKee70 + crKee67*crKee71 + crKee67*crKee72 + crKee67*crKee73 + rho*(N_v[0]*r_f(0,1) + N_v[1]*r_f(1,1) + N_v[2]*r_f(2,1) + N_v[3]*r_f(3,1) + N_v[4]*r_f(4,1) + N_v[5]*r_f(5,1)) - rho*(N_v[0]*(rData.BDF0*r_v(0,1) + rData.BDF1*r_vn(0,1) + rData.BDF2*r_vnn(0,1)) + N_v[1]*(rData.BDF0*r_v(1,1) + rData.BDF1*r_vn(1,1) + rData.BDF2*r_vnn(1,1)) + N_v[2]*(rData.BDF0*r_v(2,1) + rData.BDF1*r_vn(2,1) + rData.BDF2*r_vnn(2,1)) + N_v[3]*(rData.BDF0*r_v(3,1) + rData.BDF1*r_vn(3,1) + rData.BDF2*r_vnn(3,1)) + N_v[4]*(rData.BDF0*r_v(4,1) + rData.BDF1*r_vn(4,1) + rData.BDF2*r_vnn(4,1)) + N_v[5]*(rData.BDF0*r_v(5,1) + rData.BDF1*r_vn(5,1) + rData.BDF2*r_vnn(5,1))) - rho*(DN_e(0,1)*N_e*pow(v_e(0,1), 2) + N_e*crKee58*v_e(0,0) + crKee10*crKee80 + crKee11*v_e(0,1) + crKee46*crKee6 + crKee58*crKee8 + crKee81*v_e(0,1) + crKee9*v_e(0,1)));
const double crKee83 = crKee12 + 2*crKee5 + crKee7 + crKee81;
const double crKee84 = N_e*crKee16*crKee83;
rKee(0,0)+=gauss_weight*(-DN_e(0,0)*crKee48 + DN_e(0,0)*(C(0,0)*DN_e(0,0) + C(0,2)*DN_e(0,1)) + DN_e(0,1)*crKee1 + crKee13*crKee14 + crKee17*crKee21 + crKee18*crKee19 + crKee18*crKee20);
rKee(0,1)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_e(0,1) + crKee0) - DN_e(0,1)*crKee48 + DN_e(0,1)*crKee50 + crKee19*crKee55 + crKee20*crKee55 + crKee51*crKee53 + crKee51*crKee57);
rKee(1,0)+=gauss_weight*(DN_e(0,0)*crKee1 - DN_e(0,0)*crKee82 + DN_e(0,1)*(C(0,1)*DN_e(0,0) + crKee49) + crKee19*crKee60 + crKee20*crKee60 + crKee53*crKee59 + crKee57*crKee59);
rKee(1,1)+=gauss_weight*(DN_e(0,0)*crKee50 - DN_e(0,1)*crKee82 + DN_e(0,1)*(C(1,1)*DN_e(0,1) + C(1,2)*DN_e(0,0)) + crKee14*crKee83 + crKee19*crKee84 + crKee20*crKee84 + crKee56*crKee83);

}

template <>
void IncompressibleNavierStokesDivStable<3>::ComputeGaussPointEnrichmentContribution(
    const ElementDataContainer& rData,
    array_1d<double, 3>& rRHSee,
    BoundedMatrix<double, LocalSize, 3>& rKue,
    BoundedMatrix<double, 3, LocalSize>& rKeu,
    BoundedMatrix<double, 3, 3>& rKee)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double bdf1 = rData.BDF1;
    const double bdf2 = rData.BDF2;
    const double dt = rData.DeltaTime;

    // Get stabilization data
    const double h = rData.ElementSize;
    const double stab_c1 = rData.StabC1;
    const double stab_c2 = rData.StabC2;
    const double dyn_tau = rData.DynamicTau;

    // Get nodal data
    const auto& r_v = rData.Velocity;
    const auto& r_vn = rData.VelocityOld1;
    const auto& r_vnn = rData.VelocityOld2;
    const auto& r_vmesh = rData.MeshVelocity;
    const auto& r_f = rData.BodyForce;
    const auto& r_p = rData.Pressure;

    // Calculate convective velocity
    const BoundedMatrix<double,3,10> vconv = r_v - r_vmesh;

    // Get stress from material response
    const auto& r_stress = rData.ShearStress;
    const auto& C = rData.ConstitutiveMatrix;

    // Get shape function values
    const auto& N_p = rData.N_p;
    const auto& N_v = rData.N_v;
    const auto& DN_p = rData.DN_p;
    const auto& DN_v = rData.DN_v;
    const double N_e = rData.N_e;
    const auto& DN_e = rData.DN_e;
    const auto& DDN_v = rData.DDN_v;

    // Assemble enrichment contribution
    const double gauss_weight = rData.Weight;
    const BoundedMatrix<double,1,3> v_e = ZeroMatrix(1,3);

    //substitute_rhs_ee_3D

    //substitute_K_ue_3D

    //substitute_K_eu_3D

    //substitute_K_ee_3D
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}


template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class IncompressibleNavierStokesDivStable<2>;
template class IncompressibleNavierStokesDivStable<3>;

}