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
    CalculateKinematics(weights, velocity_N, pressure_N, velocity_DN, pressure_DN, velocity_bubble, velocity_bubble_grad);

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
    rData.StabC1 = 4.0;
    rData.StabC2 = 2.0;
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
    std::vector<BoundedMatrix<double, 1, TDim>>& rVelocityBubbleGrad)
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
            rVelocityBubble[g] = r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2];
            r_vel_bub_grad_g(0,0) = r_bar_coords_grads(0,0) * r_bar_coords[1] * r_bar_coords[2] + r_bar_coords[0] * r_bar_coords_grads(1,0) * r_bar_coords[2] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,0);
            r_vel_bub_grad_g(0,1) = r_bar_coords_grads(0,1) * r_bar_coords[1] * r_bar_coords[2] + r_bar_coords[0] * r_bar_coords_grads(1,1) * r_bar_coords[2] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,1);
        } else {
            rVelocityBubble[g] = r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords[3];
            r_vel_bub_grad_g(0,0) = r_bar_coords_grads(0,0) * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords_grads(1,0) * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,0) * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords_grads(3,0);
            r_vel_bub_grad_g(0,1) = r_bar_coords_grads(0,1) * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords_grads(1,1) * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,1) * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords_grads(3,1);
            r_vel_bub_grad_g(0,2) = r_bar_coords_grads(0,2) * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords_grads(1,2) * r_bar_coords[2] * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords_grads(2,2) * r_bar_coords[3] + r_bar_coords[0] * r_bar_coords[1] * r_bar_coords[2] * r_bar_coords_grads(3,2);
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

    // Assemble LHS contribution
    const double gauss_weight = rData.Weight;
    const double crLHS0 = C(0,0)*DN_v(0,0) + C(0,2)*DN_v(0,1);
const double crLHS1 = C(0,2)*DN_v(0,0);
const double crLHS2 = C(2,2)*DN_v(0,1) + crLHS1;
const double crLHS3 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crLHS4 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crLHS5 = DN_v(0,0)*crLHS3 + DN_v(0,1)*crLHS4;
const double crLHS6 = N_v[0]*rho;
const double crLHS7 = rData.BDF0*rho;
const double crLHS8 = N_v[0]*rData.BDF0;
const double crLHS9 = crLHS5 + crLHS8;
const double crLHS10 = pow(rho, 2);
const double crLHS11 = 1.0/(dyn_tau*rho/rData.DeltaTime + rho*stab_c2*sqrt(pow(crLHS3, 2) + pow(crLHS4, 2))/h + mu*stab_c1/pow(h, 2));
const double crLHS12 = crLHS11*crLHS5;
const double crLHS13 = crLHS10*crLHS12;
const double crLHS14 = crLHS11*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1));
const double crLHS15 = crLHS10*crLHS14;
const double crLHS16 = N_v[0]*crLHS15;
const double crLHS17 = pow(N_v[0], 2)*crLHS7 + crLHS13*crLHS9 + crLHS16*crLHS9 + crLHS5*crLHS6;
const double crLHS18 = C(0,1)*DN_v(0,1) + crLHS1;
const double crLHS19 = C(1,2)*DN_v(0,1);
const double crLHS20 = C(2,2)*DN_v(0,0) + crLHS19;
const double crLHS21 = C(0,0)*DN_v(1,0) + C(0,2)*DN_v(1,1);
const double crLHS22 = C(0,2)*DN_v(1,0);
const double crLHS23 = C(2,2)*DN_v(1,1) + crLHS22;
const double crLHS24 = DN_v(1,0)*crLHS3 + DN_v(1,1)*crLHS4;
const double crLHS25 = N_v[1]*rho;
const double crLHS26 = crLHS25*crLHS8;
const double crLHS27 = N_v[1]*rData.BDF0;
const double crLHS28 = crLHS24 + crLHS27;
const double crLHS29 = crLHS13*crLHS28 + crLHS16*crLHS28 + crLHS24*crLHS6 + crLHS26;
const double crLHS30 = C(0,1)*DN_v(1,1) + crLHS22;
const double crLHS31 = C(1,2)*DN_v(1,1);
const double crLHS32 = C(2,2)*DN_v(1,0) + crLHS31;
const double crLHS33 = C(0,0)*DN_v(2,0) + C(0,2)*DN_v(2,1);
const double crLHS34 = C(0,2)*DN_v(2,0);
const double crLHS35 = C(2,2)*DN_v(2,1) + crLHS34;
const double crLHS36 = DN_v(2,0)*crLHS3 + DN_v(2,1)*crLHS4;
const double crLHS37 = N_v[2]*rho;
const double crLHS38 = crLHS37*crLHS8;
const double crLHS39 = N_v[2]*rData.BDF0;
const double crLHS40 = crLHS36 + crLHS39;
const double crLHS41 = crLHS13*crLHS40 + crLHS16*crLHS40 + crLHS36*crLHS6 + crLHS38;
const double crLHS42 = C(0,1)*DN_v(2,1) + crLHS34;
const double crLHS43 = C(1,2)*DN_v(2,1);
const double crLHS44 = C(2,2)*DN_v(2,0) + crLHS43;
const double crLHS45 = C(0,0)*DN_v(3,0) + C(0,2)*DN_v(3,1);
const double crLHS46 = C(0,2)*DN_v(3,0);
const double crLHS47 = C(2,2)*DN_v(3,1) + crLHS46;
const double crLHS48 = DN_v(3,0)*crLHS3 + DN_v(3,1)*crLHS4;
const double crLHS49 = N_v[3]*rho;
const double crLHS50 = crLHS49*crLHS8;
const double crLHS51 = N_v[3]*rData.BDF0;
const double crLHS52 = crLHS48 + crLHS51;
const double crLHS53 = crLHS13*crLHS52 + crLHS16*crLHS52 + crLHS48*crLHS6 + crLHS50;
const double crLHS54 = C(0,1)*DN_v(3,1) + crLHS46;
const double crLHS55 = C(1,2)*DN_v(3,1);
const double crLHS56 = C(2,2)*DN_v(3,0) + crLHS55;
const double crLHS57 = C(0,0)*DN_v(4,0) + C(0,2)*DN_v(4,1);
const double crLHS58 = C(0,2)*DN_v(4,0);
const double crLHS59 = C(2,2)*DN_v(4,1) + crLHS58;
const double crLHS60 = DN_v(4,0)*crLHS3 + DN_v(4,1)*crLHS4;
const double crLHS61 = N_v[4]*rho;
const double crLHS62 = crLHS61*crLHS8;
const double crLHS63 = N_v[4]*rData.BDF0;
const double crLHS64 = crLHS60 + crLHS63;
const double crLHS65 = crLHS13*crLHS64 + crLHS16*crLHS64 + crLHS6*crLHS60 + crLHS62;
const double crLHS66 = C(0,1)*DN_v(4,1) + crLHS58;
const double crLHS67 = C(1,2)*DN_v(4,1);
const double crLHS68 = C(2,2)*DN_v(4,0) + crLHS67;
const double crLHS69 = C(0,0)*DN_v(5,0) + C(0,2)*DN_v(5,1);
const double crLHS70 = C(0,2)*DN_v(5,0);
const double crLHS71 = C(2,2)*DN_v(5,1) + crLHS70;
const double crLHS72 = DN_v(5,0)*crLHS3 + DN_v(5,1)*crLHS4;
const double crLHS73 = N_v[5]*rho;
const double crLHS74 = crLHS73*crLHS8;
const double crLHS75 = N_v[5]*rData.BDF0 + crLHS72;
const double crLHS76 = crLHS13*crLHS75 + crLHS16*crLHS75 + crLHS6*crLHS72 + crLHS74;
const double crLHS77 = C(0,1)*DN_v(5,1) + crLHS70;
const double crLHS78 = C(1,2)*DN_v(5,1);
const double crLHS79 = C(2,2)*DN_v(5,0) + crLHS78;
const double crLHS80 = DN_v(0,0)*N_p[0];
const double crLHS81 = crLHS14*crLHS6;
const double crLHS82 = crLHS12*rho;
const double crLHS83 = DN_v(0,0)*N_p[1];
const double crLHS84 = DN_v(0,0)*N_p[2];
const double crLHS85 = C(0,1)*DN_v(0,0) + crLHS19;
const double crLHS86 = C(1,1)*DN_v(0,1) + C(1,2)*DN_v(0,0);
const double crLHS87 = C(0,1)*DN_v(1,0) + crLHS31;
const double crLHS88 = C(1,1)*DN_v(1,1) + C(1,2)*DN_v(1,0);
const double crLHS89 = C(0,1)*DN_v(2,0) + crLHS43;
const double crLHS90 = C(1,1)*DN_v(2,1) + C(1,2)*DN_v(2,0);
const double crLHS91 = C(0,1)*DN_v(3,0) + crLHS55;
const double crLHS92 = C(1,1)*DN_v(3,1) + C(1,2)*DN_v(3,0);
const double crLHS93 = C(0,1)*DN_v(4,0) + crLHS67;
const double crLHS94 = C(1,1)*DN_v(4,1) + C(1,2)*DN_v(4,0);
const double crLHS95 = C(0,1)*DN_v(5,0) + crLHS78;
const double crLHS96 = C(1,1)*DN_v(5,1) + C(1,2)*DN_v(5,0);
const double crLHS97 = DN_v(0,1)*N_p[0];
const double crLHS98 = DN_v(0,1)*N_p[1];
const double crLHS99 = DN_v(0,1)*N_p[2];
const double crLHS100 = crLHS11*crLHS24;
const double crLHS101 = crLHS10*crLHS100;
const double crLHS102 = N_v[1]*crLHS15;
const double crLHS103 = crLHS101*crLHS9 + crLHS102*crLHS9 + crLHS25*crLHS5 + crLHS26;
const double crLHS104 = pow(N_v[1], 2)*crLHS7 + crLHS101*crLHS28 + crLHS102*crLHS28 + crLHS24*crLHS25;
const double crLHS105 = crLHS27*crLHS37;
const double crLHS106 = crLHS101*crLHS40 + crLHS102*crLHS40 + crLHS105 + crLHS25*crLHS36;
const double crLHS107 = crLHS27*crLHS49;
const double crLHS108 = crLHS101*crLHS52 + crLHS102*crLHS52 + crLHS107 + crLHS25*crLHS48;
const double crLHS109 = crLHS27*crLHS61;
const double crLHS110 = crLHS101*crLHS64 + crLHS102*crLHS64 + crLHS109 + crLHS25*crLHS60;
const double crLHS111 = crLHS27*crLHS73;
const double crLHS112 = crLHS101*crLHS75 + crLHS102*crLHS75 + crLHS111 + crLHS25*crLHS72;
const double crLHS113 = DN_v(1,0)*N_p[0];
const double crLHS114 = crLHS14*crLHS25;
const double crLHS115 = crLHS100*rho;
const double crLHS116 = DN_v(1,0)*N_p[1];
const double crLHS117 = DN_v(1,0)*N_p[2];
const double crLHS118 = DN_v(1,1)*N_p[0];
const double crLHS119 = DN_v(1,1)*N_p[1];
const double crLHS120 = DN_v(1,1)*N_p[2];
const double crLHS121 = crLHS11*crLHS36;
const double crLHS122 = crLHS10*crLHS121;
const double crLHS123 = N_v[2]*crLHS15;
const double crLHS124 = crLHS122*crLHS9 + crLHS123*crLHS9 + crLHS37*crLHS5 + crLHS38;
const double crLHS125 = crLHS105 + crLHS122*crLHS28 + crLHS123*crLHS28 + crLHS24*crLHS37;
const double crLHS126 = pow(N_v[2], 2)*crLHS7 + crLHS122*crLHS40 + crLHS123*crLHS40 + crLHS36*crLHS37;
const double crLHS127 = crLHS39*crLHS49;
const double crLHS128 = crLHS122*crLHS52 + crLHS123*crLHS52 + crLHS127 + crLHS37*crLHS48;
const double crLHS129 = crLHS39*crLHS61;
const double crLHS130 = crLHS122*crLHS64 + crLHS123*crLHS64 + crLHS129 + crLHS37*crLHS60;
const double crLHS131 = crLHS39*crLHS73;
const double crLHS132 = crLHS122*crLHS75 + crLHS123*crLHS75 + crLHS131 + crLHS37*crLHS72;
const double crLHS133 = DN_v(2,0)*N_p[0];
const double crLHS134 = crLHS14*crLHS37;
const double crLHS135 = crLHS121*rho;
const double crLHS136 = DN_v(2,0)*N_p[1];
const double crLHS137 = DN_v(2,0)*N_p[2];
const double crLHS138 = DN_v(2,1)*N_p[0];
const double crLHS139 = DN_v(2,1)*N_p[1];
const double crLHS140 = DN_v(2,1)*N_p[2];
const double crLHS141 = crLHS11*crLHS48;
const double crLHS142 = crLHS10*crLHS141;
const double crLHS143 = N_v[3]*crLHS15;
const double crLHS144 = crLHS142*crLHS9 + crLHS143*crLHS9 + crLHS49*crLHS5 + crLHS50;
const double crLHS145 = crLHS107 + crLHS142*crLHS28 + crLHS143*crLHS28 + crLHS24*crLHS49;
const double crLHS146 = crLHS127 + crLHS142*crLHS40 + crLHS143*crLHS40 + crLHS36*crLHS49;
const double crLHS147 = pow(N_v[3], 2)*crLHS7 + crLHS142*crLHS52 + crLHS143*crLHS52 + crLHS48*crLHS49;
const double crLHS148 = crLHS51*crLHS61;
const double crLHS149 = crLHS142*crLHS64 + crLHS143*crLHS64 + crLHS148 + crLHS49*crLHS60;
const double crLHS150 = crLHS51*crLHS73;
const double crLHS151 = crLHS142*crLHS75 + crLHS143*crLHS75 + crLHS150 + crLHS49*crLHS72;
const double crLHS152 = DN_v(3,0)*N_p[0];
const double crLHS153 = crLHS14*crLHS49;
const double crLHS154 = crLHS141*rho;
const double crLHS155 = DN_v(3,0)*N_p[1];
const double crLHS156 = DN_v(3,0)*N_p[2];
const double crLHS157 = DN_v(3,1)*N_p[0];
const double crLHS158 = DN_v(3,1)*N_p[1];
const double crLHS159 = DN_v(3,1)*N_p[2];
const double crLHS160 = crLHS11*crLHS60;
const double crLHS161 = crLHS10*crLHS160;
const double crLHS162 = N_v[4]*crLHS15;
const double crLHS163 = crLHS161*crLHS9 + crLHS162*crLHS9 + crLHS5*crLHS61 + crLHS62;
const double crLHS164 = crLHS109 + crLHS161*crLHS28 + crLHS162*crLHS28 + crLHS24*crLHS61;
const double crLHS165 = crLHS129 + crLHS161*crLHS40 + crLHS162*crLHS40 + crLHS36*crLHS61;
const double crLHS166 = crLHS148 + crLHS161*crLHS52 + crLHS162*crLHS52 + crLHS48*crLHS61;
const double crLHS167 = pow(N_v[4], 2)*crLHS7 + crLHS161*crLHS64 + crLHS162*crLHS64 + crLHS60*crLHS61;
const double crLHS168 = crLHS63*crLHS73;
const double crLHS169 = crLHS161*crLHS75 + crLHS162*crLHS75 + crLHS168 + crLHS61*crLHS72;
const double crLHS170 = DN_v(4,0)*N_p[0];
const double crLHS171 = crLHS14*crLHS61;
const double crLHS172 = crLHS160*rho;
const double crLHS173 = DN_v(4,0)*N_p[1];
const double crLHS174 = DN_v(4,0)*N_p[2];
const double crLHS175 = DN_v(4,1)*N_p[0];
const double crLHS176 = DN_v(4,1)*N_p[1];
const double crLHS177 = DN_v(4,1)*N_p[2];
const double crLHS178 = crLHS11*crLHS72;
const double crLHS179 = crLHS10*crLHS178;
const double crLHS180 = N_v[5]*crLHS15;
const double crLHS181 = crLHS179*crLHS9 + crLHS180*crLHS9 + crLHS5*crLHS73 + crLHS74;
const double crLHS182 = crLHS111 + crLHS179*crLHS28 + crLHS180*crLHS28 + crLHS24*crLHS73;
const double crLHS183 = crLHS131 + crLHS179*crLHS40 + crLHS180*crLHS40 + crLHS36*crLHS73;
const double crLHS184 = crLHS150 + crLHS179*crLHS52 + crLHS180*crLHS52 + crLHS48*crLHS73;
const double crLHS185 = crLHS168 + crLHS179*crLHS64 + crLHS180*crLHS64 + crLHS60*crLHS73;
const double crLHS186 = pow(N_v[5], 2)*crLHS7 + crLHS179*crLHS75 + crLHS180*crLHS75 + crLHS72*crLHS73;
const double crLHS187 = DN_v(5,0)*N_p[0];
const double crLHS188 = crLHS14*crLHS73;
const double crLHS189 = crLHS178*rho;
const double crLHS190 = DN_v(5,0)*N_p[1];
const double crLHS191 = DN_v(5,0)*N_p[2];
const double crLHS192 = DN_v(5,1)*N_p[0];
const double crLHS193 = DN_v(5,1)*N_p[1];
const double crLHS194 = DN_v(5,1)*N_p[2];
const double crLHS195 = crLHS11*rho;
const double crLHS196 = crLHS195*crLHS9;
const double crLHS197 = crLHS195*crLHS28;
const double crLHS198 = crLHS195*crLHS40;
const double crLHS199 = crLHS195*crLHS52;
const double crLHS200 = crLHS195*crLHS64;
const double crLHS201 = crLHS195*crLHS75;
const double crLHS202 = crLHS11*gauss_weight;
const double crLHS203 = crLHS202*(DN_p(0,0)*DN_p(1,0) + DN_p(0,1)*DN_p(1,1));
const double crLHS204 = crLHS202*(DN_p(0,0)*DN_p(2,0) + DN_p(0,1)*DN_p(2,1));
const double crLHS205 = crLHS202*(DN_p(1,0)*DN_p(2,0) + DN_p(1,1)*DN_p(2,1));
rLHS(0,0)+=gauss_weight*(DN_v(0,0)*crLHS0 + DN_v(0,1)*crLHS2 + crLHS17);
rLHS(0,1)+=gauss_weight*(DN_v(0,0)*crLHS18 + DN_v(0,1)*crLHS20);
rLHS(0,2)+=gauss_weight*(DN_v(0,0)*crLHS21 + DN_v(0,1)*crLHS23 + crLHS29);
rLHS(0,3)+=gauss_weight*(DN_v(0,0)*crLHS30 + DN_v(0,1)*crLHS32);
rLHS(0,4)+=gauss_weight*(DN_v(0,0)*crLHS33 + DN_v(0,1)*crLHS35 + crLHS41);
rLHS(0,5)+=gauss_weight*(DN_v(0,0)*crLHS42 + DN_v(0,1)*crLHS44);
rLHS(0,6)+=gauss_weight*(DN_v(0,0)*crLHS45 + DN_v(0,1)*crLHS47 + crLHS53);
rLHS(0,7)+=gauss_weight*(DN_v(0,0)*crLHS54 + DN_v(0,1)*crLHS56);
rLHS(0,8)+=gauss_weight*(DN_v(0,0)*crLHS57 + DN_v(0,1)*crLHS59 + crLHS65);
rLHS(0,9)+=gauss_weight*(DN_v(0,0)*crLHS66 + DN_v(0,1)*crLHS68);
rLHS(0,10)+=gauss_weight*(DN_v(0,0)*crLHS69 + DN_v(0,1)*crLHS71 + crLHS76);
rLHS(0,11)+=gauss_weight*(DN_v(0,0)*crLHS77 + DN_v(0,1)*crLHS79);
rLHS(0,12)+=gauss_weight*(DN_p(0,0)*crLHS81 + DN_p(0,0)*crLHS82 - crLHS80);
rLHS(0,13)+=gauss_weight*(DN_p(1,0)*crLHS81 + DN_p(1,0)*crLHS82 - crLHS83);
rLHS(0,14)+=gauss_weight*(DN_p(2,0)*crLHS81 + DN_p(2,0)*crLHS82 - crLHS84);
rLHS(1,0)+=gauss_weight*(DN_v(0,0)*crLHS2 + DN_v(0,1)*crLHS85);
rLHS(1,1)+=gauss_weight*(DN_v(0,0)*crLHS20 + DN_v(0,1)*crLHS86 + crLHS17);
rLHS(1,2)+=gauss_weight*(DN_v(0,0)*crLHS23 + DN_v(0,1)*crLHS87);
rLHS(1,3)+=gauss_weight*(DN_v(0,0)*crLHS32 + DN_v(0,1)*crLHS88 + crLHS29);
rLHS(1,4)+=gauss_weight*(DN_v(0,0)*crLHS35 + DN_v(0,1)*crLHS89);
rLHS(1,5)+=gauss_weight*(DN_v(0,0)*crLHS44 + DN_v(0,1)*crLHS90 + crLHS41);
rLHS(1,6)+=gauss_weight*(DN_v(0,0)*crLHS47 + DN_v(0,1)*crLHS91);
rLHS(1,7)+=gauss_weight*(DN_v(0,0)*crLHS56 + DN_v(0,1)*crLHS92 + crLHS53);
rLHS(1,8)+=gauss_weight*(DN_v(0,0)*crLHS59 + DN_v(0,1)*crLHS93);
rLHS(1,9)+=gauss_weight*(DN_v(0,0)*crLHS68 + DN_v(0,1)*crLHS94 + crLHS65);
rLHS(1,10)+=gauss_weight*(DN_v(0,0)*crLHS71 + DN_v(0,1)*crLHS95);
rLHS(1,11)+=gauss_weight*(DN_v(0,0)*crLHS79 + DN_v(0,1)*crLHS96 + crLHS76);
rLHS(1,12)+=gauss_weight*(DN_p(0,1)*crLHS81 + DN_p(0,1)*crLHS82 - crLHS97);
rLHS(1,13)+=gauss_weight*(DN_p(1,1)*crLHS81 + DN_p(1,1)*crLHS82 - crLHS98);
rLHS(1,14)+=gauss_weight*(DN_p(2,1)*crLHS81 + DN_p(2,1)*crLHS82 - crLHS99);
rLHS(2,0)+=gauss_weight*(DN_v(1,0)*crLHS0 + DN_v(1,1)*crLHS2 + crLHS103);
rLHS(2,1)+=gauss_weight*(DN_v(1,0)*crLHS18 + DN_v(1,1)*crLHS20);
rLHS(2,2)+=gauss_weight*(DN_v(1,0)*crLHS21 + DN_v(1,1)*crLHS23 + crLHS104);
rLHS(2,3)+=gauss_weight*(DN_v(1,0)*crLHS30 + DN_v(1,1)*crLHS32);
rLHS(2,4)+=gauss_weight*(DN_v(1,0)*crLHS33 + DN_v(1,1)*crLHS35 + crLHS106);
rLHS(2,5)+=gauss_weight*(DN_v(1,0)*crLHS42 + DN_v(1,1)*crLHS44);
rLHS(2,6)+=gauss_weight*(DN_v(1,0)*crLHS45 + DN_v(1,1)*crLHS47 + crLHS108);
rLHS(2,7)+=gauss_weight*(DN_v(1,0)*crLHS54 + DN_v(1,1)*crLHS56);
rLHS(2,8)+=gauss_weight*(DN_v(1,0)*crLHS57 + DN_v(1,1)*crLHS59 + crLHS110);
rLHS(2,9)+=gauss_weight*(DN_v(1,0)*crLHS66 + DN_v(1,1)*crLHS68);
rLHS(2,10)+=gauss_weight*(DN_v(1,0)*crLHS69 + DN_v(1,1)*crLHS71 + crLHS112);
rLHS(2,11)+=gauss_weight*(DN_v(1,0)*crLHS77 + DN_v(1,1)*crLHS79);
rLHS(2,12)+=gauss_weight*(DN_p(0,0)*crLHS114 + DN_p(0,0)*crLHS115 - crLHS113);
rLHS(2,13)+=gauss_weight*(DN_p(1,0)*crLHS114 + DN_p(1,0)*crLHS115 - crLHS116);
rLHS(2,14)+=gauss_weight*(DN_p(2,0)*crLHS114 + DN_p(2,0)*crLHS115 - crLHS117);
rLHS(3,0)+=gauss_weight*(DN_v(1,0)*crLHS2 + DN_v(1,1)*crLHS85);
rLHS(3,1)+=gauss_weight*(DN_v(1,0)*crLHS20 + DN_v(1,1)*crLHS86 + crLHS103);
rLHS(3,2)+=gauss_weight*(DN_v(1,0)*crLHS23 + DN_v(1,1)*crLHS87);
rLHS(3,3)+=gauss_weight*(DN_v(1,0)*crLHS32 + DN_v(1,1)*crLHS88 + crLHS104);
rLHS(3,4)+=gauss_weight*(DN_v(1,0)*crLHS35 + DN_v(1,1)*crLHS89);
rLHS(3,5)+=gauss_weight*(DN_v(1,0)*crLHS44 + DN_v(1,1)*crLHS90 + crLHS106);
rLHS(3,6)+=gauss_weight*(DN_v(1,0)*crLHS47 + DN_v(1,1)*crLHS91);
rLHS(3,7)+=gauss_weight*(DN_v(1,0)*crLHS56 + DN_v(1,1)*crLHS92 + crLHS108);
rLHS(3,8)+=gauss_weight*(DN_v(1,0)*crLHS59 + DN_v(1,1)*crLHS93);
rLHS(3,9)+=gauss_weight*(DN_v(1,0)*crLHS68 + DN_v(1,1)*crLHS94 + crLHS110);
rLHS(3,10)+=gauss_weight*(DN_v(1,0)*crLHS71 + DN_v(1,1)*crLHS95);
rLHS(3,11)+=gauss_weight*(DN_v(1,0)*crLHS79 + DN_v(1,1)*crLHS96 + crLHS112);
rLHS(3,12)+=gauss_weight*(DN_p(0,1)*crLHS114 + DN_p(0,1)*crLHS115 - crLHS118);
rLHS(3,13)+=gauss_weight*(DN_p(1,1)*crLHS114 + DN_p(1,1)*crLHS115 - crLHS119);
rLHS(3,14)+=gauss_weight*(DN_p(2,1)*crLHS114 + DN_p(2,1)*crLHS115 - crLHS120);
rLHS(4,0)+=gauss_weight*(DN_v(2,0)*crLHS0 + DN_v(2,1)*crLHS2 + crLHS124);
rLHS(4,1)+=gauss_weight*(DN_v(2,0)*crLHS18 + DN_v(2,1)*crLHS20);
rLHS(4,2)+=gauss_weight*(DN_v(2,0)*crLHS21 + DN_v(2,1)*crLHS23 + crLHS125);
rLHS(4,3)+=gauss_weight*(DN_v(2,0)*crLHS30 + DN_v(2,1)*crLHS32);
rLHS(4,4)+=gauss_weight*(DN_v(2,0)*crLHS33 + DN_v(2,1)*crLHS35 + crLHS126);
rLHS(4,5)+=gauss_weight*(DN_v(2,0)*crLHS42 + DN_v(2,1)*crLHS44);
rLHS(4,6)+=gauss_weight*(DN_v(2,0)*crLHS45 + DN_v(2,1)*crLHS47 + crLHS128);
rLHS(4,7)+=gauss_weight*(DN_v(2,0)*crLHS54 + DN_v(2,1)*crLHS56);
rLHS(4,8)+=gauss_weight*(DN_v(2,0)*crLHS57 + DN_v(2,1)*crLHS59 + crLHS130);
rLHS(4,9)+=gauss_weight*(DN_v(2,0)*crLHS66 + DN_v(2,1)*crLHS68);
rLHS(4,10)+=gauss_weight*(DN_v(2,0)*crLHS69 + DN_v(2,1)*crLHS71 + crLHS132);
rLHS(4,11)+=gauss_weight*(DN_v(2,0)*crLHS77 + DN_v(2,1)*crLHS79);
rLHS(4,12)+=gauss_weight*(DN_p(0,0)*crLHS134 + DN_p(0,0)*crLHS135 - crLHS133);
rLHS(4,13)+=gauss_weight*(DN_p(1,0)*crLHS134 + DN_p(1,0)*crLHS135 - crLHS136);
rLHS(4,14)+=gauss_weight*(DN_p(2,0)*crLHS134 + DN_p(2,0)*crLHS135 - crLHS137);
rLHS(5,0)+=gauss_weight*(DN_v(2,0)*crLHS2 + DN_v(2,1)*crLHS85);
rLHS(5,1)+=gauss_weight*(DN_v(2,0)*crLHS20 + DN_v(2,1)*crLHS86 + crLHS124);
rLHS(5,2)+=gauss_weight*(DN_v(2,0)*crLHS23 + DN_v(2,1)*crLHS87);
rLHS(5,3)+=gauss_weight*(DN_v(2,0)*crLHS32 + DN_v(2,1)*crLHS88 + crLHS125);
rLHS(5,4)+=gauss_weight*(DN_v(2,0)*crLHS35 + DN_v(2,1)*crLHS89);
rLHS(5,5)+=gauss_weight*(DN_v(2,0)*crLHS44 + DN_v(2,1)*crLHS90 + crLHS126);
rLHS(5,6)+=gauss_weight*(DN_v(2,0)*crLHS47 + DN_v(2,1)*crLHS91);
rLHS(5,7)+=gauss_weight*(DN_v(2,0)*crLHS56 + DN_v(2,1)*crLHS92 + crLHS128);
rLHS(5,8)+=gauss_weight*(DN_v(2,0)*crLHS59 + DN_v(2,1)*crLHS93);
rLHS(5,9)+=gauss_weight*(DN_v(2,0)*crLHS68 + DN_v(2,1)*crLHS94 + crLHS130);
rLHS(5,10)+=gauss_weight*(DN_v(2,0)*crLHS71 + DN_v(2,1)*crLHS95);
rLHS(5,11)+=gauss_weight*(DN_v(2,0)*crLHS79 + DN_v(2,1)*crLHS96 + crLHS132);
rLHS(5,12)+=gauss_weight*(DN_p(0,1)*crLHS134 + DN_p(0,1)*crLHS135 - crLHS138);
rLHS(5,13)+=gauss_weight*(DN_p(1,1)*crLHS134 + DN_p(1,1)*crLHS135 - crLHS139);
rLHS(5,14)+=gauss_weight*(DN_p(2,1)*crLHS134 + DN_p(2,1)*crLHS135 - crLHS140);
rLHS(6,0)+=gauss_weight*(DN_v(3,0)*crLHS0 + DN_v(3,1)*crLHS2 + crLHS144);
rLHS(6,1)+=gauss_weight*(DN_v(3,0)*crLHS18 + DN_v(3,1)*crLHS20);
rLHS(6,2)+=gauss_weight*(DN_v(3,0)*crLHS21 + DN_v(3,1)*crLHS23 + crLHS145);
rLHS(6,3)+=gauss_weight*(DN_v(3,0)*crLHS30 + DN_v(3,1)*crLHS32);
rLHS(6,4)+=gauss_weight*(DN_v(3,0)*crLHS33 + DN_v(3,1)*crLHS35 + crLHS146);
rLHS(6,5)+=gauss_weight*(DN_v(3,0)*crLHS42 + DN_v(3,1)*crLHS44);
rLHS(6,6)+=gauss_weight*(DN_v(3,0)*crLHS45 + DN_v(3,1)*crLHS47 + crLHS147);
rLHS(6,7)+=gauss_weight*(DN_v(3,0)*crLHS54 + DN_v(3,1)*crLHS56);
rLHS(6,8)+=gauss_weight*(DN_v(3,0)*crLHS57 + DN_v(3,1)*crLHS59 + crLHS149);
rLHS(6,9)+=gauss_weight*(DN_v(3,0)*crLHS66 + DN_v(3,1)*crLHS68);
rLHS(6,10)+=gauss_weight*(DN_v(3,0)*crLHS69 + DN_v(3,1)*crLHS71 + crLHS151);
rLHS(6,11)+=gauss_weight*(DN_v(3,0)*crLHS77 + DN_v(3,1)*crLHS79);
rLHS(6,12)+=gauss_weight*(DN_p(0,0)*crLHS153 + DN_p(0,0)*crLHS154 - crLHS152);
rLHS(6,13)+=gauss_weight*(DN_p(1,0)*crLHS153 + DN_p(1,0)*crLHS154 - crLHS155);
rLHS(6,14)+=gauss_weight*(DN_p(2,0)*crLHS153 + DN_p(2,0)*crLHS154 - crLHS156);
rLHS(7,0)+=gauss_weight*(DN_v(3,0)*crLHS2 + DN_v(3,1)*crLHS85);
rLHS(7,1)+=gauss_weight*(DN_v(3,0)*crLHS20 + DN_v(3,1)*crLHS86 + crLHS144);
rLHS(7,2)+=gauss_weight*(DN_v(3,0)*crLHS23 + DN_v(3,1)*crLHS87);
rLHS(7,3)+=gauss_weight*(DN_v(3,0)*crLHS32 + DN_v(3,1)*crLHS88 + crLHS145);
rLHS(7,4)+=gauss_weight*(DN_v(3,0)*crLHS35 + DN_v(3,1)*crLHS89);
rLHS(7,5)+=gauss_weight*(DN_v(3,0)*crLHS44 + DN_v(3,1)*crLHS90 + crLHS146);
rLHS(7,6)+=gauss_weight*(DN_v(3,0)*crLHS47 + DN_v(3,1)*crLHS91);
rLHS(7,7)+=gauss_weight*(DN_v(3,0)*crLHS56 + DN_v(3,1)*crLHS92 + crLHS147);
rLHS(7,8)+=gauss_weight*(DN_v(3,0)*crLHS59 + DN_v(3,1)*crLHS93);
rLHS(7,9)+=gauss_weight*(DN_v(3,0)*crLHS68 + DN_v(3,1)*crLHS94 + crLHS149);
rLHS(7,10)+=gauss_weight*(DN_v(3,0)*crLHS71 + DN_v(3,1)*crLHS95);
rLHS(7,11)+=gauss_weight*(DN_v(3,0)*crLHS79 + DN_v(3,1)*crLHS96 + crLHS151);
rLHS(7,12)+=gauss_weight*(DN_p(0,1)*crLHS153 + DN_p(0,1)*crLHS154 - crLHS157);
rLHS(7,13)+=gauss_weight*(DN_p(1,1)*crLHS153 + DN_p(1,1)*crLHS154 - crLHS158);
rLHS(7,14)+=gauss_weight*(DN_p(2,1)*crLHS153 + DN_p(2,1)*crLHS154 - crLHS159);
rLHS(8,0)+=gauss_weight*(DN_v(4,0)*crLHS0 + DN_v(4,1)*crLHS2 + crLHS163);
rLHS(8,1)+=gauss_weight*(DN_v(4,0)*crLHS18 + DN_v(4,1)*crLHS20);
rLHS(8,2)+=gauss_weight*(DN_v(4,0)*crLHS21 + DN_v(4,1)*crLHS23 + crLHS164);
rLHS(8,3)+=gauss_weight*(DN_v(4,0)*crLHS30 + DN_v(4,1)*crLHS32);
rLHS(8,4)+=gauss_weight*(DN_v(4,0)*crLHS33 + DN_v(4,1)*crLHS35 + crLHS165);
rLHS(8,5)+=gauss_weight*(DN_v(4,0)*crLHS42 + DN_v(4,1)*crLHS44);
rLHS(8,6)+=gauss_weight*(DN_v(4,0)*crLHS45 + DN_v(4,1)*crLHS47 + crLHS166);
rLHS(8,7)+=gauss_weight*(DN_v(4,0)*crLHS54 + DN_v(4,1)*crLHS56);
rLHS(8,8)+=gauss_weight*(DN_v(4,0)*crLHS57 + DN_v(4,1)*crLHS59 + crLHS167);
rLHS(8,9)+=gauss_weight*(DN_v(4,0)*crLHS66 + DN_v(4,1)*crLHS68);
rLHS(8,10)+=gauss_weight*(DN_v(4,0)*crLHS69 + DN_v(4,1)*crLHS71 + crLHS169);
rLHS(8,11)+=gauss_weight*(DN_v(4,0)*crLHS77 + DN_v(4,1)*crLHS79);
rLHS(8,12)+=gauss_weight*(DN_p(0,0)*crLHS171 + DN_p(0,0)*crLHS172 - crLHS170);
rLHS(8,13)+=gauss_weight*(DN_p(1,0)*crLHS171 + DN_p(1,0)*crLHS172 - crLHS173);
rLHS(8,14)+=gauss_weight*(DN_p(2,0)*crLHS171 + DN_p(2,0)*crLHS172 - crLHS174);
rLHS(9,0)+=gauss_weight*(DN_v(4,0)*crLHS2 + DN_v(4,1)*crLHS85);
rLHS(9,1)+=gauss_weight*(DN_v(4,0)*crLHS20 + DN_v(4,1)*crLHS86 + crLHS163);
rLHS(9,2)+=gauss_weight*(DN_v(4,0)*crLHS23 + DN_v(4,1)*crLHS87);
rLHS(9,3)+=gauss_weight*(DN_v(4,0)*crLHS32 + DN_v(4,1)*crLHS88 + crLHS164);
rLHS(9,4)+=gauss_weight*(DN_v(4,0)*crLHS35 + DN_v(4,1)*crLHS89);
rLHS(9,5)+=gauss_weight*(DN_v(4,0)*crLHS44 + DN_v(4,1)*crLHS90 + crLHS165);
rLHS(9,6)+=gauss_weight*(DN_v(4,0)*crLHS47 + DN_v(4,1)*crLHS91);
rLHS(9,7)+=gauss_weight*(DN_v(4,0)*crLHS56 + DN_v(4,1)*crLHS92 + crLHS166);
rLHS(9,8)+=gauss_weight*(DN_v(4,0)*crLHS59 + DN_v(4,1)*crLHS93);
rLHS(9,9)+=gauss_weight*(DN_v(4,0)*crLHS68 + DN_v(4,1)*crLHS94 + crLHS167);
rLHS(9,10)+=gauss_weight*(DN_v(4,0)*crLHS71 + DN_v(4,1)*crLHS95);
rLHS(9,11)+=gauss_weight*(DN_v(4,0)*crLHS79 + DN_v(4,1)*crLHS96 + crLHS169);
rLHS(9,12)+=gauss_weight*(DN_p(0,1)*crLHS171 + DN_p(0,1)*crLHS172 - crLHS175);
rLHS(9,13)+=gauss_weight*(DN_p(1,1)*crLHS171 + DN_p(1,1)*crLHS172 - crLHS176);
rLHS(9,14)+=gauss_weight*(DN_p(2,1)*crLHS171 + DN_p(2,1)*crLHS172 - crLHS177);
rLHS(10,0)+=gauss_weight*(DN_v(5,0)*crLHS0 + DN_v(5,1)*crLHS2 + crLHS181);
rLHS(10,1)+=gauss_weight*(DN_v(5,0)*crLHS18 + DN_v(5,1)*crLHS20);
rLHS(10,2)+=gauss_weight*(DN_v(5,0)*crLHS21 + DN_v(5,1)*crLHS23 + crLHS182);
rLHS(10,3)+=gauss_weight*(DN_v(5,0)*crLHS30 + DN_v(5,1)*crLHS32);
rLHS(10,4)+=gauss_weight*(DN_v(5,0)*crLHS33 + DN_v(5,1)*crLHS35 + crLHS183);
rLHS(10,5)+=gauss_weight*(DN_v(5,0)*crLHS42 + DN_v(5,1)*crLHS44);
rLHS(10,6)+=gauss_weight*(DN_v(5,0)*crLHS45 + DN_v(5,1)*crLHS47 + crLHS184);
rLHS(10,7)+=gauss_weight*(DN_v(5,0)*crLHS54 + DN_v(5,1)*crLHS56);
rLHS(10,8)+=gauss_weight*(DN_v(5,0)*crLHS57 + DN_v(5,1)*crLHS59 + crLHS185);
rLHS(10,9)+=gauss_weight*(DN_v(5,0)*crLHS66 + DN_v(5,1)*crLHS68);
rLHS(10,10)+=gauss_weight*(DN_v(5,0)*crLHS69 + DN_v(5,1)*crLHS71 + crLHS186);
rLHS(10,11)+=gauss_weight*(DN_v(5,0)*crLHS77 + DN_v(5,1)*crLHS79);
rLHS(10,12)+=gauss_weight*(DN_p(0,0)*crLHS188 + DN_p(0,0)*crLHS189 - crLHS187);
rLHS(10,13)+=gauss_weight*(DN_p(1,0)*crLHS188 + DN_p(1,0)*crLHS189 - crLHS190);
rLHS(10,14)+=gauss_weight*(DN_p(2,0)*crLHS188 + DN_p(2,0)*crLHS189 - crLHS191);
rLHS(11,0)+=gauss_weight*(DN_v(5,0)*crLHS2 + DN_v(5,1)*crLHS85);
rLHS(11,1)+=gauss_weight*(DN_v(5,0)*crLHS20 + DN_v(5,1)*crLHS86 + crLHS181);
rLHS(11,2)+=gauss_weight*(DN_v(5,0)*crLHS23 + DN_v(5,1)*crLHS87);
rLHS(11,3)+=gauss_weight*(DN_v(5,0)*crLHS32 + DN_v(5,1)*crLHS88 + crLHS182);
rLHS(11,4)+=gauss_weight*(DN_v(5,0)*crLHS35 + DN_v(5,1)*crLHS89);
rLHS(11,5)+=gauss_weight*(DN_v(5,0)*crLHS44 + DN_v(5,1)*crLHS90 + crLHS183);
rLHS(11,6)+=gauss_weight*(DN_v(5,0)*crLHS47 + DN_v(5,1)*crLHS91);
rLHS(11,7)+=gauss_weight*(DN_v(5,0)*crLHS56 + DN_v(5,1)*crLHS92 + crLHS184);
rLHS(11,8)+=gauss_weight*(DN_v(5,0)*crLHS59 + DN_v(5,1)*crLHS93);
rLHS(11,9)+=gauss_weight*(DN_v(5,0)*crLHS68 + DN_v(5,1)*crLHS94 + crLHS185);
rLHS(11,10)+=gauss_weight*(DN_v(5,0)*crLHS71 + DN_v(5,1)*crLHS95);
rLHS(11,11)+=gauss_weight*(DN_v(5,0)*crLHS79 + DN_v(5,1)*crLHS96 + crLHS186);
rLHS(11,12)+=gauss_weight*(DN_p(0,1)*crLHS188 + DN_p(0,1)*crLHS189 - crLHS192);
rLHS(11,13)+=gauss_weight*(DN_p(1,1)*crLHS188 + DN_p(1,1)*crLHS189 - crLHS193);
rLHS(11,14)+=gauss_weight*(DN_p(2,1)*crLHS188 + DN_p(2,1)*crLHS189 - crLHS194);
rLHS(12,0)+=gauss_weight*(DN_p(0,0)*crLHS196 + crLHS80);
rLHS(12,1)+=gauss_weight*(DN_p(0,1)*crLHS196 + crLHS97);
rLHS(12,2)+=gauss_weight*(DN_p(0,0)*crLHS197 + crLHS113);
rLHS(12,3)+=gauss_weight*(DN_p(0,1)*crLHS197 + crLHS118);
rLHS(12,4)+=gauss_weight*(DN_p(0,0)*crLHS198 + crLHS133);
rLHS(12,5)+=gauss_weight*(DN_p(0,1)*crLHS198 + crLHS138);
rLHS(12,6)+=gauss_weight*(DN_p(0,0)*crLHS199 + crLHS152);
rLHS(12,7)+=gauss_weight*(DN_p(0,1)*crLHS199 + crLHS157);
rLHS(12,8)+=gauss_weight*(DN_p(0,0)*crLHS200 + crLHS170);
rLHS(12,9)+=gauss_weight*(DN_p(0,1)*crLHS200 + crLHS175);
rLHS(12,10)+=gauss_weight*(DN_p(0,0)*crLHS201 + crLHS187);
rLHS(12,11)+=gauss_weight*(DN_p(0,1)*crLHS201 + crLHS192);
rLHS(12,12)+=crLHS202*(pow(DN_p(0,0), 2) + pow(DN_p(0,1), 2));
rLHS(12,13)+=crLHS203;
rLHS(12,14)+=crLHS204;
rLHS(13,0)+=gauss_weight*(DN_p(1,0)*crLHS196 + crLHS83);
rLHS(13,1)+=gauss_weight*(DN_p(1,1)*crLHS196 + crLHS98);
rLHS(13,2)+=gauss_weight*(DN_p(1,0)*crLHS197 + crLHS116);
rLHS(13,3)+=gauss_weight*(DN_p(1,1)*crLHS197 + crLHS119);
rLHS(13,4)+=gauss_weight*(DN_p(1,0)*crLHS198 + crLHS136);
rLHS(13,5)+=gauss_weight*(DN_p(1,1)*crLHS198 + crLHS139);
rLHS(13,6)+=gauss_weight*(DN_p(1,0)*crLHS199 + crLHS155);
rLHS(13,7)+=gauss_weight*(DN_p(1,1)*crLHS199 + crLHS158);
rLHS(13,8)+=gauss_weight*(DN_p(1,0)*crLHS200 + crLHS173);
rLHS(13,9)+=gauss_weight*(DN_p(1,1)*crLHS200 + crLHS176);
rLHS(13,10)+=gauss_weight*(DN_p(1,0)*crLHS201 + crLHS190);
rLHS(13,11)+=gauss_weight*(DN_p(1,1)*crLHS201 + crLHS193);
rLHS(13,12)+=crLHS203;
rLHS(13,13)+=crLHS202*(pow(DN_p(1,0), 2) + pow(DN_p(1,1), 2));
rLHS(13,14)+=crLHS205;
rLHS(14,0)+=gauss_weight*(DN_p(2,0)*crLHS196 + crLHS84);
rLHS(14,1)+=gauss_weight*(DN_p(2,1)*crLHS196 + crLHS99);
rLHS(14,2)+=gauss_weight*(DN_p(2,0)*crLHS197 + crLHS117);
rLHS(14,3)+=gauss_weight*(DN_p(2,1)*crLHS197 + crLHS120);
rLHS(14,4)+=gauss_weight*(DN_p(2,0)*crLHS198 + crLHS137);
rLHS(14,5)+=gauss_weight*(DN_p(2,1)*crLHS198 + crLHS140);
rLHS(14,6)+=gauss_weight*(DN_p(2,0)*crLHS199 + crLHS156);
rLHS(14,7)+=gauss_weight*(DN_p(2,1)*crLHS199 + crLHS159);
rLHS(14,8)+=gauss_weight*(DN_p(2,0)*crLHS200 + crLHS174);
rLHS(14,9)+=gauss_weight*(DN_p(2,1)*crLHS200 + crLHS177);
rLHS(14,10)+=gauss_weight*(DN_p(2,0)*crLHS201 + crLHS191);
rLHS(14,11)+=gauss_weight*(DN_p(2,1)*crLHS201 + crLHS194);
rLHS(14,12)+=crLHS204;
rLHS(14,13)+=crLHS205;
rLHS(14,14)+=crLHS202*(pow(DN_p(2,0), 2) + pow(DN_p(2,1), 2));

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

    // Get shape function values
    const auto& N_p = rData.N_p;
    const auto& N_v = rData.N_v;
    const auto& DN_p = rData.DN_p;
    const auto& DN_v = rData.DN_v;

    // Assemble RHS contribution
    const double gauss_weight = rData.Weight;
    const double crRHS0 = N_p[0]*r_p[0] + N_p[1]*r_p[1] + N_p[2]*r_p[2];
const double crRHS1 = rho*(N_v[0]*r_f(0,0) + N_v[1]*r_f(1,0) + N_v[2]*r_f(2,0) + N_v[3]*r_f(3,0) + N_v[4]*r_f(4,0) + N_v[5]*r_f(5,0));
const double crRHS2 = rho*(N_v[0]*(rData.BDF0*r_v(0,0) + rData.BDF1*r_vn(0,0) + rData.BDF2*r_vnn(0,0)) + N_v[1]*(rData.BDF0*r_v(1,0) + rData.BDF1*r_vn(1,0) + rData.BDF2*r_vnn(1,0)) + N_v[2]*(rData.BDF0*r_v(2,0) + rData.BDF1*r_vn(2,0) + rData.BDF2*r_vnn(2,0)) + N_v[3]*(rData.BDF0*r_v(3,0) + rData.BDF1*r_vn(3,0) + rData.BDF2*r_vnn(3,0)) + N_v[4]*(rData.BDF0*r_v(4,0) + rData.BDF1*r_vn(4,0) + rData.BDF2*r_vnn(4,0)) + N_v[5]*(rData.BDF0*r_v(5,0) + rData.BDF1*r_vn(5,0) + rData.BDF2*r_vnn(5,0)));
const double crRHS3 = DN_v(0,0)*r_v(0,0) + DN_v(1,0)*r_v(1,0) + DN_v(2,0)*r_v(2,0) + DN_v(3,0)*r_v(3,0) + DN_v(4,0)*r_v(4,0) + DN_v(5,0)*r_v(5,0);
const double crRHS4 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crRHS5 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crRHS6 = rho*(crRHS3*crRHS4 + crRHS5*(DN_v(0,1)*r_v(0,0) + DN_v(1,1)*r_v(1,0) + DN_v(2,1)*r_v(2,0) + DN_v(3,1)*r_v(3,0) + DN_v(4,1)*r_v(4,0) + DN_v(5,1)*r_v(5,0)));
const double crRHS7 = 1.0/(dyn_tau*rho/rData.DeltaTime + rho*stab_c2*sqrt(pow(crRHS4, 2) + pow(crRHS5, 2))/h + mu*stab_c1/pow(h, 2));
const double crRHS8 = crRHS7*(DN_p(0,0)*r_p[0] + DN_p(1,0)*r_p[1] + DN_p(2,0)*r_p[2] - crRHS1 + crRHS2 + crRHS6);
const double crRHS9 = rho*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1));
const double crRHS10 = N_v[0]*crRHS9;
const double crRHS11 = rho*(DN_v(0,0)*crRHS4 + DN_v(0,1)*crRHS5);
const double crRHS12 = rho*(N_v[0]*r_f(0,1) + N_v[1]*r_f(1,1) + N_v[2]*r_f(2,1) + N_v[3]*r_f(3,1) + N_v[4]*r_f(4,1) + N_v[5]*r_f(5,1));
const double crRHS13 = rho*(N_v[0]*(rData.BDF0*r_v(0,1) + rData.BDF1*r_vn(0,1) + rData.BDF2*r_vnn(0,1)) + N_v[1]*(rData.BDF0*r_v(1,1) + rData.BDF1*r_vn(1,1) + rData.BDF2*r_vnn(1,1)) + N_v[2]*(rData.BDF0*r_v(2,1) + rData.BDF1*r_vn(2,1) + rData.BDF2*r_vnn(2,1)) + N_v[3]*(rData.BDF0*r_v(3,1) + rData.BDF1*r_vn(3,1) + rData.BDF2*r_vnn(3,1)) + N_v[4]*(rData.BDF0*r_v(4,1) + rData.BDF1*r_vn(4,1) + rData.BDF2*r_vnn(4,1)) + N_v[5]*(rData.BDF0*r_v(5,1) + rData.BDF1*r_vn(5,1) + rData.BDF2*r_vnn(5,1)));
const double crRHS14 = DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1);
const double crRHS15 = rho*(crRHS14*crRHS5 + crRHS4*(DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1)));
const double crRHS16 = crRHS7*(DN_p(0,1)*r_p[0] + DN_p(1,1)*r_p[1] + DN_p(2,1)*r_p[2] - crRHS12 + crRHS13 + crRHS15);
const double crRHS17 = N_v[1]*crRHS9;
const double crRHS18 = rho*(DN_v(1,0)*crRHS4 + DN_v(1,1)*crRHS5);
const double crRHS19 = N_v[2]*crRHS9;
const double crRHS20 = rho*(DN_v(2,0)*crRHS4 + DN_v(2,1)*crRHS5);
const double crRHS21 = N_v[3]*crRHS9;
const double crRHS22 = rho*(DN_v(3,0)*crRHS4 + DN_v(3,1)*crRHS5);
const double crRHS23 = N_v[4]*crRHS9;
const double crRHS24 = rho*(DN_v(4,0)*crRHS4 + DN_v(4,1)*crRHS5);
const double crRHS25 = N_v[5]*crRHS9;
const double crRHS26 = rho*(DN_v(5,0)*crRHS4 + DN_v(5,1)*crRHS5);
const double crRHS27 = crRHS14 + crRHS3;
rRHS[0]+=-gauss_weight*(-DN_v(0,0)*crRHS0 + DN_v(0,0)*r_stress[0] + DN_v(0,1)*r_stress[2] - N_v[0]*crRHS1 + N_v[0]*crRHS2 + N_v[0]*crRHS6 + crRHS10*crRHS8 + crRHS11*crRHS8);
rRHS[1]+=-gauss_weight*(DN_v(0,0)*r_stress[2] - DN_v(0,1)*crRHS0 + DN_v(0,1)*r_stress[1] - N_v[0]*crRHS12 + N_v[0]*crRHS13 + N_v[0]*crRHS15 + crRHS10*crRHS16 + crRHS11*crRHS16);
rRHS[2]+=-gauss_weight*(-DN_v(1,0)*crRHS0 + DN_v(1,0)*r_stress[0] + DN_v(1,1)*r_stress[2] - N_v[1]*crRHS1 + N_v[1]*crRHS2 + N_v[1]*crRHS6 + crRHS17*crRHS8 + crRHS18*crRHS8);
rRHS[3]+=-gauss_weight*(DN_v(1,0)*r_stress[2] - DN_v(1,1)*crRHS0 + DN_v(1,1)*r_stress[1] - N_v[1]*crRHS12 + N_v[1]*crRHS13 + N_v[1]*crRHS15 + crRHS16*crRHS17 + crRHS16*crRHS18);
rRHS[4]+=-gauss_weight*(-DN_v(2,0)*crRHS0 + DN_v(2,0)*r_stress[0] + DN_v(2,1)*r_stress[2] - N_v[2]*crRHS1 + N_v[2]*crRHS2 + N_v[2]*crRHS6 + crRHS19*crRHS8 + crRHS20*crRHS8);
rRHS[5]+=-gauss_weight*(DN_v(2,0)*r_stress[2] - DN_v(2,1)*crRHS0 + DN_v(2,1)*r_stress[1] - N_v[2]*crRHS12 + N_v[2]*crRHS13 + N_v[2]*crRHS15 + crRHS16*crRHS19 + crRHS16*crRHS20);
rRHS[6]+=-gauss_weight*(-DN_v(3,0)*crRHS0 + DN_v(3,0)*r_stress[0] + DN_v(3,1)*r_stress[2] - N_v[3]*crRHS1 + N_v[3]*crRHS2 + N_v[3]*crRHS6 + crRHS21*crRHS8 + crRHS22*crRHS8);
rRHS[7]+=-gauss_weight*(DN_v(3,0)*r_stress[2] - DN_v(3,1)*crRHS0 + DN_v(3,1)*r_stress[1] - N_v[3]*crRHS12 + N_v[3]*crRHS13 + N_v[3]*crRHS15 + crRHS16*crRHS21 + crRHS16*crRHS22);
rRHS[8]+=-gauss_weight*(-DN_v(4,0)*crRHS0 + DN_v(4,0)*r_stress[0] + DN_v(4,1)*r_stress[2] - N_v[4]*crRHS1 + N_v[4]*crRHS2 + N_v[4]*crRHS6 + crRHS23*crRHS8 + crRHS24*crRHS8);
rRHS[9]+=-gauss_weight*(DN_v(4,0)*r_stress[2] - DN_v(4,1)*crRHS0 + DN_v(4,1)*r_stress[1] - N_v[4]*crRHS12 + N_v[4]*crRHS13 + N_v[4]*crRHS15 + crRHS16*crRHS23 + crRHS16*crRHS24);
rRHS[10]+=-gauss_weight*(-DN_v(5,0)*crRHS0 + DN_v(5,0)*r_stress[0] + DN_v(5,1)*r_stress[2] - N_v[5]*crRHS1 + N_v[5]*crRHS2 + N_v[5]*crRHS6 + crRHS25*crRHS8 + crRHS26*crRHS8);
rRHS[11]+=-gauss_weight*(DN_v(5,0)*r_stress[2] - DN_v(5,1)*crRHS0 + DN_v(5,1)*r_stress[1] - N_v[5]*crRHS12 + N_v[5]*crRHS13 + N_v[5]*crRHS15 + crRHS16*crRHS25 + crRHS16*crRHS26);
rRHS[12]+=-gauss_weight*(DN_p(0,0)*crRHS8 + DN_p(0,1)*crRHS16 + N_p[0]*crRHS27);
rRHS[13]+=-gauss_weight*(DN_p(1,0)*crRHS8 + DN_p(1,1)*crRHS16 + N_p[1]*crRHS27);
rRHS[14]+=-gauss_weight*(DN_p(2,0)*crRHS8 + DN_p(2,1)*crRHS16 + N_p[2]*crRHS27);

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

    // Get shape function values
    const auto& N_p = rData.N_p;
    const auto& N_v = rData.N_v;
    const auto& DN_p = rData.DN_p;
    const auto& DN_v = rData.DN_v;

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

    // Assemble enrichment contribution
    const double gauss_weight = rData.Weight;
    const BoundedMatrix<double,1,2> v_e = ZeroMatrix(1,2);

    const double crRHSee0 = N_p[0]*r_p[0] + N_p[1]*r_p[1] + N_p[2]*r_p[2];
const double crRHSee1 = N_e*rho;
const double crRHSee2 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crRHSee3 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crRHSee4 = crRHSee1*(DN_e(0,0)*crRHSee2 + DN_e(0,1)*crRHSee3);
rRHSee[0]+=-gauss_weight*(-DN_e(0,0)*crRHSee0 + DN_e(0,0)*r_stress[0] + DN_e(0,1)*r_stress[2] + crRHSee1*(crRHSee2*(DN_v(0,0)*r_v(0,0) + DN_v(1,0)*r_v(1,0) + DN_v(2,0)*r_v(2,0) + DN_v(3,0)*r_v(3,0) + DN_v(4,0)*r_v(4,0) + DN_v(5,0)*r_v(5,0)) + crRHSee3*(DN_v(0,1)*r_v(0,0) + DN_v(1,1)*r_v(1,0) + DN_v(2,1)*r_v(2,0) + DN_v(3,1)*r_v(3,0) + DN_v(4,1)*r_v(4,0) + DN_v(5,1)*r_v(5,0))) - crRHSee1*(N_v[0]*r_f(0,0) + N_v[1]*r_f(1,0) + N_v[2]*r_f(2,0) + N_v[3]*r_f(3,0) + N_v[4]*r_f(4,0) + N_v[5]*r_f(5,0)) + crRHSee1*(N_v[0]*(rData.BDF0*r_v(0,0) + rData.BDF1*r_vn(0,0) + rData.BDF2*r_vnn(0,0)) + N_v[1]*(rData.BDF0*r_v(1,0) + rData.BDF1*r_vn(1,0) + rData.BDF2*r_vnn(1,0)) + N_v[2]*(rData.BDF0*r_v(2,0) + rData.BDF1*r_vn(2,0) + rData.BDF2*r_vnn(2,0)) + N_v[3]*(rData.BDF0*r_v(3,0) + rData.BDF1*r_vn(3,0) + rData.BDF2*r_vnn(3,0)) + N_v[4]*(rData.BDF0*r_v(4,0) + rData.BDF1*r_vn(4,0) + rData.BDF2*r_vnn(4,0)) + N_v[5]*(rData.BDF0*r_v(5,0) + rData.BDF1*r_vn(5,0) + rData.BDF2*r_vnn(5,0))) - crRHSee4*v_e(0,0));
rRHSee[1]+=-gauss_weight*(DN_e(0,0)*r_stress[2] - DN_e(0,1)*crRHSee0 + DN_e(0,1)*r_stress[1] + crRHSee1*(crRHSee2*(DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1)) + crRHSee3*(DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1))) - crRHSee1*(N_v[0]*r_f(0,1) + N_v[1]*r_f(1,1) + N_v[2]*r_f(2,1) + N_v[3]*r_f(3,1) + N_v[4]*r_f(4,1) + N_v[5]*r_f(5,1)) + crRHSee1*(N_v[0]*(rData.BDF0*r_v(0,1) + rData.BDF1*r_vn(0,1) + rData.BDF2*r_vnn(0,1)) + N_v[1]*(rData.BDF0*r_v(1,1) + rData.BDF1*r_vn(1,1) + rData.BDF2*r_vnn(1,1)) + N_v[2]*(rData.BDF0*r_v(2,1) + rData.BDF1*r_vn(2,1) + rData.BDF2*r_vnn(2,1)) + N_v[3]*(rData.BDF0*r_v(3,1) + rData.BDF1*r_vn(3,1) + rData.BDF2*r_vnn(3,1)) + N_v[4]*(rData.BDF0*r_v(4,1) + rData.BDF1*r_vn(4,1) + rData.BDF2*r_vnn(4,1)) + N_v[5]*(rData.BDF0*r_v(5,1) + rData.BDF1*r_vn(5,1) + rData.BDF2*r_vnn(5,1))) - crRHSee4*v_e(0,1));


    const double crKue0 = N_p[0]*gauss_weight;
const double crKue1 = N_p[1]*gauss_weight;
const double crKue2 = N_p[2]*gauss_weight;
rKue(0,0)+=0;
rKue(0,1)+=0;
rKue(1,0)+=0;
rKue(1,1)+=0;
rKue(2,0)+=0;
rKue(2,1)+=0;
rKue(3,0)+=0;
rKue(3,1)+=0;
rKue(4,0)+=0;
rKue(4,1)+=0;
rKue(5,0)+=0;
rKue(5,1)+=0;
rKue(6,0)+=0;
rKue(6,1)+=0;
rKue(7,0)+=0;
rKue(7,1)+=0;
rKue(8,0)+=0;
rKue(8,1)+=0;
rKue(9,0)+=0;
rKue(9,1)+=0;
rKue(10,0)+=0;
rKue(10,1)+=0;
rKue(11,0)+=0;
rKue(11,1)+=0;
rKue(12,0)+=DN_e(0,0)*crKue0;
rKue(12,1)+=DN_e(0,1)*crKue0;
rKue(13,0)+=DN_e(0,0)*crKue1;
rKue(13,1)+=DN_e(0,1)*crKue1;
rKue(14,0)+=DN_e(0,0)*crKue2;
rKue(14,1)+=DN_e(0,1)*crKue2;


    const double crKeu0 = C(0,2)*DN_v(0,0);
const double crKeu1 = C(2,2)*DN_v(0,1) + crKeu0;
const double crKeu2 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crKeu3 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crKeu4 = N_e*rho;
const double crKeu5 = crKeu4*rData.BDF0;
const double crKeu6 = N_v[0]*crKeu5 + crKeu4*(DN_v(0,0)*crKeu2 + DN_v(0,1)*crKeu3);
const double crKeu7 = C(1,2)*DN_v(0,1);
const double crKeu8 = C(2,2)*DN_v(0,0) + crKeu7;
const double crKeu9 = C(0,2)*DN_v(1,0);
const double crKeu10 = C(2,2)*DN_v(1,1) + crKeu9;
const double crKeu11 = N_v[1]*crKeu5 + crKeu4*(DN_v(1,0)*crKeu2 + DN_v(1,1)*crKeu3);
const double crKeu12 = C(1,2)*DN_v(1,1);
const double crKeu13 = C(2,2)*DN_v(1,0) + crKeu12;
const double crKeu14 = C(0,2)*DN_v(2,0);
const double crKeu15 = C(2,2)*DN_v(2,1) + crKeu14;
const double crKeu16 = N_v[2]*crKeu5 + crKeu4*(DN_v(2,0)*crKeu2 + DN_v(2,1)*crKeu3);
const double crKeu17 = C(1,2)*DN_v(2,1);
const double crKeu18 = C(2,2)*DN_v(2,0) + crKeu17;
const double crKeu19 = C(0,2)*DN_v(3,0);
const double crKeu20 = C(2,2)*DN_v(3,1) + crKeu19;
const double crKeu21 = N_v[3]*crKeu5 + crKeu4*(DN_v(3,0)*crKeu2 + DN_v(3,1)*crKeu3);
const double crKeu22 = C(1,2)*DN_v(3,1);
const double crKeu23 = C(2,2)*DN_v(3,0) + crKeu22;
const double crKeu24 = C(0,2)*DN_v(4,0);
const double crKeu25 = C(2,2)*DN_v(4,1) + crKeu24;
const double crKeu26 = N_v[4]*crKeu5 + crKeu4*(DN_v(4,0)*crKeu2 + DN_v(4,1)*crKeu3);
const double crKeu27 = C(1,2)*DN_v(4,1);
const double crKeu28 = C(2,2)*DN_v(4,0) + crKeu27;
const double crKeu29 = C(0,2)*DN_v(5,0);
const double crKeu30 = C(2,2)*DN_v(5,1) + crKeu29;
const double crKeu31 = N_v[5]*crKeu5 + crKeu4*(DN_v(5,0)*crKeu2 + DN_v(5,1)*crKeu3);
const double crKeu32 = C(1,2)*DN_v(5,1);
const double crKeu33 = C(2,2)*DN_v(5,0) + crKeu32;
const double crKeu34 = DN_e(0,0)*gauss_weight;
const double crKeu35 = DN_e(0,1)*gauss_weight;
rKeu(0,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(0,0) + C(0,2)*DN_v(0,1)) + DN_e(0,1)*crKeu1 + crKeu6);
rKeu(0,1)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(0,1) + crKeu0) + DN_e(0,1)*crKeu8);
rKeu(0,2)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(1,0) + C(0,2)*DN_v(1,1)) + DN_e(0,1)*crKeu10 + crKeu11);
rKeu(0,3)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(1,1) + crKeu9) + DN_e(0,1)*crKeu13);
rKeu(0,4)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(2,0) + C(0,2)*DN_v(2,1)) + DN_e(0,1)*crKeu15 + crKeu16);
rKeu(0,5)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(2,1) + crKeu14) + DN_e(0,1)*crKeu18);
rKeu(0,6)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(3,0) + C(0,2)*DN_v(3,1)) + DN_e(0,1)*crKeu20 + crKeu21);
rKeu(0,7)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(3,1) + crKeu19) + DN_e(0,1)*crKeu23);
rKeu(0,8)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(4,0) + C(0,2)*DN_v(4,1)) + DN_e(0,1)*crKeu25 + crKeu26);
rKeu(0,9)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(4,1) + crKeu24) + DN_e(0,1)*crKeu28);
rKeu(0,10)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_v(5,0) + C(0,2)*DN_v(5,1)) + DN_e(0,1)*crKeu30 + crKeu31);
rKeu(0,11)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_v(5,1) + crKeu29) + DN_e(0,1)*crKeu33);
rKeu(0,12)+=-N_p[0]*crKeu34;
rKeu(0,13)+=-N_p[1]*crKeu34;
rKeu(0,14)+=-N_p[2]*crKeu34;
rKeu(1,0)+=gauss_weight*(DN_e(0,0)*crKeu1 + DN_e(0,1)*(C(0,1)*DN_v(0,0) + crKeu7));
rKeu(1,1)+=gauss_weight*(DN_e(0,0)*crKeu8 + DN_e(0,1)*(C(1,1)*DN_v(0,1) + C(1,2)*DN_v(0,0)) + crKeu6);
rKeu(1,2)+=gauss_weight*(DN_e(0,0)*crKeu10 + DN_e(0,1)*(C(0,1)*DN_v(1,0) + crKeu12));
rKeu(1,3)+=gauss_weight*(DN_e(0,0)*crKeu13 + DN_e(0,1)*(C(1,1)*DN_v(1,1) + C(1,2)*DN_v(1,0)) + crKeu11);
rKeu(1,4)+=gauss_weight*(DN_e(0,0)*crKeu15 + DN_e(0,1)*(C(0,1)*DN_v(2,0) + crKeu17));
rKeu(1,5)+=gauss_weight*(DN_e(0,0)*crKeu18 + DN_e(0,1)*(C(1,1)*DN_v(2,1) + C(1,2)*DN_v(2,0)) + crKeu16);
rKeu(1,6)+=gauss_weight*(DN_e(0,0)*crKeu20 + DN_e(0,1)*(C(0,1)*DN_v(3,0) + crKeu22));
rKeu(1,7)+=gauss_weight*(DN_e(0,0)*crKeu23 + DN_e(0,1)*(C(1,1)*DN_v(3,1) + C(1,2)*DN_v(3,0)) + crKeu21);
rKeu(1,8)+=gauss_weight*(DN_e(0,0)*crKeu25 + DN_e(0,1)*(C(0,1)*DN_v(4,0) + crKeu27));
rKeu(1,9)+=gauss_weight*(DN_e(0,0)*crKeu28 + DN_e(0,1)*(C(1,1)*DN_v(4,1) + C(1,2)*DN_v(4,0)) + crKeu26);
rKeu(1,10)+=gauss_weight*(DN_e(0,0)*crKeu30 + DN_e(0,1)*(C(0,1)*DN_v(5,0) + crKeu32));
rKeu(1,11)+=gauss_weight*(DN_e(0,0)*crKeu33 + DN_e(0,1)*(C(1,1)*DN_v(5,1) + C(1,2)*DN_v(5,0)) + crKeu31);
rKeu(1,12)+=-N_p[0]*crKeu35;
rKeu(1,13)+=-N_p[1]*crKeu35;
rKeu(1,14)+=-N_p[2]*crKeu35;


    const double crKee0 = C(0,2)*DN_e(0,0);
const double crKee1 = C(2,2)*DN_e(0,1) + crKee0;
const double crKee2 = -N_e*rho*(DN_e(0,0)*(N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0)) + DN_e(0,1)*(N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1)));
const double crKee3 = C(1,2)*DN_e(0,1);
const double crKee4 = C(2,2)*DN_e(0,0) + crKee3;
rKee(0,0)+=gauss_weight*(DN_e(0,0)*(C(0,0)*DN_e(0,0) + C(0,2)*DN_e(0,1)) + DN_e(0,1)*crKee1 + crKee2);
rKee(0,1)+=gauss_weight*(DN_e(0,0)*(C(0,1)*DN_e(0,1) + crKee0) + DN_e(0,1)*crKee4);
rKee(1,0)+=gauss_weight*(DN_e(0,0)*crKee1 + DN_e(0,1)*(C(0,1)*DN_e(0,0) + crKee3));
rKee(1,1)+=gauss_weight*(DN_e(0,0)*crKee4 + DN_e(0,1)*(C(1,1)*DN_e(0,1) + C(1,2)*DN_e(0,0)) + crKee2);

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