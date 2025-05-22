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
#include "incompressible_navier_stokes_p2_p1_continuous.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
IncompressibleNavierStokesP2P1Continuous<TDim>::IncompressibleNavierStokesP2P1Continuous(IndexType NewId)
    : Element(NewId)
{}

template< unsigned int TDim >
IncompressibleNavierStokesP2P1Continuous<TDim>::IncompressibleNavierStokesP2P1Continuous(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{}

template< unsigned int TDim >
IncompressibleNavierStokesP2P1Continuous<TDim>::IncompressibleNavierStokesP2P1Continuous(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{}

template< unsigned int TDim >
IncompressibleNavierStokesP2P1Continuous<TDim>::IncompressibleNavierStokesP2P1Continuous(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{}

template< unsigned int TDim >
IncompressibleNavierStokesP2P1Continuous<TDim>::~IncompressibleNavierStokesP2P1Continuous()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer IncompressibleNavierStokesP2P1Continuous<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<IncompressibleNavierStokesP2P1Continuous<TDim>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< unsigned int TDim >
Element::Pointer IncompressibleNavierStokesP2P1Continuous<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<IncompressibleNavierStokesP2P1Continuous<TDim>>(NewId, pGeom, pProperties);
}

template< unsigned int TDim >
void IncompressibleNavierStokesP2P1Continuous<TDim>::Initialize(const ProcessInfo& rCurrentProcessInfo)
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
void IncompressibleNavierStokesP2P1Continuous<TDim>::EquationIdVector(
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
void IncompressibleNavierStokesP2P1Continuous<TDim>::GetDofList(
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
void IncompressibleNavierStokesP2P1Continuous<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and initialize output
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
    DenseVector<GeometryType::ShapeFunctionsSecondDerivativesType> velocity_DDN;
    CalculateKinematics(weights, velocity_N, pressure_N, velocity_DN, pressure_DN, velocity_DDN);

    // Loop Gauss points
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(IntegrationMethod);
    for (IndexType g = 0; g < n_gauss; ++g) {
        // Set current Gauss point kinematics
        noalias(aux_data.N_v) = row(velocity_N, g);
        noalias(aux_data.N_p) = row(pressure_N, g);
        noalias(aux_data.DN_v) = velocity_DN[g];
        noalias(aux_data.DN_p) = pressure_DN[g];
        aux_data.DDN_v = velocity_DDN[g];
        aux_data.Weight = weights[g];

        // Calculate current Gauss point material response
        CalculateStrainRate(aux_data);
        mpConstitutiveLaw->CalculateMaterialResponseCauchy(cons_law_params);
        mpConstitutiveLaw->CalculateValue(cons_law_params, EFFECTIVE_VISCOSITY, aux_data.EffectiveViscosity);

        // Assemble standard Galerkin contribution
        AddGaussPointLeftHandSideContribution(aux_data, rLeftHandSideMatrix);
        AddGaussPointRightHandSideContribution(aux_data, rRightHandSideVector);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template< unsigned int TDim >
int IncompressibleNavierStokesP2P1Continuous<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
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
const Parameters IncompressibleNavierStokesP2P1Continuous<TDim>::GetSpecifications() const
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
std::string IncompressibleNavierStokesP2P1Continuous<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "IncompressibleNavierStokesP2P1Continuous" << TDim << "D" << VelocityNumNodes << "N #" << this->Id();
    return buffer.str();
}

template< unsigned int TDim >
void IncompressibleNavierStokesP2P1Continuous<TDim>::PrintInfo(std::ostream& rOStream) const
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
void IncompressibleNavierStokesP2P1Continuous<TDim>::SetElementData(
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
    rData.StabC1 = 12.0;
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
void IncompressibleNavierStokesP2P1Continuous<TDim>::CalculateKinematics(
    Vector& rGaussWeights,
    Matrix& rVelocityN,
    Matrix& rPressureN,
    GeometryType::ShapeFunctionsGradientsType& rVelocityDNDX,
    GeometryType::ShapeFunctionsGradientsType& rPressureDNDX,
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

    // Calculate integration points weight
    if (rGaussWeights.size() != n_gauss) {
        rGaussWeights.resize(n_gauss, false);
    }
    for (IndexType g = 0; g < n_gauss; ++g) {
        rGaussWeights[g] = det_J_vect[g] * integration_points[g].Weight();
    }
}

template< unsigned int TDim >
void IncompressibleNavierStokesP2P1Continuous<TDim>::CalculateStrainRate(ElementDataContainer& rData)
{
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
void IncompressibleNavierStokesP2P1Continuous<2>::AddGaussPointLeftHandSideContribution(
    const ElementDataContainer& rData,
    MatrixType& rLHS)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

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
    const double crLHS3 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
    const double crLHS4 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
    const double crLHS5 = rho*stab_c2*std::sqrt(std::pow(crLHS3, 2) + std::pow(crLHS4, 2));
    const double crLHS6 = crLHS5*h/stab_c1 + mu;
    const double crLHS7 = 1.0*C(0,0);
    const double crLHS8 = C(0,2)*DDN_v[0](0,0);
    const double crLHS9 = 1.0*DDN_v[0](0,1);
    const double crLHS10 = rho*(DN_v(0,0)*crLHS3 + DN_v(0,1)*crLHS4);
    const double crLHS11 = rData.BDF0*rho;
    const double crLHS12 = N_v[0]*crLHS11;
    const double crLHS13 = -crLHS10 - crLHS12;
    const double crLHS14 = C(0,2)*crLHS9 + C(2,2)*crLHS9 + DDN_v[0](0,0)*crLHS7 + crLHS13 + 1.0*crLHS8;
    const double crLHS15 = 1.0/(crLHS5/h + dyn_tau*rho/rData.DeltaTime + mu*stab_c1/std::pow(h, 2));
    const double crLHS16 = crLHS10*crLHS15;
    const double crLHS17 = rho*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1));
    const double crLHS18 = N_v[0]*crLHS17;
    const double crLHS19 = crLHS15*crLHS18;
    const double crLHS20 = std::pow(N_v[0], 2)*crLHS11 + N_v[0]*crLHS10;
    const double crLHS21 = C(0,1)*DN_v(0,1) + crLHS1;
    const double crLHS22 = C(1,2)*DN_v(0,1);
    const double crLHS23 = C(2,2)*DN_v(0,0) + crLHS22;
    const double crLHS24 = DN_v(0,0)*crLHS6;
    const double crLHS25 = DN_v(0,1)*crLHS24;
    const double crLHS26 = C(0,1)*DDN_v[0](0,1) + C(1,2)*DDN_v[0](0,1) + C(2,2)*DDN_v[0](0,0) + crLHS8;
    const double crLHS27 = C(0,0)*DN_v(1,0) + C(0,2)*DN_v(1,1);
    const double crLHS28 = C(0,2)*DN_v(1,0);
    const double crLHS29 = C(2,2)*DN_v(1,1) + crLHS28;
    const double crLHS30 = DN_v(1,0)*crLHS24;
    const double crLHS31 = C(0,2)*DDN_v[1](0,0);
    const double crLHS32 = 1.0*DDN_v[1](0,1);
    const double crLHS33 = rho*(DN_v(1,0)*crLHS3 + DN_v(1,1)*crLHS4);
    const double crLHS34 = N_v[1]*crLHS11;
    const double crLHS35 = -crLHS33 - crLHS34;
    const double crLHS36 = C(0,2)*crLHS32 + C(2,2)*crLHS32 + DDN_v[1](0,0)*crLHS7 + 1.0*crLHS31 + crLHS35;
    const double crLHS37 = N_v[1]*crLHS12;
    const double crLHS38 = N_v[0]*crLHS33 + crLHS37;
    const double crLHS39 = C(0,1)*DN_v(1,1) + crLHS28;
    const double crLHS40 = C(1,2)*DN_v(1,1);
    const double crLHS41 = C(2,2)*DN_v(1,0) + crLHS40;
    const double crLHS42 = DN_v(1,1)*crLHS24;
    const double crLHS43 = C(0,1)*DDN_v[1](0,1) + C(1,2)*DDN_v[1](0,1) + C(2,2)*DDN_v[1](0,0) + crLHS31;
    const double crLHS44 = C(0,0)*DN_v(2,0) + C(0,2)*DN_v(2,1);
    const double crLHS45 = C(0,2)*DN_v(2,0);
    const double crLHS46 = C(2,2)*DN_v(2,1) + crLHS45;
    const double crLHS47 = DN_v(2,0)*crLHS24;
    const double crLHS48 = C(0,2)*DDN_v[2](0,0);
    const double crLHS49 = 1.0*DDN_v[2](0,1);
    const double crLHS50 = rho*(DN_v(2,0)*crLHS3 + DN_v(2,1)*crLHS4);
    const double crLHS51 = N_v[2]*crLHS11;
    const double crLHS52 = -crLHS50 - crLHS51;
    const double crLHS53 = C(0,2)*crLHS49 + C(2,2)*crLHS49 + DDN_v[2](0,0)*crLHS7 + 1.0*crLHS48 + crLHS52;
    const double crLHS54 = N_v[2]*crLHS12;
    const double crLHS55 = N_v[0]*crLHS50 + crLHS54;
    const double crLHS56 = C(0,1)*DN_v(2,1) + crLHS45;
    const double crLHS57 = C(1,2)*DN_v(2,1);
    const double crLHS58 = C(2,2)*DN_v(2,0) + crLHS57;
    const double crLHS59 = DN_v(2,1)*crLHS24;
    const double crLHS60 = C(0,1)*DDN_v[2](0,1) + C(1,2)*DDN_v[2](0,1) + C(2,2)*DDN_v[2](0,0) + crLHS48;
    const double crLHS61 = C(0,0)*DN_v(3,0) + C(0,2)*DN_v(3,1);
    const double crLHS62 = C(0,2)*DN_v(3,0);
    const double crLHS63 = C(2,2)*DN_v(3,1) + crLHS62;
    const double crLHS64 = DN_v(3,0)*crLHS24;
    const double crLHS65 = C(0,2)*DDN_v[3](0,0);
    const double crLHS66 = 1.0*DDN_v[3](0,1);
    const double crLHS67 = rho*(DN_v(3,0)*crLHS3 + DN_v(3,1)*crLHS4);
    const double crLHS68 = N_v[3]*crLHS11;
    const double crLHS69 = -crLHS67 - crLHS68;
    const double crLHS70 = C(0,2)*crLHS66 + C(2,2)*crLHS66 + DDN_v[3](0,0)*crLHS7 + 1.0*crLHS65 + crLHS69;
    const double crLHS71 = N_v[3]*crLHS12;
    const double crLHS72 = N_v[0]*crLHS67 + crLHS71;
    const double crLHS73 = C(0,1)*DN_v(3,1) + crLHS62;
    const double crLHS74 = C(1,2)*DN_v(3,1);
    const double crLHS75 = C(2,2)*DN_v(3,0) + crLHS74;
    const double crLHS76 = DN_v(3,1)*crLHS24;
    const double crLHS77 = C(0,1)*DDN_v[3](0,1) + C(1,2)*DDN_v[3](0,1) + C(2,2)*DDN_v[3](0,0) + crLHS65;
    const double crLHS78 = C(0,0)*DN_v(4,0) + C(0,2)*DN_v(4,1);
    const double crLHS79 = C(0,2)*DN_v(4,0);
    const double crLHS80 = C(2,2)*DN_v(4,1) + crLHS79;
    const double crLHS81 = DN_v(4,0)*crLHS24;
    const double crLHS82 = C(0,2)*DDN_v[4](0,0);
    const double crLHS83 = 1.0*DDN_v[4](0,1);
    const double crLHS84 = rho*(DN_v(4,0)*crLHS3 + DN_v(4,1)*crLHS4);
    const double crLHS85 = N_v[4]*crLHS11;
    const double crLHS86 = -crLHS84 - crLHS85;
    const double crLHS87 = C(0,2)*crLHS83 + C(2,2)*crLHS83 + DDN_v[4](0,0)*crLHS7 + 1.0*crLHS82 + crLHS86;
    const double crLHS88 = N_v[4]*crLHS12;
    const double crLHS89 = N_v[0]*crLHS84 + crLHS88;
    const double crLHS90 = C(0,1)*DN_v(4,1) + crLHS79;
    const double crLHS91 = C(1,2)*DN_v(4,1);
    const double crLHS92 = C(2,2)*DN_v(4,0) + crLHS91;
    const double crLHS93 = DN_v(4,1)*crLHS24;
    const double crLHS94 = C(0,1)*DDN_v[4](0,1) + C(1,2)*DDN_v[4](0,1) + C(2,2)*DDN_v[4](0,0) + crLHS82;
    const double crLHS95 = C(0,0)*DN_v(5,0) + C(0,2)*DN_v(5,1);
    const double crLHS96 = C(0,2)*DN_v(5,0);
    const double crLHS97 = C(2,2)*DN_v(5,1) + crLHS96;
    const double crLHS98 = DN_v(5,0)*crLHS24;
    const double crLHS99 = C(0,2)*DDN_v[5](0,0);
    const double crLHS100 = 1.0*DDN_v[5](0,1);
    const double crLHS101 = rho*(DN_v(5,0)*crLHS3 + DN_v(5,1)*crLHS4);
    const double crLHS102 = -N_v[5]*crLHS11 - crLHS101;
    const double crLHS103 = C(0,2)*crLHS100 + C(2,2)*crLHS100 + DDN_v[5](0,0)*crLHS7 + crLHS102 + 1.0*crLHS99;
    const double crLHS104 = N_v[5]*crLHS12;
    const double crLHS105 = N_v[0]*crLHS101 + crLHS104;
    const double crLHS106 = C(0,1)*DN_v(5,1) + crLHS96;
    const double crLHS107 = C(1,2)*DN_v(5,1);
    const double crLHS108 = C(2,2)*DN_v(5,0) + crLHS107;
    const double crLHS109 = DN_v(5,1)*crLHS24;
    const double crLHS110 = C(0,1)*DDN_v[5](0,1) + C(1,2)*DDN_v[5](0,1) + C(2,2)*DDN_v[5](0,0) + crLHS99;
    const double crLHS111 = -DN_v(0,0)*N_p[0];
    const double crLHS112 = DN_p(0,0)*crLHS15;
    const double crLHS113 = -DN_v(0,0)*N_p[1];
    const double crLHS114 = DN_p(1,0)*crLHS15;
    const double crLHS115 = -DN_v(0,0)*N_p[2];
    const double crLHS116 = DN_p(2,0)*crLHS15;
    const double crLHS117 = C(0,1)*DN_v(0,0) + crLHS22;
    const double crLHS118 = C(1,2)*DDN_v[0](1,1);
    const double crLHS119 = C(0,1)*DDN_v[0](1,0) + C(0,2)*DDN_v[0](1,0) + C(2,2)*DDN_v[0](1,1) + crLHS118;
    const double crLHS120 = C(1,1)*DN_v(0,1) + C(1,2)*DN_v(0,0);
    const double crLHS121 = 1.0*C(1,1);
    const double crLHS122 = 1.0*DDN_v[0](1,0);
    const double crLHS123 = C(1,2)*crLHS122 + C(2,2)*crLHS122 + DDN_v[0](1,1)*crLHS121 + 1.0*crLHS118 + crLHS13;
    const double crLHS124 = C(0,1)*DN_v(1,0) + crLHS40;
    const double crLHS125 = DN_v(0,1)*crLHS6;
    const double crLHS126 = DN_v(1,0)*crLHS125;
    const double crLHS127 = C(1,2)*DDN_v[1](1,1);
    const double crLHS128 = C(0,1)*DDN_v[1](1,0) + C(0,2)*DDN_v[1](1,0) + C(2,2)*DDN_v[1](1,1) + crLHS127;
    const double crLHS129 = C(1,1)*DN_v(1,1) + C(1,2)*DN_v(1,0);
    const double crLHS130 = DN_v(1,1)*crLHS125;
    const double crLHS131 = 1.0*DDN_v[1](1,0);
    const double crLHS132 = C(1,2)*crLHS131 + C(2,2)*crLHS131 + DDN_v[1](1,1)*crLHS121 + 1.0*crLHS127 + crLHS35;
    const double crLHS133 = C(0,1)*DN_v(2,0) + crLHS57;
    const double crLHS134 = DN_v(2,0)*crLHS125;
    const double crLHS135 = C(1,2)*DDN_v[2](1,1);
    const double crLHS136 = C(0,1)*DDN_v[2](1,0) + C(0,2)*DDN_v[2](1,0) + C(2,2)*DDN_v[2](1,1) + crLHS135;
    const double crLHS137 = C(1,1)*DN_v(2,1) + C(1,2)*DN_v(2,0);
    const double crLHS138 = DN_v(2,1)*crLHS125;
    const double crLHS139 = 1.0*DDN_v[2](1,0);
    const double crLHS140 = C(1,2)*crLHS139 + C(2,2)*crLHS139 + DDN_v[2](1,1)*crLHS121 + 1.0*crLHS135 + crLHS52;
    const double crLHS141 = C(0,1)*DN_v(3,0) + crLHS74;
    const double crLHS142 = DN_v(3,0)*crLHS125;
    const double crLHS143 = C(1,2)*DDN_v[3](1,1);
    const double crLHS144 = C(0,1)*DDN_v[3](1,0) + C(0,2)*DDN_v[3](1,0) + C(2,2)*DDN_v[3](1,1) + crLHS143;
    const double crLHS145 = C(1,1)*DN_v(3,1) + C(1,2)*DN_v(3,0);
    const double crLHS146 = DN_v(3,1)*crLHS125;
    const double crLHS147 = 1.0*DDN_v[3](1,0);
    const double crLHS148 = C(1,2)*crLHS147 + C(2,2)*crLHS147 + DDN_v[3](1,1)*crLHS121 + 1.0*crLHS143 + crLHS69;
    const double crLHS149 = C(0,1)*DN_v(4,0) + crLHS91;
    const double crLHS150 = DN_v(4,0)*crLHS125;
    const double crLHS151 = C(1,2)*DDN_v[4](1,1);
    const double crLHS152 = C(0,1)*DDN_v[4](1,0) + C(0,2)*DDN_v[4](1,0) + C(2,2)*DDN_v[4](1,1) + crLHS151;
    const double crLHS153 = C(1,1)*DN_v(4,1) + C(1,2)*DN_v(4,0);
    const double crLHS154 = DN_v(4,1)*crLHS125;
    const double crLHS155 = 1.0*DDN_v[4](1,0);
    const double crLHS156 = C(1,2)*crLHS155 + C(2,2)*crLHS155 + DDN_v[4](1,1)*crLHS121 + 1.0*crLHS151 + crLHS86;
    const double crLHS157 = C(0,1)*DN_v(5,0) + crLHS107;
    const double crLHS158 = DN_v(5,0)*crLHS125;
    const double crLHS159 = C(1,2)*DDN_v[5](1,1);
    const double crLHS160 = C(0,1)*DDN_v[5](1,0) + C(0,2)*DDN_v[5](1,0) + C(2,2)*DDN_v[5](1,1) + crLHS159;
    const double crLHS161 = C(1,1)*DN_v(5,1) + C(1,2)*DN_v(5,0);
    const double crLHS162 = DN_v(5,1)*crLHS125;
    const double crLHS163 = 1.0*DDN_v[5](1,0);
    const double crLHS164 = C(1,2)*crLHS163 + C(2,2)*crLHS163 + DDN_v[5](1,1)*crLHS121 + crLHS102 + 1.0*crLHS159;
    const double crLHS165 = -DN_v(0,1)*N_p[0];
    const double crLHS166 = DN_p(0,1)*crLHS15;
    const double crLHS167 = -DN_v(0,1)*N_p[1];
    const double crLHS168 = DN_p(1,1)*crLHS15;
    const double crLHS169 = -DN_v(0,1)*N_p[2];
    const double crLHS170 = DN_p(2,1)*crLHS15;
    const double crLHS171 = crLHS15*crLHS33;
    const double crLHS172 = N_v[1]*crLHS17;
    const double crLHS173 = crLHS15*crLHS172;
    const double crLHS174 = N_v[1]*crLHS10 + crLHS37;
    const double crLHS175 = std::pow(N_v[1], 2)*crLHS11 + N_v[1]*crLHS33;
    const double crLHS176 = DN_v(1,0)*crLHS6;
    const double crLHS177 = DN_v(1,1)*crLHS176;
    const double crLHS178 = DN_v(2,0)*crLHS176;
    const double crLHS179 = N_v[2]*crLHS34;
    const double crLHS180 = N_v[1]*crLHS50 + crLHS179;
    const double crLHS181 = DN_v(2,1)*crLHS176;
    const double crLHS182 = DN_v(3,0)*crLHS176;
    const double crLHS183 = N_v[3]*crLHS34;
    const double crLHS184 = N_v[1]*crLHS67 + crLHS183;
    const double crLHS185 = DN_v(3,1)*crLHS176;
    const double crLHS186 = DN_v(4,0)*crLHS176;
    const double crLHS187 = N_v[4]*crLHS34;
    const double crLHS188 = N_v[1]*crLHS84 + crLHS187;
    const double crLHS189 = DN_v(4,1)*crLHS176;
    const double crLHS190 = DN_v(5,0)*crLHS176;
    const double crLHS191 = N_v[5]*crLHS34;
    const double crLHS192 = N_v[1]*crLHS101 + crLHS191;
    const double crLHS193 = DN_v(5,1)*crLHS176;
    const double crLHS194 = -DN_v(1,0)*N_p[0];
    const double crLHS195 = -DN_v(1,0)*N_p[1];
    const double crLHS196 = -DN_v(1,0)*N_p[2];
    const double crLHS197 = DN_v(1,1)*crLHS6;
    const double crLHS198 = DN_v(2,0)*crLHS197;
    const double crLHS199 = DN_v(2,1)*crLHS197;
    const double crLHS200 = DN_v(3,0)*crLHS197;
    const double crLHS201 = DN_v(3,1)*crLHS197;
    const double crLHS202 = DN_v(4,0)*crLHS197;
    const double crLHS203 = DN_v(4,1)*crLHS197;
    const double crLHS204 = DN_v(5,0)*crLHS197;
    const double crLHS205 = DN_v(5,1)*crLHS197;
    const double crLHS206 = -DN_v(1,1)*N_p[0];
    const double crLHS207 = -DN_v(1,1)*N_p[1];
    const double crLHS208 = -DN_v(1,1)*N_p[2];
    const double crLHS209 = crLHS15*crLHS50;
    const double crLHS210 = N_v[2]*crLHS17;
    const double crLHS211 = crLHS15*crLHS210;
    const double crLHS212 = N_v[2]*crLHS10 + crLHS54;
    const double crLHS213 = N_v[2]*crLHS33 + crLHS179;
    const double crLHS214 = std::pow(N_v[2], 2)*crLHS11 + N_v[2]*crLHS50;
    const double crLHS215 = DN_v(2,0)*crLHS6;
    const double crLHS216 = DN_v(2,1)*crLHS215;
    const double crLHS217 = DN_v(3,0)*crLHS215;
    const double crLHS218 = N_v[3]*crLHS51;
    const double crLHS219 = N_v[2]*crLHS67 + crLHS218;
    const double crLHS220 = DN_v(3,1)*crLHS215;
    const double crLHS221 = DN_v(4,0)*crLHS215;
    const double crLHS222 = N_v[4]*crLHS51;
    const double crLHS223 = N_v[2]*crLHS84 + crLHS222;
    const double crLHS224 = DN_v(4,1)*crLHS215;
    const double crLHS225 = DN_v(5,0)*crLHS215;
    const double crLHS226 = N_v[5]*crLHS51;
    const double crLHS227 = N_v[2]*crLHS101 + crLHS226;
    const double crLHS228 = DN_v(5,1)*crLHS215;
    const double crLHS229 = -DN_v(2,0)*N_p[0];
    const double crLHS230 = -DN_v(2,0)*N_p[1];
    const double crLHS231 = -DN_v(2,0)*N_p[2];
    const double crLHS232 = DN_v(2,1)*crLHS6;
    const double crLHS233 = DN_v(3,0)*crLHS232;
    const double crLHS234 = DN_v(3,1)*crLHS232;
    const double crLHS235 = DN_v(4,0)*crLHS232;
    const double crLHS236 = DN_v(4,1)*crLHS232;
    const double crLHS237 = DN_v(5,0)*crLHS232;
    const double crLHS238 = DN_v(5,1)*crLHS232;
    const double crLHS239 = -DN_v(2,1)*N_p[0];
    const double crLHS240 = -DN_v(2,1)*N_p[1];
    const double crLHS241 = -DN_v(2,1)*N_p[2];
    const double crLHS242 = crLHS15*crLHS67;
    const double crLHS243 = N_v[3]*crLHS17;
    const double crLHS244 = crLHS15*crLHS243;
    const double crLHS245 = N_v[3]*crLHS10 + crLHS71;
    const double crLHS246 = N_v[3]*crLHS33 + crLHS183;
    const double crLHS247 = N_v[3]*crLHS50 + crLHS218;
    const double crLHS248 = std::pow(N_v[3], 2)*crLHS11 + N_v[3]*crLHS67;
    const double crLHS249 = DN_v(3,0)*crLHS6;
    const double crLHS250 = DN_v(3,1)*crLHS249;
    const double crLHS251 = DN_v(4,0)*crLHS249;
    const double crLHS252 = N_v[4]*crLHS68;
    const double crLHS253 = N_v[3]*crLHS84 + crLHS252;
    const double crLHS254 = DN_v(4,1)*crLHS249;
    const double crLHS255 = DN_v(5,0)*crLHS249;
    const double crLHS256 = N_v[5]*crLHS68;
    const double crLHS257 = N_v[3]*crLHS101 + crLHS256;
    const double crLHS258 = DN_v(5,1)*crLHS249;
    const double crLHS259 = -DN_v(3,0)*N_p[0];
    const double crLHS260 = -DN_v(3,0)*N_p[1];
    const double crLHS261 = -DN_v(3,0)*N_p[2];
    const double crLHS262 = DN_v(3,1)*crLHS6;
    const double crLHS263 = DN_v(4,0)*crLHS262;
    const double crLHS264 = DN_v(4,1)*crLHS262;
    const double crLHS265 = DN_v(5,0)*crLHS262;
    const double crLHS266 = DN_v(5,1)*crLHS262;
    const double crLHS267 = -DN_v(3,1)*N_p[0];
    const double crLHS268 = -DN_v(3,1)*N_p[1];
    const double crLHS269 = -DN_v(3,1)*N_p[2];
    const double crLHS270 = crLHS15*crLHS84;
    const double crLHS271 = N_v[4]*crLHS17;
    const double crLHS272 = crLHS15*crLHS271;
    const double crLHS273 = N_v[4]*crLHS10 + crLHS88;
    const double crLHS274 = N_v[4]*crLHS33 + crLHS187;
    const double crLHS275 = N_v[4]*crLHS50 + crLHS222;
    const double crLHS276 = N_v[4]*crLHS67 + crLHS252;
    const double crLHS277 = std::pow(N_v[4], 2)*crLHS11 + N_v[4]*crLHS84;
    const double crLHS278 = DN_v(4,0)*crLHS6;
    const double crLHS279 = DN_v(4,1)*crLHS278;
    const double crLHS280 = DN_v(5,0)*crLHS278;
    const double crLHS281 = N_v[5]*crLHS85;
    const double crLHS282 = N_v[4]*crLHS101 + crLHS281;
    const double crLHS283 = DN_v(5,1)*crLHS278;
    const double crLHS284 = -DN_v(4,0)*N_p[0];
    const double crLHS285 = -DN_v(4,0)*N_p[1];
    const double crLHS286 = -DN_v(4,0)*N_p[2];
    const double crLHS287 = DN_v(4,1)*crLHS6;
    const double crLHS288 = DN_v(5,0)*crLHS287;
    const double crLHS289 = DN_v(5,1)*crLHS287;
    const double crLHS290 = -DN_v(4,1)*N_p[0];
    const double crLHS291 = -DN_v(4,1)*N_p[1];
    const double crLHS292 = -DN_v(4,1)*N_p[2];
    const double crLHS293 = crLHS101*crLHS15;
    const double crLHS294 = N_v[5]*crLHS17;
    const double crLHS295 = crLHS15*crLHS294;
    const double crLHS296 = N_v[5]*crLHS10 + crLHS104;
    const double crLHS297 = N_v[5]*crLHS33 + crLHS191;
    const double crLHS298 = N_v[5]*crLHS50 + crLHS226;
    const double crLHS299 = N_v[5]*crLHS67 + crLHS256;
    const double crLHS300 = N_v[5]*crLHS84 + crLHS281;
    const double crLHS301 = std::pow(N_v[5], 2)*crLHS11 + N_v[5]*crLHS101;
    const double crLHS302 = DN_v(5,0)*DN_v(5,1)*crLHS6;
    const double crLHS303 = -DN_v(5,0)*N_p[0];
    const double crLHS304 = -DN_v(5,0)*N_p[1];
    const double crLHS305 = -DN_v(5,0)*N_p[2];
    const double crLHS306 = -DN_v(5,1)*N_p[0];
    const double crLHS307 = -DN_v(5,1)*N_p[1];
    const double crLHS308 = -DN_v(5,1)*N_p[2];
    const double crLHS309 = crLHS15*gauss_weight;
    const double crLHS310 = crLHS309*(DN_p(0,0)*DN_p(1,0) + DN_p(0,1)*DN_p(1,1));
    const double crLHS311 = crLHS309*(DN_p(0,0)*DN_p(2,0) + DN_p(0,1)*DN_p(2,1));
    const double crLHS312 = crLHS309*(DN_p(1,0)*DN_p(2,0) + DN_p(1,1)*DN_p(2,1));
    rLHS(0,0)+=gauss_weight*(std::pow(DN_v(0,0), 2)*crLHS6 + DN_v(0,0)*crLHS0 + DN_v(0,1)*crLHS2 - crLHS14*crLHS16 - crLHS14*crLHS19 + crLHS20);
    rLHS(0,1)+=gauss_weight*(DN_v(0,0)*crLHS21 + DN_v(0,1)*crLHS23 - crLHS16*crLHS26 - crLHS19*crLHS26 + crLHS25);
    rLHS(0,2)+=gauss_weight*(DN_v(0,0)*crLHS27 + DN_v(0,1)*crLHS29 - crLHS16*crLHS36 - crLHS19*crLHS36 + crLHS30 + crLHS38);
    rLHS(0,3)+=gauss_weight*(DN_v(0,0)*crLHS39 + DN_v(0,1)*crLHS41 - crLHS16*crLHS43 - crLHS19*crLHS43 + crLHS42);
    rLHS(0,4)+=gauss_weight*(DN_v(0,0)*crLHS44 + DN_v(0,1)*crLHS46 - crLHS16*crLHS53 - crLHS19*crLHS53 + crLHS47 + crLHS55);
    rLHS(0,5)+=gauss_weight*(DN_v(0,0)*crLHS56 + DN_v(0,1)*crLHS58 - crLHS16*crLHS60 - crLHS19*crLHS60 + crLHS59);
    rLHS(0,6)+=gauss_weight*(DN_v(0,0)*crLHS61 + DN_v(0,1)*crLHS63 - crLHS16*crLHS70 - crLHS19*crLHS70 + crLHS64 + crLHS72);
    rLHS(0,7)+=gauss_weight*(DN_v(0,0)*crLHS73 + DN_v(0,1)*crLHS75 - crLHS16*crLHS77 - crLHS19*crLHS77 + crLHS76);
    rLHS(0,8)+=gauss_weight*(DN_v(0,0)*crLHS78 + DN_v(0,1)*crLHS80 - crLHS16*crLHS87 - crLHS19*crLHS87 + crLHS81 + crLHS89);
    rLHS(0,9)+=gauss_weight*(DN_v(0,0)*crLHS90 + DN_v(0,1)*crLHS92 - crLHS16*crLHS94 - crLHS19*crLHS94 + crLHS93);
    rLHS(0,10)+=gauss_weight*(DN_v(0,0)*crLHS95 + DN_v(0,1)*crLHS97 - crLHS103*crLHS16 - crLHS103*crLHS19 + crLHS105 + crLHS98);
    rLHS(0,11)+=gauss_weight*(DN_v(0,0)*crLHS106 + DN_v(0,1)*crLHS108 + crLHS109 - crLHS110*crLHS16 - crLHS110*crLHS19);
    rLHS(0,12)+=gauss_weight*(crLHS10*crLHS112 + crLHS111 + crLHS112*crLHS18);
    rLHS(0,13)+=gauss_weight*(crLHS10*crLHS114 + crLHS113 + crLHS114*crLHS18);
    rLHS(0,14)+=gauss_weight*(crLHS10*crLHS116 + crLHS115 + crLHS116*crLHS18);
    rLHS(1,0)+=gauss_weight*(DN_v(0,0)*crLHS2 + DN_v(0,1)*crLHS117 - crLHS119*crLHS16 - crLHS119*crLHS19 + crLHS25);
    rLHS(1,1)+=gauss_weight*(DN_v(0,0)*crLHS23 + std::pow(DN_v(0,1), 2)*crLHS6 + DN_v(0,1)*crLHS120 - crLHS123*crLHS16 - crLHS123*crLHS19 + crLHS20);
    rLHS(1,2)+=gauss_weight*(DN_v(0,0)*crLHS29 + DN_v(0,1)*crLHS124 + crLHS126 - crLHS128*crLHS16 - crLHS128*crLHS19);
    rLHS(1,3)+=gauss_weight*(DN_v(0,0)*crLHS41 + DN_v(0,1)*crLHS129 + crLHS130 - crLHS132*crLHS16 - crLHS132*crLHS19 + crLHS38);
    rLHS(1,4)+=gauss_weight*(DN_v(0,0)*crLHS46 + DN_v(0,1)*crLHS133 + crLHS134 - crLHS136*crLHS16 - crLHS136*crLHS19);
    rLHS(1,5)+=gauss_weight*(DN_v(0,0)*crLHS58 + DN_v(0,1)*crLHS137 + crLHS138 - crLHS140*crLHS16 - crLHS140*crLHS19 + crLHS55);
    rLHS(1,6)+=gauss_weight*(DN_v(0,0)*crLHS63 + DN_v(0,1)*crLHS141 + crLHS142 - crLHS144*crLHS16 - crLHS144*crLHS19);
    rLHS(1,7)+=gauss_weight*(DN_v(0,0)*crLHS75 + DN_v(0,1)*crLHS145 + crLHS146 - crLHS148*crLHS16 - crLHS148*crLHS19 + crLHS72);
    rLHS(1,8)+=gauss_weight*(DN_v(0,0)*crLHS80 + DN_v(0,1)*crLHS149 + crLHS150 - crLHS152*crLHS16 - crLHS152*crLHS19);
    rLHS(1,9)+=gauss_weight*(DN_v(0,0)*crLHS92 + DN_v(0,1)*crLHS153 + crLHS154 - crLHS156*crLHS16 - crLHS156*crLHS19 + crLHS89);
    rLHS(1,10)+=gauss_weight*(DN_v(0,0)*crLHS97 + DN_v(0,1)*crLHS157 + crLHS158 - crLHS16*crLHS160 - crLHS160*crLHS19);
    rLHS(1,11)+=gauss_weight*(DN_v(0,0)*crLHS108 + DN_v(0,1)*crLHS161 + crLHS105 - crLHS16*crLHS164 + crLHS162 - crLHS164*crLHS19);
    rLHS(1,12)+=gauss_weight*(crLHS10*crLHS166 + crLHS165 + crLHS166*crLHS18);
    rLHS(1,13)+=gauss_weight*(crLHS10*crLHS168 + crLHS167 + crLHS168*crLHS18);
    rLHS(1,14)+=gauss_weight*(crLHS10*crLHS170 + crLHS169 + crLHS170*crLHS18);
    rLHS(2,0)+=gauss_weight*(DN_v(1,0)*crLHS0 + DN_v(1,1)*crLHS2 - crLHS14*crLHS171 - crLHS14*crLHS173 + crLHS174 + crLHS30);
    rLHS(2,1)+=gauss_weight*(DN_v(1,0)*crLHS21 + DN_v(1,1)*crLHS23 + crLHS126 - crLHS171*crLHS26 - crLHS173*crLHS26);
    rLHS(2,2)+=gauss_weight*(std::pow(DN_v(1,0), 2)*crLHS6 + DN_v(1,0)*crLHS27 + DN_v(1,1)*crLHS29 - crLHS171*crLHS36 - crLHS173*crLHS36 + crLHS175);
    rLHS(2,3)+=gauss_weight*(DN_v(1,0)*crLHS39 + DN_v(1,1)*crLHS41 - crLHS171*crLHS43 - crLHS173*crLHS43 + crLHS177);
    rLHS(2,4)+=gauss_weight*(DN_v(1,0)*crLHS44 + DN_v(1,1)*crLHS46 - crLHS171*crLHS53 - crLHS173*crLHS53 + crLHS178 + crLHS180);
    rLHS(2,5)+=gauss_weight*(DN_v(1,0)*crLHS56 + DN_v(1,1)*crLHS58 - crLHS171*crLHS60 - crLHS173*crLHS60 + crLHS181);
    rLHS(2,6)+=gauss_weight*(DN_v(1,0)*crLHS61 + DN_v(1,1)*crLHS63 - crLHS171*crLHS70 - crLHS173*crLHS70 + crLHS182 + crLHS184);
    rLHS(2,7)+=gauss_weight*(DN_v(1,0)*crLHS73 + DN_v(1,1)*crLHS75 - crLHS171*crLHS77 - crLHS173*crLHS77 + crLHS185);
    rLHS(2,8)+=gauss_weight*(DN_v(1,0)*crLHS78 + DN_v(1,1)*crLHS80 - crLHS171*crLHS87 - crLHS173*crLHS87 + crLHS186 + crLHS188);
    rLHS(2,9)+=gauss_weight*(DN_v(1,0)*crLHS90 + DN_v(1,1)*crLHS92 - crLHS171*crLHS94 - crLHS173*crLHS94 + crLHS189);
    rLHS(2,10)+=gauss_weight*(DN_v(1,0)*crLHS95 + DN_v(1,1)*crLHS97 - crLHS103*crLHS171 - crLHS103*crLHS173 + crLHS190 + crLHS192);
    rLHS(2,11)+=gauss_weight*(DN_v(1,0)*crLHS106 + DN_v(1,1)*crLHS108 - crLHS110*crLHS171 - crLHS110*crLHS173 + crLHS193);
    rLHS(2,12)+=gauss_weight*(crLHS112*crLHS172 + crLHS112*crLHS33 + crLHS194);
    rLHS(2,13)+=gauss_weight*(crLHS114*crLHS172 + crLHS114*crLHS33 + crLHS195);
    rLHS(2,14)+=gauss_weight*(crLHS116*crLHS172 + crLHS116*crLHS33 + crLHS196);
    rLHS(3,0)+=gauss_weight*(DN_v(1,0)*crLHS2 + DN_v(1,1)*crLHS117 - crLHS119*crLHS171 - crLHS119*crLHS173 + crLHS42);
    rLHS(3,1)+=gauss_weight*(DN_v(1,0)*crLHS23 + DN_v(1,1)*crLHS120 - crLHS123*crLHS171 - crLHS123*crLHS173 + crLHS130 + crLHS174);
    rLHS(3,2)+=gauss_weight*(DN_v(1,0)*crLHS29 + DN_v(1,1)*crLHS124 - crLHS128*crLHS171 - crLHS128*crLHS173 + crLHS177);
    rLHS(3,3)+=gauss_weight*(DN_v(1,0)*crLHS41 + std::pow(DN_v(1,1), 2)*crLHS6 + DN_v(1,1)*crLHS129 - crLHS132*crLHS171 - crLHS132*crLHS173 + crLHS175);
    rLHS(3,4)+=gauss_weight*(DN_v(1,0)*crLHS46 + DN_v(1,1)*crLHS133 - crLHS136*crLHS171 - crLHS136*crLHS173 + crLHS198);
    rLHS(3,5)+=gauss_weight*(DN_v(1,0)*crLHS58 + DN_v(1,1)*crLHS137 - crLHS140*crLHS171 - crLHS140*crLHS173 + crLHS180 + crLHS199);
    rLHS(3,6)+=gauss_weight*(DN_v(1,0)*crLHS63 + DN_v(1,1)*crLHS141 - crLHS144*crLHS171 - crLHS144*crLHS173 + crLHS200);
    rLHS(3,7)+=gauss_weight*(DN_v(1,0)*crLHS75 + DN_v(1,1)*crLHS145 - crLHS148*crLHS171 - crLHS148*crLHS173 + crLHS184 + crLHS201);
    rLHS(3,8)+=gauss_weight*(DN_v(1,0)*crLHS80 + DN_v(1,1)*crLHS149 - crLHS152*crLHS171 - crLHS152*crLHS173 + crLHS202);
    rLHS(3,9)+=gauss_weight*(DN_v(1,0)*crLHS92 + DN_v(1,1)*crLHS153 - crLHS156*crLHS171 - crLHS156*crLHS173 + crLHS188 + crLHS203);
    rLHS(3,10)+=gauss_weight*(DN_v(1,0)*crLHS97 + DN_v(1,1)*crLHS157 - crLHS160*crLHS171 - crLHS160*crLHS173 + crLHS204);
    rLHS(3,11)+=gauss_weight*(DN_v(1,0)*crLHS108 + DN_v(1,1)*crLHS161 - crLHS164*crLHS171 - crLHS164*crLHS173 + crLHS192 + crLHS205);
    rLHS(3,12)+=gauss_weight*(crLHS166*crLHS172 + crLHS166*crLHS33 + crLHS206);
    rLHS(3,13)+=gauss_weight*(crLHS168*crLHS172 + crLHS168*crLHS33 + crLHS207);
    rLHS(3,14)+=gauss_weight*(crLHS170*crLHS172 + crLHS170*crLHS33 + crLHS208);
    rLHS(4,0)+=gauss_weight*(DN_v(2,0)*crLHS0 + DN_v(2,1)*crLHS2 - crLHS14*crLHS209 - crLHS14*crLHS211 + crLHS212 + crLHS47);
    rLHS(4,1)+=gauss_weight*(DN_v(2,0)*crLHS21 + DN_v(2,1)*crLHS23 + crLHS134 - crLHS209*crLHS26 - crLHS211*crLHS26);
    rLHS(4,2)+=gauss_weight*(DN_v(2,0)*crLHS27 + DN_v(2,1)*crLHS29 + crLHS178 - crLHS209*crLHS36 - crLHS211*crLHS36 + crLHS213);
    rLHS(4,3)+=gauss_weight*(DN_v(2,0)*crLHS39 + DN_v(2,1)*crLHS41 + crLHS198 - crLHS209*crLHS43 - crLHS211*crLHS43);
    rLHS(4,4)+=gauss_weight*(std::pow(DN_v(2,0), 2)*crLHS6 + DN_v(2,0)*crLHS44 + DN_v(2,1)*crLHS46 - crLHS209*crLHS53 - crLHS211*crLHS53 + crLHS214);
    rLHS(4,5)+=gauss_weight*(DN_v(2,0)*crLHS56 + DN_v(2,1)*crLHS58 - crLHS209*crLHS60 - crLHS211*crLHS60 + crLHS216);
    rLHS(4,6)+=gauss_weight*(DN_v(2,0)*crLHS61 + DN_v(2,1)*crLHS63 - crLHS209*crLHS70 - crLHS211*crLHS70 + crLHS217 + crLHS219);
    rLHS(4,7)+=gauss_weight*(DN_v(2,0)*crLHS73 + DN_v(2,1)*crLHS75 - crLHS209*crLHS77 - crLHS211*crLHS77 + crLHS220);
    rLHS(4,8)+=gauss_weight*(DN_v(2,0)*crLHS78 + DN_v(2,1)*crLHS80 - crLHS209*crLHS87 - crLHS211*crLHS87 + crLHS221 + crLHS223);
    rLHS(4,9)+=gauss_weight*(DN_v(2,0)*crLHS90 + DN_v(2,1)*crLHS92 - crLHS209*crLHS94 - crLHS211*crLHS94 + crLHS224);
    rLHS(4,10)+=gauss_weight*(DN_v(2,0)*crLHS95 + DN_v(2,1)*crLHS97 - crLHS103*crLHS209 - crLHS103*crLHS211 + crLHS225 + crLHS227);
    rLHS(4,11)+=gauss_weight*(DN_v(2,0)*crLHS106 + DN_v(2,1)*crLHS108 - crLHS110*crLHS209 - crLHS110*crLHS211 + crLHS228);
    rLHS(4,12)+=gauss_weight*(crLHS112*crLHS210 + crLHS112*crLHS50 + crLHS229);
    rLHS(4,13)+=gauss_weight*(crLHS114*crLHS210 + crLHS114*crLHS50 + crLHS230);
    rLHS(4,14)+=gauss_weight*(crLHS116*crLHS210 + crLHS116*crLHS50 + crLHS231);
    rLHS(5,0)+=gauss_weight*(DN_v(2,0)*crLHS2 + DN_v(2,1)*crLHS117 - crLHS119*crLHS209 - crLHS119*crLHS211 + crLHS59);
    rLHS(5,1)+=gauss_weight*(DN_v(2,0)*crLHS23 + DN_v(2,1)*crLHS120 - crLHS123*crLHS209 - crLHS123*crLHS211 + crLHS138 + crLHS212);
    rLHS(5,2)+=gauss_weight*(DN_v(2,0)*crLHS29 + DN_v(2,1)*crLHS124 - crLHS128*crLHS209 - crLHS128*crLHS211 + crLHS181);
    rLHS(5,3)+=gauss_weight*(DN_v(2,0)*crLHS41 + DN_v(2,1)*crLHS129 - crLHS132*crLHS209 - crLHS132*crLHS211 + crLHS199 + crLHS213);
    rLHS(5,4)+=gauss_weight*(DN_v(2,0)*crLHS46 + DN_v(2,1)*crLHS133 - crLHS136*crLHS209 - crLHS136*crLHS211 + crLHS216);
    rLHS(5,5)+=gauss_weight*(DN_v(2,0)*crLHS58 + std::pow(DN_v(2,1), 2)*crLHS6 + DN_v(2,1)*crLHS137 - crLHS140*crLHS209 - crLHS140*crLHS211 + crLHS214);
    rLHS(5,6)+=gauss_weight*(DN_v(2,0)*crLHS63 + DN_v(2,1)*crLHS141 - crLHS144*crLHS209 - crLHS144*crLHS211 + crLHS233);
    rLHS(5,7)+=gauss_weight*(DN_v(2,0)*crLHS75 + DN_v(2,1)*crLHS145 - crLHS148*crLHS209 - crLHS148*crLHS211 + crLHS219 + crLHS234);
    rLHS(5,8)+=gauss_weight*(DN_v(2,0)*crLHS80 + DN_v(2,1)*crLHS149 - crLHS152*crLHS209 - crLHS152*crLHS211 + crLHS235);
    rLHS(5,9)+=gauss_weight*(DN_v(2,0)*crLHS92 + DN_v(2,1)*crLHS153 - crLHS156*crLHS209 - crLHS156*crLHS211 + crLHS223 + crLHS236);
    rLHS(5,10)+=gauss_weight*(DN_v(2,0)*crLHS97 + DN_v(2,1)*crLHS157 - crLHS160*crLHS209 - crLHS160*crLHS211 + crLHS237);
    rLHS(5,11)+=gauss_weight*(DN_v(2,0)*crLHS108 + DN_v(2,1)*crLHS161 - crLHS164*crLHS209 - crLHS164*crLHS211 + crLHS227 + crLHS238);
    rLHS(5,12)+=gauss_weight*(crLHS166*crLHS210 + crLHS166*crLHS50 + crLHS239);
    rLHS(5,13)+=gauss_weight*(crLHS168*crLHS210 + crLHS168*crLHS50 + crLHS240);
    rLHS(5,14)+=gauss_weight*(crLHS170*crLHS210 + crLHS170*crLHS50 + crLHS241);
    rLHS(6,0)+=gauss_weight*(DN_v(3,0)*crLHS0 + DN_v(3,1)*crLHS2 - crLHS14*crLHS242 - crLHS14*crLHS244 + crLHS245 + crLHS64);
    rLHS(6,1)+=gauss_weight*(DN_v(3,0)*crLHS21 + DN_v(3,1)*crLHS23 + crLHS142 - crLHS242*crLHS26 - crLHS244*crLHS26);
    rLHS(6,2)+=gauss_weight*(DN_v(3,0)*crLHS27 + DN_v(3,1)*crLHS29 + crLHS182 - crLHS242*crLHS36 - crLHS244*crLHS36 + crLHS246);
    rLHS(6,3)+=gauss_weight*(DN_v(3,0)*crLHS39 + DN_v(3,1)*crLHS41 + crLHS200 - crLHS242*crLHS43 - crLHS244*crLHS43);
    rLHS(6,4)+=gauss_weight*(DN_v(3,0)*crLHS44 + DN_v(3,1)*crLHS46 + crLHS217 - crLHS242*crLHS53 - crLHS244*crLHS53 + crLHS247);
    rLHS(6,5)+=gauss_weight*(DN_v(3,0)*crLHS56 + DN_v(3,1)*crLHS58 + crLHS233 - crLHS242*crLHS60 - crLHS244*crLHS60);
    rLHS(6,6)+=gauss_weight*(std::pow(DN_v(3,0), 2)*crLHS6 + DN_v(3,0)*crLHS61 + DN_v(3,1)*crLHS63 - crLHS242*crLHS70 - crLHS244*crLHS70 + crLHS248);
    rLHS(6,7)+=gauss_weight*(DN_v(3,0)*crLHS73 + DN_v(3,1)*crLHS75 - crLHS242*crLHS77 - crLHS244*crLHS77 + crLHS250);
    rLHS(6,8)+=gauss_weight*(DN_v(3,0)*crLHS78 + DN_v(3,1)*crLHS80 - crLHS242*crLHS87 - crLHS244*crLHS87 + crLHS251 + crLHS253);
    rLHS(6,9)+=gauss_weight*(DN_v(3,0)*crLHS90 + DN_v(3,1)*crLHS92 - crLHS242*crLHS94 - crLHS244*crLHS94 + crLHS254);
    rLHS(6,10)+=gauss_weight*(DN_v(3,0)*crLHS95 + DN_v(3,1)*crLHS97 - crLHS103*crLHS242 - crLHS103*crLHS244 + crLHS255 + crLHS257);
    rLHS(6,11)+=gauss_weight*(DN_v(3,0)*crLHS106 + DN_v(3,1)*crLHS108 - crLHS110*crLHS242 - crLHS110*crLHS244 + crLHS258);
    rLHS(6,12)+=gauss_weight*(crLHS112*crLHS243 + crLHS112*crLHS67 + crLHS259);
    rLHS(6,13)+=gauss_weight*(crLHS114*crLHS243 + crLHS114*crLHS67 + crLHS260);
    rLHS(6,14)+=gauss_weight*(crLHS116*crLHS243 + crLHS116*crLHS67 + crLHS261);
    rLHS(7,0)+=gauss_weight*(DN_v(3,0)*crLHS2 + DN_v(3,1)*crLHS117 - crLHS119*crLHS242 - crLHS119*crLHS244 + crLHS76);
    rLHS(7,1)+=gauss_weight*(DN_v(3,0)*crLHS23 + DN_v(3,1)*crLHS120 - crLHS123*crLHS242 - crLHS123*crLHS244 + crLHS146 + crLHS245);
    rLHS(7,2)+=gauss_weight*(DN_v(3,0)*crLHS29 + DN_v(3,1)*crLHS124 - crLHS128*crLHS242 - crLHS128*crLHS244 + crLHS185);
    rLHS(7,3)+=gauss_weight*(DN_v(3,0)*crLHS41 + DN_v(3,1)*crLHS129 - crLHS132*crLHS242 - crLHS132*crLHS244 + crLHS201 + crLHS246);
    rLHS(7,4)+=gauss_weight*(DN_v(3,0)*crLHS46 + DN_v(3,1)*crLHS133 - crLHS136*crLHS242 - crLHS136*crLHS244 + crLHS220);
    rLHS(7,5)+=gauss_weight*(DN_v(3,0)*crLHS58 + DN_v(3,1)*crLHS137 - crLHS140*crLHS242 - crLHS140*crLHS244 + crLHS234 + crLHS247);
    rLHS(7,6)+=gauss_weight*(DN_v(3,0)*crLHS63 + DN_v(3,1)*crLHS141 - crLHS144*crLHS242 - crLHS144*crLHS244 + crLHS250);
    rLHS(7,7)+=gauss_weight*(DN_v(3,0)*crLHS75 + std::pow(DN_v(3,1), 2)*crLHS6 + DN_v(3,1)*crLHS145 - crLHS148*crLHS242 - crLHS148*crLHS244 + crLHS248);
    rLHS(7,8)+=gauss_weight*(DN_v(3,0)*crLHS80 + DN_v(3,1)*crLHS149 - crLHS152*crLHS242 - crLHS152*crLHS244 + crLHS263);
    rLHS(7,9)+=gauss_weight*(DN_v(3,0)*crLHS92 + DN_v(3,1)*crLHS153 - crLHS156*crLHS242 - crLHS156*crLHS244 + crLHS253 + crLHS264);
    rLHS(7,10)+=gauss_weight*(DN_v(3,0)*crLHS97 + DN_v(3,1)*crLHS157 - crLHS160*crLHS242 - crLHS160*crLHS244 + crLHS265);
    rLHS(7,11)+=gauss_weight*(DN_v(3,0)*crLHS108 + DN_v(3,1)*crLHS161 - crLHS164*crLHS242 - crLHS164*crLHS244 + crLHS257 + crLHS266);
    rLHS(7,12)+=gauss_weight*(crLHS166*crLHS243 + crLHS166*crLHS67 + crLHS267);
    rLHS(7,13)+=gauss_weight*(crLHS168*crLHS243 + crLHS168*crLHS67 + crLHS268);
    rLHS(7,14)+=gauss_weight*(crLHS170*crLHS243 + crLHS170*crLHS67 + crLHS269);
    rLHS(8,0)+=gauss_weight*(DN_v(4,0)*crLHS0 + DN_v(4,1)*crLHS2 - crLHS14*crLHS270 - crLHS14*crLHS272 + crLHS273 + crLHS81);
    rLHS(8,1)+=gauss_weight*(DN_v(4,0)*crLHS21 + DN_v(4,1)*crLHS23 + crLHS150 - crLHS26*crLHS270 - crLHS26*crLHS272);
    rLHS(8,2)+=gauss_weight*(DN_v(4,0)*crLHS27 + DN_v(4,1)*crLHS29 + crLHS186 - crLHS270*crLHS36 - crLHS272*crLHS36 + crLHS274);
    rLHS(8,3)+=gauss_weight*(DN_v(4,0)*crLHS39 + DN_v(4,1)*crLHS41 + crLHS202 - crLHS270*crLHS43 - crLHS272*crLHS43);
    rLHS(8,4)+=gauss_weight*(DN_v(4,0)*crLHS44 + DN_v(4,1)*crLHS46 + crLHS221 - crLHS270*crLHS53 - crLHS272*crLHS53 + crLHS275);
    rLHS(8,5)+=gauss_weight*(DN_v(4,0)*crLHS56 + DN_v(4,1)*crLHS58 + crLHS235 - crLHS270*crLHS60 - crLHS272*crLHS60);
    rLHS(8,6)+=gauss_weight*(DN_v(4,0)*crLHS61 + DN_v(4,1)*crLHS63 + crLHS251 - crLHS270*crLHS70 - crLHS272*crLHS70 + crLHS276);
    rLHS(8,7)+=gauss_weight*(DN_v(4,0)*crLHS73 + DN_v(4,1)*crLHS75 + crLHS263 - crLHS270*crLHS77 - crLHS272*crLHS77);
    rLHS(8,8)+=gauss_weight*(std::pow(DN_v(4,0), 2)*crLHS6 + DN_v(4,0)*crLHS78 + DN_v(4,1)*crLHS80 - crLHS270*crLHS87 - crLHS272*crLHS87 + crLHS277);
    rLHS(8,9)+=gauss_weight*(DN_v(4,0)*crLHS90 + DN_v(4,1)*crLHS92 - crLHS270*crLHS94 - crLHS272*crLHS94 + crLHS279);
    rLHS(8,10)+=gauss_weight*(DN_v(4,0)*crLHS95 + DN_v(4,1)*crLHS97 - crLHS103*crLHS270 - crLHS103*crLHS272 + crLHS280 + crLHS282);
    rLHS(8,11)+=gauss_weight*(DN_v(4,0)*crLHS106 + DN_v(4,1)*crLHS108 - crLHS110*crLHS270 - crLHS110*crLHS272 + crLHS283);
    rLHS(8,12)+=gauss_weight*(crLHS112*crLHS271 + crLHS112*crLHS84 + crLHS284);
    rLHS(8,13)+=gauss_weight*(crLHS114*crLHS271 + crLHS114*crLHS84 + crLHS285);
    rLHS(8,14)+=gauss_weight*(crLHS116*crLHS271 + crLHS116*crLHS84 + crLHS286);
    rLHS(9,0)+=gauss_weight*(DN_v(4,0)*crLHS2 + DN_v(4,1)*crLHS117 - crLHS119*crLHS270 - crLHS119*crLHS272 + crLHS93);
    rLHS(9,1)+=gauss_weight*(DN_v(4,0)*crLHS23 + DN_v(4,1)*crLHS120 - crLHS123*crLHS270 - crLHS123*crLHS272 + crLHS154 + crLHS273);
    rLHS(9,2)+=gauss_weight*(DN_v(4,0)*crLHS29 + DN_v(4,1)*crLHS124 - crLHS128*crLHS270 - crLHS128*crLHS272 + crLHS189);
    rLHS(9,3)+=gauss_weight*(DN_v(4,0)*crLHS41 + DN_v(4,1)*crLHS129 - crLHS132*crLHS270 - crLHS132*crLHS272 + crLHS203 + crLHS274);
    rLHS(9,4)+=gauss_weight*(DN_v(4,0)*crLHS46 + DN_v(4,1)*crLHS133 - crLHS136*crLHS270 - crLHS136*crLHS272 + crLHS224);
    rLHS(9,5)+=gauss_weight*(DN_v(4,0)*crLHS58 + DN_v(4,1)*crLHS137 - crLHS140*crLHS270 - crLHS140*crLHS272 + crLHS236 + crLHS275);
    rLHS(9,6)+=gauss_weight*(DN_v(4,0)*crLHS63 + DN_v(4,1)*crLHS141 - crLHS144*crLHS270 - crLHS144*crLHS272 + crLHS254);
    rLHS(9,7)+=gauss_weight*(DN_v(4,0)*crLHS75 + DN_v(4,1)*crLHS145 - crLHS148*crLHS270 - crLHS148*crLHS272 + crLHS264 + crLHS276);
    rLHS(9,8)+=gauss_weight*(DN_v(4,0)*crLHS80 + DN_v(4,1)*crLHS149 - crLHS152*crLHS270 - crLHS152*crLHS272 + crLHS279);
    rLHS(9,9)+=gauss_weight*(DN_v(4,0)*crLHS92 + std::pow(DN_v(4,1), 2)*crLHS6 + DN_v(4,1)*crLHS153 - crLHS156*crLHS270 - crLHS156*crLHS272 + crLHS277);
    rLHS(9,10)+=gauss_weight*(DN_v(4,0)*crLHS97 + DN_v(4,1)*crLHS157 - crLHS160*crLHS270 - crLHS160*crLHS272 + crLHS288);
    rLHS(9,11)+=gauss_weight*(DN_v(4,0)*crLHS108 + DN_v(4,1)*crLHS161 - crLHS164*crLHS270 - crLHS164*crLHS272 + crLHS282 + crLHS289);
    rLHS(9,12)+=gauss_weight*(crLHS166*crLHS271 + crLHS166*crLHS84 + crLHS290);
    rLHS(9,13)+=gauss_weight*(crLHS168*crLHS271 + crLHS168*crLHS84 + crLHS291);
    rLHS(9,14)+=gauss_weight*(crLHS170*crLHS271 + crLHS170*crLHS84 + crLHS292);
    rLHS(10,0)+=gauss_weight*(DN_v(5,0)*crLHS0 + DN_v(5,1)*crLHS2 - crLHS14*crLHS293 - crLHS14*crLHS295 + crLHS296 + crLHS98);
    rLHS(10,1)+=gauss_weight*(DN_v(5,0)*crLHS21 + DN_v(5,1)*crLHS23 + crLHS158 - crLHS26*crLHS293 - crLHS26*crLHS295);
    rLHS(10,2)+=gauss_weight*(DN_v(5,0)*crLHS27 + DN_v(5,1)*crLHS29 + crLHS190 - crLHS293*crLHS36 - crLHS295*crLHS36 + crLHS297);
    rLHS(10,3)+=gauss_weight*(DN_v(5,0)*crLHS39 + DN_v(5,1)*crLHS41 + crLHS204 - crLHS293*crLHS43 - crLHS295*crLHS43);
    rLHS(10,4)+=gauss_weight*(DN_v(5,0)*crLHS44 + DN_v(5,1)*crLHS46 + crLHS225 - crLHS293*crLHS53 - crLHS295*crLHS53 + crLHS298);
    rLHS(10,5)+=gauss_weight*(DN_v(5,0)*crLHS56 + DN_v(5,1)*crLHS58 + crLHS237 - crLHS293*crLHS60 - crLHS295*crLHS60);
    rLHS(10,6)+=gauss_weight*(DN_v(5,0)*crLHS61 + DN_v(5,1)*crLHS63 + crLHS255 - crLHS293*crLHS70 - crLHS295*crLHS70 + crLHS299);
    rLHS(10,7)+=gauss_weight*(DN_v(5,0)*crLHS73 + DN_v(5,1)*crLHS75 + crLHS265 - crLHS293*crLHS77 - crLHS295*crLHS77);
    rLHS(10,8)+=gauss_weight*(DN_v(5,0)*crLHS78 + DN_v(5,1)*crLHS80 + crLHS280 - crLHS293*crLHS87 - crLHS295*crLHS87 + crLHS300);
    rLHS(10,9)+=gauss_weight*(DN_v(5,0)*crLHS90 + DN_v(5,1)*crLHS92 + crLHS288 - crLHS293*crLHS94 - crLHS295*crLHS94);
    rLHS(10,10)+=gauss_weight*(std::pow(DN_v(5,0), 2)*crLHS6 + DN_v(5,0)*crLHS95 + DN_v(5,1)*crLHS97 - crLHS103*crLHS293 - crLHS103*crLHS295 + crLHS301);
    rLHS(10,11)+=gauss_weight*(DN_v(5,0)*crLHS106 + DN_v(5,1)*crLHS108 - crLHS110*crLHS293 - crLHS110*crLHS295 + crLHS302);
    rLHS(10,12)+=gauss_weight*(crLHS101*crLHS112 + crLHS112*crLHS294 + crLHS303);
    rLHS(10,13)+=gauss_weight*(crLHS101*crLHS114 + crLHS114*crLHS294 + crLHS304);
    rLHS(10,14)+=gauss_weight*(crLHS101*crLHS116 + crLHS116*crLHS294 + crLHS305);
    rLHS(11,0)+=gauss_weight*(DN_v(5,0)*crLHS2 + DN_v(5,1)*crLHS117 + crLHS109 - crLHS119*crLHS293 - crLHS119*crLHS295);
    rLHS(11,1)+=gauss_weight*(DN_v(5,0)*crLHS23 + DN_v(5,1)*crLHS120 - crLHS123*crLHS293 - crLHS123*crLHS295 + crLHS162 + crLHS296);
    rLHS(11,2)+=gauss_weight*(DN_v(5,0)*crLHS29 + DN_v(5,1)*crLHS124 - crLHS128*crLHS293 - crLHS128*crLHS295 + crLHS193);
    rLHS(11,3)+=gauss_weight*(DN_v(5,0)*crLHS41 + DN_v(5,1)*crLHS129 - crLHS132*crLHS293 - crLHS132*crLHS295 + crLHS205 + crLHS297);
    rLHS(11,4)+=gauss_weight*(DN_v(5,0)*crLHS46 + DN_v(5,1)*crLHS133 - crLHS136*crLHS293 - crLHS136*crLHS295 + crLHS228);
    rLHS(11,5)+=gauss_weight*(DN_v(5,0)*crLHS58 + DN_v(5,1)*crLHS137 - crLHS140*crLHS293 - crLHS140*crLHS295 + crLHS238 + crLHS298);
    rLHS(11,6)+=gauss_weight*(DN_v(5,0)*crLHS63 + DN_v(5,1)*crLHS141 - crLHS144*crLHS293 - crLHS144*crLHS295 + crLHS258);
    rLHS(11,7)+=gauss_weight*(DN_v(5,0)*crLHS75 + DN_v(5,1)*crLHS145 - crLHS148*crLHS293 - crLHS148*crLHS295 + crLHS266 + crLHS299);
    rLHS(11,8)+=gauss_weight*(DN_v(5,0)*crLHS80 + DN_v(5,1)*crLHS149 - crLHS152*crLHS293 - crLHS152*crLHS295 + crLHS283);
    rLHS(11,9)+=gauss_weight*(DN_v(5,0)*crLHS92 + DN_v(5,1)*crLHS153 - crLHS156*crLHS293 - crLHS156*crLHS295 + crLHS289 + crLHS300);
    rLHS(11,10)+=gauss_weight*(DN_v(5,0)*crLHS97 + DN_v(5,1)*crLHS157 - crLHS160*crLHS293 - crLHS160*crLHS295 + crLHS302);
    rLHS(11,11)+=gauss_weight*(DN_v(5,0)*crLHS108 + std::pow(DN_v(5,1), 2)*crLHS6 + DN_v(5,1)*crLHS161 - crLHS164*crLHS293 - crLHS164*crLHS295 + crLHS301);
    rLHS(11,12)+=gauss_weight*(crLHS101*crLHS166 + crLHS166*crLHS294 + crLHS306);
    rLHS(11,13)+=gauss_weight*(crLHS101*crLHS168 + crLHS168*crLHS294 + crLHS307);
    rLHS(11,14)+=gauss_weight*(crLHS101*crLHS170 + crLHS170*crLHS294 + crLHS308);
    rLHS(12,0)+=-gauss_weight*(crLHS111 + crLHS112*crLHS14 + crLHS119*crLHS166);
    rLHS(12,1)+=-gauss_weight*(crLHS112*crLHS26 + crLHS123*crLHS166 + crLHS165);
    rLHS(12,2)+=-gauss_weight*(crLHS112*crLHS36 + crLHS128*crLHS166 + crLHS194);
    rLHS(12,3)+=-gauss_weight*(crLHS112*crLHS43 + crLHS132*crLHS166 + crLHS206);
    rLHS(12,4)+=-gauss_weight*(crLHS112*crLHS53 + crLHS136*crLHS166 + crLHS229);
    rLHS(12,5)+=-gauss_weight*(crLHS112*crLHS60 + crLHS140*crLHS166 + crLHS239);
    rLHS(12,6)+=-gauss_weight*(crLHS112*crLHS70 + crLHS144*crLHS166 + crLHS259);
    rLHS(12,7)+=-gauss_weight*(crLHS112*crLHS77 + crLHS148*crLHS166 + crLHS267);
    rLHS(12,8)+=-gauss_weight*(crLHS112*crLHS87 + crLHS152*crLHS166 + crLHS284);
    rLHS(12,9)+=-gauss_weight*(crLHS112*crLHS94 + crLHS156*crLHS166 + crLHS290);
    rLHS(12,10)+=-gauss_weight*(crLHS103*crLHS112 + crLHS160*crLHS166 + crLHS303);
    rLHS(12,11)+=-gauss_weight*(crLHS110*crLHS112 + crLHS164*crLHS166 + crLHS306);
    rLHS(12,12)+=crLHS309*(std::pow(DN_p(0,0), 2) + std::pow(DN_p(0,1), 2));
    rLHS(12,13)+=crLHS310;
    rLHS(12,14)+=crLHS311;
    rLHS(13,0)+=-gauss_weight*(crLHS113 + crLHS114*crLHS14 + crLHS119*crLHS168);
    rLHS(13,1)+=-gauss_weight*(crLHS114*crLHS26 + crLHS123*crLHS168 + crLHS167);
    rLHS(13,2)+=-gauss_weight*(crLHS114*crLHS36 + crLHS128*crLHS168 + crLHS195);
    rLHS(13,3)+=-gauss_weight*(crLHS114*crLHS43 + crLHS132*crLHS168 + crLHS207);
    rLHS(13,4)+=-gauss_weight*(crLHS114*crLHS53 + crLHS136*crLHS168 + crLHS230);
    rLHS(13,5)+=-gauss_weight*(crLHS114*crLHS60 + crLHS140*crLHS168 + crLHS240);
    rLHS(13,6)+=-gauss_weight*(crLHS114*crLHS70 + crLHS144*crLHS168 + crLHS260);
    rLHS(13,7)+=-gauss_weight*(crLHS114*crLHS77 + crLHS148*crLHS168 + crLHS268);
    rLHS(13,8)+=-gauss_weight*(crLHS114*crLHS87 + crLHS152*crLHS168 + crLHS285);
    rLHS(13,9)+=-gauss_weight*(crLHS114*crLHS94 + crLHS156*crLHS168 + crLHS291);
    rLHS(13,10)+=-gauss_weight*(crLHS103*crLHS114 + crLHS160*crLHS168 + crLHS304);
    rLHS(13,11)+=-gauss_weight*(crLHS110*crLHS114 + crLHS164*crLHS168 + crLHS307);
    rLHS(13,12)+=crLHS310;
    rLHS(13,13)+=crLHS309*(std::pow(DN_p(1,0), 2) + std::pow(DN_p(1,1), 2));
    rLHS(13,14)+=crLHS312;
    rLHS(14,0)+=-gauss_weight*(crLHS115 + crLHS116*crLHS14 + crLHS119*crLHS170);
    rLHS(14,1)+=-gauss_weight*(crLHS116*crLHS26 + crLHS123*crLHS170 + crLHS169);
    rLHS(14,2)+=-gauss_weight*(crLHS116*crLHS36 + crLHS128*crLHS170 + crLHS196);
    rLHS(14,3)+=-gauss_weight*(crLHS116*crLHS43 + crLHS132*crLHS170 + crLHS208);
    rLHS(14,4)+=-gauss_weight*(crLHS116*crLHS53 + crLHS136*crLHS170 + crLHS231);
    rLHS(14,5)+=-gauss_weight*(crLHS116*crLHS60 + crLHS140*crLHS170 + crLHS241);
    rLHS(14,6)+=-gauss_weight*(crLHS116*crLHS70 + crLHS144*crLHS170 + crLHS261);
    rLHS(14,7)+=-gauss_weight*(crLHS116*crLHS77 + crLHS148*crLHS170 + crLHS269);
    rLHS(14,8)+=-gauss_weight*(crLHS116*crLHS87 + crLHS152*crLHS170 + crLHS286);
    rLHS(14,9)+=-gauss_weight*(crLHS116*crLHS94 + crLHS156*crLHS170 + crLHS292);
    rLHS(14,10)+=-gauss_weight*(crLHS103*crLHS116 + crLHS160*crLHS170 + crLHS305);
    rLHS(14,11)+=-gauss_weight*(crLHS110*crLHS116 + crLHS164*crLHS170 + crLHS308);
    rLHS(14,12)+=crLHS311;
    rLHS(14,13)+=crLHS312;
    rLHS(14,14)+=crLHS309*(std::pow(DN_p(2,0), 2) + std::pow(DN_p(2,1), 2));

}

template <>
void IncompressibleNavierStokesP2P1Continuous<3>::AddGaussPointLeftHandSideContribution(
    const ElementDataContainer& rData,
    MatrixType& rLHS)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

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
    const double crLHS0 = C(0,0)*DN_v(0,0) + C(0,3)*DN_v(0,1) + C(0,5)*DN_v(0,2);
    const double crLHS1 = C(0,3)*DN_v(0,0);
    const double crLHS2 = C(3,3)*DN_v(0,1) + C(3,5)*DN_v(0,2) + crLHS1;
    const double crLHS3 = C(0,5)*DN_v(0,0);
    const double crLHS4 = C(3,5)*DN_v(0,1) + C(5,5)*DN_v(0,2) + crLHS3;
    const double crLHS5 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0) + N_v[6]*vconv(6,0) + N_v[7]*vconv(7,0) + N_v[8]*vconv(8,0) + N_v[9]*vconv(9,0);
    const double crLHS6 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1) + N_v[6]*vconv(6,1) + N_v[7]*vconv(7,1) + N_v[8]*vconv(8,1) + N_v[9]*vconv(9,1);
    const double crLHS7 = N_v[0]*vconv(0,2) + N_v[1]*vconv(1,2) + N_v[2]*vconv(2,2) + N_v[3]*vconv(3,2) + N_v[4]*vconv(4,2) + N_v[5]*vconv(5,2) + N_v[6]*vconv(6,2) + N_v[7]*vconv(7,2) + N_v[8]*vconv(8,2) + N_v[9]*vconv(9,2);
    const double crLHS8 = rho*stab_c2*std::sqrt(std::pow(crLHS5, 2) + std::pow(crLHS6, 2) + std::pow(crLHS7, 2));
    const double crLHS9 = crLHS8*h/stab_c1 + mu;
    const double crLHS10 = 1.0*C(0,0);
    const double crLHS11 = C(0,3)*DDN_v[0](0,0);
    const double crLHS12 = 1.0*DDN_v[0](0,1);
    const double crLHS13 = C(0,5)*DDN_v[0](0,0);
    const double crLHS14 = 1.0*DDN_v[0](0,2);
    const double crLHS15 = rho*(DN_v(0,0)*crLHS5 + DN_v(0,1)*crLHS6 + DN_v(0,2)*crLHS7);
    const double crLHS16 = rData.BDF0*rho;
    const double crLHS17 = N_v[0]*crLHS16;
    const double crLHS18 = -crLHS15 - crLHS17;
    const double crLHS19 = C(0,3)*crLHS12 + C(0,5)*crLHS14 + C(3,3)*crLHS12 + C(3,5)*crLHS12 + C(3,5)*crLHS14 + C(5,5)*crLHS14 + DDN_v[0](0,0)*crLHS10 + 1.0*crLHS11 + 1.0*crLHS13 + crLHS18;
    const double crLHS20 = 1.0/(crLHS8/h + dyn_tau*rho/rData.DeltaTime + mu*stab_c1/std::pow(h, 2));
    const double crLHS21 = crLHS15*crLHS20;
    const double crLHS22 = rho*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(0,2)*vconv(0,2) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(1,2)*vconv(1,2) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(2,2)*vconv(2,2) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(3,2)*vconv(3,2) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(4,2)*vconv(4,2) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1) + DN_v(5,2)*vconv(5,2) + DN_v(6,0)*vconv(6,0) + DN_v(6,1)*vconv(6,1) + DN_v(6,2)*vconv(6,2) + DN_v(7,0)*vconv(7,0) + DN_v(7,1)*vconv(7,1) + DN_v(7,2)*vconv(7,2) + DN_v(8,0)*vconv(8,0) + DN_v(8,1)*vconv(8,1) + DN_v(8,2)*vconv(8,2) + DN_v(9,0)*vconv(9,0) + DN_v(9,1)*vconv(9,1) + DN_v(9,2)*vconv(9,2));
    const double crLHS23 = N_v[0]*crLHS22;
    const double crLHS24 = crLHS20*crLHS23;
    const double crLHS25 = std::pow(N_v[0], 2)*crLHS16 + N_v[0]*crLHS15;
    const double crLHS26 = C(0,1)*DN_v(0,1) + C(0,4)*DN_v(0,2) + crLHS1;
    const double crLHS27 = C(1,3)*DN_v(0,1);
    const double crLHS28 = C(3,3)*DN_v(0,0) + C(3,4)*DN_v(0,2) + crLHS27;
    const double crLHS29 = C(3,5)*DN_v(0,0);
    const double crLHS30 = C(4,5)*DN_v(0,2);
    const double crLHS31 = C(1,5)*DN_v(0,1) + crLHS29 + crLHS30;
    const double crLHS32 = DN_v(0,0)*crLHS9;
    const double crLHS33 = DN_v(0,1)*crLHS32;
    const double crLHS34 = C(3,5)*DDN_v[0](0,0);
    const double crLHS35 = C(0,1)*DDN_v[0](0,1) + C(0,4)*DDN_v[0](0,2) + C(1,3)*DDN_v[0](0,1) + C(1,5)*DDN_v[0](0,1) + C(3,3)*DDN_v[0](0,0) + C(3,4)*DDN_v[0](0,2) + C(4,5)*DDN_v[0](0,2) + crLHS11 + crLHS34;
    const double crLHS36 = C(0,2)*DN_v(0,2) + C(0,4)*DN_v(0,1) + crLHS3;
    const double crLHS37 = C(3,4)*DN_v(0,1);
    const double crLHS38 = C(2,3)*DN_v(0,2) + crLHS29 + crLHS37;
    const double crLHS39 = C(2,5)*DN_v(0,2);
    const double crLHS40 = C(4,5)*DN_v(0,1) + C(5,5)*DN_v(0,0) + crLHS39;
    const double crLHS41 = DN_v(0,2)*crLHS32;
    const double crLHS42 = C(0,2)*DDN_v[0](0,2) + C(0,4)*DDN_v[0](0,1) + C(2,3)*DDN_v[0](0,2) + C(2,5)*DDN_v[0](0,2) + C(3,4)*DDN_v[0](0,1) + C(4,5)*DDN_v[0](0,1) + C(5,5)*DDN_v[0](0,0) + crLHS13 + crLHS34;
    const double crLHS43 = C(0,0)*DN_v(1,0) + C(0,3)*DN_v(1,1) + C(0,5)*DN_v(1,2);
    const double crLHS44 = C(0,3)*DN_v(1,0);
    const double crLHS45 = C(3,3)*DN_v(1,1) + C(3,5)*DN_v(1,2) + crLHS44;
    const double crLHS46 = C(0,5)*DN_v(1,0);
    const double crLHS47 = C(3,5)*DN_v(1,1) + C(5,5)*DN_v(1,2) + crLHS46;
    const double crLHS48 = DN_v(1,0)*crLHS32;
    const double crLHS49 = C(0,3)*DDN_v[1](0,0);
    const double crLHS50 = 1.0*DDN_v[1](0,1);
    const double crLHS51 = C(0,5)*DDN_v[1](0,0);
    const double crLHS52 = 1.0*DDN_v[1](0,2);
    const double crLHS53 = rho*(DN_v(1,0)*crLHS5 + DN_v(1,1)*crLHS6 + DN_v(1,2)*crLHS7);
    const double crLHS54 = N_v[1]*crLHS16;
    const double crLHS55 = -crLHS53 - crLHS54;
    const double crLHS56 = C(0,3)*crLHS50 + C(0,5)*crLHS52 + C(3,3)*crLHS50 + C(3,5)*crLHS50 + C(3,5)*crLHS52 + C(5,5)*crLHS52 + DDN_v[1](0,0)*crLHS10 + 1.0*crLHS49 + 1.0*crLHS51 + crLHS55;
    const double crLHS57 = N_v[1]*crLHS17;
    const double crLHS58 = N_v[0]*crLHS53 + crLHS57;
    const double crLHS59 = C(0,1)*DN_v(1,1) + C(0,4)*DN_v(1,2) + crLHS44;
    const double crLHS60 = C(1,3)*DN_v(1,1);
    const double crLHS61 = C(3,3)*DN_v(1,0) + C(3,4)*DN_v(1,2) + crLHS60;
    const double crLHS62 = C(3,5)*DN_v(1,0);
    const double crLHS63 = C(4,5)*DN_v(1,2);
    const double crLHS64 = C(1,5)*DN_v(1,1) + crLHS62 + crLHS63;
    const double crLHS65 = DN_v(1,1)*crLHS32;
    const double crLHS66 = C(3,5)*DDN_v[1](0,0);
    const double crLHS67 = C(0,1)*DDN_v[1](0,1) + C(0,4)*DDN_v[1](0,2) + C(1,3)*DDN_v[1](0,1) + C(1,5)*DDN_v[1](0,1) + C(3,3)*DDN_v[1](0,0) + C(3,4)*DDN_v[1](0,2) + C(4,5)*DDN_v[1](0,2) + crLHS49 + crLHS66;
    const double crLHS68 = C(0,2)*DN_v(1,2) + C(0,4)*DN_v(1,1) + crLHS46;
    const double crLHS69 = C(3,4)*DN_v(1,1);
    const double crLHS70 = C(2,3)*DN_v(1,2) + crLHS62 + crLHS69;
    const double crLHS71 = C(2,5)*DN_v(1,2);
    const double crLHS72 = C(4,5)*DN_v(1,1) + C(5,5)*DN_v(1,0) + crLHS71;
    const double crLHS73 = DN_v(1,2)*crLHS32;
    const double crLHS74 = C(0,2)*DDN_v[1](0,2) + C(0,4)*DDN_v[1](0,1) + C(2,3)*DDN_v[1](0,2) + C(2,5)*DDN_v[1](0,2) + C(3,4)*DDN_v[1](0,1) + C(4,5)*DDN_v[1](0,1) + C(5,5)*DDN_v[1](0,0) + crLHS51 + crLHS66;
    const double crLHS75 = C(0,0)*DN_v(2,0) + C(0,3)*DN_v(2,1) + C(0,5)*DN_v(2,2);
    const double crLHS76 = C(0,3)*DN_v(2,0);
    const double crLHS77 = C(3,3)*DN_v(2,1) + C(3,5)*DN_v(2,2) + crLHS76;
    const double crLHS78 = C(0,5)*DN_v(2,0);
    const double crLHS79 = C(3,5)*DN_v(2,1) + C(5,5)*DN_v(2,2) + crLHS78;
    const double crLHS80 = DN_v(2,0)*crLHS32;
    const double crLHS81 = C(0,3)*DDN_v[2](0,0);
    const double crLHS82 = 1.0*DDN_v[2](0,1);
    const double crLHS83 = C(0,5)*DDN_v[2](0,0);
    const double crLHS84 = 1.0*DDN_v[2](0,2);
    const double crLHS85 = rho*(DN_v(2,0)*crLHS5 + DN_v(2,1)*crLHS6 + DN_v(2,2)*crLHS7);
    const double crLHS86 = N_v[2]*crLHS16;
    const double crLHS87 = -crLHS85 - crLHS86;
    const double crLHS88 = C(0,3)*crLHS82 + C(0,5)*crLHS84 + C(3,3)*crLHS82 + C(3,5)*crLHS82 + C(3,5)*crLHS84 + C(5,5)*crLHS84 + DDN_v[2](0,0)*crLHS10 + 1.0*crLHS81 + 1.0*crLHS83 + crLHS87;
    const double crLHS89 = N_v[2]*crLHS17;
    const double crLHS90 = N_v[0]*crLHS85 + crLHS89;
    const double crLHS91 = C(0,1)*DN_v(2,1) + C(0,4)*DN_v(2,2) + crLHS76;
    const double crLHS92 = C(1,3)*DN_v(2,1);
    const double crLHS93 = C(3,3)*DN_v(2,0) + C(3,4)*DN_v(2,2) + crLHS92;
    const double crLHS94 = C(3,5)*DN_v(2,0);
    const double crLHS95 = C(4,5)*DN_v(2,2);
    const double crLHS96 = C(1,5)*DN_v(2,1) + crLHS94 + crLHS95;
    const double crLHS97 = DN_v(2,1)*crLHS32;
    const double crLHS98 = C(3,5)*DDN_v[2](0,0);
    const double crLHS99 = C(0,1)*DDN_v[2](0,1) + C(0,4)*DDN_v[2](0,2) + C(1,3)*DDN_v[2](0,1) + C(1,5)*DDN_v[2](0,1) + C(3,3)*DDN_v[2](0,0) + C(3,4)*DDN_v[2](0,2) + C(4,5)*DDN_v[2](0,2) + crLHS81 + crLHS98;
    const double crLHS100 = C(0,2)*DN_v(2,2) + C(0,4)*DN_v(2,1) + crLHS78;
    const double crLHS101 = C(3,4)*DN_v(2,1);
    const double crLHS102 = C(2,3)*DN_v(2,2) + crLHS101 + crLHS94;
    const double crLHS103 = C(2,5)*DN_v(2,2);
    const double crLHS104 = C(4,5)*DN_v(2,1) + C(5,5)*DN_v(2,0) + crLHS103;
    const double crLHS105 = DN_v(2,2)*crLHS32;
    const double crLHS106 = C(0,2)*DDN_v[2](0,2) + C(0,4)*DDN_v[2](0,1) + C(2,3)*DDN_v[2](0,2) + C(2,5)*DDN_v[2](0,2) + C(3,4)*DDN_v[2](0,1) + C(4,5)*DDN_v[2](0,1) + C(5,5)*DDN_v[2](0,0) + crLHS83 + crLHS98;
    const double crLHS107 = C(0,0)*DN_v(3,0) + C(0,3)*DN_v(3,1) + C(0,5)*DN_v(3,2);
    const double crLHS108 = C(0,3)*DN_v(3,0);
    const double crLHS109 = C(3,3)*DN_v(3,1) + C(3,5)*DN_v(3,2) + crLHS108;
    const double crLHS110 = C(0,5)*DN_v(3,0);
    const double crLHS111 = C(3,5)*DN_v(3,1) + C(5,5)*DN_v(3,2) + crLHS110;
    const double crLHS112 = DN_v(3,0)*crLHS32;
    const double crLHS113 = C(0,3)*DDN_v[3](0,0);
    const double crLHS114 = 1.0*DDN_v[3](0,1);
    const double crLHS115 = C(0,5)*DDN_v[3](0,0);
    const double crLHS116 = 1.0*DDN_v[3](0,2);
    const double crLHS117 = rho*(DN_v(3,0)*crLHS5 + DN_v(3,1)*crLHS6 + DN_v(3,2)*crLHS7);
    const double crLHS118 = N_v[3]*crLHS16;
    const double crLHS119 = -crLHS117 - crLHS118;
    const double crLHS120 = C(0,3)*crLHS114 + C(0,5)*crLHS116 + C(3,3)*crLHS114 + C(3,5)*crLHS114 + C(3,5)*crLHS116 + C(5,5)*crLHS116 + DDN_v[3](0,0)*crLHS10 + 1.0*crLHS113 + 1.0*crLHS115 + crLHS119;
    const double crLHS121 = N_v[3]*crLHS17;
    const double crLHS122 = N_v[0]*crLHS117 + crLHS121;
    const double crLHS123 = C(0,1)*DN_v(3,1) + C(0,4)*DN_v(3,2) + crLHS108;
    const double crLHS124 = C(1,3)*DN_v(3,1);
    const double crLHS125 = C(3,3)*DN_v(3,0) + C(3,4)*DN_v(3,2) + crLHS124;
    const double crLHS126 = C(3,5)*DN_v(3,0);
    const double crLHS127 = C(4,5)*DN_v(3,2);
    const double crLHS128 = C(1,5)*DN_v(3,1) + crLHS126 + crLHS127;
    const double crLHS129 = DN_v(3,1)*crLHS32;
    const double crLHS130 = C(3,5)*DDN_v[3](0,0);
    const double crLHS131 = C(0,1)*DDN_v[3](0,1) + C(0,4)*DDN_v[3](0,2) + C(1,3)*DDN_v[3](0,1) + C(1,5)*DDN_v[3](0,1) + C(3,3)*DDN_v[3](0,0) + C(3,4)*DDN_v[3](0,2) + C(4,5)*DDN_v[3](0,2) + crLHS113 + crLHS130;
    const double crLHS132 = C(0,2)*DN_v(3,2) + C(0,4)*DN_v(3,1) + crLHS110;
    const double crLHS133 = C(3,4)*DN_v(3,1);
    const double crLHS134 = C(2,3)*DN_v(3,2) + crLHS126 + crLHS133;
    const double crLHS135 = C(2,5)*DN_v(3,2);
    const double crLHS136 = C(4,5)*DN_v(3,1) + C(5,5)*DN_v(3,0) + crLHS135;
    const double crLHS137 = DN_v(3,2)*crLHS32;
    const double crLHS138 = C(0,2)*DDN_v[3](0,2) + C(0,4)*DDN_v[3](0,1) + C(2,3)*DDN_v[3](0,2) + C(2,5)*DDN_v[3](0,2) + C(3,4)*DDN_v[3](0,1) + C(4,5)*DDN_v[3](0,1) + C(5,5)*DDN_v[3](0,0) + crLHS115 + crLHS130;
    const double crLHS139 = C(0,0)*DN_v(4,0) + C(0,3)*DN_v(4,1) + C(0,5)*DN_v(4,2);
    const double crLHS140 = C(0,3)*DN_v(4,0);
    const double crLHS141 = C(3,3)*DN_v(4,1) + C(3,5)*DN_v(4,2) + crLHS140;
    const double crLHS142 = C(0,5)*DN_v(4,0);
    const double crLHS143 = C(3,5)*DN_v(4,1) + C(5,5)*DN_v(4,2) + crLHS142;
    const double crLHS144 = DN_v(4,0)*crLHS32;
    const double crLHS145 = C(0,3)*DDN_v[4](0,0);
    const double crLHS146 = 1.0*DDN_v[4](0,1);
    const double crLHS147 = C(0,5)*DDN_v[4](0,0);
    const double crLHS148 = 1.0*DDN_v[4](0,2);
    const double crLHS149 = rho*(DN_v(4,0)*crLHS5 + DN_v(4,1)*crLHS6 + DN_v(4,2)*crLHS7);
    const double crLHS150 = N_v[4]*crLHS16;
    const double crLHS151 = -crLHS149 - crLHS150;
    const double crLHS152 = C(0,3)*crLHS146 + C(0,5)*crLHS148 + C(3,3)*crLHS146 + C(3,5)*crLHS146 + C(3,5)*crLHS148 + C(5,5)*crLHS148 + DDN_v[4](0,0)*crLHS10 + 1.0*crLHS145 + 1.0*crLHS147 + crLHS151;
    const double crLHS153 = N_v[4]*crLHS17;
    const double crLHS154 = N_v[0]*crLHS149 + crLHS153;
    const double crLHS155 = C(0,1)*DN_v(4,1) + C(0,4)*DN_v(4,2) + crLHS140;
    const double crLHS156 = C(1,3)*DN_v(4,1);
    const double crLHS157 = C(3,3)*DN_v(4,0) + C(3,4)*DN_v(4,2) + crLHS156;
    const double crLHS158 = C(3,5)*DN_v(4,0);
    const double crLHS159 = C(4,5)*DN_v(4,2);
    const double crLHS160 = C(1,5)*DN_v(4,1) + crLHS158 + crLHS159;
    const double crLHS161 = DN_v(4,1)*crLHS32;
    const double crLHS162 = C(3,5)*DDN_v[4](0,0);
    const double crLHS163 = C(0,1)*DDN_v[4](0,1) + C(0,4)*DDN_v[4](0,2) + C(1,3)*DDN_v[4](0,1) + C(1,5)*DDN_v[4](0,1) + C(3,3)*DDN_v[4](0,0) + C(3,4)*DDN_v[4](0,2) + C(4,5)*DDN_v[4](0,2) + crLHS145 + crLHS162;
    const double crLHS164 = C(0,2)*DN_v(4,2) + C(0,4)*DN_v(4,1) + crLHS142;
    const double crLHS165 = C(3,4)*DN_v(4,1);
    const double crLHS166 = C(2,3)*DN_v(4,2) + crLHS158 + crLHS165;
    const double crLHS167 = C(2,5)*DN_v(4,2);
    const double crLHS168 = C(4,5)*DN_v(4,1) + C(5,5)*DN_v(4,0) + crLHS167;
    const double crLHS169 = DN_v(4,2)*crLHS32;
    const double crLHS170 = C(0,2)*DDN_v[4](0,2) + C(0,4)*DDN_v[4](0,1) + C(2,3)*DDN_v[4](0,2) + C(2,5)*DDN_v[4](0,2) + C(3,4)*DDN_v[4](0,1) + C(4,5)*DDN_v[4](0,1) + C(5,5)*DDN_v[4](0,0) + crLHS147 + crLHS162;
    const double crLHS171 = C(0,0)*DN_v(5,0) + C(0,3)*DN_v(5,1) + C(0,5)*DN_v(5,2);
    const double crLHS172 = C(0,3)*DN_v(5,0);
    const double crLHS173 = C(3,3)*DN_v(5,1) + C(3,5)*DN_v(5,2) + crLHS172;
    const double crLHS174 = C(0,5)*DN_v(5,0);
    const double crLHS175 = C(3,5)*DN_v(5,1) + C(5,5)*DN_v(5,2) + crLHS174;
    const double crLHS176 = DN_v(5,0)*crLHS32;
    const double crLHS177 = C(0,3)*DDN_v[5](0,0);
    const double crLHS178 = 1.0*DDN_v[5](0,1);
    const double crLHS179 = C(0,5)*DDN_v[5](0,0);
    const double crLHS180 = 1.0*DDN_v[5](0,2);
    const double crLHS181 = rho*(DN_v(5,0)*crLHS5 + DN_v(5,1)*crLHS6 + DN_v(5,2)*crLHS7);
    const double crLHS182 = N_v[5]*crLHS16;
    const double crLHS183 = -crLHS181 - crLHS182;
    const double crLHS184 = C(0,3)*crLHS178 + C(0,5)*crLHS180 + C(3,3)*crLHS178 + C(3,5)*crLHS178 + C(3,5)*crLHS180 + C(5,5)*crLHS180 + DDN_v[5](0,0)*crLHS10 + 1.0*crLHS177 + 1.0*crLHS179 + crLHS183;
    const double crLHS185 = N_v[5]*crLHS17;
    const double crLHS186 = N_v[0]*crLHS181 + crLHS185;
    const double crLHS187 = C(0,1)*DN_v(5,1) + C(0,4)*DN_v(5,2) + crLHS172;
    const double crLHS188 = C(1,3)*DN_v(5,1);
    const double crLHS189 = C(3,3)*DN_v(5,0) + C(3,4)*DN_v(5,2) + crLHS188;
    const double crLHS190 = C(3,5)*DN_v(5,0);
    const double crLHS191 = C(4,5)*DN_v(5,2);
    const double crLHS192 = C(1,5)*DN_v(5,1) + crLHS190 + crLHS191;
    const double crLHS193 = DN_v(5,1)*crLHS32;
    const double crLHS194 = C(3,5)*DDN_v[5](0,0);
    const double crLHS195 = C(0,1)*DDN_v[5](0,1) + C(0,4)*DDN_v[5](0,2) + C(1,3)*DDN_v[5](0,1) + C(1,5)*DDN_v[5](0,1) + C(3,3)*DDN_v[5](0,0) + C(3,4)*DDN_v[5](0,2) + C(4,5)*DDN_v[5](0,2) + crLHS177 + crLHS194;
    const double crLHS196 = C(0,2)*DN_v(5,2) + C(0,4)*DN_v(5,1) + crLHS174;
    const double crLHS197 = C(3,4)*DN_v(5,1);
    const double crLHS198 = C(2,3)*DN_v(5,2) + crLHS190 + crLHS197;
    const double crLHS199 = C(2,5)*DN_v(5,2);
    const double crLHS200 = C(4,5)*DN_v(5,1) + C(5,5)*DN_v(5,0) + crLHS199;
    const double crLHS201 = DN_v(5,2)*crLHS32;
    const double crLHS202 = C(0,2)*DDN_v[5](0,2) + C(0,4)*DDN_v[5](0,1) + C(2,3)*DDN_v[5](0,2) + C(2,5)*DDN_v[5](0,2) + C(3,4)*DDN_v[5](0,1) + C(4,5)*DDN_v[5](0,1) + C(5,5)*DDN_v[5](0,0) + crLHS179 + crLHS194;
    const double crLHS203 = C(0,0)*DN_v(6,0) + C(0,3)*DN_v(6,1) + C(0,5)*DN_v(6,2);
    const double crLHS204 = C(0,3)*DN_v(6,0);
    const double crLHS205 = C(3,3)*DN_v(6,1) + C(3,5)*DN_v(6,2) + crLHS204;
    const double crLHS206 = C(0,5)*DN_v(6,0);
    const double crLHS207 = C(3,5)*DN_v(6,1) + C(5,5)*DN_v(6,2) + crLHS206;
    const double crLHS208 = DN_v(6,0)*crLHS32;
    const double crLHS209 = C(0,3)*DDN_v[6](0,0);
    const double crLHS210 = 1.0*DDN_v[6](0,1);
    const double crLHS211 = C(0,5)*DDN_v[6](0,0);
    const double crLHS212 = 1.0*DDN_v[6](0,2);
    const double crLHS213 = rho*(DN_v(6,0)*crLHS5 + DN_v(6,1)*crLHS6 + DN_v(6,2)*crLHS7);
    const double crLHS214 = N_v[6]*crLHS16;
    const double crLHS215 = -crLHS213 - crLHS214;
    const double crLHS216 = C(0,3)*crLHS210 + C(0,5)*crLHS212 + C(3,3)*crLHS210 + C(3,5)*crLHS210 + C(3,5)*crLHS212 + C(5,5)*crLHS212 + DDN_v[6](0,0)*crLHS10 + 1.0*crLHS209 + 1.0*crLHS211 + crLHS215;
    const double crLHS217 = N_v[6]*crLHS17;
    const double crLHS218 = N_v[0]*crLHS213 + crLHS217;
    const double crLHS219 = C(0,1)*DN_v(6,1) + C(0,4)*DN_v(6,2) + crLHS204;
    const double crLHS220 = C(1,3)*DN_v(6,1);
    const double crLHS221 = C(3,3)*DN_v(6,0) + C(3,4)*DN_v(6,2) + crLHS220;
    const double crLHS222 = C(3,5)*DN_v(6,0);
    const double crLHS223 = C(4,5)*DN_v(6,2);
    const double crLHS224 = C(1,5)*DN_v(6,1) + crLHS222 + crLHS223;
    const double crLHS225 = DN_v(6,1)*crLHS32;
    const double crLHS226 = C(3,5)*DDN_v[6](0,0);
    const double crLHS227 = C(0,1)*DDN_v[6](0,1) + C(0,4)*DDN_v[6](0,2) + C(1,3)*DDN_v[6](0,1) + C(1,5)*DDN_v[6](0,1) + C(3,3)*DDN_v[6](0,0) + C(3,4)*DDN_v[6](0,2) + C(4,5)*DDN_v[6](0,2) + crLHS209 + crLHS226;
    const double crLHS228 = C(0,2)*DN_v(6,2) + C(0,4)*DN_v(6,1) + crLHS206;
    const double crLHS229 = C(3,4)*DN_v(6,1);
    const double crLHS230 = C(2,3)*DN_v(6,2) + crLHS222 + crLHS229;
    const double crLHS231 = C(2,5)*DN_v(6,2);
    const double crLHS232 = C(4,5)*DN_v(6,1) + C(5,5)*DN_v(6,0) + crLHS231;
    const double crLHS233 = DN_v(6,2)*crLHS32;
    const double crLHS234 = C(0,2)*DDN_v[6](0,2) + C(0,4)*DDN_v[6](0,1) + C(2,3)*DDN_v[6](0,2) + C(2,5)*DDN_v[6](0,2) + C(3,4)*DDN_v[6](0,1) + C(4,5)*DDN_v[6](0,1) + C(5,5)*DDN_v[6](0,0) + crLHS211 + crLHS226;
    const double crLHS235 = C(0,0)*DN_v(7,0) + C(0,3)*DN_v(7,1) + C(0,5)*DN_v(7,2);
    const double crLHS236 = C(0,3)*DN_v(7,0);
    const double crLHS237 = C(3,3)*DN_v(7,1) + C(3,5)*DN_v(7,2) + crLHS236;
    const double crLHS238 = C(0,5)*DN_v(7,0);
    const double crLHS239 = C(3,5)*DN_v(7,1) + C(5,5)*DN_v(7,2) + crLHS238;
    const double crLHS240 = DN_v(7,0)*crLHS32;
    const double crLHS241 = C(0,3)*DDN_v[7](0,0);
    const double crLHS242 = 1.0*DDN_v[7](0,1);
    const double crLHS243 = C(0,5)*DDN_v[7](0,0);
    const double crLHS244 = 1.0*DDN_v[7](0,2);
    const double crLHS245 = rho*(DN_v(7,0)*crLHS5 + DN_v(7,1)*crLHS6 + DN_v(7,2)*crLHS7);
    const double crLHS246 = N_v[7]*crLHS16;
    const double crLHS247 = -crLHS245 - crLHS246;
    const double crLHS248 = C(0,3)*crLHS242 + C(0,5)*crLHS244 + C(3,3)*crLHS242 + C(3,5)*crLHS242 + C(3,5)*crLHS244 + C(5,5)*crLHS244 + DDN_v[7](0,0)*crLHS10 + 1.0*crLHS241 + 1.0*crLHS243 + crLHS247;
    const double crLHS249 = N_v[7]*crLHS17;
    const double crLHS250 = N_v[0]*crLHS245 + crLHS249;
    const double crLHS251 = C(0,1)*DN_v(7,1) + C(0,4)*DN_v(7,2) + crLHS236;
    const double crLHS252 = C(1,3)*DN_v(7,1);
    const double crLHS253 = C(3,3)*DN_v(7,0) + C(3,4)*DN_v(7,2) + crLHS252;
    const double crLHS254 = C(3,5)*DN_v(7,0);
    const double crLHS255 = C(4,5)*DN_v(7,2);
    const double crLHS256 = C(1,5)*DN_v(7,1) + crLHS254 + crLHS255;
    const double crLHS257 = DN_v(7,1)*crLHS32;
    const double crLHS258 = C(3,5)*DDN_v[7](0,0);
    const double crLHS259 = C(0,1)*DDN_v[7](0,1) + C(0,4)*DDN_v[7](0,2) + C(1,3)*DDN_v[7](0,1) + C(1,5)*DDN_v[7](0,1) + C(3,3)*DDN_v[7](0,0) + C(3,4)*DDN_v[7](0,2) + C(4,5)*DDN_v[7](0,2) + crLHS241 + crLHS258;
    const double crLHS260 = C(0,2)*DN_v(7,2) + C(0,4)*DN_v(7,1) + crLHS238;
    const double crLHS261 = C(3,4)*DN_v(7,1);
    const double crLHS262 = C(2,3)*DN_v(7,2) + crLHS254 + crLHS261;
    const double crLHS263 = C(2,5)*DN_v(7,2);
    const double crLHS264 = C(4,5)*DN_v(7,1) + C(5,5)*DN_v(7,0) + crLHS263;
    const double crLHS265 = DN_v(7,2)*crLHS32;
    const double crLHS266 = C(0,2)*DDN_v[7](0,2) + C(0,4)*DDN_v[7](0,1) + C(2,3)*DDN_v[7](0,2) + C(2,5)*DDN_v[7](0,2) + C(3,4)*DDN_v[7](0,1) + C(4,5)*DDN_v[7](0,1) + C(5,5)*DDN_v[7](0,0) + crLHS243 + crLHS258;
    const double crLHS267 = C(0,0)*DN_v(8,0) + C(0,3)*DN_v(8,1) + C(0,5)*DN_v(8,2);
    const double crLHS268 = C(0,3)*DN_v(8,0);
    const double crLHS269 = C(3,3)*DN_v(8,1) + C(3,5)*DN_v(8,2) + crLHS268;
    const double crLHS270 = C(0,5)*DN_v(8,0);
    const double crLHS271 = C(3,5)*DN_v(8,1) + C(5,5)*DN_v(8,2) + crLHS270;
    const double crLHS272 = DN_v(8,0)*crLHS32;
    const double crLHS273 = C(0,3)*DDN_v[8](0,0);
    const double crLHS274 = 1.0*DDN_v[8](0,1);
    const double crLHS275 = C(0,5)*DDN_v[8](0,0);
    const double crLHS276 = 1.0*DDN_v[8](0,2);
    const double crLHS277 = rho*(DN_v(8,0)*crLHS5 + DN_v(8,1)*crLHS6 + DN_v(8,2)*crLHS7);
    const double crLHS278 = N_v[8]*crLHS16;
    const double crLHS279 = -crLHS277 - crLHS278;
    const double crLHS280 = C(0,3)*crLHS274 + C(0,5)*crLHS276 + C(3,3)*crLHS274 + C(3,5)*crLHS274 + C(3,5)*crLHS276 + C(5,5)*crLHS276 + DDN_v[8](0,0)*crLHS10 + 1.0*crLHS273 + 1.0*crLHS275 + crLHS279;
    const double crLHS281 = N_v[8]*crLHS17;
    const double crLHS282 = N_v[0]*crLHS277 + crLHS281;
    const double crLHS283 = C(0,1)*DN_v(8,1) + C(0,4)*DN_v(8,2) + crLHS268;
    const double crLHS284 = C(1,3)*DN_v(8,1);
    const double crLHS285 = C(3,3)*DN_v(8,0) + C(3,4)*DN_v(8,2) + crLHS284;
    const double crLHS286 = C(3,5)*DN_v(8,0);
    const double crLHS287 = C(4,5)*DN_v(8,2);
    const double crLHS288 = C(1,5)*DN_v(8,1) + crLHS286 + crLHS287;
    const double crLHS289 = DN_v(8,1)*crLHS32;
    const double crLHS290 = C(3,5)*DDN_v[8](0,0);
    const double crLHS291 = C(0,1)*DDN_v[8](0,1) + C(0,4)*DDN_v[8](0,2) + C(1,3)*DDN_v[8](0,1) + C(1,5)*DDN_v[8](0,1) + C(3,3)*DDN_v[8](0,0) + C(3,4)*DDN_v[8](0,2) + C(4,5)*DDN_v[8](0,2) + crLHS273 + crLHS290;
    const double crLHS292 = C(0,2)*DN_v(8,2) + C(0,4)*DN_v(8,1) + crLHS270;
    const double crLHS293 = C(3,4)*DN_v(8,1);
    const double crLHS294 = C(2,3)*DN_v(8,2) + crLHS286 + crLHS293;
    const double crLHS295 = C(2,5)*DN_v(8,2);
    const double crLHS296 = C(4,5)*DN_v(8,1) + C(5,5)*DN_v(8,0) + crLHS295;
    const double crLHS297 = DN_v(8,2)*crLHS32;
    const double crLHS298 = C(0,2)*DDN_v[8](0,2) + C(0,4)*DDN_v[8](0,1) + C(2,3)*DDN_v[8](0,2) + C(2,5)*DDN_v[8](0,2) + C(3,4)*DDN_v[8](0,1) + C(4,5)*DDN_v[8](0,1) + C(5,5)*DDN_v[8](0,0) + crLHS275 + crLHS290;
    const double crLHS299 = C(0,0)*DN_v(9,0) + C(0,3)*DN_v(9,1) + C(0,5)*DN_v(9,2);
    const double crLHS300 = C(0,3)*DN_v(9,0);
    const double crLHS301 = C(3,3)*DN_v(9,1) + C(3,5)*DN_v(9,2) + crLHS300;
    const double crLHS302 = C(0,5)*DN_v(9,0);
    const double crLHS303 = C(3,5)*DN_v(9,1) + C(5,5)*DN_v(9,2) + crLHS302;
    const double crLHS304 = DN_v(9,0)*crLHS32;
    const double crLHS305 = C(0,3)*DDN_v[9](0,0);
    const double crLHS306 = 1.0*DDN_v[9](0,1);
    const double crLHS307 = C(0,5)*DDN_v[9](0,0);
    const double crLHS308 = 1.0*DDN_v[9](0,2);
    const double crLHS309 = rho*(DN_v(9,0)*crLHS5 + DN_v(9,1)*crLHS6 + DN_v(9,2)*crLHS7);
    const double crLHS310 = -N_v[9]*crLHS16 - crLHS309;
    const double crLHS311 = C(0,3)*crLHS306 + C(0,5)*crLHS308 + C(3,3)*crLHS306 + C(3,5)*crLHS306 + C(3,5)*crLHS308 + C(5,5)*crLHS308 + DDN_v[9](0,0)*crLHS10 + 1.0*crLHS305 + 1.0*crLHS307 + crLHS310;
    const double crLHS312 = N_v[9]*crLHS17;
    const double crLHS313 = N_v[0]*crLHS309 + crLHS312;
    const double crLHS314 = C(0,1)*DN_v(9,1) + C(0,4)*DN_v(9,2) + crLHS300;
    const double crLHS315 = C(1,3)*DN_v(9,1);
    const double crLHS316 = C(3,3)*DN_v(9,0) + C(3,4)*DN_v(9,2) + crLHS315;
    const double crLHS317 = C(3,5)*DN_v(9,0);
    const double crLHS318 = C(4,5)*DN_v(9,2);
    const double crLHS319 = C(1,5)*DN_v(9,1) + crLHS317 + crLHS318;
    const double crLHS320 = DN_v(9,1)*crLHS32;
    const double crLHS321 = C(3,5)*DDN_v[9](0,0);
    const double crLHS322 = C(0,1)*DDN_v[9](0,1) + C(0,4)*DDN_v[9](0,2) + C(1,3)*DDN_v[9](0,1) + C(1,5)*DDN_v[9](0,1) + C(3,3)*DDN_v[9](0,0) + C(3,4)*DDN_v[9](0,2) + C(4,5)*DDN_v[9](0,2) + crLHS305 + crLHS321;
    const double crLHS323 = C(0,2)*DN_v(9,2) + C(0,4)*DN_v(9,1) + crLHS302;
    const double crLHS324 = C(3,4)*DN_v(9,1);
    const double crLHS325 = C(2,3)*DN_v(9,2) + crLHS317 + crLHS324;
    const double crLHS326 = C(2,5)*DN_v(9,2);
    const double crLHS327 = C(4,5)*DN_v(9,1) + C(5,5)*DN_v(9,0) + crLHS326;
    const double crLHS328 = DN_v(9,2)*crLHS32;
    const double crLHS329 = C(0,2)*DDN_v[9](0,2) + C(0,4)*DDN_v[9](0,1) + C(2,3)*DDN_v[9](0,2) + C(2,5)*DDN_v[9](0,2) + C(3,4)*DDN_v[9](0,1) + C(4,5)*DDN_v[9](0,1) + C(5,5)*DDN_v[9](0,0) + crLHS307 + crLHS321;
    const double crLHS330 = -DN_v(0,0)*N_p[0];
    const double crLHS331 = DN_p(0,0)*crLHS20;
    const double crLHS332 = -DN_v(0,0)*N_p[1];
    const double crLHS333 = DN_p(1,0)*crLHS20;
    const double crLHS334 = -DN_v(0,0)*N_p[2];
    const double crLHS335 = DN_p(2,0)*crLHS20;
    const double crLHS336 = -DN_v(0,0)*N_p[3];
    const double crLHS337 = DN_p(3,0)*crLHS20;
    const double crLHS338 = C(0,1)*DN_v(0,0) + C(1,5)*DN_v(0,2) + crLHS27;
    const double crLHS339 = C(0,4)*DN_v(0,0) + crLHS30 + crLHS37;
    const double crLHS340 = C(1,3)*DDN_v[0](1,1);
    const double crLHS341 = C(3,4)*DDN_v[0](1,1);
    const double crLHS342 = C(0,1)*DDN_v[0](1,0) + C(0,3)*DDN_v[0](1,0) + C(0,4)*DDN_v[0](1,0) + C(1,5)*DDN_v[0](1,2) + C(3,3)*DDN_v[0](1,1) + C(3,5)*DDN_v[0](1,2) + C(4,5)*DDN_v[0](1,2) + crLHS340 + crLHS341;
    const double crLHS343 = C(1,1)*DN_v(0,1) + C(1,3)*DN_v(0,0) + C(1,4)*DN_v(0,2);
    const double crLHS344 = C(1,4)*DN_v(0,1);
    const double crLHS345 = C(3,4)*DN_v(0,0) + C(4,4)*DN_v(0,2) + crLHS344;
    const double crLHS346 = 1.0*C(1,1);
    const double crLHS347 = 1.0*DDN_v[0](1,0);
    const double crLHS348 = C(1,4)*DDN_v[0](1,1);
    const double crLHS349 = 1.0*DDN_v[0](1,2);
    const double crLHS350 = C(1,3)*crLHS347 + C(1,4)*crLHS349 + C(3,3)*crLHS347 + C(3,4)*crLHS347 + C(3,4)*crLHS349 + C(4,4)*crLHS349 + DDN_v[0](1,1)*crLHS346 + crLHS18 + 1.0*crLHS340 + 1.0*crLHS348;
    const double crLHS351 = C(1,2)*DN_v(0,2) + C(1,5)*DN_v(0,0) + crLHS344;
    const double crLHS352 = C(2,4)*DN_v(0,2);
    const double crLHS353 = C(4,4)*DN_v(0,1) + C(4,5)*DN_v(0,0) + crLHS352;
    const double crLHS354 = DN_v(0,1)*crLHS9;
    const double crLHS355 = DN_v(0,2)*crLHS354;
    const double crLHS356 = C(1,2)*DDN_v[0](1,2) + C(1,5)*DDN_v[0](1,0) + C(2,3)*DDN_v[0](1,2) + C(2,4)*DDN_v[0](1,2) + C(3,5)*DDN_v[0](1,0) + C(4,4)*DDN_v[0](1,1) + C(4,5)*DDN_v[0](1,0) + crLHS341 + crLHS348;
    const double crLHS357 = C(0,1)*DN_v(1,0) + C(1,5)*DN_v(1,2) + crLHS60;
    const double crLHS358 = C(0,4)*DN_v(1,0) + crLHS63 + crLHS69;
    const double crLHS359 = DN_v(1,0)*crLHS354;
    const double crLHS360 = C(1,3)*DDN_v[1](1,1);
    const double crLHS361 = C(3,4)*DDN_v[1](1,1);
    const double crLHS362 = C(0,1)*DDN_v[1](1,0) + C(0,3)*DDN_v[1](1,0) + C(0,4)*DDN_v[1](1,0) + C(1,5)*DDN_v[1](1,2) + C(3,3)*DDN_v[1](1,1) + C(3,5)*DDN_v[1](1,2) + C(4,5)*DDN_v[1](1,2) + crLHS360 + crLHS361;
    const double crLHS363 = C(1,1)*DN_v(1,1) + C(1,3)*DN_v(1,0) + C(1,4)*DN_v(1,2);
    const double crLHS364 = C(1,4)*DN_v(1,1);
    const double crLHS365 = C(3,4)*DN_v(1,0) + C(4,4)*DN_v(1,2) + crLHS364;
    const double crLHS366 = DN_v(1,1)*crLHS354;
    const double crLHS367 = 1.0*DDN_v[1](1,0);
    const double crLHS368 = C(1,4)*DDN_v[1](1,1);
    const double crLHS369 = 1.0*DDN_v[1](1,2);
    const double crLHS370 = C(1,3)*crLHS367 + C(1,4)*crLHS369 + C(3,3)*crLHS367 + C(3,4)*crLHS367 + C(3,4)*crLHS369 + C(4,4)*crLHS369 + DDN_v[1](1,1)*crLHS346 + 1.0*crLHS360 + 1.0*crLHS368 + crLHS55;
    const double crLHS371 = C(1,2)*DN_v(1,2) + C(1,5)*DN_v(1,0) + crLHS364;
    const double crLHS372 = C(2,4)*DN_v(1,2);
    const double crLHS373 = C(4,4)*DN_v(1,1) + C(4,5)*DN_v(1,0) + crLHS372;
    const double crLHS374 = DN_v(1,2)*crLHS354;
    const double crLHS375 = C(1,2)*DDN_v[1](1,2) + C(1,5)*DDN_v[1](1,0) + C(2,3)*DDN_v[1](1,2) + C(2,4)*DDN_v[1](1,2) + C(3,5)*DDN_v[1](1,0) + C(4,4)*DDN_v[1](1,1) + C(4,5)*DDN_v[1](1,0) + crLHS361 + crLHS368;
    const double crLHS376 = C(0,1)*DN_v(2,0) + C(1,5)*DN_v(2,2) + crLHS92;
    const double crLHS377 = C(0,4)*DN_v(2,0) + crLHS101 + crLHS95;
    const double crLHS378 = DN_v(2,0)*crLHS354;
    const double crLHS379 = C(1,3)*DDN_v[2](1,1);
    const double crLHS380 = C(3,4)*DDN_v[2](1,1);
    const double crLHS381 = C(0,1)*DDN_v[2](1,0) + C(0,3)*DDN_v[2](1,0) + C(0,4)*DDN_v[2](1,0) + C(1,5)*DDN_v[2](1,2) + C(3,3)*DDN_v[2](1,1) + C(3,5)*DDN_v[2](1,2) + C(4,5)*DDN_v[2](1,2) + crLHS379 + crLHS380;
    const double crLHS382 = C(1,1)*DN_v(2,1) + C(1,3)*DN_v(2,0) + C(1,4)*DN_v(2,2);
    const double crLHS383 = C(1,4)*DN_v(2,1);
    const double crLHS384 = C(3,4)*DN_v(2,0) + C(4,4)*DN_v(2,2) + crLHS383;
    const double crLHS385 = DN_v(2,1)*crLHS354;
    const double crLHS386 = 1.0*DDN_v[2](1,0);
    const double crLHS387 = C(1,4)*DDN_v[2](1,1);
    const double crLHS388 = 1.0*DDN_v[2](1,2);
    const double crLHS389 = C(1,3)*crLHS386 + C(1,4)*crLHS388 + C(3,3)*crLHS386 + C(3,4)*crLHS386 + C(3,4)*crLHS388 + C(4,4)*crLHS388 + DDN_v[2](1,1)*crLHS346 + 1.0*crLHS379 + 1.0*crLHS387 + crLHS87;
    const double crLHS390 = C(1,2)*DN_v(2,2) + C(1,5)*DN_v(2,0) + crLHS383;
    const double crLHS391 = C(2,4)*DN_v(2,2);
    const double crLHS392 = C(4,4)*DN_v(2,1) + C(4,5)*DN_v(2,0) + crLHS391;
    const double crLHS393 = DN_v(2,2)*crLHS354;
    const double crLHS394 = C(1,2)*DDN_v[2](1,2) + C(1,5)*DDN_v[2](1,0) + C(2,3)*DDN_v[2](1,2) + C(2,4)*DDN_v[2](1,2) + C(3,5)*DDN_v[2](1,0) + C(4,4)*DDN_v[2](1,1) + C(4,5)*DDN_v[2](1,0) + crLHS380 + crLHS387;
    const double crLHS395 = C(0,1)*DN_v(3,0) + C(1,5)*DN_v(3,2) + crLHS124;
    const double crLHS396 = C(0,4)*DN_v(3,0) + crLHS127 + crLHS133;
    const double crLHS397 = DN_v(3,0)*crLHS354;
    const double crLHS398 = C(1,3)*DDN_v[3](1,1);
    const double crLHS399 = C(3,4)*DDN_v[3](1,1);
    const double crLHS400 = C(0,1)*DDN_v[3](1,0) + C(0,3)*DDN_v[3](1,0) + C(0,4)*DDN_v[3](1,0) + C(1,5)*DDN_v[3](1,2) + C(3,3)*DDN_v[3](1,1) + C(3,5)*DDN_v[3](1,2) + C(4,5)*DDN_v[3](1,2) + crLHS398 + crLHS399;
    const double crLHS401 = C(1,1)*DN_v(3,1) + C(1,3)*DN_v(3,0) + C(1,4)*DN_v(3,2);
    const double crLHS402 = C(1,4)*DN_v(3,1);
    const double crLHS403 = C(3,4)*DN_v(3,0) + C(4,4)*DN_v(3,2) + crLHS402;
    const double crLHS404 = DN_v(3,1)*crLHS354;
    const double crLHS405 = 1.0*DDN_v[3](1,0);
    const double crLHS406 = C(1,4)*DDN_v[3](1,1);
    const double crLHS407 = 1.0*DDN_v[3](1,2);
    const double crLHS408 = C(1,3)*crLHS405 + C(1,4)*crLHS407 + C(3,3)*crLHS405 + C(3,4)*crLHS405 + C(3,4)*crLHS407 + C(4,4)*crLHS407 + DDN_v[3](1,1)*crLHS346 + crLHS119 + 1.0*crLHS398 + 1.0*crLHS406;
    const double crLHS409 = C(1,2)*DN_v(3,2) + C(1,5)*DN_v(3,0) + crLHS402;
    const double crLHS410 = C(2,4)*DN_v(3,2);
    const double crLHS411 = C(4,4)*DN_v(3,1) + C(4,5)*DN_v(3,0) + crLHS410;
    const double crLHS412 = DN_v(3,2)*crLHS354;
    const double crLHS413 = C(1,2)*DDN_v[3](1,2) + C(1,5)*DDN_v[3](1,0) + C(2,3)*DDN_v[3](1,2) + C(2,4)*DDN_v[3](1,2) + C(3,5)*DDN_v[3](1,0) + C(4,4)*DDN_v[3](1,1) + C(4,5)*DDN_v[3](1,0) + crLHS399 + crLHS406;
    const double crLHS414 = C(0,1)*DN_v(4,0) + C(1,5)*DN_v(4,2) + crLHS156;
    const double crLHS415 = C(0,4)*DN_v(4,0) + crLHS159 + crLHS165;
    const double crLHS416 = DN_v(4,0)*crLHS354;
    const double crLHS417 = C(1,3)*DDN_v[4](1,1);
    const double crLHS418 = C(3,4)*DDN_v[4](1,1);
    const double crLHS419 = C(0,1)*DDN_v[4](1,0) + C(0,3)*DDN_v[4](1,0) + C(0,4)*DDN_v[4](1,0) + C(1,5)*DDN_v[4](1,2) + C(3,3)*DDN_v[4](1,1) + C(3,5)*DDN_v[4](1,2) + C(4,5)*DDN_v[4](1,2) + crLHS417 + crLHS418;
    const double crLHS420 = C(1,1)*DN_v(4,1) + C(1,3)*DN_v(4,0) + C(1,4)*DN_v(4,2);
    const double crLHS421 = C(1,4)*DN_v(4,1);
    const double crLHS422 = C(3,4)*DN_v(4,0) + C(4,4)*DN_v(4,2) + crLHS421;
    const double crLHS423 = DN_v(4,1)*crLHS354;
    const double crLHS424 = 1.0*DDN_v[4](1,0);
    const double crLHS425 = C(1,4)*DDN_v[4](1,1);
    const double crLHS426 = 1.0*DDN_v[4](1,2);
    const double crLHS427 = C(1,3)*crLHS424 + C(1,4)*crLHS426 + C(3,3)*crLHS424 + C(3,4)*crLHS424 + C(3,4)*crLHS426 + C(4,4)*crLHS426 + DDN_v[4](1,1)*crLHS346 + crLHS151 + 1.0*crLHS417 + 1.0*crLHS425;
    const double crLHS428 = C(1,2)*DN_v(4,2) + C(1,5)*DN_v(4,0) + crLHS421;
    const double crLHS429 = C(2,4)*DN_v(4,2);
    const double crLHS430 = C(4,4)*DN_v(4,1) + C(4,5)*DN_v(4,0) + crLHS429;
    const double crLHS431 = DN_v(4,2)*crLHS354;
    const double crLHS432 = C(1,2)*DDN_v[4](1,2) + C(1,5)*DDN_v[4](1,0) + C(2,3)*DDN_v[4](1,2) + C(2,4)*DDN_v[4](1,2) + C(3,5)*DDN_v[4](1,0) + C(4,4)*DDN_v[4](1,1) + C(4,5)*DDN_v[4](1,0) + crLHS418 + crLHS425;
    const double crLHS433 = C(0,1)*DN_v(5,0) + C(1,5)*DN_v(5,2) + crLHS188;
    const double crLHS434 = C(0,4)*DN_v(5,0) + crLHS191 + crLHS197;
    const double crLHS435 = DN_v(5,0)*crLHS354;
    const double crLHS436 = C(1,3)*DDN_v[5](1,1);
    const double crLHS437 = C(3,4)*DDN_v[5](1,1);
    const double crLHS438 = C(0,1)*DDN_v[5](1,0) + C(0,3)*DDN_v[5](1,0) + C(0,4)*DDN_v[5](1,0) + C(1,5)*DDN_v[5](1,2) + C(3,3)*DDN_v[5](1,1) + C(3,5)*DDN_v[5](1,2) + C(4,5)*DDN_v[5](1,2) + crLHS436 + crLHS437;
    const double crLHS439 = C(1,1)*DN_v(5,1) + C(1,3)*DN_v(5,0) + C(1,4)*DN_v(5,2);
    const double crLHS440 = C(1,4)*DN_v(5,1);
    const double crLHS441 = C(3,4)*DN_v(5,0) + C(4,4)*DN_v(5,2) + crLHS440;
    const double crLHS442 = DN_v(5,1)*crLHS354;
    const double crLHS443 = 1.0*DDN_v[5](1,0);
    const double crLHS444 = C(1,4)*DDN_v[5](1,1);
    const double crLHS445 = 1.0*DDN_v[5](1,2);
    const double crLHS446 = C(1,3)*crLHS443 + C(1,4)*crLHS445 + C(3,3)*crLHS443 + C(3,4)*crLHS443 + C(3,4)*crLHS445 + C(4,4)*crLHS445 + DDN_v[5](1,1)*crLHS346 + crLHS183 + 1.0*crLHS436 + 1.0*crLHS444;
    const double crLHS447 = C(1,2)*DN_v(5,2) + C(1,5)*DN_v(5,0) + crLHS440;
    const double crLHS448 = C(2,4)*DN_v(5,2);
    const double crLHS449 = C(4,4)*DN_v(5,1) + C(4,5)*DN_v(5,0) + crLHS448;
    const double crLHS450 = DN_v(5,2)*crLHS354;
    const double crLHS451 = C(1,2)*DDN_v[5](1,2) + C(1,5)*DDN_v[5](1,0) + C(2,3)*DDN_v[5](1,2) + C(2,4)*DDN_v[5](1,2) + C(3,5)*DDN_v[5](1,0) + C(4,4)*DDN_v[5](1,1) + C(4,5)*DDN_v[5](1,0) + crLHS437 + crLHS444;
    const double crLHS452 = C(0,1)*DN_v(6,0) + C(1,5)*DN_v(6,2) + crLHS220;
    const double crLHS453 = C(0,4)*DN_v(6,0) + crLHS223 + crLHS229;
    const double crLHS454 = DN_v(6,0)*crLHS354;
    const double crLHS455 = C(1,3)*DDN_v[6](1,1);
    const double crLHS456 = C(3,4)*DDN_v[6](1,1);
    const double crLHS457 = C(0,1)*DDN_v[6](1,0) + C(0,3)*DDN_v[6](1,0) + C(0,4)*DDN_v[6](1,0) + C(1,5)*DDN_v[6](1,2) + C(3,3)*DDN_v[6](1,1) + C(3,5)*DDN_v[6](1,2) + C(4,5)*DDN_v[6](1,2) + crLHS455 + crLHS456;
    const double crLHS458 = C(1,1)*DN_v(6,1) + C(1,3)*DN_v(6,0) + C(1,4)*DN_v(6,2);
    const double crLHS459 = C(1,4)*DN_v(6,1);
    const double crLHS460 = C(3,4)*DN_v(6,0) + C(4,4)*DN_v(6,2) + crLHS459;
    const double crLHS461 = DN_v(6,1)*crLHS354;
    const double crLHS462 = 1.0*DDN_v[6](1,0);
    const double crLHS463 = C(1,4)*DDN_v[6](1,1);
    const double crLHS464 = 1.0*DDN_v[6](1,2);
    const double crLHS465 = C(1,3)*crLHS462 + C(1,4)*crLHS464 + C(3,3)*crLHS462 + C(3,4)*crLHS462 + C(3,4)*crLHS464 + C(4,4)*crLHS464 + DDN_v[6](1,1)*crLHS346 + crLHS215 + 1.0*crLHS455 + 1.0*crLHS463;
    const double crLHS466 = C(1,2)*DN_v(6,2) + C(1,5)*DN_v(6,0) + crLHS459;
    const double crLHS467 = C(2,4)*DN_v(6,2);
    const double crLHS468 = C(4,4)*DN_v(6,1) + C(4,5)*DN_v(6,0) + crLHS467;
    const double crLHS469 = DN_v(6,2)*crLHS354;
    const double crLHS470 = C(1,2)*DDN_v[6](1,2) + C(1,5)*DDN_v[6](1,0) + C(2,3)*DDN_v[6](1,2) + C(2,4)*DDN_v[6](1,2) + C(3,5)*DDN_v[6](1,0) + C(4,4)*DDN_v[6](1,1) + C(4,5)*DDN_v[6](1,0) + crLHS456 + crLHS463;
    const double crLHS471 = C(0,1)*DN_v(7,0) + C(1,5)*DN_v(7,2) + crLHS252;
    const double crLHS472 = C(0,4)*DN_v(7,0) + crLHS255 + crLHS261;
    const double crLHS473 = DN_v(7,0)*crLHS354;
    const double crLHS474 = C(1,3)*DDN_v[7](1,1);
    const double crLHS475 = C(3,4)*DDN_v[7](1,1);
    const double crLHS476 = C(0,1)*DDN_v[7](1,0) + C(0,3)*DDN_v[7](1,0) + C(0,4)*DDN_v[7](1,0) + C(1,5)*DDN_v[7](1,2) + C(3,3)*DDN_v[7](1,1) + C(3,5)*DDN_v[7](1,2) + C(4,5)*DDN_v[7](1,2) + crLHS474 + crLHS475;
    const double crLHS477 = C(1,1)*DN_v(7,1) + C(1,3)*DN_v(7,0) + C(1,4)*DN_v(7,2);
    const double crLHS478 = C(1,4)*DN_v(7,1);
    const double crLHS479 = C(3,4)*DN_v(7,0) + C(4,4)*DN_v(7,2) + crLHS478;
    const double crLHS480 = DN_v(7,1)*crLHS354;
    const double crLHS481 = 1.0*DDN_v[7](1,0);
    const double crLHS482 = C(1,4)*DDN_v[7](1,1);
    const double crLHS483 = 1.0*DDN_v[7](1,2);
    const double crLHS484 = C(1,3)*crLHS481 + C(1,4)*crLHS483 + C(3,3)*crLHS481 + C(3,4)*crLHS481 + C(3,4)*crLHS483 + C(4,4)*crLHS483 + DDN_v[7](1,1)*crLHS346 + crLHS247 + 1.0*crLHS474 + 1.0*crLHS482;
    const double crLHS485 = C(1,2)*DN_v(7,2) + C(1,5)*DN_v(7,0) + crLHS478;
    const double crLHS486 = C(2,4)*DN_v(7,2);
    const double crLHS487 = C(4,4)*DN_v(7,1) + C(4,5)*DN_v(7,0) + crLHS486;
    const double crLHS488 = DN_v(7,2)*crLHS354;
    const double crLHS489 = C(1,2)*DDN_v[7](1,2) + C(1,5)*DDN_v[7](1,0) + C(2,3)*DDN_v[7](1,2) + C(2,4)*DDN_v[7](1,2) + C(3,5)*DDN_v[7](1,0) + C(4,4)*DDN_v[7](1,1) + C(4,5)*DDN_v[7](1,0) + crLHS475 + crLHS482;
    const double crLHS490 = C(0,1)*DN_v(8,0) + C(1,5)*DN_v(8,2) + crLHS284;
    const double crLHS491 = C(0,4)*DN_v(8,0) + crLHS287 + crLHS293;
    const double crLHS492 = DN_v(8,0)*crLHS354;
    const double crLHS493 = C(1,3)*DDN_v[8](1,1);
    const double crLHS494 = C(3,4)*DDN_v[8](1,1);
    const double crLHS495 = C(0,1)*DDN_v[8](1,0) + C(0,3)*DDN_v[8](1,0) + C(0,4)*DDN_v[8](1,0) + C(1,5)*DDN_v[8](1,2) + C(3,3)*DDN_v[8](1,1) + C(3,5)*DDN_v[8](1,2) + C(4,5)*DDN_v[8](1,2) + crLHS493 + crLHS494;
    const double crLHS496 = C(1,1)*DN_v(8,1) + C(1,3)*DN_v(8,0) + C(1,4)*DN_v(8,2);
    const double crLHS497 = C(1,4)*DN_v(8,1);
    const double crLHS498 = C(3,4)*DN_v(8,0) + C(4,4)*DN_v(8,2) + crLHS497;
    const double crLHS499 = DN_v(8,1)*crLHS354;
    const double crLHS500 = 1.0*DDN_v[8](1,0);
    const double crLHS501 = C(1,4)*DDN_v[8](1,1);
    const double crLHS502 = 1.0*DDN_v[8](1,2);
    const double crLHS503 = C(1,3)*crLHS500 + C(1,4)*crLHS502 + C(3,3)*crLHS500 + C(3,4)*crLHS500 + C(3,4)*crLHS502 + C(4,4)*crLHS502 + DDN_v[8](1,1)*crLHS346 + crLHS279 + 1.0*crLHS493 + 1.0*crLHS501;
    const double crLHS504 = C(1,2)*DN_v(8,2) + C(1,5)*DN_v(8,0) + crLHS497;
    const double crLHS505 = C(2,4)*DN_v(8,2);
    const double crLHS506 = C(4,4)*DN_v(8,1) + C(4,5)*DN_v(8,0) + crLHS505;
    const double crLHS507 = DN_v(8,2)*crLHS354;
    const double crLHS508 = C(1,2)*DDN_v[8](1,2) + C(1,5)*DDN_v[8](1,0) + C(2,3)*DDN_v[8](1,2) + C(2,4)*DDN_v[8](1,2) + C(3,5)*DDN_v[8](1,0) + C(4,4)*DDN_v[8](1,1) + C(4,5)*DDN_v[8](1,0) + crLHS494 + crLHS501;
    const double crLHS509 = C(0,1)*DN_v(9,0) + C(1,5)*DN_v(9,2) + crLHS315;
    const double crLHS510 = C(0,4)*DN_v(9,0) + crLHS318 + crLHS324;
    const double crLHS511 = DN_v(9,0)*crLHS354;
    const double crLHS512 = C(1,3)*DDN_v[9](1,1);
    const double crLHS513 = C(3,4)*DDN_v[9](1,1);
    const double crLHS514 = C(0,1)*DDN_v[9](1,0) + C(0,3)*DDN_v[9](1,0) + C(0,4)*DDN_v[9](1,0) + C(1,5)*DDN_v[9](1,2) + C(3,3)*DDN_v[9](1,1) + C(3,5)*DDN_v[9](1,2) + C(4,5)*DDN_v[9](1,2) + crLHS512 + crLHS513;
    const double crLHS515 = C(1,1)*DN_v(9,1) + C(1,3)*DN_v(9,0) + C(1,4)*DN_v(9,2);
    const double crLHS516 = C(1,4)*DN_v(9,1);
    const double crLHS517 = C(3,4)*DN_v(9,0) + C(4,4)*DN_v(9,2) + crLHS516;
    const double crLHS518 = DN_v(9,1)*crLHS354;
    const double crLHS519 = 1.0*DDN_v[9](1,0);
    const double crLHS520 = C(1,4)*DDN_v[9](1,1);
    const double crLHS521 = 1.0*DDN_v[9](1,2);
    const double crLHS522 = C(1,3)*crLHS519 + C(1,4)*crLHS521 + C(3,3)*crLHS519 + C(3,4)*crLHS519 + C(3,4)*crLHS521 + C(4,4)*crLHS521 + DDN_v[9](1,1)*crLHS346 + crLHS310 + 1.0*crLHS512 + 1.0*crLHS520;
    const double crLHS523 = C(1,2)*DN_v(9,2) + C(1,5)*DN_v(9,0) + crLHS516;
    const double crLHS524 = C(2,4)*DN_v(9,2);
    const double crLHS525 = C(4,4)*DN_v(9,1) + C(4,5)*DN_v(9,0) + crLHS524;
    const double crLHS526 = DN_v(9,2)*crLHS354;
    const double crLHS527 = C(1,2)*DDN_v[9](1,2) + C(1,5)*DDN_v[9](1,0) + C(2,3)*DDN_v[9](1,2) + C(2,4)*DDN_v[9](1,2) + C(3,5)*DDN_v[9](1,0) + C(4,4)*DDN_v[9](1,1) + C(4,5)*DDN_v[9](1,0) + crLHS513 + crLHS520;
    const double crLHS528 = -DN_v(0,1)*N_p[0];
    const double crLHS529 = DN_p(0,1)*crLHS20;
    const double crLHS530 = -DN_v(0,1)*N_p[1];
    const double crLHS531 = DN_p(1,1)*crLHS20;
    const double crLHS532 = -DN_v(0,1)*N_p[2];
    const double crLHS533 = DN_p(2,1)*crLHS20;
    const double crLHS534 = -DN_v(0,1)*N_p[3];
    const double crLHS535 = DN_p(3,1)*crLHS20;
    const double crLHS536 = C(0,2)*DN_v(0,0) + C(2,3)*DN_v(0,1) + crLHS39;
    const double crLHS537 = C(2,5)*DDN_v[0](2,2);
    const double crLHS538 = C(4,5)*DDN_v[0](2,2);
    const double crLHS539 = C(0,2)*DDN_v[0](2,0) + C(0,4)*DDN_v[0](2,0) + C(0,5)*DDN_v[0](2,0) + C(2,3)*DDN_v[0](2,1) + C(3,4)*DDN_v[0](2,1) + C(3,5)*DDN_v[0](2,1) + C(5,5)*DDN_v[0](2,2) + crLHS537 + crLHS538;
    const double crLHS540 = C(1,2)*DN_v(0,1) + C(2,3)*DN_v(0,0) + crLHS352;
    const double crLHS541 = C(2,4)*DDN_v[0](2,2);
    const double crLHS542 = C(1,2)*DDN_v[0](2,1) + C(1,4)*DDN_v[0](2,1) + C(1,5)*DDN_v[0](2,1) + C(2,3)*DDN_v[0](2,0) + C(3,4)*DDN_v[0](2,0) + C(3,5)*DDN_v[0](2,0) + C(4,4)*DDN_v[0](2,2) + crLHS538 + crLHS541;
    const double crLHS543 = C(2,2)*DN_v(0,2) + C(2,4)*DN_v(0,1) + C(2,5)*DN_v(0,0);
    const double crLHS544 = 1.0*C(2,2);
    const double crLHS545 = 1.0*DDN_v[0](2,1);
    const double crLHS546 = 1.0*DDN_v[0](2,0);
    const double crLHS547 = C(2,4)*crLHS545 + C(2,5)*crLHS546 + C(4,4)*crLHS545 + C(4,5)*crLHS545 + C(4,5)*crLHS546 + C(5,5)*crLHS546 + DDN_v[0](2,2)*crLHS544 + crLHS18 + 1.0*crLHS537 + 1.0*crLHS541;
    const double crLHS548 = C(0,2)*DN_v(1,0) + C(2,3)*DN_v(1,1) + crLHS71;
    const double crLHS549 = DN_v(0,2)*crLHS9;
    const double crLHS550 = DN_v(1,0)*crLHS549;
    const double crLHS551 = C(2,5)*DDN_v[1](2,2);
    const double crLHS552 = C(4,5)*DDN_v[1](2,2);
    const double crLHS553 = C(0,2)*DDN_v[1](2,0) + C(0,4)*DDN_v[1](2,0) + C(0,5)*DDN_v[1](2,0) + C(2,3)*DDN_v[1](2,1) + C(3,4)*DDN_v[1](2,1) + C(3,5)*DDN_v[1](2,1) + C(5,5)*DDN_v[1](2,2) + crLHS551 + crLHS552;
    const double crLHS554 = C(1,2)*DN_v(1,1) + C(2,3)*DN_v(1,0) + crLHS372;
    const double crLHS555 = DN_v(1,1)*crLHS549;
    const double crLHS556 = C(2,4)*DDN_v[1](2,2);
    const double crLHS557 = C(1,2)*DDN_v[1](2,1) + C(1,4)*DDN_v[1](2,1) + C(1,5)*DDN_v[1](2,1) + C(2,3)*DDN_v[1](2,0) + C(3,4)*DDN_v[1](2,0) + C(3,5)*DDN_v[1](2,0) + C(4,4)*DDN_v[1](2,2) + crLHS552 + crLHS556;
    const double crLHS558 = C(2,2)*DN_v(1,2) + C(2,4)*DN_v(1,1) + C(2,5)*DN_v(1,0);
    const double crLHS559 = DN_v(1,2)*crLHS549;
    const double crLHS560 = 1.0*DDN_v[1](2,1);
    const double crLHS561 = 1.0*DDN_v[1](2,0);
    const double crLHS562 = C(2,4)*crLHS560 + C(2,5)*crLHS561 + C(4,4)*crLHS560 + C(4,5)*crLHS560 + C(4,5)*crLHS561 + C(5,5)*crLHS561 + DDN_v[1](2,2)*crLHS544 + crLHS55 + 1.0*crLHS551 + 1.0*crLHS556;
    const double crLHS563 = C(0,2)*DN_v(2,0) + C(2,3)*DN_v(2,1) + crLHS103;
    const double crLHS564 = DN_v(2,0)*crLHS549;
    const double crLHS565 = C(2,5)*DDN_v[2](2,2);
    const double crLHS566 = C(4,5)*DDN_v[2](2,2);
    const double crLHS567 = C(0,2)*DDN_v[2](2,0) + C(0,4)*DDN_v[2](2,0) + C(0,5)*DDN_v[2](2,0) + C(2,3)*DDN_v[2](2,1) + C(3,4)*DDN_v[2](2,1) + C(3,5)*DDN_v[2](2,1) + C(5,5)*DDN_v[2](2,2) + crLHS565 + crLHS566;
    const double crLHS568 = C(1,2)*DN_v(2,1) + C(2,3)*DN_v(2,0) + crLHS391;
    const double crLHS569 = DN_v(2,1)*crLHS549;
    const double crLHS570 = C(2,4)*DDN_v[2](2,2);
    const double crLHS571 = C(1,2)*DDN_v[2](2,1) + C(1,4)*DDN_v[2](2,1) + C(1,5)*DDN_v[2](2,1) + C(2,3)*DDN_v[2](2,0) + C(3,4)*DDN_v[2](2,0) + C(3,5)*DDN_v[2](2,0) + C(4,4)*DDN_v[2](2,2) + crLHS566 + crLHS570;
    const double crLHS572 = C(2,2)*DN_v(2,2) + C(2,4)*DN_v(2,1) + C(2,5)*DN_v(2,0);
    const double crLHS573 = DN_v(2,2)*crLHS549;
    const double crLHS574 = 1.0*DDN_v[2](2,1);
    const double crLHS575 = 1.0*DDN_v[2](2,0);
    const double crLHS576 = C(2,4)*crLHS574 + C(2,5)*crLHS575 + C(4,4)*crLHS574 + C(4,5)*crLHS574 + C(4,5)*crLHS575 + C(5,5)*crLHS575 + DDN_v[2](2,2)*crLHS544 + 1.0*crLHS565 + 1.0*crLHS570 + crLHS87;
    const double crLHS577 = C(0,2)*DN_v(3,0) + C(2,3)*DN_v(3,1) + crLHS135;
    const double crLHS578 = DN_v(3,0)*crLHS549;
    const double crLHS579 = C(2,5)*DDN_v[3](2,2);
    const double crLHS580 = C(4,5)*DDN_v[3](2,2);
    const double crLHS581 = C(0,2)*DDN_v[3](2,0) + C(0,4)*DDN_v[3](2,0) + C(0,5)*DDN_v[3](2,0) + C(2,3)*DDN_v[3](2,1) + C(3,4)*DDN_v[3](2,1) + C(3,5)*DDN_v[3](2,1) + C(5,5)*DDN_v[3](2,2) + crLHS579 + crLHS580;
    const double crLHS582 = C(1,2)*DN_v(3,1) + C(2,3)*DN_v(3,0) + crLHS410;
    const double crLHS583 = DN_v(3,1)*crLHS549;
    const double crLHS584 = C(2,4)*DDN_v[3](2,2);
    const double crLHS585 = C(1,2)*DDN_v[3](2,1) + C(1,4)*DDN_v[3](2,1) + C(1,5)*DDN_v[3](2,1) + C(2,3)*DDN_v[3](2,0) + C(3,4)*DDN_v[3](2,0) + C(3,5)*DDN_v[3](2,0) + C(4,4)*DDN_v[3](2,2) + crLHS580 + crLHS584;
    const double crLHS586 = C(2,2)*DN_v(3,2) + C(2,4)*DN_v(3,1) + C(2,5)*DN_v(3,0);
    const double crLHS587 = DN_v(3,2)*crLHS549;
    const double crLHS588 = 1.0*DDN_v[3](2,1);
    const double crLHS589 = 1.0*DDN_v[3](2,0);
    const double crLHS590 = C(2,4)*crLHS588 + C(2,5)*crLHS589 + C(4,4)*crLHS588 + C(4,5)*crLHS588 + C(4,5)*crLHS589 + C(5,5)*crLHS589 + DDN_v[3](2,2)*crLHS544 + crLHS119 + 1.0*crLHS579 + 1.0*crLHS584;
    const double crLHS591 = C(0,2)*DN_v(4,0) + C(2,3)*DN_v(4,1) + crLHS167;
    const double crLHS592 = DN_v(4,0)*crLHS549;
    const double crLHS593 = C(2,5)*DDN_v[4](2,2);
    const double crLHS594 = C(4,5)*DDN_v[4](2,2);
    const double crLHS595 = C(0,2)*DDN_v[4](2,0) + C(0,4)*DDN_v[4](2,0) + C(0,5)*DDN_v[4](2,0) + C(2,3)*DDN_v[4](2,1) + C(3,4)*DDN_v[4](2,1) + C(3,5)*DDN_v[4](2,1) + C(5,5)*DDN_v[4](2,2) + crLHS593 + crLHS594;
    const double crLHS596 = C(1,2)*DN_v(4,1) + C(2,3)*DN_v(4,0) + crLHS429;
    const double crLHS597 = DN_v(4,1)*crLHS549;
    const double crLHS598 = C(2,4)*DDN_v[4](2,2);
    const double crLHS599 = C(1,2)*DDN_v[4](2,1) + C(1,4)*DDN_v[4](2,1) + C(1,5)*DDN_v[4](2,1) + C(2,3)*DDN_v[4](2,0) + C(3,4)*DDN_v[4](2,0) + C(3,5)*DDN_v[4](2,0) + C(4,4)*DDN_v[4](2,2) + crLHS594 + crLHS598;
    const double crLHS600 = C(2,2)*DN_v(4,2) + C(2,4)*DN_v(4,1) + C(2,5)*DN_v(4,0);
    const double crLHS601 = DN_v(4,2)*crLHS549;
    const double crLHS602 = 1.0*DDN_v[4](2,1);
    const double crLHS603 = 1.0*DDN_v[4](2,0);
    const double crLHS604 = C(2,4)*crLHS602 + C(2,5)*crLHS603 + C(4,4)*crLHS602 + C(4,5)*crLHS602 + C(4,5)*crLHS603 + C(5,5)*crLHS603 + DDN_v[4](2,2)*crLHS544 + crLHS151 + 1.0*crLHS593 + 1.0*crLHS598;
    const double crLHS605 = C(0,2)*DN_v(5,0) + C(2,3)*DN_v(5,1) + crLHS199;
    const double crLHS606 = DN_v(5,0)*crLHS549;
    const double crLHS607 = C(2,5)*DDN_v[5](2,2);
    const double crLHS608 = C(4,5)*DDN_v[5](2,2);
    const double crLHS609 = C(0,2)*DDN_v[5](2,0) + C(0,4)*DDN_v[5](2,0) + C(0,5)*DDN_v[5](2,0) + C(2,3)*DDN_v[5](2,1) + C(3,4)*DDN_v[5](2,1) + C(3,5)*DDN_v[5](2,1) + C(5,5)*DDN_v[5](2,2) + crLHS607 + crLHS608;
    const double crLHS610 = C(1,2)*DN_v(5,1) + C(2,3)*DN_v(5,0) + crLHS448;
    const double crLHS611 = DN_v(5,1)*crLHS549;
    const double crLHS612 = C(2,4)*DDN_v[5](2,2);
    const double crLHS613 = C(1,2)*DDN_v[5](2,1) + C(1,4)*DDN_v[5](2,1) + C(1,5)*DDN_v[5](2,1) + C(2,3)*DDN_v[5](2,0) + C(3,4)*DDN_v[5](2,0) + C(3,5)*DDN_v[5](2,0) + C(4,4)*DDN_v[5](2,2) + crLHS608 + crLHS612;
    const double crLHS614 = C(2,2)*DN_v(5,2) + C(2,4)*DN_v(5,1) + C(2,5)*DN_v(5,0);
    const double crLHS615 = DN_v(5,2)*crLHS549;
    const double crLHS616 = 1.0*DDN_v[5](2,1);
    const double crLHS617 = 1.0*DDN_v[5](2,0);
    const double crLHS618 = C(2,4)*crLHS616 + C(2,5)*crLHS617 + C(4,4)*crLHS616 + C(4,5)*crLHS616 + C(4,5)*crLHS617 + C(5,5)*crLHS617 + DDN_v[5](2,2)*crLHS544 + crLHS183 + 1.0*crLHS607 + 1.0*crLHS612;
    const double crLHS619 = C(0,2)*DN_v(6,0) + C(2,3)*DN_v(6,1) + crLHS231;
    const double crLHS620 = DN_v(6,0)*crLHS549;
    const double crLHS621 = C(2,5)*DDN_v[6](2,2);
    const double crLHS622 = C(4,5)*DDN_v[6](2,2);
    const double crLHS623 = C(0,2)*DDN_v[6](2,0) + C(0,4)*DDN_v[6](2,0) + C(0,5)*DDN_v[6](2,0) + C(2,3)*DDN_v[6](2,1) + C(3,4)*DDN_v[6](2,1) + C(3,5)*DDN_v[6](2,1) + C(5,5)*DDN_v[6](2,2) + crLHS621 + crLHS622;
    const double crLHS624 = C(1,2)*DN_v(6,1) + C(2,3)*DN_v(6,0) + crLHS467;
    const double crLHS625 = DN_v(6,1)*crLHS549;
    const double crLHS626 = C(2,4)*DDN_v[6](2,2);
    const double crLHS627 = C(1,2)*DDN_v[6](2,1) + C(1,4)*DDN_v[6](2,1) + C(1,5)*DDN_v[6](2,1) + C(2,3)*DDN_v[6](2,0) + C(3,4)*DDN_v[6](2,0) + C(3,5)*DDN_v[6](2,0) + C(4,4)*DDN_v[6](2,2) + crLHS622 + crLHS626;
    const double crLHS628 = C(2,2)*DN_v(6,2) + C(2,4)*DN_v(6,1) + C(2,5)*DN_v(6,0);
    const double crLHS629 = DN_v(6,2)*crLHS549;
    const double crLHS630 = 1.0*DDN_v[6](2,1);
    const double crLHS631 = 1.0*DDN_v[6](2,0);
    const double crLHS632 = C(2,4)*crLHS630 + C(2,5)*crLHS631 + C(4,4)*crLHS630 + C(4,5)*crLHS630 + C(4,5)*crLHS631 + C(5,5)*crLHS631 + DDN_v[6](2,2)*crLHS544 + crLHS215 + 1.0*crLHS621 + 1.0*crLHS626;
    const double crLHS633 = C(0,2)*DN_v(7,0) + C(2,3)*DN_v(7,1) + crLHS263;
    const double crLHS634 = DN_v(7,0)*crLHS549;
    const double crLHS635 = C(2,5)*DDN_v[7](2,2);
    const double crLHS636 = C(4,5)*DDN_v[7](2,2);
    const double crLHS637 = C(0,2)*DDN_v[7](2,0) + C(0,4)*DDN_v[7](2,0) + C(0,5)*DDN_v[7](2,0) + C(2,3)*DDN_v[7](2,1) + C(3,4)*DDN_v[7](2,1) + C(3,5)*DDN_v[7](2,1) + C(5,5)*DDN_v[7](2,2) + crLHS635 + crLHS636;
    const double crLHS638 = C(1,2)*DN_v(7,1) + C(2,3)*DN_v(7,0) + crLHS486;
    const double crLHS639 = DN_v(7,1)*crLHS549;
    const double crLHS640 = C(2,4)*DDN_v[7](2,2);
    const double crLHS641 = C(1,2)*DDN_v[7](2,1) + C(1,4)*DDN_v[7](2,1) + C(1,5)*DDN_v[7](2,1) + C(2,3)*DDN_v[7](2,0) + C(3,4)*DDN_v[7](2,0) + C(3,5)*DDN_v[7](2,0) + C(4,4)*DDN_v[7](2,2) + crLHS636 + crLHS640;
    const double crLHS642 = C(2,2)*DN_v(7,2) + C(2,4)*DN_v(7,1) + C(2,5)*DN_v(7,0);
    const double crLHS643 = DN_v(7,2)*crLHS549;
    const double crLHS644 = 1.0*DDN_v[7](2,1);
    const double crLHS645 = 1.0*DDN_v[7](2,0);
    const double crLHS646 = C(2,4)*crLHS644 + C(2,5)*crLHS645 + C(4,4)*crLHS644 + C(4,5)*crLHS644 + C(4,5)*crLHS645 + C(5,5)*crLHS645 + DDN_v[7](2,2)*crLHS544 + crLHS247 + 1.0*crLHS635 + 1.0*crLHS640;
    const double crLHS647 = C(0,2)*DN_v(8,0) + C(2,3)*DN_v(8,1) + crLHS295;
    const double crLHS648 = DN_v(8,0)*crLHS549;
    const double crLHS649 = C(2,5)*DDN_v[8](2,2);
    const double crLHS650 = C(4,5)*DDN_v[8](2,2);
    const double crLHS651 = C(0,2)*DDN_v[8](2,0) + C(0,4)*DDN_v[8](2,0) + C(0,5)*DDN_v[8](2,0) + C(2,3)*DDN_v[8](2,1) + C(3,4)*DDN_v[8](2,1) + C(3,5)*DDN_v[8](2,1) + C(5,5)*DDN_v[8](2,2) + crLHS649 + crLHS650;
    const double crLHS652 = C(1,2)*DN_v(8,1) + C(2,3)*DN_v(8,0) + crLHS505;
    const double crLHS653 = DN_v(8,1)*crLHS549;
    const double crLHS654 = C(2,4)*DDN_v[8](2,2);
    const double crLHS655 = C(1,2)*DDN_v[8](2,1) + C(1,4)*DDN_v[8](2,1) + C(1,5)*DDN_v[8](2,1) + C(2,3)*DDN_v[8](2,0) + C(3,4)*DDN_v[8](2,0) + C(3,5)*DDN_v[8](2,0) + C(4,4)*DDN_v[8](2,2) + crLHS650 + crLHS654;
    const double crLHS656 = C(2,2)*DN_v(8,2) + C(2,4)*DN_v(8,1) + C(2,5)*DN_v(8,0);
    const double crLHS657 = DN_v(8,2)*crLHS549;
    const double crLHS658 = 1.0*DDN_v[8](2,1);
    const double crLHS659 = 1.0*DDN_v[8](2,0);
    const double crLHS660 = C(2,4)*crLHS658 + C(2,5)*crLHS659 + C(4,4)*crLHS658 + C(4,5)*crLHS658 + C(4,5)*crLHS659 + C(5,5)*crLHS659 + DDN_v[8](2,2)*crLHS544 + crLHS279 + 1.0*crLHS649 + 1.0*crLHS654;
    const double crLHS661 = C(0,2)*DN_v(9,0) + C(2,3)*DN_v(9,1) + crLHS326;
    const double crLHS662 = DN_v(9,0)*crLHS549;
    const double crLHS663 = C(2,5)*DDN_v[9](2,2);
    const double crLHS664 = C(4,5)*DDN_v[9](2,2);
    const double crLHS665 = C(0,2)*DDN_v[9](2,0) + C(0,4)*DDN_v[9](2,0) + C(0,5)*DDN_v[9](2,0) + C(2,3)*DDN_v[9](2,1) + C(3,4)*DDN_v[9](2,1) + C(3,5)*DDN_v[9](2,1) + C(5,5)*DDN_v[9](2,2) + crLHS663 + crLHS664;
    const double crLHS666 = C(1,2)*DN_v(9,1) + C(2,3)*DN_v(9,0) + crLHS524;
    const double crLHS667 = DN_v(9,1)*crLHS549;
    const double crLHS668 = C(2,4)*DDN_v[9](2,2);
    const double crLHS669 = C(1,2)*DDN_v[9](2,1) + C(1,4)*DDN_v[9](2,1) + C(1,5)*DDN_v[9](2,1) + C(2,3)*DDN_v[9](2,0) + C(3,4)*DDN_v[9](2,0) + C(3,5)*DDN_v[9](2,0) + C(4,4)*DDN_v[9](2,2) + crLHS664 + crLHS668;
    const double crLHS670 = C(2,2)*DN_v(9,2) + C(2,4)*DN_v(9,1) + C(2,5)*DN_v(9,0);
    const double crLHS671 = DN_v(9,2)*crLHS549;
    const double crLHS672 = 1.0*DDN_v[9](2,1);
    const double crLHS673 = 1.0*DDN_v[9](2,0);
    const double crLHS674 = C(2,4)*crLHS672 + C(2,5)*crLHS673 + C(4,4)*crLHS672 + C(4,5)*crLHS672 + C(4,5)*crLHS673 + C(5,5)*crLHS673 + DDN_v[9](2,2)*crLHS544 + crLHS310 + 1.0*crLHS663 + 1.0*crLHS668;
    const double crLHS675 = -DN_v(0,2)*N_p[0];
    const double crLHS676 = DN_p(0,2)*crLHS20;
    const double crLHS677 = -DN_v(0,2)*N_p[1];
    const double crLHS678 = DN_p(1,2)*crLHS20;
    const double crLHS679 = -DN_v(0,2)*N_p[2];
    const double crLHS680 = DN_p(2,2)*crLHS20;
    const double crLHS681 = -DN_v(0,2)*N_p[3];
    const double crLHS682 = DN_p(3,2)*crLHS20;
    const double crLHS683 = crLHS20*crLHS53;
    const double crLHS684 = N_v[1]*crLHS22;
    const double crLHS685 = crLHS20*crLHS684;
    const double crLHS686 = N_v[1]*crLHS15 + crLHS57;
    const double crLHS687 = std::pow(N_v[1], 2)*crLHS16 + N_v[1]*crLHS53;
    const double crLHS688 = DN_v(1,0)*crLHS9;
    const double crLHS689 = DN_v(1,1)*crLHS688;
    const double crLHS690 = DN_v(1,2)*crLHS688;
    const double crLHS691 = DN_v(2,0)*crLHS688;
    const double crLHS692 = N_v[2]*crLHS54;
    const double crLHS693 = N_v[1]*crLHS85 + crLHS692;
    const double crLHS694 = DN_v(2,1)*crLHS688;
    const double crLHS695 = DN_v(2,2)*crLHS688;
    const double crLHS696 = DN_v(3,0)*crLHS688;
    const double crLHS697 = N_v[3]*crLHS54;
    const double crLHS698 = N_v[1]*crLHS117 + crLHS697;
    const double crLHS699 = DN_v(3,1)*crLHS688;
    const double crLHS700 = DN_v(3,2)*crLHS688;
    const double crLHS701 = DN_v(4,0)*crLHS688;
    const double crLHS702 = N_v[4]*crLHS54;
    const double crLHS703 = N_v[1]*crLHS149 + crLHS702;
    const double crLHS704 = DN_v(4,1)*crLHS688;
    const double crLHS705 = DN_v(4,2)*crLHS688;
    const double crLHS706 = DN_v(5,0)*crLHS688;
    const double crLHS707 = N_v[5]*crLHS54;
    const double crLHS708 = N_v[1]*crLHS181 + crLHS707;
    const double crLHS709 = DN_v(5,1)*crLHS688;
    const double crLHS710 = DN_v(5,2)*crLHS688;
    const double crLHS711 = DN_v(6,0)*crLHS688;
    const double crLHS712 = N_v[6]*crLHS54;
    const double crLHS713 = N_v[1]*crLHS213 + crLHS712;
    const double crLHS714 = DN_v(6,1)*crLHS688;
    const double crLHS715 = DN_v(6,2)*crLHS688;
    const double crLHS716 = DN_v(7,0)*crLHS688;
    const double crLHS717 = N_v[7]*crLHS54;
    const double crLHS718 = N_v[1]*crLHS245 + crLHS717;
    const double crLHS719 = DN_v(7,1)*crLHS688;
    const double crLHS720 = DN_v(7,2)*crLHS688;
    const double crLHS721 = DN_v(8,0)*crLHS688;
    const double crLHS722 = N_v[8]*crLHS54;
    const double crLHS723 = N_v[1]*crLHS277 + crLHS722;
    const double crLHS724 = DN_v(8,1)*crLHS688;
    const double crLHS725 = DN_v(8,2)*crLHS688;
    const double crLHS726 = DN_v(9,0)*crLHS688;
    const double crLHS727 = N_v[9]*crLHS54;
    const double crLHS728 = N_v[1]*crLHS309 + crLHS727;
    const double crLHS729 = DN_v(9,1)*crLHS688;
    const double crLHS730 = DN_v(9,2)*crLHS688;
    const double crLHS731 = -DN_v(1,0)*N_p[0];
    const double crLHS732 = -DN_v(1,0)*N_p[1];
    const double crLHS733 = -DN_v(1,0)*N_p[2];
    const double crLHS734 = -DN_v(1,0)*N_p[3];
    const double crLHS735 = DN_v(1,1)*crLHS9;
    const double crLHS736 = DN_v(1,2)*crLHS735;
    const double crLHS737 = DN_v(2,0)*crLHS735;
    const double crLHS738 = DN_v(2,1)*crLHS735;
    const double crLHS739 = DN_v(2,2)*crLHS735;
    const double crLHS740 = DN_v(3,0)*crLHS735;
    const double crLHS741 = DN_v(3,1)*crLHS735;
    const double crLHS742 = DN_v(3,2)*crLHS735;
    const double crLHS743 = DN_v(4,0)*crLHS735;
    const double crLHS744 = DN_v(4,1)*crLHS735;
    const double crLHS745 = DN_v(4,2)*crLHS735;
    const double crLHS746 = DN_v(5,0)*crLHS735;
    const double crLHS747 = DN_v(5,1)*crLHS735;
    const double crLHS748 = DN_v(5,2)*crLHS735;
    const double crLHS749 = DN_v(6,0)*crLHS735;
    const double crLHS750 = DN_v(6,1)*crLHS735;
    const double crLHS751 = DN_v(6,2)*crLHS735;
    const double crLHS752 = DN_v(7,0)*crLHS735;
    const double crLHS753 = DN_v(7,1)*crLHS735;
    const double crLHS754 = DN_v(7,2)*crLHS735;
    const double crLHS755 = DN_v(8,0)*crLHS735;
    const double crLHS756 = DN_v(8,1)*crLHS735;
    const double crLHS757 = DN_v(8,2)*crLHS735;
    const double crLHS758 = DN_v(9,0)*crLHS735;
    const double crLHS759 = DN_v(9,1)*crLHS735;
    const double crLHS760 = DN_v(9,2)*crLHS735;
    const double crLHS761 = -DN_v(1,1)*N_p[0];
    const double crLHS762 = -DN_v(1,1)*N_p[1];
    const double crLHS763 = -DN_v(1,1)*N_p[2];
    const double crLHS764 = -DN_v(1,1)*N_p[3];
    const double crLHS765 = DN_v(1,2)*crLHS9;
    const double crLHS766 = DN_v(2,0)*crLHS765;
    const double crLHS767 = DN_v(2,1)*crLHS765;
    const double crLHS768 = DN_v(2,2)*crLHS765;
    const double crLHS769 = DN_v(3,0)*crLHS765;
    const double crLHS770 = DN_v(3,1)*crLHS765;
    const double crLHS771 = DN_v(3,2)*crLHS765;
    const double crLHS772 = DN_v(4,0)*crLHS765;
    const double crLHS773 = DN_v(4,1)*crLHS765;
    const double crLHS774 = DN_v(4,2)*crLHS765;
    const double crLHS775 = DN_v(5,0)*crLHS765;
    const double crLHS776 = DN_v(5,1)*crLHS765;
    const double crLHS777 = DN_v(5,2)*crLHS765;
    const double crLHS778 = DN_v(6,0)*crLHS765;
    const double crLHS779 = DN_v(6,1)*crLHS765;
    const double crLHS780 = DN_v(6,2)*crLHS765;
    const double crLHS781 = DN_v(7,0)*crLHS765;
    const double crLHS782 = DN_v(7,1)*crLHS765;
    const double crLHS783 = DN_v(7,2)*crLHS765;
    const double crLHS784 = DN_v(8,0)*crLHS765;
    const double crLHS785 = DN_v(8,1)*crLHS765;
    const double crLHS786 = DN_v(8,2)*crLHS765;
    const double crLHS787 = DN_v(9,0)*crLHS765;
    const double crLHS788 = DN_v(9,1)*crLHS765;
    const double crLHS789 = DN_v(9,2)*crLHS765;
    const double crLHS790 = -DN_v(1,2)*N_p[0];
    const double crLHS791 = -DN_v(1,2)*N_p[1];
    const double crLHS792 = -DN_v(1,2)*N_p[2];
    const double crLHS793 = -DN_v(1,2)*N_p[3];
    const double crLHS794 = crLHS20*crLHS85;
    const double crLHS795 = N_v[2]*crLHS22;
    const double crLHS796 = crLHS20*crLHS795;
    const double crLHS797 = N_v[2]*crLHS15 + crLHS89;
    const double crLHS798 = N_v[2]*crLHS53 + crLHS692;
    const double crLHS799 = std::pow(N_v[2], 2)*crLHS16 + N_v[2]*crLHS85;
    const double crLHS800 = DN_v(2,0)*crLHS9;
    const double crLHS801 = DN_v(2,1)*crLHS800;
    const double crLHS802 = DN_v(2,2)*crLHS800;
    const double crLHS803 = DN_v(3,0)*crLHS800;
    const double crLHS804 = N_v[3]*crLHS86;
    const double crLHS805 = N_v[2]*crLHS117 + crLHS804;
    const double crLHS806 = DN_v(3,1)*crLHS800;
    const double crLHS807 = DN_v(3,2)*crLHS800;
    const double crLHS808 = DN_v(4,0)*crLHS800;
    const double crLHS809 = N_v[4]*crLHS86;
    const double crLHS810 = N_v[2]*crLHS149 + crLHS809;
    const double crLHS811 = DN_v(4,1)*crLHS800;
    const double crLHS812 = DN_v(4,2)*crLHS800;
    const double crLHS813 = DN_v(5,0)*crLHS800;
    const double crLHS814 = N_v[5]*crLHS86;
    const double crLHS815 = N_v[2]*crLHS181 + crLHS814;
    const double crLHS816 = DN_v(5,1)*crLHS800;
    const double crLHS817 = DN_v(5,2)*crLHS800;
    const double crLHS818 = DN_v(6,0)*crLHS800;
    const double crLHS819 = N_v[6]*crLHS86;
    const double crLHS820 = N_v[2]*crLHS213 + crLHS819;
    const double crLHS821 = DN_v(6,1)*crLHS800;
    const double crLHS822 = DN_v(6,2)*crLHS800;
    const double crLHS823 = DN_v(7,0)*crLHS800;
    const double crLHS824 = N_v[7]*crLHS86;
    const double crLHS825 = N_v[2]*crLHS245 + crLHS824;
    const double crLHS826 = DN_v(7,1)*crLHS800;
    const double crLHS827 = DN_v(7,2)*crLHS800;
    const double crLHS828 = DN_v(8,0)*crLHS800;
    const double crLHS829 = N_v[8]*crLHS86;
    const double crLHS830 = N_v[2]*crLHS277 + crLHS829;
    const double crLHS831 = DN_v(8,1)*crLHS800;
    const double crLHS832 = DN_v(8,2)*crLHS800;
    const double crLHS833 = DN_v(9,0)*crLHS800;
    const double crLHS834 = N_v[9]*crLHS86;
    const double crLHS835 = N_v[2]*crLHS309 + crLHS834;
    const double crLHS836 = DN_v(9,1)*crLHS800;
    const double crLHS837 = DN_v(9,2)*crLHS800;
    const double crLHS838 = -DN_v(2,0)*N_p[0];
    const double crLHS839 = -DN_v(2,0)*N_p[1];
    const double crLHS840 = -DN_v(2,0)*N_p[2];
    const double crLHS841 = -DN_v(2,0)*N_p[3];
    const double crLHS842 = DN_v(2,1)*crLHS9;
    const double crLHS843 = DN_v(2,2)*crLHS842;
    const double crLHS844 = DN_v(3,0)*crLHS842;
    const double crLHS845 = DN_v(3,1)*crLHS842;
    const double crLHS846 = DN_v(3,2)*crLHS842;
    const double crLHS847 = DN_v(4,0)*crLHS842;
    const double crLHS848 = DN_v(4,1)*crLHS842;
    const double crLHS849 = DN_v(4,2)*crLHS842;
    const double crLHS850 = DN_v(5,0)*crLHS842;
    const double crLHS851 = DN_v(5,1)*crLHS842;
    const double crLHS852 = DN_v(5,2)*crLHS842;
    const double crLHS853 = DN_v(6,0)*crLHS842;
    const double crLHS854 = DN_v(6,1)*crLHS842;
    const double crLHS855 = DN_v(6,2)*crLHS842;
    const double crLHS856 = DN_v(7,0)*crLHS842;
    const double crLHS857 = DN_v(7,1)*crLHS842;
    const double crLHS858 = DN_v(7,2)*crLHS842;
    const double crLHS859 = DN_v(8,0)*crLHS842;
    const double crLHS860 = DN_v(8,1)*crLHS842;
    const double crLHS861 = DN_v(8,2)*crLHS842;
    const double crLHS862 = DN_v(9,0)*crLHS842;
    const double crLHS863 = DN_v(9,1)*crLHS842;
    const double crLHS864 = DN_v(9,2)*crLHS842;
    const double crLHS865 = -DN_v(2,1)*N_p[0];
    const double crLHS866 = -DN_v(2,1)*N_p[1];
    const double crLHS867 = -DN_v(2,1)*N_p[2];
    const double crLHS868 = -DN_v(2,1)*N_p[3];
    const double crLHS869 = DN_v(2,2)*crLHS9;
    const double crLHS870 = DN_v(3,0)*crLHS869;
    const double crLHS871 = DN_v(3,1)*crLHS869;
    const double crLHS872 = DN_v(3,2)*crLHS869;
    const double crLHS873 = DN_v(4,0)*crLHS869;
    const double crLHS874 = DN_v(4,1)*crLHS869;
    const double crLHS875 = DN_v(4,2)*crLHS869;
    const double crLHS876 = DN_v(5,0)*crLHS869;
    const double crLHS877 = DN_v(5,1)*crLHS869;
    const double crLHS878 = DN_v(5,2)*crLHS869;
    const double crLHS879 = DN_v(6,0)*crLHS869;
    const double crLHS880 = DN_v(6,1)*crLHS869;
    const double crLHS881 = DN_v(6,2)*crLHS869;
    const double crLHS882 = DN_v(7,0)*crLHS869;
    const double crLHS883 = DN_v(7,1)*crLHS869;
    const double crLHS884 = DN_v(7,2)*crLHS869;
    const double crLHS885 = DN_v(8,0)*crLHS869;
    const double crLHS886 = DN_v(8,1)*crLHS869;
    const double crLHS887 = DN_v(8,2)*crLHS869;
    const double crLHS888 = DN_v(9,0)*crLHS869;
    const double crLHS889 = DN_v(9,1)*crLHS869;
    const double crLHS890 = DN_v(9,2)*crLHS869;
    const double crLHS891 = -DN_v(2,2)*N_p[0];
    const double crLHS892 = -DN_v(2,2)*N_p[1];
    const double crLHS893 = -DN_v(2,2)*N_p[2];
    const double crLHS894 = -DN_v(2,2)*N_p[3];
    const double crLHS895 = crLHS117*crLHS20;
    const double crLHS896 = N_v[3]*crLHS22;
    const double crLHS897 = crLHS20*crLHS896;
    const double crLHS898 = N_v[3]*crLHS15 + crLHS121;
    const double crLHS899 = N_v[3]*crLHS53 + crLHS697;
    const double crLHS900 = N_v[3]*crLHS85 + crLHS804;
    const double crLHS901 = std::pow(N_v[3], 2)*crLHS16 + N_v[3]*crLHS117;
    const double crLHS902 = DN_v(3,0)*crLHS9;
    const double crLHS903 = DN_v(3,1)*crLHS902;
    const double crLHS904 = DN_v(3,2)*crLHS902;
    const double crLHS905 = DN_v(4,0)*crLHS902;
    const double crLHS906 = N_v[4]*crLHS118;
    const double crLHS907 = N_v[3]*crLHS149 + crLHS906;
    const double crLHS908 = DN_v(4,1)*crLHS902;
    const double crLHS909 = DN_v(4,2)*crLHS902;
    const double crLHS910 = DN_v(5,0)*crLHS902;
    const double crLHS911 = N_v[5]*crLHS118;
    const double crLHS912 = N_v[3]*crLHS181 + crLHS911;
    const double crLHS913 = DN_v(5,1)*crLHS902;
    const double crLHS914 = DN_v(5,2)*crLHS902;
    const double crLHS915 = DN_v(6,0)*crLHS902;
    const double crLHS916 = N_v[6]*crLHS118;
    const double crLHS917 = N_v[3]*crLHS213 + crLHS916;
    const double crLHS918 = DN_v(6,1)*crLHS902;
    const double crLHS919 = DN_v(6,2)*crLHS902;
    const double crLHS920 = DN_v(7,0)*crLHS902;
    const double crLHS921 = N_v[7]*crLHS118;
    const double crLHS922 = N_v[3]*crLHS245 + crLHS921;
    const double crLHS923 = DN_v(7,1)*crLHS902;
    const double crLHS924 = DN_v(7,2)*crLHS902;
    const double crLHS925 = DN_v(8,0)*crLHS902;
    const double crLHS926 = N_v[8]*crLHS118;
    const double crLHS927 = N_v[3]*crLHS277 + crLHS926;
    const double crLHS928 = DN_v(8,1)*crLHS902;
    const double crLHS929 = DN_v(8,2)*crLHS902;
    const double crLHS930 = DN_v(9,0)*crLHS902;
    const double crLHS931 = N_v[9]*crLHS118;
    const double crLHS932 = N_v[3]*crLHS309 + crLHS931;
    const double crLHS933 = DN_v(9,1)*crLHS902;
    const double crLHS934 = DN_v(9,2)*crLHS902;
    const double crLHS935 = -DN_v(3,0)*N_p[0];
    const double crLHS936 = -DN_v(3,0)*N_p[1];
    const double crLHS937 = -DN_v(3,0)*N_p[2];
    const double crLHS938 = -DN_v(3,0)*N_p[3];
    const double crLHS939 = DN_v(3,1)*crLHS9;
    const double crLHS940 = DN_v(3,2)*crLHS939;
    const double crLHS941 = DN_v(4,0)*crLHS939;
    const double crLHS942 = DN_v(4,1)*crLHS939;
    const double crLHS943 = DN_v(4,2)*crLHS939;
    const double crLHS944 = DN_v(5,0)*crLHS939;
    const double crLHS945 = DN_v(5,1)*crLHS939;
    const double crLHS946 = DN_v(5,2)*crLHS939;
    const double crLHS947 = DN_v(6,0)*crLHS939;
    const double crLHS948 = DN_v(6,1)*crLHS939;
    const double crLHS949 = DN_v(6,2)*crLHS939;
    const double crLHS950 = DN_v(7,0)*crLHS939;
    const double crLHS951 = DN_v(7,1)*crLHS939;
    const double crLHS952 = DN_v(7,2)*crLHS939;
    const double crLHS953 = DN_v(8,0)*crLHS939;
    const double crLHS954 = DN_v(8,1)*crLHS939;
    const double crLHS955 = DN_v(8,2)*crLHS939;
    const double crLHS956 = DN_v(9,0)*crLHS939;
    const double crLHS957 = DN_v(9,1)*crLHS939;
    const double crLHS958 = DN_v(9,2)*crLHS939;
    const double crLHS959 = -DN_v(3,1)*N_p[0];
    const double crLHS960 = -DN_v(3,1)*N_p[1];
    const double crLHS961 = -DN_v(3,1)*N_p[2];
    const double crLHS962 = -DN_v(3,1)*N_p[3];
    const double crLHS963 = DN_v(3,2)*crLHS9;
    const double crLHS964 = DN_v(4,0)*crLHS963;
    const double crLHS965 = DN_v(4,1)*crLHS963;
    const double crLHS966 = DN_v(4,2)*crLHS963;
    const double crLHS967 = DN_v(5,0)*crLHS963;
    const double crLHS968 = DN_v(5,1)*crLHS963;
    const double crLHS969 = DN_v(5,2)*crLHS963;
    const double crLHS970 = DN_v(6,0)*crLHS963;
    const double crLHS971 = DN_v(6,1)*crLHS963;
    const double crLHS972 = DN_v(6,2)*crLHS963;
    const double crLHS973 = DN_v(7,0)*crLHS963;
    const double crLHS974 = DN_v(7,1)*crLHS963;
    const double crLHS975 = DN_v(7,2)*crLHS963;
    const double crLHS976 = DN_v(8,0)*crLHS963;
    const double crLHS977 = DN_v(8,1)*crLHS963;
    const double crLHS978 = DN_v(8,2)*crLHS963;
    const double crLHS979 = DN_v(9,0)*crLHS963;
    const double crLHS980 = DN_v(9,1)*crLHS963;
    const double crLHS981 = DN_v(9,2)*crLHS963;
    const double crLHS982 = -DN_v(3,2)*N_p[0];
    const double crLHS983 = -DN_v(3,2)*N_p[1];
    const double crLHS984 = -DN_v(3,2)*N_p[2];
    const double crLHS985 = -DN_v(3,2)*N_p[3];
    const double crLHS986 = crLHS149*crLHS20;
    const double crLHS987 = N_v[4]*crLHS22;
    const double crLHS988 = crLHS20*crLHS987;
    const double crLHS989 = N_v[4]*crLHS15 + crLHS153;
    const double crLHS990 = N_v[4]*crLHS53 + crLHS702;
    const double crLHS991 = N_v[4]*crLHS85 + crLHS809;
    const double crLHS992 = N_v[4]*crLHS117 + crLHS906;
    const double crLHS993 = std::pow(N_v[4], 2)*crLHS16 + N_v[4]*crLHS149;
    const double crLHS994 = DN_v(4,0)*crLHS9;
    const double crLHS995 = DN_v(4,1)*crLHS994;
    const double crLHS996 = DN_v(4,2)*crLHS994;
    const double crLHS997 = DN_v(5,0)*crLHS994;
    const double crLHS998 = N_v[5]*crLHS150;
    const double crLHS999 = N_v[4]*crLHS181 + crLHS998;
    const double crLHS1000 = DN_v(5,1)*crLHS994;
    const double crLHS1001 = DN_v(5,2)*crLHS994;
    const double crLHS1002 = DN_v(6,0)*crLHS994;
    const double crLHS1003 = N_v[6]*crLHS150;
    const double crLHS1004 = N_v[4]*crLHS213 + crLHS1003;
    const double crLHS1005 = DN_v(6,1)*crLHS994;
    const double crLHS1006 = DN_v(6,2)*crLHS994;
    const double crLHS1007 = DN_v(7,0)*crLHS994;
    const double crLHS1008 = N_v[7]*crLHS150;
    const double crLHS1009 = N_v[4]*crLHS245 + crLHS1008;
    const double crLHS1010 = DN_v(7,1)*crLHS994;
    const double crLHS1011 = DN_v(7,2)*crLHS994;
    const double crLHS1012 = DN_v(8,0)*crLHS994;
    const double crLHS1013 = N_v[8]*crLHS150;
    const double crLHS1014 = N_v[4]*crLHS277 + crLHS1013;
    const double crLHS1015 = DN_v(8,1)*crLHS994;
    const double crLHS1016 = DN_v(8,2)*crLHS994;
    const double crLHS1017 = DN_v(9,0)*crLHS994;
    const double crLHS1018 = N_v[9]*crLHS150;
    const double crLHS1019 = N_v[4]*crLHS309 + crLHS1018;
    const double crLHS1020 = DN_v(9,1)*crLHS994;
    const double crLHS1021 = DN_v(9,2)*crLHS994;
    const double crLHS1022 = -DN_v(4,0)*N_p[0];
    const double crLHS1023 = -DN_v(4,0)*N_p[1];
    const double crLHS1024 = -DN_v(4,0)*N_p[2];
    const double crLHS1025 = -DN_v(4,0)*N_p[3];
    const double crLHS1026 = DN_v(4,1)*crLHS9;
    const double crLHS1027 = DN_v(4,2)*crLHS1026;
    const double crLHS1028 = DN_v(5,0)*crLHS1026;
    const double crLHS1029 = DN_v(5,1)*crLHS1026;
    const double crLHS1030 = DN_v(5,2)*crLHS1026;
    const double crLHS1031 = DN_v(6,0)*crLHS1026;
    const double crLHS1032 = DN_v(6,1)*crLHS1026;
    const double crLHS1033 = DN_v(6,2)*crLHS1026;
    const double crLHS1034 = DN_v(7,0)*crLHS1026;
    const double crLHS1035 = DN_v(7,1)*crLHS1026;
    const double crLHS1036 = DN_v(7,2)*crLHS1026;
    const double crLHS1037 = DN_v(8,0)*crLHS1026;
    const double crLHS1038 = DN_v(8,1)*crLHS1026;
    const double crLHS1039 = DN_v(8,2)*crLHS1026;
    const double crLHS1040 = DN_v(9,0)*crLHS1026;
    const double crLHS1041 = DN_v(9,1)*crLHS1026;
    const double crLHS1042 = DN_v(9,2)*crLHS1026;
    const double crLHS1043 = -DN_v(4,1)*N_p[0];
    const double crLHS1044 = -DN_v(4,1)*N_p[1];
    const double crLHS1045 = -DN_v(4,1)*N_p[2];
    const double crLHS1046 = -DN_v(4,1)*N_p[3];
    const double crLHS1047 = DN_v(4,2)*crLHS9;
    const double crLHS1048 = DN_v(5,0)*crLHS1047;
    const double crLHS1049 = DN_v(5,1)*crLHS1047;
    const double crLHS1050 = DN_v(5,2)*crLHS1047;
    const double crLHS1051 = DN_v(6,0)*crLHS1047;
    const double crLHS1052 = DN_v(6,1)*crLHS1047;
    const double crLHS1053 = DN_v(6,2)*crLHS1047;
    const double crLHS1054 = DN_v(7,0)*crLHS1047;
    const double crLHS1055 = DN_v(7,1)*crLHS1047;
    const double crLHS1056 = DN_v(7,2)*crLHS1047;
    const double crLHS1057 = DN_v(8,0)*crLHS1047;
    const double crLHS1058 = DN_v(8,1)*crLHS1047;
    const double crLHS1059 = DN_v(8,2)*crLHS1047;
    const double crLHS1060 = DN_v(9,0)*crLHS1047;
    const double crLHS1061 = DN_v(9,1)*crLHS1047;
    const double crLHS1062 = DN_v(9,2)*crLHS1047;
    const double crLHS1063 = -DN_v(4,2)*N_p[0];
    const double crLHS1064 = -DN_v(4,2)*N_p[1];
    const double crLHS1065 = -DN_v(4,2)*N_p[2];
    const double crLHS1066 = -DN_v(4,2)*N_p[3];
    const double crLHS1067 = crLHS181*crLHS20;
    const double crLHS1068 = N_v[5]*crLHS22;
    const double crLHS1069 = crLHS1068*crLHS20;
    const double crLHS1070 = N_v[5]*crLHS15 + crLHS185;
    const double crLHS1071 = N_v[5]*crLHS53 + crLHS707;
    const double crLHS1072 = N_v[5]*crLHS85 + crLHS814;
    const double crLHS1073 = N_v[5]*crLHS117 + crLHS911;
    const double crLHS1074 = N_v[5]*crLHS149 + crLHS998;
    const double crLHS1075 = std::pow(N_v[5], 2)*crLHS16 + N_v[5]*crLHS181;
    const double crLHS1076 = DN_v(5,0)*crLHS9;
    const double crLHS1077 = DN_v(5,1)*crLHS1076;
    const double crLHS1078 = DN_v(5,2)*crLHS1076;
    const double crLHS1079 = DN_v(6,0)*crLHS1076;
    const double crLHS1080 = N_v[6]*crLHS182;
    const double crLHS1081 = N_v[5]*crLHS213 + crLHS1080;
    const double crLHS1082 = DN_v(6,1)*crLHS1076;
    const double crLHS1083 = DN_v(6,2)*crLHS1076;
    const double crLHS1084 = DN_v(7,0)*crLHS1076;
    const double crLHS1085 = N_v[7]*crLHS182;
    const double crLHS1086 = N_v[5]*crLHS245 + crLHS1085;
    const double crLHS1087 = DN_v(7,1)*crLHS1076;
    const double crLHS1088 = DN_v(7,2)*crLHS1076;
    const double crLHS1089 = DN_v(8,0)*crLHS1076;
    const double crLHS1090 = N_v[8]*crLHS182;
    const double crLHS1091 = N_v[5]*crLHS277 + crLHS1090;
    const double crLHS1092 = DN_v(8,1)*crLHS1076;
    const double crLHS1093 = DN_v(8,2)*crLHS1076;
    const double crLHS1094 = DN_v(9,0)*crLHS1076;
    const double crLHS1095 = N_v[9]*crLHS182;
    const double crLHS1096 = N_v[5]*crLHS309 + crLHS1095;
    const double crLHS1097 = DN_v(9,1)*crLHS1076;
    const double crLHS1098 = DN_v(9,2)*crLHS1076;
    const double crLHS1099 = -DN_v(5,0)*N_p[0];
    const double crLHS1100 = -DN_v(5,0)*N_p[1];
    const double crLHS1101 = -DN_v(5,0)*N_p[2];
    const double crLHS1102 = -DN_v(5,0)*N_p[3];
    const double crLHS1103 = DN_v(5,1)*crLHS9;
    const double crLHS1104 = DN_v(5,2)*crLHS1103;
    const double crLHS1105 = DN_v(6,0)*crLHS1103;
    const double crLHS1106 = DN_v(6,1)*crLHS1103;
    const double crLHS1107 = DN_v(6,2)*crLHS1103;
    const double crLHS1108 = DN_v(7,0)*crLHS1103;
    const double crLHS1109 = DN_v(7,1)*crLHS1103;
    const double crLHS1110 = DN_v(7,2)*crLHS1103;
    const double crLHS1111 = DN_v(8,0)*crLHS1103;
    const double crLHS1112 = DN_v(8,1)*crLHS1103;
    const double crLHS1113 = DN_v(8,2)*crLHS1103;
    const double crLHS1114 = DN_v(9,0)*crLHS1103;
    const double crLHS1115 = DN_v(9,1)*crLHS1103;
    const double crLHS1116 = DN_v(9,2)*crLHS1103;
    const double crLHS1117 = -DN_v(5,1)*N_p[0];
    const double crLHS1118 = -DN_v(5,1)*N_p[1];
    const double crLHS1119 = -DN_v(5,1)*N_p[2];
    const double crLHS1120 = -DN_v(5,1)*N_p[3];
    const double crLHS1121 = DN_v(5,2)*crLHS9;
    const double crLHS1122 = DN_v(6,0)*crLHS1121;
    const double crLHS1123 = DN_v(6,1)*crLHS1121;
    const double crLHS1124 = DN_v(6,2)*crLHS1121;
    const double crLHS1125 = DN_v(7,0)*crLHS1121;
    const double crLHS1126 = DN_v(7,1)*crLHS1121;
    const double crLHS1127 = DN_v(7,2)*crLHS1121;
    const double crLHS1128 = DN_v(8,0)*crLHS1121;
    const double crLHS1129 = DN_v(8,1)*crLHS1121;
    const double crLHS1130 = DN_v(8,2)*crLHS1121;
    const double crLHS1131 = DN_v(9,0)*crLHS1121;
    const double crLHS1132 = DN_v(9,1)*crLHS1121;
    const double crLHS1133 = DN_v(9,2)*crLHS1121;
    const double crLHS1134 = -DN_v(5,2)*N_p[0];
    const double crLHS1135 = -DN_v(5,2)*N_p[1];
    const double crLHS1136 = -DN_v(5,2)*N_p[2];
    const double crLHS1137 = -DN_v(5,2)*N_p[3];
    const double crLHS1138 = crLHS20*crLHS213;
    const double crLHS1139 = N_v[6]*crLHS22;
    const double crLHS1140 = crLHS1139*crLHS20;
    const double crLHS1141 = N_v[6]*crLHS15 + crLHS217;
    const double crLHS1142 = N_v[6]*crLHS53 + crLHS712;
    const double crLHS1143 = N_v[6]*crLHS85 + crLHS819;
    const double crLHS1144 = N_v[6]*crLHS117 + crLHS916;
    const double crLHS1145 = N_v[6]*crLHS149 + crLHS1003;
    const double crLHS1146 = N_v[6]*crLHS181 + crLHS1080;
    const double crLHS1147 = std::pow(N_v[6], 2)*crLHS16 + N_v[6]*crLHS213;
    const double crLHS1148 = DN_v(6,0)*crLHS9;
    const double crLHS1149 = DN_v(6,1)*crLHS1148;
    const double crLHS1150 = DN_v(6,2)*crLHS1148;
    const double crLHS1151 = DN_v(7,0)*crLHS1148;
    const double crLHS1152 = N_v[7]*crLHS214;
    const double crLHS1153 = N_v[6]*crLHS245 + crLHS1152;
    const double crLHS1154 = DN_v(7,1)*crLHS1148;
    const double crLHS1155 = DN_v(7,2)*crLHS1148;
    const double crLHS1156 = DN_v(8,0)*crLHS1148;
    const double crLHS1157 = N_v[8]*crLHS214;
    const double crLHS1158 = N_v[6]*crLHS277 + crLHS1157;
    const double crLHS1159 = DN_v(8,1)*crLHS1148;
    const double crLHS1160 = DN_v(8,2)*crLHS1148;
    const double crLHS1161 = DN_v(9,0)*crLHS1148;
    const double crLHS1162 = N_v[9]*crLHS214;
    const double crLHS1163 = N_v[6]*crLHS309 + crLHS1162;
    const double crLHS1164 = DN_v(9,1)*crLHS1148;
    const double crLHS1165 = DN_v(9,2)*crLHS1148;
    const double crLHS1166 = -DN_v(6,0)*N_p[0];
    const double crLHS1167 = -DN_v(6,0)*N_p[1];
    const double crLHS1168 = -DN_v(6,0)*N_p[2];
    const double crLHS1169 = -DN_v(6,0)*N_p[3];
    const double crLHS1170 = DN_v(6,1)*crLHS9;
    const double crLHS1171 = DN_v(6,2)*crLHS1170;
    const double crLHS1172 = DN_v(7,0)*crLHS1170;
    const double crLHS1173 = DN_v(7,1)*crLHS1170;
    const double crLHS1174 = DN_v(7,2)*crLHS1170;
    const double crLHS1175 = DN_v(8,0)*crLHS1170;
    const double crLHS1176 = DN_v(8,1)*crLHS1170;
    const double crLHS1177 = DN_v(8,2)*crLHS1170;
    const double crLHS1178 = DN_v(9,0)*crLHS1170;
    const double crLHS1179 = DN_v(9,1)*crLHS1170;
    const double crLHS1180 = DN_v(9,2)*crLHS1170;
    const double crLHS1181 = -DN_v(6,1)*N_p[0];
    const double crLHS1182 = -DN_v(6,1)*N_p[1];
    const double crLHS1183 = -DN_v(6,1)*N_p[2];
    const double crLHS1184 = -DN_v(6,1)*N_p[3];
    const double crLHS1185 = DN_v(6,2)*crLHS9;
    const double crLHS1186 = DN_v(7,0)*crLHS1185;
    const double crLHS1187 = DN_v(7,1)*crLHS1185;
    const double crLHS1188 = DN_v(7,2)*crLHS1185;
    const double crLHS1189 = DN_v(8,0)*crLHS1185;
    const double crLHS1190 = DN_v(8,1)*crLHS1185;
    const double crLHS1191 = DN_v(8,2)*crLHS1185;
    const double crLHS1192 = DN_v(9,0)*crLHS1185;
    const double crLHS1193 = DN_v(9,1)*crLHS1185;
    const double crLHS1194 = DN_v(9,2)*crLHS1185;
    const double crLHS1195 = -DN_v(6,2)*N_p[0];
    const double crLHS1196 = -DN_v(6,2)*N_p[1];
    const double crLHS1197 = -DN_v(6,2)*N_p[2];
    const double crLHS1198 = -DN_v(6,2)*N_p[3];
    const double crLHS1199 = crLHS20*crLHS245;
    const double crLHS1200 = N_v[7]*crLHS22;
    const double crLHS1201 = crLHS1200*crLHS20;
    const double crLHS1202 = N_v[7]*crLHS15 + crLHS249;
    const double crLHS1203 = N_v[7]*crLHS53 + crLHS717;
    const double crLHS1204 = N_v[7]*crLHS85 + crLHS824;
    const double crLHS1205 = N_v[7]*crLHS117 + crLHS921;
    const double crLHS1206 = N_v[7]*crLHS149 + crLHS1008;
    const double crLHS1207 = N_v[7]*crLHS181 + crLHS1085;
    const double crLHS1208 = N_v[7]*crLHS213 + crLHS1152;
    const double crLHS1209 = std::pow(N_v[7], 2)*crLHS16 + N_v[7]*crLHS245;
    const double crLHS1210 = DN_v(7,0)*crLHS9;
    const double crLHS1211 = DN_v(7,1)*crLHS1210;
    const double crLHS1212 = DN_v(7,2)*crLHS1210;
    const double crLHS1213 = DN_v(8,0)*crLHS1210;
    const double crLHS1214 = N_v[8]*crLHS246;
    const double crLHS1215 = N_v[7]*crLHS277 + crLHS1214;
    const double crLHS1216 = DN_v(8,1)*crLHS1210;
    const double crLHS1217 = DN_v(8,2)*crLHS1210;
    const double crLHS1218 = DN_v(9,0)*crLHS1210;
    const double crLHS1219 = N_v[9]*crLHS246;
    const double crLHS1220 = N_v[7]*crLHS309 + crLHS1219;
    const double crLHS1221 = DN_v(9,1)*crLHS1210;
    const double crLHS1222 = DN_v(9,2)*crLHS1210;
    const double crLHS1223 = -DN_v(7,0)*N_p[0];
    const double crLHS1224 = -DN_v(7,0)*N_p[1];
    const double crLHS1225 = -DN_v(7,0)*N_p[2];
    const double crLHS1226 = -DN_v(7,0)*N_p[3];
    const double crLHS1227 = DN_v(7,1)*crLHS9;
    const double crLHS1228 = DN_v(7,2)*crLHS1227;
    const double crLHS1229 = DN_v(8,0)*crLHS1227;
    const double crLHS1230 = DN_v(8,1)*crLHS1227;
    const double crLHS1231 = DN_v(8,2)*crLHS1227;
    const double crLHS1232 = DN_v(9,0)*crLHS1227;
    const double crLHS1233 = DN_v(9,1)*crLHS1227;
    const double crLHS1234 = DN_v(9,2)*crLHS1227;
    const double crLHS1235 = -DN_v(7,1)*N_p[0];
    const double crLHS1236 = -DN_v(7,1)*N_p[1];
    const double crLHS1237 = -DN_v(7,1)*N_p[2];
    const double crLHS1238 = -DN_v(7,1)*N_p[3];
    const double crLHS1239 = DN_v(7,2)*crLHS9;
    const double crLHS1240 = DN_v(8,0)*crLHS1239;
    const double crLHS1241 = DN_v(8,1)*crLHS1239;
    const double crLHS1242 = DN_v(8,2)*crLHS1239;
    const double crLHS1243 = DN_v(9,0)*crLHS1239;
    const double crLHS1244 = DN_v(9,1)*crLHS1239;
    const double crLHS1245 = DN_v(9,2)*crLHS1239;
    const double crLHS1246 = -DN_v(7,2)*N_p[0];
    const double crLHS1247 = -DN_v(7,2)*N_p[1];
    const double crLHS1248 = -DN_v(7,2)*N_p[2];
    const double crLHS1249 = -DN_v(7,2)*N_p[3];
    const double crLHS1250 = crLHS20*crLHS277;
    const double crLHS1251 = N_v[8]*crLHS22;
    const double crLHS1252 = crLHS1251*crLHS20;
    const double crLHS1253 = N_v[8]*crLHS15 + crLHS281;
    const double crLHS1254 = N_v[8]*crLHS53 + crLHS722;
    const double crLHS1255 = N_v[8]*crLHS85 + crLHS829;
    const double crLHS1256 = N_v[8]*crLHS117 + crLHS926;
    const double crLHS1257 = N_v[8]*crLHS149 + crLHS1013;
    const double crLHS1258 = N_v[8]*crLHS181 + crLHS1090;
    const double crLHS1259 = N_v[8]*crLHS213 + crLHS1157;
    const double crLHS1260 = N_v[8]*crLHS245 + crLHS1214;
    const double crLHS1261 = std::pow(N_v[8], 2)*crLHS16 + N_v[8]*crLHS277;
    const double crLHS1262 = DN_v(8,0)*crLHS9;
    const double crLHS1263 = DN_v(8,1)*crLHS1262;
    const double crLHS1264 = DN_v(8,2)*crLHS1262;
    const double crLHS1265 = DN_v(9,0)*crLHS1262;
    const double crLHS1266 = N_v[9]*crLHS278;
    const double crLHS1267 = N_v[8]*crLHS309 + crLHS1266;
    const double crLHS1268 = DN_v(9,1)*crLHS1262;
    const double crLHS1269 = DN_v(9,2)*crLHS1262;
    const double crLHS1270 = -DN_v(8,0)*N_p[0];
    const double crLHS1271 = -DN_v(8,0)*N_p[1];
    const double crLHS1272 = -DN_v(8,0)*N_p[2];
    const double crLHS1273 = -DN_v(8,0)*N_p[3];
    const double crLHS1274 = DN_v(8,1)*crLHS9;
    const double crLHS1275 = DN_v(8,2)*crLHS1274;
    const double crLHS1276 = DN_v(9,0)*crLHS1274;
    const double crLHS1277 = DN_v(9,1)*crLHS1274;
    const double crLHS1278 = DN_v(9,2)*crLHS1274;
    const double crLHS1279 = -DN_v(8,1)*N_p[0];
    const double crLHS1280 = -DN_v(8,1)*N_p[1];
    const double crLHS1281 = -DN_v(8,1)*N_p[2];
    const double crLHS1282 = -DN_v(8,1)*N_p[3];
    const double crLHS1283 = DN_v(8,2)*crLHS9;
    const double crLHS1284 = DN_v(9,0)*crLHS1283;
    const double crLHS1285 = DN_v(9,1)*crLHS1283;
    const double crLHS1286 = DN_v(9,2)*crLHS1283;
    const double crLHS1287 = -DN_v(8,2)*N_p[0];
    const double crLHS1288 = -DN_v(8,2)*N_p[1];
    const double crLHS1289 = -DN_v(8,2)*N_p[2];
    const double crLHS1290 = -DN_v(8,2)*N_p[3];
    const double crLHS1291 = crLHS20*crLHS309;
    const double crLHS1292 = N_v[9]*crLHS22;
    const double crLHS1293 = crLHS1292*crLHS20;
    const double crLHS1294 = N_v[9]*crLHS15 + crLHS312;
    const double crLHS1295 = N_v[9]*crLHS53 + crLHS727;
    const double crLHS1296 = N_v[9]*crLHS85 + crLHS834;
    const double crLHS1297 = N_v[9]*crLHS117 + crLHS931;
    const double crLHS1298 = N_v[9]*crLHS149 + crLHS1018;
    const double crLHS1299 = N_v[9]*crLHS181 + crLHS1095;
    const double crLHS1300 = N_v[9]*crLHS213 + crLHS1162;
    const double crLHS1301 = N_v[9]*crLHS245 + crLHS1219;
    const double crLHS1302 = N_v[9]*crLHS277 + crLHS1266;
    const double crLHS1303 = std::pow(N_v[9], 2)*crLHS16 + N_v[9]*crLHS309;
    const double crLHS1304 = DN_v(9,0)*crLHS9;
    const double crLHS1305 = DN_v(9,1)*crLHS1304;
    const double crLHS1306 = DN_v(9,2)*crLHS1304;
    const double crLHS1307 = -DN_v(9,0)*N_p[0];
    const double crLHS1308 = -DN_v(9,0)*N_p[1];
    const double crLHS1309 = -DN_v(9,0)*N_p[2];
    const double crLHS1310 = -DN_v(9,0)*N_p[3];
    const double crLHS1311 = DN_v(9,1)*DN_v(9,2)*crLHS9;
    const double crLHS1312 = -DN_v(9,1)*N_p[0];
    const double crLHS1313 = -DN_v(9,1)*N_p[1];
    const double crLHS1314 = -DN_v(9,1)*N_p[2];
    const double crLHS1315 = -DN_v(9,1)*N_p[3];
    const double crLHS1316 = -DN_v(9,2)*N_p[0];
    const double crLHS1317 = -DN_v(9,2)*N_p[1];
    const double crLHS1318 = -DN_v(9,2)*N_p[2];
    const double crLHS1319 = -DN_v(9,2)*N_p[3];
    const double crLHS1320 = crLHS20*gauss_weight;
    const double crLHS1321 = crLHS1320*(DN_p(0,0)*DN_p(1,0) + DN_p(0,1)*DN_p(1,1) + DN_p(0,2)*DN_p(1,2));
    const double crLHS1322 = crLHS1320*(DN_p(0,0)*DN_p(2,0) + DN_p(0,1)*DN_p(2,1) + DN_p(0,2)*DN_p(2,2));
    const double crLHS1323 = crLHS1320*(DN_p(0,0)*DN_p(3,0) + DN_p(0,1)*DN_p(3,1) + DN_p(0,2)*DN_p(3,2));
    const double crLHS1324 = crLHS1320*(DN_p(1,0)*DN_p(2,0) + DN_p(1,1)*DN_p(2,1) + DN_p(1,2)*DN_p(2,2));
    const double crLHS1325 = crLHS1320*(DN_p(1,0)*DN_p(3,0) + DN_p(1,1)*DN_p(3,1) + DN_p(1,2)*DN_p(3,2));
    const double crLHS1326 = crLHS1320*(DN_p(2,0)*DN_p(3,0) + DN_p(2,1)*DN_p(3,1) + DN_p(2,2)*DN_p(3,2));
    rLHS(0,0)+=gauss_weight*(std::pow(DN_v(0,0), 2)*crLHS9 + DN_v(0,0)*crLHS0 + DN_v(0,1)*crLHS2 + DN_v(0,2)*crLHS4 - crLHS19*crLHS21 - crLHS19*crLHS24 + crLHS25);
    rLHS(0,1)+=gauss_weight*(DN_v(0,0)*crLHS26 + DN_v(0,1)*crLHS28 + DN_v(0,2)*crLHS31 - crLHS21*crLHS35 - crLHS24*crLHS35 + crLHS33);
    rLHS(0,2)+=gauss_weight*(DN_v(0,0)*crLHS36 + DN_v(0,1)*crLHS38 + DN_v(0,2)*crLHS40 - crLHS21*crLHS42 - crLHS24*crLHS42 + crLHS41);
    rLHS(0,3)+=gauss_weight*(DN_v(0,0)*crLHS43 + DN_v(0,1)*crLHS45 + DN_v(0,2)*crLHS47 - crLHS21*crLHS56 - crLHS24*crLHS56 + crLHS48 + crLHS58);
    rLHS(0,4)+=gauss_weight*(DN_v(0,0)*crLHS59 + DN_v(0,1)*crLHS61 + DN_v(0,2)*crLHS64 - crLHS21*crLHS67 - crLHS24*crLHS67 + crLHS65);
    rLHS(0,5)+=gauss_weight*(DN_v(0,0)*crLHS68 + DN_v(0,1)*crLHS70 + DN_v(0,2)*crLHS72 - crLHS21*crLHS74 - crLHS24*crLHS74 + crLHS73);
    rLHS(0,6)+=gauss_weight*(DN_v(0,0)*crLHS75 + DN_v(0,1)*crLHS77 + DN_v(0,2)*crLHS79 - crLHS21*crLHS88 - crLHS24*crLHS88 + crLHS80 + crLHS90);
    rLHS(0,7)+=gauss_weight*(DN_v(0,0)*crLHS91 + DN_v(0,1)*crLHS93 + DN_v(0,2)*crLHS96 - crLHS21*crLHS99 - crLHS24*crLHS99 + crLHS97);
    rLHS(0,8)+=gauss_weight*(DN_v(0,0)*crLHS100 + DN_v(0,1)*crLHS102 + DN_v(0,2)*crLHS104 + crLHS105 - crLHS106*crLHS21 - crLHS106*crLHS24);
    rLHS(0,9)+=gauss_weight*(DN_v(0,0)*crLHS107 + DN_v(0,1)*crLHS109 + DN_v(0,2)*crLHS111 + crLHS112 - crLHS120*crLHS21 - crLHS120*crLHS24 + crLHS122);
    rLHS(0,10)+=gauss_weight*(DN_v(0,0)*crLHS123 + DN_v(0,1)*crLHS125 + DN_v(0,2)*crLHS128 + crLHS129 - crLHS131*crLHS21 - crLHS131*crLHS24);
    rLHS(0,11)+=gauss_weight*(DN_v(0,0)*crLHS132 + DN_v(0,1)*crLHS134 + DN_v(0,2)*crLHS136 + crLHS137 - crLHS138*crLHS21 - crLHS138*crLHS24);
    rLHS(0,12)+=gauss_weight*(DN_v(0,0)*crLHS139 + DN_v(0,1)*crLHS141 + DN_v(0,2)*crLHS143 + crLHS144 - crLHS152*crLHS21 - crLHS152*crLHS24 + crLHS154);
    rLHS(0,13)+=gauss_weight*(DN_v(0,0)*crLHS155 + DN_v(0,1)*crLHS157 + DN_v(0,2)*crLHS160 + crLHS161 - crLHS163*crLHS21 - crLHS163*crLHS24);
    rLHS(0,14)+=gauss_weight*(DN_v(0,0)*crLHS164 + DN_v(0,1)*crLHS166 + DN_v(0,2)*crLHS168 + crLHS169 - crLHS170*crLHS21 - crLHS170*crLHS24);
    rLHS(0,15)+=gauss_weight*(DN_v(0,0)*crLHS171 + DN_v(0,1)*crLHS173 + DN_v(0,2)*crLHS175 + crLHS176 - crLHS184*crLHS21 - crLHS184*crLHS24 + crLHS186);
    rLHS(0,16)+=gauss_weight*(DN_v(0,0)*crLHS187 + DN_v(0,1)*crLHS189 + DN_v(0,2)*crLHS192 + crLHS193 - crLHS195*crLHS21 - crLHS195*crLHS24);
    rLHS(0,17)+=gauss_weight*(DN_v(0,0)*crLHS196 + DN_v(0,1)*crLHS198 + DN_v(0,2)*crLHS200 + crLHS201 - crLHS202*crLHS21 - crLHS202*crLHS24);
    rLHS(0,18)+=gauss_weight*(DN_v(0,0)*crLHS203 + DN_v(0,1)*crLHS205 + DN_v(0,2)*crLHS207 + crLHS208 - crLHS21*crLHS216 - crLHS216*crLHS24 + crLHS218);
    rLHS(0,19)+=gauss_weight*(DN_v(0,0)*crLHS219 + DN_v(0,1)*crLHS221 + DN_v(0,2)*crLHS224 - crLHS21*crLHS227 + crLHS225 - crLHS227*crLHS24);
    rLHS(0,20)+=gauss_weight*(DN_v(0,0)*crLHS228 + DN_v(0,1)*crLHS230 + DN_v(0,2)*crLHS232 - crLHS21*crLHS234 + crLHS233 - crLHS234*crLHS24);
    rLHS(0,21)+=gauss_weight*(DN_v(0,0)*crLHS235 + DN_v(0,1)*crLHS237 + DN_v(0,2)*crLHS239 - crLHS21*crLHS248 - crLHS24*crLHS248 + crLHS240 + crLHS250);
    rLHS(0,22)+=gauss_weight*(DN_v(0,0)*crLHS251 + DN_v(0,1)*crLHS253 + DN_v(0,2)*crLHS256 - crLHS21*crLHS259 - crLHS24*crLHS259 + crLHS257);
    rLHS(0,23)+=gauss_weight*(DN_v(0,0)*crLHS260 + DN_v(0,1)*crLHS262 + DN_v(0,2)*crLHS264 - crLHS21*crLHS266 - crLHS24*crLHS266 + crLHS265);
    rLHS(0,24)+=gauss_weight*(DN_v(0,0)*crLHS267 + DN_v(0,1)*crLHS269 + DN_v(0,2)*crLHS271 - crLHS21*crLHS280 - crLHS24*crLHS280 + crLHS272 + crLHS282);
    rLHS(0,25)+=gauss_weight*(DN_v(0,0)*crLHS283 + DN_v(0,1)*crLHS285 + DN_v(0,2)*crLHS288 - crLHS21*crLHS291 - crLHS24*crLHS291 + crLHS289);
    rLHS(0,26)+=gauss_weight*(DN_v(0,0)*crLHS292 + DN_v(0,1)*crLHS294 + DN_v(0,2)*crLHS296 - crLHS21*crLHS298 - crLHS24*crLHS298 + crLHS297);
    rLHS(0,27)+=gauss_weight*(DN_v(0,0)*crLHS299 + DN_v(0,1)*crLHS301 + DN_v(0,2)*crLHS303 - crLHS21*crLHS311 - crLHS24*crLHS311 + crLHS304 + crLHS313);
    rLHS(0,28)+=gauss_weight*(DN_v(0,0)*crLHS314 + DN_v(0,1)*crLHS316 + DN_v(0,2)*crLHS319 - crLHS21*crLHS322 - crLHS24*crLHS322 + crLHS320);
    rLHS(0,29)+=gauss_weight*(DN_v(0,0)*crLHS323 + DN_v(0,1)*crLHS325 + DN_v(0,2)*crLHS327 - crLHS21*crLHS329 - crLHS24*crLHS329 + crLHS328);
    rLHS(0,30)+=gauss_weight*(crLHS15*crLHS331 + crLHS23*crLHS331 + crLHS330);
    rLHS(0,31)+=gauss_weight*(crLHS15*crLHS333 + crLHS23*crLHS333 + crLHS332);
    rLHS(0,32)+=gauss_weight*(crLHS15*crLHS335 + crLHS23*crLHS335 + crLHS334);
    rLHS(0,33)+=gauss_weight*(crLHS15*crLHS337 + crLHS23*crLHS337 + crLHS336);
    rLHS(1,0)+=gauss_weight*(DN_v(0,0)*crLHS2 + DN_v(0,1)*crLHS338 + DN_v(0,2)*crLHS339 - crLHS21*crLHS342 - crLHS24*crLHS342 + crLHS33);
    rLHS(1,1)+=gauss_weight*(DN_v(0,0)*crLHS28 + std::pow(DN_v(0,1), 2)*crLHS9 + DN_v(0,1)*crLHS343 + DN_v(0,2)*crLHS345 - crLHS21*crLHS350 - crLHS24*crLHS350 + crLHS25);
    rLHS(1,2)+=gauss_weight*(DN_v(0,0)*crLHS38 + DN_v(0,1)*crLHS351 + DN_v(0,2)*crLHS353 - crLHS21*crLHS356 - crLHS24*crLHS356 + crLHS355);
    rLHS(1,3)+=gauss_weight*(DN_v(0,0)*crLHS45 + DN_v(0,1)*crLHS357 + DN_v(0,2)*crLHS358 - crLHS21*crLHS362 - crLHS24*crLHS362 + crLHS359);
    rLHS(1,4)+=gauss_weight*(DN_v(0,0)*crLHS61 + DN_v(0,1)*crLHS363 + DN_v(0,2)*crLHS365 - crLHS21*crLHS370 - crLHS24*crLHS370 + crLHS366 + crLHS58);
    rLHS(1,5)+=gauss_weight*(DN_v(0,0)*crLHS70 + DN_v(0,1)*crLHS371 + DN_v(0,2)*crLHS373 - crLHS21*crLHS375 - crLHS24*crLHS375 + crLHS374);
    rLHS(1,6)+=gauss_weight*(DN_v(0,0)*crLHS77 + DN_v(0,1)*crLHS376 + DN_v(0,2)*crLHS377 - crLHS21*crLHS381 - crLHS24*crLHS381 + crLHS378);
    rLHS(1,7)+=gauss_weight*(DN_v(0,0)*crLHS93 + DN_v(0,1)*crLHS382 + DN_v(0,2)*crLHS384 - crLHS21*crLHS389 - crLHS24*crLHS389 + crLHS385 + crLHS90);
    rLHS(1,8)+=gauss_weight*(DN_v(0,0)*crLHS102 + DN_v(0,1)*crLHS390 + DN_v(0,2)*crLHS392 - crLHS21*crLHS394 - crLHS24*crLHS394 + crLHS393);
    rLHS(1,9)+=gauss_weight*(DN_v(0,0)*crLHS109 + DN_v(0,1)*crLHS395 + DN_v(0,2)*crLHS396 - crLHS21*crLHS400 - crLHS24*crLHS400 + crLHS397);
    rLHS(1,10)+=gauss_weight*(DN_v(0,0)*crLHS125 + DN_v(0,1)*crLHS401 + DN_v(0,2)*crLHS403 + crLHS122 - crLHS21*crLHS408 - crLHS24*crLHS408 + crLHS404);
    rLHS(1,11)+=gauss_weight*(DN_v(0,0)*crLHS134 + DN_v(0,1)*crLHS409 + DN_v(0,2)*crLHS411 - crLHS21*crLHS413 - crLHS24*crLHS413 + crLHS412);
    rLHS(1,12)+=gauss_weight*(DN_v(0,0)*crLHS141 + DN_v(0,1)*crLHS414 + DN_v(0,2)*crLHS415 - crLHS21*crLHS419 - crLHS24*crLHS419 + crLHS416);
    rLHS(1,13)+=gauss_weight*(DN_v(0,0)*crLHS157 + DN_v(0,1)*crLHS420 + DN_v(0,2)*crLHS422 + crLHS154 - crLHS21*crLHS427 - crLHS24*crLHS427 + crLHS423);
    rLHS(1,14)+=gauss_weight*(DN_v(0,0)*crLHS166 + DN_v(0,1)*crLHS428 + DN_v(0,2)*crLHS430 - crLHS21*crLHS432 - crLHS24*crLHS432 + crLHS431);
    rLHS(1,15)+=gauss_weight*(DN_v(0,0)*crLHS173 + DN_v(0,1)*crLHS433 + DN_v(0,2)*crLHS434 - crLHS21*crLHS438 - crLHS24*crLHS438 + crLHS435);
    rLHS(1,16)+=gauss_weight*(DN_v(0,0)*crLHS189 + DN_v(0,1)*crLHS439 + DN_v(0,2)*crLHS441 + crLHS186 - crLHS21*crLHS446 - crLHS24*crLHS446 + crLHS442);
    rLHS(1,17)+=gauss_weight*(DN_v(0,0)*crLHS198 + DN_v(0,1)*crLHS447 + DN_v(0,2)*crLHS449 - crLHS21*crLHS451 - crLHS24*crLHS451 + crLHS450);
    rLHS(1,18)+=gauss_weight*(DN_v(0,0)*crLHS205 + DN_v(0,1)*crLHS452 + DN_v(0,2)*crLHS453 - crLHS21*crLHS457 - crLHS24*crLHS457 + crLHS454);
    rLHS(1,19)+=gauss_weight*(DN_v(0,0)*crLHS221 + DN_v(0,1)*crLHS458 + DN_v(0,2)*crLHS460 - crLHS21*crLHS465 + crLHS218 - crLHS24*crLHS465 + crLHS461);
    rLHS(1,20)+=gauss_weight*(DN_v(0,0)*crLHS230 + DN_v(0,1)*crLHS466 + DN_v(0,2)*crLHS468 - crLHS21*crLHS470 - crLHS24*crLHS470 + crLHS469);
    rLHS(1,21)+=gauss_weight*(DN_v(0,0)*crLHS237 + DN_v(0,1)*crLHS471 + DN_v(0,2)*crLHS472 - crLHS21*crLHS476 - crLHS24*crLHS476 + crLHS473);
    rLHS(1,22)+=gauss_weight*(DN_v(0,0)*crLHS253 + DN_v(0,1)*crLHS477 + DN_v(0,2)*crLHS479 - crLHS21*crLHS484 - crLHS24*crLHS484 + crLHS250 + crLHS480);
    rLHS(1,23)+=gauss_weight*(DN_v(0,0)*crLHS262 + DN_v(0,1)*crLHS485 + DN_v(0,2)*crLHS487 - crLHS21*crLHS489 - crLHS24*crLHS489 + crLHS488);
    rLHS(1,24)+=gauss_weight*(DN_v(0,0)*crLHS269 + DN_v(0,1)*crLHS490 + DN_v(0,2)*crLHS491 - crLHS21*crLHS495 - crLHS24*crLHS495 + crLHS492);
    rLHS(1,25)+=gauss_weight*(DN_v(0,0)*crLHS285 + DN_v(0,1)*crLHS496 + DN_v(0,2)*crLHS498 - crLHS21*crLHS503 - crLHS24*crLHS503 + crLHS282 + crLHS499);
    rLHS(1,26)+=gauss_weight*(DN_v(0,0)*crLHS294 + DN_v(0,1)*crLHS504 + DN_v(0,2)*crLHS506 - crLHS21*crLHS508 - crLHS24*crLHS508 + crLHS507);
    rLHS(1,27)+=gauss_weight*(DN_v(0,0)*crLHS301 + DN_v(0,1)*crLHS509 + DN_v(0,2)*crLHS510 - crLHS21*crLHS514 - crLHS24*crLHS514 + crLHS511);
    rLHS(1,28)+=gauss_weight*(DN_v(0,0)*crLHS316 + DN_v(0,1)*crLHS515 + DN_v(0,2)*crLHS517 - crLHS21*crLHS522 - crLHS24*crLHS522 + crLHS313 + crLHS518);
    rLHS(1,29)+=gauss_weight*(DN_v(0,0)*crLHS325 + DN_v(0,1)*crLHS523 + DN_v(0,2)*crLHS525 - crLHS21*crLHS527 - crLHS24*crLHS527 + crLHS526);
    rLHS(1,30)+=gauss_weight*(crLHS15*crLHS529 + crLHS23*crLHS529 + crLHS528);
    rLHS(1,31)+=gauss_weight*(crLHS15*crLHS531 + crLHS23*crLHS531 + crLHS530);
    rLHS(1,32)+=gauss_weight*(crLHS15*crLHS533 + crLHS23*crLHS533 + crLHS532);
    rLHS(1,33)+=gauss_weight*(crLHS15*crLHS535 + crLHS23*crLHS535 + crLHS534);
    rLHS(2,0)+=gauss_weight*(DN_v(0,0)*crLHS4 + DN_v(0,1)*crLHS339 + DN_v(0,2)*crLHS536 - crLHS21*crLHS539 - crLHS24*crLHS539 + crLHS41);
    rLHS(2,1)+=gauss_weight*(DN_v(0,0)*crLHS31 + DN_v(0,1)*crLHS345 + DN_v(0,2)*crLHS540 - crLHS21*crLHS542 - crLHS24*crLHS542 + crLHS355);
    rLHS(2,2)+=gauss_weight*(DN_v(0,0)*crLHS40 + DN_v(0,1)*crLHS353 + std::pow(DN_v(0,2), 2)*crLHS9 + DN_v(0,2)*crLHS543 - crLHS21*crLHS547 - crLHS24*crLHS547 + crLHS25);
    rLHS(2,3)+=gauss_weight*(DN_v(0,0)*crLHS47 + DN_v(0,1)*crLHS358 + DN_v(0,2)*crLHS548 - crLHS21*crLHS553 - crLHS24*crLHS553 + crLHS550);
    rLHS(2,4)+=gauss_weight*(DN_v(0,0)*crLHS64 + DN_v(0,1)*crLHS365 + DN_v(0,2)*crLHS554 - crLHS21*crLHS557 - crLHS24*crLHS557 + crLHS555);
    rLHS(2,5)+=gauss_weight*(DN_v(0,0)*crLHS72 + DN_v(0,1)*crLHS373 + DN_v(0,2)*crLHS558 - crLHS21*crLHS562 - crLHS24*crLHS562 + crLHS559 + crLHS58);
    rLHS(2,6)+=gauss_weight*(DN_v(0,0)*crLHS79 + DN_v(0,1)*crLHS377 + DN_v(0,2)*crLHS563 - crLHS21*crLHS567 - crLHS24*crLHS567 + crLHS564);
    rLHS(2,7)+=gauss_weight*(DN_v(0,0)*crLHS96 + DN_v(0,1)*crLHS384 + DN_v(0,2)*crLHS568 - crLHS21*crLHS571 - crLHS24*crLHS571 + crLHS569);
    rLHS(2,8)+=gauss_weight*(DN_v(0,0)*crLHS104 + DN_v(0,1)*crLHS392 + DN_v(0,2)*crLHS572 - crLHS21*crLHS576 - crLHS24*crLHS576 + crLHS573 + crLHS90);
    rLHS(2,9)+=gauss_weight*(DN_v(0,0)*crLHS111 + DN_v(0,1)*crLHS396 + DN_v(0,2)*crLHS577 - crLHS21*crLHS581 - crLHS24*crLHS581 + crLHS578);
    rLHS(2,10)+=gauss_weight*(DN_v(0,0)*crLHS128 + DN_v(0,1)*crLHS403 + DN_v(0,2)*crLHS582 - crLHS21*crLHS585 - crLHS24*crLHS585 + crLHS583);
    rLHS(2,11)+=gauss_weight*(DN_v(0,0)*crLHS136 + DN_v(0,1)*crLHS411 + DN_v(0,2)*crLHS586 + crLHS122 - crLHS21*crLHS590 - crLHS24*crLHS590 + crLHS587);
    rLHS(2,12)+=gauss_weight*(DN_v(0,0)*crLHS143 + DN_v(0,1)*crLHS415 + DN_v(0,2)*crLHS591 - crLHS21*crLHS595 - crLHS24*crLHS595 + crLHS592);
    rLHS(2,13)+=gauss_weight*(DN_v(0,0)*crLHS160 + DN_v(0,1)*crLHS422 + DN_v(0,2)*crLHS596 - crLHS21*crLHS599 - crLHS24*crLHS599 + crLHS597);
    rLHS(2,14)+=gauss_weight*(DN_v(0,0)*crLHS168 + DN_v(0,1)*crLHS430 + DN_v(0,2)*crLHS600 + crLHS154 - crLHS21*crLHS604 - crLHS24*crLHS604 + crLHS601);
    rLHS(2,15)+=gauss_weight*(DN_v(0,0)*crLHS175 + DN_v(0,1)*crLHS434 + DN_v(0,2)*crLHS605 - crLHS21*crLHS609 - crLHS24*crLHS609 + crLHS606);
    rLHS(2,16)+=gauss_weight*(DN_v(0,0)*crLHS192 + DN_v(0,1)*crLHS441 + DN_v(0,2)*crLHS610 - crLHS21*crLHS613 - crLHS24*crLHS613 + crLHS611);
    rLHS(2,17)+=gauss_weight*(DN_v(0,0)*crLHS200 + DN_v(0,1)*crLHS449 + DN_v(0,2)*crLHS614 + crLHS186 - crLHS21*crLHS618 - crLHS24*crLHS618 + crLHS615);
    rLHS(2,18)+=gauss_weight*(DN_v(0,0)*crLHS207 + DN_v(0,1)*crLHS453 + DN_v(0,2)*crLHS619 - crLHS21*crLHS623 - crLHS24*crLHS623 + crLHS620);
    rLHS(2,19)+=gauss_weight*(DN_v(0,0)*crLHS224 + DN_v(0,1)*crLHS460 + DN_v(0,2)*crLHS624 - crLHS21*crLHS627 - crLHS24*crLHS627 + crLHS625);
    rLHS(2,20)+=gauss_weight*(DN_v(0,0)*crLHS232 + DN_v(0,1)*crLHS468 + DN_v(0,2)*crLHS628 - crLHS21*crLHS632 + crLHS218 - crLHS24*crLHS632 + crLHS629);
    rLHS(2,21)+=gauss_weight*(DN_v(0,0)*crLHS239 + DN_v(0,1)*crLHS472 + DN_v(0,2)*crLHS633 - crLHS21*crLHS637 - crLHS24*crLHS637 + crLHS634);
    rLHS(2,22)+=gauss_weight*(DN_v(0,0)*crLHS256 + DN_v(0,1)*crLHS479 + DN_v(0,2)*crLHS638 - crLHS21*crLHS641 - crLHS24*crLHS641 + crLHS639);
    rLHS(2,23)+=gauss_weight*(DN_v(0,0)*crLHS264 + DN_v(0,1)*crLHS487 + DN_v(0,2)*crLHS642 - crLHS21*crLHS646 - crLHS24*crLHS646 + crLHS250 + crLHS643);
    rLHS(2,24)+=gauss_weight*(DN_v(0,0)*crLHS271 + DN_v(0,1)*crLHS491 + DN_v(0,2)*crLHS647 - crLHS21*crLHS651 - crLHS24*crLHS651 + crLHS648);
    rLHS(2,25)+=gauss_weight*(DN_v(0,0)*crLHS288 + DN_v(0,1)*crLHS498 + DN_v(0,2)*crLHS652 - crLHS21*crLHS655 - crLHS24*crLHS655 + crLHS653);
    rLHS(2,26)+=gauss_weight*(DN_v(0,0)*crLHS296 + DN_v(0,1)*crLHS506 + DN_v(0,2)*crLHS656 - crLHS21*crLHS660 - crLHS24*crLHS660 + crLHS282 + crLHS657);
    rLHS(2,27)+=gauss_weight*(DN_v(0,0)*crLHS303 + DN_v(0,1)*crLHS510 + DN_v(0,2)*crLHS661 - crLHS21*crLHS665 - crLHS24*crLHS665 + crLHS662);
    rLHS(2,28)+=gauss_weight*(DN_v(0,0)*crLHS319 + DN_v(0,1)*crLHS517 + DN_v(0,2)*crLHS666 - crLHS21*crLHS669 - crLHS24*crLHS669 + crLHS667);
    rLHS(2,29)+=gauss_weight*(DN_v(0,0)*crLHS327 + DN_v(0,1)*crLHS525 + DN_v(0,2)*crLHS670 - crLHS21*crLHS674 - crLHS24*crLHS674 + crLHS313 + crLHS671);
    rLHS(2,30)+=gauss_weight*(crLHS15*crLHS676 + crLHS23*crLHS676 + crLHS675);
    rLHS(2,31)+=gauss_weight*(crLHS15*crLHS678 + crLHS23*crLHS678 + crLHS677);
    rLHS(2,32)+=gauss_weight*(crLHS15*crLHS680 + crLHS23*crLHS680 + crLHS679);
    rLHS(2,33)+=gauss_weight*(crLHS15*crLHS682 + crLHS23*crLHS682 + crLHS681);
    rLHS(3,0)+=gauss_weight*(DN_v(1,0)*crLHS0 + DN_v(1,1)*crLHS2 + DN_v(1,2)*crLHS4 - crLHS19*crLHS683 - crLHS19*crLHS685 + crLHS48 + crLHS686);
    rLHS(3,1)+=gauss_weight*(DN_v(1,0)*crLHS26 + DN_v(1,1)*crLHS28 + DN_v(1,2)*crLHS31 - crLHS35*crLHS683 - crLHS35*crLHS685 + crLHS359);
    rLHS(3,2)+=gauss_weight*(DN_v(1,0)*crLHS36 + DN_v(1,1)*crLHS38 + DN_v(1,2)*crLHS40 - crLHS42*crLHS683 - crLHS42*crLHS685 + crLHS550);
    rLHS(3,3)+=gauss_weight*(std::pow(DN_v(1,0), 2)*crLHS9 + DN_v(1,0)*crLHS43 + DN_v(1,1)*crLHS45 + DN_v(1,2)*crLHS47 - crLHS56*crLHS683 - crLHS56*crLHS685 + crLHS687);
    rLHS(3,4)+=gauss_weight*(DN_v(1,0)*crLHS59 + DN_v(1,1)*crLHS61 + DN_v(1,2)*crLHS64 - crLHS67*crLHS683 - crLHS67*crLHS685 + crLHS689);
    rLHS(3,5)+=gauss_weight*(DN_v(1,0)*crLHS68 + DN_v(1,1)*crLHS70 + DN_v(1,2)*crLHS72 - crLHS683*crLHS74 - crLHS685*crLHS74 + crLHS690);
    rLHS(3,6)+=gauss_weight*(DN_v(1,0)*crLHS75 + DN_v(1,1)*crLHS77 + DN_v(1,2)*crLHS79 - crLHS683*crLHS88 - crLHS685*crLHS88 + crLHS691 + crLHS693);
    rLHS(3,7)+=gauss_weight*(DN_v(1,0)*crLHS91 + DN_v(1,1)*crLHS93 + DN_v(1,2)*crLHS96 - crLHS683*crLHS99 - crLHS685*crLHS99 + crLHS694);
    rLHS(3,8)+=gauss_weight*(DN_v(1,0)*crLHS100 + DN_v(1,1)*crLHS102 + DN_v(1,2)*crLHS104 - crLHS106*crLHS683 - crLHS106*crLHS685 + crLHS695);
    rLHS(3,9)+=gauss_weight*(DN_v(1,0)*crLHS107 + DN_v(1,1)*crLHS109 + DN_v(1,2)*crLHS111 - crLHS120*crLHS683 - crLHS120*crLHS685 + crLHS696 + crLHS698);
    rLHS(3,10)+=gauss_weight*(DN_v(1,0)*crLHS123 + DN_v(1,1)*crLHS125 + DN_v(1,2)*crLHS128 - crLHS131*crLHS683 - crLHS131*crLHS685 + crLHS699);
    rLHS(3,11)+=gauss_weight*(DN_v(1,0)*crLHS132 + DN_v(1,1)*crLHS134 + DN_v(1,2)*crLHS136 - crLHS138*crLHS683 - crLHS138*crLHS685 + crLHS700);
    rLHS(3,12)+=gauss_weight*(DN_v(1,0)*crLHS139 + DN_v(1,1)*crLHS141 + DN_v(1,2)*crLHS143 - crLHS152*crLHS683 - crLHS152*crLHS685 + crLHS701 + crLHS703);
    rLHS(3,13)+=gauss_weight*(DN_v(1,0)*crLHS155 + DN_v(1,1)*crLHS157 + DN_v(1,2)*crLHS160 - crLHS163*crLHS683 - crLHS163*crLHS685 + crLHS704);
    rLHS(3,14)+=gauss_weight*(DN_v(1,0)*crLHS164 + DN_v(1,1)*crLHS166 + DN_v(1,2)*crLHS168 - crLHS170*crLHS683 - crLHS170*crLHS685 + crLHS705);
    rLHS(3,15)+=gauss_weight*(DN_v(1,0)*crLHS171 + DN_v(1,1)*crLHS173 + DN_v(1,2)*crLHS175 - crLHS184*crLHS683 - crLHS184*crLHS685 + crLHS706 + crLHS708);
    rLHS(3,16)+=gauss_weight*(DN_v(1,0)*crLHS187 + DN_v(1,1)*crLHS189 + DN_v(1,2)*crLHS192 - crLHS195*crLHS683 - crLHS195*crLHS685 + crLHS709);
    rLHS(3,17)+=gauss_weight*(DN_v(1,0)*crLHS196 + DN_v(1,1)*crLHS198 + DN_v(1,2)*crLHS200 - crLHS202*crLHS683 - crLHS202*crLHS685 + crLHS710);
    rLHS(3,18)+=gauss_weight*(DN_v(1,0)*crLHS203 + DN_v(1,1)*crLHS205 + DN_v(1,2)*crLHS207 - crLHS216*crLHS683 - crLHS216*crLHS685 + crLHS711 + crLHS713);
    rLHS(3,19)+=gauss_weight*(DN_v(1,0)*crLHS219 + DN_v(1,1)*crLHS221 + DN_v(1,2)*crLHS224 - crLHS227*crLHS683 - crLHS227*crLHS685 + crLHS714);
    rLHS(3,20)+=gauss_weight*(DN_v(1,0)*crLHS228 + DN_v(1,1)*crLHS230 + DN_v(1,2)*crLHS232 - crLHS234*crLHS683 - crLHS234*crLHS685 + crLHS715);
    rLHS(3,21)+=gauss_weight*(DN_v(1,0)*crLHS235 + DN_v(1,1)*crLHS237 + DN_v(1,2)*crLHS239 - crLHS248*crLHS683 - crLHS248*crLHS685 + crLHS716 + crLHS718);
    rLHS(3,22)+=gauss_weight*(DN_v(1,0)*crLHS251 + DN_v(1,1)*crLHS253 + DN_v(1,2)*crLHS256 - crLHS259*crLHS683 - crLHS259*crLHS685 + crLHS719);
    rLHS(3,23)+=gauss_weight*(DN_v(1,0)*crLHS260 + DN_v(1,1)*crLHS262 + DN_v(1,2)*crLHS264 - crLHS266*crLHS683 - crLHS266*crLHS685 + crLHS720);
    rLHS(3,24)+=gauss_weight*(DN_v(1,0)*crLHS267 + DN_v(1,1)*crLHS269 + DN_v(1,2)*crLHS271 - crLHS280*crLHS683 - crLHS280*crLHS685 + crLHS721 + crLHS723);
    rLHS(3,25)+=gauss_weight*(DN_v(1,0)*crLHS283 + DN_v(1,1)*crLHS285 + DN_v(1,2)*crLHS288 - crLHS291*crLHS683 - crLHS291*crLHS685 + crLHS724);
    rLHS(3,26)+=gauss_weight*(DN_v(1,0)*crLHS292 + DN_v(1,1)*crLHS294 + DN_v(1,2)*crLHS296 - crLHS298*crLHS683 - crLHS298*crLHS685 + crLHS725);
    rLHS(3,27)+=gauss_weight*(DN_v(1,0)*crLHS299 + DN_v(1,1)*crLHS301 + DN_v(1,2)*crLHS303 - crLHS311*crLHS683 - crLHS311*crLHS685 + crLHS726 + crLHS728);
    rLHS(3,28)+=gauss_weight*(DN_v(1,0)*crLHS314 + DN_v(1,1)*crLHS316 + DN_v(1,2)*crLHS319 - crLHS322*crLHS683 - crLHS322*crLHS685 + crLHS729);
    rLHS(3,29)+=gauss_weight*(DN_v(1,0)*crLHS323 + DN_v(1,1)*crLHS325 + DN_v(1,2)*crLHS327 - crLHS329*crLHS683 - crLHS329*crLHS685 + crLHS730);
    rLHS(3,30)+=gauss_weight*(crLHS331*crLHS53 + crLHS331*crLHS684 + crLHS731);
    rLHS(3,31)+=gauss_weight*(crLHS333*crLHS53 + crLHS333*crLHS684 + crLHS732);
    rLHS(3,32)+=gauss_weight*(crLHS335*crLHS53 + crLHS335*crLHS684 + crLHS733);
    rLHS(3,33)+=gauss_weight*(crLHS337*crLHS53 + crLHS337*crLHS684 + crLHS734);
    rLHS(4,0)+=gauss_weight*(DN_v(1,0)*crLHS2 + DN_v(1,1)*crLHS338 + DN_v(1,2)*crLHS339 - crLHS342*crLHS683 - crLHS342*crLHS685 + crLHS65);
    rLHS(4,1)+=gauss_weight*(DN_v(1,0)*crLHS28 + DN_v(1,1)*crLHS343 + DN_v(1,2)*crLHS345 - crLHS350*crLHS683 - crLHS350*crLHS685 + crLHS366 + crLHS686);
    rLHS(4,2)+=gauss_weight*(DN_v(1,0)*crLHS38 + DN_v(1,1)*crLHS351 + DN_v(1,2)*crLHS353 - crLHS356*crLHS683 - crLHS356*crLHS685 + crLHS555);
    rLHS(4,3)+=gauss_weight*(DN_v(1,0)*crLHS45 + DN_v(1,1)*crLHS357 + DN_v(1,2)*crLHS358 - crLHS362*crLHS683 - crLHS362*crLHS685 + crLHS689);
    rLHS(4,4)+=gauss_weight*(DN_v(1,0)*crLHS61 + std::pow(DN_v(1,1), 2)*crLHS9 + DN_v(1,1)*crLHS363 + DN_v(1,2)*crLHS365 - crLHS370*crLHS683 - crLHS370*crLHS685 + crLHS687);
    rLHS(4,5)+=gauss_weight*(DN_v(1,0)*crLHS70 + DN_v(1,1)*crLHS371 + DN_v(1,2)*crLHS373 - crLHS375*crLHS683 - crLHS375*crLHS685 + crLHS736);
    rLHS(4,6)+=gauss_weight*(DN_v(1,0)*crLHS77 + DN_v(1,1)*crLHS376 + DN_v(1,2)*crLHS377 - crLHS381*crLHS683 - crLHS381*crLHS685 + crLHS737);
    rLHS(4,7)+=gauss_weight*(DN_v(1,0)*crLHS93 + DN_v(1,1)*crLHS382 + DN_v(1,2)*crLHS384 - crLHS389*crLHS683 - crLHS389*crLHS685 + crLHS693 + crLHS738);
    rLHS(4,8)+=gauss_weight*(DN_v(1,0)*crLHS102 + DN_v(1,1)*crLHS390 + DN_v(1,2)*crLHS392 - crLHS394*crLHS683 - crLHS394*crLHS685 + crLHS739);
    rLHS(4,9)+=gauss_weight*(DN_v(1,0)*crLHS109 + DN_v(1,1)*crLHS395 + DN_v(1,2)*crLHS396 - crLHS400*crLHS683 - crLHS400*crLHS685 + crLHS740);
    rLHS(4,10)+=gauss_weight*(DN_v(1,0)*crLHS125 + DN_v(1,1)*crLHS401 + DN_v(1,2)*crLHS403 - crLHS408*crLHS683 - crLHS408*crLHS685 + crLHS698 + crLHS741);
    rLHS(4,11)+=gauss_weight*(DN_v(1,0)*crLHS134 + DN_v(1,1)*crLHS409 + DN_v(1,2)*crLHS411 - crLHS413*crLHS683 - crLHS413*crLHS685 + crLHS742);
    rLHS(4,12)+=gauss_weight*(DN_v(1,0)*crLHS141 + DN_v(1,1)*crLHS414 + DN_v(1,2)*crLHS415 - crLHS419*crLHS683 - crLHS419*crLHS685 + crLHS743);
    rLHS(4,13)+=gauss_weight*(DN_v(1,0)*crLHS157 + DN_v(1,1)*crLHS420 + DN_v(1,2)*crLHS422 - crLHS427*crLHS683 - crLHS427*crLHS685 + crLHS703 + crLHS744);
    rLHS(4,14)+=gauss_weight*(DN_v(1,0)*crLHS166 + DN_v(1,1)*crLHS428 + DN_v(1,2)*crLHS430 - crLHS432*crLHS683 - crLHS432*crLHS685 + crLHS745);
    rLHS(4,15)+=gauss_weight*(DN_v(1,0)*crLHS173 + DN_v(1,1)*crLHS433 + DN_v(1,2)*crLHS434 - crLHS438*crLHS683 - crLHS438*crLHS685 + crLHS746);
    rLHS(4,16)+=gauss_weight*(DN_v(1,0)*crLHS189 + DN_v(1,1)*crLHS439 + DN_v(1,2)*crLHS441 - crLHS446*crLHS683 - crLHS446*crLHS685 + crLHS708 + crLHS747);
    rLHS(4,17)+=gauss_weight*(DN_v(1,0)*crLHS198 + DN_v(1,1)*crLHS447 + DN_v(1,2)*crLHS449 - crLHS451*crLHS683 - crLHS451*crLHS685 + crLHS748);
    rLHS(4,18)+=gauss_weight*(DN_v(1,0)*crLHS205 + DN_v(1,1)*crLHS452 + DN_v(1,2)*crLHS453 - crLHS457*crLHS683 - crLHS457*crLHS685 + crLHS749);
    rLHS(4,19)+=gauss_weight*(DN_v(1,0)*crLHS221 + DN_v(1,1)*crLHS458 + DN_v(1,2)*crLHS460 - crLHS465*crLHS683 - crLHS465*crLHS685 + crLHS713 + crLHS750);
    rLHS(4,20)+=gauss_weight*(DN_v(1,0)*crLHS230 + DN_v(1,1)*crLHS466 + DN_v(1,2)*crLHS468 - crLHS470*crLHS683 - crLHS470*crLHS685 + crLHS751);
    rLHS(4,21)+=gauss_weight*(DN_v(1,0)*crLHS237 + DN_v(1,1)*crLHS471 + DN_v(1,2)*crLHS472 - crLHS476*crLHS683 - crLHS476*crLHS685 + crLHS752);
    rLHS(4,22)+=gauss_weight*(DN_v(1,0)*crLHS253 + DN_v(1,1)*crLHS477 + DN_v(1,2)*crLHS479 - crLHS484*crLHS683 - crLHS484*crLHS685 + crLHS718 + crLHS753);
    rLHS(4,23)+=gauss_weight*(DN_v(1,0)*crLHS262 + DN_v(1,1)*crLHS485 + DN_v(1,2)*crLHS487 - crLHS489*crLHS683 - crLHS489*crLHS685 + crLHS754);
    rLHS(4,24)+=gauss_weight*(DN_v(1,0)*crLHS269 + DN_v(1,1)*crLHS490 + DN_v(1,2)*crLHS491 - crLHS495*crLHS683 - crLHS495*crLHS685 + crLHS755);
    rLHS(4,25)+=gauss_weight*(DN_v(1,0)*crLHS285 + DN_v(1,1)*crLHS496 + DN_v(1,2)*crLHS498 - crLHS503*crLHS683 - crLHS503*crLHS685 + crLHS723 + crLHS756);
    rLHS(4,26)+=gauss_weight*(DN_v(1,0)*crLHS294 + DN_v(1,1)*crLHS504 + DN_v(1,2)*crLHS506 - crLHS508*crLHS683 - crLHS508*crLHS685 + crLHS757);
    rLHS(4,27)+=gauss_weight*(DN_v(1,0)*crLHS301 + DN_v(1,1)*crLHS509 + DN_v(1,2)*crLHS510 - crLHS514*crLHS683 - crLHS514*crLHS685 + crLHS758);
    rLHS(4,28)+=gauss_weight*(DN_v(1,0)*crLHS316 + DN_v(1,1)*crLHS515 + DN_v(1,2)*crLHS517 - crLHS522*crLHS683 - crLHS522*crLHS685 + crLHS728 + crLHS759);
    rLHS(4,29)+=gauss_weight*(DN_v(1,0)*crLHS325 + DN_v(1,1)*crLHS523 + DN_v(1,2)*crLHS525 - crLHS527*crLHS683 - crLHS527*crLHS685 + crLHS760);
    rLHS(4,30)+=gauss_weight*(crLHS529*crLHS53 + crLHS529*crLHS684 + crLHS761);
    rLHS(4,31)+=gauss_weight*(crLHS53*crLHS531 + crLHS531*crLHS684 + crLHS762);
    rLHS(4,32)+=gauss_weight*(crLHS53*crLHS533 + crLHS533*crLHS684 + crLHS763);
    rLHS(4,33)+=gauss_weight*(crLHS53*crLHS535 + crLHS535*crLHS684 + crLHS764);
    rLHS(5,0)+=gauss_weight*(DN_v(1,0)*crLHS4 + DN_v(1,1)*crLHS339 + DN_v(1,2)*crLHS536 - crLHS539*crLHS683 - crLHS539*crLHS685 + crLHS73);
    rLHS(5,1)+=gauss_weight*(DN_v(1,0)*crLHS31 + DN_v(1,1)*crLHS345 + DN_v(1,2)*crLHS540 + crLHS374 - crLHS542*crLHS683 - crLHS542*crLHS685);
    rLHS(5,2)+=gauss_weight*(DN_v(1,0)*crLHS40 + DN_v(1,1)*crLHS353 + DN_v(1,2)*crLHS543 - crLHS547*crLHS683 - crLHS547*crLHS685 + crLHS559 + crLHS686);
    rLHS(5,3)+=gauss_weight*(DN_v(1,0)*crLHS47 + DN_v(1,1)*crLHS358 + DN_v(1,2)*crLHS548 - crLHS553*crLHS683 - crLHS553*crLHS685 + crLHS690);
    rLHS(5,4)+=gauss_weight*(DN_v(1,0)*crLHS64 + DN_v(1,1)*crLHS365 + DN_v(1,2)*crLHS554 - crLHS557*crLHS683 - crLHS557*crLHS685 + crLHS736);
    rLHS(5,5)+=gauss_weight*(DN_v(1,0)*crLHS72 + DN_v(1,1)*crLHS373 + std::pow(DN_v(1,2), 2)*crLHS9 + DN_v(1,2)*crLHS558 - crLHS562*crLHS683 - crLHS562*crLHS685 + crLHS687);
    rLHS(5,6)+=gauss_weight*(DN_v(1,0)*crLHS79 + DN_v(1,1)*crLHS377 + DN_v(1,2)*crLHS563 - crLHS567*crLHS683 - crLHS567*crLHS685 + crLHS766);
    rLHS(5,7)+=gauss_weight*(DN_v(1,0)*crLHS96 + DN_v(1,1)*crLHS384 + DN_v(1,2)*crLHS568 - crLHS571*crLHS683 - crLHS571*crLHS685 + crLHS767);
    rLHS(5,8)+=gauss_weight*(DN_v(1,0)*crLHS104 + DN_v(1,1)*crLHS392 + DN_v(1,2)*crLHS572 - crLHS576*crLHS683 - crLHS576*crLHS685 + crLHS693 + crLHS768);
    rLHS(5,9)+=gauss_weight*(DN_v(1,0)*crLHS111 + DN_v(1,1)*crLHS396 + DN_v(1,2)*crLHS577 - crLHS581*crLHS683 - crLHS581*crLHS685 + crLHS769);
    rLHS(5,10)+=gauss_weight*(DN_v(1,0)*crLHS128 + DN_v(1,1)*crLHS403 + DN_v(1,2)*crLHS582 - crLHS585*crLHS683 - crLHS585*crLHS685 + crLHS770);
    rLHS(5,11)+=gauss_weight*(DN_v(1,0)*crLHS136 + DN_v(1,1)*crLHS411 + DN_v(1,2)*crLHS586 - crLHS590*crLHS683 - crLHS590*crLHS685 + crLHS698 + crLHS771);
    rLHS(5,12)+=gauss_weight*(DN_v(1,0)*crLHS143 + DN_v(1,1)*crLHS415 + DN_v(1,2)*crLHS591 - crLHS595*crLHS683 - crLHS595*crLHS685 + crLHS772);
    rLHS(5,13)+=gauss_weight*(DN_v(1,0)*crLHS160 + DN_v(1,1)*crLHS422 + DN_v(1,2)*crLHS596 - crLHS599*crLHS683 - crLHS599*crLHS685 + crLHS773);
    rLHS(5,14)+=gauss_weight*(DN_v(1,0)*crLHS168 + DN_v(1,1)*crLHS430 + DN_v(1,2)*crLHS600 - crLHS604*crLHS683 - crLHS604*crLHS685 + crLHS703 + crLHS774);
    rLHS(5,15)+=gauss_weight*(DN_v(1,0)*crLHS175 + DN_v(1,1)*crLHS434 + DN_v(1,2)*crLHS605 - crLHS609*crLHS683 - crLHS609*crLHS685 + crLHS775);
    rLHS(5,16)+=gauss_weight*(DN_v(1,0)*crLHS192 + DN_v(1,1)*crLHS441 + DN_v(1,2)*crLHS610 - crLHS613*crLHS683 - crLHS613*crLHS685 + crLHS776);
    rLHS(5,17)+=gauss_weight*(DN_v(1,0)*crLHS200 + DN_v(1,1)*crLHS449 + DN_v(1,2)*crLHS614 - crLHS618*crLHS683 - crLHS618*crLHS685 + crLHS708 + crLHS777);
    rLHS(5,18)+=gauss_weight*(DN_v(1,0)*crLHS207 + DN_v(1,1)*crLHS453 + DN_v(1,2)*crLHS619 - crLHS623*crLHS683 - crLHS623*crLHS685 + crLHS778);
    rLHS(5,19)+=gauss_weight*(DN_v(1,0)*crLHS224 + DN_v(1,1)*crLHS460 + DN_v(1,2)*crLHS624 - crLHS627*crLHS683 - crLHS627*crLHS685 + crLHS779);
    rLHS(5,20)+=gauss_weight*(DN_v(1,0)*crLHS232 + DN_v(1,1)*crLHS468 + DN_v(1,2)*crLHS628 - crLHS632*crLHS683 - crLHS632*crLHS685 + crLHS713 + crLHS780);
    rLHS(5,21)+=gauss_weight*(DN_v(1,0)*crLHS239 + DN_v(1,1)*crLHS472 + DN_v(1,2)*crLHS633 - crLHS637*crLHS683 - crLHS637*crLHS685 + crLHS781);
    rLHS(5,22)+=gauss_weight*(DN_v(1,0)*crLHS256 + DN_v(1,1)*crLHS479 + DN_v(1,2)*crLHS638 - crLHS641*crLHS683 - crLHS641*crLHS685 + crLHS782);
    rLHS(5,23)+=gauss_weight*(DN_v(1,0)*crLHS264 + DN_v(1,1)*crLHS487 + DN_v(1,2)*crLHS642 - crLHS646*crLHS683 - crLHS646*crLHS685 + crLHS718 + crLHS783);
    rLHS(5,24)+=gauss_weight*(DN_v(1,0)*crLHS271 + DN_v(1,1)*crLHS491 + DN_v(1,2)*crLHS647 - crLHS651*crLHS683 - crLHS651*crLHS685 + crLHS784);
    rLHS(5,25)+=gauss_weight*(DN_v(1,0)*crLHS288 + DN_v(1,1)*crLHS498 + DN_v(1,2)*crLHS652 - crLHS655*crLHS683 - crLHS655*crLHS685 + crLHS785);
    rLHS(5,26)+=gauss_weight*(DN_v(1,0)*crLHS296 + DN_v(1,1)*crLHS506 + DN_v(1,2)*crLHS656 - crLHS660*crLHS683 - crLHS660*crLHS685 + crLHS723 + crLHS786);
    rLHS(5,27)+=gauss_weight*(DN_v(1,0)*crLHS303 + DN_v(1,1)*crLHS510 + DN_v(1,2)*crLHS661 - crLHS665*crLHS683 - crLHS665*crLHS685 + crLHS787);
    rLHS(5,28)+=gauss_weight*(DN_v(1,0)*crLHS319 + DN_v(1,1)*crLHS517 + DN_v(1,2)*crLHS666 - crLHS669*crLHS683 - crLHS669*crLHS685 + crLHS788);
    rLHS(5,29)+=gauss_weight*(DN_v(1,0)*crLHS327 + DN_v(1,1)*crLHS525 + DN_v(1,2)*crLHS670 - crLHS674*crLHS683 - crLHS674*crLHS685 + crLHS728 + crLHS789);
    rLHS(5,30)+=gauss_weight*(crLHS53*crLHS676 + crLHS676*crLHS684 + crLHS790);
    rLHS(5,31)+=gauss_weight*(crLHS53*crLHS678 + crLHS678*crLHS684 + crLHS791);
    rLHS(5,32)+=gauss_weight*(crLHS53*crLHS680 + crLHS680*crLHS684 + crLHS792);
    rLHS(5,33)+=gauss_weight*(crLHS53*crLHS682 + crLHS682*crLHS684 + crLHS793);
    rLHS(6,0)+=gauss_weight*(DN_v(2,0)*crLHS0 + DN_v(2,1)*crLHS2 + DN_v(2,2)*crLHS4 - crLHS19*crLHS794 - crLHS19*crLHS796 + crLHS797 + crLHS80);
    rLHS(6,1)+=gauss_weight*(DN_v(2,0)*crLHS26 + DN_v(2,1)*crLHS28 + DN_v(2,2)*crLHS31 - crLHS35*crLHS794 - crLHS35*crLHS796 + crLHS378);
    rLHS(6,2)+=gauss_weight*(DN_v(2,0)*crLHS36 + DN_v(2,1)*crLHS38 + DN_v(2,2)*crLHS40 - crLHS42*crLHS794 - crLHS42*crLHS796 + crLHS564);
    rLHS(6,3)+=gauss_weight*(DN_v(2,0)*crLHS43 + DN_v(2,1)*crLHS45 + DN_v(2,2)*crLHS47 - crLHS56*crLHS794 - crLHS56*crLHS796 + crLHS691 + crLHS798);
    rLHS(6,4)+=gauss_weight*(DN_v(2,0)*crLHS59 + DN_v(2,1)*crLHS61 + DN_v(2,2)*crLHS64 - crLHS67*crLHS794 - crLHS67*crLHS796 + crLHS737);
    rLHS(6,5)+=gauss_weight*(DN_v(2,0)*crLHS68 + DN_v(2,1)*crLHS70 + DN_v(2,2)*crLHS72 - crLHS74*crLHS794 - crLHS74*crLHS796 + crLHS766);
    rLHS(6,6)+=gauss_weight*(std::pow(DN_v(2,0), 2)*crLHS9 + DN_v(2,0)*crLHS75 + DN_v(2,1)*crLHS77 + DN_v(2,2)*crLHS79 - crLHS794*crLHS88 - crLHS796*crLHS88 + crLHS799);
    rLHS(6,7)+=gauss_weight*(DN_v(2,0)*crLHS91 + DN_v(2,1)*crLHS93 + DN_v(2,2)*crLHS96 - crLHS794*crLHS99 - crLHS796*crLHS99 + crLHS801);
    rLHS(6,8)+=gauss_weight*(DN_v(2,0)*crLHS100 + DN_v(2,1)*crLHS102 + DN_v(2,2)*crLHS104 - crLHS106*crLHS794 - crLHS106*crLHS796 + crLHS802);
    rLHS(6,9)+=gauss_weight*(DN_v(2,0)*crLHS107 + DN_v(2,1)*crLHS109 + DN_v(2,2)*crLHS111 - crLHS120*crLHS794 - crLHS120*crLHS796 + crLHS803 + crLHS805);
    rLHS(6,10)+=gauss_weight*(DN_v(2,0)*crLHS123 + DN_v(2,1)*crLHS125 + DN_v(2,2)*crLHS128 - crLHS131*crLHS794 - crLHS131*crLHS796 + crLHS806);
    rLHS(6,11)+=gauss_weight*(DN_v(2,0)*crLHS132 + DN_v(2,1)*crLHS134 + DN_v(2,2)*crLHS136 - crLHS138*crLHS794 - crLHS138*crLHS796 + crLHS807);
    rLHS(6,12)+=gauss_weight*(DN_v(2,0)*crLHS139 + DN_v(2,1)*crLHS141 + DN_v(2,2)*crLHS143 - crLHS152*crLHS794 - crLHS152*crLHS796 + crLHS808 + crLHS810);
    rLHS(6,13)+=gauss_weight*(DN_v(2,0)*crLHS155 + DN_v(2,1)*crLHS157 + DN_v(2,2)*crLHS160 - crLHS163*crLHS794 - crLHS163*crLHS796 + crLHS811);
    rLHS(6,14)+=gauss_weight*(DN_v(2,0)*crLHS164 + DN_v(2,1)*crLHS166 + DN_v(2,2)*crLHS168 - crLHS170*crLHS794 - crLHS170*crLHS796 + crLHS812);
    rLHS(6,15)+=gauss_weight*(DN_v(2,0)*crLHS171 + DN_v(2,1)*crLHS173 + DN_v(2,2)*crLHS175 - crLHS184*crLHS794 - crLHS184*crLHS796 + crLHS813 + crLHS815);
    rLHS(6,16)+=gauss_weight*(DN_v(2,0)*crLHS187 + DN_v(2,1)*crLHS189 + DN_v(2,2)*crLHS192 - crLHS195*crLHS794 - crLHS195*crLHS796 + crLHS816);
    rLHS(6,17)+=gauss_weight*(DN_v(2,0)*crLHS196 + DN_v(2,1)*crLHS198 + DN_v(2,2)*crLHS200 - crLHS202*crLHS794 - crLHS202*crLHS796 + crLHS817);
    rLHS(6,18)+=gauss_weight*(DN_v(2,0)*crLHS203 + DN_v(2,1)*crLHS205 + DN_v(2,2)*crLHS207 - crLHS216*crLHS794 - crLHS216*crLHS796 + crLHS818 + crLHS820);
    rLHS(6,19)+=gauss_weight*(DN_v(2,0)*crLHS219 + DN_v(2,1)*crLHS221 + DN_v(2,2)*crLHS224 - crLHS227*crLHS794 - crLHS227*crLHS796 + crLHS821);
    rLHS(6,20)+=gauss_weight*(DN_v(2,0)*crLHS228 + DN_v(2,1)*crLHS230 + DN_v(2,2)*crLHS232 - crLHS234*crLHS794 - crLHS234*crLHS796 + crLHS822);
    rLHS(6,21)+=gauss_weight*(DN_v(2,0)*crLHS235 + DN_v(2,1)*crLHS237 + DN_v(2,2)*crLHS239 - crLHS248*crLHS794 - crLHS248*crLHS796 + crLHS823 + crLHS825);
    rLHS(6,22)+=gauss_weight*(DN_v(2,0)*crLHS251 + DN_v(2,1)*crLHS253 + DN_v(2,2)*crLHS256 - crLHS259*crLHS794 - crLHS259*crLHS796 + crLHS826);
    rLHS(6,23)+=gauss_weight*(DN_v(2,0)*crLHS260 + DN_v(2,1)*crLHS262 + DN_v(2,2)*crLHS264 - crLHS266*crLHS794 - crLHS266*crLHS796 + crLHS827);
    rLHS(6,24)+=gauss_weight*(DN_v(2,0)*crLHS267 + DN_v(2,1)*crLHS269 + DN_v(2,2)*crLHS271 - crLHS280*crLHS794 - crLHS280*crLHS796 + crLHS828 + crLHS830);
    rLHS(6,25)+=gauss_weight*(DN_v(2,0)*crLHS283 + DN_v(2,1)*crLHS285 + DN_v(2,2)*crLHS288 - crLHS291*crLHS794 - crLHS291*crLHS796 + crLHS831);
    rLHS(6,26)+=gauss_weight*(DN_v(2,0)*crLHS292 + DN_v(2,1)*crLHS294 + DN_v(2,2)*crLHS296 - crLHS298*crLHS794 - crLHS298*crLHS796 + crLHS832);
    rLHS(6,27)+=gauss_weight*(DN_v(2,0)*crLHS299 + DN_v(2,1)*crLHS301 + DN_v(2,2)*crLHS303 - crLHS311*crLHS794 - crLHS311*crLHS796 + crLHS833 + crLHS835);
    rLHS(6,28)+=gauss_weight*(DN_v(2,0)*crLHS314 + DN_v(2,1)*crLHS316 + DN_v(2,2)*crLHS319 - crLHS322*crLHS794 - crLHS322*crLHS796 + crLHS836);
    rLHS(6,29)+=gauss_weight*(DN_v(2,0)*crLHS323 + DN_v(2,1)*crLHS325 + DN_v(2,2)*crLHS327 - crLHS329*crLHS794 - crLHS329*crLHS796 + crLHS837);
    rLHS(6,30)+=gauss_weight*(crLHS331*crLHS795 + crLHS331*crLHS85 + crLHS838);
    rLHS(6,31)+=gauss_weight*(crLHS333*crLHS795 + crLHS333*crLHS85 + crLHS839);
    rLHS(6,32)+=gauss_weight*(crLHS335*crLHS795 + crLHS335*crLHS85 + crLHS840);
    rLHS(6,33)+=gauss_weight*(crLHS337*crLHS795 + crLHS337*crLHS85 + crLHS841);
    rLHS(7,0)+=gauss_weight*(DN_v(2,0)*crLHS2 + DN_v(2,1)*crLHS338 + DN_v(2,2)*crLHS339 - crLHS342*crLHS794 - crLHS342*crLHS796 + crLHS97);
    rLHS(7,1)+=gauss_weight*(DN_v(2,0)*crLHS28 + DN_v(2,1)*crLHS343 + DN_v(2,2)*crLHS345 - crLHS350*crLHS794 - crLHS350*crLHS796 + crLHS385 + crLHS797);
    rLHS(7,2)+=gauss_weight*(DN_v(2,0)*crLHS38 + DN_v(2,1)*crLHS351 + DN_v(2,2)*crLHS353 - crLHS356*crLHS794 - crLHS356*crLHS796 + crLHS569);
    rLHS(7,3)+=gauss_weight*(DN_v(2,0)*crLHS45 + DN_v(2,1)*crLHS357 + DN_v(2,2)*crLHS358 - crLHS362*crLHS794 - crLHS362*crLHS796 + crLHS694);
    rLHS(7,4)+=gauss_weight*(DN_v(2,0)*crLHS61 + DN_v(2,1)*crLHS363 + DN_v(2,2)*crLHS365 - crLHS370*crLHS794 - crLHS370*crLHS796 + crLHS738 + crLHS798);
    rLHS(7,5)+=gauss_weight*(DN_v(2,0)*crLHS70 + DN_v(2,1)*crLHS371 + DN_v(2,2)*crLHS373 - crLHS375*crLHS794 - crLHS375*crLHS796 + crLHS767);
    rLHS(7,6)+=gauss_weight*(DN_v(2,0)*crLHS77 + DN_v(2,1)*crLHS376 + DN_v(2,2)*crLHS377 - crLHS381*crLHS794 - crLHS381*crLHS796 + crLHS801);
    rLHS(7,7)+=gauss_weight*(DN_v(2,0)*crLHS93 + std::pow(DN_v(2,1), 2)*crLHS9 + DN_v(2,1)*crLHS382 + DN_v(2,2)*crLHS384 - crLHS389*crLHS794 - crLHS389*crLHS796 + crLHS799);
    rLHS(7,8)+=gauss_weight*(DN_v(2,0)*crLHS102 + DN_v(2,1)*crLHS390 + DN_v(2,2)*crLHS392 - crLHS394*crLHS794 - crLHS394*crLHS796 + crLHS843);
    rLHS(7,9)+=gauss_weight*(DN_v(2,0)*crLHS109 + DN_v(2,1)*crLHS395 + DN_v(2,2)*crLHS396 - crLHS400*crLHS794 - crLHS400*crLHS796 + crLHS844);
    rLHS(7,10)+=gauss_weight*(DN_v(2,0)*crLHS125 + DN_v(2,1)*crLHS401 + DN_v(2,2)*crLHS403 - crLHS408*crLHS794 - crLHS408*crLHS796 + crLHS805 + crLHS845);
    rLHS(7,11)+=gauss_weight*(DN_v(2,0)*crLHS134 + DN_v(2,1)*crLHS409 + DN_v(2,2)*crLHS411 - crLHS413*crLHS794 - crLHS413*crLHS796 + crLHS846);
    rLHS(7,12)+=gauss_weight*(DN_v(2,0)*crLHS141 + DN_v(2,1)*crLHS414 + DN_v(2,2)*crLHS415 - crLHS419*crLHS794 - crLHS419*crLHS796 + crLHS847);
    rLHS(7,13)+=gauss_weight*(DN_v(2,0)*crLHS157 + DN_v(2,1)*crLHS420 + DN_v(2,2)*crLHS422 - crLHS427*crLHS794 - crLHS427*crLHS796 + crLHS810 + crLHS848);
    rLHS(7,14)+=gauss_weight*(DN_v(2,0)*crLHS166 + DN_v(2,1)*crLHS428 + DN_v(2,2)*crLHS430 - crLHS432*crLHS794 - crLHS432*crLHS796 + crLHS849);
    rLHS(7,15)+=gauss_weight*(DN_v(2,0)*crLHS173 + DN_v(2,1)*crLHS433 + DN_v(2,2)*crLHS434 - crLHS438*crLHS794 - crLHS438*crLHS796 + crLHS850);
    rLHS(7,16)+=gauss_weight*(DN_v(2,0)*crLHS189 + DN_v(2,1)*crLHS439 + DN_v(2,2)*crLHS441 - crLHS446*crLHS794 - crLHS446*crLHS796 + crLHS815 + crLHS851);
    rLHS(7,17)+=gauss_weight*(DN_v(2,0)*crLHS198 + DN_v(2,1)*crLHS447 + DN_v(2,2)*crLHS449 - crLHS451*crLHS794 - crLHS451*crLHS796 + crLHS852);
    rLHS(7,18)+=gauss_weight*(DN_v(2,0)*crLHS205 + DN_v(2,1)*crLHS452 + DN_v(2,2)*crLHS453 - crLHS457*crLHS794 - crLHS457*crLHS796 + crLHS853);
    rLHS(7,19)+=gauss_weight*(DN_v(2,0)*crLHS221 + DN_v(2,1)*crLHS458 + DN_v(2,2)*crLHS460 - crLHS465*crLHS794 - crLHS465*crLHS796 + crLHS820 + crLHS854);
    rLHS(7,20)+=gauss_weight*(DN_v(2,0)*crLHS230 + DN_v(2,1)*crLHS466 + DN_v(2,2)*crLHS468 - crLHS470*crLHS794 - crLHS470*crLHS796 + crLHS855);
    rLHS(7,21)+=gauss_weight*(DN_v(2,0)*crLHS237 + DN_v(2,1)*crLHS471 + DN_v(2,2)*crLHS472 - crLHS476*crLHS794 - crLHS476*crLHS796 + crLHS856);
    rLHS(7,22)+=gauss_weight*(DN_v(2,0)*crLHS253 + DN_v(2,1)*crLHS477 + DN_v(2,2)*crLHS479 - crLHS484*crLHS794 - crLHS484*crLHS796 + crLHS825 + crLHS857);
    rLHS(7,23)+=gauss_weight*(DN_v(2,0)*crLHS262 + DN_v(2,1)*crLHS485 + DN_v(2,2)*crLHS487 - crLHS489*crLHS794 - crLHS489*crLHS796 + crLHS858);
    rLHS(7,24)+=gauss_weight*(DN_v(2,0)*crLHS269 + DN_v(2,1)*crLHS490 + DN_v(2,2)*crLHS491 - crLHS495*crLHS794 - crLHS495*crLHS796 + crLHS859);
    rLHS(7,25)+=gauss_weight*(DN_v(2,0)*crLHS285 + DN_v(2,1)*crLHS496 + DN_v(2,2)*crLHS498 - crLHS503*crLHS794 - crLHS503*crLHS796 + crLHS830 + crLHS860);
    rLHS(7,26)+=gauss_weight*(DN_v(2,0)*crLHS294 + DN_v(2,1)*crLHS504 + DN_v(2,2)*crLHS506 - crLHS508*crLHS794 - crLHS508*crLHS796 + crLHS861);
    rLHS(7,27)+=gauss_weight*(DN_v(2,0)*crLHS301 + DN_v(2,1)*crLHS509 + DN_v(2,2)*crLHS510 - crLHS514*crLHS794 - crLHS514*crLHS796 + crLHS862);
    rLHS(7,28)+=gauss_weight*(DN_v(2,0)*crLHS316 + DN_v(2,1)*crLHS515 + DN_v(2,2)*crLHS517 - crLHS522*crLHS794 - crLHS522*crLHS796 + crLHS835 + crLHS863);
    rLHS(7,29)+=gauss_weight*(DN_v(2,0)*crLHS325 + DN_v(2,1)*crLHS523 + DN_v(2,2)*crLHS525 - crLHS527*crLHS794 - crLHS527*crLHS796 + crLHS864);
    rLHS(7,30)+=gauss_weight*(crLHS529*crLHS795 + crLHS529*crLHS85 + crLHS865);
    rLHS(7,31)+=gauss_weight*(crLHS531*crLHS795 + crLHS531*crLHS85 + crLHS866);
    rLHS(7,32)+=gauss_weight*(crLHS533*crLHS795 + crLHS533*crLHS85 + crLHS867);
    rLHS(7,33)+=gauss_weight*(crLHS535*crLHS795 + crLHS535*crLHS85 + crLHS868);
    rLHS(8,0)+=gauss_weight*(DN_v(2,0)*crLHS4 + DN_v(2,1)*crLHS339 + DN_v(2,2)*crLHS536 + crLHS105 - crLHS539*crLHS794 - crLHS539*crLHS796);
    rLHS(8,1)+=gauss_weight*(DN_v(2,0)*crLHS31 + DN_v(2,1)*crLHS345 + DN_v(2,2)*crLHS540 + crLHS393 - crLHS542*crLHS794 - crLHS542*crLHS796);
    rLHS(8,2)+=gauss_weight*(DN_v(2,0)*crLHS40 + DN_v(2,1)*crLHS353 + DN_v(2,2)*crLHS543 - crLHS547*crLHS794 - crLHS547*crLHS796 + crLHS573 + crLHS797);
    rLHS(8,3)+=gauss_weight*(DN_v(2,0)*crLHS47 + DN_v(2,1)*crLHS358 + DN_v(2,2)*crLHS548 - crLHS553*crLHS794 - crLHS553*crLHS796 + crLHS695);
    rLHS(8,4)+=gauss_weight*(DN_v(2,0)*crLHS64 + DN_v(2,1)*crLHS365 + DN_v(2,2)*crLHS554 - crLHS557*crLHS794 - crLHS557*crLHS796 + crLHS739);
    rLHS(8,5)+=gauss_weight*(DN_v(2,0)*crLHS72 + DN_v(2,1)*crLHS373 + DN_v(2,2)*crLHS558 - crLHS562*crLHS794 - crLHS562*crLHS796 + crLHS768 + crLHS798);
    rLHS(8,6)+=gauss_weight*(DN_v(2,0)*crLHS79 + DN_v(2,1)*crLHS377 + DN_v(2,2)*crLHS563 - crLHS567*crLHS794 - crLHS567*crLHS796 + crLHS802);
    rLHS(8,7)+=gauss_weight*(DN_v(2,0)*crLHS96 + DN_v(2,1)*crLHS384 + DN_v(2,2)*crLHS568 - crLHS571*crLHS794 - crLHS571*crLHS796 + crLHS843);
    rLHS(8,8)+=gauss_weight*(DN_v(2,0)*crLHS104 + DN_v(2,1)*crLHS392 + std::pow(DN_v(2,2), 2)*crLHS9 + DN_v(2,2)*crLHS572 - crLHS576*crLHS794 - crLHS576*crLHS796 + crLHS799);
    rLHS(8,9)+=gauss_weight*(DN_v(2,0)*crLHS111 + DN_v(2,1)*crLHS396 + DN_v(2,2)*crLHS577 - crLHS581*crLHS794 - crLHS581*crLHS796 + crLHS870);
    rLHS(8,10)+=gauss_weight*(DN_v(2,0)*crLHS128 + DN_v(2,1)*crLHS403 + DN_v(2,2)*crLHS582 - crLHS585*crLHS794 - crLHS585*crLHS796 + crLHS871);
    rLHS(8,11)+=gauss_weight*(DN_v(2,0)*crLHS136 + DN_v(2,1)*crLHS411 + DN_v(2,2)*crLHS586 - crLHS590*crLHS794 - crLHS590*crLHS796 + crLHS805 + crLHS872);
    rLHS(8,12)+=gauss_weight*(DN_v(2,0)*crLHS143 + DN_v(2,1)*crLHS415 + DN_v(2,2)*crLHS591 - crLHS595*crLHS794 - crLHS595*crLHS796 + crLHS873);
    rLHS(8,13)+=gauss_weight*(DN_v(2,0)*crLHS160 + DN_v(2,1)*crLHS422 + DN_v(2,2)*crLHS596 - crLHS599*crLHS794 - crLHS599*crLHS796 + crLHS874);
    rLHS(8,14)+=gauss_weight*(DN_v(2,0)*crLHS168 + DN_v(2,1)*crLHS430 + DN_v(2,2)*crLHS600 - crLHS604*crLHS794 - crLHS604*crLHS796 + crLHS810 + crLHS875);
    rLHS(8,15)+=gauss_weight*(DN_v(2,0)*crLHS175 + DN_v(2,1)*crLHS434 + DN_v(2,2)*crLHS605 - crLHS609*crLHS794 - crLHS609*crLHS796 + crLHS876);
    rLHS(8,16)+=gauss_weight*(DN_v(2,0)*crLHS192 + DN_v(2,1)*crLHS441 + DN_v(2,2)*crLHS610 - crLHS613*crLHS794 - crLHS613*crLHS796 + crLHS877);
    rLHS(8,17)+=gauss_weight*(DN_v(2,0)*crLHS200 + DN_v(2,1)*crLHS449 + DN_v(2,2)*crLHS614 - crLHS618*crLHS794 - crLHS618*crLHS796 + crLHS815 + crLHS878);
    rLHS(8,18)+=gauss_weight*(DN_v(2,0)*crLHS207 + DN_v(2,1)*crLHS453 + DN_v(2,2)*crLHS619 - crLHS623*crLHS794 - crLHS623*crLHS796 + crLHS879);
    rLHS(8,19)+=gauss_weight*(DN_v(2,0)*crLHS224 + DN_v(2,1)*crLHS460 + DN_v(2,2)*crLHS624 - crLHS627*crLHS794 - crLHS627*crLHS796 + crLHS880);
    rLHS(8,20)+=gauss_weight*(DN_v(2,0)*crLHS232 + DN_v(2,1)*crLHS468 + DN_v(2,2)*crLHS628 - crLHS632*crLHS794 - crLHS632*crLHS796 + crLHS820 + crLHS881);
    rLHS(8,21)+=gauss_weight*(DN_v(2,0)*crLHS239 + DN_v(2,1)*crLHS472 + DN_v(2,2)*crLHS633 - crLHS637*crLHS794 - crLHS637*crLHS796 + crLHS882);
    rLHS(8,22)+=gauss_weight*(DN_v(2,0)*crLHS256 + DN_v(2,1)*crLHS479 + DN_v(2,2)*crLHS638 - crLHS641*crLHS794 - crLHS641*crLHS796 + crLHS883);
    rLHS(8,23)+=gauss_weight*(DN_v(2,0)*crLHS264 + DN_v(2,1)*crLHS487 + DN_v(2,2)*crLHS642 - crLHS646*crLHS794 - crLHS646*crLHS796 + crLHS825 + crLHS884);
    rLHS(8,24)+=gauss_weight*(DN_v(2,0)*crLHS271 + DN_v(2,1)*crLHS491 + DN_v(2,2)*crLHS647 - crLHS651*crLHS794 - crLHS651*crLHS796 + crLHS885);
    rLHS(8,25)+=gauss_weight*(DN_v(2,0)*crLHS288 + DN_v(2,1)*crLHS498 + DN_v(2,2)*crLHS652 - crLHS655*crLHS794 - crLHS655*crLHS796 + crLHS886);
    rLHS(8,26)+=gauss_weight*(DN_v(2,0)*crLHS296 + DN_v(2,1)*crLHS506 + DN_v(2,2)*crLHS656 - crLHS660*crLHS794 - crLHS660*crLHS796 + crLHS830 + crLHS887);
    rLHS(8,27)+=gauss_weight*(DN_v(2,0)*crLHS303 + DN_v(2,1)*crLHS510 + DN_v(2,2)*crLHS661 - crLHS665*crLHS794 - crLHS665*crLHS796 + crLHS888);
    rLHS(8,28)+=gauss_weight*(DN_v(2,0)*crLHS319 + DN_v(2,1)*crLHS517 + DN_v(2,2)*crLHS666 - crLHS669*crLHS794 - crLHS669*crLHS796 + crLHS889);
    rLHS(8,29)+=gauss_weight*(DN_v(2,0)*crLHS327 + DN_v(2,1)*crLHS525 + DN_v(2,2)*crLHS670 - crLHS674*crLHS794 - crLHS674*crLHS796 + crLHS835 + crLHS890);
    rLHS(8,30)+=gauss_weight*(crLHS676*crLHS795 + crLHS676*crLHS85 + crLHS891);
    rLHS(8,31)+=gauss_weight*(crLHS678*crLHS795 + crLHS678*crLHS85 + crLHS892);
    rLHS(8,32)+=gauss_weight*(crLHS680*crLHS795 + crLHS680*crLHS85 + crLHS893);
    rLHS(8,33)+=gauss_weight*(crLHS682*crLHS795 + crLHS682*crLHS85 + crLHS894);
    rLHS(9,0)+=gauss_weight*(DN_v(3,0)*crLHS0 + DN_v(3,1)*crLHS2 + DN_v(3,2)*crLHS4 + crLHS112 - crLHS19*crLHS895 - crLHS19*crLHS897 + crLHS898);
    rLHS(9,1)+=gauss_weight*(DN_v(3,0)*crLHS26 + DN_v(3,1)*crLHS28 + DN_v(3,2)*crLHS31 - crLHS35*crLHS895 - crLHS35*crLHS897 + crLHS397);
    rLHS(9,2)+=gauss_weight*(DN_v(3,0)*crLHS36 + DN_v(3,1)*crLHS38 + DN_v(3,2)*crLHS40 - crLHS42*crLHS895 - crLHS42*crLHS897 + crLHS578);
    rLHS(9,3)+=gauss_weight*(DN_v(3,0)*crLHS43 + DN_v(3,1)*crLHS45 + DN_v(3,2)*crLHS47 - crLHS56*crLHS895 - crLHS56*crLHS897 + crLHS696 + crLHS899);
    rLHS(9,4)+=gauss_weight*(DN_v(3,0)*crLHS59 + DN_v(3,1)*crLHS61 + DN_v(3,2)*crLHS64 - crLHS67*crLHS895 - crLHS67*crLHS897 + crLHS740);
    rLHS(9,5)+=gauss_weight*(DN_v(3,0)*crLHS68 + DN_v(3,1)*crLHS70 + DN_v(3,2)*crLHS72 - crLHS74*crLHS895 - crLHS74*crLHS897 + crLHS769);
    rLHS(9,6)+=gauss_weight*(DN_v(3,0)*crLHS75 + DN_v(3,1)*crLHS77 + DN_v(3,2)*crLHS79 + crLHS803 - crLHS88*crLHS895 - crLHS88*crLHS897 + crLHS900);
    rLHS(9,7)+=gauss_weight*(DN_v(3,0)*crLHS91 + DN_v(3,1)*crLHS93 + DN_v(3,2)*crLHS96 + crLHS844 - crLHS895*crLHS99 - crLHS897*crLHS99);
    rLHS(9,8)+=gauss_weight*(DN_v(3,0)*crLHS100 + DN_v(3,1)*crLHS102 + DN_v(3,2)*crLHS104 - crLHS106*crLHS895 - crLHS106*crLHS897 + crLHS870);
    rLHS(9,9)+=gauss_weight*(std::pow(DN_v(3,0), 2)*crLHS9 + DN_v(3,0)*crLHS107 + DN_v(3,1)*crLHS109 + DN_v(3,2)*crLHS111 - crLHS120*crLHS895 - crLHS120*crLHS897 + crLHS901);
    rLHS(9,10)+=gauss_weight*(DN_v(3,0)*crLHS123 + DN_v(3,1)*crLHS125 + DN_v(3,2)*crLHS128 - crLHS131*crLHS895 - crLHS131*crLHS897 + crLHS903);
    rLHS(9,11)+=gauss_weight*(DN_v(3,0)*crLHS132 + DN_v(3,1)*crLHS134 + DN_v(3,2)*crLHS136 - crLHS138*crLHS895 - crLHS138*crLHS897 + crLHS904);
    rLHS(9,12)+=gauss_weight*(DN_v(3,0)*crLHS139 + DN_v(3,1)*crLHS141 + DN_v(3,2)*crLHS143 - crLHS152*crLHS895 - crLHS152*crLHS897 + crLHS905 + crLHS907);
    rLHS(9,13)+=gauss_weight*(DN_v(3,0)*crLHS155 + DN_v(3,1)*crLHS157 + DN_v(3,2)*crLHS160 - crLHS163*crLHS895 - crLHS163*crLHS897 + crLHS908);
    rLHS(9,14)+=gauss_weight*(DN_v(3,0)*crLHS164 + DN_v(3,1)*crLHS166 + DN_v(3,2)*crLHS168 - crLHS170*crLHS895 - crLHS170*crLHS897 + crLHS909);
    rLHS(9,15)+=gauss_weight*(DN_v(3,0)*crLHS171 + DN_v(3,1)*crLHS173 + DN_v(3,2)*crLHS175 - crLHS184*crLHS895 - crLHS184*crLHS897 + crLHS910 + crLHS912);
    rLHS(9,16)+=gauss_weight*(DN_v(3,0)*crLHS187 + DN_v(3,1)*crLHS189 + DN_v(3,2)*crLHS192 - crLHS195*crLHS895 - crLHS195*crLHS897 + crLHS913);
    rLHS(9,17)+=gauss_weight*(DN_v(3,0)*crLHS196 + DN_v(3,1)*crLHS198 + DN_v(3,2)*crLHS200 - crLHS202*crLHS895 - crLHS202*crLHS897 + crLHS914);
    rLHS(9,18)+=gauss_weight*(DN_v(3,0)*crLHS203 + DN_v(3,1)*crLHS205 + DN_v(3,2)*crLHS207 - crLHS216*crLHS895 - crLHS216*crLHS897 + crLHS915 + crLHS917);
    rLHS(9,19)+=gauss_weight*(DN_v(3,0)*crLHS219 + DN_v(3,1)*crLHS221 + DN_v(3,2)*crLHS224 - crLHS227*crLHS895 - crLHS227*crLHS897 + crLHS918);
    rLHS(9,20)+=gauss_weight*(DN_v(3,0)*crLHS228 + DN_v(3,1)*crLHS230 + DN_v(3,2)*crLHS232 - crLHS234*crLHS895 - crLHS234*crLHS897 + crLHS919);
    rLHS(9,21)+=gauss_weight*(DN_v(3,0)*crLHS235 + DN_v(3,1)*crLHS237 + DN_v(3,2)*crLHS239 - crLHS248*crLHS895 - crLHS248*crLHS897 + crLHS920 + crLHS922);
    rLHS(9,22)+=gauss_weight*(DN_v(3,0)*crLHS251 + DN_v(3,1)*crLHS253 + DN_v(3,2)*crLHS256 - crLHS259*crLHS895 - crLHS259*crLHS897 + crLHS923);
    rLHS(9,23)+=gauss_weight*(DN_v(3,0)*crLHS260 + DN_v(3,1)*crLHS262 + DN_v(3,2)*crLHS264 - crLHS266*crLHS895 - crLHS266*crLHS897 + crLHS924);
    rLHS(9,24)+=gauss_weight*(DN_v(3,0)*crLHS267 + DN_v(3,1)*crLHS269 + DN_v(3,2)*crLHS271 - crLHS280*crLHS895 - crLHS280*crLHS897 + crLHS925 + crLHS927);
    rLHS(9,25)+=gauss_weight*(DN_v(3,0)*crLHS283 + DN_v(3,1)*crLHS285 + DN_v(3,2)*crLHS288 - crLHS291*crLHS895 - crLHS291*crLHS897 + crLHS928);
    rLHS(9,26)+=gauss_weight*(DN_v(3,0)*crLHS292 + DN_v(3,1)*crLHS294 + DN_v(3,2)*crLHS296 - crLHS298*crLHS895 - crLHS298*crLHS897 + crLHS929);
    rLHS(9,27)+=gauss_weight*(DN_v(3,0)*crLHS299 + DN_v(3,1)*crLHS301 + DN_v(3,2)*crLHS303 - crLHS311*crLHS895 - crLHS311*crLHS897 + crLHS930 + crLHS932);
    rLHS(9,28)+=gauss_weight*(DN_v(3,0)*crLHS314 + DN_v(3,1)*crLHS316 + DN_v(3,2)*crLHS319 - crLHS322*crLHS895 - crLHS322*crLHS897 + crLHS933);
    rLHS(9,29)+=gauss_weight*(DN_v(3,0)*crLHS323 + DN_v(3,1)*crLHS325 + DN_v(3,2)*crLHS327 - crLHS329*crLHS895 - crLHS329*crLHS897 + crLHS934);
    rLHS(9,30)+=gauss_weight*(crLHS117*crLHS331 + crLHS331*crLHS896 + crLHS935);
    rLHS(9,31)+=gauss_weight*(crLHS117*crLHS333 + crLHS333*crLHS896 + crLHS936);
    rLHS(9,32)+=gauss_weight*(crLHS117*crLHS335 + crLHS335*crLHS896 + crLHS937);
    rLHS(9,33)+=gauss_weight*(crLHS117*crLHS337 + crLHS337*crLHS896 + crLHS938);
    rLHS(10,0)+=gauss_weight*(DN_v(3,0)*crLHS2 + DN_v(3,1)*crLHS338 + DN_v(3,2)*crLHS339 + crLHS129 - crLHS342*crLHS895 - crLHS342*crLHS897);
    rLHS(10,1)+=gauss_weight*(DN_v(3,0)*crLHS28 + DN_v(3,1)*crLHS343 + DN_v(3,2)*crLHS345 - crLHS350*crLHS895 - crLHS350*crLHS897 + crLHS404 + crLHS898);
    rLHS(10,2)+=gauss_weight*(DN_v(3,0)*crLHS38 + DN_v(3,1)*crLHS351 + DN_v(3,2)*crLHS353 - crLHS356*crLHS895 - crLHS356*crLHS897 + crLHS583);
    rLHS(10,3)+=gauss_weight*(DN_v(3,0)*crLHS45 + DN_v(3,1)*crLHS357 + DN_v(3,2)*crLHS358 - crLHS362*crLHS895 - crLHS362*crLHS897 + crLHS699);
    rLHS(10,4)+=gauss_weight*(DN_v(3,0)*crLHS61 + DN_v(3,1)*crLHS363 + DN_v(3,2)*crLHS365 - crLHS370*crLHS895 - crLHS370*crLHS897 + crLHS741 + crLHS899);
    rLHS(10,5)+=gauss_weight*(DN_v(3,0)*crLHS70 + DN_v(3,1)*crLHS371 + DN_v(3,2)*crLHS373 - crLHS375*crLHS895 - crLHS375*crLHS897 + crLHS770);
    rLHS(10,6)+=gauss_weight*(DN_v(3,0)*crLHS77 + DN_v(3,1)*crLHS376 + DN_v(3,2)*crLHS377 - crLHS381*crLHS895 - crLHS381*crLHS897 + crLHS806);
    rLHS(10,7)+=gauss_weight*(DN_v(3,0)*crLHS93 + DN_v(3,1)*crLHS382 + DN_v(3,2)*crLHS384 - crLHS389*crLHS895 - crLHS389*crLHS897 + crLHS845 + crLHS900);
    rLHS(10,8)+=gauss_weight*(DN_v(3,0)*crLHS102 + DN_v(3,1)*crLHS390 + DN_v(3,2)*crLHS392 - crLHS394*crLHS895 - crLHS394*crLHS897 + crLHS871);
    rLHS(10,9)+=gauss_weight*(DN_v(3,0)*crLHS109 + DN_v(3,1)*crLHS395 + DN_v(3,2)*crLHS396 - crLHS400*crLHS895 - crLHS400*crLHS897 + crLHS903);
    rLHS(10,10)+=gauss_weight*(DN_v(3,0)*crLHS125 + std::pow(DN_v(3,1), 2)*crLHS9 + DN_v(3,1)*crLHS401 + DN_v(3,2)*crLHS403 - crLHS408*crLHS895 - crLHS408*crLHS897 + crLHS901);
    rLHS(10,11)+=gauss_weight*(DN_v(3,0)*crLHS134 + DN_v(3,1)*crLHS409 + DN_v(3,2)*crLHS411 - crLHS413*crLHS895 - crLHS413*crLHS897 + crLHS940);
    rLHS(10,12)+=gauss_weight*(DN_v(3,0)*crLHS141 + DN_v(3,1)*crLHS414 + DN_v(3,2)*crLHS415 - crLHS419*crLHS895 - crLHS419*crLHS897 + crLHS941);
    rLHS(10,13)+=gauss_weight*(DN_v(3,0)*crLHS157 + DN_v(3,1)*crLHS420 + DN_v(3,2)*crLHS422 - crLHS427*crLHS895 - crLHS427*crLHS897 + crLHS907 + crLHS942);
    rLHS(10,14)+=gauss_weight*(DN_v(3,0)*crLHS166 + DN_v(3,1)*crLHS428 + DN_v(3,2)*crLHS430 - crLHS432*crLHS895 - crLHS432*crLHS897 + crLHS943);
    rLHS(10,15)+=gauss_weight*(DN_v(3,0)*crLHS173 + DN_v(3,1)*crLHS433 + DN_v(3,2)*crLHS434 - crLHS438*crLHS895 - crLHS438*crLHS897 + crLHS944);
    rLHS(10,16)+=gauss_weight*(DN_v(3,0)*crLHS189 + DN_v(3,1)*crLHS439 + DN_v(3,2)*crLHS441 - crLHS446*crLHS895 - crLHS446*crLHS897 + crLHS912 + crLHS945);
    rLHS(10,17)+=gauss_weight*(DN_v(3,0)*crLHS198 + DN_v(3,1)*crLHS447 + DN_v(3,2)*crLHS449 - crLHS451*crLHS895 - crLHS451*crLHS897 + crLHS946);
    rLHS(10,18)+=gauss_weight*(DN_v(3,0)*crLHS205 + DN_v(3,1)*crLHS452 + DN_v(3,2)*crLHS453 - crLHS457*crLHS895 - crLHS457*crLHS897 + crLHS947);
    rLHS(10,19)+=gauss_weight*(DN_v(3,0)*crLHS221 + DN_v(3,1)*crLHS458 + DN_v(3,2)*crLHS460 - crLHS465*crLHS895 - crLHS465*crLHS897 + crLHS917 + crLHS948);
    rLHS(10,20)+=gauss_weight*(DN_v(3,0)*crLHS230 + DN_v(3,1)*crLHS466 + DN_v(3,2)*crLHS468 - crLHS470*crLHS895 - crLHS470*crLHS897 + crLHS949);
    rLHS(10,21)+=gauss_weight*(DN_v(3,0)*crLHS237 + DN_v(3,1)*crLHS471 + DN_v(3,2)*crLHS472 - crLHS476*crLHS895 - crLHS476*crLHS897 + crLHS950);
    rLHS(10,22)+=gauss_weight*(DN_v(3,0)*crLHS253 + DN_v(3,1)*crLHS477 + DN_v(3,2)*crLHS479 - crLHS484*crLHS895 - crLHS484*crLHS897 + crLHS922 + crLHS951);
    rLHS(10,23)+=gauss_weight*(DN_v(3,0)*crLHS262 + DN_v(3,1)*crLHS485 + DN_v(3,2)*crLHS487 - crLHS489*crLHS895 - crLHS489*crLHS897 + crLHS952);
    rLHS(10,24)+=gauss_weight*(DN_v(3,0)*crLHS269 + DN_v(3,1)*crLHS490 + DN_v(3,2)*crLHS491 - crLHS495*crLHS895 - crLHS495*crLHS897 + crLHS953);
    rLHS(10,25)+=gauss_weight*(DN_v(3,0)*crLHS285 + DN_v(3,1)*crLHS496 + DN_v(3,2)*crLHS498 - crLHS503*crLHS895 - crLHS503*crLHS897 + crLHS927 + crLHS954);
    rLHS(10,26)+=gauss_weight*(DN_v(3,0)*crLHS294 + DN_v(3,1)*crLHS504 + DN_v(3,2)*crLHS506 - crLHS508*crLHS895 - crLHS508*crLHS897 + crLHS955);
    rLHS(10,27)+=gauss_weight*(DN_v(3,0)*crLHS301 + DN_v(3,1)*crLHS509 + DN_v(3,2)*crLHS510 - crLHS514*crLHS895 - crLHS514*crLHS897 + crLHS956);
    rLHS(10,28)+=gauss_weight*(DN_v(3,0)*crLHS316 + DN_v(3,1)*crLHS515 + DN_v(3,2)*crLHS517 - crLHS522*crLHS895 - crLHS522*crLHS897 + crLHS932 + crLHS957);
    rLHS(10,29)+=gauss_weight*(DN_v(3,0)*crLHS325 + DN_v(3,1)*crLHS523 + DN_v(3,2)*crLHS525 - crLHS527*crLHS895 - crLHS527*crLHS897 + crLHS958);
    rLHS(10,30)+=gauss_weight*(crLHS117*crLHS529 + crLHS529*crLHS896 + crLHS959);
    rLHS(10,31)+=gauss_weight*(crLHS117*crLHS531 + crLHS531*crLHS896 + crLHS960);
    rLHS(10,32)+=gauss_weight*(crLHS117*crLHS533 + crLHS533*crLHS896 + crLHS961);
    rLHS(10,33)+=gauss_weight*(crLHS117*crLHS535 + crLHS535*crLHS896 + crLHS962);
    rLHS(11,0)+=gauss_weight*(DN_v(3,0)*crLHS4 + DN_v(3,1)*crLHS339 + DN_v(3,2)*crLHS536 + crLHS137 - crLHS539*crLHS895 - crLHS539*crLHS897);
    rLHS(11,1)+=gauss_weight*(DN_v(3,0)*crLHS31 + DN_v(3,1)*crLHS345 + DN_v(3,2)*crLHS540 + crLHS412 - crLHS542*crLHS895 - crLHS542*crLHS897);
    rLHS(11,2)+=gauss_weight*(DN_v(3,0)*crLHS40 + DN_v(3,1)*crLHS353 + DN_v(3,2)*crLHS543 - crLHS547*crLHS895 - crLHS547*crLHS897 + crLHS587 + crLHS898);
    rLHS(11,3)+=gauss_weight*(DN_v(3,0)*crLHS47 + DN_v(3,1)*crLHS358 + DN_v(3,2)*crLHS548 - crLHS553*crLHS895 - crLHS553*crLHS897 + crLHS700);
    rLHS(11,4)+=gauss_weight*(DN_v(3,0)*crLHS64 + DN_v(3,1)*crLHS365 + DN_v(3,2)*crLHS554 - crLHS557*crLHS895 - crLHS557*crLHS897 + crLHS742);
    rLHS(11,5)+=gauss_weight*(DN_v(3,0)*crLHS72 + DN_v(3,1)*crLHS373 + DN_v(3,2)*crLHS558 - crLHS562*crLHS895 - crLHS562*crLHS897 + crLHS771 + crLHS899);
    rLHS(11,6)+=gauss_weight*(DN_v(3,0)*crLHS79 + DN_v(3,1)*crLHS377 + DN_v(3,2)*crLHS563 - crLHS567*crLHS895 - crLHS567*crLHS897 + crLHS807);
    rLHS(11,7)+=gauss_weight*(DN_v(3,0)*crLHS96 + DN_v(3,1)*crLHS384 + DN_v(3,2)*crLHS568 - crLHS571*crLHS895 - crLHS571*crLHS897 + crLHS846);
    rLHS(11,8)+=gauss_weight*(DN_v(3,0)*crLHS104 + DN_v(3,1)*crLHS392 + DN_v(3,2)*crLHS572 - crLHS576*crLHS895 - crLHS576*crLHS897 + crLHS872 + crLHS900);
    rLHS(11,9)+=gauss_weight*(DN_v(3,0)*crLHS111 + DN_v(3,1)*crLHS396 + DN_v(3,2)*crLHS577 - crLHS581*crLHS895 - crLHS581*crLHS897 + crLHS904);
    rLHS(11,10)+=gauss_weight*(DN_v(3,0)*crLHS128 + DN_v(3,1)*crLHS403 + DN_v(3,2)*crLHS582 - crLHS585*crLHS895 - crLHS585*crLHS897 + crLHS940);
    rLHS(11,11)+=gauss_weight*(DN_v(3,0)*crLHS136 + DN_v(3,1)*crLHS411 + std::pow(DN_v(3,2), 2)*crLHS9 + DN_v(3,2)*crLHS586 - crLHS590*crLHS895 - crLHS590*crLHS897 + crLHS901);
    rLHS(11,12)+=gauss_weight*(DN_v(3,0)*crLHS143 + DN_v(3,1)*crLHS415 + DN_v(3,2)*crLHS591 - crLHS595*crLHS895 - crLHS595*crLHS897 + crLHS964);
    rLHS(11,13)+=gauss_weight*(DN_v(3,0)*crLHS160 + DN_v(3,1)*crLHS422 + DN_v(3,2)*crLHS596 - crLHS599*crLHS895 - crLHS599*crLHS897 + crLHS965);
    rLHS(11,14)+=gauss_weight*(DN_v(3,0)*crLHS168 + DN_v(3,1)*crLHS430 + DN_v(3,2)*crLHS600 - crLHS604*crLHS895 - crLHS604*crLHS897 + crLHS907 + crLHS966);
    rLHS(11,15)+=gauss_weight*(DN_v(3,0)*crLHS175 + DN_v(3,1)*crLHS434 + DN_v(3,2)*crLHS605 - crLHS609*crLHS895 - crLHS609*crLHS897 + crLHS967);
    rLHS(11,16)+=gauss_weight*(DN_v(3,0)*crLHS192 + DN_v(3,1)*crLHS441 + DN_v(3,2)*crLHS610 - crLHS613*crLHS895 - crLHS613*crLHS897 + crLHS968);
    rLHS(11,17)+=gauss_weight*(DN_v(3,0)*crLHS200 + DN_v(3,1)*crLHS449 + DN_v(3,2)*crLHS614 - crLHS618*crLHS895 - crLHS618*crLHS897 + crLHS912 + crLHS969);
    rLHS(11,18)+=gauss_weight*(DN_v(3,0)*crLHS207 + DN_v(3,1)*crLHS453 + DN_v(3,2)*crLHS619 - crLHS623*crLHS895 - crLHS623*crLHS897 + crLHS970);
    rLHS(11,19)+=gauss_weight*(DN_v(3,0)*crLHS224 + DN_v(3,1)*crLHS460 + DN_v(3,2)*crLHS624 - crLHS627*crLHS895 - crLHS627*crLHS897 + crLHS971);
    rLHS(11,20)+=gauss_weight*(DN_v(3,0)*crLHS232 + DN_v(3,1)*crLHS468 + DN_v(3,2)*crLHS628 - crLHS632*crLHS895 - crLHS632*crLHS897 + crLHS917 + crLHS972);
    rLHS(11,21)+=gauss_weight*(DN_v(3,0)*crLHS239 + DN_v(3,1)*crLHS472 + DN_v(3,2)*crLHS633 - crLHS637*crLHS895 - crLHS637*crLHS897 + crLHS973);
    rLHS(11,22)+=gauss_weight*(DN_v(3,0)*crLHS256 + DN_v(3,1)*crLHS479 + DN_v(3,2)*crLHS638 - crLHS641*crLHS895 - crLHS641*crLHS897 + crLHS974);
    rLHS(11,23)+=gauss_weight*(DN_v(3,0)*crLHS264 + DN_v(3,1)*crLHS487 + DN_v(3,2)*crLHS642 - crLHS646*crLHS895 - crLHS646*crLHS897 + crLHS922 + crLHS975);
    rLHS(11,24)+=gauss_weight*(DN_v(3,0)*crLHS271 + DN_v(3,1)*crLHS491 + DN_v(3,2)*crLHS647 - crLHS651*crLHS895 - crLHS651*crLHS897 + crLHS976);
    rLHS(11,25)+=gauss_weight*(DN_v(3,0)*crLHS288 + DN_v(3,1)*crLHS498 + DN_v(3,2)*crLHS652 - crLHS655*crLHS895 - crLHS655*crLHS897 + crLHS977);
    rLHS(11,26)+=gauss_weight*(DN_v(3,0)*crLHS296 + DN_v(3,1)*crLHS506 + DN_v(3,2)*crLHS656 - crLHS660*crLHS895 - crLHS660*crLHS897 + crLHS927 + crLHS978);
    rLHS(11,27)+=gauss_weight*(DN_v(3,0)*crLHS303 + DN_v(3,1)*crLHS510 + DN_v(3,2)*crLHS661 - crLHS665*crLHS895 - crLHS665*crLHS897 + crLHS979);
    rLHS(11,28)+=gauss_weight*(DN_v(3,0)*crLHS319 + DN_v(3,1)*crLHS517 + DN_v(3,2)*crLHS666 - crLHS669*crLHS895 - crLHS669*crLHS897 + crLHS980);
    rLHS(11,29)+=gauss_weight*(DN_v(3,0)*crLHS327 + DN_v(3,1)*crLHS525 + DN_v(3,2)*crLHS670 - crLHS674*crLHS895 - crLHS674*crLHS897 + crLHS932 + crLHS981);
    rLHS(11,30)+=gauss_weight*(crLHS117*crLHS676 + crLHS676*crLHS896 + crLHS982);
    rLHS(11,31)+=gauss_weight*(crLHS117*crLHS678 + crLHS678*crLHS896 + crLHS983);
    rLHS(11,32)+=gauss_weight*(crLHS117*crLHS680 + crLHS680*crLHS896 + crLHS984);
    rLHS(11,33)+=gauss_weight*(crLHS117*crLHS682 + crLHS682*crLHS896 + crLHS985);
    rLHS(12,0)+=gauss_weight*(DN_v(4,0)*crLHS0 + DN_v(4,1)*crLHS2 + DN_v(4,2)*crLHS4 + crLHS144 - crLHS19*crLHS986 - crLHS19*crLHS988 + crLHS989);
    rLHS(12,1)+=gauss_weight*(DN_v(4,0)*crLHS26 + DN_v(4,1)*crLHS28 + DN_v(4,2)*crLHS31 - crLHS35*crLHS986 - crLHS35*crLHS988 + crLHS416);
    rLHS(12,2)+=gauss_weight*(DN_v(4,0)*crLHS36 + DN_v(4,1)*crLHS38 + DN_v(4,2)*crLHS40 - crLHS42*crLHS986 - crLHS42*crLHS988 + crLHS592);
    rLHS(12,3)+=gauss_weight*(DN_v(4,0)*crLHS43 + DN_v(4,1)*crLHS45 + DN_v(4,2)*crLHS47 - crLHS56*crLHS986 - crLHS56*crLHS988 + crLHS701 + crLHS990);
    rLHS(12,4)+=gauss_weight*(DN_v(4,0)*crLHS59 + DN_v(4,1)*crLHS61 + DN_v(4,2)*crLHS64 - crLHS67*crLHS986 - crLHS67*crLHS988 + crLHS743);
    rLHS(12,5)+=gauss_weight*(DN_v(4,0)*crLHS68 + DN_v(4,1)*crLHS70 + DN_v(4,2)*crLHS72 - crLHS74*crLHS986 - crLHS74*crLHS988 + crLHS772);
    rLHS(12,6)+=gauss_weight*(DN_v(4,0)*crLHS75 + DN_v(4,1)*crLHS77 + DN_v(4,2)*crLHS79 + crLHS808 - crLHS88*crLHS986 - crLHS88*crLHS988 + crLHS991);
    rLHS(12,7)+=gauss_weight*(DN_v(4,0)*crLHS91 + DN_v(4,1)*crLHS93 + DN_v(4,2)*crLHS96 + crLHS847 - crLHS986*crLHS99 - crLHS988*crLHS99);
    rLHS(12,8)+=gauss_weight*(DN_v(4,0)*crLHS100 + DN_v(4,1)*crLHS102 + DN_v(4,2)*crLHS104 - crLHS106*crLHS986 - crLHS106*crLHS988 + crLHS873);
    rLHS(12,9)+=gauss_weight*(DN_v(4,0)*crLHS107 + DN_v(4,1)*crLHS109 + DN_v(4,2)*crLHS111 - crLHS120*crLHS986 - crLHS120*crLHS988 + crLHS905 + crLHS992);
    rLHS(12,10)+=gauss_weight*(DN_v(4,0)*crLHS123 + DN_v(4,1)*crLHS125 + DN_v(4,2)*crLHS128 - crLHS131*crLHS986 - crLHS131*crLHS988 + crLHS941);
    rLHS(12,11)+=gauss_weight*(DN_v(4,0)*crLHS132 + DN_v(4,1)*crLHS134 + DN_v(4,2)*crLHS136 - crLHS138*crLHS986 - crLHS138*crLHS988 + crLHS964);
    rLHS(12,12)+=gauss_weight*(std::pow(DN_v(4,0), 2)*crLHS9 + DN_v(4,0)*crLHS139 + DN_v(4,1)*crLHS141 + DN_v(4,2)*crLHS143 - crLHS152*crLHS986 - crLHS152*crLHS988 + crLHS993);
    rLHS(12,13)+=gauss_weight*(DN_v(4,0)*crLHS155 + DN_v(4,1)*crLHS157 + DN_v(4,2)*crLHS160 - crLHS163*crLHS986 - crLHS163*crLHS988 + crLHS995);
    rLHS(12,14)+=gauss_weight*(DN_v(4,0)*crLHS164 + DN_v(4,1)*crLHS166 + DN_v(4,2)*crLHS168 - crLHS170*crLHS986 - crLHS170*crLHS988 + crLHS996);
    rLHS(12,15)+=gauss_weight*(DN_v(4,0)*crLHS171 + DN_v(4,1)*crLHS173 + DN_v(4,2)*crLHS175 - crLHS184*crLHS986 - crLHS184*crLHS988 + crLHS997 + crLHS999);
    rLHS(12,16)+=gauss_weight*(DN_v(4,0)*crLHS187 + DN_v(4,1)*crLHS189 + DN_v(4,2)*crLHS192 + crLHS1000 - crLHS195*crLHS986 - crLHS195*crLHS988);
    rLHS(12,17)+=gauss_weight*(DN_v(4,0)*crLHS196 + DN_v(4,1)*crLHS198 + DN_v(4,2)*crLHS200 + crLHS1001 - crLHS202*crLHS986 - crLHS202*crLHS988);
    rLHS(12,18)+=gauss_weight*(DN_v(4,0)*crLHS203 + DN_v(4,1)*crLHS205 + DN_v(4,2)*crLHS207 + crLHS1002 + crLHS1004 - crLHS216*crLHS986 - crLHS216*crLHS988);
    rLHS(12,19)+=gauss_weight*(DN_v(4,0)*crLHS219 + DN_v(4,1)*crLHS221 + DN_v(4,2)*crLHS224 + crLHS1005 - crLHS227*crLHS986 - crLHS227*crLHS988);
    rLHS(12,20)+=gauss_weight*(DN_v(4,0)*crLHS228 + DN_v(4,1)*crLHS230 + DN_v(4,2)*crLHS232 + crLHS1006 - crLHS234*crLHS986 - crLHS234*crLHS988);
    rLHS(12,21)+=gauss_weight*(DN_v(4,0)*crLHS235 + DN_v(4,1)*crLHS237 + DN_v(4,2)*crLHS239 + crLHS1007 + crLHS1009 - crLHS248*crLHS986 - crLHS248*crLHS988);
    rLHS(12,22)+=gauss_weight*(DN_v(4,0)*crLHS251 + DN_v(4,1)*crLHS253 + DN_v(4,2)*crLHS256 + crLHS1010 - crLHS259*crLHS986 - crLHS259*crLHS988);
    rLHS(12,23)+=gauss_weight*(DN_v(4,0)*crLHS260 + DN_v(4,1)*crLHS262 + DN_v(4,2)*crLHS264 + crLHS1011 - crLHS266*crLHS986 - crLHS266*crLHS988);
    rLHS(12,24)+=gauss_weight*(DN_v(4,0)*crLHS267 + DN_v(4,1)*crLHS269 + DN_v(4,2)*crLHS271 + crLHS1012 + crLHS1014 - crLHS280*crLHS986 - crLHS280*crLHS988);
    rLHS(12,25)+=gauss_weight*(DN_v(4,0)*crLHS283 + DN_v(4,1)*crLHS285 + DN_v(4,2)*crLHS288 + crLHS1015 - crLHS291*crLHS986 - crLHS291*crLHS988);
    rLHS(12,26)+=gauss_weight*(DN_v(4,0)*crLHS292 + DN_v(4,1)*crLHS294 + DN_v(4,2)*crLHS296 + crLHS1016 - crLHS298*crLHS986 - crLHS298*crLHS988);
    rLHS(12,27)+=gauss_weight*(DN_v(4,0)*crLHS299 + DN_v(4,1)*crLHS301 + DN_v(4,2)*crLHS303 + crLHS1017 + crLHS1019 - crLHS311*crLHS986 - crLHS311*crLHS988);
    rLHS(12,28)+=gauss_weight*(DN_v(4,0)*crLHS314 + DN_v(4,1)*crLHS316 + DN_v(4,2)*crLHS319 + crLHS1020 - crLHS322*crLHS986 - crLHS322*crLHS988);
    rLHS(12,29)+=gauss_weight*(DN_v(4,0)*crLHS323 + DN_v(4,1)*crLHS325 + DN_v(4,2)*crLHS327 + crLHS1021 - crLHS329*crLHS986 - crLHS329*crLHS988);
    rLHS(12,30)+=gauss_weight*(crLHS1022 + crLHS149*crLHS331 + crLHS331*crLHS987);
    rLHS(12,31)+=gauss_weight*(crLHS1023 + crLHS149*crLHS333 + crLHS333*crLHS987);
    rLHS(12,32)+=gauss_weight*(crLHS1024 + crLHS149*crLHS335 + crLHS335*crLHS987);
    rLHS(12,33)+=gauss_weight*(crLHS1025 + crLHS149*crLHS337 + crLHS337*crLHS987);
    rLHS(13,0)+=gauss_weight*(DN_v(4,0)*crLHS2 + DN_v(4,1)*crLHS338 + DN_v(4,2)*crLHS339 + crLHS161 - crLHS342*crLHS986 - crLHS342*crLHS988);
    rLHS(13,1)+=gauss_weight*(DN_v(4,0)*crLHS28 + DN_v(4,1)*crLHS343 + DN_v(4,2)*crLHS345 - crLHS350*crLHS986 - crLHS350*crLHS988 + crLHS423 + crLHS989);
    rLHS(13,2)+=gauss_weight*(DN_v(4,0)*crLHS38 + DN_v(4,1)*crLHS351 + DN_v(4,2)*crLHS353 - crLHS356*crLHS986 - crLHS356*crLHS988 + crLHS597);
    rLHS(13,3)+=gauss_weight*(DN_v(4,0)*crLHS45 + DN_v(4,1)*crLHS357 + DN_v(4,2)*crLHS358 - crLHS362*crLHS986 - crLHS362*crLHS988 + crLHS704);
    rLHS(13,4)+=gauss_weight*(DN_v(4,0)*crLHS61 + DN_v(4,1)*crLHS363 + DN_v(4,2)*crLHS365 - crLHS370*crLHS986 - crLHS370*crLHS988 + crLHS744 + crLHS990);
    rLHS(13,5)+=gauss_weight*(DN_v(4,0)*crLHS70 + DN_v(4,1)*crLHS371 + DN_v(4,2)*crLHS373 - crLHS375*crLHS986 - crLHS375*crLHS988 + crLHS773);
    rLHS(13,6)+=gauss_weight*(DN_v(4,0)*crLHS77 + DN_v(4,1)*crLHS376 + DN_v(4,2)*crLHS377 - crLHS381*crLHS986 - crLHS381*crLHS988 + crLHS811);
    rLHS(13,7)+=gauss_weight*(DN_v(4,0)*crLHS93 + DN_v(4,1)*crLHS382 + DN_v(4,2)*crLHS384 - crLHS389*crLHS986 - crLHS389*crLHS988 + crLHS848 + crLHS991);
    rLHS(13,8)+=gauss_weight*(DN_v(4,0)*crLHS102 + DN_v(4,1)*crLHS390 + DN_v(4,2)*crLHS392 - crLHS394*crLHS986 - crLHS394*crLHS988 + crLHS874);
    rLHS(13,9)+=gauss_weight*(DN_v(4,0)*crLHS109 + DN_v(4,1)*crLHS395 + DN_v(4,2)*crLHS396 - crLHS400*crLHS986 - crLHS400*crLHS988 + crLHS908);
    rLHS(13,10)+=gauss_weight*(DN_v(4,0)*crLHS125 + DN_v(4,1)*crLHS401 + DN_v(4,2)*crLHS403 - crLHS408*crLHS986 - crLHS408*crLHS988 + crLHS942 + crLHS992);
    rLHS(13,11)+=gauss_weight*(DN_v(4,0)*crLHS134 + DN_v(4,1)*crLHS409 + DN_v(4,2)*crLHS411 - crLHS413*crLHS986 - crLHS413*crLHS988 + crLHS965);
    rLHS(13,12)+=gauss_weight*(DN_v(4,0)*crLHS141 + DN_v(4,1)*crLHS414 + DN_v(4,2)*crLHS415 - crLHS419*crLHS986 - crLHS419*crLHS988 + crLHS995);
    rLHS(13,13)+=gauss_weight*(DN_v(4,0)*crLHS157 + std::pow(DN_v(4,1), 2)*crLHS9 + DN_v(4,1)*crLHS420 + DN_v(4,2)*crLHS422 - crLHS427*crLHS986 - crLHS427*crLHS988 + crLHS993);
    rLHS(13,14)+=gauss_weight*(DN_v(4,0)*crLHS166 + DN_v(4,1)*crLHS428 + DN_v(4,2)*crLHS430 + crLHS1027 - crLHS432*crLHS986 - crLHS432*crLHS988);
    rLHS(13,15)+=gauss_weight*(DN_v(4,0)*crLHS173 + DN_v(4,1)*crLHS433 + DN_v(4,2)*crLHS434 + crLHS1028 - crLHS438*crLHS986 - crLHS438*crLHS988);
    rLHS(13,16)+=gauss_weight*(DN_v(4,0)*crLHS189 + DN_v(4,1)*crLHS439 + DN_v(4,2)*crLHS441 + crLHS1029 - crLHS446*crLHS986 - crLHS446*crLHS988 + crLHS999);
    rLHS(13,17)+=gauss_weight*(DN_v(4,0)*crLHS198 + DN_v(4,1)*crLHS447 + DN_v(4,2)*crLHS449 + crLHS1030 - crLHS451*crLHS986 - crLHS451*crLHS988);
    rLHS(13,18)+=gauss_weight*(DN_v(4,0)*crLHS205 + DN_v(4,1)*crLHS452 + DN_v(4,2)*crLHS453 + crLHS1031 - crLHS457*crLHS986 - crLHS457*crLHS988);
    rLHS(13,19)+=gauss_weight*(DN_v(4,0)*crLHS221 + DN_v(4,1)*crLHS458 + DN_v(4,2)*crLHS460 + crLHS1004 + crLHS1032 - crLHS465*crLHS986 - crLHS465*crLHS988);
    rLHS(13,20)+=gauss_weight*(DN_v(4,0)*crLHS230 + DN_v(4,1)*crLHS466 + DN_v(4,2)*crLHS468 + crLHS1033 - crLHS470*crLHS986 - crLHS470*crLHS988);
    rLHS(13,21)+=gauss_weight*(DN_v(4,0)*crLHS237 + DN_v(4,1)*crLHS471 + DN_v(4,2)*crLHS472 + crLHS1034 - crLHS476*crLHS986 - crLHS476*crLHS988);
    rLHS(13,22)+=gauss_weight*(DN_v(4,0)*crLHS253 + DN_v(4,1)*crLHS477 + DN_v(4,2)*crLHS479 + crLHS1009 + crLHS1035 - crLHS484*crLHS986 - crLHS484*crLHS988);
    rLHS(13,23)+=gauss_weight*(DN_v(4,0)*crLHS262 + DN_v(4,1)*crLHS485 + DN_v(4,2)*crLHS487 + crLHS1036 - crLHS489*crLHS986 - crLHS489*crLHS988);
    rLHS(13,24)+=gauss_weight*(DN_v(4,0)*crLHS269 + DN_v(4,1)*crLHS490 + DN_v(4,2)*crLHS491 + crLHS1037 - crLHS495*crLHS986 - crLHS495*crLHS988);
    rLHS(13,25)+=gauss_weight*(DN_v(4,0)*crLHS285 + DN_v(4,1)*crLHS496 + DN_v(4,2)*crLHS498 + crLHS1014 + crLHS1038 - crLHS503*crLHS986 - crLHS503*crLHS988);
    rLHS(13,26)+=gauss_weight*(DN_v(4,0)*crLHS294 + DN_v(4,1)*crLHS504 + DN_v(4,2)*crLHS506 + crLHS1039 - crLHS508*crLHS986 - crLHS508*crLHS988);
    rLHS(13,27)+=gauss_weight*(DN_v(4,0)*crLHS301 + DN_v(4,1)*crLHS509 + DN_v(4,2)*crLHS510 + crLHS1040 - crLHS514*crLHS986 - crLHS514*crLHS988);
    rLHS(13,28)+=gauss_weight*(DN_v(4,0)*crLHS316 + DN_v(4,1)*crLHS515 + DN_v(4,2)*crLHS517 + crLHS1019 + crLHS1041 - crLHS522*crLHS986 - crLHS522*crLHS988);
    rLHS(13,29)+=gauss_weight*(DN_v(4,0)*crLHS325 + DN_v(4,1)*crLHS523 + DN_v(4,2)*crLHS525 + crLHS1042 - crLHS527*crLHS986 - crLHS527*crLHS988);
    rLHS(13,30)+=gauss_weight*(crLHS1043 + crLHS149*crLHS529 + crLHS529*crLHS987);
    rLHS(13,31)+=gauss_weight*(crLHS1044 + crLHS149*crLHS531 + crLHS531*crLHS987);
    rLHS(13,32)+=gauss_weight*(crLHS1045 + crLHS149*crLHS533 + crLHS533*crLHS987);
    rLHS(13,33)+=gauss_weight*(crLHS1046 + crLHS149*crLHS535 + crLHS535*crLHS987);
    rLHS(14,0)+=gauss_weight*(DN_v(4,0)*crLHS4 + DN_v(4,1)*crLHS339 + DN_v(4,2)*crLHS536 + crLHS169 - crLHS539*crLHS986 - crLHS539*crLHS988);
    rLHS(14,1)+=gauss_weight*(DN_v(4,0)*crLHS31 + DN_v(4,1)*crLHS345 + DN_v(4,2)*crLHS540 + crLHS431 - crLHS542*crLHS986 - crLHS542*crLHS988);
    rLHS(14,2)+=gauss_weight*(DN_v(4,0)*crLHS40 + DN_v(4,1)*crLHS353 + DN_v(4,2)*crLHS543 - crLHS547*crLHS986 - crLHS547*crLHS988 + crLHS601 + crLHS989);
    rLHS(14,3)+=gauss_weight*(DN_v(4,0)*crLHS47 + DN_v(4,1)*crLHS358 + DN_v(4,2)*crLHS548 - crLHS553*crLHS986 - crLHS553*crLHS988 + crLHS705);
    rLHS(14,4)+=gauss_weight*(DN_v(4,0)*crLHS64 + DN_v(4,1)*crLHS365 + DN_v(4,2)*crLHS554 - crLHS557*crLHS986 - crLHS557*crLHS988 + crLHS745);
    rLHS(14,5)+=gauss_weight*(DN_v(4,0)*crLHS72 + DN_v(4,1)*crLHS373 + DN_v(4,2)*crLHS558 - crLHS562*crLHS986 - crLHS562*crLHS988 + crLHS774 + crLHS990);
    rLHS(14,6)+=gauss_weight*(DN_v(4,0)*crLHS79 + DN_v(4,1)*crLHS377 + DN_v(4,2)*crLHS563 - crLHS567*crLHS986 - crLHS567*crLHS988 + crLHS812);
    rLHS(14,7)+=gauss_weight*(DN_v(4,0)*crLHS96 + DN_v(4,1)*crLHS384 + DN_v(4,2)*crLHS568 - crLHS571*crLHS986 - crLHS571*crLHS988 + crLHS849);
    rLHS(14,8)+=gauss_weight*(DN_v(4,0)*crLHS104 + DN_v(4,1)*crLHS392 + DN_v(4,2)*crLHS572 - crLHS576*crLHS986 - crLHS576*crLHS988 + crLHS875 + crLHS991);
    rLHS(14,9)+=gauss_weight*(DN_v(4,0)*crLHS111 + DN_v(4,1)*crLHS396 + DN_v(4,2)*crLHS577 - crLHS581*crLHS986 - crLHS581*crLHS988 + crLHS909);
    rLHS(14,10)+=gauss_weight*(DN_v(4,0)*crLHS128 + DN_v(4,1)*crLHS403 + DN_v(4,2)*crLHS582 - crLHS585*crLHS986 - crLHS585*crLHS988 + crLHS943);
    rLHS(14,11)+=gauss_weight*(DN_v(4,0)*crLHS136 + DN_v(4,1)*crLHS411 + DN_v(4,2)*crLHS586 - crLHS590*crLHS986 - crLHS590*crLHS988 + crLHS966 + crLHS992);
    rLHS(14,12)+=gauss_weight*(DN_v(4,0)*crLHS143 + DN_v(4,1)*crLHS415 + DN_v(4,2)*crLHS591 - crLHS595*crLHS986 - crLHS595*crLHS988 + crLHS996);
    rLHS(14,13)+=gauss_weight*(DN_v(4,0)*crLHS160 + DN_v(4,1)*crLHS422 + DN_v(4,2)*crLHS596 + crLHS1027 - crLHS599*crLHS986 - crLHS599*crLHS988);
    rLHS(14,14)+=gauss_weight*(DN_v(4,0)*crLHS168 + DN_v(4,1)*crLHS430 + std::pow(DN_v(4,2), 2)*crLHS9 + DN_v(4,2)*crLHS600 - crLHS604*crLHS986 - crLHS604*crLHS988 + crLHS993);
    rLHS(14,15)+=gauss_weight*(DN_v(4,0)*crLHS175 + DN_v(4,1)*crLHS434 + DN_v(4,2)*crLHS605 + crLHS1048 - crLHS609*crLHS986 - crLHS609*crLHS988);
    rLHS(14,16)+=gauss_weight*(DN_v(4,0)*crLHS192 + DN_v(4,1)*crLHS441 + DN_v(4,2)*crLHS610 + crLHS1049 - crLHS613*crLHS986 - crLHS613*crLHS988);
    rLHS(14,17)+=gauss_weight*(DN_v(4,0)*crLHS200 + DN_v(4,1)*crLHS449 + DN_v(4,2)*crLHS614 + crLHS1050 - crLHS618*crLHS986 - crLHS618*crLHS988 + crLHS999);
    rLHS(14,18)+=gauss_weight*(DN_v(4,0)*crLHS207 + DN_v(4,1)*crLHS453 + DN_v(4,2)*crLHS619 + crLHS1051 - crLHS623*crLHS986 - crLHS623*crLHS988);
    rLHS(14,19)+=gauss_weight*(DN_v(4,0)*crLHS224 + DN_v(4,1)*crLHS460 + DN_v(4,2)*crLHS624 + crLHS1052 - crLHS627*crLHS986 - crLHS627*crLHS988);
    rLHS(14,20)+=gauss_weight*(DN_v(4,0)*crLHS232 + DN_v(4,1)*crLHS468 + DN_v(4,2)*crLHS628 + crLHS1004 + crLHS1053 - crLHS632*crLHS986 - crLHS632*crLHS988);
    rLHS(14,21)+=gauss_weight*(DN_v(4,0)*crLHS239 + DN_v(4,1)*crLHS472 + DN_v(4,2)*crLHS633 + crLHS1054 - crLHS637*crLHS986 - crLHS637*crLHS988);
    rLHS(14,22)+=gauss_weight*(DN_v(4,0)*crLHS256 + DN_v(4,1)*crLHS479 + DN_v(4,2)*crLHS638 + crLHS1055 - crLHS641*crLHS986 - crLHS641*crLHS988);
    rLHS(14,23)+=gauss_weight*(DN_v(4,0)*crLHS264 + DN_v(4,1)*crLHS487 + DN_v(4,2)*crLHS642 + crLHS1009 + crLHS1056 - crLHS646*crLHS986 - crLHS646*crLHS988);
    rLHS(14,24)+=gauss_weight*(DN_v(4,0)*crLHS271 + DN_v(4,1)*crLHS491 + DN_v(4,2)*crLHS647 + crLHS1057 - crLHS651*crLHS986 - crLHS651*crLHS988);
    rLHS(14,25)+=gauss_weight*(DN_v(4,0)*crLHS288 + DN_v(4,1)*crLHS498 + DN_v(4,2)*crLHS652 + crLHS1058 - crLHS655*crLHS986 - crLHS655*crLHS988);
    rLHS(14,26)+=gauss_weight*(DN_v(4,0)*crLHS296 + DN_v(4,1)*crLHS506 + DN_v(4,2)*crLHS656 + crLHS1014 + crLHS1059 - crLHS660*crLHS986 - crLHS660*crLHS988);
    rLHS(14,27)+=gauss_weight*(DN_v(4,0)*crLHS303 + DN_v(4,1)*crLHS510 + DN_v(4,2)*crLHS661 + crLHS1060 - crLHS665*crLHS986 - crLHS665*crLHS988);
    rLHS(14,28)+=gauss_weight*(DN_v(4,0)*crLHS319 + DN_v(4,1)*crLHS517 + DN_v(4,2)*crLHS666 + crLHS1061 - crLHS669*crLHS986 - crLHS669*crLHS988);
    rLHS(14,29)+=gauss_weight*(DN_v(4,0)*crLHS327 + DN_v(4,1)*crLHS525 + DN_v(4,2)*crLHS670 + crLHS1019 + crLHS1062 - crLHS674*crLHS986 - crLHS674*crLHS988);
    rLHS(14,30)+=gauss_weight*(crLHS1063 + crLHS149*crLHS676 + crLHS676*crLHS987);
    rLHS(14,31)+=gauss_weight*(crLHS1064 + crLHS149*crLHS678 + crLHS678*crLHS987);
    rLHS(14,32)+=gauss_weight*(crLHS1065 + crLHS149*crLHS680 + crLHS680*crLHS987);
    rLHS(14,33)+=gauss_weight*(crLHS1066 + crLHS149*crLHS682 + crLHS682*crLHS987);
    rLHS(15,0)+=gauss_weight*(DN_v(5,0)*crLHS0 + DN_v(5,1)*crLHS2 + DN_v(5,2)*crLHS4 - crLHS1067*crLHS19 - crLHS1069*crLHS19 + crLHS1070 + crLHS176);
    rLHS(15,1)+=gauss_weight*(DN_v(5,0)*crLHS26 + DN_v(5,1)*crLHS28 + DN_v(5,2)*crLHS31 - crLHS1067*crLHS35 - crLHS1069*crLHS35 + crLHS435);
    rLHS(15,2)+=gauss_weight*(DN_v(5,0)*crLHS36 + DN_v(5,1)*crLHS38 + DN_v(5,2)*crLHS40 - crLHS1067*crLHS42 - crLHS1069*crLHS42 + crLHS606);
    rLHS(15,3)+=gauss_weight*(DN_v(5,0)*crLHS43 + DN_v(5,1)*crLHS45 + DN_v(5,2)*crLHS47 - crLHS1067*crLHS56 - crLHS1069*crLHS56 + crLHS1071 + crLHS706);
    rLHS(15,4)+=gauss_weight*(DN_v(5,0)*crLHS59 + DN_v(5,1)*crLHS61 + DN_v(5,2)*crLHS64 - crLHS1067*crLHS67 - crLHS1069*crLHS67 + crLHS746);
    rLHS(15,5)+=gauss_weight*(DN_v(5,0)*crLHS68 + DN_v(5,1)*crLHS70 + DN_v(5,2)*crLHS72 - crLHS1067*crLHS74 - crLHS1069*crLHS74 + crLHS775);
    rLHS(15,6)+=gauss_weight*(DN_v(5,0)*crLHS75 + DN_v(5,1)*crLHS77 + DN_v(5,2)*crLHS79 - crLHS1067*crLHS88 - crLHS1069*crLHS88 + crLHS1072 + crLHS813);
    rLHS(15,7)+=gauss_weight*(DN_v(5,0)*crLHS91 + DN_v(5,1)*crLHS93 + DN_v(5,2)*crLHS96 - crLHS1067*crLHS99 - crLHS1069*crLHS99 + crLHS850);
    rLHS(15,8)+=gauss_weight*(DN_v(5,0)*crLHS100 + DN_v(5,1)*crLHS102 + DN_v(5,2)*crLHS104 - crLHS106*crLHS1067 - crLHS106*crLHS1069 + crLHS876);
    rLHS(15,9)+=gauss_weight*(DN_v(5,0)*crLHS107 + DN_v(5,1)*crLHS109 + DN_v(5,2)*crLHS111 - crLHS1067*crLHS120 - crLHS1069*crLHS120 + crLHS1073 + crLHS910);
    rLHS(15,10)+=gauss_weight*(DN_v(5,0)*crLHS123 + DN_v(5,1)*crLHS125 + DN_v(5,2)*crLHS128 - crLHS1067*crLHS131 - crLHS1069*crLHS131 + crLHS944);
    rLHS(15,11)+=gauss_weight*(DN_v(5,0)*crLHS132 + DN_v(5,1)*crLHS134 + DN_v(5,2)*crLHS136 - crLHS1067*crLHS138 - crLHS1069*crLHS138 + crLHS967);
    rLHS(15,12)+=gauss_weight*(DN_v(5,0)*crLHS139 + DN_v(5,1)*crLHS141 + DN_v(5,2)*crLHS143 - crLHS1067*crLHS152 - crLHS1069*crLHS152 + crLHS1074 + crLHS997);
    rLHS(15,13)+=gauss_weight*(DN_v(5,0)*crLHS155 + DN_v(5,1)*crLHS157 + DN_v(5,2)*crLHS160 + crLHS1028 - crLHS1067*crLHS163 - crLHS1069*crLHS163);
    rLHS(15,14)+=gauss_weight*(DN_v(5,0)*crLHS164 + DN_v(5,1)*crLHS166 + DN_v(5,2)*crLHS168 + crLHS1048 - crLHS1067*crLHS170 - crLHS1069*crLHS170);
    rLHS(15,15)+=gauss_weight*(std::pow(DN_v(5,0), 2)*crLHS9 + DN_v(5,0)*crLHS171 + DN_v(5,1)*crLHS173 + DN_v(5,2)*crLHS175 - crLHS1067*crLHS184 - crLHS1069*crLHS184 + crLHS1075);
    rLHS(15,16)+=gauss_weight*(DN_v(5,0)*crLHS187 + DN_v(5,1)*crLHS189 + DN_v(5,2)*crLHS192 - crLHS1067*crLHS195 - crLHS1069*crLHS195 + crLHS1077);
    rLHS(15,17)+=gauss_weight*(DN_v(5,0)*crLHS196 + DN_v(5,1)*crLHS198 + DN_v(5,2)*crLHS200 - crLHS1067*crLHS202 - crLHS1069*crLHS202 + crLHS1078);
    rLHS(15,18)+=gauss_weight*(DN_v(5,0)*crLHS203 + DN_v(5,1)*crLHS205 + DN_v(5,2)*crLHS207 - crLHS1067*crLHS216 - crLHS1069*crLHS216 + crLHS1079 + crLHS1081);
    rLHS(15,19)+=gauss_weight*(DN_v(5,0)*crLHS219 + DN_v(5,1)*crLHS221 + DN_v(5,2)*crLHS224 - crLHS1067*crLHS227 - crLHS1069*crLHS227 + crLHS1082);
    rLHS(15,20)+=gauss_weight*(DN_v(5,0)*crLHS228 + DN_v(5,1)*crLHS230 + DN_v(5,2)*crLHS232 - crLHS1067*crLHS234 - crLHS1069*crLHS234 + crLHS1083);
    rLHS(15,21)+=gauss_weight*(DN_v(5,0)*crLHS235 + DN_v(5,1)*crLHS237 + DN_v(5,2)*crLHS239 - crLHS1067*crLHS248 - crLHS1069*crLHS248 + crLHS1084 + crLHS1086);
    rLHS(15,22)+=gauss_weight*(DN_v(5,0)*crLHS251 + DN_v(5,1)*crLHS253 + DN_v(5,2)*crLHS256 - crLHS1067*crLHS259 - crLHS1069*crLHS259 + crLHS1087);
    rLHS(15,23)+=gauss_weight*(DN_v(5,0)*crLHS260 + DN_v(5,1)*crLHS262 + DN_v(5,2)*crLHS264 - crLHS1067*crLHS266 - crLHS1069*crLHS266 + crLHS1088);
    rLHS(15,24)+=gauss_weight*(DN_v(5,0)*crLHS267 + DN_v(5,1)*crLHS269 + DN_v(5,2)*crLHS271 - crLHS1067*crLHS280 - crLHS1069*crLHS280 + crLHS1089 + crLHS1091);
    rLHS(15,25)+=gauss_weight*(DN_v(5,0)*crLHS283 + DN_v(5,1)*crLHS285 + DN_v(5,2)*crLHS288 - crLHS1067*crLHS291 - crLHS1069*crLHS291 + crLHS1092);
    rLHS(15,26)+=gauss_weight*(DN_v(5,0)*crLHS292 + DN_v(5,1)*crLHS294 + DN_v(5,2)*crLHS296 - crLHS1067*crLHS298 - crLHS1069*crLHS298 + crLHS1093);
    rLHS(15,27)+=gauss_weight*(DN_v(5,0)*crLHS299 + DN_v(5,1)*crLHS301 + DN_v(5,2)*crLHS303 - crLHS1067*crLHS311 - crLHS1069*crLHS311 + crLHS1094 + crLHS1096);
    rLHS(15,28)+=gauss_weight*(DN_v(5,0)*crLHS314 + DN_v(5,1)*crLHS316 + DN_v(5,2)*crLHS319 - crLHS1067*crLHS322 - crLHS1069*crLHS322 + crLHS1097);
    rLHS(15,29)+=gauss_weight*(DN_v(5,0)*crLHS323 + DN_v(5,1)*crLHS325 + DN_v(5,2)*crLHS327 - crLHS1067*crLHS329 - crLHS1069*crLHS329 + crLHS1098);
    rLHS(15,30)+=gauss_weight*(crLHS1068*crLHS331 + crLHS1099 + crLHS181*crLHS331);
    rLHS(15,31)+=gauss_weight*(crLHS1068*crLHS333 + crLHS1100 + crLHS181*crLHS333);
    rLHS(15,32)+=gauss_weight*(crLHS1068*crLHS335 + crLHS1101 + crLHS181*crLHS335);
    rLHS(15,33)+=gauss_weight*(crLHS1068*crLHS337 + crLHS1102 + crLHS181*crLHS337);
    rLHS(16,0)+=gauss_weight*(DN_v(5,0)*crLHS2 + DN_v(5,1)*crLHS338 + DN_v(5,2)*crLHS339 - crLHS1067*crLHS342 - crLHS1069*crLHS342 + crLHS193);
    rLHS(16,1)+=gauss_weight*(DN_v(5,0)*crLHS28 + DN_v(5,1)*crLHS343 + DN_v(5,2)*crLHS345 - crLHS1067*crLHS350 - crLHS1069*crLHS350 + crLHS1070 + crLHS442);
    rLHS(16,2)+=gauss_weight*(DN_v(5,0)*crLHS38 + DN_v(5,1)*crLHS351 + DN_v(5,2)*crLHS353 - crLHS1067*crLHS356 - crLHS1069*crLHS356 + crLHS611);
    rLHS(16,3)+=gauss_weight*(DN_v(5,0)*crLHS45 + DN_v(5,1)*crLHS357 + DN_v(5,2)*crLHS358 - crLHS1067*crLHS362 - crLHS1069*crLHS362 + crLHS709);
    rLHS(16,4)+=gauss_weight*(DN_v(5,0)*crLHS61 + DN_v(5,1)*crLHS363 + DN_v(5,2)*crLHS365 - crLHS1067*crLHS370 - crLHS1069*crLHS370 + crLHS1071 + crLHS747);
    rLHS(16,5)+=gauss_weight*(DN_v(5,0)*crLHS70 + DN_v(5,1)*crLHS371 + DN_v(5,2)*crLHS373 - crLHS1067*crLHS375 - crLHS1069*crLHS375 + crLHS776);
    rLHS(16,6)+=gauss_weight*(DN_v(5,0)*crLHS77 + DN_v(5,1)*crLHS376 + DN_v(5,2)*crLHS377 - crLHS1067*crLHS381 - crLHS1069*crLHS381 + crLHS816);
    rLHS(16,7)+=gauss_weight*(DN_v(5,0)*crLHS93 + DN_v(5,1)*crLHS382 + DN_v(5,2)*crLHS384 - crLHS1067*crLHS389 - crLHS1069*crLHS389 + crLHS1072 + crLHS851);
    rLHS(16,8)+=gauss_weight*(DN_v(5,0)*crLHS102 + DN_v(5,1)*crLHS390 + DN_v(5,2)*crLHS392 - crLHS1067*crLHS394 - crLHS1069*crLHS394 + crLHS877);
    rLHS(16,9)+=gauss_weight*(DN_v(5,0)*crLHS109 + DN_v(5,1)*crLHS395 + DN_v(5,2)*crLHS396 - crLHS1067*crLHS400 - crLHS1069*crLHS400 + crLHS913);
    rLHS(16,10)+=gauss_weight*(DN_v(5,0)*crLHS125 + DN_v(5,1)*crLHS401 + DN_v(5,2)*crLHS403 - crLHS1067*crLHS408 - crLHS1069*crLHS408 + crLHS1073 + crLHS945);
    rLHS(16,11)+=gauss_weight*(DN_v(5,0)*crLHS134 + DN_v(5,1)*crLHS409 + DN_v(5,2)*crLHS411 - crLHS1067*crLHS413 - crLHS1069*crLHS413 + crLHS968);
    rLHS(16,12)+=gauss_weight*(DN_v(5,0)*crLHS141 + DN_v(5,1)*crLHS414 + DN_v(5,2)*crLHS415 + crLHS1000 - crLHS1067*crLHS419 - crLHS1069*crLHS419);
    rLHS(16,13)+=gauss_weight*(DN_v(5,0)*crLHS157 + DN_v(5,1)*crLHS420 + DN_v(5,2)*crLHS422 + crLHS1029 - crLHS1067*crLHS427 - crLHS1069*crLHS427 + crLHS1074);
    rLHS(16,14)+=gauss_weight*(DN_v(5,0)*crLHS166 + DN_v(5,1)*crLHS428 + DN_v(5,2)*crLHS430 + crLHS1049 - crLHS1067*crLHS432 - crLHS1069*crLHS432);
    rLHS(16,15)+=gauss_weight*(DN_v(5,0)*crLHS173 + DN_v(5,1)*crLHS433 + DN_v(5,2)*crLHS434 - crLHS1067*crLHS438 - crLHS1069*crLHS438 + crLHS1077);
    rLHS(16,16)+=gauss_weight*(DN_v(5,0)*crLHS189 + std::pow(DN_v(5,1), 2)*crLHS9 + DN_v(5,1)*crLHS439 + DN_v(5,2)*crLHS441 - crLHS1067*crLHS446 - crLHS1069*crLHS446 + crLHS1075);
    rLHS(16,17)+=gauss_weight*(DN_v(5,0)*crLHS198 + DN_v(5,1)*crLHS447 + DN_v(5,2)*crLHS449 - crLHS1067*crLHS451 - crLHS1069*crLHS451 + crLHS1104);
    rLHS(16,18)+=gauss_weight*(DN_v(5,0)*crLHS205 + DN_v(5,1)*crLHS452 + DN_v(5,2)*crLHS453 - crLHS1067*crLHS457 - crLHS1069*crLHS457 + crLHS1105);
    rLHS(16,19)+=gauss_weight*(DN_v(5,0)*crLHS221 + DN_v(5,1)*crLHS458 + DN_v(5,2)*crLHS460 - crLHS1067*crLHS465 - crLHS1069*crLHS465 + crLHS1081 + crLHS1106);
    rLHS(16,20)+=gauss_weight*(DN_v(5,0)*crLHS230 + DN_v(5,1)*crLHS466 + DN_v(5,2)*crLHS468 - crLHS1067*crLHS470 - crLHS1069*crLHS470 + crLHS1107);
    rLHS(16,21)+=gauss_weight*(DN_v(5,0)*crLHS237 + DN_v(5,1)*crLHS471 + DN_v(5,2)*crLHS472 - crLHS1067*crLHS476 - crLHS1069*crLHS476 + crLHS1108);
    rLHS(16,22)+=gauss_weight*(DN_v(5,0)*crLHS253 + DN_v(5,1)*crLHS477 + DN_v(5,2)*crLHS479 - crLHS1067*crLHS484 - crLHS1069*crLHS484 + crLHS1086 + crLHS1109);
    rLHS(16,23)+=gauss_weight*(DN_v(5,0)*crLHS262 + DN_v(5,1)*crLHS485 + DN_v(5,2)*crLHS487 - crLHS1067*crLHS489 - crLHS1069*crLHS489 + crLHS1110);
    rLHS(16,24)+=gauss_weight*(DN_v(5,0)*crLHS269 + DN_v(5,1)*crLHS490 + DN_v(5,2)*crLHS491 - crLHS1067*crLHS495 - crLHS1069*crLHS495 + crLHS1111);
    rLHS(16,25)+=gauss_weight*(DN_v(5,0)*crLHS285 + DN_v(5,1)*crLHS496 + DN_v(5,2)*crLHS498 - crLHS1067*crLHS503 - crLHS1069*crLHS503 + crLHS1091 + crLHS1112);
    rLHS(16,26)+=gauss_weight*(DN_v(5,0)*crLHS294 + DN_v(5,1)*crLHS504 + DN_v(5,2)*crLHS506 - crLHS1067*crLHS508 - crLHS1069*crLHS508 + crLHS1113);
    rLHS(16,27)+=gauss_weight*(DN_v(5,0)*crLHS301 + DN_v(5,1)*crLHS509 + DN_v(5,2)*crLHS510 - crLHS1067*crLHS514 - crLHS1069*crLHS514 + crLHS1114);
    rLHS(16,28)+=gauss_weight*(DN_v(5,0)*crLHS316 + DN_v(5,1)*crLHS515 + DN_v(5,2)*crLHS517 - crLHS1067*crLHS522 - crLHS1069*crLHS522 + crLHS1096 + crLHS1115);
    rLHS(16,29)+=gauss_weight*(DN_v(5,0)*crLHS325 + DN_v(5,1)*crLHS523 + DN_v(5,2)*crLHS525 - crLHS1067*crLHS527 - crLHS1069*crLHS527 + crLHS1116);
    rLHS(16,30)+=gauss_weight*(crLHS1068*crLHS529 + crLHS1117 + crLHS181*crLHS529);
    rLHS(16,31)+=gauss_weight*(crLHS1068*crLHS531 + crLHS1118 + crLHS181*crLHS531);
    rLHS(16,32)+=gauss_weight*(crLHS1068*crLHS533 + crLHS1119 + crLHS181*crLHS533);
    rLHS(16,33)+=gauss_weight*(crLHS1068*crLHS535 + crLHS1120 + crLHS181*crLHS535);
    rLHS(17,0)+=gauss_weight*(DN_v(5,0)*crLHS4 + DN_v(5,1)*crLHS339 + DN_v(5,2)*crLHS536 - crLHS1067*crLHS539 - crLHS1069*crLHS539 + crLHS201);
    rLHS(17,1)+=gauss_weight*(DN_v(5,0)*crLHS31 + DN_v(5,1)*crLHS345 + DN_v(5,2)*crLHS540 - crLHS1067*crLHS542 - crLHS1069*crLHS542 + crLHS450);
    rLHS(17,2)+=gauss_weight*(DN_v(5,0)*crLHS40 + DN_v(5,1)*crLHS353 + DN_v(5,2)*crLHS543 - crLHS1067*crLHS547 - crLHS1069*crLHS547 + crLHS1070 + crLHS615);
    rLHS(17,3)+=gauss_weight*(DN_v(5,0)*crLHS47 + DN_v(5,1)*crLHS358 + DN_v(5,2)*crLHS548 - crLHS1067*crLHS553 - crLHS1069*crLHS553 + crLHS710);
    rLHS(17,4)+=gauss_weight*(DN_v(5,0)*crLHS64 + DN_v(5,1)*crLHS365 + DN_v(5,2)*crLHS554 - crLHS1067*crLHS557 - crLHS1069*crLHS557 + crLHS748);
    rLHS(17,5)+=gauss_weight*(DN_v(5,0)*crLHS72 + DN_v(5,1)*crLHS373 + DN_v(5,2)*crLHS558 - crLHS1067*crLHS562 - crLHS1069*crLHS562 + crLHS1071 + crLHS777);
    rLHS(17,6)+=gauss_weight*(DN_v(5,0)*crLHS79 + DN_v(5,1)*crLHS377 + DN_v(5,2)*crLHS563 - crLHS1067*crLHS567 - crLHS1069*crLHS567 + crLHS817);
    rLHS(17,7)+=gauss_weight*(DN_v(5,0)*crLHS96 + DN_v(5,1)*crLHS384 + DN_v(5,2)*crLHS568 - crLHS1067*crLHS571 - crLHS1069*crLHS571 + crLHS852);
    rLHS(17,8)+=gauss_weight*(DN_v(5,0)*crLHS104 + DN_v(5,1)*crLHS392 + DN_v(5,2)*crLHS572 - crLHS1067*crLHS576 - crLHS1069*crLHS576 + crLHS1072 + crLHS878);
    rLHS(17,9)+=gauss_weight*(DN_v(5,0)*crLHS111 + DN_v(5,1)*crLHS396 + DN_v(5,2)*crLHS577 - crLHS1067*crLHS581 - crLHS1069*crLHS581 + crLHS914);
    rLHS(17,10)+=gauss_weight*(DN_v(5,0)*crLHS128 + DN_v(5,1)*crLHS403 + DN_v(5,2)*crLHS582 - crLHS1067*crLHS585 - crLHS1069*crLHS585 + crLHS946);
    rLHS(17,11)+=gauss_weight*(DN_v(5,0)*crLHS136 + DN_v(5,1)*crLHS411 + DN_v(5,2)*crLHS586 - crLHS1067*crLHS590 - crLHS1069*crLHS590 + crLHS1073 + crLHS969);
    rLHS(17,12)+=gauss_weight*(DN_v(5,0)*crLHS143 + DN_v(5,1)*crLHS415 + DN_v(5,2)*crLHS591 + crLHS1001 - crLHS1067*crLHS595 - crLHS1069*crLHS595);
    rLHS(17,13)+=gauss_weight*(DN_v(5,0)*crLHS160 + DN_v(5,1)*crLHS422 + DN_v(5,2)*crLHS596 + crLHS1030 - crLHS1067*crLHS599 - crLHS1069*crLHS599);
    rLHS(17,14)+=gauss_weight*(DN_v(5,0)*crLHS168 + DN_v(5,1)*crLHS430 + DN_v(5,2)*crLHS600 + crLHS1050 - crLHS1067*crLHS604 - crLHS1069*crLHS604 + crLHS1074);
    rLHS(17,15)+=gauss_weight*(DN_v(5,0)*crLHS175 + DN_v(5,1)*crLHS434 + DN_v(5,2)*crLHS605 - crLHS1067*crLHS609 - crLHS1069*crLHS609 + crLHS1078);
    rLHS(17,16)+=gauss_weight*(DN_v(5,0)*crLHS192 + DN_v(5,1)*crLHS441 + DN_v(5,2)*crLHS610 - crLHS1067*crLHS613 - crLHS1069*crLHS613 + crLHS1104);
    rLHS(17,17)+=gauss_weight*(DN_v(5,0)*crLHS200 + DN_v(5,1)*crLHS449 + std::pow(DN_v(5,2), 2)*crLHS9 + DN_v(5,2)*crLHS614 - crLHS1067*crLHS618 - crLHS1069*crLHS618 + crLHS1075);
    rLHS(17,18)+=gauss_weight*(DN_v(5,0)*crLHS207 + DN_v(5,1)*crLHS453 + DN_v(5,2)*crLHS619 - crLHS1067*crLHS623 - crLHS1069*crLHS623 + crLHS1122);
    rLHS(17,19)+=gauss_weight*(DN_v(5,0)*crLHS224 + DN_v(5,1)*crLHS460 + DN_v(5,2)*crLHS624 - crLHS1067*crLHS627 - crLHS1069*crLHS627 + crLHS1123);
    rLHS(17,20)+=gauss_weight*(DN_v(5,0)*crLHS232 + DN_v(5,1)*crLHS468 + DN_v(5,2)*crLHS628 - crLHS1067*crLHS632 - crLHS1069*crLHS632 + crLHS1081 + crLHS1124);
    rLHS(17,21)+=gauss_weight*(DN_v(5,0)*crLHS239 + DN_v(5,1)*crLHS472 + DN_v(5,2)*crLHS633 - crLHS1067*crLHS637 - crLHS1069*crLHS637 + crLHS1125);
    rLHS(17,22)+=gauss_weight*(DN_v(5,0)*crLHS256 + DN_v(5,1)*crLHS479 + DN_v(5,2)*crLHS638 - crLHS1067*crLHS641 - crLHS1069*crLHS641 + crLHS1126);
    rLHS(17,23)+=gauss_weight*(DN_v(5,0)*crLHS264 + DN_v(5,1)*crLHS487 + DN_v(5,2)*crLHS642 - crLHS1067*crLHS646 - crLHS1069*crLHS646 + crLHS1086 + crLHS1127);
    rLHS(17,24)+=gauss_weight*(DN_v(5,0)*crLHS271 + DN_v(5,1)*crLHS491 + DN_v(5,2)*crLHS647 - crLHS1067*crLHS651 - crLHS1069*crLHS651 + crLHS1128);
    rLHS(17,25)+=gauss_weight*(DN_v(5,0)*crLHS288 + DN_v(5,1)*crLHS498 + DN_v(5,2)*crLHS652 - crLHS1067*crLHS655 - crLHS1069*crLHS655 + crLHS1129);
    rLHS(17,26)+=gauss_weight*(DN_v(5,0)*crLHS296 + DN_v(5,1)*crLHS506 + DN_v(5,2)*crLHS656 - crLHS1067*crLHS660 - crLHS1069*crLHS660 + crLHS1091 + crLHS1130);
    rLHS(17,27)+=gauss_weight*(DN_v(5,0)*crLHS303 + DN_v(5,1)*crLHS510 + DN_v(5,2)*crLHS661 - crLHS1067*crLHS665 - crLHS1069*crLHS665 + crLHS1131);
    rLHS(17,28)+=gauss_weight*(DN_v(5,0)*crLHS319 + DN_v(5,1)*crLHS517 + DN_v(5,2)*crLHS666 - crLHS1067*crLHS669 - crLHS1069*crLHS669 + crLHS1132);
    rLHS(17,29)+=gauss_weight*(DN_v(5,0)*crLHS327 + DN_v(5,1)*crLHS525 + DN_v(5,2)*crLHS670 - crLHS1067*crLHS674 - crLHS1069*crLHS674 + crLHS1096 + crLHS1133);
    rLHS(17,30)+=gauss_weight*(crLHS1068*crLHS676 + crLHS1134 + crLHS181*crLHS676);
    rLHS(17,31)+=gauss_weight*(crLHS1068*crLHS678 + crLHS1135 + crLHS181*crLHS678);
    rLHS(17,32)+=gauss_weight*(crLHS1068*crLHS680 + crLHS1136 + crLHS181*crLHS680);
    rLHS(17,33)+=gauss_weight*(crLHS1068*crLHS682 + crLHS1137 + crLHS181*crLHS682);
    rLHS(18,0)+=gauss_weight*(DN_v(6,0)*crLHS0 + DN_v(6,1)*crLHS2 + DN_v(6,2)*crLHS4 - crLHS1138*crLHS19 - crLHS1140*crLHS19 + crLHS1141 + crLHS208);
    rLHS(18,1)+=gauss_weight*(DN_v(6,0)*crLHS26 + DN_v(6,1)*crLHS28 + DN_v(6,2)*crLHS31 - crLHS1138*crLHS35 - crLHS1140*crLHS35 + crLHS454);
    rLHS(18,2)+=gauss_weight*(DN_v(6,0)*crLHS36 + DN_v(6,1)*crLHS38 + DN_v(6,2)*crLHS40 - crLHS1138*crLHS42 - crLHS1140*crLHS42 + crLHS620);
    rLHS(18,3)+=gauss_weight*(DN_v(6,0)*crLHS43 + DN_v(6,1)*crLHS45 + DN_v(6,2)*crLHS47 - crLHS1138*crLHS56 - crLHS1140*crLHS56 + crLHS1142 + crLHS711);
    rLHS(18,4)+=gauss_weight*(DN_v(6,0)*crLHS59 + DN_v(6,1)*crLHS61 + DN_v(6,2)*crLHS64 - crLHS1138*crLHS67 - crLHS1140*crLHS67 + crLHS749);
    rLHS(18,5)+=gauss_weight*(DN_v(6,0)*crLHS68 + DN_v(6,1)*crLHS70 + DN_v(6,2)*crLHS72 - crLHS1138*crLHS74 - crLHS1140*crLHS74 + crLHS778);
    rLHS(18,6)+=gauss_weight*(DN_v(6,0)*crLHS75 + DN_v(6,1)*crLHS77 + DN_v(6,2)*crLHS79 - crLHS1138*crLHS88 - crLHS1140*crLHS88 + crLHS1143 + crLHS818);
    rLHS(18,7)+=gauss_weight*(DN_v(6,0)*crLHS91 + DN_v(6,1)*crLHS93 + DN_v(6,2)*crLHS96 - crLHS1138*crLHS99 - crLHS1140*crLHS99 + crLHS853);
    rLHS(18,8)+=gauss_weight*(DN_v(6,0)*crLHS100 + DN_v(6,1)*crLHS102 + DN_v(6,2)*crLHS104 - crLHS106*crLHS1138 - crLHS106*crLHS1140 + crLHS879);
    rLHS(18,9)+=gauss_weight*(DN_v(6,0)*crLHS107 + DN_v(6,1)*crLHS109 + DN_v(6,2)*crLHS111 - crLHS1138*crLHS120 - crLHS1140*crLHS120 + crLHS1144 + crLHS915);
    rLHS(18,10)+=gauss_weight*(DN_v(6,0)*crLHS123 + DN_v(6,1)*crLHS125 + DN_v(6,2)*crLHS128 - crLHS1138*crLHS131 - crLHS1140*crLHS131 + crLHS947);
    rLHS(18,11)+=gauss_weight*(DN_v(6,0)*crLHS132 + DN_v(6,1)*crLHS134 + DN_v(6,2)*crLHS136 - crLHS1138*crLHS138 - crLHS1140*crLHS138 + crLHS970);
    rLHS(18,12)+=gauss_weight*(DN_v(6,0)*crLHS139 + DN_v(6,1)*crLHS141 + DN_v(6,2)*crLHS143 + crLHS1002 - crLHS1138*crLHS152 - crLHS1140*crLHS152 + crLHS1145);
    rLHS(18,13)+=gauss_weight*(DN_v(6,0)*crLHS155 + DN_v(6,1)*crLHS157 + DN_v(6,2)*crLHS160 + crLHS1031 - crLHS1138*crLHS163 - crLHS1140*crLHS163);
    rLHS(18,14)+=gauss_weight*(DN_v(6,0)*crLHS164 + DN_v(6,1)*crLHS166 + DN_v(6,2)*crLHS168 + crLHS1051 - crLHS1138*crLHS170 - crLHS1140*crLHS170);
    rLHS(18,15)+=gauss_weight*(DN_v(6,0)*crLHS171 + DN_v(6,1)*crLHS173 + DN_v(6,2)*crLHS175 + crLHS1079 - crLHS1138*crLHS184 - crLHS1140*crLHS184 + crLHS1146);
    rLHS(18,16)+=gauss_weight*(DN_v(6,0)*crLHS187 + DN_v(6,1)*crLHS189 + DN_v(6,2)*crLHS192 + crLHS1105 - crLHS1138*crLHS195 - crLHS1140*crLHS195);
    rLHS(18,17)+=gauss_weight*(DN_v(6,0)*crLHS196 + DN_v(6,1)*crLHS198 + DN_v(6,2)*crLHS200 + crLHS1122 - crLHS1138*crLHS202 - crLHS1140*crLHS202);
    rLHS(18,18)+=gauss_weight*(std::pow(DN_v(6,0), 2)*crLHS9 + DN_v(6,0)*crLHS203 + DN_v(6,1)*crLHS205 + DN_v(6,2)*crLHS207 - crLHS1138*crLHS216 - crLHS1140*crLHS216 + crLHS1147);
    rLHS(18,19)+=gauss_weight*(DN_v(6,0)*crLHS219 + DN_v(6,1)*crLHS221 + DN_v(6,2)*crLHS224 - crLHS1138*crLHS227 - crLHS1140*crLHS227 + crLHS1149);
    rLHS(18,20)+=gauss_weight*(DN_v(6,0)*crLHS228 + DN_v(6,1)*crLHS230 + DN_v(6,2)*crLHS232 - crLHS1138*crLHS234 - crLHS1140*crLHS234 + crLHS1150);
    rLHS(18,21)+=gauss_weight*(DN_v(6,0)*crLHS235 + DN_v(6,1)*crLHS237 + DN_v(6,2)*crLHS239 - crLHS1138*crLHS248 - crLHS1140*crLHS248 + crLHS1151 + crLHS1153);
    rLHS(18,22)+=gauss_weight*(DN_v(6,0)*crLHS251 + DN_v(6,1)*crLHS253 + DN_v(6,2)*crLHS256 - crLHS1138*crLHS259 - crLHS1140*crLHS259 + crLHS1154);
    rLHS(18,23)+=gauss_weight*(DN_v(6,0)*crLHS260 + DN_v(6,1)*crLHS262 + DN_v(6,2)*crLHS264 - crLHS1138*crLHS266 - crLHS1140*crLHS266 + crLHS1155);
    rLHS(18,24)+=gauss_weight*(DN_v(6,0)*crLHS267 + DN_v(6,1)*crLHS269 + DN_v(6,2)*crLHS271 - crLHS1138*crLHS280 - crLHS1140*crLHS280 + crLHS1156 + crLHS1158);
    rLHS(18,25)+=gauss_weight*(DN_v(6,0)*crLHS283 + DN_v(6,1)*crLHS285 + DN_v(6,2)*crLHS288 - crLHS1138*crLHS291 - crLHS1140*crLHS291 + crLHS1159);
    rLHS(18,26)+=gauss_weight*(DN_v(6,0)*crLHS292 + DN_v(6,1)*crLHS294 + DN_v(6,2)*crLHS296 - crLHS1138*crLHS298 - crLHS1140*crLHS298 + crLHS1160);
    rLHS(18,27)+=gauss_weight*(DN_v(6,0)*crLHS299 + DN_v(6,1)*crLHS301 + DN_v(6,2)*crLHS303 - crLHS1138*crLHS311 - crLHS1140*crLHS311 + crLHS1161 + crLHS1163);
    rLHS(18,28)+=gauss_weight*(DN_v(6,0)*crLHS314 + DN_v(6,1)*crLHS316 + DN_v(6,2)*crLHS319 - crLHS1138*crLHS322 - crLHS1140*crLHS322 + crLHS1164);
    rLHS(18,29)+=gauss_weight*(DN_v(6,0)*crLHS323 + DN_v(6,1)*crLHS325 + DN_v(6,2)*crLHS327 - crLHS1138*crLHS329 - crLHS1140*crLHS329 + crLHS1165);
    rLHS(18,30)+=gauss_weight*(crLHS1139*crLHS331 + crLHS1166 + crLHS213*crLHS331);
    rLHS(18,31)+=gauss_weight*(crLHS1139*crLHS333 + crLHS1167 + crLHS213*crLHS333);
    rLHS(18,32)+=gauss_weight*(crLHS1139*crLHS335 + crLHS1168 + crLHS213*crLHS335);
    rLHS(18,33)+=gauss_weight*(crLHS1139*crLHS337 + crLHS1169 + crLHS213*crLHS337);
    rLHS(19,0)+=gauss_weight*(DN_v(6,0)*crLHS2 + DN_v(6,1)*crLHS338 + DN_v(6,2)*crLHS339 - crLHS1138*crLHS342 - crLHS1140*crLHS342 + crLHS225);
    rLHS(19,1)+=gauss_weight*(DN_v(6,0)*crLHS28 + DN_v(6,1)*crLHS343 + DN_v(6,2)*crLHS345 - crLHS1138*crLHS350 - crLHS1140*crLHS350 + crLHS1141 + crLHS461);
    rLHS(19,2)+=gauss_weight*(DN_v(6,0)*crLHS38 + DN_v(6,1)*crLHS351 + DN_v(6,2)*crLHS353 - crLHS1138*crLHS356 - crLHS1140*crLHS356 + crLHS625);
    rLHS(19,3)+=gauss_weight*(DN_v(6,0)*crLHS45 + DN_v(6,1)*crLHS357 + DN_v(6,2)*crLHS358 - crLHS1138*crLHS362 - crLHS1140*crLHS362 + crLHS714);
    rLHS(19,4)+=gauss_weight*(DN_v(6,0)*crLHS61 + DN_v(6,1)*crLHS363 + DN_v(6,2)*crLHS365 - crLHS1138*crLHS370 - crLHS1140*crLHS370 + crLHS1142 + crLHS750);
    rLHS(19,5)+=gauss_weight*(DN_v(6,0)*crLHS70 + DN_v(6,1)*crLHS371 + DN_v(6,2)*crLHS373 - crLHS1138*crLHS375 - crLHS1140*crLHS375 + crLHS779);
    rLHS(19,6)+=gauss_weight*(DN_v(6,0)*crLHS77 + DN_v(6,1)*crLHS376 + DN_v(6,2)*crLHS377 - crLHS1138*crLHS381 - crLHS1140*crLHS381 + crLHS821);
    rLHS(19,7)+=gauss_weight*(DN_v(6,0)*crLHS93 + DN_v(6,1)*crLHS382 + DN_v(6,2)*crLHS384 - crLHS1138*crLHS389 - crLHS1140*crLHS389 + crLHS1143 + crLHS854);
    rLHS(19,8)+=gauss_weight*(DN_v(6,0)*crLHS102 + DN_v(6,1)*crLHS390 + DN_v(6,2)*crLHS392 - crLHS1138*crLHS394 - crLHS1140*crLHS394 + crLHS880);
    rLHS(19,9)+=gauss_weight*(DN_v(6,0)*crLHS109 + DN_v(6,1)*crLHS395 + DN_v(6,2)*crLHS396 - crLHS1138*crLHS400 - crLHS1140*crLHS400 + crLHS918);
    rLHS(19,10)+=gauss_weight*(DN_v(6,0)*crLHS125 + DN_v(6,1)*crLHS401 + DN_v(6,2)*crLHS403 - crLHS1138*crLHS408 - crLHS1140*crLHS408 + crLHS1144 + crLHS948);
    rLHS(19,11)+=gauss_weight*(DN_v(6,0)*crLHS134 + DN_v(6,1)*crLHS409 + DN_v(6,2)*crLHS411 - crLHS1138*crLHS413 - crLHS1140*crLHS413 + crLHS971);
    rLHS(19,12)+=gauss_weight*(DN_v(6,0)*crLHS141 + DN_v(6,1)*crLHS414 + DN_v(6,2)*crLHS415 + crLHS1005 - crLHS1138*crLHS419 - crLHS1140*crLHS419);
    rLHS(19,13)+=gauss_weight*(DN_v(6,0)*crLHS157 + DN_v(6,1)*crLHS420 + DN_v(6,2)*crLHS422 + crLHS1032 - crLHS1138*crLHS427 - crLHS1140*crLHS427 + crLHS1145);
    rLHS(19,14)+=gauss_weight*(DN_v(6,0)*crLHS166 + DN_v(6,1)*crLHS428 + DN_v(6,2)*crLHS430 + crLHS1052 - crLHS1138*crLHS432 - crLHS1140*crLHS432);
    rLHS(19,15)+=gauss_weight*(DN_v(6,0)*crLHS173 + DN_v(6,1)*crLHS433 + DN_v(6,2)*crLHS434 + crLHS1082 - crLHS1138*crLHS438 - crLHS1140*crLHS438);
    rLHS(19,16)+=gauss_weight*(DN_v(6,0)*crLHS189 + DN_v(6,1)*crLHS439 + DN_v(6,2)*crLHS441 + crLHS1106 - crLHS1138*crLHS446 - crLHS1140*crLHS446 + crLHS1146);
    rLHS(19,17)+=gauss_weight*(DN_v(6,0)*crLHS198 + DN_v(6,1)*crLHS447 + DN_v(6,2)*crLHS449 + crLHS1123 - crLHS1138*crLHS451 - crLHS1140*crLHS451);
    rLHS(19,18)+=gauss_weight*(DN_v(6,0)*crLHS205 + DN_v(6,1)*crLHS452 + DN_v(6,2)*crLHS453 - crLHS1138*crLHS457 - crLHS1140*crLHS457 + crLHS1149);
    rLHS(19,19)+=gauss_weight*(DN_v(6,0)*crLHS221 + std::pow(DN_v(6,1), 2)*crLHS9 + DN_v(6,1)*crLHS458 + DN_v(6,2)*crLHS460 - crLHS1138*crLHS465 - crLHS1140*crLHS465 + crLHS1147);
    rLHS(19,20)+=gauss_weight*(DN_v(6,0)*crLHS230 + DN_v(6,1)*crLHS466 + DN_v(6,2)*crLHS468 - crLHS1138*crLHS470 - crLHS1140*crLHS470 + crLHS1171);
    rLHS(19,21)+=gauss_weight*(DN_v(6,0)*crLHS237 + DN_v(6,1)*crLHS471 + DN_v(6,2)*crLHS472 - crLHS1138*crLHS476 - crLHS1140*crLHS476 + crLHS1172);
    rLHS(19,22)+=gauss_weight*(DN_v(6,0)*crLHS253 + DN_v(6,1)*crLHS477 + DN_v(6,2)*crLHS479 - crLHS1138*crLHS484 - crLHS1140*crLHS484 + crLHS1153 + crLHS1173);
    rLHS(19,23)+=gauss_weight*(DN_v(6,0)*crLHS262 + DN_v(6,1)*crLHS485 + DN_v(6,2)*crLHS487 - crLHS1138*crLHS489 - crLHS1140*crLHS489 + crLHS1174);
    rLHS(19,24)+=gauss_weight*(DN_v(6,0)*crLHS269 + DN_v(6,1)*crLHS490 + DN_v(6,2)*crLHS491 - crLHS1138*crLHS495 - crLHS1140*crLHS495 + crLHS1175);
    rLHS(19,25)+=gauss_weight*(DN_v(6,0)*crLHS285 + DN_v(6,1)*crLHS496 + DN_v(6,2)*crLHS498 - crLHS1138*crLHS503 - crLHS1140*crLHS503 + crLHS1158 + crLHS1176);
    rLHS(19,26)+=gauss_weight*(DN_v(6,0)*crLHS294 + DN_v(6,1)*crLHS504 + DN_v(6,2)*crLHS506 - crLHS1138*crLHS508 - crLHS1140*crLHS508 + crLHS1177);
    rLHS(19,27)+=gauss_weight*(DN_v(6,0)*crLHS301 + DN_v(6,1)*crLHS509 + DN_v(6,2)*crLHS510 - crLHS1138*crLHS514 - crLHS1140*crLHS514 + crLHS1178);
    rLHS(19,28)+=gauss_weight*(DN_v(6,0)*crLHS316 + DN_v(6,1)*crLHS515 + DN_v(6,2)*crLHS517 - crLHS1138*crLHS522 - crLHS1140*crLHS522 + crLHS1163 + crLHS1179);
    rLHS(19,29)+=gauss_weight*(DN_v(6,0)*crLHS325 + DN_v(6,1)*crLHS523 + DN_v(6,2)*crLHS525 - crLHS1138*crLHS527 - crLHS1140*crLHS527 + crLHS1180);
    rLHS(19,30)+=gauss_weight*(crLHS1139*crLHS529 + crLHS1181 + crLHS213*crLHS529);
    rLHS(19,31)+=gauss_weight*(crLHS1139*crLHS531 + crLHS1182 + crLHS213*crLHS531);
    rLHS(19,32)+=gauss_weight*(crLHS1139*crLHS533 + crLHS1183 + crLHS213*crLHS533);
    rLHS(19,33)+=gauss_weight*(crLHS1139*crLHS535 + crLHS1184 + crLHS213*crLHS535);
    rLHS(20,0)+=gauss_weight*(DN_v(6,0)*crLHS4 + DN_v(6,1)*crLHS339 + DN_v(6,2)*crLHS536 - crLHS1138*crLHS539 - crLHS1140*crLHS539 + crLHS233);
    rLHS(20,1)+=gauss_weight*(DN_v(6,0)*crLHS31 + DN_v(6,1)*crLHS345 + DN_v(6,2)*crLHS540 - crLHS1138*crLHS542 - crLHS1140*crLHS542 + crLHS469);
    rLHS(20,2)+=gauss_weight*(DN_v(6,0)*crLHS40 + DN_v(6,1)*crLHS353 + DN_v(6,2)*crLHS543 - crLHS1138*crLHS547 - crLHS1140*crLHS547 + crLHS1141 + crLHS629);
    rLHS(20,3)+=gauss_weight*(DN_v(6,0)*crLHS47 + DN_v(6,1)*crLHS358 + DN_v(6,2)*crLHS548 - crLHS1138*crLHS553 - crLHS1140*crLHS553 + crLHS715);
    rLHS(20,4)+=gauss_weight*(DN_v(6,0)*crLHS64 + DN_v(6,1)*crLHS365 + DN_v(6,2)*crLHS554 - crLHS1138*crLHS557 - crLHS1140*crLHS557 + crLHS751);
    rLHS(20,5)+=gauss_weight*(DN_v(6,0)*crLHS72 + DN_v(6,1)*crLHS373 + DN_v(6,2)*crLHS558 - crLHS1138*crLHS562 - crLHS1140*crLHS562 + crLHS1142 + crLHS780);
    rLHS(20,6)+=gauss_weight*(DN_v(6,0)*crLHS79 + DN_v(6,1)*crLHS377 + DN_v(6,2)*crLHS563 - crLHS1138*crLHS567 - crLHS1140*crLHS567 + crLHS822);
    rLHS(20,7)+=gauss_weight*(DN_v(6,0)*crLHS96 + DN_v(6,1)*crLHS384 + DN_v(6,2)*crLHS568 - crLHS1138*crLHS571 - crLHS1140*crLHS571 + crLHS855);
    rLHS(20,8)+=gauss_weight*(DN_v(6,0)*crLHS104 + DN_v(6,1)*crLHS392 + DN_v(6,2)*crLHS572 - crLHS1138*crLHS576 - crLHS1140*crLHS576 + crLHS1143 + crLHS881);
    rLHS(20,9)+=gauss_weight*(DN_v(6,0)*crLHS111 + DN_v(6,1)*crLHS396 + DN_v(6,2)*crLHS577 - crLHS1138*crLHS581 - crLHS1140*crLHS581 + crLHS919);
    rLHS(20,10)+=gauss_weight*(DN_v(6,0)*crLHS128 + DN_v(6,1)*crLHS403 + DN_v(6,2)*crLHS582 - crLHS1138*crLHS585 - crLHS1140*crLHS585 + crLHS949);
    rLHS(20,11)+=gauss_weight*(DN_v(6,0)*crLHS136 + DN_v(6,1)*crLHS411 + DN_v(6,2)*crLHS586 - crLHS1138*crLHS590 - crLHS1140*crLHS590 + crLHS1144 + crLHS972);
    rLHS(20,12)+=gauss_weight*(DN_v(6,0)*crLHS143 + DN_v(6,1)*crLHS415 + DN_v(6,2)*crLHS591 + crLHS1006 - crLHS1138*crLHS595 - crLHS1140*crLHS595);
    rLHS(20,13)+=gauss_weight*(DN_v(6,0)*crLHS160 + DN_v(6,1)*crLHS422 + DN_v(6,2)*crLHS596 + crLHS1033 - crLHS1138*crLHS599 - crLHS1140*crLHS599);
    rLHS(20,14)+=gauss_weight*(DN_v(6,0)*crLHS168 + DN_v(6,1)*crLHS430 + DN_v(6,2)*crLHS600 + crLHS1053 - crLHS1138*crLHS604 - crLHS1140*crLHS604 + crLHS1145);
    rLHS(20,15)+=gauss_weight*(DN_v(6,0)*crLHS175 + DN_v(6,1)*crLHS434 + DN_v(6,2)*crLHS605 + crLHS1083 - crLHS1138*crLHS609 - crLHS1140*crLHS609);
    rLHS(20,16)+=gauss_weight*(DN_v(6,0)*crLHS192 + DN_v(6,1)*crLHS441 + DN_v(6,2)*crLHS610 + crLHS1107 - crLHS1138*crLHS613 - crLHS1140*crLHS613);
    rLHS(20,17)+=gauss_weight*(DN_v(6,0)*crLHS200 + DN_v(6,1)*crLHS449 + DN_v(6,2)*crLHS614 + crLHS1124 - crLHS1138*crLHS618 - crLHS1140*crLHS618 + crLHS1146);
    rLHS(20,18)+=gauss_weight*(DN_v(6,0)*crLHS207 + DN_v(6,1)*crLHS453 + DN_v(6,2)*crLHS619 - crLHS1138*crLHS623 - crLHS1140*crLHS623 + crLHS1150);
    rLHS(20,19)+=gauss_weight*(DN_v(6,0)*crLHS224 + DN_v(6,1)*crLHS460 + DN_v(6,2)*crLHS624 - crLHS1138*crLHS627 - crLHS1140*crLHS627 + crLHS1171);
    rLHS(20,20)+=gauss_weight*(DN_v(6,0)*crLHS232 + DN_v(6,1)*crLHS468 + std::pow(DN_v(6,2), 2)*crLHS9 + DN_v(6,2)*crLHS628 - crLHS1138*crLHS632 - crLHS1140*crLHS632 + crLHS1147);
    rLHS(20,21)+=gauss_weight*(DN_v(6,0)*crLHS239 + DN_v(6,1)*crLHS472 + DN_v(6,2)*crLHS633 - crLHS1138*crLHS637 - crLHS1140*crLHS637 + crLHS1186);
    rLHS(20,22)+=gauss_weight*(DN_v(6,0)*crLHS256 + DN_v(6,1)*crLHS479 + DN_v(6,2)*crLHS638 - crLHS1138*crLHS641 - crLHS1140*crLHS641 + crLHS1187);
    rLHS(20,23)+=gauss_weight*(DN_v(6,0)*crLHS264 + DN_v(6,1)*crLHS487 + DN_v(6,2)*crLHS642 - crLHS1138*crLHS646 - crLHS1140*crLHS646 + crLHS1153 + crLHS1188);
    rLHS(20,24)+=gauss_weight*(DN_v(6,0)*crLHS271 + DN_v(6,1)*crLHS491 + DN_v(6,2)*crLHS647 - crLHS1138*crLHS651 - crLHS1140*crLHS651 + crLHS1189);
    rLHS(20,25)+=gauss_weight*(DN_v(6,0)*crLHS288 + DN_v(6,1)*crLHS498 + DN_v(6,2)*crLHS652 - crLHS1138*crLHS655 - crLHS1140*crLHS655 + crLHS1190);
    rLHS(20,26)+=gauss_weight*(DN_v(6,0)*crLHS296 + DN_v(6,1)*crLHS506 + DN_v(6,2)*crLHS656 - crLHS1138*crLHS660 - crLHS1140*crLHS660 + crLHS1158 + crLHS1191);
    rLHS(20,27)+=gauss_weight*(DN_v(6,0)*crLHS303 + DN_v(6,1)*crLHS510 + DN_v(6,2)*crLHS661 - crLHS1138*crLHS665 - crLHS1140*crLHS665 + crLHS1192);
    rLHS(20,28)+=gauss_weight*(DN_v(6,0)*crLHS319 + DN_v(6,1)*crLHS517 + DN_v(6,2)*crLHS666 - crLHS1138*crLHS669 - crLHS1140*crLHS669 + crLHS1193);
    rLHS(20,29)+=gauss_weight*(DN_v(6,0)*crLHS327 + DN_v(6,1)*crLHS525 + DN_v(6,2)*crLHS670 - crLHS1138*crLHS674 - crLHS1140*crLHS674 + crLHS1163 + crLHS1194);
    rLHS(20,30)+=gauss_weight*(crLHS1139*crLHS676 + crLHS1195 + crLHS213*crLHS676);
    rLHS(20,31)+=gauss_weight*(crLHS1139*crLHS678 + crLHS1196 + crLHS213*crLHS678);
    rLHS(20,32)+=gauss_weight*(crLHS1139*crLHS680 + crLHS1197 + crLHS213*crLHS680);
    rLHS(20,33)+=gauss_weight*(crLHS1139*crLHS682 + crLHS1198 + crLHS213*crLHS682);
    rLHS(21,0)+=gauss_weight*(DN_v(7,0)*crLHS0 + DN_v(7,1)*crLHS2 + DN_v(7,2)*crLHS4 - crLHS1199*crLHS19 - crLHS1201*crLHS19 + crLHS1202 + crLHS240);
    rLHS(21,1)+=gauss_weight*(DN_v(7,0)*crLHS26 + DN_v(7,1)*crLHS28 + DN_v(7,2)*crLHS31 - crLHS1199*crLHS35 - crLHS1201*crLHS35 + crLHS473);
    rLHS(21,2)+=gauss_weight*(DN_v(7,0)*crLHS36 + DN_v(7,1)*crLHS38 + DN_v(7,2)*crLHS40 - crLHS1199*crLHS42 - crLHS1201*crLHS42 + crLHS634);
    rLHS(21,3)+=gauss_weight*(DN_v(7,0)*crLHS43 + DN_v(7,1)*crLHS45 + DN_v(7,2)*crLHS47 - crLHS1199*crLHS56 - crLHS1201*crLHS56 + crLHS1203 + crLHS716);
    rLHS(21,4)+=gauss_weight*(DN_v(7,0)*crLHS59 + DN_v(7,1)*crLHS61 + DN_v(7,2)*crLHS64 - crLHS1199*crLHS67 - crLHS1201*crLHS67 + crLHS752);
    rLHS(21,5)+=gauss_weight*(DN_v(7,0)*crLHS68 + DN_v(7,1)*crLHS70 + DN_v(7,2)*crLHS72 - crLHS1199*crLHS74 - crLHS1201*crLHS74 + crLHS781);
    rLHS(21,6)+=gauss_weight*(DN_v(7,0)*crLHS75 + DN_v(7,1)*crLHS77 + DN_v(7,2)*crLHS79 - crLHS1199*crLHS88 - crLHS1201*crLHS88 + crLHS1204 + crLHS823);
    rLHS(21,7)+=gauss_weight*(DN_v(7,0)*crLHS91 + DN_v(7,1)*crLHS93 + DN_v(7,2)*crLHS96 - crLHS1199*crLHS99 - crLHS1201*crLHS99 + crLHS856);
    rLHS(21,8)+=gauss_weight*(DN_v(7,0)*crLHS100 + DN_v(7,1)*crLHS102 + DN_v(7,2)*crLHS104 - crLHS106*crLHS1199 - crLHS106*crLHS1201 + crLHS882);
    rLHS(21,9)+=gauss_weight*(DN_v(7,0)*crLHS107 + DN_v(7,1)*crLHS109 + DN_v(7,2)*crLHS111 - crLHS1199*crLHS120 - crLHS120*crLHS1201 + crLHS1205 + crLHS920);
    rLHS(21,10)+=gauss_weight*(DN_v(7,0)*crLHS123 + DN_v(7,1)*crLHS125 + DN_v(7,2)*crLHS128 - crLHS1199*crLHS131 - crLHS1201*crLHS131 + crLHS950);
    rLHS(21,11)+=gauss_weight*(DN_v(7,0)*crLHS132 + DN_v(7,1)*crLHS134 + DN_v(7,2)*crLHS136 - crLHS1199*crLHS138 - crLHS1201*crLHS138 + crLHS973);
    rLHS(21,12)+=gauss_weight*(DN_v(7,0)*crLHS139 + DN_v(7,1)*crLHS141 + DN_v(7,2)*crLHS143 + crLHS1007 - crLHS1199*crLHS152 - crLHS1201*crLHS152 + crLHS1206);
    rLHS(21,13)+=gauss_weight*(DN_v(7,0)*crLHS155 + DN_v(7,1)*crLHS157 + DN_v(7,2)*crLHS160 + crLHS1034 - crLHS1199*crLHS163 - crLHS1201*crLHS163);
    rLHS(21,14)+=gauss_weight*(DN_v(7,0)*crLHS164 + DN_v(7,1)*crLHS166 + DN_v(7,2)*crLHS168 + crLHS1054 - crLHS1199*crLHS170 - crLHS1201*crLHS170);
    rLHS(21,15)+=gauss_weight*(DN_v(7,0)*crLHS171 + DN_v(7,1)*crLHS173 + DN_v(7,2)*crLHS175 + crLHS1084 - crLHS1199*crLHS184 - crLHS1201*crLHS184 + crLHS1207);
    rLHS(21,16)+=gauss_weight*(DN_v(7,0)*crLHS187 + DN_v(7,1)*crLHS189 + DN_v(7,2)*crLHS192 + crLHS1108 - crLHS1199*crLHS195 - crLHS1201*crLHS195);
    rLHS(21,17)+=gauss_weight*(DN_v(7,0)*crLHS196 + DN_v(7,1)*crLHS198 + DN_v(7,2)*crLHS200 + crLHS1125 - crLHS1199*crLHS202 - crLHS1201*crLHS202);
    rLHS(21,18)+=gauss_weight*(DN_v(7,0)*crLHS203 + DN_v(7,1)*crLHS205 + DN_v(7,2)*crLHS207 + crLHS1151 - crLHS1199*crLHS216 - crLHS1201*crLHS216 + crLHS1208);
    rLHS(21,19)+=gauss_weight*(DN_v(7,0)*crLHS219 + DN_v(7,1)*crLHS221 + DN_v(7,2)*crLHS224 + crLHS1172 - crLHS1199*crLHS227 - crLHS1201*crLHS227);
    rLHS(21,20)+=gauss_weight*(DN_v(7,0)*crLHS228 + DN_v(7,1)*crLHS230 + DN_v(7,2)*crLHS232 + crLHS1186 - crLHS1199*crLHS234 - crLHS1201*crLHS234);
    rLHS(21,21)+=gauss_weight*(std::pow(DN_v(7,0), 2)*crLHS9 + DN_v(7,0)*crLHS235 + DN_v(7,1)*crLHS237 + DN_v(7,2)*crLHS239 - crLHS1199*crLHS248 - crLHS1201*crLHS248 + crLHS1209);
    rLHS(21,22)+=gauss_weight*(DN_v(7,0)*crLHS251 + DN_v(7,1)*crLHS253 + DN_v(7,2)*crLHS256 - crLHS1199*crLHS259 - crLHS1201*crLHS259 + crLHS1211);
    rLHS(21,23)+=gauss_weight*(DN_v(7,0)*crLHS260 + DN_v(7,1)*crLHS262 + DN_v(7,2)*crLHS264 - crLHS1199*crLHS266 - crLHS1201*crLHS266 + crLHS1212);
    rLHS(21,24)+=gauss_weight*(DN_v(7,0)*crLHS267 + DN_v(7,1)*crLHS269 + DN_v(7,2)*crLHS271 - crLHS1199*crLHS280 - crLHS1201*crLHS280 + crLHS1213 + crLHS1215);
    rLHS(21,25)+=gauss_weight*(DN_v(7,0)*crLHS283 + DN_v(7,1)*crLHS285 + DN_v(7,2)*crLHS288 - crLHS1199*crLHS291 - crLHS1201*crLHS291 + crLHS1216);
    rLHS(21,26)+=gauss_weight*(DN_v(7,0)*crLHS292 + DN_v(7,1)*crLHS294 + DN_v(7,2)*crLHS296 - crLHS1199*crLHS298 - crLHS1201*crLHS298 + crLHS1217);
    rLHS(21,27)+=gauss_weight*(DN_v(7,0)*crLHS299 + DN_v(7,1)*crLHS301 + DN_v(7,2)*crLHS303 - crLHS1199*crLHS311 - crLHS1201*crLHS311 + crLHS1218 + crLHS1220);
    rLHS(21,28)+=gauss_weight*(DN_v(7,0)*crLHS314 + DN_v(7,1)*crLHS316 + DN_v(7,2)*crLHS319 - crLHS1199*crLHS322 - crLHS1201*crLHS322 + crLHS1221);
    rLHS(21,29)+=gauss_weight*(DN_v(7,0)*crLHS323 + DN_v(7,1)*crLHS325 + DN_v(7,2)*crLHS327 - crLHS1199*crLHS329 - crLHS1201*crLHS329 + crLHS1222);
    rLHS(21,30)+=gauss_weight*(crLHS1200*crLHS331 + crLHS1223 + crLHS245*crLHS331);
    rLHS(21,31)+=gauss_weight*(crLHS1200*crLHS333 + crLHS1224 + crLHS245*crLHS333);
    rLHS(21,32)+=gauss_weight*(crLHS1200*crLHS335 + crLHS1225 + crLHS245*crLHS335);
    rLHS(21,33)+=gauss_weight*(crLHS1200*crLHS337 + crLHS1226 + crLHS245*crLHS337);
    rLHS(22,0)+=gauss_weight*(DN_v(7,0)*crLHS2 + DN_v(7,1)*crLHS338 + DN_v(7,2)*crLHS339 - crLHS1199*crLHS342 - crLHS1201*crLHS342 + crLHS257);
    rLHS(22,1)+=gauss_weight*(DN_v(7,0)*crLHS28 + DN_v(7,1)*crLHS343 + DN_v(7,2)*crLHS345 - crLHS1199*crLHS350 - crLHS1201*crLHS350 + crLHS1202 + crLHS480);
    rLHS(22,2)+=gauss_weight*(DN_v(7,0)*crLHS38 + DN_v(7,1)*crLHS351 + DN_v(7,2)*crLHS353 - crLHS1199*crLHS356 - crLHS1201*crLHS356 + crLHS639);
    rLHS(22,3)+=gauss_weight*(DN_v(7,0)*crLHS45 + DN_v(7,1)*crLHS357 + DN_v(7,2)*crLHS358 - crLHS1199*crLHS362 - crLHS1201*crLHS362 + crLHS719);
    rLHS(22,4)+=gauss_weight*(DN_v(7,0)*crLHS61 + DN_v(7,1)*crLHS363 + DN_v(7,2)*crLHS365 - crLHS1199*crLHS370 - crLHS1201*crLHS370 + crLHS1203 + crLHS753);
    rLHS(22,5)+=gauss_weight*(DN_v(7,0)*crLHS70 + DN_v(7,1)*crLHS371 + DN_v(7,2)*crLHS373 - crLHS1199*crLHS375 - crLHS1201*crLHS375 + crLHS782);
    rLHS(22,6)+=gauss_weight*(DN_v(7,0)*crLHS77 + DN_v(7,1)*crLHS376 + DN_v(7,2)*crLHS377 - crLHS1199*crLHS381 - crLHS1201*crLHS381 + crLHS826);
    rLHS(22,7)+=gauss_weight*(DN_v(7,0)*crLHS93 + DN_v(7,1)*crLHS382 + DN_v(7,2)*crLHS384 - crLHS1199*crLHS389 - crLHS1201*crLHS389 + crLHS1204 + crLHS857);
    rLHS(22,8)+=gauss_weight*(DN_v(7,0)*crLHS102 + DN_v(7,1)*crLHS390 + DN_v(7,2)*crLHS392 - crLHS1199*crLHS394 - crLHS1201*crLHS394 + crLHS883);
    rLHS(22,9)+=gauss_weight*(DN_v(7,0)*crLHS109 + DN_v(7,1)*crLHS395 + DN_v(7,2)*crLHS396 - crLHS1199*crLHS400 - crLHS1201*crLHS400 + crLHS923);
    rLHS(22,10)+=gauss_weight*(DN_v(7,0)*crLHS125 + DN_v(7,1)*crLHS401 + DN_v(7,2)*crLHS403 - crLHS1199*crLHS408 - crLHS1201*crLHS408 + crLHS1205 + crLHS951);
    rLHS(22,11)+=gauss_weight*(DN_v(7,0)*crLHS134 + DN_v(7,1)*crLHS409 + DN_v(7,2)*crLHS411 - crLHS1199*crLHS413 - crLHS1201*crLHS413 + crLHS974);
    rLHS(22,12)+=gauss_weight*(DN_v(7,0)*crLHS141 + DN_v(7,1)*crLHS414 + DN_v(7,2)*crLHS415 + crLHS1010 - crLHS1199*crLHS419 - crLHS1201*crLHS419);
    rLHS(22,13)+=gauss_weight*(DN_v(7,0)*crLHS157 + DN_v(7,1)*crLHS420 + DN_v(7,2)*crLHS422 + crLHS1035 - crLHS1199*crLHS427 - crLHS1201*crLHS427 + crLHS1206);
    rLHS(22,14)+=gauss_weight*(DN_v(7,0)*crLHS166 + DN_v(7,1)*crLHS428 + DN_v(7,2)*crLHS430 + crLHS1055 - crLHS1199*crLHS432 - crLHS1201*crLHS432);
    rLHS(22,15)+=gauss_weight*(DN_v(7,0)*crLHS173 + DN_v(7,1)*crLHS433 + DN_v(7,2)*crLHS434 + crLHS1087 - crLHS1199*crLHS438 - crLHS1201*crLHS438);
    rLHS(22,16)+=gauss_weight*(DN_v(7,0)*crLHS189 + DN_v(7,1)*crLHS439 + DN_v(7,2)*crLHS441 + crLHS1109 - crLHS1199*crLHS446 - crLHS1201*crLHS446 + crLHS1207);
    rLHS(22,17)+=gauss_weight*(DN_v(7,0)*crLHS198 + DN_v(7,1)*crLHS447 + DN_v(7,2)*crLHS449 + crLHS1126 - crLHS1199*crLHS451 - crLHS1201*crLHS451);
    rLHS(22,18)+=gauss_weight*(DN_v(7,0)*crLHS205 + DN_v(7,1)*crLHS452 + DN_v(7,2)*crLHS453 + crLHS1154 - crLHS1199*crLHS457 - crLHS1201*crLHS457);
    rLHS(22,19)+=gauss_weight*(DN_v(7,0)*crLHS221 + DN_v(7,1)*crLHS458 + DN_v(7,2)*crLHS460 + crLHS1173 - crLHS1199*crLHS465 - crLHS1201*crLHS465 + crLHS1208);
    rLHS(22,20)+=gauss_weight*(DN_v(7,0)*crLHS230 + DN_v(7,1)*crLHS466 + DN_v(7,2)*crLHS468 + crLHS1187 - crLHS1199*crLHS470 - crLHS1201*crLHS470);
    rLHS(22,21)+=gauss_weight*(DN_v(7,0)*crLHS237 + DN_v(7,1)*crLHS471 + DN_v(7,2)*crLHS472 - crLHS1199*crLHS476 - crLHS1201*crLHS476 + crLHS1211);
    rLHS(22,22)+=gauss_weight*(DN_v(7,0)*crLHS253 + std::pow(DN_v(7,1), 2)*crLHS9 + DN_v(7,1)*crLHS477 + DN_v(7,2)*crLHS479 - crLHS1199*crLHS484 - crLHS1201*crLHS484 + crLHS1209);
    rLHS(22,23)+=gauss_weight*(DN_v(7,0)*crLHS262 + DN_v(7,1)*crLHS485 + DN_v(7,2)*crLHS487 - crLHS1199*crLHS489 - crLHS1201*crLHS489 + crLHS1228);
    rLHS(22,24)+=gauss_weight*(DN_v(7,0)*crLHS269 + DN_v(7,1)*crLHS490 + DN_v(7,2)*crLHS491 - crLHS1199*crLHS495 - crLHS1201*crLHS495 + crLHS1229);
    rLHS(22,25)+=gauss_weight*(DN_v(7,0)*crLHS285 + DN_v(7,1)*crLHS496 + DN_v(7,2)*crLHS498 - crLHS1199*crLHS503 - crLHS1201*crLHS503 + crLHS1215 + crLHS1230);
    rLHS(22,26)+=gauss_weight*(DN_v(7,0)*crLHS294 + DN_v(7,1)*crLHS504 + DN_v(7,2)*crLHS506 - crLHS1199*crLHS508 - crLHS1201*crLHS508 + crLHS1231);
    rLHS(22,27)+=gauss_weight*(DN_v(7,0)*crLHS301 + DN_v(7,1)*crLHS509 + DN_v(7,2)*crLHS510 - crLHS1199*crLHS514 - crLHS1201*crLHS514 + crLHS1232);
    rLHS(22,28)+=gauss_weight*(DN_v(7,0)*crLHS316 + DN_v(7,1)*crLHS515 + DN_v(7,2)*crLHS517 - crLHS1199*crLHS522 - crLHS1201*crLHS522 + crLHS1220 + crLHS1233);
    rLHS(22,29)+=gauss_weight*(DN_v(7,0)*crLHS325 + DN_v(7,1)*crLHS523 + DN_v(7,2)*crLHS525 - crLHS1199*crLHS527 - crLHS1201*crLHS527 + crLHS1234);
    rLHS(22,30)+=gauss_weight*(crLHS1200*crLHS529 + crLHS1235 + crLHS245*crLHS529);
    rLHS(22,31)+=gauss_weight*(crLHS1200*crLHS531 + crLHS1236 + crLHS245*crLHS531);
    rLHS(22,32)+=gauss_weight*(crLHS1200*crLHS533 + crLHS1237 + crLHS245*crLHS533);
    rLHS(22,33)+=gauss_weight*(crLHS1200*crLHS535 + crLHS1238 + crLHS245*crLHS535);
    rLHS(23,0)+=gauss_weight*(DN_v(7,0)*crLHS4 + DN_v(7,1)*crLHS339 + DN_v(7,2)*crLHS536 - crLHS1199*crLHS539 - crLHS1201*crLHS539 + crLHS265);
    rLHS(23,1)+=gauss_weight*(DN_v(7,0)*crLHS31 + DN_v(7,1)*crLHS345 + DN_v(7,2)*crLHS540 - crLHS1199*crLHS542 - crLHS1201*crLHS542 + crLHS488);
    rLHS(23,2)+=gauss_weight*(DN_v(7,0)*crLHS40 + DN_v(7,1)*crLHS353 + DN_v(7,2)*crLHS543 - crLHS1199*crLHS547 - crLHS1201*crLHS547 + crLHS1202 + crLHS643);
    rLHS(23,3)+=gauss_weight*(DN_v(7,0)*crLHS47 + DN_v(7,1)*crLHS358 + DN_v(7,2)*crLHS548 - crLHS1199*crLHS553 - crLHS1201*crLHS553 + crLHS720);
    rLHS(23,4)+=gauss_weight*(DN_v(7,0)*crLHS64 + DN_v(7,1)*crLHS365 + DN_v(7,2)*crLHS554 - crLHS1199*crLHS557 - crLHS1201*crLHS557 + crLHS754);
    rLHS(23,5)+=gauss_weight*(DN_v(7,0)*crLHS72 + DN_v(7,1)*crLHS373 + DN_v(7,2)*crLHS558 - crLHS1199*crLHS562 - crLHS1201*crLHS562 + crLHS1203 + crLHS783);
    rLHS(23,6)+=gauss_weight*(DN_v(7,0)*crLHS79 + DN_v(7,1)*crLHS377 + DN_v(7,2)*crLHS563 - crLHS1199*crLHS567 - crLHS1201*crLHS567 + crLHS827);
    rLHS(23,7)+=gauss_weight*(DN_v(7,0)*crLHS96 + DN_v(7,1)*crLHS384 + DN_v(7,2)*crLHS568 - crLHS1199*crLHS571 - crLHS1201*crLHS571 + crLHS858);
    rLHS(23,8)+=gauss_weight*(DN_v(7,0)*crLHS104 + DN_v(7,1)*crLHS392 + DN_v(7,2)*crLHS572 - crLHS1199*crLHS576 - crLHS1201*crLHS576 + crLHS1204 + crLHS884);
    rLHS(23,9)+=gauss_weight*(DN_v(7,0)*crLHS111 + DN_v(7,1)*crLHS396 + DN_v(7,2)*crLHS577 - crLHS1199*crLHS581 - crLHS1201*crLHS581 + crLHS924);
    rLHS(23,10)+=gauss_weight*(DN_v(7,0)*crLHS128 + DN_v(7,1)*crLHS403 + DN_v(7,2)*crLHS582 - crLHS1199*crLHS585 - crLHS1201*crLHS585 + crLHS952);
    rLHS(23,11)+=gauss_weight*(DN_v(7,0)*crLHS136 + DN_v(7,1)*crLHS411 + DN_v(7,2)*crLHS586 - crLHS1199*crLHS590 - crLHS1201*crLHS590 + crLHS1205 + crLHS975);
    rLHS(23,12)+=gauss_weight*(DN_v(7,0)*crLHS143 + DN_v(7,1)*crLHS415 + DN_v(7,2)*crLHS591 + crLHS1011 - crLHS1199*crLHS595 - crLHS1201*crLHS595);
    rLHS(23,13)+=gauss_weight*(DN_v(7,0)*crLHS160 + DN_v(7,1)*crLHS422 + DN_v(7,2)*crLHS596 + crLHS1036 - crLHS1199*crLHS599 - crLHS1201*crLHS599);
    rLHS(23,14)+=gauss_weight*(DN_v(7,0)*crLHS168 + DN_v(7,1)*crLHS430 + DN_v(7,2)*crLHS600 + crLHS1056 - crLHS1199*crLHS604 - crLHS1201*crLHS604 + crLHS1206);
    rLHS(23,15)+=gauss_weight*(DN_v(7,0)*crLHS175 + DN_v(7,1)*crLHS434 + DN_v(7,2)*crLHS605 + crLHS1088 - crLHS1199*crLHS609 - crLHS1201*crLHS609);
    rLHS(23,16)+=gauss_weight*(DN_v(7,0)*crLHS192 + DN_v(7,1)*crLHS441 + DN_v(7,2)*crLHS610 + crLHS1110 - crLHS1199*crLHS613 - crLHS1201*crLHS613);
    rLHS(23,17)+=gauss_weight*(DN_v(7,0)*crLHS200 + DN_v(7,1)*crLHS449 + DN_v(7,2)*crLHS614 + crLHS1127 - crLHS1199*crLHS618 - crLHS1201*crLHS618 + crLHS1207);
    rLHS(23,18)+=gauss_weight*(DN_v(7,0)*crLHS207 + DN_v(7,1)*crLHS453 + DN_v(7,2)*crLHS619 + crLHS1155 - crLHS1199*crLHS623 - crLHS1201*crLHS623);
    rLHS(23,19)+=gauss_weight*(DN_v(7,0)*crLHS224 + DN_v(7,1)*crLHS460 + DN_v(7,2)*crLHS624 + crLHS1174 - crLHS1199*crLHS627 - crLHS1201*crLHS627);
    rLHS(23,20)+=gauss_weight*(DN_v(7,0)*crLHS232 + DN_v(7,1)*crLHS468 + DN_v(7,2)*crLHS628 + crLHS1188 - crLHS1199*crLHS632 - crLHS1201*crLHS632 + crLHS1208);
    rLHS(23,21)+=gauss_weight*(DN_v(7,0)*crLHS239 + DN_v(7,1)*crLHS472 + DN_v(7,2)*crLHS633 - crLHS1199*crLHS637 - crLHS1201*crLHS637 + crLHS1212);
    rLHS(23,22)+=gauss_weight*(DN_v(7,0)*crLHS256 + DN_v(7,1)*crLHS479 + DN_v(7,2)*crLHS638 - crLHS1199*crLHS641 - crLHS1201*crLHS641 + crLHS1228);
    rLHS(23,23)+=gauss_weight*(DN_v(7,0)*crLHS264 + DN_v(7,1)*crLHS487 + std::pow(DN_v(7,2), 2)*crLHS9 + DN_v(7,2)*crLHS642 - crLHS1199*crLHS646 - crLHS1201*crLHS646 + crLHS1209);
    rLHS(23,24)+=gauss_weight*(DN_v(7,0)*crLHS271 + DN_v(7,1)*crLHS491 + DN_v(7,2)*crLHS647 - crLHS1199*crLHS651 - crLHS1201*crLHS651 + crLHS1240);
    rLHS(23,25)+=gauss_weight*(DN_v(7,0)*crLHS288 + DN_v(7,1)*crLHS498 + DN_v(7,2)*crLHS652 - crLHS1199*crLHS655 - crLHS1201*crLHS655 + crLHS1241);
    rLHS(23,26)+=gauss_weight*(DN_v(7,0)*crLHS296 + DN_v(7,1)*crLHS506 + DN_v(7,2)*crLHS656 - crLHS1199*crLHS660 - crLHS1201*crLHS660 + crLHS1215 + crLHS1242);
    rLHS(23,27)+=gauss_weight*(DN_v(7,0)*crLHS303 + DN_v(7,1)*crLHS510 + DN_v(7,2)*crLHS661 - crLHS1199*crLHS665 - crLHS1201*crLHS665 + crLHS1243);
    rLHS(23,28)+=gauss_weight*(DN_v(7,0)*crLHS319 + DN_v(7,1)*crLHS517 + DN_v(7,2)*crLHS666 - crLHS1199*crLHS669 - crLHS1201*crLHS669 + crLHS1244);
    rLHS(23,29)+=gauss_weight*(DN_v(7,0)*crLHS327 + DN_v(7,1)*crLHS525 + DN_v(7,2)*crLHS670 - crLHS1199*crLHS674 - crLHS1201*crLHS674 + crLHS1220 + crLHS1245);
    rLHS(23,30)+=gauss_weight*(crLHS1200*crLHS676 + crLHS1246 + crLHS245*crLHS676);
    rLHS(23,31)+=gauss_weight*(crLHS1200*crLHS678 + crLHS1247 + crLHS245*crLHS678);
    rLHS(23,32)+=gauss_weight*(crLHS1200*crLHS680 + crLHS1248 + crLHS245*crLHS680);
    rLHS(23,33)+=gauss_weight*(crLHS1200*crLHS682 + crLHS1249 + crLHS245*crLHS682);
    rLHS(24,0)+=gauss_weight*(DN_v(8,0)*crLHS0 + DN_v(8,1)*crLHS2 + DN_v(8,2)*crLHS4 - crLHS1250*crLHS19 - crLHS1252*crLHS19 + crLHS1253 + crLHS272);
    rLHS(24,1)+=gauss_weight*(DN_v(8,0)*crLHS26 + DN_v(8,1)*crLHS28 + DN_v(8,2)*crLHS31 - crLHS1250*crLHS35 - crLHS1252*crLHS35 + crLHS492);
    rLHS(24,2)+=gauss_weight*(DN_v(8,0)*crLHS36 + DN_v(8,1)*crLHS38 + DN_v(8,2)*crLHS40 - crLHS1250*crLHS42 - crLHS1252*crLHS42 + crLHS648);
    rLHS(24,3)+=gauss_weight*(DN_v(8,0)*crLHS43 + DN_v(8,1)*crLHS45 + DN_v(8,2)*crLHS47 - crLHS1250*crLHS56 - crLHS1252*crLHS56 + crLHS1254 + crLHS721);
    rLHS(24,4)+=gauss_weight*(DN_v(8,0)*crLHS59 + DN_v(8,1)*crLHS61 + DN_v(8,2)*crLHS64 - crLHS1250*crLHS67 - crLHS1252*crLHS67 + crLHS755);
    rLHS(24,5)+=gauss_weight*(DN_v(8,0)*crLHS68 + DN_v(8,1)*crLHS70 + DN_v(8,2)*crLHS72 - crLHS1250*crLHS74 - crLHS1252*crLHS74 + crLHS784);
    rLHS(24,6)+=gauss_weight*(DN_v(8,0)*crLHS75 + DN_v(8,1)*crLHS77 + DN_v(8,2)*crLHS79 - crLHS1250*crLHS88 - crLHS1252*crLHS88 + crLHS1255 + crLHS828);
    rLHS(24,7)+=gauss_weight*(DN_v(8,0)*crLHS91 + DN_v(8,1)*crLHS93 + DN_v(8,2)*crLHS96 - crLHS1250*crLHS99 - crLHS1252*crLHS99 + crLHS859);
    rLHS(24,8)+=gauss_weight*(DN_v(8,0)*crLHS100 + DN_v(8,1)*crLHS102 + DN_v(8,2)*crLHS104 - crLHS106*crLHS1250 - crLHS106*crLHS1252 + crLHS885);
    rLHS(24,9)+=gauss_weight*(DN_v(8,0)*crLHS107 + DN_v(8,1)*crLHS109 + DN_v(8,2)*crLHS111 - crLHS120*crLHS1250 - crLHS120*crLHS1252 + crLHS1256 + crLHS925);
    rLHS(24,10)+=gauss_weight*(DN_v(8,0)*crLHS123 + DN_v(8,1)*crLHS125 + DN_v(8,2)*crLHS128 - crLHS1250*crLHS131 - crLHS1252*crLHS131 + crLHS953);
    rLHS(24,11)+=gauss_weight*(DN_v(8,0)*crLHS132 + DN_v(8,1)*crLHS134 + DN_v(8,2)*crLHS136 - crLHS1250*crLHS138 - crLHS1252*crLHS138 + crLHS976);
    rLHS(24,12)+=gauss_weight*(DN_v(8,0)*crLHS139 + DN_v(8,1)*crLHS141 + DN_v(8,2)*crLHS143 + crLHS1012 - crLHS1250*crLHS152 - crLHS1252*crLHS152 + crLHS1257);
    rLHS(24,13)+=gauss_weight*(DN_v(8,0)*crLHS155 + DN_v(8,1)*crLHS157 + DN_v(8,2)*crLHS160 + crLHS1037 - crLHS1250*crLHS163 - crLHS1252*crLHS163);
    rLHS(24,14)+=gauss_weight*(DN_v(8,0)*crLHS164 + DN_v(8,1)*crLHS166 + DN_v(8,2)*crLHS168 + crLHS1057 - crLHS1250*crLHS170 - crLHS1252*crLHS170);
    rLHS(24,15)+=gauss_weight*(DN_v(8,0)*crLHS171 + DN_v(8,1)*crLHS173 + DN_v(8,2)*crLHS175 + crLHS1089 - crLHS1250*crLHS184 - crLHS1252*crLHS184 + crLHS1258);
    rLHS(24,16)+=gauss_weight*(DN_v(8,0)*crLHS187 + DN_v(8,1)*crLHS189 + DN_v(8,2)*crLHS192 + crLHS1111 - crLHS1250*crLHS195 - crLHS1252*crLHS195);
    rLHS(24,17)+=gauss_weight*(DN_v(8,0)*crLHS196 + DN_v(8,1)*crLHS198 + DN_v(8,2)*crLHS200 + crLHS1128 - crLHS1250*crLHS202 - crLHS1252*crLHS202);
    rLHS(24,18)+=gauss_weight*(DN_v(8,0)*crLHS203 + DN_v(8,1)*crLHS205 + DN_v(8,2)*crLHS207 + crLHS1156 - crLHS1250*crLHS216 - crLHS1252*crLHS216 + crLHS1259);
    rLHS(24,19)+=gauss_weight*(DN_v(8,0)*crLHS219 + DN_v(8,1)*crLHS221 + DN_v(8,2)*crLHS224 + crLHS1175 - crLHS1250*crLHS227 - crLHS1252*crLHS227);
    rLHS(24,20)+=gauss_weight*(DN_v(8,0)*crLHS228 + DN_v(8,1)*crLHS230 + DN_v(8,2)*crLHS232 + crLHS1189 - crLHS1250*crLHS234 - crLHS1252*crLHS234);
    rLHS(24,21)+=gauss_weight*(DN_v(8,0)*crLHS235 + DN_v(8,1)*crLHS237 + DN_v(8,2)*crLHS239 + crLHS1213 - crLHS1250*crLHS248 - crLHS1252*crLHS248 + crLHS1260);
    rLHS(24,22)+=gauss_weight*(DN_v(8,0)*crLHS251 + DN_v(8,1)*crLHS253 + DN_v(8,2)*crLHS256 + crLHS1229 - crLHS1250*crLHS259 - crLHS1252*crLHS259);
    rLHS(24,23)+=gauss_weight*(DN_v(8,0)*crLHS260 + DN_v(8,1)*crLHS262 + DN_v(8,2)*crLHS264 + crLHS1240 - crLHS1250*crLHS266 - crLHS1252*crLHS266);
    rLHS(24,24)+=gauss_weight*(std::pow(DN_v(8,0), 2)*crLHS9 + DN_v(8,0)*crLHS267 + DN_v(8,1)*crLHS269 + DN_v(8,2)*crLHS271 - crLHS1250*crLHS280 - crLHS1252*crLHS280 + crLHS1261);
    rLHS(24,25)+=gauss_weight*(DN_v(8,0)*crLHS283 + DN_v(8,1)*crLHS285 + DN_v(8,2)*crLHS288 - crLHS1250*crLHS291 - crLHS1252*crLHS291 + crLHS1263);
    rLHS(24,26)+=gauss_weight*(DN_v(8,0)*crLHS292 + DN_v(8,1)*crLHS294 + DN_v(8,2)*crLHS296 - crLHS1250*crLHS298 - crLHS1252*crLHS298 + crLHS1264);
    rLHS(24,27)+=gauss_weight*(DN_v(8,0)*crLHS299 + DN_v(8,1)*crLHS301 + DN_v(8,2)*crLHS303 - crLHS1250*crLHS311 - crLHS1252*crLHS311 + crLHS1265 + crLHS1267);
    rLHS(24,28)+=gauss_weight*(DN_v(8,0)*crLHS314 + DN_v(8,1)*crLHS316 + DN_v(8,2)*crLHS319 - crLHS1250*crLHS322 - crLHS1252*crLHS322 + crLHS1268);
    rLHS(24,29)+=gauss_weight*(DN_v(8,0)*crLHS323 + DN_v(8,1)*crLHS325 + DN_v(8,2)*crLHS327 - crLHS1250*crLHS329 - crLHS1252*crLHS329 + crLHS1269);
    rLHS(24,30)+=gauss_weight*(crLHS1251*crLHS331 + crLHS1270 + crLHS277*crLHS331);
    rLHS(24,31)+=gauss_weight*(crLHS1251*crLHS333 + crLHS1271 + crLHS277*crLHS333);
    rLHS(24,32)+=gauss_weight*(crLHS1251*crLHS335 + crLHS1272 + crLHS277*crLHS335);
    rLHS(24,33)+=gauss_weight*(crLHS1251*crLHS337 + crLHS1273 + crLHS277*crLHS337);
    rLHS(25,0)+=gauss_weight*(DN_v(8,0)*crLHS2 + DN_v(8,1)*crLHS338 + DN_v(8,2)*crLHS339 - crLHS1250*crLHS342 - crLHS1252*crLHS342 + crLHS289);
    rLHS(25,1)+=gauss_weight*(DN_v(8,0)*crLHS28 + DN_v(8,1)*crLHS343 + DN_v(8,2)*crLHS345 - crLHS1250*crLHS350 - crLHS1252*crLHS350 + crLHS1253 + crLHS499);
    rLHS(25,2)+=gauss_weight*(DN_v(8,0)*crLHS38 + DN_v(8,1)*crLHS351 + DN_v(8,2)*crLHS353 - crLHS1250*crLHS356 - crLHS1252*crLHS356 + crLHS653);
    rLHS(25,3)+=gauss_weight*(DN_v(8,0)*crLHS45 + DN_v(8,1)*crLHS357 + DN_v(8,2)*crLHS358 - crLHS1250*crLHS362 - crLHS1252*crLHS362 + crLHS724);
    rLHS(25,4)+=gauss_weight*(DN_v(8,0)*crLHS61 + DN_v(8,1)*crLHS363 + DN_v(8,2)*crLHS365 - crLHS1250*crLHS370 - crLHS1252*crLHS370 + crLHS1254 + crLHS756);
    rLHS(25,5)+=gauss_weight*(DN_v(8,0)*crLHS70 + DN_v(8,1)*crLHS371 + DN_v(8,2)*crLHS373 - crLHS1250*crLHS375 - crLHS1252*crLHS375 + crLHS785);
    rLHS(25,6)+=gauss_weight*(DN_v(8,0)*crLHS77 + DN_v(8,1)*crLHS376 + DN_v(8,2)*crLHS377 - crLHS1250*crLHS381 - crLHS1252*crLHS381 + crLHS831);
    rLHS(25,7)+=gauss_weight*(DN_v(8,0)*crLHS93 + DN_v(8,1)*crLHS382 + DN_v(8,2)*crLHS384 - crLHS1250*crLHS389 - crLHS1252*crLHS389 + crLHS1255 + crLHS860);
    rLHS(25,8)+=gauss_weight*(DN_v(8,0)*crLHS102 + DN_v(8,1)*crLHS390 + DN_v(8,2)*crLHS392 - crLHS1250*crLHS394 - crLHS1252*crLHS394 + crLHS886);
    rLHS(25,9)+=gauss_weight*(DN_v(8,0)*crLHS109 + DN_v(8,1)*crLHS395 + DN_v(8,2)*crLHS396 - crLHS1250*crLHS400 - crLHS1252*crLHS400 + crLHS928);
    rLHS(25,10)+=gauss_weight*(DN_v(8,0)*crLHS125 + DN_v(8,1)*crLHS401 + DN_v(8,2)*crLHS403 - crLHS1250*crLHS408 - crLHS1252*crLHS408 + crLHS1256 + crLHS954);
    rLHS(25,11)+=gauss_weight*(DN_v(8,0)*crLHS134 + DN_v(8,1)*crLHS409 + DN_v(8,2)*crLHS411 - crLHS1250*crLHS413 - crLHS1252*crLHS413 + crLHS977);
    rLHS(25,12)+=gauss_weight*(DN_v(8,0)*crLHS141 + DN_v(8,1)*crLHS414 + DN_v(8,2)*crLHS415 + crLHS1015 - crLHS1250*crLHS419 - crLHS1252*crLHS419);
    rLHS(25,13)+=gauss_weight*(DN_v(8,0)*crLHS157 + DN_v(8,1)*crLHS420 + DN_v(8,2)*crLHS422 + crLHS1038 - crLHS1250*crLHS427 - crLHS1252*crLHS427 + crLHS1257);
    rLHS(25,14)+=gauss_weight*(DN_v(8,0)*crLHS166 + DN_v(8,1)*crLHS428 + DN_v(8,2)*crLHS430 + crLHS1058 - crLHS1250*crLHS432 - crLHS1252*crLHS432);
    rLHS(25,15)+=gauss_weight*(DN_v(8,0)*crLHS173 + DN_v(8,1)*crLHS433 + DN_v(8,2)*crLHS434 + crLHS1092 - crLHS1250*crLHS438 - crLHS1252*crLHS438);
    rLHS(25,16)+=gauss_weight*(DN_v(8,0)*crLHS189 + DN_v(8,1)*crLHS439 + DN_v(8,2)*crLHS441 + crLHS1112 - crLHS1250*crLHS446 - crLHS1252*crLHS446 + crLHS1258);
    rLHS(25,17)+=gauss_weight*(DN_v(8,0)*crLHS198 + DN_v(8,1)*crLHS447 + DN_v(8,2)*crLHS449 + crLHS1129 - crLHS1250*crLHS451 - crLHS1252*crLHS451);
    rLHS(25,18)+=gauss_weight*(DN_v(8,0)*crLHS205 + DN_v(8,1)*crLHS452 + DN_v(8,2)*crLHS453 + crLHS1159 - crLHS1250*crLHS457 - crLHS1252*crLHS457);
    rLHS(25,19)+=gauss_weight*(DN_v(8,0)*crLHS221 + DN_v(8,1)*crLHS458 + DN_v(8,2)*crLHS460 + crLHS1176 - crLHS1250*crLHS465 - crLHS1252*crLHS465 + crLHS1259);
    rLHS(25,20)+=gauss_weight*(DN_v(8,0)*crLHS230 + DN_v(8,1)*crLHS466 + DN_v(8,2)*crLHS468 + crLHS1190 - crLHS1250*crLHS470 - crLHS1252*crLHS470);
    rLHS(25,21)+=gauss_weight*(DN_v(8,0)*crLHS237 + DN_v(8,1)*crLHS471 + DN_v(8,2)*crLHS472 + crLHS1216 - crLHS1250*crLHS476 - crLHS1252*crLHS476);
    rLHS(25,22)+=gauss_weight*(DN_v(8,0)*crLHS253 + DN_v(8,1)*crLHS477 + DN_v(8,2)*crLHS479 + crLHS1230 - crLHS1250*crLHS484 - crLHS1252*crLHS484 + crLHS1260);
    rLHS(25,23)+=gauss_weight*(DN_v(8,0)*crLHS262 + DN_v(8,1)*crLHS485 + DN_v(8,2)*crLHS487 + crLHS1241 - crLHS1250*crLHS489 - crLHS1252*crLHS489);
    rLHS(25,24)+=gauss_weight*(DN_v(8,0)*crLHS269 + DN_v(8,1)*crLHS490 + DN_v(8,2)*crLHS491 - crLHS1250*crLHS495 - crLHS1252*crLHS495 + crLHS1263);
    rLHS(25,25)+=gauss_weight*(DN_v(8,0)*crLHS285 + std::pow(DN_v(8,1), 2)*crLHS9 + DN_v(8,1)*crLHS496 + DN_v(8,2)*crLHS498 - crLHS1250*crLHS503 - crLHS1252*crLHS503 + crLHS1261);
    rLHS(25,26)+=gauss_weight*(DN_v(8,0)*crLHS294 + DN_v(8,1)*crLHS504 + DN_v(8,2)*crLHS506 - crLHS1250*crLHS508 - crLHS1252*crLHS508 + crLHS1275);
    rLHS(25,27)+=gauss_weight*(DN_v(8,0)*crLHS301 + DN_v(8,1)*crLHS509 + DN_v(8,2)*crLHS510 - crLHS1250*crLHS514 - crLHS1252*crLHS514 + crLHS1276);
    rLHS(25,28)+=gauss_weight*(DN_v(8,0)*crLHS316 + DN_v(8,1)*crLHS515 + DN_v(8,2)*crLHS517 - crLHS1250*crLHS522 - crLHS1252*crLHS522 + crLHS1267 + crLHS1277);
    rLHS(25,29)+=gauss_weight*(DN_v(8,0)*crLHS325 + DN_v(8,1)*crLHS523 + DN_v(8,2)*crLHS525 - crLHS1250*crLHS527 - crLHS1252*crLHS527 + crLHS1278);
    rLHS(25,30)+=gauss_weight*(crLHS1251*crLHS529 + crLHS1279 + crLHS277*crLHS529);
    rLHS(25,31)+=gauss_weight*(crLHS1251*crLHS531 + crLHS1280 + crLHS277*crLHS531);
    rLHS(25,32)+=gauss_weight*(crLHS1251*crLHS533 + crLHS1281 + crLHS277*crLHS533);
    rLHS(25,33)+=gauss_weight*(crLHS1251*crLHS535 + crLHS1282 + crLHS277*crLHS535);
    rLHS(26,0)+=gauss_weight*(DN_v(8,0)*crLHS4 + DN_v(8,1)*crLHS339 + DN_v(8,2)*crLHS536 - crLHS1250*crLHS539 - crLHS1252*crLHS539 + crLHS297);
    rLHS(26,1)+=gauss_weight*(DN_v(8,0)*crLHS31 + DN_v(8,1)*crLHS345 + DN_v(8,2)*crLHS540 - crLHS1250*crLHS542 - crLHS1252*crLHS542 + crLHS507);
    rLHS(26,2)+=gauss_weight*(DN_v(8,0)*crLHS40 + DN_v(8,1)*crLHS353 + DN_v(8,2)*crLHS543 - crLHS1250*crLHS547 - crLHS1252*crLHS547 + crLHS1253 + crLHS657);
    rLHS(26,3)+=gauss_weight*(DN_v(8,0)*crLHS47 + DN_v(8,1)*crLHS358 + DN_v(8,2)*crLHS548 - crLHS1250*crLHS553 - crLHS1252*crLHS553 + crLHS725);
    rLHS(26,4)+=gauss_weight*(DN_v(8,0)*crLHS64 + DN_v(8,1)*crLHS365 + DN_v(8,2)*crLHS554 - crLHS1250*crLHS557 - crLHS1252*crLHS557 + crLHS757);
    rLHS(26,5)+=gauss_weight*(DN_v(8,0)*crLHS72 + DN_v(8,1)*crLHS373 + DN_v(8,2)*crLHS558 - crLHS1250*crLHS562 - crLHS1252*crLHS562 + crLHS1254 + crLHS786);
    rLHS(26,6)+=gauss_weight*(DN_v(8,0)*crLHS79 + DN_v(8,1)*crLHS377 + DN_v(8,2)*crLHS563 - crLHS1250*crLHS567 - crLHS1252*crLHS567 + crLHS832);
    rLHS(26,7)+=gauss_weight*(DN_v(8,0)*crLHS96 + DN_v(8,1)*crLHS384 + DN_v(8,2)*crLHS568 - crLHS1250*crLHS571 - crLHS1252*crLHS571 + crLHS861);
    rLHS(26,8)+=gauss_weight*(DN_v(8,0)*crLHS104 + DN_v(8,1)*crLHS392 + DN_v(8,2)*crLHS572 - crLHS1250*crLHS576 - crLHS1252*crLHS576 + crLHS1255 + crLHS887);
    rLHS(26,9)+=gauss_weight*(DN_v(8,0)*crLHS111 + DN_v(8,1)*crLHS396 + DN_v(8,2)*crLHS577 - crLHS1250*crLHS581 - crLHS1252*crLHS581 + crLHS929);
    rLHS(26,10)+=gauss_weight*(DN_v(8,0)*crLHS128 + DN_v(8,1)*crLHS403 + DN_v(8,2)*crLHS582 - crLHS1250*crLHS585 - crLHS1252*crLHS585 + crLHS955);
    rLHS(26,11)+=gauss_weight*(DN_v(8,0)*crLHS136 + DN_v(8,1)*crLHS411 + DN_v(8,2)*crLHS586 - crLHS1250*crLHS590 - crLHS1252*crLHS590 + crLHS1256 + crLHS978);
    rLHS(26,12)+=gauss_weight*(DN_v(8,0)*crLHS143 + DN_v(8,1)*crLHS415 + DN_v(8,2)*crLHS591 + crLHS1016 - crLHS1250*crLHS595 - crLHS1252*crLHS595);
    rLHS(26,13)+=gauss_weight*(DN_v(8,0)*crLHS160 + DN_v(8,1)*crLHS422 + DN_v(8,2)*crLHS596 + crLHS1039 - crLHS1250*crLHS599 - crLHS1252*crLHS599);
    rLHS(26,14)+=gauss_weight*(DN_v(8,0)*crLHS168 + DN_v(8,1)*crLHS430 + DN_v(8,2)*crLHS600 + crLHS1059 - crLHS1250*crLHS604 - crLHS1252*crLHS604 + crLHS1257);
    rLHS(26,15)+=gauss_weight*(DN_v(8,0)*crLHS175 + DN_v(8,1)*crLHS434 + DN_v(8,2)*crLHS605 + crLHS1093 - crLHS1250*crLHS609 - crLHS1252*crLHS609);
    rLHS(26,16)+=gauss_weight*(DN_v(8,0)*crLHS192 + DN_v(8,1)*crLHS441 + DN_v(8,2)*crLHS610 + crLHS1113 - crLHS1250*crLHS613 - crLHS1252*crLHS613);
    rLHS(26,17)+=gauss_weight*(DN_v(8,0)*crLHS200 + DN_v(8,1)*crLHS449 + DN_v(8,2)*crLHS614 + crLHS1130 - crLHS1250*crLHS618 - crLHS1252*crLHS618 + crLHS1258);
    rLHS(26,18)+=gauss_weight*(DN_v(8,0)*crLHS207 + DN_v(8,1)*crLHS453 + DN_v(8,2)*crLHS619 + crLHS1160 - crLHS1250*crLHS623 - crLHS1252*crLHS623);
    rLHS(26,19)+=gauss_weight*(DN_v(8,0)*crLHS224 + DN_v(8,1)*crLHS460 + DN_v(8,2)*crLHS624 + crLHS1177 - crLHS1250*crLHS627 - crLHS1252*crLHS627);
    rLHS(26,20)+=gauss_weight*(DN_v(8,0)*crLHS232 + DN_v(8,1)*crLHS468 + DN_v(8,2)*crLHS628 + crLHS1191 - crLHS1250*crLHS632 - crLHS1252*crLHS632 + crLHS1259);
    rLHS(26,21)+=gauss_weight*(DN_v(8,0)*crLHS239 + DN_v(8,1)*crLHS472 + DN_v(8,2)*crLHS633 + crLHS1217 - crLHS1250*crLHS637 - crLHS1252*crLHS637);
    rLHS(26,22)+=gauss_weight*(DN_v(8,0)*crLHS256 + DN_v(8,1)*crLHS479 + DN_v(8,2)*crLHS638 + crLHS1231 - crLHS1250*crLHS641 - crLHS1252*crLHS641);
    rLHS(26,23)+=gauss_weight*(DN_v(8,0)*crLHS264 + DN_v(8,1)*crLHS487 + DN_v(8,2)*crLHS642 + crLHS1242 - crLHS1250*crLHS646 - crLHS1252*crLHS646 + crLHS1260);
    rLHS(26,24)+=gauss_weight*(DN_v(8,0)*crLHS271 + DN_v(8,1)*crLHS491 + DN_v(8,2)*crLHS647 - crLHS1250*crLHS651 - crLHS1252*crLHS651 + crLHS1264);
    rLHS(26,25)+=gauss_weight*(DN_v(8,0)*crLHS288 + DN_v(8,1)*crLHS498 + DN_v(8,2)*crLHS652 - crLHS1250*crLHS655 - crLHS1252*crLHS655 + crLHS1275);
    rLHS(26,26)+=gauss_weight*(DN_v(8,0)*crLHS296 + DN_v(8,1)*crLHS506 + std::pow(DN_v(8,2), 2)*crLHS9 + DN_v(8,2)*crLHS656 - crLHS1250*crLHS660 - crLHS1252*crLHS660 + crLHS1261);
    rLHS(26,27)+=gauss_weight*(DN_v(8,0)*crLHS303 + DN_v(8,1)*crLHS510 + DN_v(8,2)*crLHS661 - crLHS1250*crLHS665 - crLHS1252*crLHS665 + crLHS1284);
    rLHS(26,28)+=gauss_weight*(DN_v(8,0)*crLHS319 + DN_v(8,1)*crLHS517 + DN_v(8,2)*crLHS666 - crLHS1250*crLHS669 - crLHS1252*crLHS669 + crLHS1285);
    rLHS(26,29)+=gauss_weight*(DN_v(8,0)*crLHS327 + DN_v(8,1)*crLHS525 + DN_v(8,2)*crLHS670 - crLHS1250*crLHS674 - crLHS1252*crLHS674 + crLHS1267 + crLHS1286);
    rLHS(26,30)+=gauss_weight*(crLHS1251*crLHS676 + crLHS1287 + crLHS277*crLHS676);
    rLHS(26,31)+=gauss_weight*(crLHS1251*crLHS678 + crLHS1288 + crLHS277*crLHS678);
    rLHS(26,32)+=gauss_weight*(crLHS1251*crLHS680 + crLHS1289 + crLHS277*crLHS680);
    rLHS(26,33)+=gauss_weight*(crLHS1251*crLHS682 + crLHS1290 + crLHS277*crLHS682);
    rLHS(27,0)+=gauss_weight*(DN_v(9,0)*crLHS0 + DN_v(9,1)*crLHS2 + DN_v(9,2)*crLHS4 - crLHS1291*crLHS19 - crLHS1293*crLHS19 + crLHS1294 + crLHS304);
    rLHS(27,1)+=gauss_weight*(DN_v(9,0)*crLHS26 + DN_v(9,1)*crLHS28 + DN_v(9,2)*crLHS31 - crLHS1291*crLHS35 - crLHS1293*crLHS35 + crLHS511);
    rLHS(27,2)+=gauss_weight*(DN_v(9,0)*crLHS36 + DN_v(9,1)*crLHS38 + DN_v(9,2)*crLHS40 - crLHS1291*crLHS42 - crLHS1293*crLHS42 + crLHS662);
    rLHS(27,3)+=gauss_weight*(DN_v(9,0)*crLHS43 + DN_v(9,1)*crLHS45 + DN_v(9,2)*crLHS47 - crLHS1291*crLHS56 - crLHS1293*crLHS56 + crLHS1295 + crLHS726);
    rLHS(27,4)+=gauss_weight*(DN_v(9,0)*crLHS59 + DN_v(9,1)*crLHS61 + DN_v(9,2)*crLHS64 - crLHS1291*crLHS67 - crLHS1293*crLHS67 + crLHS758);
    rLHS(27,5)+=gauss_weight*(DN_v(9,0)*crLHS68 + DN_v(9,1)*crLHS70 + DN_v(9,2)*crLHS72 - crLHS1291*crLHS74 - crLHS1293*crLHS74 + crLHS787);
    rLHS(27,6)+=gauss_weight*(DN_v(9,0)*crLHS75 + DN_v(9,1)*crLHS77 + DN_v(9,2)*crLHS79 - crLHS1291*crLHS88 - crLHS1293*crLHS88 + crLHS1296 + crLHS833);
    rLHS(27,7)+=gauss_weight*(DN_v(9,0)*crLHS91 + DN_v(9,1)*crLHS93 + DN_v(9,2)*crLHS96 - crLHS1291*crLHS99 - crLHS1293*crLHS99 + crLHS862);
    rLHS(27,8)+=gauss_weight*(DN_v(9,0)*crLHS100 + DN_v(9,1)*crLHS102 + DN_v(9,2)*crLHS104 - crLHS106*crLHS1291 - crLHS106*crLHS1293 + crLHS888);
    rLHS(27,9)+=gauss_weight*(DN_v(9,0)*crLHS107 + DN_v(9,1)*crLHS109 + DN_v(9,2)*crLHS111 - crLHS120*crLHS1291 - crLHS120*crLHS1293 + crLHS1297 + crLHS930);
    rLHS(27,10)+=gauss_weight*(DN_v(9,0)*crLHS123 + DN_v(9,1)*crLHS125 + DN_v(9,2)*crLHS128 - crLHS1291*crLHS131 - crLHS1293*crLHS131 + crLHS956);
    rLHS(27,11)+=gauss_weight*(DN_v(9,0)*crLHS132 + DN_v(9,1)*crLHS134 + DN_v(9,2)*crLHS136 - crLHS1291*crLHS138 - crLHS1293*crLHS138 + crLHS979);
    rLHS(27,12)+=gauss_weight*(DN_v(9,0)*crLHS139 + DN_v(9,1)*crLHS141 + DN_v(9,2)*crLHS143 + crLHS1017 - crLHS1291*crLHS152 - crLHS1293*crLHS152 + crLHS1298);
    rLHS(27,13)+=gauss_weight*(DN_v(9,0)*crLHS155 + DN_v(9,1)*crLHS157 + DN_v(9,2)*crLHS160 + crLHS1040 - crLHS1291*crLHS163 - crLHS1293*crLHS163);
    rLHS(27,14)+=gauss_weight*(DN_v(9,0)*crLHS164 + DN_v(9,1)*crLHS166 + DN_v(9,2)*crLHS168 + crLHS1060 - crLHS1291*crLHS170 - crLHS1293*crLHS170);
    rLHS(27,15)+=gauss_weight*(DN_v(9,0)*crLHS171 + DN_v(9,1)*crLHS173 + DN_v(9,2)*crLHS175 + crLHS1094 - crLHS1291*crLHS184 - crLHS1293*crLHS184 + crLHS1299);
    rLHS(27,16)+=gauss_weight*(DN_v(9,0)*crLHS187 + DN_v(9,1)*crLHS189 + DN_v(9,2)*crLHS192 + crLHS1114 - crLHS1291*crLHS195 - crLHS1293*crLHS195);
    rLHS(27,17)+=gauss_weight*(DN_v(9,0)*crLHS196 + DN_v(9,1)*crLHS198 + DN_v(9,2)*crLHS200 + crLHS1131 - crLHS1291*crLHS202 - crLHS1293*crLHS202);
    rLHS(27,18)+=gauss_weight*(DN_v(9,0)*crLHS203 + DN_v(9,1)*crLHS205 + DN_v(9,2)*crLHS207 + crLHS1161 - crLHS1291*crLHS216 - crLHS1293*crLHS216 + crLHS1300);
    rLHS(27,19)+=gauss_weight*(DN_v(9,0)*crLHS219 + DN_v(9,1)*crLHS221 + DN_v(9,2)*crLHS224 + crLHS1178 - crLHS1291*crLHS227 - crLHS1293*crLHS227);
    rLHS(27,20)+=gauss_weight*(DN_v(9,0)*crLHS228 + DN_v(9,1)*crLHS230 + DN_v(9,2)*crLHS232 + crLHS1192 - crLHS1291*crLHS234 - crLHS1293*crLHS234);
    rLHS(27,21)+=gauss_weight*(DN_v(9,0)*crLHS235 + DN_v(9,1)*crLHS237 + DN_v(9,2)*crLHS239 + crLHS1218 - crLHS1291*crLHS248 - crLHS1293*crLHS248 + crLHS1301);
    rLHS(27,22)+=gauss_weight*(DN_v(9,0)*crLHS251 + DN_v(9,1)*crLHS253 + DN_v(9,2)*crLHS256 + crLHS1232 - crLHS1291*crLHS259 - crLHS1293*crLHS259);
    rLHS(27,23)+=gauss_weight*(DN_v(9,0)*crLHS260 + DN_v(9,1)*crLHS262 + DN_v(9,2)*crLHS264 + crLHS1243 - crLHS1291*crLHS266 - crLHS1293*crLHS266);
    rLHS(27,24)+=gauss_weight*(DN_v(9,0)*crLHS267 + DN_v(9,1)*crLHS269 + DN_v(9,2)*crLHS271 + crLHS1265 - crLHS1291*crLHS280 - crLHS1293*crLHS280 + crLHS1302);
    rLHS(27,25)+=gauss_weight*(DN_v(9,0)*crLHS283 + DN_v(9,1)*crLHS285 + DN_v(9,2)*crLHS288 + crLHS1276 - crLHS1291*crLHS291 - crLHS1293*crLHS291);
    rLHS(27,26)+=gauss_weight*(DN_v(9,0)*crLHS292 + DN_v(9,1)*crLHS294 + DN_v(9,2)*crLHS296 + crLHS1284 - crLHS1291*crLHS298 - crLHS1293*crLHS298);
    rLHS(27,27)+=gauss_weight*(std::pow(DN_v(9,0), 2)*crLHS9 + DN_v(9,0)*crLHS299 + DN_v(9,1)*crLHS301 + DN_v(9,2)*crLHS303 - crLHS1291*crLHS311 - crLHS1293*crLHS311 + crLHS1303);
    rLHS(27,28)+=gauss_weight*(DN_v(9,0)*crLHS314 + DN_v(9,1)*crLHS316 + DN_v(9,2)*crLHS319 - crLHS1291*crLHS322 - crLHS1293*crLHS322 + crLHS1305);
    rLHS(27,29)+=gauss_weight*(DN_v(9,0)*crLHS323 + DN_v(9,1)*crLHS325 + DN_v(9,2)*crLHS327 - crLHS1291*crLHS329 - crLHS1293*crLHS329 + crLHS1306);
    rLHS(27,30)+=gauss_weight*(crLHS1292*crLHS331 + crLHS1307 + crLHS309*crLHS331);
    rLHS(27,31)+=gauss_weight*(crLHS1292*crLHS333 + crLHS1308 + crLHS309*crLHS333);
    rLHS(27,32)+=gauss_weight*(crLHS1292*crLHS335 + crLHS1309 + crLHS309*crLHS335);
    rLHS(27,33)+=gauss_weight*(crLHS1292*crLHS337 + crLHS1310 + crLHS309*crLHS337);
    rLHS(28,0)+=gauss_weight*(DN_v(9,0)*crLHS2 + DN_v(9,1)*crLHS338 + DN_v(9,2)*crLHS339 - crLHS1291*crLHS342 - crLHS1293*crLHS342 + crLHS320);
    rLHS(28,1)+=gauss_weight*(DN_v(9,0)*crLHS28 + DN_v(9,1)*crLHS343 + DN_v(9,2)*crLHS345 - crLHS1291*crLHS350 - crLHS1293*crLHS350 + crLHS1294 + crLHS518);
    rLHS(28,2)+=gauss_weight*(DN_v(9,0)*crLHS38 + DN_v(9,1)*crLHS351 + DN_v(9,2)*crLHS353 - crLHS1291*crLHS356 - crLHS1293*crLHS356 + crLHS667);
    rLHS(28,3)+=gauss_weight*(DN_v(9,0)*crLHS45 + DN_v(9,1)*crLHS357 + DN_v(9,2)*crLHS358 - crLHS1291*crLHS362 - crLHS1293*crLHS362 + crLHS729);
    rLHS(28,4)+=gauss_weight*(DN_v(9,0)*crLHS61 + DN_v(9,1)*crLHS363 + DN_v(9,2)*crLHS365 - crLHS1291*crLHS370 - crLHS1293*crLHS370 + crLHS1295 + crLHS759);
    rLHS(28,5)+=gauss_weight*(DN_v(9,0)*crLHS70 + DN_v(9,1)*crLHS371 + DN_v(9,2)*crLHS373 - crLHS1291*crLHS375 - crLHS1293*crLHS375 + crLHS788);
    rLHS(28,6)+=gauss_weight*(DN_v(9,0)*crLHS77 + DN_v(9,1)*crLHS376 + DN_v(9,2)*crLHS377 - crLHS1291*crLHS381 - crLHS1293*crLHS381 + crLHS836);
    rLHS(28,7)+=gauss_weight*(DN_v(9,0)*crLHS93 + DN_v(9,1)*crLHS382 + DN_v(9,2)*crLHS384 - crLHS1291*crLHS389 - crLHS1293*crLHS389 + crLHS1296 + crLHS863);
    rLHS(28,8)+=gauss_weight*(DN_v(9,0)*crLHS102 + DN_v(9,1)*crLHS390 + DN_v(9,2)*crLHS392 - crLHS1291*crLHS394 - crLHS1293*crLHS394 + crLHS889);
    rLHS(28,9)+=gauss_weight*(DN_v(9,0)*crLHS109 + DN_v(9,1)*crLHS395 + DN_v(9,2)*crLHS396 - crLHS1291*crLHS400 - crLHS1293*crLHS400 + crLHS933);
    rLHS(28,10)+=gauss_weight*(DN_v(9,0)*crLHS125 + DN_v(9,1)*crLHS401 + DN_v(9,2)*crLHS403 - crLHS1291*crLHS408 - crLHS1293*crLHS408 + crLHS1297 + crLHS957);
    rLHS(28,11)+=gauss_weight*(DN_v(9,0)*crLHS134 + DN_v(9,1)*crLHS409 + DN_v(9,2)*crLHS411 - crLHS1291*crLHS413 - crLHS1293*crLHS413 + crLHS980);
    rLHS(28,12)+=gauss_weight*(DN_v(9,0)*crLHS141 + DN_v(9,1)*crLHS414 + DN_v(9,2)*crLHS415 + crLHS1020 - crLHS1291*crLHS419 - crLHS1293*crLHS419);
    rLHS(28,13)+=gauss_weight*(DN_v(9,0)*crLHS157 + DN_v(9,1)*crLHS420 + DN_v(9,2)*crLHS422 + crLHS1041 - crLHS1291*crLHS427 - crLHS1293*crLHS427 + crLHS1298);
    rLHS(28,14)+=gauss_weight*(DN_v(9,0)*crLHS166 + DN_v(9,1)*crLHS428 + DN_v(9,2)*crLHS430 + crLHS1061 - crLHS1291*crLHS432 - crLHS1293*crLHS432);
    rLHS(28,15)+=gauss_weight*(DN_v(9,0)*crLHS173 + DN_v(9,1)*crLHS433 + DN_v(9,2)*crLHS434 + crLHS1097 - crLHS1291*crLHS438 - crLHS1293*crLHS438);
    rLHS(28,16)+=gauss_weight*(DN_v(9,0)*crLHS189 + DN_v(9,1)*crLHS439 + DN_v(9,2)*crLHS441 + crLHS1115 - crLHS1291*crLHS446 - crLHS1293*crLHS446 + crLHS1299);
    rLHS(28,17)+=gauss_weight*(DN_v(9,0)*crLHS198 + DN_v(9,1)*crLHS447 + DN_v(9,2)*crLHS449 + crLHS1132 - crLHS1291*crLHS451 - crLHS1293*crLHS451);
    rLHS(28,18)+=gauss_weight*(DN_v(9,0)*crLHS205 + DN_v(9,1)*crLHS452 + DN_v(9,2)*crLHS453 + crLHS1164 - crLHS1291*crLHS457 - crLHS1293*crLHS457);
    rLHS(28,19)+=gauss_weight*(DN_v(9,0)*crLHS221 + DN_v(9,1)*crLHS458 + DN_v(9,2)*crLHS460 + crLHS1179 - crLHS1291*crLHS465 - crLHS1293*crLHS465 + crLHS1300);
    rLHS(28,20)+=gauss_weight*(DN_v(9,0)*crLHS230 + DN_v(9,1)*crLHS466 + DN_v(9,2)*crLHS468 + crLHS1193 - crLHS1291*crLHS470 - crLHS1293*crLHS470);
    rLHS(28,21)+=gauss_weight*(DN_v(9,0)*crLHS237 + DN_v(9,1)*crLHS471 + DN_v(9,2)*crLHS472 + crLHS1221 - crLHS1291*crLHS476 - crLHS1293*crLHS476);
    rLHS(28,22)+=gauss_weight*(DN_v(9,0)*crLHS253 + DN_v(9,1)*crLHS477 + DN_v(9,2)*crLHS479 + crLHS1233 - crLHS1291*crLHS484 - crLHS1293*crLHS484 + crLHS1301);
    rLHS(28,23)+=gauss_weight*(DN_v(9,0)*crLHS262 + DN_v(9,1)*crLHS485 + DN_v(9,2)*crLHS487 + crLHS1244 - crLHS1291*crLHS489 - crLHS1293*crLHS489);
    rLHS(28,24)+=gauss_weight*(DN_v(9,0)*crLHS269 + DN_v(9,1)*crLHS490 + DN_v(9,2)*crLHS491 + crLHS1268 - crLHS1291*crLHS495 - crLHS1293*crLHS495);
    rLHS(28,25)+=gauss_weight*(DN_v(9,0)*crLHS285 + DN_v(9,1)*crLHS496 + DN_v(9,2)*crLHS498 + crLHS1277 - crLHS1291*crLHS503 - crLHS1293*crLHS503 + crLHS1302);
    rLHS(28,26)+=gauss_weight*(DN_v(9,0)*crLHS294 + DN_v(9,1)*crLHS504 + DN_v(9,2)*crLHS506 + crLHS1285 - crLHS1291*crLHS508 - crLHS1293*crLHS508);
    rLHS(28,27)+=gauss_weight*(DN_v(9,0)*crLHS301 + DN_v(9,1)*crLHS509 + DN_v(9,2)*crLHS510 - crLHS1291*crLHS514 - crLHS1293*crLHS514 + crLHS1305);
    rLHS(28,28)+=gauss_weight*(DN_v(9,0)*crLHS316 + std::pow(DN_v(9,1), 2)*crLHS9 + DN_v(9,1)*crLHS515 + DN_v(9,2)*crLHS517 - crLHS1291*crLHS522 - crLHS1293*crLHS522 + crLHS1303);
    rLHS(28,29)+=gauss_weight*(DN_v(9,0)*crLHS325 + DN_v(9,1)*crLHS523 + DN_v(9,2)*crLHS525 - crLHS1291*crLHS527 - crLHS1293*crLHS527 + crLHS1311);
    rLHS(28,30)+=gauss_weight*(crLHS1292*crLHS529 + crLHS1312 + crLHS309*crLHS529);
    rLHS(28,31)+=gauss_weight*(crLHS1292*crLHS531 + crLHS1313 + crLHS309*crLHS531);
    rLHS(28,32)+=gauss_weight*(crLHS1292*crLHS533 + crLHS1314 + crLHS309*crLHS533);
    rLHS(28,33)+=gauss_weight*(crLHS1292*crLHS535 + crLHS1315 + crLHS309*crLHS535);
    rLHS(29,0)+=gauss_weight*(DN_v(9,0)*crLHS4 + DN_v(9,1)*crLHS339 + DN_v(9,2)*crLHS536 - crLHS1291*crLHS539 - crLHS1293*crLHS539 + crLHS328);
    rLHS(29,1)+=gauss_weight*(DN_v(9,0)*crLHS31 + DN_v(9,1)*crLHS345 + DN_v(9,2)*crLHS540 - crLHS1291*crLHS542 - crLHS1293*crLHS542 + crLHS526);
    rLHS(29,2)+=gauss_weight*(DN_v(9,0)*crLHS40 + DN_v(9,1)*crLHS353 + DN_v(9,2)*crLHS543 - crLHS1291*crLHS547 - crLHS1293*crLHS547 + crLHS1294 + crLHS671);
    rLHS(29,3)+=gauss_weight*(DN_v(9,0)*crLHS47 + DN_v(9,1)*crLHS358 + DN_v(9,2)*crLHS548 - crLHS1291*crLHS553 - crLHS1293*crLHS553 + crLHS730);
    rLHS(29,4)+=gauss_weight*(DN_v(9,0)*crLHS64 + DN_v(9,1)*crLHS365 + DN_v(9,2)*crLHS554 - crLHS1291*crLHS557 - crLHS1293*crLHS557 + crLHS760);
    rLHS(29,5)+=gauss_weight*(DN_v(9,0)*crLHS72 + DN_v(9,1)*crLHS373 + DN_v(9,2)*crLHS558 - crLHS1291*crLHS562 - crLHS1293*crLHS562 + crLHS1295 + crLHS789);
    rLHS(29,6)+=gauss_weight*(DN_v(9,0)*crLHS79 + DN_v(9,1)*crLHS377 + DN_v(9,2)*crLHS563 - crLHS1291*crLHS567 - crLHS1293*crLHS567 + crLHS837);
    rLHS(29,7)+=gauss_weight*(DN_v(9,0)*crLHS96 + DN_v(9,1)*crLHS384 + DN_v(9,2)*crLHS568 - crLHS1291*crLHS571 - crLHS1293*crLHS571 + crLHS864);
    rLHS(29,8)+=gauss_weight*(DN_v(9,0)*crLHS104 + DN_v(9,1)*crLHS392 + DN_v(9,2)*crLHS572 - crLHS1291*crLHS576 - crLHS1293*crLHS576 + crLHS1296 + crLHS890);
    rLHS(29,9)+=gauss_weight*(DN_v(9,0)*crLHS111 + DN_v(9,1)*crLHS396 + DN_v(9,2)*crLHS577 - crLHS1291*crLHS581 - crLHS1293*crLHS581 + crLHS934);
    rLHS(29,10)+=gauss_weight*(DN_v(9,0)*crLHS128 + DN_v(9,1)*crLHS403 + DN_v(9,2)*crLHS582 - crLHS1291*crLHS585 - crLHS1293*crLHS585 + crLHS958);
    rLHS(29,11)+=gauss_weight*(DN_v(9,0)*crLHS136 + DN_v(9,1)*crLHS411 + DN_v(9,2)*crLHS586 - crLHS1291*crLHS590 - crLHS1293*crLHS590 + crLHS1297 + crLHS981);
    rLHS(29,12)+=gauss_weight*(DN_v(9,0)*crLHS143 + DN_v(9,1)*crLHS415 + DN_v(9,2)*crLHS591 + crLHS1021 - crLHS1291*crLHS595 - crLHS1293*crLHS595);
    rLHS(29,13)+=gauss_weight*(DN_v(9,0)*crLHS160 + DN_v(9,1)*crLHS422 + DN_v(9,2)*crLHS596 + crLHS1042 - crLHS1291*crLHS599 - crLHS1293*crLHS599);
    rLHS(29,14)+=gauss_weight*(DN_v(9,0)*crLHS168 + DN_v(9,1)*crLHS430 + DN_v(9,2)*crLHS600 + crLHS1062 - crLHS1291*crLHS604 - crLHS1293*crLHS604 + crLHS1298);
    rLHS(29,15)+=gauss_weight*(DN_v(9,0)*crLHS175 + DN_v(9,1)*crLHS434 + DN_v(9,2)*crLHS605 + crLHS1098 - crLHS1291*crLHS609 - crLHS1293*crLHS609);
    rLHS(29,16)+=gauss_weight*(DN_v(9,0)*crLHS192 + DN_v(9,1)*crLHS441 + DN_v(9,2)*crLHS610 + crLHS1116 - crLHS1291*crLHS613 - crLHS1293*crLHS613);
    rLHS(29,17)+=gauss_weight*(DN_v(9,0)*crLHS200 + DN_v(9,1)*crLHS449 + DN_v(9,2)*crLHS614 + crLHS1133 - crLHS1291*crLHS618 - crLHS1293*crLHS618 + crLHS1299);
    rLHS(29,18)+=gauss_weight*(DN_v(9,0)*crLHS207 + DN_v(9,1)*crLHS453 + DN_v(9,2)*crLHS619 + crLHS1165 - crLHS1291*crLHS623 - crLHS1293*crLHS623);
    rLHS(29,19)+=gauss_weight*(DN_v(9,0)*crLHS224 + DN_v(9,1)*crLHS460 + DN_v(9,2)*crLHS624 + crLHS1180 - crLHS1291*crLHS627 - crLHS1293*crLHS627);
    rLHS(29,20)+=gauss_weight*(DN_v(9,0)*crLHS232 + DN_v(9,1)*crLHS468 + DN_v(9,2)*crLHS628 + crLHS1194 - crLHS1291*crLHS632 - crLHS1293*crLHS632 + crLHS1300);
    rLHS(29,21)+=gauss_weight*(DN_v(9,0)*crLHS239 + DN_v(9,1)*crLHS472 + DN_v(9,2)*crLHS633 + crLHS1222 - crLHS1291*crLHS637 - crLHS1293*crLHS637);
    rLHS(29,22)+=gauss_weight*(DN_v(9,0)*crLHS256 + DN_v(9,1)*crLHS479 + DN_v(9,2)*crLHS638 + crLHS1234 - crLHS1291*crLHS641 - crLHS1293*crLHS641);
    rLHS(29,23)+=gauss_weight*(DN_v(9,0)*crLHS264 + DN_v(9,1)*crLHS487 + DN_v(9,2)*crLHS642 + crLHS1245 - crLHS1291*crLHS646 - crLHS1293*crLHS646 + crLHS1301);
    rLHS(29,24)+=gauss_weight*(DN_v(9,0)*crLHS271 + DN_v(9,1)*crLHS491 + DN_v(9,2)*crLHS647 + crLHS1269 - crLHS1291*crLHS651 - crLHS1293*crLHS651);
    rLHS(29,25)+=gauss_weight*(DN_v(9,0)*crLHS288 + DN_v(9,1)*crLHS498 + DN_v(9,2)*crLHS652 + crLHS1278 - crLHS1291*crLHS655 - crLHS1293*crLHS655);
    rLHS(29,26)+=gauss_weight*(DN_v(9,0)*crLHS296 + DN_v(9,1)*crLHS506 + DN_v(9,2)*crLHS656 + crLHS1286 - crLHS1291*crLHS660 - crLHS1293*crLHS660 + crLHS1302);
    rLHS(29,27)+=gauss_weight*(DN_v(9,0)*crLHS303 + DN_v(9,1)*crLHS510 + DN_v(9,2)*crLHS661 - crLHS1291*crLHS665 - crLHS1293*crLHS665 + crLHS1306);
    rLHS(29,28)+=gauss_weight*(DN_v(9,0)*crLHS319 + DN_v(9,1)*crLHS517 + DN_v(9,2)*crLHS666 - crLHS1291*crLHS669 - crLHS1293*crLHS669 + crLHS1311);
    rLHS(29,29)+=gauss_weight*(DN_v(9,0)*crLHS327 + DN_v(9,1)*crLHS525 + std::pow(DN_v(9,2), 2)*crLHS9 + DN_v(9,2)*crLHS670 - crLHS1291*crLHS674 - crLHS1293*crLHS674 + crLHS1303);
    rLHS(29,30)+=gauss_weight*(crLHS1292*crLHS676 + crLHS1316 + crLHS309*crLHS676);
    rLHS(29,31)+=gauss_weight*(crLHS1292*crLHS678 + crLHS1317 + crLHS309*crLHS678);
    rLHS(29,32)+=gauss_weight*(crLHS1292*crLHS680 + crLHS1318 + crLHS309*crLHS680);
    rLHS(29,33)+=gauss_weight*(crLHS1292*crLHS682 + crLHS1319 + crLHS309*crLHS682);
    rLHS(30,0)+=-gauss_weight*(crLHS19*crLHS331 + crLHS330 + crLHS342*crLHS529 + crLHS539*crLHS676);
    rLHS(30,1)+=-gauss_weight*(crLHS331*crLHS35 + crLHS350*crLHS529 + crLHS528 + crLHS542*crLHS676);
    rLHS(30,2)+=-gauss_weight*(crLHS331*crLHS42 + crLHS356*crLHS529 + crLHS547*crLHS676 + crLHS675);
    rLHS(30,3)+=-gauss_weight*(crLHS331*crLHS56 + crLHS362*crLHS529 + crLHS553*crLHS676 + crLHS731);
    rLHS(30,4)+=-gauss_weight*(crLHS331*crLHS67 + crLHS370*crLHS529 + crLHS557*crLHS676 + crLHS761);
    rLHS(30,5)+=-gauss_weight*(crLHS331*crLHS74 + crLHS375*crLHS529 + crLHS562*crLHS676 + crLHS790);
    rLHS(30,6)+=-gauss_weight*(crLHS331*crLHS88 + crLHS381*crLHS529 + crLHS567*crLHS676 + crLHS838);
    rLHS(30,7)+=-gauss_weight*(crLHS331*crLHS99 + crLHS389*crLHS529 + crLHS571*crLHS676 + crLHS865);
    rLHS(30,8)+=-gauss_weight*(crLHS106*crLHS331 + crLHS394*crLHS529 + crLHS576*crLHS676 + crLHS891);
    rLHS(30,9)+=-gauss_weight*(crLHS120*crLHS331 + crLHS400*crLHS529 + crLHS581*crLHS676 + crLHS935);
    rLHS(30,10)+=-gauss_weight*(crLHS131*crLHS331 + crLHS408*crLHS529 + crLHS585*crLHS676 + crLHS959);
    rLHS(30,11)+=-gauss_weight*(crLHS138*crLHS331 + crLHS413*crLHS529 + crLHS590*crLHS676 + crLHS982);
    rLHS(30,12)+=-gauss_weight*(crLHS1022 + crLHS152*crLHS331 + crLHS419*crLHS529 + crLHS595*crLHS676);
    rLHS(30,13)+=-gauss_weight*(crLHS1043 + crLHS163*crLHS331 + crLHS427*crLHS529 + crLHS599*crLHS676);
    rLHS(30,14)+=-gauss_weight*(crLHS1063 + crLHS170*crLHS331 + crLHS432*crLHS529 + crLHS604*crLHS676);
    rLHS(30,15)+=-gauss_weight*(crLHS1099 + crLHS184*crLHS331 + crLHS438*crLHS529 + crLHS609*crLHS676);
    rLHS(30,16)+=-gauss_weight*(crLHS1117 + crLHS195*crLHS331 + crLHS446*crLHS529 + crLHS613*crLHS676);
    rLHS(30,17)+=-gauss_weight*(crLHS1134 + crLHS202*crLHS331 + crLHS451*crLHS529 + crLHS618*crLHS676);
    rLHS(30,18)+=-gauss_weight*(crLHS1166 + crLHS216*crLHS331 + crLHS457*crLHS529 + crLHS623*crLHS676);
    rLHS(30,19)+=-gauss_weight*(crLHS1181 + crLHS227*crLHS331 + crLHS465*crLHS529 + crLHS627*crLHS676);
    rLHS(30,20)+=-gauss_weight*(crLHS1195 + crLHS234*crLHS331 + crLHS470*crLHS529 + crLHS632*crLHS676);
    rLHS(30,21)+=-gauss_weight*(crLHS1223 + crLHS248*crLHS331 + crLHS476*crLHS529 + crLHS637*crLHS676);
    rLHS(30,22)+=-gauss_weight*(crLHS1235 + crLHS259*crLHS331 + crLHS484*crLHS529 + crLHS641*crLHS676);
    rLHS(30,23)+=-gauss_weight*(crLHS1246 + crLHS266*crLHS331 + crLHS489*crLHS529 + crLHS646*crLHS676);
    rLHS(30,24)+=-gauss_weight*(crLHS1270 + crLHS280*crLHS331 + crLHS495*crLHS529 + crLHS651*crLHS676);
    rLHS(30,25)+=-gauss_weight*(crLHS1279 + crLHS291*crLHS331 + crLHS503*crLHS529 + crLHS655*crLHS676);
    rLHS(30,26)+=-gauss_weight*(crLHS1287 + crLHS298*crLHS331 + crLHS508*crLHS529 + crLHS660*crLHS676);
    rLHS(30,27)+=-gauss_weight*(crLHS1307 + crLHS311*crLHS331 + crLHS514*crLHS529 + crLHS665*crLHS676);
    rLHS(30,28)+=-gauss_weight*(crLHS1312 + crLHS322*crLHS331 + crLHS522*crLHS529 + crLHS669*crLHS676);
    rLHS(30,29)+=-gauss_weight*(crLHS1316 + crLHS329*crLHS331 + crLHS527*crLHS529 + crLHS674*crLHS676);
    rLHS(30,30)+=crLHS1320*(std::pow(DN_p(0,0), 2) + std::pow(DN_p(0,1), 2) + std::pow(DN_p(0,2), 2));
    rLHS(30,31)+=crLHS1321;
    rLHS(30,32)+=crLHS1322;
    rLHS(30,33)+=crLHS1323;
    rLHS(31,0)+=-gauss_weight*(crLHS19*crLHS333 + crLHS332 + crLHS342*crLHS531 + crLHS539*crLHS678);
    rLHS(31,1)+=-gauss_weight*(crLHS333*crLHS35 + crLHS350*crLHS531 + crLHS530 + crLHS542*crLHS678);
    rLHS(31,2)+=-gauss_weight*(crLHS333*crLHS42 + crLHS356*crLHS531 + crLHS547*crLHS678 + crLHS677);
    rLHS(31,3)+=-gauss_weight*(crLHS333*crLHS56 + crLHS362*crLHS531 + crLHS553*crLHS678 + crLHS732);
    rLHS(31,4)+=-gauss_weight*(crLHS333*crLHS67 + crLHS370*crLHS531 + crLHS557*crLHS678 + crLHS762);
    rLHS(31,5)+=-gauss_weight*(crLHS333*crLHS74 + crLHS375*crLHS531 + crLHS562*crLHS678 + crLHS791);
    rLHS(31,6)+=-gauss_weight*(crLHS333*crLHS88 + crLHS381*crLHS531 + crLHS567*crLHS678 + crLHS839);
    rLHS(31,7)+=-gauss_weight*(crLHS333*crLHS99 + crLHS389*crLHS531 + crLHS571*crLHS678 + crLHS866);
    rLHS(31,8)+=-gauss_weight*(crLHS106*crLHS333 + crLHS394*crLHS531 + crLHS576*crLHS678 + crLHS892);
    rLHS(31,9)+=-gauss_weight*(crLHS120*crLHS333 + crLHS400*crLHS531 + crLHS581*crLHS678 + crLHS936);
    rLHS(31,10)+=-gauss_weight*(crLHS131*crLHS333 + crLHS408*crLHS531 + crLHS585*crLHS678 + crLHS960);
    rLHS(31,11)+=-gauss_weight*(crLHS138*crLHS333 + crLHS413*crLHS531 + crLHS590*crLHS678 + crLHS983);
    rLHS(31,12)+=-gauss_weight*(crLHS1023 + crLHS152*crLHS333 + crLHS419*crLHS531 + crLHS595*crLHS678);
    rLHS(31,13)+=-gauss_weight*(crLHS1044 + crLHS163*crLHS333 + crLHS427*crLHS531 + crLHS599*crLHS678);
    rLHS(31,14)+=-gauss_weight*(crLHS1064 + crLHS170*crLHS333 + crLHS432*crLHS531 + crLHS604*crLHS678);
    rLHS(31,15)+=-gauss_weight*(crLHS1100 + crLHS184*crLHS333 + crLHS438*crLHS531 + crLHS609*crLHS678);
    rLHS(31,16)+=-gauss_weight*(crLHS1118 + crLHS195*crLHS333 + crLHS446*crLHS531 + crLHS613*crLHS678);
    rLHS(31,17)+=-gauss_weight*(crLHS1135 + crLHS202*crLHS333 + crLHS451*crLHS531 + crLHS618*crLHS678);
    rLHS(31,18)+=-gauss_weight*(crLHS1167 + crLHS216*crLHS333 + crLHS457*crLHS531 + crLHS623*crLHS678);
    rLHS(31,19)+=-gauss_weight*(crLHS1182 + crLHS227*crLHS333 + crLHS465*crLHS531 + crLHS627*crLHS678);
    rLHS(31,20)+=-gauss_weight*(crLHS1196 + crLHS234*crLHS333 + crLHS470*crLHS531 + crLHS632*crLHS678);
    rLHS(31,21)+=-gauss_weight*(crLHS1224 + crLHS248*crLHS333 + crLHS476*crLHS531 + crLHS637*crLHS678);
    rLHS(31,22)+=-gauss_weight*(crLHS1236 + crLHS259*crLHS333 + crLHS484*crLHS531 + crLHS641*crLHS678);
    rLHS(31,23)+=-gauss_weight*(crLHS1247 + crLHS266*crLHS333 + crLHS489*crLHS531 + crLHS646*crLHS678);
    rLHS(31,24)+=-gauss_weight*(crLHS1271 + crLHS280*crLHS333 + crLHS495*crLHS531 + crLHS651*crLHS678);
    rLHS(31,25)+=-gauss_weight*(crLHS1280 + crLHS291*crLHS333 + crLHS503*crLHS531 + crLHS655*crLHS678);
    rLHS(31,26)+=-gauss_weight*(crLHS1288 + crLHS298*crLHS333 + crLHS508*crLHS531 + crLHS660*crLHS678);
    rLHS(31,27)+=-gauss_weight*(crLHS1308 + crLHS311*crLHS333 + crLHS514*crLHS531 + crLHS665*crLHS678);
    rLHS(31,28)+=-gauss_weight*(crLHS1313 + crLHS322*crLHS333 + crLHS522*crLHS531 + crLHS669*crLHS678);
    rLHS(31,29)+=-gauss_weight*(crLHS1317 + crLHS329*crLHS333 + crLHS527*crLHS531 + crLHS674*crLHS678);
    rLHS(31,30)+=crLHS1321;
    rLHS(31,31)+=crLHS1320*(std::pow(DN_p(1,0), 2) + std::pow(DN_p(1,1), 2) + std::pow(DN_p(1,2), 2));
    rLHS(31,32)+=crLHS1324;
    rLHS(31,33)+=crLHS1325;
    rLHS(32,0)+=-gauss_weight*(crLHS19*crLHS335 + crLHS334 + crLHS342*crLHS533 + crLHS539*crLHS680);
    rLHS(32,1)+=-gauss_weight*(crLHS335*crLHS35 + crLHS350*crLHS533 + crLHS532 + crLHS542*crLHS680);
    rLHS(32,2)+=-gauss_weight*(crLHS335*crLHS42 + crLHS356*crLHS533 + crLHS547*crLHS680 + crLHS679);
    rLHS(32,3)+=-gauss_weight*(crLHS335*crLHS56 + crLHS362*crLHS533 + crLHS553*crLHS680 + crLHS733);
    rLHS(32,4)+=-gauss_weight*(crLHS335*crLHS67 + crLHS370*crLHS533 + crLHS557*crLHS680 + crLHS763);
    rLHS(32,5)+=-gauss_weight*(crLHS335*crLHS74 + crLHS375*crLHS533 + crLHS562*crLHS680 + crLHS792);
    rLHS(32,6)+=-gauss_weight*(crLHS335*crLHS88 + crLHS381*crLHS533 + crLHS567*crLHS680 + crLHS840);
    rLHS(32,7)+=-gauss_weight*(crLHS335*crLHS99 + crLHS389*crLHS533 + crLHS571*crLHS680 + crLHS867);
    rLHS(32,8)+=-gauss_weight*(crLHS106*crLHS335 + crLHS394*crLHS533 + crLHS576*crLHS680 + crLHS893);
    rLHS(32,9)+=-gauss_weight*(crLHS120*crLHS335 + crLHS400*crLHS533 + crLHS581*crLHS680 + crLHS937);
    rLHS(32,10)+=-gauss_weight*(crLHS131*crLHS335 + crLHS408*crLHS533 + crLHS585*crLHS680 + crLHS961);
    rLHS(32,11)+=-gauss_weight*(crLHS138*crLHS335 + crLHS413*crLHS533 + crLHS590*crLHS680 + crLHS984);
    rLHS(32,12)+=-gauss_weight*(crLHS1024 + crLHS152*crLHS335 + crLHS419*crLHS533 + crLHS595*crLHS680);
    rLHS(32,13)+=-gauss_weight*(crLHS1045 + crLHS163*crLHS335 + crLHS427*crLHS533 + crLHS599*crLHS680);
    rLHS(32,14)+=-gauss_weight*(crLHS1065 + crLHS170*crLHS335 + crLHS432*crLHS533 + crLHS604*crLHS680);
    rLHS(32,15)+=-gauss_weight*(crLHS1101 + crLHS184*crLHS335 + crLHS438*crLHS533 + crLHS609*crLHS680);
    rLHS(32,16)+=-gauss_weight*(crLHS1119 + crLHS195*crLHS335 + crLHS446*crLHS533 + crLHS613*crLHS680);
    rLHS(32,17)+=-gauss_weight*(crLHS1136 + crLHS202*crLHS335 + crLHS451*crLHS533 + crLHS618*crLHS680);
    rLHS(32,18)+=-gauss_weight*(crLHS1168 + crLHS216*crLHS335 + crLHS457*crLHS533 + crLHS623*crLHS680);
    rLHS(32,19)+=-gauss_weight*(crLHS1183 + crLHS227*crLHS335 + crLHS465*crLHS533 + crLHS627*crLHS680);
    rLHS(32,20)+=-gauss_weight*(crLHS1197 + crLHS234*crLHS335 + crLHS470*crLHS533 + crLHS632*crLHS680);
    rLHS(32,21)+=-gauss_weight*(crLHS1225 + crLHS248*crLHS335 + crLHS476*crLHS533 + crLHS637*crLHS680);
    rLHS(32,22)+=-gauss_weight*(crLHS1237 + crLHS259*crLHS335 + crLHS484*crLHS533 + crLHS641*crLHS680);
    rLHS(32,23)+=-gauss_weight*(crLHS1248 + crLHS266*crLHS335 + crLHS489*crLHS533 + crLHS646*crLHS680);
    rLHS(32,24)+=-gauss_weight*(crLHS1272 + crLHS280*crLHS335 + crLHS495*crLHS533 + crLHS651*crLHS680);
    rLHS(32,25)+=-gauss_weight*(crLHS1281 + crLHS291*crLHS335 + crLHS503*crLHS533 + crLHS655*crLHS680);
    rLHS(32,26)+=-gauss_weight*(crLHS1289 + crLHS298*crLHS335 + crLHS508*crLHS533 + crLHS660*crLHS680);
    rLHS(32,27)+=-gauss_weight*(crLHS1309 + crLHS311*crLHS335 + crLHS514*crLHS533 + crLHS665*crLHS680);
    rLHS(32,28)+=-gauss_weight*(crLHS1314 + crLHS322*crLHS335 + crLHS522*crLHS533 + crLHS669*crLHS680);
    rLHS(32,29)+=-gauss_weight*(crLHS1318 + crLHS329*crLHS335 + crLHS527*crLHS533 + crLHS674*crLHS680);
    rLHS(32,30)+=crLHS1322;
    rLHS(32,31)+=crLHS1324;
    rLHS(32,32)+=crLHS1320*(std::pow(DN_p(2,0), 2) + std::pow(DN_p(2,1), 2) + std::pow(DN_p(2,2), 2));
    rLHS(32,33)+=crLHS1326;
    rLHS(33,0)+=-gauss_weight*(crLHS19*crLHS337 + crLHS336 + crLHS342*crLHS535 + crLHS539*crLHS682);
    rLHS(33,1)+=-gauss_weight*(crLHS337*crLHS35 + crLHS350*crLHS535 + crLHS534 + crLHS542*crLHS682);
    rLHS(33,2)+=-gauss_weight*(crLHS337*crLHS42 + crLHS356*crLHS535 + crLHS547*crLHS682 + crLHS681);
    rLHS(33,3)+=-gauss_weight*(crLHS337*crLHS56 + crLHS362*crLHS535 + crLHS553*crLHS682 + crLHS734);
    rLHS(33,4)+=-gauss_weight*(crLHS337*crLHS67 + crLHS370*crLHS535 + crLHS557*crLHS682 + crLHS764);
    rLHS(33,5)+=-gauss_weight*(crLHS337*crLHS74 + crLHS375*crLHS535 + crLHS562*crLHS682 + crLHS793);
    rLHS(33,6)+=-gauss_weight*(crLHS337*crLHS88 + crLHS381*crLHS535 + crLHS567*crLHS682 + crLHS841);
    rLHS(33,7)+=-gauss_weight*(crLHS337*crLHS99 + crLHS389*crLHS535 + crLHS571*crLHS682 + crLHS868);
    rLHS(33,8)+=-gauss_weight*(crLHS106*crLHS337 + crLHS394*crLHS535 + crLHS576*crLHS682 + crLHS894);
    rLHS(33,9)+=-gauss_weight*(crLHS120*crLHS337 + crLHS400*crLHS535 + crLHS581*crLHS682 + crLHS938);
    rLHS(33,10)+=-gauss_weight*(crLHS131*crLHS337 + crLHS408*crLHS535 + crLHS585*crLHS682 + crLHS962);
    rLHS(33,11)+=-gauss_weight*(crLHS138*crLHS337 + crLHS413*crLHS535 + crLHS590*crLHS682 + crLHS985);
    rLHS(33,12)+=-gauss_weight*(crLHS1025 + crLHS152*crLHS337 + crLHS419*crLHS535 + crLHS595*crLHS682);
    rLHS(33,13)+=-gauss_weight*(crLHS1046 + crLHS163*crLHS337 + crLHS427*crLHS535 + crLHS599*crLHS682);
    rLHS(33,14)+=-gauss_weight*(crLHS1066 + crLHS170*crLHS337 + crLHS432*crLHS535 + crLHS604*crLHS682);
    rLHS(33,15)+=-gauss_weight*(crLHS1102 + crLHS184*crLHS337 + crLHS438*crLHS535 + crLHS609*crLHS682);
    rLHS(33,16)+=-gauss_weight*(crLHS1120 + crLHS195*crLHS337 + crLHS446*crLHS535 + crLHS613*crLHS682);
    rLHS(33,17)+=-gauss_weight*(crLHS1137 + crLHS202*crLHS337 + crLHS451*crLHS535 + crLHS618*crLHS682);
    rLHS(33,18)+=-gauss_weight*(crLHS1169 + crLHS216*crLHS337 + crLHS457*crLHS535 + crLHS623*crLHS682);
    rLHS(33,19)+=-gauss_weight*(crLHS1184 + crLHS227*crLHS337 + crLHS465*crLHS535 + crLHS627*crLHS682);
    rLHS(33,20)+=-gauss_weight*(crLHS1198 + crLHS234*crLHS337 + crLHS470*crLHS535 + crLHS632*crLHS682);
    rLHS(33,21)+=-gauss_weight*(crLHS1226 + crLHS248*crLHS337 + crLHS476*crLHS535 + crLHS637*crLHS682);
    rLHS(33,22)+=-gauss_weight*(crLHS1238 + crLHS259*crLHS337 + crLHS484*crLHS535 + crLHS641*crLHS682);
    rLHS(33,23)+=-gauss_weight*(crLHS1249 + crLHS266*crLHS337 + crLHS489*crLHS535 + crLHS646*crLHS682);
    rLHS(33,24)+=-gauss_weight*(crLHS1273 + crLHS280*crLHS337 + crLHS495*crLHS535 + crLHS651*crLHS682);
    rLHS(33,25)+=-gauss_weight*(crLHS1282 + crLHS291*crLHS337 + crLHS503*crLHS535 + crLHS655*crLHS682);
    rLHS(33,26)+=-gauss_weight*(crLHS1290 + crLHS298*crLHS337 + crLHS508*crLHS535 + crLHS660*crLHS682);
    rLHS(33,27)+=-gauss_weight*(crLHS1310 + crLHS311*crLHS337 + crLHS514*crLHS535 + crLHS665*crLHS682);
    rLHS(33,28)+=-gauss_weight*(crLHS1315 + crLHS322*crLHS337 + crLHS522*crLHS535 + crLHS669*crLHS682);
    rLHS(33,29)+=-gauss_weight*(crLHS1319 + crLHS329*crLHS337 + crLHS527*crLHS535 + crLHS674*crLHS682);
    rLHS(33,30)+=crLHS1323;
    rLHS(33,31)+=crLHS1325;
    rLHS(33,32)+=crLHS1326;
    rLHS(33,33)+=crLHS1320*(std::pow(DN_p(3,0), 2) + std::pow(DN_p(3,1), 2) + std::pow(DN_p(3,2), 2));

}

template <>
void IncompressibleNavierStokesP2P1Continuous<2>::AddGaussPointRightHandSideContribution(
    const ElementDataContainer& rData,
    VectorType& rRHS)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

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
    const double crRHS7 = DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1);
    const double crRHS8 = crRHS3 + crRHS7;
    const double crRHS9 = rho*stab_c2*std::sqrt(std::pow(crRHS4, 2) + std::pow(crRHS5, 2));
    const double crRHS10 = crRHS8*(crRHS9*h/stab_c1 + mu);
    const double crRHS11 = 1.0*C(0,0);
    const double crRHS12 = DDN_v[0](0,0)*r_v(0,0);
    const double crRHS13 = DDN_v[1](0,0)*r_v(1,0);
    const double crRHS14 = DDN_v[2](0,0)*r_v(2,0);
    const double crRHS15 = DDN_v[3](0,0)*r_v(3,0);
    const double crRHS16 = DDN_v[4](0,0)*r_v(4,0);
    const double crRHS17 = DDN_v[5](0,0)*r_v(5,0);
    const double crRHS18 = 1.0*C(0,1);
    const double crRHS19 = DDN_v[0](0,1)*r_v(0,1);
    const double crRHS20 = DDN_v[1](0,1)*r_v(1,1);
    const double crRHS21 = DDN_v[2](0,1)*r_v(2,1);
    const double crRHS22 = DDN_v[3](0,1)*r_v(3,1);
    const double crRHS23 = DDN_v[4](0,1)*r_v(4,1);
    const double crRHS24 = DDN_v[5](0,1)*r_v(5,1);
    const double crRHS25 = 1.0*C(0,2);
    const double crRHS26 = 1.0*C(1,2);
    const double crRHS27 = DDN_v[0](0,0)*r_v(0,1) + DDN_v[0](0,1)*r_v(0,0);
    const double crRHS28 = DDN_v[1](0,0)*r_v(1,1) + DDN_v[1](0,1)*r_v(1,0);
    const double crRHS29 = DDN_v[2](0,0)*r_v(2,1) + DDN_v[2](0,1)*r_v(2,0);
    const double crRHS30 = DDN_v[3](0,0)*r_v(3,1) + DDN_v[3](0,1)*r_v(3,0);
    const double crRHS31 = DDN_v[4](0,0)*r_v(4,1) + DDN_v[4](0,1)*r_v(4,0);
    const double crRHS32 = DDN_v[5](0,0)*r_v(5,1) + DDN_v[5](0,1)*r_v(5,0);
    const double crRHS33 = 1.0*C(2,2);
    const double crRHS34 = 1.0/(crRHS9/h + dyn_tau*rho/rData.DeltaTime + mu*stab_c1/std::pow(h, 2));
    const double crRHS35 = crRHS34*(-DN_p(0,0)*r_p[0] - DN_p(1,0)*r_p[1] - DN_p(2,0)*r_p[2] + crRHS1 + crRHS11*crRHS12 + crRHS11*crRHS13 + crRHS11*crRHS14 + crRHS11*crRHS15 + crRHS11*crRHS16 + crRHS11*crRHS17 + crRHS12*crRHS25 + crRHS13*crRHS25 + crRHS14*crRHS25 + crRHS15*crRHS25 + crRHS16*crRHS25 + crRHS17*crRHS25 + crRHS18*crRHS19 + crRHS18*crRHS20 + crRHS18*crRHS21 + crRHS18*crRHS22 + crRHS18*crRHS23 + crRHS18*crRHS24 + crRHS19*crRHS26 - crRHS2 + crRHS20*crRHS26 + crRHS21*crRHS26 + crRHS22*crRHS26 + crRHS23*crRHS26 + crRHS24*crRHS26 + crRHS25*crRHS27 + crRHS25*crRHS28 + crRHS25*crRHS29 + crRHS25*crRHS30 + crRHS25*crRHS31 + crRHS25*crRHS32 + crRHS27*crRHS33 + crRHS28*crRHS33 + crRHS29*crRHS33 + crRHS30*crRHS33 + crRHS31*crRHS33 + crRHS32*crRHS33 - crRHS6);
    const double crRHS36 = rho*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1));
    const double crRHS37 = N_v[0]*crRHS36;
    const double crRHS38 = rho*(DN_v(0,0)*crRHS4 + DN_v(0,1)*crRHS5);
    const double crRHS39 = rho*(N_v[0]*r_f(0,1) + N_v[1]*r_f(1,1) + N_v[2]*r_f(2,1) + N_v[3]*r_f(3,1) + N_v[4]*r_f(4,1) + N_v[5]*r_f(5,1));
    const double crRHS40 = rho*(N_v[0]*(rData.BDF0*r_v(0,1) + rData.BDF1*r_vn(0,1) + rData.BDF2*r_vnn(0,1)) + N_v[1]*(rData.BDF0*r_v(1,1) + rData.BDF1*r_vn(1,1) + rData.BDF2*r_vnn(1,1)) + N_v[2]*(rData.BDF0*r_v(2,1) + rData.BDF1*r_vn(2,1) + rData.BDF2*r_vnn(2,1)) + N_v[3]*(rData.BDF0*r_v(3,1) + rData.BDF1*r_vn(3,1) + rData.BDF2*r_vnn(3,1)) + N_v[4]*(rData.BDF0*r_v(4,1) + rData.BDF1*r_vn(4,1) + rData.BDF2*r_vnn(4,1)) + N_v[5]*(rData.BDF0*r_v(5,1) + rData.BDF1*r_vn(5,1) + rData.BDF2*r_vnn(5,1)));
    const double crRHS41 = rho*(crRHS4*(DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1)) + crRHS5*crRHS7);
    const double crRHS42 = DDN_v[0](1,0)*r_v(0,0);
    const double crRHS43 = DDN_v[1](1,0)*r_v(1,0);
    const double crRHS44 = DDN_v[2](1,0)*r_v(2,0);
    const double crRHS45 = DDN_v[3](1,0)*r_v(3,0);
    const double crRHS46 = DDN_v[4](1,0)*r_v(4,0);
    const double crRHS47 = DDN_v[5](1,0)*r_v(5,0);
    const double crRHS48 = 1.0*C(1,1);
    const double crRHS49 = DDN_v[0](1,1)*r_v(0,1);
    const double crRHS50 = DDN_v[1](1,1)*r_v(1,1);
    const double crRHS51 = DDN_v[2](1,1)*r_v(2,1);
    const double crRHS52 = DDN_v[3](1,1)*r_v(3,1);
    const double crRHS53 = DDN_v[4](1,1)*r_v(4,1);
    const double crRHS54 = DDN_v[5](1,1)*r_v(5,1);
    const double crRHS55 = DDN_v[0](1,0)*r_v(0,1) + DDN_v[0](1,1)*r_v(0,0);
    const double crRHS56 = DDN_v[1](1,0)*r_v(1,1) + DDN_v[1](1,1)*r_v(1,0);
    const double crRHS57 = DDN_v[2](1,0)*r_v(2,1) + DDN_v[2](1,1)*r_v(2,0);
    const double crRHS58 = DDN_v[3](1,0)*r_v(3,1) + DDN_v[3](1,1)*r_v(3,0);
    const double crRHS59 = DDN_v[4](1,0)*r_v(4,1) + DDN_v[4](1,1)*r_v(4,0);
    const double crRHS60 = DDN_v[5](1,0)*r_v(5,1) + DDN_v[5](1,1)*r_v(5,0);
    const double crRHS61 = crRHS34*(-DN_p(0,1)*r_p[0] - DN_p(1,1)*r_p[1] - DN_p(2,1)*r_p[2] + crRHS18*crRHS42 + crRHS18*crRHS43 + crRHS18*crRHS44 + crRHS18*crRHS45 + crRHS18*crRHS46 + crRHS18*crRHS47 + crRHS25*crRHS42 + crRHS25*crRHS43 + crRHS25*crRHS44 + crRHS25*crRHS45 + crRHS25*crRHS46 + crRHS25*crRHS47 + crRHS26*crRHS49 + crRHS26*crRHS50 + crRHS26*crRHS51 + crRHS26*crRHS52 + crRHS26*crRHS53 + crRHS26*crRHS54 + crRHS26*crRHS55 + crRHS26*crRHS56 + crRHS26*crRHS57 + crRHS26*crRHS58 + crRHS26*crRHS59 + crRHS26*crRHS60 + crRHS33*crRHS55 + crRHS33*crRHS56 + crRHS33*crRHS57 + crRHS33*crRHS58 + crRHS33*crRHS59 + crRHS33*crRHS60 + crRHS39 - crRHS40 - crRHS41 + crRHS48*crRHS49 + crRHS48*crRHS50 + crRHS48*crRHS51 + crRHS48*crRHS52 + crRHS48*crRHS53 + crRHS48*crRHS54);
    const double crRHS62 = N_v[1]*crRHS36;
    const double crRHS63 = rho*(DN_v(1,0)*crRHS4 + DN_v(1,1)*crRHS5);
    const double crRHS64 = N_v[2]*crRHS36;
    const double crRHS65 = rho*(DN_v(2,0)*crRHS4 + DN_v(2,1)*crRHS5);
    const double crRHS66 = N_v[3]*crRHS36;
    const double crRHS67 = rho*(DN_v(3,0)*crRHS4 + DN_v(3,1)*crRHS5);
    const double crRHS68 = N_v[4]*crRHS36;
    const double crRHS69 = rho*(DN_v(4,0)*crRHS4 + DN_v(4,1)*crRHS5);
    const double crRHS70 = N_v[5]*crRHS36;
    const double crRHS71 = rho*(DN_v(5,0)*crRHS4 + DN_v(5,1)*crRHS5);
    rRHS[0]+=-gauss_weight*(-DN_v(0,0)*crRHS0 + DN_v(0,0)*crRHS10 + DN_v(0,0)*r_stress[0] + DN_v(0,1)*r_stress[2] - N_v[0]*crRHS1 + N_v[0]*crRHS2 + N_v[0]*crRHS6 - crRHS35*crRHS37 - crRHS35*crRHS38);
    rRHS[1]+=-gauss_weight*(DN_v(0,0)*r_stress[2] - DN_v(0,1)*crRHS0 + DN_v(0,1)*crRHS10 + DN_v(0,1)*r_stress[1] - N_v[0]*crRHS39 + N_v[0]*crRHS40 + N_v[0]*crRHS41 - crRHS37*crRHS61 - crRHS38*crRHS61);
    rRHS[2]+=-gauss_weight*(-DN_v(1,0)*crRHS0 + DN_v(1,0)*crRHS10 + DN_v(1,0)*r_stress[0] + DN_v(1,1)*r_stress[2] - N_v[1]*crRHS1 + N_v[1]*crRHS2 + N_v[1]*crRHS6 - crRHS35*crRHS62 - crRHS35*crRHS63);
    rRHS[3]+=-gauss_weight*(DN_v(1,0)*r_stress[2] - DN_v(1,1)*crRHS0 + DN_v(1,1)*crRHS10 + DN_v(1,1)*r_stress[1] - N_v[1]*crRHS39 + N_v[1]*crRHS40 + N_v[1]*crRHS41 - crRHS61*crRHS62 - crRHS61*crRHS63);
    rRHS[4]+=-gauss_weight*(-DN_v(2,0)*crRHS0 + DN_v(2,0)*crRHS10 + DN_v(2,0)*r_stress[0] + DN_v(2,1)*r_stress[2] - N_v[2]*crRHS1 + N_v[2]*crRHS2 + N_v[2]*crRHS6 - crRHS35*crRHS64 - crRHS35*crRHS65);
    rRHS[5]+=-gauss_weight*(DN_v(2,0)*r_stress[2] - DN_v(2,1)*crRHS0 + DN_v(2,1)*crRHS10 + DN_v(2,1)*r_stress[1] - N_v[2]*crRHS39 + N_v[2]*crRHS40 + N_v[2]*crRHS41 - crRHS61*crRHS64 - crRHS61*crRHS65);
    rRHS[6]+=-gauss_weight*(-DN_v(3,0)*crRHS0 + DN_v(3,0)*crRHS10 + DN_v(3,0)*r_stress[0] + DN_v(3,1)*r_stress[2] - N_v[3]*crRHS1 + N_v[3]*crRHS2 + N_v[3]*crRHS6 - crRHS35*crRHS66 - crRHS35*crRHS67);
    rRHS[7]+=-gauss_weight*(DN_v(3,0)*r_stress[2] - DN_v(3,1)*crRHS0 + DN_v(3,1)*crRHS10 + DN_v(3,1)*r_stress[1] - N_v[3]*crRHS39 + N_v[3]*crRHS40 + N_v[3]*crRHS41 - crRHS61*crRHS66 - crRHS61*crRHS67);
    rRHS[8]+=-gauss_weight*(-DN_v(4,0)*crRHS0 + DN_v(4,0)*crRHS10 + DN_v(4,0)*r_stress[0] + DN_v(4,1)*r_stress[2] - N_v[4]*crRHS1 + N_v[4]*crRHS2 + N_v[4]*crRHS6 - crRHS35*crRHS68 - crRHS35*crRHS69);
    rRHS[9]+=-gauss_weight*(DN_v(4,0)*r_stress[2] - DN_v(4,1)*crRHS0 + DN_v(4,1)*crRHS10 + DN_v(4,1)*r_stress[1] - N_v[4]*crRHS39 + N_v[4]*crRHS40 + N_v[4]*crRHS41 - crRHS61*crRHS68 - crRHS61*crRHS69);
    rRHS[10]+=-gauss_weight*(-DN_v(5,0)*crRHS0 + DN_v(5,0)*crRHS10 + DN_v(5,0)*r_stress[0] + DN_v(5,1)*r_stress[2] - N_v[5]*crRHS1 + N_v[5]*crRHS2 + N_v[5]*crRHS6 - crRHS35*crRHS70 - crRHS35*crRHS71);
    rRHS[11]+=-gauss_weight*(DN_v(5,0)*r_stress[2] - DN_v(5,1)*crRHS0 + DN_v(5,1)*crRHS10 + DN_v(5,1)*r_stress[1] - N_v[5]*crRHS39 + N_v[5]*crRHS40 + N_v[5]*crRHS41 - crRHS61*crRHS70 - crRHS61*crRHS71);
    rRHS[12]+=gauss_weight*(DN_p(0,0)*crRHS35 + DN_p(0,1)*crRHS61 - N_p[0]*crRHS8);
    rRHS[13]+=gauss_weight*(DN_p(1,0)*crRHS35 + DN_p(1,1)*crRHS61 - N_p[1]*crRHS8);
    rRHS[14]+=gauss_weight*(DN_p(2,0)*crRHS35 + DN_p(2,1)*crRHS61 - N_p[2]*crRHS8);

}

template <>
void IncompressibleNavierStokesP2P1Continuous<3>::AddGaussPointRightHandSideContribution(
    const ElementDataContainer& rData,
    VectorType& rRHS)
{
    // Get material data
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

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
    const double crRHS0 = N_p[0]*r_p[0] + N_p[1]*r_p[1] + N_p[2]*r_p[2] + N_p[3]*r_p[3];
    const double crRHS1 = rho*(N_v[0]*r_f(0,0) + N_v[1]*r_f(1,0) + N_v[2]*r_f(2,0) + N_v[3]*r_f(3,0) + N_v[4]*r_f(4,0) + N_v[5]*r_f(5,0) + N_v[6]*r_f(6,0) + N_v[7]*r_f(7,0) + N_v[8]*r_f(8,0) + N_v[9]*r_f(9,0));
    const double crRHS2 = rho*(N_v[0]*(rData.BDF0*r_v(0,0) + rData.BDF1*r_vn(0,0) + rData.BDF2*r_vnn(0,0)) + N_v[1]*(rData.BDF0*r_v(1,0) + rData.BDF1*r_vn(1,0) + rData.BDF2*r_vnn(1,0)) + N_v[2]*(rData.BDF0*r_v(2,0) + rData.BDF1*r_vn(2,0) + rData.BDF2*r_vnn(2,0)) + N_v[3]*(rData.BDF0*r_v(3,0) + rData.BDF1*r_vn(3,0) + rData.BDF2*r_vnn(3,0)) + N_v[4]*(rData.BDF0*r_v(4,0) + rData.BDF1*r_vn(4,0) + rData.BDF2*r_vnn(4,0)) + N_v[5]*(rData.BDF0*r_v(5,0) + rData.BDF1*r_vn(5,0) + rData.BDF2*r_vnn(5,0)) + N_v[6]*(rData.BDF0*r_v(6,0) + rData.BDF1*r_vn(6,0) + rData.BDF2*r_vnn(6,0)) + N_v[7]*(rData.BDF0*r_v(7,0) + rData.BDF1*r_vn(7,0) + rData.BDF2*r_vnn(7,0)) + N_v[8]*(rData.BDF0*r_v(8,0) + rData.BDF1*r_vn(8,0) + rData.BDF2*r_vnn(8,0)) + N_v[9]*(rData.BDF0*r_v(9,0) + rData.BDF1*r_vn(9,0) + rData.BDF2*r_vnn(9,0)));
    const double crRHS3 = DN_v(0,0)*r_v(0,0) + DN_v(1,0)*r_v(1,0) + DN_v(2,0)*r_v(2,0) + DN_v(3,0)*r_v(3,0) + DN_v(4,0)*r_v(4,0) + DN_v(5,0)*r_v(5,0) + DN_v(6,0)*r_v(6,0) + DN_v(7,0)*r_v(7,0) + DN_v(8,0)*r_v(8,0) + DN_v(9,0)*r_v(9,0);
    const double crRHS4 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0) + N_v[6]*vconv(6,0) + N_v[7]*vconv(7,0) + N_v[8]*vconv(8,0) + N_v[9]*vconv(9,0);
    const double crRHS5 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1) + N_v[6]*vconv(6,1) + N_v[7]*vconv(7,1) + N_v[8]*vconv(8,1) + N_v[9]*vconv(9,1);
    const double crRHS6 = N_v[0]*vconv(0,2) + N_v[1]*vconv(1,2) + N_v[2]*vconv(2,2) + N_v[3]*vconv(3,2) + N_v[4]*vconv(4,2) + N_v[5]*vconv(5,2) + N_v[6]*vconv(6,2) + N_v[7]*vconv(7,2) + N_v[8]*vconv(8,2) + N_v[9]*vconv(9,2);
    const double crRHS7 = rho*(crRHS3*crRHS4 + crRHS5*(DN_v(0,1)*r_v(0,0) + DN_v(1,1)*r_v(1,0) + DN_v(2,1)*r_v(2,0) + DN_v(3,1)*r_v(3,0) + DN_v(4,1)*r_v(4,0) + DN_v(5,1)*r_v(5,0) + DN_v(6,1)*r_v(6,0) + DN_v(7,1)*r_v(7,0) + DN_v(8,1)*r_v(8,0) + DN_v(9,1)*r_v(9,0)) + crRHS6*(DN_v(0,2)*r_v(0,0) + DN_v(1,2)*r_v(1,0) + DN_v(2,2)*r_v(2,0) + DN_v(3,2)*r_v(3,0) + DN_v(4,2)*r_v(4,0) + DN_v(5,2)*r_v(5,0) + DN_v(6,2)*r_v(6,0) + DN_v(7,2)*r_v(7,0) + DN_v(8,2)*r_v(8,0) + DN_v(9,2)*r_v(9,0)));
    const double crRHS8 = DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1) + DN_v(6,1)*r_v(6,1) + DN_v(7,1)*r_v(7,1) + DN_v(8,1)*r_v(8,1) + DN_v(9,1)*r_v(9,1);
    const double crRHS9 = DN_v(0,2)*r_v(0,2) + DN_v(1,2)*r_v(1,2) + DN_v(2,2)*r_v(2,2) + DN_v(3,2)*r_v(3,2) + DN_v(4,2)*r_v(4,2) + DN_v(5,2)*r_v(5,2) + DN_v(6,2)*r_v(6,2) + DN_v(7,2)*r_v(7,2) + DN_v(8,2)*r_v(8,2) + DN_v(9,2)*r_v(9,2);
    const double crRHS10 = crRHS3 + crRHS8 + crRHS9;
    const double crRHS11 = rho*stab_c2*std::sqrt(std::pow(crRHS4, 2) + std::pow(crRHS5, 2) + std::pow(crRHS6, 2));
    const double crRHS12 = crRHS10*(crRHS11*h/stab_c1 + mu);
    const double crRHS13 = 1.0*C(0,0);
    const double crRHS14 = DDN_v[0](0,0)*r_v(0,0);
    const double crRHS15 = DDN_v[1](0,0)*r_v(1,0);
    const double crRHS16 = DDN_v[2](0,0)*r_v(2,0);
    const double crRHS17 = DDN_v[3](0,0)*r_v(3,0);
    const double crRHS18 = DDN_v[4](0,0)*r_v(4,0);
    const double crRHS19 = DDN_v[5](0,0)*r_v(5,0);
    const double crRHS20 = DDN_v[6](0,0)*r_v(6,0);
    const double crRHS21 = DDN_v[7](0,0)*r_v(7,0);
    const double crRHS22 = DDN_v[8](0,0)*r_v(8,0);
    const double crRHS23 = DDN_v[9](0,0)*r_v(9,0);
    const double crRHS24 = 1.0*C(0,1);
    const double crRHS25 = DDN_v[0](0,1)*r_v(0,1);
    const double crRHS26 = DDN_v[1](0,1)*r_v(1,1);
    const double crRHS27 = DDN_v[2](0,1)*r_v(2,1);
    const double crRHS28 = DDN_v[3](0,1)*r_v(3,1);
    const double crRHS29 = DDN_v[4](0,1)*r_v(4,1);
    const double crRHS30 = DDN_v[5](0,1)*r_v(5,1);
    const double crRHS31 = DDN_v[6](0,1)*r_v(6,1);
    const double crRHS32 = DDN_v[7](0,1)*r_v(7,1);
    const double crRHS33 = DDN_v[8](0,1)*r_v(8,1);
    const double crRHS34 = DDN_v[9](0,1)*r_v(9,1);
    const double crRHS35 = 1.0*C(0,2);
    const double crRHS36 = DDN_v[0](0,2)*r_v(0,2);
    const double crRHS37 = DDN_v[1](0,2)*r_v(1,2);
    const double crRHS38 = DDN_v[2](0,2)*r_v(2,2);
    const double crRHS39 = DDN_v[3](0,2)*r_v(3,2);
    const double crRHS40 = DDN_v[4](0,2)*r_v(4,2);
    const double crRHS41 = DDN_v[5](0,2)*r_v(5,2);
    const double crRHS42 = DDN_v[6](0,2)*r_v(6,2);
    const double crRHS43 = DDN_v[7](0,2)*r_v(7,2);
    const double crRHS44 = DDN_v[8](0,2)*r_v(8,2);
    const double crRHS45 = DDN_v[9](0,2)*r_v(9,2);
    const double crRHS46 = 1.0*C(0,3);
    const double crRHS47 = 1.0*C(0,5);
    const double crRHS48 = 1.0*C(1,3);
    const double crRHS49 = 1.0*C(1,5);
    const double crRHS50 = 1.0*C(2,3);
    const double crRHS51 = 1.0*C(2,5);
    const double crRHS52 = DDN_v[0](0,0)*r_v(0,1) + DDN_v[0](0,1)*r_v(0,0);
    const double crRHS53 = DDN_v[1](0,0)*r_v(1,1) + DDN_v[1](0,1)*r_v(1,0);
    const double crRHS54 = DDN_v[2](0,0)*r_v(2,1) + DDN_v[2](0,1)*r_v(2,0);
    const double crRHS55 = DDN_v[3](0,0)*r_v(3,1) + DDN_v[3](0,1)*r_v(3,0);
    const double crRHS56 = DDN_v[4](0,0)*r_v(4,1) + DDN_v[4](0,1)*r_v(4,0);
    const double crRHS57 = DDN_v[5](0,0)*r_v(5,1) + DDN_v[5](0,1)*r_v(5,0);
    const double crRHS58 = DDN_v[6](0,0)*r_v(6,1) + DDN_v[6](0,1)*r_v(6,0);
    const double crRHS59 = DDN_v[7](0,0)*r_v(7,1) + DDN_v[7](0,1)*r_v(7,0);
    const double crRHS60 = DDN_v[8](0,0)*r_v(8,1) + DDN_v[8](0,1)*r_v(8,0);
    const double crRHS61 = DDN_v[9](0,0)*r_v(9,1) + DDN_v[9](0,1)*r_v(9,0);
    const double crRHS62 = DDN_v[0](0,1)*r_v(0,2) + DDN_v[0](0,2)*r_v(0,1);
    const double crRHS63 = 1.0*C(0,4);
    const double crRHS64 = DDN_v[1](0,1)*r_v(1,2) + DDN_v[1](0,2)*r_v(1,1);
    const double crRHS65 = DDN_v[2](0,1)*r_v(2,2) + DDN_v[2](0,2)*r_v(2,1);
    const double crRHS66 = DDN_v[3](0,1)*r_v(3,2) + DDN_v[3](0,2)*r_v(3,1);
    const double crRHS67 = DDN_v[4](0,1)*r_v(4,2) + DDN_v[4](0,2)*r_v(4,1);
    const double crRHS68 = DDN_v[5](0,1)*r_v(5,2) + DDN_v[5](0,2)*r_v(5,1);
    const double crRHS69 = DDN_v[6](0,1)*r_v(6,2) + DDN_v[6](0,2)*r_v(6,1);
    const double crRHS70 = DDN_v[7](0,1)*r_v(7,2) + DDN_v[7](0,2)*r_v(7,1);
    const double crRHS71 = DDN_v[8](0,1)*r_v(8,2) + DDN_v[8](0,2)*r_v(8,1);
    const double crRHS72 = DDN_v[9](0,1)*r_v(9,2) + DDN_v[9](0,2)*r_v(9,1);
    const double crRHS73 = DDN_v[0](0,0)*r_v(0,2) + DDN_v[0](0,2)*r_v(0,0);
    const double crRHS74 = DDN_v[1](0,0)*r_v(1,2) + DDN_v[1](0,2)*r_v(1,0);
    const double crRHS75 = DDN_v[2](0,0)*r_v(2,2) + DDN_v[2](0,2)*r_v(2,0);
    const double crRHS76 = DDN_v[3](0,0)*r_v(3,2) + DDN_v[3](0,2)*r_v(3,0);
    const double crRHS77 = DDN_v[4](0,0)*r_v(4,2) + DDN_v[4](0,2)*r_v(4,0);
    const double crRHS78 = DDN_v[5](0,0)*r_v(5,2) + DDN_v[5](0,2)*r_v(5,0);
    const double crRHS79 = DDN_v[6](0,0)*r_v(6,2) + DDN_v[6](0,2)*r_v(6,0);
    const double crRHS80 = DDN_v[7](0,0)*r_v(7,2) + DDN_v[7](0,2)*r_v(7,0);
    const double crRHS81 = DDN_v[8](0,0)*r_v(8,2) + DDN_v[8](0,2)*r_v(8,0);
    const double crRHS82 = DDN_v[9](0,0)*r_v(9,2) + DDN_v[9](0,2)*r_v(9,0);
    const double crRHS83 = 1.0*C(3,3);
    const double crRHS84 = 1.0*C(3,4);
    const double crRHS85 = 1.0*C(3,5);
    const double crRHS86 = 1.0*C(4,5);
    const double crRHS87 = 1.0*C(5,5);
    const double crRHS88 = 1.0/(crRHS11/h + dyn_tau*rho/rData.DeltaTime + mu*stab_c1/std::pow(h, 2));
    const double crRHS89 = crRHS88*(-DN_p(0,0)*r_p[0] - DN_p(1,0)*r_p[1] - DN_p(2,0)*r_p[2] - DN_p(3,0)*r_p[3] + crRHS1 + crRHS13*crRHS14 + crRHS13*crRHS15 + crRHS13*crRHS16 + crRHS13*crRHS17 + crRHS13*crRHS18 + crRHS13*crRHS19 + crRHS13*crRHS20 + crRHS13*crRHS21 + crRHS13*crRHS22 + crRHS13*crRHS23 + crRHS14*crRHS46 + crRHS14*crRHS47 + crRHS15*crRHS46 + crRHS15*crRHS47 + crRHS16*crRHS46 + crRHS16*crRHS47 + crRHS17*crRHS46 + crRHS17*crRHS47 + crRHS18*crRHS46 + crRHS18*crRHS47 + crRHS19*crRHS46 + crRHS19*crRHS47 - crRHS2 + crRHS20*crRHS46 + crRHS20*crRHS47 + crRHS21*crRHS46 + crRHS21*crRHS47 + crRHS22*crRHS46 + crRHS22*crRHS47 + crRHS23*crRHS46 + crRHS23*crRHS47 + crRHS24*crRHS25 + crRHS24*crRHS26 + crRHS24*crRHS27 + crRHS24*crRHS28 + crRHS24*crRHS29 + crRHS24*crRHS30 + crRHS24*crRHS31 + crRHS24*crRHS32 + crRHS24*crRHS33 + crRHS24*crRHS34 + crRHS25*crRHS48 + crRHS25*crRHS49 + crRHS26*crRHS48 + crRHS26*crRHS49 + crRHS27*crRHS48 + crRHS27*crRHS49 + crRHS28*crRHS48 + crRHS28*crRHS49 + crRHS29*crRHS48 + crRHS29*crRHS49 + crRHS30*crRHS48 + crRHS30*crRHS49 + crRHS31*crRHS48 + crRHS31*crRHS49 + crRHS32*crRHS48 + crRHS32*crRHS49 + crRHS33*crRHS48 + crRHS33*crRHS49 + crRHS34*crRHS48 + crRHS34*crRHS49 + crRHS35*crRHS36 + crRHS35*crRHS37 + crRHS35*crRHS38 + crRHS35*crRHS39 + crRHS35*crRHS40 + crRHS35*crRHS41 + crRHS35*crRHS42 + crRHS35*crRHS43 + crRHS35*crRHS44 + crRHS35*crRHS45 + crRHS36*crRHS50 + crRHS36*crRHS51 + crRHS37*crRHS50 + crRHS37*crRHS51 + crRHS38*crRHS50 + crRHS38*crRHS51 + crRHS39*crRHS50 + crRHS39*crRHS51 + crRHS40*crRHS50 + crRHS40*crRHS51 + crRHS41*crRHS50 + crRHS41*crRHS51 + crRHS42*crRHS50 + crRHS42*crRHS51 + crRHS43*crRHS50 + crRHS43*crRHS51 + crRHS44*crRHS50 + crRHS44*crRHS51 + crRHS45*crRHS50 + crRHS45*crRHS51 + crRHS46*crRHS52 + crRHS46*crRHS53 + crRHS46*crRHS54 + crRHS46*crRHS55 + crRHS46*crRHS56 + crRHS46*crRHS57 + crRHS46*crRHS58 + crRHS46*crRHS59 + crRHS46*crRHS60 + crRHS46*crRHS61 + crRHS47*crRHS73 + crRHS47*crRHS74 + crRHS47*crRHS75 + crRHS47*crRHS76 + crRHS47*crRHS77 + crRHS47*crRHS78 + crRHS47*crRHS79 + crRHS47*crRHS80 + crRHS47*crRHS81 + crRHS47*crRHS82 + crRHS52*crRHS83 + crRHS52*crRHS85 + crRHS53*crRHS83 + crRHS53*crRHS85 + crRHS54*crRHS83 + crRHS54*crRHS85 + crRHS55*crRHS83 + crRHS55*crRHS85 + crRHS56*crRHS83 + crRHS56*crRHS85 + crRHS57*crRHS83 + crRHS57*crRHS85 + crRHS58*crRHS83 + crRHS58*crRHS85 + crRHS59*crRHS83 + crRHS59*crRHS85 + crRHS60*crRHS83 + crRHS60*crRHS85 + crRHS61*crRHS83 + crRHS61*crRHS85 + crRHS62*crRHS63 + crRHS62*crRHS84 + crRHS62*crRHS86 + crRHS63*crRHS64 + crRHS63*crRHS65 + crRHS63*crRHS66 + crRHS63*crRHS67 + crRHS63*crRHS68 + crRHS63*crRHS69 + crRHS63*crRHS70 + crRHS63*crRHS71 + crRHS63*crRHS72 + crRHS64*crRHS84 + crRHS64*crRHS86 + crRHS65*crRHS84 + crRHS65*crRHS86 + crRHS66*crRHS84 + crRHS66*crRHS86 + crRHS67*crRHS84 + crRHS67*crRHS86 + crRHS68*crRHS84 + crRHS68*crRHS86 + crRHS69*crRHS84 + crRHS69*crRHS86 - crRHS7 + crRHS70*crRHS84 + crRHS70*crRHS86 + crRHS71*crRHS84 + crRHS71*crRHS86 + crRHS72*crRHS84 + crRHS72*crRHS86 + crRHS73*crRHS85 + crRHS73*crRHS87 + crRHS74*crRHS85 + crRHS74*crRHS87 + crRHS75*crRHS85 + crRHS75*crRHS87 + crRHS76*crRHS85 + crRHS76*crRHS87 + crRHS77*crRHS85 + crRHS77*crRHS87 + crRHS78*crRHS85 + crRHS78*crRHS87 + crRHS79*crRHS85 + crRHS79*crRHS87 + crRHS80*crRHS85 + crRHS80*crRHS87 + crRHS81*crRHS85 + crRHS81*crRHS87 + crRHS82*crRHS85 + crRHS82*crRHS87);
    const double crRHS90 = rho*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(0,2)*vconv(0,2) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(1,2)*vconv(1,2) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(2,2)*vconv(2,2) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(3,2)*vconv(3,2) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(4,2)*vconv(4,2) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1) + DN_v(5,2)*vconv(5,2) + DN_v(6,0)*vconv(6,0) + DN_v(6,1)*vconv(6,1) + DN_v(6,2)*vconv(6,2) + DN_v(7,0)*vconv(7,0) + DN_v(7,1)*vconv(7,1) + DN_v(7,2)*vconv(7,2) + DN_v(8,0)*vconv(8,0) + DN_v(8,1)*vconv(8,1) + DN_v(8,2)*vconv(8,2) + DN_v(9,0)*vconv(9,0) + DN_v(9,1)*vconv(9,1) + DN_v(9,2)*vconv(9,2));
    const double crRHS91 = N_v[0]*crRHS90;
    const double crRHS92 = rho*(DN_v(0,0)*crRHS4 + DN_v(0,1)*crRHS5 + DN_v(0,2)*crRHS6);
    const double crRHS93 = rho*(N_v[0]*r_f(0,1) + N_v[1]*r_f(1,1) + N_v[2]*r_f(2,1) + N_v[3]*r_f(3,1) + N_v[4]*r_f(4,1) + N_v[5]*r_f(5,1) + N_v[6]*r_f(6,1) + N_v[7]*r_f(7,1) + N_v[8]*r_f(8,1) + N_v[9]*r_f(9,1));
    const double crRHS94 = rho*(N_v[0]*(rData.BDF0*r_v(0,1) + rData.BDF1*r_vn(0,1) + rData.BDF2*r_vnn(0,1)) + N_v[1]*(rData.BDF0*r_v(1,1) + rData.BDF1*r_vn(1,1) + rData.BDF2*r_vnn(1,1)) + N_v[2]*(rData.BDF0*r_v(2,1) + rData.BDF1*r_vn(2,1) + rData.BDF2*r_vnn(2,1)) + N_v[3]*(rData.BDF0*r_v(3,1) + rData.BDF1*r_vn(3,1) + rData.BDF2*r_vnn(3,1)) + N_v[4]*(rData.BDF0*r_v(4,1) + rData.BDF1*r_vn(4,1) + rData.BDF2*r_vnn(4,1)) + N_v[5]*(rData.BDF0*r_v(5,1) + rData.BDF1*r_vn(5,1) + rData.BDF2*r_vnn(5,1)) + N_v[6]*(rData.BDF0*r_v(6,1) + rData.BDF1*r_vn(6,1) + rData.BDF2*r_vnn(6,1)) + N_v[7]*(rData.BDF0*r_v(7,1) + rData.BDF1*r_vn(7,1) + rData.BDF2*r_vnn(7,1)) + N_v[8]*(rData.BDF0*r_v(8,1) + rData.BDF1*r_vn(8,1) + rData.BDF2*r_vnn(8,1)) + N_v[9]*(rData.BDF0*r_v(9,1) + rData.BDF1*r_vn(9,1) + rData.BDF2*r_vnn(9,1)));
    const double crRHS95 = rho*(crRHS4*(DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1) + DN_v(6,0)*r_v(6,1) + DN_v(7,0)*r_v(7,1) + DN_v(8,0)*r_v(8,1) + DN_v(9,0)*r_v(9,1)) + crRHS5*crRHS8 + crRHS6*(DN_v(0,2)*r_v(0,1) + DN_v(1,2)*r_v(1,1) + DN_v(2,2)*r_v(2,1) + DN_v(3,2)*r_v(3,1) + DN_v(4,2)*r_v(4,1) + DN_v(5,2)*r_v(5,1) + DN_v(6,2)*r_v(6,1) + DN_v(7,2)*r_v(7,1) + DN_v(8,2)*r_v(8,1) + DN_v(9,2)*r_v(9,1)));
    const double crRHS96 = DDN_v[0](1,0)*r_v(0,0);
    const double crRHS97 = DDN_v[1](1,0)*r_v(1,0);
    const double crRHS98 = DDN_v[2](1,0)*r_v(2,0);
    const double crRHS99 = DDN_v[3](1,0)*r_v(3,0);
    const double crRHS100 = DDN_v[4](1,0)*r_v(4,0);
    const double crRHS101 = DDN_v[5](1,0)*r_v(5,0);
    const double crRHS102 = DDN_v[6](1,0)*r_v(6,0);
    const double crRHS103 = DDN_v[7](1,0)*r_v(7,0);
    const double crRHS104 = DDN_v[8](1,0)*r_v(8,0);
    const double crRHS105 = DDN_v[9](1,0)*r_v(9,0);
    const double crRHS106 = 1.0*C(1,1);
    const double crRHS107 = DDN_v[0](1,1)*r_v(0,1);
    const double crRHS108 = DDN_v[1](1,1)*r_v(1,1);
    const double crRHS109 = DDN_v[2](1,1)*r_v(2,1);
    const double crRHS110 = DDN_v[3](1,1)*r_v(3,1);
    const double crRHS111 = DDN_v[4](1,1)*r_v(4,1);
    const double crRHS112 = DDN_v[5](1,1)*r_v(5,1);
    const double crRHS113 = DDN_v[6](1,1)*r_v(6,1);
    const double crRHS114 = DDN_v[7](1,1)*r_v(7,1);
    const double crRHS115 = DDN_v[8](1,1)*r_v(8,1);
    const double crRHS116 = DDN_v[9](1,1)*r_v(9,1);
    const double crRHS117 = 1.0*C(1,2);
    const double crRHS118 = DDN_v[0](1,2)*r_v(0,2);
    const double crRHS119 = DDN_v[1](1,2)*r_v(1,2);
    const double crRHS120 = DDN_v[2](1,2)*r_v(2,2);
    const double crRHS121 = DDN_v[3](1,2)*r_v(3,2);
    const double crRHS122 = DDN_v[4](1,2)*r_v(4,2);
    const double crRHS123 = DDN_v[5](1,2)*r_v(5,2);
    const double crRHS124 = DDN_v[6](1,2)*r_v(6,2);
    const double crRHS125 = DDN_v[7](1,2)*r_v(7,2);
    const double crRHS126 = DDN_v[8](1,2)*r_v(8,2);
    const double crRHS127 = DDN_v[9](1,2)*r_v(9,2);
    const double crRHS128 = 1.0*C(1,4);
    const double crRHS129 = 1.0*C(2,4);
    const double crRHS130 = DDN_v[0](1,0)*r_v(0,1) + DDN_v[0](1,1)*r_v(0,0);
    const double crRHS131 = DDN_v[1](1,0)*r_v(1,1) + DDN_v[1](1,1)*r_v(1,0);
    const double crRHS132 = DDN_v[2](1,0)*r_v(2,1) + DDN_v[2](1,1)*r_v(2,0);
    const double crRHS133 = DDN_v[3](1,0)*r_v(3,1) + DDN_v[3](1,1)*r_v(3,0);
    const double crRHS134 = DDN_v[4](1,0)*r_v(4,1) + DDN_v[4](1,1)*r_v(4,0);
    const double crRHS135 = DDN_v[5](1,0)*r_v(5,1) + DDN_v[5](1,1)*r_v(5,0);
    const double crRHS136 = DDN_v[6](1,0)*r_v(6,1) + DDN_v[6](1,1)*r_v(6,0);
    const double crRHS137 = DDN_v[7](1,0)*r_v(7,1) + DDN_v[7](1,1)*r_v(7,0);
    const double crRHS138 = DDN_v[8](1,0)*r_v(8,1) + DDN_v[8](1,1)*r_v(8,0);
    const double crRHS139 = DDN_v[9](1,0)*r_v(9,1) + DDN_v[9](1,1)*r_v(9,0);
    const double crRHS140 = DDN_v[0](1,1)*r_v(0,2) + DDN_v[0](1,2)*r_v(0,1);
    const double crRHS141 = DDN_v[1](1,1)*r_v(1,2) + DDN_v[1](1,2)*r_v(1,1);
    const double crRHS142 = DDN_v[2](1,1)*r_v(2,2) + DDN_v[2](1,2)*r_v(2,1);
    const double crRHS143 = DDN_v[3](1,1)*r_v(3,2) + DDN_v[3](1,2)*r_v(3,1);
    const double crRHS144 = DDN_v[4](1,1)*r_v(4,2) + DDN_v[4](1,2)*r_v(4,1);
    const double crRHS145 = DDN_v[5](1,1)*r_v(5,2) + DDN_v[5](1,2)*r_v(5,1);
    const double crRHS146 = DDN_v[6](1,1)*r_v(6,2) + DDN_v[6](1,2)*r_v(6,1);
    const double crRHS147 = DDN_v[7](1,1)*r_v(7,2) + DDN_v[7](1,2)*r_v(7,1);
    const double crRHS148 = DDN_v[8](1,1)*r_v(8,2) + DDN_v[8](1,2)*r_v(8,1);
    const double crRHS149 = DDN_v[9](1,1)*r_v(9,2) + DDN_v[9](1,2)*r_v(9,1);
    const double crRHS150 = DDN_v[0](1,0)*r_v(0,2) + DDN_v[0](1,2)*r_v(0,0);
    const double crRHS151 = DDN_v[1](1,0)*r_v(1,2) + DDN_v[1](1,2)*r_v(1,0);
    const double crRHS152 = DDN_v[2](1,0)*r_v(2,2) + DDN_v[2](1,2)*r_v(2,0);
    const double crRHS153 = DDN_v[3](1,0)*r_v(3,2) + DDN_v[3](1,2)*r_v(3,0);
    const double crRHS154 = DDN_v[4](1,0)*r_v(4,2) + DDN_v[4](1,2)*r_v(4,0);
    const double crRHS155 = DDN_v[5](1,0)*r_v(5,2) + DDN_v[5](1,2)*r_v(5,0);
    const double crRHS156 = DDN_v[6](1,0)*r_v(6,2) + DDN_v[6](1,2)*r_v(6,0);
    const double crRHS157 = DDN_v[7](1,0)*r_v(7,2) + DDN_v[7](1,2)*r_v(7,0);
    const double crRHS158 = DDN_v[8](1,0)*r_v(8,2) + DDN_v[8](1,2)*r_v(8,0);
    const double crRHS159 = DDN_v[9](1,0)*r_v(9,2) + DDN_v[9](1,2)*r_v(9,0);
    const double crRHS160 = 1.0*C(4,4);
    const double crRHS161 = crRHS88*(-DN_p(0,1)*r_p[0] - DN_p(1,1)*r_p[1] - DN_p(2,1)*r_p[2] - DN_p(3,1)*r_p[3] + crRHS100*crRHS24 + crRHS100*crRHS46 + crRHS100*crRHS63 + crRHS101*crRHS24 + crRHS101*crRHS46 + crRHS101*crRHS63 + crRHS102*crRHS24 + crRHS102*crRHS46 + crRHS102*crRHS63 + crRHS103*crRHS24 + crRHS103*crRHS46 + crRHS103*crRHS63 + crRHS104*crRHS24 + crRHS104*crRHS46 + crRHS104*crRHS63 + crRHS105*crRHS24 + crRHS105*crRHS46 + crRHS105*crRHS63 + crRHS106*crRHS107 + crRHS106*crRHS108 + crRHS106*crRHS109 + crRHS106*crRHS110 + crRHS106*crRHS111 + crRHS106*crRHS112 + crRHS106*crRHS113 + crRHS106*crRHS114 + crRHS106*crRHS115 + crRHS106*crRHS116 + crRHS107*crRHS128 + crRHS107*crRHS48 + crRHS108*crRHS128 + crRHS108*crRHS48 + crRHS109*crRHS128 + crRHS109*crRHS48 + crRHS110*crRHS128 + crRHS110*crRHS48 + crRHS111*crRHS128 + crRHS111*crRHS48 + crRHS112*crRHS128 + crRHS112*crRHS48 + crRHS113*crRHS128 + crRHS113*crRHS48 + crRHS114*crRHS128 + crRHS114*crRHS48 + crRHS115*crRHS128 + crRHS115*crRHS48 + crRHS116*crRHS128 + crRHS116*crRHS48 + crRHS117*crRHS118 + crRHS117*crRHS119 + crRHS117*crRHS120 + crRHS117*crRHS121 + crRHS117*crRHS122 + crRHS117*crRHS123 + crRHS117*crRHS124 + crRHS117*crRHS125 + crRHS117*crRHS126 + crRHS117*crRHS127 + crRHS118*crRHS129 + crRHS118*crRHS50 + crRHS119*crRHS129 + crRHS119*crRHS50 + crRHS120*crRHS129 + crRHS120*crRHS50 + crRHS121*crRHS129 + crRHS121*crRHS50 + crRHS122*crRHS129 + crRHS122*crRHS50 + crRHS123*crRHS129 + crRHS123*crRHS50 + crRHS124*crRHS129 + crRHS124*crRHS50 + crRHS125*crRHS129 + crRHS125*crRHS50 + crRHS126*crRHS129 + crRHS126*crRHS50 + crRHS127*crRHS129 + crRHS127*crRHS50 + crRHS128*crRHS140 + crRHS128*crRHS141 + crRHS128*crRHS142 + crRHS128*crRHS143 + crRHS128*crRHS144 + crRHS128*crRHS145 + crRHS128*crRHS146 + crRHS128*crRHS147 + crRHS128*crRHS148 + crRHS128*crRHS149 + crRHS130*crRHS48 + crRHS130*crRHS83 + crRHS130*crRHS84 + crRHS131*crRHS48 + crRHS131*crRHS83 + crRHS131*crRHS84 + crRHS132*crRHS48 + crRHS132*crRHS83 + crRHS132*crRHS84 + crRHS133*crRHS48 + crRHS133*crRHS83 + crRHS133*crRHS84 + crRHS134*crRHS48 + crRHS134*crRHS83 + crRHS134*crRHS84 + crRHS135*crRHS48 + crRHS135*crRHS83 + crRHS135*crRHS84 + crRHS136*crRHS48 + crRHS136*crRHS83 + crRHS136*crRHS84 + crRHS137*crRHS48 + crRHS137*crRHS83 + crRHS137*crRHS84 + crRHS138*crRHS48 + crRHS138*crRHS83 + crRHS138*crRHS84 + crRHS139*crRHS48 + crRHS139*crRHS83 + crRHS139*crRHS84 + crRHS140*crRHS160 + crRHS140*crRHS84 + crRHS141*crRHS160 + crRHS141*crRHS84 + crRHS142*crRHS160 + crRHS142*crRHS84 + crRHS143*crRHS160 + crRHS143*crRHS84 + crRHS144*crRHS160 + crRHS144*crRHS84 + crRHS145*crRHS160 + crRHS145*crRHS84 + crRHS146*crRHS160 + crRHS146*crRHS84 + crRHS147*crRHS160 + crRHS147*crRHS84 + crRHS148*crRHS160 + crRHS148*crRHS84 + crRHS149*crRHS160 + crRHS149*crRHS84 + crRHS150*crRHS49 + crRHS150*crRHS85 + crRHS150*crRHS86 + crRHS151*crRHS49 + crRHS151*crRHS85 + crRHS151*crRHS86 + crRHS152*crRHS49 + crRHS152*crRHS85 + crRHS152*crRHS86 + crRHS153*crRHS49 + crRHS153*crRHS85 + crRHS153*crRHS86 + crRHS154*crRHS49 + crRHS154*crRHS85 + crRHS154*crRHS86 + crRHS155*crRHS49 + crRHS155*crRHS85 + crRHS155*crRHS86 + crRHS156*crRHS49 + crRHS156*crRHS85 + crRHS156*crRHS86 + crRHS157*crRHS49 + crRHS157*crRHS85 + crRHS157*crRHS86 + crRHS158*crRHS49 + crRHS158*crRHS85 + crRHS158*crRHS86 + crRHS159*crRHS49 + crRHS159*crRHS85 + crRHS159*crRHS86 + crRHS24*crRHS96 + crRHS24*crRHS97 + crRHS24*crRHS98 + crRHS24*crRHS99 + crRHS46*crRHS96 + crRHS46*crRHS97 + crRHS46*crRHS98 + crRHS46*crRHS99 + crRHS63*crRHS96 + crRHS63*crRHS97 + crRHS63*crRHS98 + crRHS63*crRHS99 + crRHS93 - crRHS94 - crRHS95);
    const double crRHS162 = rho*(N_v[0]*r_f(0,2) + N_v[1]*r_f(1,2) + N_v[2]*r_f(2,2) + N_v[3]*r_f(3,2) + N_v[4]*r_f(4,2) + N_v[5]*r_f(5,2) + N_v[6]*r_f(6,2) + N_v[7]*r_f(7,2) + N_v[8]*r_f(8,2) + N_v[9]*r_f(9,2));
    const double crRHS163 = rho*(N_v[0]*(rData.BDF0*r_v(0,2) + rData.BDF1*r_vn(0,2) + rData.BDF2*r_vnn(0,2)) + N_v[1]*(rData.BDF0*r_v(1,2) + rData.BDF1*r_vn(1,2) + rData.BDF2*r_vnn(1,2)) + N_v[2]*(rData.BDF0*r_v(2,2) + rData.BDF1*r_vn(2,2) + rData.BDF2*r_vnn(2,2)) + N_v[3]*(rData.BDF0*r_v(3,2) + rData.BDF1*r_vn(3,2) + rData.BDF2*r_vnn(3,2)) + N_v[4]*(rData.BDF0*r_v(4,2) + rData.BDF1*r_vn(4,2) + rData.BDF2*r_vnn(4,2)) + N_v[5]*(rData.BDF0*r_v(5,2) + rData.BDF1*r_vn(5,2) + rData.BDF2*r_vnn(5,2)) + N_v[6]*(rData.BDF0*r_v(6,2) + rData.BDF1*r_vn(6,2) + rData.BDF2*r_vnn(6,2)) + N_v[7]*(rData.BDF0*r_v(7,2) + rData.BDF1*r_vn(7,2) + rData.BDF2*r_vnn(7,2)) + N_v[8]*(rData.BDF0*r_v(8,2) + rData.BDF1*r_vn(8,2) + rData.BDF2*r_vnn(8,2)) + N_v[9]*(rData.BDF0*r_v(9,2) + rData.BDF1*r_vn(9,2) + rData.BDF2*r_vnn(9,2)));
    const double crRHS164 = rho*(crRHS4*(DN_v(0,0)*r_v(0,2) + DN_v(1,0)*r_v(1,2) + DN_v(2,0)*r_v(2,2) + DN_v(3,0)*r_v(3,2) + DN_v(4,0)*r_v(4,2) + DN_v(5,0)*r_v(5,2) + DN_v(6,0)*r_v(6,2) + DN_v(7,0)*r_v(7,2) + DN_v(8,0)*r_v(8,2) + DN_v(9,0)*r_v(9,2)) + crRHS5*(DN_v(0,1)*r_v(0,2) + DN_v(1,1)*r_v(1,2) + DN_v(2,1)*r_v(2,2) + DN_v(3,1)*r_v(3,2) + DN_v(4,1)*r_v(4,2) + DN_v(5,1)*r_v(5,2) + DN_v(6,1)*r_v(6,2) + DN_v(7,1)*r_v(7,2) + DN_v(8,1)*r_v(8,2) + DN_v(9,1)*r_v(9,2)) + crRHS6*crRHS9);
    const double crRHS165 = DDN_v[0](2,0)*r_v(0,0);
    const double crRHS166 = DDN_v[1](2,0)*r_v(1,0);
    const double crRHS167 = DDN_v[2](2,0)*r_v(2,0);
    const double crRHS168 = DDN_v[3](2,0)*r_v(3,0);
    const double crRHS169 = DDN_v[4](2,0)*r_v(4,0);
    const double crRHS170 = DDN_v[5](2,0)*r_v(5,0);
    const double crRHS171 = DDN_v[6](2,0)*r_v(6,0);
    const double crRHS172 = DDN_v[7](2,0)*r_v(7,0);
    const double crRHS173 = DDN_v[8](2,0)*r_v(8,0);
    const double crRHS174 = DDN_v[9](2,0)*r_v(9,0);
    const double crRHS175 = DDN_v[0](2,1)*r_v(0,1);
    const double crRHS176 = DDN_v[1](2,1)*r_v(1,1);
    const double crRHS177 = DDN_v[2](2,1)*r_v(2,1);
    const double crRHS178 = DDN_v[3](2,1)*r_v(3,1);
    const double crRHS179 = DDN_v[4](2,1)*r_v(4,1);
    const double crRHS180 = DDN_v[5](2,1)*r_v(5,1);
    const double crRHS181 = DDN_v[6](2,1)*r_v(6,1);
    const double crRHS182 = DDN_v[7](2,1)*r_v(7,1);
    const double crRHS183 = DDN_v[8](2,1)*r_v(8,1);
    const double crRHS184 = DDN_v[9](2,1)*r_v(9,1);
    const double crRHS185 = 1.0*C(2,2);
    const double crRHS186 = DDN_v[0](2,2)*r_v(0,2);
    const double crRHS187 = DDN_v[1](2,2)*r_v(1,2);
    const double crRHS188 = DDN_v[2](2,2)*r_v(2,2);
    const double crRHS189 = DDN_v[3](2,2)*r_v(3,2);
    const double crRHS190 = DDN_v[4](2,2)*r_v(4,2);
    const double crRHS191 = DDN_v[5](2,2)*r_v(5,2);
    const double crRHS192 = DDN_v[6](2,2)*r_v(6,2);
    const double crRHS193 = DDN_v[7](2,2)*r_v(7,2);
    const double crRHS194 = DDN_v[8](2,2)*r_v(8,2);
    const double crRHS195 = DDN_v[9](2,2)*r_v(9,2);
    const double crRHS196 = DDN_v[0](2,0)*r_v(0,1) + DDN_v[0](2,1)*r_v(0,0);
    const double crRHS197 = DDN_v[1](2,0)*r_v(1,1) + DDN_v[1](2,1)*r_v(1,0);
    const double crRHS198 = DDN_v[2](2,0)*r_v(2,1) + DDN_v[2](2,1)*r_v(2,0);
    const double crRHS199 = DDN_v[3](2,0)*r_v(3,1) + DDN_v[3](2,1)*r_v(3,0);
    const double crRHS200 = DDN_v[4](2,0)*r_v(4,1) + DDN_v[4](2,1)*r_v(4,0);
    const double crRHS201 = DDN_v[5](2,0)*r_v(5,1) + DDN_v[5](2,1)*r_v(5,0);
    const double crRHS202 = DDN_v[6](2,0)*r_v(6,1) + DDN_v[6](2,1)*r_v(6,0);
    const double crRHS203 = DDN_v[7](2,0)*r_v(7,1) + DDN_v[7](2,1)*r_v(7,0);
    const double crRHS204 = DDN_v[8](2,0)*r_v(8,1) + DDN_v[8](2,1)*r_v(8,0);
    const double crRHS205 = DDN_v[9](2,0)*r_v(9,1) + DDN_v[9](2,1)*r_v(9,0);
    const double crRHS206 = DDN_v[0](2,1)*r_v(0,2) + DDN_v[0](2,2)*r_v(0,1);
    const double crRHS207 = DDN_v[1](2,1)*r_v(1,2) + DDN_v[1](2,2)*r_v(1,1);
    const double crRHS208 = DDN_v[2](2,1)*r_v(2,2) + DDN_v[2](2,2)*r_v(2,1);
    const double crRHS209 = DDN_v[3](2,1)*r_v(3,2) + DDN_v[3](2,2)*r_v(3,1);
    const double crRHS210 = DDN_v[4](2,1)*r_v(4,2) + DDN_v[4](2,2)*r_v(4,1);
    const double crRHS211 = DDN_v[5](2,1)*r_v(5,2) + DDN_v[5](2,2)*r_v(5,1);
    const double crRHS212 = DDN_v[6](2,1)*r_v(6,2) + DDN_v[6](2,2)*r_v(6,1);
    const double crRHS213 = DDN_v[7](2,1)*r_v(7,2) + DDN_v[7](2,2)*r_v(7,1);
    const double crRHS214 = DDN_v[8](2,1)*r_v(8,2) + DDN_v[8](2,2)*r_v(8,1);
    const double crRHS215 = DDN_v[9](2,1)*r_v(9,2) + DDN_v[9](2,2)*r_v(9,1);
    const double crRHS216 = DDN_v[0](2,0)*r_v(0,2) + DDN_v[0](2,2)*r_v(0,0);
    const double crRHS217 = DDN_v[1](2,0)*r_v(1,2) + DDN_v[1](2,2)*r_v(1,0);
    const double crRHS218 = DDN_v[2](2,0)*r_v(2,2) + DDN_v[2](2,2)*r_v(2,0);
    const double crRHS219 = DDN_v[3](2,0)*r_v(3,2) + DDN_v[3](2,2)*r_v(3,0);
    const double crRHS220 = DDN_v[4](2,0)*r_v(4,2) + DDN_v[4](2,2)*r_v(4,0);
    const double crRHS221 = DDN_v[5](2,0)*r_v(5,2) + DDN_v[5](2,2)*r_v(5,0);
    const double crRHS222 = DDN_v[6](2,0)*r_v(6,2) + DDN_v[6](2,2)*r_v(6,0);
    const double crRHS223 = DDN_v[7](2,0)*r_v(7,2) + DDN_v[7](2,2)*r_v(7,0);
    const double crRHS224 = DDN_v[8](2,0)*r_v(8,2) + DDN_v[8](2,2)*r_v(8,0);
    const double crRHS225 = DDN_v[9](2,0)*r_v(9,2) + DDN_v[9](2,2)*r_v(9,0);
    const double crRHS226 = crRHS88*(-DN_p(0,2)*r_p[0] - DN_p(1,2)*r_p[1] - DN_p(2,2)*r_p[2] - DN_p(3,2)*r_p[3] + crRHS117*crRHS175 + crRHS117*crRHS176 + crRHS117*crRHS177 + crRHS117*crRHS178 + crRHS117*crRHS179 + crRHS117*crRHS180 + crRHS117*crRHS181 + crRHS117*crRHS182 + crRHS117*crRHS183 + crRHS117*crRHS184 + crRHS128*crRHS175 + crRHS128*crRHS176 + crRHS128*crRHS177 + crRHS128*crRHS178 + crRHS128*crRHS179 + crRHS128*crRHS180 + crRHS128*crRHS181 + crRHS128*crRHS182 + crRHS128*crRHS183 + crRHS128*crRHS184 + crRHS129*crRHS186 + crRHS129*crRHS187 + crRHS129*crRHS188 + crRHS129*crRHS189 + crRHS129*crRHS190 + crRHS129*crRHS191 + crRHS129*crRHS192 + crRHS129*crRHS193 + crRHS129*crRHS194 + crRHS129*crRHS195 + crRHS129*crRHS206 + crRHS129*crRHS207 + crRHS129*crRHS208 + crRHS129*crRHS209 + crRHS129*crRHS210 + crRHS129*crRHS211 + crRHS129*crRHS212 + crRHS129*crRHS213 + crRHS129*crRHS214 + crRHS129*crRHS215 + crRHS160*crRHS206 + crRHS160*crRHS207 + crRHS160*crRHS208 + crRHS160*crRHS209 + crRHS160*crRHS210 + crRHS160*crRHS211 + crRHS160*crRHS212 + crRHS160*crRHS213 + crRHS160*crRHS214 + crRHS160*crRHS215 + crRHS162 - crRHS163 - crRHS164 + crRHS165*crRHS35 + crRHS165*crRHS47 + crRHS165*crRHS63 + crRHS166*crRHS35 + crRHS166*crRHS47 + crRHS166*crRHS63 + crRHS167*crRHS35 + crRHS167*crRHS47 + crRHS167*crRHS63 + crRHS168*crRHS35 + crRHS168*crRHS47 + crRHS168*crRHS63 + crRHS169*crRHS35 + crRHS169*crRHS47 + crRHS169*crRHS63 + crRHS170*crRHS35 + crRHS170*crRHS47 + crRHS170*crRHS63 + crRHS171*crRHS35 + crRHS171*crRHS47 + crRHS171*crRHS63 + crRHS172*crRHS35 + crRHS172*crRHS47 + crRHS172*crRHS63 + crRHS173*crRHS35 + crRHS173*crRHS47 + crRHS173*crRHS63 + crRHS174*crRHS35 + crRHS174*crRHS47 + crRHS174*crRHS63 + crRHS175*crRHS49 + crRHS176*crRHS49 + crRHS177*crRHS49 + crRHS178*crRHS49 + crRHS179*crRHS49 + crRHS180*crRHS49 + crRHS181*crRHS49 + crRHS182*crRHS49 + crRHS183*crRHS49 + crRHS184*crRHS49 + crRHS185*crRHS186 + crRHS185*crRHS187 + crRHS185*crRHS188 + crRHS185*crRHS189 + crRHS185*crRHS190 + crRHS185*crRHS191 + crRHS185*crRHS192 + crRHS185*crRHS193 + crRHS185*crRHS194 + crRHS185*crRHS195 + crRHS186*crRHS51 + crRHS187*crRHS51 + crRHS188*crRHS51 + crRHS189*crRHS51 + crRHS190*crRHS51 + crRHS191*crRHS51 + crRHS192*crRHS51 + crRHS193*crRHS51 + crRHS194*crRHS51 + crRHS195*crRHS51 + crRHS196*crRHS50 + crRHS196*crRHS84 + crRHS196*crRHS85 + crRHS197*crRHS50 + crRHS197*crRHS84 + crRHS197*crRHS85 + crRHS198*crRHS50 + crRHS198*crRHS84 + crRHS198*crRHS85 + crRHS199*crRHS50 + crRHS199*crRHS84 + crRHS199*crRHS85 + crRHS200*crRHS50 + crRHS200*crRHS84 + crRHS200*crRHS85 + crRHS201*crRHS50 + crRHS201*crRHS84 + crRHS201*crRHS85 + crRHS202*crRHS50 + crRHS202*crRHS84 + crRHS202*crRHS85 + crRHS203*crRHS50 + crRHS203*crRHS84 + crRHS203*crRHS85 + crRHS204*crRHS50 + crRHS204*crRHS84 + crRHS204*crRHS85 + crRHS205*crRHS50 + crRHS205*crRHS84 + crRHS205*crRHS85 + crRHS206*crRHS86 + crRHS207*crRHS86 + crRHS208*crRHS86 + crRHS209*crRHS86 + crRHS210*crRHS86 + crRHS211*crRHS86 + crRHS212*crRHS86 + crRHS213*crRHS86 + crRHS214*crRHS86 + crRHS215*crRHS86 + crRHS216*crRHS51 + crRHS216*crRHS86 + crRHS216*crRHS87 + crRHS217*crRHS51 + crRHS217*crRHS86 + crRHS217*crRHS87 + crRHS218*crRHS51 + crRHS218*crRHS86 + crRHS218*crRHS87 + crRHS219*crRHS51 + crRHS219*crRHS86 + crRHS219*crRHS87 + crRHS220*crRHS51 + crRHS220*crRHS86 + crRHS220*crRHS87 + crRHS221*crRHS51 + crRHS221*crRHS86 + crRHS221*crRHS87 + crRHS222*crRHS51 + crRHS222*crRHS86 + crRHS222*crRHS87 + crRHS223*crRHS51 + crRHS223*crRHS86 + crRHS223*crRHS87 + crRHS224*crRHS51 + crRHS224*crRHS86 + crRHS224*crRHS87 + crRHS225*crRHS51 + crRHS225*crRHS86 + crRHS225*crRHS87);
    const double crRHS227 = N_v[1]*crRHS90;
    const double crRHS228 = rho*(DN_v(1,0)*crRHS4 + DN_v(1,1)*crRHS5 + DN_v(1,2)*crRHS6);
    const double crRHS229 = N_v[2]*crRHS90;
    const double crRHS230 = rho*(DN_v(2,0)*crRHS4 + DN_v(2,1)*crRHS5 + DN_v(2,2)*crRHS6);
    const double crRHS231 = N_v[3]*crRHS90;
    const double crRHS232 = rho*(DN_v(3,0)*crRHS4 + DN_v(3,1)*crRHS5 + DN_v(3,2)*crRHS6);
    const double crRHS233 = N_v[4]*crRHS90;
    const double crRHS234 = rho*(DN_v(4,0)*crRHS4 + DN_v(4,1)*crRHS5 + DN_v(4,2)*crRHS6);
    const double crRHS235 = N_v[5]*crRHS90;
    const double crRHS236 = rho*(DN_v(5,0)*crRHS4 + DN_v(5,1)*crRHS5 + DN_v(5,2)*crRHS6);
    const double crRHS237 = N_v[6]*crRHS90;
    const double crRHS238 = rho*(DN_v(6,0)*crRHS4 + DN_v(6,1)*crRHS5 + DN_v(6,2)*crRHS6);
    const double crRHS239 = N_v[7]*crRHS90;
    const double crRHS240 = rho*(DN_v(7,0)*crRHS4 + DN_v(7,1)*crRHS5 + DN_v(7,2)*crRHS6);
    const double crRHS241 = N_v[8]*crRHS90;
    const double crRHS242 = rho*(DN_v(8,0)*crRHS4 + DN_v(8,1)*crRHS5 + DN_v(8,2)*crRHS6);
    const double crRHS243 = N_v[9]*crRHS90;
    const double crRHS244 = rho*(DN_v(9,0)*crRHS4 + DN_v(9,1)*crRHS5 + DN_v(9,2)*crRHS6);
    rRHS[0]+=-gauss_weight*(-DN_v(0,0)*crRHS0 + DN_v(0,0)*crRHS12 + DN_v(0,0)*r_stress[0] + DN_v(0,1)*r_stress[3] + DN_v(0,2)*r_stress[5] - N_v[0]*crRHS1 + N_v[0]*crRHS2 + N_v[0]*crRHS7 - crRHS89*crRHS91 - crRHS89*crRHS92);
    rRHS[1]+=-gauss_weight*(DN_v(0,0)*r_stress[3] - DN_v(0,1)*crRHS0 + DN_v(0,1)*crRHS12 + DN_v(0,1)*r_stress[1] + DN_v(0,2)*r_stress[4] - N_v[0]*crRHS93 + N_v[0]*crRHS94 + N_v[0]*crRHS95 - crRHS161*crRHS91 - crRHS161*crRHS92);
    rRHS[2]+=-gauss_weight*(DN_v(0,0)*r_stress[5] + DN_v(0,1)*r_stress[4] - DN_v(0,2)*crRHS0 + DN_v(0,2)*crRHS12 + DN_v(0,2)*r_stress[2] - N_v[0]*crRHS162 + N_v[0]*crRHS163 + N_v[0]*crRHS164 - crRHS226*crRHS91 - crRHS226*crRHS92);
    rRHS[3]+=-gauss_weight*(-DN_v(1,0)*crRHS0 + DN_v(1,0)*crRHS12 + DN_v(1,0)*r_stress[0] + DN_v(1,1)*r_stress[3] + DN_v(1,2)*r_stress[5] - N_v[1]*crRHS1 + N_v[1]*crRHS2 + N_v[1]*crRHS7 - crRHS227*crRHS89 - crRHS228*crRHS89);
    rRHS[4]+=-gauss_weight*(DN_v(1,0)*r_stress[3] - DN_v(1,1)*crRHS0 + DN_v(1,1)*crRHS12 + DN_v(1,1)*r_stress[1] + DN_v(1,2)*r_stress[4] - N_v[1]*crRHS93 + N_v[1]*crRHS94 + N_v[1]*crRHS95 - crRHS161*crRHS227 - crRHS161*crRHS228);
    rRHS[5]+=-gauss_weight*(DN_v(1,0)*r_stress[5] + DN_v(1,1)*r_stress[4] - DN_v(1,2)*crRHS0 + DN_v(1,2)*crRHS12 + DN_v(1,2)*r_stress[2] - N_v[1]*crRHS162 + N_v[1]*crRHS163 + N_v[1]*crRHS164 - crRHS226*crRHS227 - crRHS226*crRHS228);
    rRHS[6]+=-gauss_weight*(-DN_v(2,0)*crRHS0 + DN_v(2,0)*crRHS12 + DN_v(2,0)*r_stress[0] + DN_v(2,1)*r_stress[3] + DN_v(2,2)*r_stress[5] - N_v[2]*crRHS1 + N_v[2]*crRHS2 + N_v[2]*crRHS7 - crRHS229*crRHS89 - crRHS230*crRHS89);
    rRHS[7]+=-gauss_weight*(DN_v(2,0)*r_stress[3] - DN_v(2,1)*crRHS0 + DN_v(2,1)*crRHS12 + DN_v(2,1)*r_stress[1] + DN_v(2,2)*r_stress[4] - N_v[2]*crRHS93 + N_v[2]*crRHS94 + N_v[2]*crRHS95 - crRHS161*crRHS229 - crRHS161*crRHS230);
    rRHS[8]+=-gauss_weight*(DN_v(2,0)*r_stress[5] + DN_v(2,1)*r_stress[4] - DN_v(2,2)*crRHS0 + DN_v(2,2)*crRHS12 + DN_v(2,2)*r_stress[2] - N_v[2]*crRHS162 + N_v[2]*crRHS163 + N_v[2]*crRHS164 - crRHS226*crRHS229 - crRHS226*crRHS230);
    rRHS[9]+=-gauss_weight*(-DN_v(3,0)*crRHS0 + DN_v(3,0)*crRHS12 + DN_v(3,0)*r_stress[0] + DN_v(3,1)*r_stress[3] + DN_v(3,2)*r_stress[5] - N_v[3]*crRHS1 + N_v[3]*crRHS2 + N_v[3]*crRHS7 - crRHS231*crRHS89 - crRHS232*crRHS89);
    rRHS[10]+=-gauss_weight*(DN_v(3,0)*r_stress[3] - DN_v(3,1)*crRHS0 + DN_v(3,1)*crRHS12 + DN_v(3,1)*r_stress[1] + DN_v(3,2)*r_stress[4] - N_v[3]*crRHS93 + N_v[3]*crRHS94 + N_v[3]*crRHS95 - crRHS161*crRHS231 - crRHS161*crRHS232);
    rRHS[11]+=-gauss_weight*(DN_v(3,0)*r_stress[5] + DN_v(3,1)*r_stress[4] - DN_v(3,2)*crRHS0 + DN_v(3,2)*crRHS12 + DN_v(3,2)*r_stress[2] - N_v[3]*crRHS162 + N_v[3]*crRHS163 + N_v[3]*crRHS164 - crRHS226*crRHS231 - crRHS226*crRHS232);
    rRHS[12]+=-gauss_weight*(-DN_v(4,0)*crRHS0 + DN_v(4,0)*crRHS12 + DN_v(4,0)*r_stress[0] + DN_v(4,1)*r_stress[3] + DN_v(4,2)*r_stress[5] - N_v[4]*crRHS1 + N_v[4]*crRHS2 + N_v[4]*crRHS7 - crRHS233*crRHS89 - crRHS234*crRHS89);
    rRHS[13]+=-gauss_weight*(DN_v(4,0)*r_stress[3] - DN_v(4,1)*crRHS0 + DN_v(4,1)*crRHS12 + DN_v(4,1)*r_stress[1] + DN_v(4,2)*r_stress[4] - N_v[4]*crRHS93 + N_v[4]*crRHS94 + N_v[4]*crRHS95 - crRHS161*crRHS233 - crRHS161*crRHS234);
    rRHS[14]+=-gauss_weight*(DN_v(4,0)*r_stress[5] + DN_v(4,1)*r_stress[4] - DN_v(4,2)*crRHS0 + DN_v(4,2)*crRHS12 + DN_v(4,2)*r_stress[2] - N_v[4]*crRHS162 + N_v[4]*crRHS163 + N_v[4]*crRHS164 - crRHS226*crRHS233 - crRHS226*crRHS234);
    rRHS[15]+=-gauss_weight*(-DN_v(5,0)*crRHS0 + DN_v(5,0)*crRHS12 + DN_v(5,0)*r_stress[0] + DN_v(5,1)*r_stress[3] + DN_v(5,2)*r_stress[5] - N_v[5]*crRHS1 + N_v[5]*crRHS2 + N_v[5]*crRHS7 - crRHS235*crRHS89 - crRHS236*crRHS89);
    rRHS[16]+=-gauss_weight*(DN_v(5,0)*r_stress[3] - DN_v(5,1)*crRHS0 + DN_v(5,1)*crRHS12 + DN_v(5,1)*r_stress[1] + DN_v(5,2)*r_stress[4] - N_v[5]*crRHS93 + N_v[5]*crRHS94 + N_v[5]*crRHS95 - crRHS161*crRHS235 - crRHS161*crRHS236);
    rRHS[17]+=-gauss_weight*(DN_v(5,0)*r_stress[5] + DN_v(5,1)*r_stress[4] - DN_v(5,2)*crRHS0 + DN_v(5,2)*crRHS12 + DN_v(5,2)*r_stress[2] - N_v[5]*crRHS162 + N_v[5]*crRHS163 + N_v[5]*crRHS164 - crRHS226*crRHS235 - crRHS226*crRHS236);
    rRHS[18]+=-gauss_weight*(-DN_v(6,0)*crRHS0 + DN_v(6,0)*crRHS12 + DN_v(6,0)*r_stress[0] + DN_v(6,1)*r_stress[3] + DN_v(6,2)*r_stress[5] - N_v[6]*crRHS1 + N_v[6]*crRHS2 + N_v[6]*crRHS7 - crRHS237*crRHS89 - crRHS238*crRHS89);
    rRHS[19]+=-gauss_weight*(DN_v(6,0)*r_stress[3] - DN_v(6,1)*crRHS0 + DN_v(6,1)*crRHS12 + DN_v(6,1)*r_stress[1] + DN_v(6,2)*r_stress[4] - N_v[6]*crRHS93 + N_v[6]*crRHS94 + N_v[6]*crRHS95 - crRHS161*crRHS237 - crRHS161*crRHS238);
    rRHS[20]+=-gauss_weight*(DN_v(6,0)*r_stress[5] + DN_v(6,1)*r_stress[4] - DN_v(6,2)*crRHS0 + DN_v(6,2)*crRHS12 + DN_v(6,2)*r_stress[2] - N_v[6]*crRHS162 + N_v[6]*crRHS163 + N_v[6]*crRHS164 - crRHS226*crRHS237 - crRHS226*crRHS238);
    rRHS[21]+=-gauss_weight*(-DN_v(7,0)*crRHS0 + DN_v(7,0)*crRHS12 + DN_v(7,0)*r_stress[0] + DN_v(7,1)*r_stress[3] + DN_v(7,2)*r_stress[5] - N_v[7]*crRHS1 + N_v[7]*crRHS2 + N_v[7]*crRHS7 - crRHS239*crRHS89 - crRHS240*crRHS89);
    rRHS[22]+=-gauss_weight*(DN_v(7,0)*r_stress[3] - DN_v(7,1)*crRHS0 + DN_v(7,1)*crRHS12 + DN_v(7,1)*r_stress[1] + DN_v(7,2)*r_stress[4] - N_v[7]*crRHS93 + N_v[7]*crRHS94 + N_v[7]*crRHS95 - crRHS161*crRHS239 - crRHS161*crRHS240);
    rRHS[23]+=-gauss_weight*(DN_v(7,0)*r_stress[5] + DN_v(7,1)*r_stress[4] - DN_v(7,2)*crRHS0 + DN_v(7,2)*crRHS12 + DN_v(7,2)*r_stress[2] - N_v[7]*crRHS162 + N_v[7]*crRHS163 + N_v[7]*crRHS164 - crRHS226*crRHS239 - crRHS226*crRHS240);
    rRHS[24]+=-gauss_weight*(-DN_v(8,0)*crRHS0 + DN_v(8,0)*crRHS12 + DN_v(8,0)*r_stress[0] + DN_v(8,1)*r_stress[3] + DN_v(8,2)*r_stress[5] - N_v[8]*crRHS1 + N_v[8]*crRHS2 + N_v[8]*crRHS7 - crRHS241*crRHS89 - crRHS242*crRHS89);
    rRHS[25]+=-gauss_weight*(DN_v(8,0)*r_stress[3] - DN_v(8,1)*crRHS0 + DN_v(8,1)*crRHS12 + DN_v(8,1)*r_stress[1] + DN_v(8,2)*r_stress[4] - N_v[8]*crRHS93 + N_v[8]*crRHS94 + N_v[8]*crRHS95 - crRHS161*crRHS241 - crRHS161*crRHS242);
    rRHS[26]+=-gauss_weight*(DN_v(8,0)*r_stress[5] + DN_v(8,1)*r_stress[4] - DN_v(8,2)*crRHS0 + DN_v(8,2)*crRHS12 + DN_v(8,2)*r_stress[2] - N_v[8]*crRHS162 + N_v[8]*crRHS163 + N_v[8]*crRHS164 - crRHS226*crRHS241 - crRHS226*crRHS242);
    rRHS[27]+=-gauss_weight*(-DN_v(9,0)*crRHS0 + DN_v(9,0)*crRHS12 + DN_v(9,0)*r_stress[0] + DN_v(9,1)*r_stress[3] + DN_v(9,2)*r_stress[5] - N_v[9]*crRHS1 + N_v[9]*crRHS2 + N_v[9]*crRHS7 - crRHS243*crRHS89 - crRHS244*crRHS89);
    rRHS[28]+=-gauss_weight*(DN_v(9,0)*r_stress[3] - DN_v(9,1)*crRHS0 + DN_v(9,1)*crRHS12 + DN_v(9,1)*r_stress[1] + DN_v(9,2)*r_stress[4] - N_v[9]*crRHS93 + N_v[9]*crRHS94 + N_v[9]*crRHS95 - crRHS161*crRHS243 - crRHS161*crRHS244);
    rRHS[29]+=-gauss_weight*(DN_v(9,0)*r_stress[5] + DN_v(9,1)*r_stress[4] - DN_v(9,2)*crRHS0 + DN_v(9,2)*crRHS12 + DN_v(9,2)*r_stress[2] - N_v[9]*crRHS162 + N_v[9]*crRHS163 + N_v[9]*crRHS164 - crRHS226*crRHS243 - crRHS226*crRHS244);
    rRHS[30]+=gauss_weight*(DN_p(0,0)*crRHS89 + DN_p(0,1)*crRHS161 + DN_p(0,2)*crRHS226 - N_p[0]*crRHS10);
    rRHS[31]+=gauss_weight*(DN_p(1,0)*crRHS89 + DN_p(1,1)*crRHS161 + DN_p(1,2)*crRHS226 - N_p[1]*crRHS10);
    rRHS[32]+=gauss_weight*(DN_p(2,0)*crRHS89 + DN_p(2,1)*crRHS161 + DN_p(2,2)*crRHS226 - N_p[2]*crRHS10);
    rRHS[33]+=gauss_weight*(DN_p(3,0)*crRHS89 + DN_p(3,1)*crRHS161 + DN_p(3,2)*crRHS226 - N_p[3]*crRHS10);

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< unsigned int TDim >
void IncompressibleNavierStokesP2P1Continuous<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}


template< unsigned int TDim >
void IncompressibleNavierStokesP2P1Continuous<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class IncompressibleNavierStokesP2P1Continuous<2>;
template class IncompressibleNavierStokesP2P1Continuous<3>;

}