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
const double crLHS3 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crLHS4 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crLHS5 = rho*stab_c2*sqrt(pow(crLHS3, 2) + pow(crLHS4, 2));
const double crLHS6 = crLHS5*h/stab_c1 + mu;
const double crLHS7 = 1.0*C(0,0);
const double crLHS8 = C(0,2)*DDN_v[0](0,0);
const double crLHS9 = 1.0*DDN_v[0](0,1);
const double crLHS10 = rho*(DN_v(0,0)*crLHS3 + DN_v(0,1)*crLHS4);
const double crLHS11 = rData.BDF0*rho;
const double crLHS12 = N_v[0]*crLHS11;
const double crLHS13 = -crLHS10 - crLHS12;
const double crLHS14 = C(0,2)*crLHS9 + C(2,2)*crLHS9 + DDN_v[0](0,0)*crLHS7 + crLHS13 + 1.0*crLHS8;
const double crLHS15 = 1.0/(crLHS5/h + dyn_tau*rho/rData.DeltaTime + mu*stab_c1/pow(h, 2));
const double crLHS16 = crLHS10*crLHS15;
const double crLHS17 = rho*(DN_v(0,0)*vconv(0,0) + DN_v(0,1)*vconv(0,1) + DN_v(1,0)*vconv(1,0) + DN_v(1,1)*vconv(1,1) + DN_v(2,0)*vconv(2,0) + DN_v(2,1)*vconv(2,1) + DN_v(3,0)*vconv(3,0) + DN_v(3,1)*vconv(3,1) + DN_v(4,0)*vconv(4,0) + DN_v(4,1)*vconv(4,1) + DN_v(5,0)*vconv(5,0) + DN_v(5,1)*vconv(5,1));
const double crLHS18 = N_v[0]*crLHS17;
const double crLHS19 = crLHS15*crLHS18;
const double crLHS20 = pow(N_v[0], 2)*crLHS11 + N_v[0]*crLHS10;
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
const double crLHS175 = pow(N_v[1], 2)*crLHS11 + N_v[1]*crLHS33;
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
const double crLHS214 = pow(N_v[2], 2)*crLHS11 + N_v[2]*crLHS50;
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
const double crLHS248 = pow(N_v[3], 2)*crLHS11 + N_v[3]*crLHS67;
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
const double crLHS277 = pow(N_v[4], 2)*crLHS11 + N_v[4]*crLHS84;
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
const double crLHS301 = pow(N_v[5], 2)*crLHS11 + N_v[5]*crLHS101;
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
rLHS(0,0)+=gauss_weight*(pow(DN_v(0,0), 2)*crLHS6 + DN_v(0,0)*crLHS0 + DN_v(0,1)*crLHS2 - crLHS14*crLHS16 - crLHS14*crLHS19 + crLHS20);
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
rLHS(1,1)+=gauss_weight*(DN_v(0,0)*crLHS23 + pow(DN_v(0,1), 2)*crLHS6 + DN_v(0,1)*crLHS120 - crLHS123*crLHS16 - crLHS123*crLHS19 + crLHS20);
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
rLHS(2,2)+=gauss_weight*(pow(DN_v(1,0), 2)*crLHS6 + DN_v(1,0)*crLHS27 + DN_v(1,1)*crLHS29 - crLHS171*crLHS36 - crLHS173*crLHS36 + crLHS175);
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
rLHS(3,3)+=gauss_weight*(DN_v(1,0)*crLHS41 + pow(DN_v(1,1), 2)*crLHS6 + DN_v(1,1)*crLHS129 - crLHS132*crLHS171 - crLHS132*crLHS173 + crLHS175);
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
rLHS(4,4)+=gauss_weight*(pow(DN_v(2,0), 2)*crLHS6 + DN_v(2,0)*crLHS44 + DN_v(2,1)*crLHS46 - crLHS209*crLHS53 - crLHS211*crLHS53 + crLHS214);
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
rLHS(5,5)+=gauss_weight*(DN_v(2,0)*crLHS58 + pow(DN_v(2,1), 2)*crLHS6 + DN_v(2,1)*crLHS137 - crLHS140*crLHS209 - crLHS140*crLHS211 + crLHS214);
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
rLHS(6,6)+=gauss_weight*(pow(DN_v(3,0), 2)*crLHS6 + DN_v(3,0)*crLHS61 + DN_v(3,1)*crLHS63 - crLHS242*crLHS70 - crLHS244*crLHS70 + crLHS248);
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
rLHS(7,7)+=gauss_weight*(DN_v(3,0)*crLHS75 + pow(DN_v(3,1), 2)*crLHS6 + DN_v(3,1)*crLHS145 - crLHS148*crLHS242 - crLHS148*crLHS244 + crLHS248);
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
rLHS(8,8)+=gauss_weight*(pow(DN_v(4,0), 2)*crLHS6 + DN_v(4,0)*crLHS78 + DN_v(4,1)*crLHS80 - crLHS270*crLHS87 - crLHS272*crLHS87 + crLHS277);
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
rLHS(9,9)+=gauss_weight*(DN_v(4,0)*crLHS92 + pow(DN_v(4,1), 2)*crLHS6 + DN_v(4,1)*crLHS153 - crLHS156*crLHS270 - crLHS156*crLHS272 + crLHS277);
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
rLHS(10,10)+=gauss_weight*(pow(DN_v(5,0), 2)*crLHS6 + DN_v(5,0)*crLHS95 + DN_v(5,1)*crLHS97 - crLHS103*crLHS293 - crLHS103*crLHS295 + crLHS301);
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
rLHS(11,11)+=gauss_weight*(DN_v(5,0)*crLHS108 + pow(DN_v(5,1), 2)*crLHS6 + DN_v(5,1)*crLHS161 - crLHS164*crLHS293 - crLHS164*crLHS295 + crLHS301);
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
rLHS(12,12)+=crLHS309*(pow(DN_p(0,0), 2) + pow(DN_p(0,1), 2));
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
rLHS(13,13)+=crLHS309*(pow(DN_p(1,0), 2) + pow(DN_p(1,1), 2));
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
rLHS(14,14)+=crLHS309*(pow(DN_p(2,0), 2) + pow(DN_p(2,1), 2));

}

template <>
void IncompressibleNavierStokesP2P1Continuous<3>::AddGaussPointLeftHandSideContribution(
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
void IncompressibleNavierStokesP2P1Continuous<2>::AddGaussPointRightHandSideContribution(
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
const double crRHS7 = DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1);
const double crRHS8 = crRHS3 + crRHS7;
const double crRHS9 = rho*stab_c2*sqrt(pow(crRHS4, 2) + pow(crRHS5, 2));
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
const double crRHS34 = 1.0/(crRHS9/h + dyn_tau*rho/rData.DeltaTime + mu*stab_c1/pow(h, 2));
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