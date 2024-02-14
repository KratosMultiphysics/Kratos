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
        ComputeGaussPointLHSContribution(aux_data, rLeftHandSideMatrix);
        ComputeGaussPointRHSContribution(aux_data, rRightHandSideVector);
    }
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
void IncompressibleNavierStokesDivStable<TDim>::CalculateStrainRate(ElementDataContainer& rData)
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