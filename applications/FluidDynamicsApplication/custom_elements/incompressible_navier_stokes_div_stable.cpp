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
    SetElementData(aux_data);

    // Initialize constitutive law parameters
    const auto& r_geom = this->GetGeometry();
    const auto p_prop = this->GetProperties();
    ConstitutiveLaw::Parameters cons_law_params(r_geom, p_prop, rCurrentProcessInfo);
    cons_law_params.SetStrainVector(aux_data.StrainRate);
    cons_law_params.SetStressVector(aux_data.ShearStress);
    cons_law_params.SetConstitutiveMatrix(aux_data.ConstitutiveMatrix);

    // Calculate kinematics
    //TODO: Bubble is missing here
    Vector weights;
    Matrix velocity_N;
    Matrix pressure_N;
    GeometryType::ShapeFunctionsGradientsType velocity_DN;
    GeometryType::ShapeFunctionsGradientsType pressure_DN;
    CalculateKinematics(weights, velocity_N, pressure_N, velocity_DN, pressure_DN);

    // Loop Gauss points
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(IntegrationMethod);
    for (IndexType g = 0; g < n_gauss; ++g) {
        // set current Gauss point kinematics
        noalias(aux_data.N_v) = row(velocity_N, g);
        noalias(aux_data.N_p) = row(pressure_N, g);
        noalias(aux_data.DN_v) = velocity_DN[g];
        noalias(aux_data.DN_p) = pressure_DN[g];
        aux_data.Weight = weights[g];

        // Calculate current Gauss point material response
        CalculateStrainRate(aux_data);
        mpConstitutiveLaw->CalculateMaterialResponseCauchy(cons_law_params);

        // Assemble standard Galerkin contribution
        ComputeGaussPointLHSContribution(aux_data, rLeftHandSideMatrix);
        ComputeGaussPointRHSContribution(aux_data, rRightHandSideVector);

        // Assemble bubble function contributions (to be condensed)
        //TODO:
    }

    // Condense bubble function contribution
    //TODO:
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
void IncompressibleNavierStokesDivStable<TDim>::SetElementData(ElementDataContainer &rData)
{
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
}

template< unsigned int TDim >
void IncompressibleNavierStokesDivStable<TDim>::CalculateKinematics(
    Vector& rGaussWeights,
    Matrix& rVelocityN,
    Matrix& rPressureN,
    GeometryType::ShapeFunctionsGradientsType& rVelocityDNDX,
    GeometryType::ShapeFunctionsGradientsType& rPressureDNDX)
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

    // Calculate enrichment bubble kinematics
    //TODO:

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
    ElementDataContainer& rData,
    MatrixType& rLHS)
{
    // Get material data
    const double rho = this->GetProperties().GetValue(DENSITY);

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double dt = rData.DeltaTime;

    // Calculate convective velocity
    const BoundedMatrix<double,2,6> vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const BoundedMatrix<double,3,3>& C = rData.ConstitutiveMatrix;

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
const double crLHS8 = pow(N_v[0], 2)*crLHS7 + crLHS5*crLHS6;
const double crLHS9 = C(0,1)*DN_v(0,1) + crLHS1;
const double crLHS10 = C(1,2)*DN_v(0,1);
const double crLHS11 = C(2,2)*DN_v(0,0) + crLHS10;
const double crLHS12 = C(0,0)*DN_v(1,0) + C(0,2)*DN_v(1,1);
const double crLHS13 = C(0,2)*DN_v(1,0);
const double crLHS14 = C(2,2)*DN_v(1,1) + crLHS13;
const double crLHS15 = DN_v(1,0)*crLHS3 + DN_v(1,1)*crLHS4;
const double crLHS16 = crLHS6*rData.BDF0;
const double crLHS17 = N_v[1]*crLHS16;
const double crLHS18 = crLHS15*crLHS6 + crLHS17;
const double crLHS19 = C(0,1)*DN_v(1,1) + crLHS13;
const double crLHS20 = C(1,2)*DN_v(1,1);
const double crLHS21 = C(2,2)*DN_v(1,0) + crLHS20;
const double crLHS22 = C(0,0)*DN_v(2,0) + C(0,2)*DN_v(2,1);
const double crLHS23 = C(0,2)*DN_v(2,0);
const double crLHS24 = C(2,2)*DN_v(2,1) + crLHS23;
const double crLHS25 = DN_v(2,0)*crLHS3 + DN_v(2,1)*crLHS4;
const double crLHS26 = N_v[2]*crLHS16;
const double crLHS27 = crLHS25*crLHS6 + crLHS26;
const double crLHS28 = C(0,1)*DN_v(2,1) + crLHS23;
const double crLHS29 = C(1,2)*DN_v(2,1);
const double crLHS30 = C(2,2)*DN_v(2,0) + crLHS29;
const double crLHS31 = C(0,0)*DN_v(3,0) + C(0,2)*DN_v(3,1);
const double crLHS32 = C(0,2)*DN_v(3,0);
const double crLHS33 = C(2,2)*DN_v(3,1) + crLHS32;
const double crLHS34 = DN_v(3,0)*crLHS3 + DN_v(3,1)*crLHS4;
const double crLHS35 = N_v[3]*crLHS16;
const double crLHS36 = crLHS34*crLHS6 + crLHS35;
const double crLHS37 = C(0,1)*DN_v(3,1) + crLHS32;
const double crLHS38 = C(1,2)*DN_v(3,1);
const double crLHS39 = C(2,2)*DN_v(3,0) + crLHS38;
const double crLHS40 = C(0,0)*DN_v(4,0) + C(0,2)*DN_v(4,1);
const double crLHS41 = C(0,2)*DN_v(4,0);
const double crLHS42 = C(2,2)*DN_v(4,1) + crLHS41;
const double crLHS43 = DN_v(4,0)*crLHS3 + DN_v(4,1)*crLHS4;
const double crLHS44 = N_v[4]*crLHS16;
const double crLHS45 = crLHS43*crLHS6 + crLHS44;
const double crLHS46 = C(0,1)*DN_v(4,1) + crLHS41;
const double crLHS47 = C(1,2)*DN_v(4,1);
const double crLHS48 = C(2,2)*DN_v(4,0) + crLHS47;
const double crLHS49 = C(0,0)*DN_v(5,0) + C(0,2)*DN_v(5,1);
const double crLHS50 = C(0,2)*DN_v(5,0);
const double crLHS51 = C(2,2)*DN_v(5,1) + crLHS50;
const double crLHS52 = DN_v(5,0)*crLHS3 + DN_v(5,1)*crLHS4;
const double crLHS53 = N_v[5]*crLHS16;
const double crLHS54 = crLHS52*crLHS6 + crLHS53;
const double crLHS55 = C(0,1)*DN_v(5,1) + crLHS50;
const double crLHS56 = C(1,2)*DN_v(5,1);
const double crLHS57 = C(2,2)*DN_v(5,0) + crLHS56;
const double crLHS58 = DN_v(0,0)*gauss_weight;
const double crLHS59 = N_p[0]*crLHS58;
const double crLHS60 = N_p[1]*crLHS58;
const double crLHS61 = N_p[2]*crLHS58;
const double crLHS62 = C(0,1)*DN_v(0,0) + crLHS10;
const double crLHS63 = C(1,1)*DN_v(0,1) + C(1,2)*DN_v(0,0);
const double crLHS64 = C(0,1)*DN_v(1,0) + crLHS20;
const double crLHS65 = C(1,1)*DN_v(1,1) + C(1,2)*DN_v(1,0);
const double crLHS66 = C(0,1)*DN_v(2,0) + crLHS29;
const double crLHS67 = C(1,1)*DN_v(2,1) + C(1,2)*DN_v(2,0);
const double crLHS68 = C(0,1)*DN_v(3,0) + crLHS38;
const double crLHS69 = C(1,1)*DN_v(3,1) + C(1,2)*DN_v(3,0);
const double crLHS70 = C(0,1)*DN_v(4,0) + crLHS47;
const double crLHS71 = C(1,1)*DN_v(4,1) + C(1,2)*DN_v(4,0);
const double crLHS72 = C(0,1)*DN_v(5,0) + crLHS56;
const double crLHS73 = C(1,1)*DN_v(5,1) + C(1,2)*DN_v(5,0);
const double crLHS74 = DN_v(0,1)*gauss_weight;
const double crLHS75 = N_p[0]*crLHS74;
const double crLHS76 = N_p[1]*crLHS74;
const double crLHS77 = N_p[2]*crLHS74;
const double crLHS78 = N_v[1]*rho;
const double crLHS79 = crLHS17 + crLHS5*crLHS78;
const double crLHS80 = pow(N_v[1], 2)*crLHS7 + crLHS15*crLHS78;
const double crLHS81 = N_v[1]*crLHS7;
const double crLHS82 = N_v[2]*crLHS81;
const double crLHS83 = crLHS25*crLHS78 + crLHS82;
const double crLHS84 = N_v[3]*crLHS81;
const double crLHS85 = crLHS34*crLHS78 + crLHS84;
const double crLHS86 = N_v[4]*crLHS81;
const double crLHS87 = crLHS43*crLHS78 + crLHS86;
const double crLHS88 = N_v[5]*crLHS81;
const double crLHS89 = crLHS52*crLHS78 + crLHS88;
const double crLHS90 = DN_v(1,0)*gauss_weight;
const double crLHS91 = N_p[0]*crLHS90;
const double crLHS92 = N_p[1]*crLHS90;
const double crLHS93 = N_p[2]*crLHS90;
const double crLHS94 = DN_v(1,1)*gauss_weight;
const double crLHS95 = N_p[0]*crLHS94;
const double crLHS96 = N_p[1]*crLHS94;
const double crLHS97 = N_p[2]*crLHS94;
const double crLHS98 = N_v[2]*rho;
const double crLHS99 = crLHS26 + crLHS5*crLHS98;
const double crLHS100 = crLHS15*crLHS98 + crLHS82;
const double crLHS101 = pow(N_v[2], 2)*crLHS7 + crLHS25*crLHS98;
const double crLHS102 = N_v[2]*crLHS7;
const double crLHS103 = N_v[3]*crLHS102;
const double crLHS104 = crLHS103 + crLHS34*crLHS98;
const double crLHS105 = N_v[4]*crLHS102;
const double crLHS106 = crLHS105 + crLHS43*crLHS98;
const double crLHS107 = N_v[5]*crLHS102;
const double crLHS108 = crLHS107 + crLHS52*crLHS98;
const double crLHS109 = DN_v(2,0)*gauss_weight;
const double crLHS110 = N_p[0]*crLHS109;
const double crLHS111 = N_p[1]*crLHS109;
const double crLHS112 = N_p[2]*crLHS109;
const double crLHS113 = DN_v(2,1)*gauss_weight;
const double crLHS114 = N_p[0]*crLHS113;
const double crLHS115 = N_p[1]*crLHS113;
const double crLHS116 = N_p[2]*crLHS113;
const double crLHS117 = N_v[3]*rho;
const double crLHS118 = crLHS117*crLHS5 + crLHS35;
const double crLHS119 = crLHS117*crLHS15 + crLHS84;
const double crLHS120 = crLHS103 + crLHS117*crLHS25;
const double crLHS121 = pow(N_v[3], 2)*crLHS7 + crLHS117*crLHS34;
const double crLHS122 = N_v[3]*crLHS7;
const double crLHS123 = N_v[4]*crLHS122;
const double crLHS124 = crLHS117*crLHS43 + crLHS123;
const double crLHS125 = N_v[5]*crLHS122;
const double crLHS126 = crLHS117*crLHS52 + crLHS125;
const double crLHS127 = DN_v(3,0)*gauss_weight;
const double crLHS128 = N_p[0]*crLHS127;
const double crLHS129 = N_p[1]*crLHS127;
const double crLHS130 = N_p[2]*crLHS127;
const double crLHS131 = DN_v(3,1)*gauss_weight;
const double crLHS132 = N_p[0]*crLHS131;
const double crLHS133 = N_p[1]*crLHS131;
const double crLHS134 = N_p[2]*crLHS131;
const double crLHS135 = N_v[4]*rho;
const double crLHS136 = crLHS135*crLHS5 + crLHS44;
const double crLHS137 = crLHS135*crLHS15 + crLHS86;
const double crLHS138 = crLHS105 + crLHS135*crLHS25;
const double crLHS139 = crLHS123 + crLHS135*crLHS34;
const double crLHS140 = pow(N_v[4], 2)*crLHS7 + crLHS135*crLHS43;
const double crLHS141 = N_v[4]*N_v[5]*crLHS7;
const double crLHS142 = crLHS135*crLHS52 + crLHS141;
const double crLHS143 = DN_v(4,0)*gauss_weight;
const double crLHS144 = N_p[0]*crLHS143;
const double crLHS145 = N_p[1]*crLHS143;
const double crLHS146 = N_p[2]*crLHS143;
const double crLHS147 = DN_v(4,1)*gauss_weight;
const double crLHS148 = N_p[0]*crLHS147;
const double crLHS149 = N_p[1]*crLHS147;
const double crLHS150 = N_p[2]*crLHS147;
const double crLHS151 = N_v[5]*rho;
const double crLHS152 = crLHS151*crLHS5 + crLHS53;
const double crLHS153 = crLHS15*crLHS151 + crLHS88;
const double crLHS154 = crLHS107 + crLHS151*crLHS25;
const double crLHS155 = crLHS125 + crLHS151*crLHS34;
const double crLHS156 = crLHS141 + crLHS151*crLHS43;
const double crLHS157 = pow(N_v[5], 2)*crLHS7 + crLHS151*crLHS52;
const double crLHS158 = DN_v(5,0)*gauss_weight;
const double crLHS159 = N_p[0]*crLHS158;
const double crLHS160 = N_p[1]*crLHS158;
const double crLHS161 = N_p[2]*crLHS158;
const double crLHS162 = DN_v(5,1)*gauss_weight;
const double crLHS163 = N_p[0]*crLHS162;
const double crLHS164 = N_p[1]*crLHS162;
const double crLHS165 = N_p[2]*crLHS162;
rLHS(0,0)+=gauss_weight*(DN_v(0,0)*crLHS0 + DN_v(0,1)*crLHS2 + crLHS8);
rLHS(0,1)+=gauss_weight*(DN_v(0,0)*crLHS9 + DN_v(0,1)*crLHS11);
rLHS(0,2)+=gauss_weight*(DN_v(0,0)*crLHS12 + DN_v(0,1)*crLHS14 + crLHS18);
rLHS(0,3)+=gauss_weight*(DN_v(0,0)*crLHS19 + DN_v(0,1)*crLHS21);
rLHS(0,4)+=gauss_weight*(DN_v(0,0)*crLHS22 + DN_v(0,1)*crLHS24 + crLHS27);
rLHS(0,5)+=gauss_weight*(DN_v(0,0)*crLHS28 + DN_v(0,1)*crLHS30);
rLHS(0,6)+=gauss_weight*(DN_v(0,0)*crLHS31 + DN_v(0,1)*crLHS33 + crLHS36);
rLHS(0,7)+=gauss_weight*(DN_v(0,0)*crLHS37 + DN_v(0,1)*crLHS39);
rLHS(0,8)+=gauss_weight*(DN_v(0,0)*crLHS40 + DN_v(0,1)*crLHS42 + crLHS45);
rLHS(0,9)+=gauss_weight*(DN_v(0,0)*crLHS46 + DN_v(0,1)*crLHS48);
rLHS(0,10)+=gauss_weight*(DN_v(0,0)*crLHS49 + DN_v(0,1)*crLHS51 + crLHS54);
rLHS(0,11)+=gauss_weight*(DN_v(0,0)*crLHS55 + DN_v(0,1)*crLHS57);
rLHS(0,12)+=-crLHS59;
rLHS(0,13)+=-crLHS60;
rLHS(0,14)+=-crLHS61;
rLHS(1,0)+=gauss_weight*(DN_v(0,0)*crLHS2 + DN_v(0,1)*crLHS62);
rLHS(1,1)+=gauss_weight*(DN_v(0,0)*crLHS11 + DN_v(0,1)*crLHS63 + crLHS8);
rLHS(1,2)+=gauss_weight*(DN_v(0,0)*crLHS14 + DN_v(0,1)*crLHS64);
rLHS(1,3)+=gauss_weight*(DN_v(0,0)*crLHS21 + DN_v(0,1)*crLHS65 + crLHS18);
rLHS(1,4)+=gauss_weight*(DN_v(0,0)*crLHS24 + DN_v(0,1)*crLHS66);
rLHS(1,5)+=gauss_weight*(DN_v(0,0)*crLHS30 + DN_v(0,1)*crLHS67 + crLHS27);
rLHS(1,6)+=gauss_weight*(DN_v(0,0)*crLHS33 + DN_v(0,1)*crLHS68);
rLHS(1,7)+=gauss_weight*(DN_v(0,0)*crLHS39 + DN_v(0,1)*crLHS69 + crLHS36);
rLHS(1,8)+=gauss_weight*(DN_v(0,0)*crLHS42 + DN_v(0,1)*crLHS70);
rLHS(1,9)+=gauss_weight*(DN_v(0,0)*crLHS48 + DN_v(0,1)*crLHS71 + crLHS45);
rLHS(1,10)+=gauss_weight*(DN_v(0,0)*crLHS51 + DN_v(0,1)*crLHS72);
rLHS(1,11)+=gauss_weight*(DN_v(0,0)*crLHS57 + DN_v(0,1)*crLHS73 + crLHS54);
rLHS(1,12)+=-crLHS75;
rLHS(1,13)+=-crLHS76;
rLHS(1,14)+=-crLHS77;
rLHS(2,0)+=gauss_weight*(DN_v(1,0)*crLHS0 + DN_v(1,1)*crLHS2 + crLHS79);
rLHS(2,1)+=gauss_weight*(DN_v(1,0)*crLHS9 + DN_v(1,1)*crLHS11);
rLHS(2,2)+=gauss_weight*(DN_v(1,0)*crLHS12 + DN_v(1,1)*crLHS14 + crLHS80);
rLHS(2,3)+=gauss_weight*(DN_v(1,0)*crLHS19 + DN_v(1,1)*crLHS21);
rLHS(2,4)+=gauss_weight*(DN_v(1,0)*crLHS22 + DN_v(1,1)*crLHS24 + crLHS83);
rLHS(2,5)+=gauss_weight*(DN_v(1,0)*crLHS28 + DN_v(1,1)*crLHS30);
rLHS(2,6)+=gauss_weight*(DN_v(1,0)*crLHS31 + DN_v(1,1)*crLHS33 + crLHS85);
rLHS(2,7)+=gauss_weight*(DN_v(1,0)*crLHS37 + DN_v(1,1)*crLHS39);
rLHS(2,8)+=gauss_weight*(DN_v(1,0)*crLHS40 + DN_v(1,1)*crLHS42 + crLHS87);
rLHS(2,9)+=gauss_weight*(DN_v(1,0)*crLHS46 + DN_v(1,1)*crLHS48);
rLHS(2,10)+=gauss_weight*(DN_v(1,0)*crLHS49 + DN_v(1,1)*crLHS51 + crLHS89);
rLHS(2,11)+=gauss_weight*(DN_v(1,0)*crLHS55 + DN_v(1,1)*crLHS57);
rLHS(2,12)+=-crLHS91;
rLHS(2,13)+=-crLHS92;
rLHS(2,14)+=-crLHS93;
rLHS(3,0)+=gauss_weight*(DN_v(1,0)*crLHS2 + DN_v(1,1)*crLHS62);
rLHS(3,1)+=gauss_weight*(DN_v(1,0)*crLHS11 + DN_v(1,1)*crLHS63 + crLHS79);
rLHS(3,2)+=gauss_weight*(DN_v(1,0)*crLHS14 + DN_v(1,1)*crLHS64);
rLHS(3,3)+=gauss_weight*(DN_v(1,0)*crLHS21 + DN_v(1,1)*crLHS65 + crLHS80);
rLHS(3,4)+=gauss_weight*(DN_v(1,0)*crLHS24 + DN_v(1,1)*crLHS66);
rLHS(3,5)+=gauss_weight*(DN_v(1,0)*crLHS30 + DN_v(1,1)*crLHS67 + crLHS83);
rLHS(3,6)+=gauss_weight*(DN_v(1,0)*crLHS33 + DN_v(1,1)*crLHS68);
rLHS(3,7)+=gauss_weight*(DN_v(1,0)*crLHS39 + DN_v(1,1)*crLHS69 + crLHS85);
rLHS(3,8)+=gauss_weight*(DN_v(1,0)*crLHS42 + DN_v(1,1)*crLHS70);
rLHS(3,9)+=gauss_weight*(DN_v(1,0)*crLHS48 + DN_v(1,1)*crLHS71 + crLHS87);
rLHS(3,10)+=gauss_weight*(DN_v(1,0)*crLHS51 + DN_v(1,1)*crLHS72);
rLHS(3,11)+=gauss_weight*(DN_v(1,0)*crLHS57 + DN_v(1,1)*crLHS73 + crLHS89);
rLHS(3,12)+=-crLHS95;
rLHS(3,13)+=-crLHS96;
rLHS(3,14)+=-crLHS97;
rLHS(4,0)+=gauss_weight*(DN_v(2,0)*crLHS0 + DN_v(2,1)*crLHS2 + crLHS99);
rLHS(4,1)+=gauss_weight*(DN_v(2,0)*crLHS9 + DN_v(2,1)*crLHS11);
rLHS(4,2)+=gauss_weight*(DN_v(2,0)*crLHS12 + DN_v(2,1)*crLHS14 + crLHS100);
rLHS(4,3)+=gauss_weight*(DN_v(2,0)*crLHS19 + DN_v(2,1)*crLHS21);
rLHS(4,4)+=gauss_weight*(DN_v(2,0)*crLHS22 + DN_v(2,1)*crLHS24 + crLHS101);
rLHS(4,5)+=gauss_weight*(DN_v(2,0)*crLHS28 + DN_v(2,1)*crLHS30);
rLHS(4,6)+=gauss_weight*(DN_v(2,0)*crLHS31 + DN_v(2,1)*crLHS33 + crLHS104);
rLHS(4,7)+=gauss_weight*(DN_v(2,0)*crLHS37 + DN_v(2,1)*crLHS39);
rLHS(4,8)+=gauss_weight*(DN_v(2,0)*crLHS40 + DN_v(2,1)*crLHS42 + crLHS106);
rLHS(4,9)+=gauss_weight*(DN_v(2,0)*crLHS46 + DN_v(2,1)*crLHS48);
rLHS(4,10)+=gauss_weight*(DN_v(2,0)*crLHS49 + DN_v(2,1)*crLHS51 + crLHS108);
rLHS(4,11)+=gauss_weight*(DN_v(2,0)*crLHS55 + DN_v(2,1)*crLHS57);
rLHS(4,12)+=-crLHS110;
rLHS(4,13)+=-crLHS111;
rLHS(4,14)+=-crLHS112;
rLHS(5,0)+=gauss_weight*(DN_v(2,0)*crLHS2 + DN_v(2,1)*crLHS62);
rLHS(5,1)+=gauss_weight*(DN_v(2,0)*crLHS11 + DN_v(2,1)*crLHS63 + crLHS99);
rLHS(5,2)+=gauss_weight*(DN_v(2,0)*crLHS14 + DN_v(2,1)*crLHS64);
rLHS(5,3)+=gauss_weight*(DN_v(2,0)*crLHS21 + DN_v(2,1)*crLHS65 + crLHS100);
rLHS(5,4)+=gauss_weight*(DN_v(2,0)*crLHS24 + DN_v(2,1)*crLHS66);
rLHS(5,5)+=gauss_weight*(DN_v(2,0)*crLHS30 + DN_v(2,1)*crLHS67 + crLHS101);
rLHS(5,6)+=gauss_weight*(DN_v(2,0)*crLHS33 + DN_v(2,1)*crLHS68);
rLHS(5,7)+=gauss_weight*(DN_v(2,0)*crLHS39 + DN_v(2,1)*crLHS69 + crLHS104);
rLHS(5,8)+=gauss_weight*(DN_v(2,0)*crLHS42 + DN_v(2,1)*crLHS70);
rLHS(5,9)+=gauss_weight*(DN_v(2,0)*crLHS48 + DN_v(2,1)*crLHS71 + crLHS106);
rLHS(5,10)+=gauss_weight*(DN_v(2,0)*crLHS51 + DN_v(2,1)*crLHS72);
rLHS(5,11)+=gauss_weight*(DN_v(2,0)*crLHS57 + DN_v(2,1)*crLHS73 + crLHS108);
rLHS(5,12)+=-crLHS114;
rLHS(5,13)+=-crLHS115;
rLHS(5,14)+=-crLHS116;
rLHS(6,0)+=gauss_weight*(DN_v(3,0)*crLHS0 + DN_v(3,1)*crLHS2 + crLHS118);
rLHS(6,1)+=gauss_weight*(DN_v(3,0)*crLHS9 + DN_v(3,1)*crLHS11);
rLHS(6,2)+=gauss_weight*(DN_v(3,0)*crLHS12 + DN_v(3,1)*crLHS14 + crLHS119);
rLHS(6,3)+=gauss_weight*(DN_v(3,0)*crLHS19 + DN_v(3,1)*crLHS21);
rLHS(6,4)+=gauss_weight*(DN_v(3,0)*crLHS22 + DN_v(3,1)*crLHS24 + crLHS120);
rLHS(6,5)+=gauss_weight*(DN_v(3,0)*crLHS28 + DN_v(3,1)*crLHS30);
rLHS(6,6)+=gauss_weight*(DN_v(3,0)*crLHS31 + DN_v(3,1)*crLHS33 + crLHS121);
rLHS(6,7)+=gauss_weight*(DN_v(3,0)*crLHS37 + DN_v(3,1)*crLHS39);
rLHS(6,8)+=gauss_weight*(DN_v(3,0)*crLHS40 + DN_v(3,1)*crLHS42 + crLHS124);
rLHS(6,9)+=gauss_weight*(DN_v(3,0)*crLHS46 + DN_v(3,1)*crLHS48);
rLHS(6,10)+=gauss_weight*(DN_v(3,0)*crLHS49 + DN_v(3,1)*crLHS51 + crLHS126);
rLHS(6,11)+=gauss_weight*(DN_v(3,0)*crLHS55 + DN_v(3,1)*crLHS57);
rLHS(6,12)+=-crLHS128;
rLHS(6,13)+=-crLHS129;
rLHS(6,14)+=-crLHS130;
rLHS(7,0)+=gauss_weight*(DN_v(3,0)*crLHS2 + DN_v(3,1)*crLHS62);
rLHS(7,1)+=gauss_weight*(DN_v(3,0)*crLHS11 + DN_v(3,1)*crLHS63 + crLHS118);
rLHS(7,2)+=gauss_weight*(DN_v(3,0)*crLHS14 + DN_v(3,1)*crLHS64);
rLHS(7,3)+=gauss_weight*(DN_v(3,0)*crLHS21 + DN_v(3,1)*crLHS65 + crLHS119);
rLHS(7,4)+=gauss_weight*(DN_v(3,0)*crLHS24 + DN_v(3,1)*crLHS66);
rLHS(7,5)+=gauss_weight*(DN_v(3,0)*crLHS30 + DN_v(3,1)*crLHS67 + crLHS120);
rLHS(7,6)+=gauss_weight*(DN_v(3,0)*crLHS33 + DN_v(3,1)*crLHS68);
rLHS(7,7)+=gauss_weight*(DN_v(3,0)*crLHS39 + DN_v(3,1)*crLHS69 + crLHS121);
rLHS(7,8)+=gauss_weight*(DN_v(3,0)*crLHS42 + DN_v(3,1)*crLHS70);
rLHS(7,9)+=gauss_weight*(DN_v(3,0)*crLHS48 + DN_v(3,1)*crLHS71 + crLHS124);
rLHS(7,10)+=gauss_weight*(DN_v(3,0)*crLHS51 + DN_v(3,1)*crLHS72);
rLHS(7,11)+=gauss_weight*(DN_v(3,0)*crLHS57 + DN_v(3,1)*crLHS73 + crLHS126);
rLHS(7,12)+=-crLHS132;
rLHS(7,13)+=-crLHS133;
rLHS(7,14)+=-crLHS134;
rLHS(8,0)+=gauss_weight*(DN_v(4,0)*crLHS0 + DN_v(4,1)*crLHS2 + crLHS136);
rLHS(8,1)+=gauss_weight*(DN_v(4,0)*crLHS9 + DN_v(4,1)*crLHS11);
rLHS(8,2)+=gauss_weight*(DN_v(4,0)*crLHS12 + DN_v(4,1)*crLHS14 + crLHS137);
rLHS(8,3)+=gauss_weight*(DN_v(4,0)*crLHS19 + DN_v(4,1)*crLHS21);
rLHS(8,4)+=gauss_weight*(DN_v(4,0)*crLHS22 + DN_v(4,1)*crLHS24 + crLHS138);
rLHS(8,5)+=gauss_weight*(DN_v(4,0)*crLHS28 + DN_v(4,1)*crLHS30);
rLHS(8,6)+=gauss_weight*(DN_v(4,0)*crLHS31 + DN_v(4,1)*crLHS33 + crLHS139);
rLHS(8,7)+=gauss_weight*(DN_v(4,0)*crLHS37 + DN_v(4,1)*crLHS39);
rLHS(8,8)+=gauss_weight*(DN_v(4,0)*crLHS40 + DN_v(4,1)*crLHS42 + crLHS140);
rLHS(8,9)+=gauss_weight*(DN_v(4,0)*crLHS46 + DN_v(4,1)*crLHS48);
rLHS(8,10)+=gauss_weight*(DN_v(4,0)*crLHS49 + DN_v(4,1)*crLHS51 + crLHS142);
rLHS(8,11)+=gauss_weight*(DN_v(4,0)*crLHS55 + DN_v(4,1)*crLHS57);
rLHS(8,12)+=-crLHS144;
rLHS(8,13)+=-crLHS145;
rLHS(8,14)+=-crLHS146;
rLHS(9,0)+=gauss_weight*(DN_v(4,0)*crLHS2 + DN_v(4,1)*crLHS62);
rLHS(9,1)+=gauss_weight*(DN_v(4,0)*crLHS11 + DN_v(4,1)*crLHS63 + crLHS136);
rLHS(9,2)+=gauss_weight*(DN_v(4,0)*crLHS14 + DN_v(4,1)*crLHS64);
rLHS(9,3)+=gauss_weight*(DN_v(4,0)*crLHS21 + DN_v(4,1)*crLHS65 + crLHS137);
rLHS(9,4)+=gauss_weight*(DN_v(4,0)*crLHS24 + DN_v(4,1)*crLHS66);
rLHS(9,5)+=gauss_weight*(DN_v(4,0)*crLHS30 + DN_v(4,1)*crLHS67 + crLHS138);
rLHS(9,6)+=gauss_weight*(DN_v(4,0)*crLHS33 + DN_v(4,1)*crLHS68);
rLHS(9,7)+=gauss_weight*(DN_v(4,0)*crLHS39 + DN_v(4,1)*crLHS69 + crLHS139);
rLHS(9,8)+=gauss_weight*(DN_v(4,0)*crLHS42 + DN_v(4,1)*crLHS70);
rLHS(9,9)+=gauss_weight*(DN_v(4,0)*crLHS48 + DN_v(4,1)*crLHS71 + crLHS140);
rLHS(9,10)+=gauss_weight*(DN_v(4,0)*crLHS51 + DN_v(4,1)*crLHS72);
rLHS(9,11)+=gauss_weight*(DN_v(4,0)*crLHS57 + DN_v(4,1)*crLHS73 + crLHS142);
rLHS(9,12)+=-crLHS148;
rLHS(9,13)+=-crLHS149;
rLHS(9,14)+=-crLHS150;
rLHS(10,0)+=gauss_weight*(DN_v(5,0)*crLHS0 + DN_v(5,1)*crLHS2 + crLHS152);
rLHS(10,1)+=gauss_weight*(DN_v(5,0)*crLHS9 + DN_v(5,1)*crLHS11);
rLHS(10,2)+=gauss_weight*(DN_v(5,0)*crLHS12 + DN_v(5,1)*crLHS14 + crLHS153);
rLHS(10,3)+=gauss_weight*(DN_v(5,0)*crLHS19 + DN_v(5,1)*crLHS21);
rLHS(10,4)+=gauss_weight*(DN_v(5,0)*crLHS22 + DN_v(5,1)*crLHS24 + crLHS154);
rLHS(10,5)+=gauss_weight*(DN_v(5,0)*crLHS28 + DN_v(5,1)*crLHS30);
rLHS(10,6)+=gauss_weight*(DN_v(5,0)*crLHS31 + DN_v(5,1)*crLHS33 + crLHS155);
rLHS(10,7)+=gauss_weight*(DN_v(5,0)*crLHS37 + DN_v(5,1)*crLHS39);
rLHS(10,8)+=gauss_weight*(DN_v(5,0)*crLHS40 + DN_v(5,1)*crLHS42 + crLHS156);
rLHS(10,9)+=gauss_weight*(DN_v(5,0)*crLHS46 + DN_v(5,1)*crLHS48);
rLHS(10,10)+=gauss_weight*(DN_v(5,0)*crLHS49 + DN_v(5,1)*crLHS51 + crLHS157);
rLHS(10,11)+=gauss_weight*(DN_v(5,0)*crLHS55 + DN_v(5,1)*crLHS57);
rLHS(10,12)+=-crLHS159;
rLHS(10,13)+=-crLHS160;
rLHS(10,14)+=-crLHS161;
rLHS(11,0)+=gauss_weight*(DN_v(5,0)*crLHS2 + DN_v(5,1)*crLHS62);
rLHS(11,1)+=gauss_weight*(DN_v(5,0)*crLHS11 + DN_v(5,1)*crLHS63 + crLHS152);
rLHS(11,2)+=gauss_weight*(DN_v(5,0)*crLHS14 + DN_v(5,1)*crLHS64);
rLHS(11,3)+=gauss_weight*(DN_v(5,0)*crLHS21 + DN_v(5,1)*crLHS65 + crLHS153);
rLHS(11,4)+=gauss_weight*(DN_v(5,0)*crLHS24 + DN_v(5,1)*crLHS66);
rLHS(11,5)+=gauss_weight*(DN_v(5,0)*crLHS30 + DN_v(5,1)*crLHS67 + crLHS154);
rLHS(11,6)+=gauss_weight*(DN_v(5,0)*crLHS33 + DN_v(5,1)*crLHS68);
rLHS(11,7)+=gauss_weight*(DN_v(5,0)*crLHS39 + DN_v(5,1)*crLHS69 + crLHS155);
rLHS(11,8)+=gauss_weight*(DN_v(5,0)*crLHS42 + DN_v(5,1)*crLHS70);
rLHS(11,9)+=gauss_weight*(DN_v(5,0)*crLHS48 + DN_v(5,1)*crLHS71 + crLHS156);
rLHS(11,10)+=gauss_weight*(DN_v(5,0)*crLHS51 + DN_v(5,1)*crLHS72);
rLHS(11,11)+=gauss_weight*(DN_v(5,0)*crLHS57 + DN_v(5,1)*crLHS73 + crLHS157);
rLHS(11,12)+=-crLHS163;
rLHS(11,13)+=-crLHS164;
rLHS(11,14)+=-crLHS165;
rLHS(12,0)+=crLHS59;
rLHS(12,1)+=crLHS75;
rLHS(12,2)+=crLHS91;
rLHS(12,3)+=crLHS95;
rLHS(12,4)+=crLHS110;
rLHS(12,5)+=crLHS114;
rLHS(12,6)+=crLHS128;
rLHS(12,7)+=crLHS132;
rLHS(12,8)+=crLHS144;
rLHS(12,9)+=crLHS148;
rLHS(12,10)+=crLHS159;
rLHS(12,11)+=crLHS163;
rLHS(12,12)+=0;
rLHS(12,13)+=0;
rLHS(12,14)+=0;
rLHS(13,0)+=crLHS60;
rLHS(13,1)+=crLHS76;
rLHS(13,2)+=crLHS92;
rLHS(13,3)+=crLHS96;
rLHS(13,4)+=crLHS111;
rLHS(13,5)+=crLHS115;
rLHS(13,6)+=crLHS129;
rLHS(13,7)+=crLHS133;
rLHS(13,8)+=crLHS145;
rLHS(13,9)+=crLHS149;
rLHS(13,10)+=crLHS160;
rLHS(13,11)+=crLHS164;
rLHS(13,12)+=0;
rLHS(13,13)+=0;
rLHS(13,14)+=0;
rLHS(14,0)+=crLHS61;
rLHS(14,1)+=crLHS77;
rLHS(14,2)+=crLHS93;
rLHS(14,3)+=crLHS97;
rLHS(14,4)+=crLHS112;
rLHS(14,5)+=crLHS116;
rLHS(14,6)+=crLHS130;
rLHS(14,7)+=crLHS134;
rLHS(14,8)+=crLHS146;
rLHS(14,9)+=crLHS150;
rLHS(14,10)+=crLHS161;
rLHS(14,11)+=crLHS165;
rLHS(14,12)+=0;
rLHS(14,13)+=0;
rLHS(14,14)+=0;

}

template <>
void IncompressibleNavierStokesDivStable<3>::ComputeGaussPointLHSContribution(
    ElementDataContainer& rData,
    MatrixType& rLHS)
{
    // Get material data
    const double rho = this->GetProperties().GetValue(DENSITY);

    // Get auxiliary data
    const double dt = rData.DeltaTime;
    const double bdf0 = rData.BDF0;

    // Calculate convective velocity
    const BoundedMatrix<double,3,10> vconv = rData.Velocity - rData.MeshVelocity;

    // Get constitutive matrix
    const BoundedMatrix<double,6,6>& C = rData.ConstitutiveMatrix;

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
    ElementDataContainer& rData,
    VectorType& rRHS)
{
    // Get material data
    const double rho = this->GetProperties().GetValue(DENSITY);

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double bdf1 = rData.BDF1;
    const double bdf2 = rData.BDF2;
    const double dt = rData.DeltaTime;

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
const double crRHS1 = N_v[0]*r_f(0,0) + N_v[1]*r_f(1,0) + N_v[2]*r_f(2,0) + N_v[3]*r_f(3,0) + N_v[4]*r_f(4,0) + N_v[5]*r_f(5,0);
const double crRHS2 = N_v[0]*rho;
const double crRHS3 = N_v[0]*(rData.BDF0*r_v(0,0) + rData.BDF1*r_vn(0,0) + rData.BDF2*r_vnn(0,0)) + N_v[1]*(rData.BDF0*r_v(1,0) + rData.BDF1*r_vn(1,0) + rData.BDF2*r_vnn(1,0)) + N_v[2]*(rData.BDF0*r_v(2,0) + rData.BDF1*r_vn(2,0) + rData.BDF2*r_vnn(2,0)) + N_v[3]*(rData.BDF0*r_v(3,0) + rData.BDF1*r_vn(3,0) + rData.BDF2*r_vnn(3,0)) + N_v[4]*(rData.BDF0*r_v(4,0) + rData.BDF1*r_vn(4,0) + rData.BDF2*r_vnn(4,0)) + N_v[5]*(rData.BDF0*r_v(5,0) + rData.BDF1*r_vn(5,0) + rData.BDF2*r_vnn(5,0));
const double crRHS4 = DN_v(0,0)*r_v(0,0) + DN_v(1,0)*r_v(1,0) + DN_v(2,0)*r_v(2,0) + DN_v(3,0)*r_v(3,0) + DN_v(4,0)*r_v(4,0) + DN_v(5,0)*r_v(5,0);
const double crRHS5 = N_v[0]*vconv(0,0) + N_v[1]*vconv(1,0) + N_v[2]*vconv(2,0) + N_v[3]*vconv(3,0) + N_v[4]*vconv(4,0) + N_v[5]*vconv(5,0);
const double crRHS6 = N_v[0]*vconv(0,1) + N_v[1]*vconv(1,1) + N_v[2]*vconv(2,1) + N_v[3]*vconv(3,1) + N_v[4]*vconv(4,1) + N_v[5]*vconv(5,1);
const double crRHS7 = crRHS4*crRHS5 + crRHS6*(DN_v(0,1)*r_v(0,0) + DN_v(1,1)*r_v(1,0) + DN_v(2,1)*r_v(2,0) + DN_v(3,1)*r_v(3,0) + DN_v(4,1)*r_v(4,0) + DN_v(5,1)*r_v(5,0));
const double crRHS8 = N_v[0]*r_f(0,1) + N_v[1]*r_f(1,1) + N_v[2]*r_f(2,1) + N_v[3]*r_f(3,1) + N_v[4]*r_f(4,1) + N_v[5]*r_f(5,1);
const double crRHS9 = N_v[0]*(rData.BDF0*r_v(0,1) + rData.BDF1*r_vn(0,1) + rData.BDF2*r_vnn(0,1)) + N_v[1]*(rData.BDF0*r_v(1,1) + rData.BDF1*r_vn(1,1) + rData.BDF2*r_vnn(1,1)) + N_v[2]*(rData.BDF0*r_v(2,1) + rData.BDF1*r_vn(2,1) + rData.BDF2*r_vnn(2,1)) + N_v[3]*(rData.BDF0*r_v(3,1) + rData.BDF1*r_vn(3,1) + rData.BDF2*r_vnn(3,1)) + N_v[4]*(rData.BDF0*r_v(4,1) + rData.BDF1*r_vn(4,1) + rData.BDF2*r_vnn(4,1)) + N_v[5]*(rData.BDF0*r_v(5,1) + rData.BDF1*r_vn(5,1) + rData.BDF2*r_vnn(5,1));
const double crRHS10 = DN_v(0,1)*r_v(0,1) + DN_v(1,1)*r_v(1,1) + DN_v(2,1)*r_v(2,1) + DN_v(3,1)*r_v(3,1) + DN_v(4,1)*r_v(4,1) + DN_v(5,1)*r_v(5,1);
const double crRHS11 = crRHS10*crRHS6 + crRHS5*(DN_v(0,0)*r_v(0,1) + DN_v(1,0)*r_v(1,1) + DN_v(2,0)*r_v(2,1) + DN_v(3,0)*r_v(3,1) + DN_v(4,0)*r_v(4,1) + DN_v(5,0)*r_v(5,1));
const double crRHS12 = N_v[1]*rho;
const double crRHS13 = N_v[2]*rho;
const double crRHS14 = N_v[3]*rho;
const double crRHS15 = N_v[4]*rho;
const double crRHS16 = N_v[5]*rho;
const double crRHS17 = gauss_weight*(crRHS10 + crRHS4);
rRHS[0]+=-gauss_weight*(-DN_v(0,0)*crRHS0 + DN_v(0,0)*r_stress[0] + DN_v(0,1)*r_stress[2] - crRHS1*crRHS2 + crRHS2*crRHS3 + crRHS2*crRHS7);
rRHS[1]+=-gauss_weight*(DN_v(0,0)*r_stress[2] - DN_v(0,1)*crRHS0 + DN_v(0,1)*r_stress[1] + crRHS11*crRHS2 - crRHS2*crRHS8 + crRHS2*crRHS9);
rRHS[2]+=-gauss_weight*(-DN_v(1,0)*crRHS0 + DN_v(1,0)*r_stress[0] + DN_v(1,1)*r_stress[2] - crRHS1*crRHS12 + crRHS12*crRHS3 + crRHS12*crRHS7);
rRHS[3]+=-gauss_weight*(DN_v(1,0)*r_stress[2] - DN_v(1,1)*crRHS0 + DN_v(1,1)*r_stress[1] + crRHS11*crRHS12 - crRHS12*crRHS8 + crRHS12*crRHS9);
rRHS[4]+=-gauss_weight*(-DN_v(2,0)*crRHS0 + DN_v(2,0)*r_stress[0] + DN_v(2,1)*r_stress[2] - crRHS1*crRHS13 + crRHS13*crRHS3 + crRHS13*crRHS7);
rRHS[5]+=-gauss_weight*(DN_v(2,0)*r_stress[2] - DN_v(2,1)*crRHS0 + DN_v(2,1)*r_stress[1] + crRHS11*crRHS13 - crRHS13*crRHS8 + crRHS13*crRHS9);
rRHS[6]+=-gauss_weight*(-DN_v(3,0)*crRHS0 + DN_v(3,0)*r_stress[0] + DN_v(3,1)*r_stress[2] - crRHS1*crRHS14 + crRHS14*crRHS3 + crRHS14*crRHS7);
rRHS[7]+=-gauss_weight*(DN_v(3,0)*r_stress[2] - DN_v(3,1)*crRHS0 + DN_v(3,1)*r_stress[1] + crRHS11*crRHS14 - crRHS14*crRHS8 + crRHS14*crRHS9);
rRHS[8]+=-gauss_weight*(-DN_v(4,0)*crRHS0 + DN_v(4,0)*r_stress[0] + DN_v(4,1)*r_stress[2] - crRHS1*crRHS15 + crRHS15*crRHS3 + crRHS15*crRHS7);
rRHS[9]+=-gauss_weight*(DN_v(4,0)*r_stress[2] - DN_v(4,1)*crRHS0 + DN_v(4,1)*r_stress[1] + crRHS11*crRHS15 - crRHS15*crRHS8 + crRHS15*crRHS9);
rRHS[10]+=-gauss_weight*(-DN_v(5,0)*crRHS0 + DN_v(5,0)*r_stress[0] + DN_v(5,1)*r_stress[2] - crRHS1*crRHS16 + crRHS16*crRHS3 + crRHS16*crRHS7);
rRHS[11]+=-gauss_weight*(DN_v(5,0)*r_stress[2] - DN_v(5,1)*crRHS0 + DN_v(5,1)*r_stress[1] + crRHS11*crRHS16 - crRHS16*crRHS8 + crRHS16*crRHS9);
rRHS[12]+=-N_p[0]*crRHS17;
rRHS[13]+=-N_p[1]*crRHS17;
rRHS[14]+=-N_p[2]*crRHS17;

}

template <>
void IncompressibleNavierStokesDivStable<3>::ComputeGaussPointRHSContribution(
    ElementDataContainer& rData,
    VectorType& rRHS)
{
    // Get material data
    const double rho = this->GetProperties().GetValue(DENSITY);

    // Get auxiliary data
    const double bdf0 = rData.BDF0;
    const double bdf1 = rData.BDF1;
    const double bdf2 = rData.BDF2;
    const double dt = rData.DeltaTime;

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