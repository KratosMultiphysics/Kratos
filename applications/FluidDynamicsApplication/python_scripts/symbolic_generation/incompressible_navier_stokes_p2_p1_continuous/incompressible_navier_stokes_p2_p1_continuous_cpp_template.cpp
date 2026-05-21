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
    //substitute_lhs_2D
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
    //substitute_rhs_2D
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