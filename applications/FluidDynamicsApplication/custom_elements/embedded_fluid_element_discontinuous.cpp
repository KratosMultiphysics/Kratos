#include "includes/kratos_flags.h"

#include "custom_elements/embedded_fluid_element_discontinuous.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/symbolic_navier_stokes.h"

#include "custom_utilities/element_size_calculator.h"
#include "custom_utilities/embedded_discontinuous_data.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"

#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"


namespace Kratos {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TBaseElement >
EmbeddedFluidElementDiscontinuous<TBaseElement>::EmbeddedFluidElementDiscontinuous(IndexType NewId):
    TBaseElement(NewId)
{}

template< class TBaseElement >
EmbeddedFluidElementDiscontinuous<TBaseElement>::EmbeddedFluidElementDiscontinuous(IndexType NewId, const NodesArrayType& ThisNodes):
    TBaseElement(NewId,ThisNodes)
{}


template< class TBaseElement >
EmbeddedFluidElementDiscontinuous<TBaseElement>::EmbeddedFluidElementDiscontinuous(IndexType NewId, Geometry<NodeType>::Pointer pGeometry):
    TBaseElement(NewId,pGeometry)
{}

template< class TBaseElement >
EmbeddedFluidElementDiscontinuous<TBaseElement>::EmbeddedFluidElementDiscontinuous(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties):
    TBaseElement(NewId,pGeometry,pProperties)
{}


template< class TBaseElement >
EmbeddedFluidElementDiscontinuous<TBaseElement>::~EmbeddedFluidElementDiscontinuous()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TBaseElement >
Element::Pointer EmbeddedFluidElementDiscontinuous<TBaseElement>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_shared<EmbeddedFluidElementDiscontinuous>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TBaseElement >
Element::Pointer EmbeddedFluidElementDiscontinuous<TBaseElement>::Create(
    IndexType NewId,
    Geometry<NodeType>::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_shared<EmbeddedFluidElementDiscontinuous>(NewId, pGeom, pProperties);
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize){
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    if (rRightHandSideVector.size() != LocalSize){
        rRightHandSideVector.resize(LocalSize, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    EmbeddedDiscontinuousElementData data;
    data.Initialize(*this, rCurrentProcessInfo);
    this->InitializeGeometryData(data);

    // Iterate over the positive side volume integration points
    const unsigned int number_of_positive_gauss_points = data.PositiveSideWeights.size();
    for (unsigned int g = 0; g < number_of_positive_gauss_points; ++g){
        const size_t gauss_pt_index = g;
        data.UpdateGeometryValues(gauss_pt_index, data.PositiveSideWeights[g], row(data.PositiveSideN, g), data.PositiveSideDNDX[g]);
        this->CalculateMaterialResponse(data);
        this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
    }

    // Iterate over the negative side volume integration points
    const unsigned int number_of_negative_gauss_points = data.NegativeSideWeights.size();
    for (unsigned int g = 0; g < number_of_negative_gauss_points; ++g){
        const size_t gauss_pt_index = g + number_of_negative_gauss_points;
        data.UpdateGeometryValues(gauss_pt_index, data.NegativeSideWeights[g], row(data.NegativeSideN, g), data.NegativeSideDNDX[g]);
        this->CalculateMaterialResponse(data);
        this->AddTimeIntegratedSystem(data, rLeftHandSideMatrix, rRightHandSideVector);
    }

    // If the element is cut, add the interface contributions
    if ( data.IsCut() ) {

        // Add the boundary term together with the interface equilibrium imposition. Note that the interface 
        // equilibrium imposition and boundary term addition yields minus the base element boundary term. 
        // Therefore, two auxiliar arrays are used (aux_LHS and aux_RHS) to store the base element boundary 
        // term contribution. Finally, the opposite of these arrays is added to the local system.
        const size_t volume_gauss_points = number_of_positive_gauss_points + number_of_negative_gauss_points;
        VectorType aux_RHS = ZeroVector(LocalSize); 
        MatrixType aux_LHS = ZeroMatrix(LocalSize, LocalSize);

        // Add the base element boundary contribution on the positive interface
        const unsigned int number_of_positive_interface_gauss_points = data.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < number_of_positive_interface_gauss_points; ++g){
            const size_t gauss_pt_index = g + volume_gauss_points;
            data.UpdateGeometryValues(gauss_pt_index, data.PositiveInterfaceWeights[g], row(data.PositiveInterfaceN, g), data.PositiveInterfaceDNDX[g]);
            this->CalculateMaterialResponse(data);
            this-> AddBoundaryTraction(data, data.PositiveInterfaceUnitNormals[g], aux_LHS, aux_RHS);
        }

        // Add the base element boundary contribution on the negative interface
        const unsigned int number_of_negative_interface_gauss_points = data.NegativeInterfaceWeights.size();
        for (unsigned int g = 0; g < number_of_negative_interface_gauss_points; ++g){
            const size_t gauss_pt_index = g + volume_gauss_points + number_of_positive_interface_gauss_points;
            data.UpdateGeometryValues(gauss_pt_index, data.NegativeInterfaceWeights[g], row(data.NegativeInterfaceN, g), data.NegativeInterfaceDNDX[g]);
            this->CalculateMaterialResponse(data);
            this-> AddBoundaryTraction(data, data.NegativeInterfaceUnitNormals[g], aux_LHS, aux_RHS);
        }

        // Recall to swap the boundary term signs because of the interface equilibrium (Neumann) imposition
        rLeftHandSideMatrix -= aux_LHS;
        rRightHandSideVector -= aux_RHS;

        // Add the Nitsche Navier boundary condition implementation (Winter, 2018)
        data.InitializeBoundaryConditionData(rCurrentProcessInfo);
        AddNormalPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
        // AddNormalSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, data); // NOTE: IMPLEMENT THE SKEW-SYMMETRIC ADJOINT IF IT IS NEEDED IN THE FUTURE. CREATE A IS_SKEW_SYMMETRIC ELEMENTAL FLAG.
        // AddTangentialPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
        // AddTangentialSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, data); // NOTE: IMPLEMENT THE SKEW-SYMMETRIC ADJOINT IF IT IS NEEDED IN THE FUTURE. CREATE A IS_SKEW_SYMMETRIC ELEMENTAL FLAG.

        // if (this->Id() == 628) {
        //     KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@")
        //     KRATOS_WATCH("BOUNDARY CONDITION CONTRIBUTION")
        //     KRATOS_WATCH("@@@@@@@@@@@@@@@@@@@@@@")
        //     for (unsigned int i = 0; i < rLeftHandSideMatrix.size1(); ++i) {
        //         std::cout << "[" << i << "] ";
        //         for (unsigned int j = 0; j < rLeftHandSideMatrix.size2(); ++j) {
        //             std::cout << rLeftHandSideMatrix(i,j) << " ";
        //         }
        //         std::cout << std::endl;
        //     }
        //     KRATOS_WATCH(rRightHandSideVector)
        // }
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::Calculate(
    const Variable<double> &rVariable,
    double& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = 0.0;
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);

    // If the element is split, integrate sigma.n over the interface
    // Note that in the ausas formulation, both interface sides need to be integrated
    if (rVariable == DRAG_FORCE) {

        EmbeddedDiscontinuousElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        const unsigned int number_of_positive_gauss_points = data.PositiveSideWeights.size();

        if ( data.IsCut() ){
            // Integrate positive interface side drag
            const unsigned int n_int_pos_gauss = data.PositiveInterfaceWeights.size();
            for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {

                // Update the Gauss pt. data
                data.UpdateGeometryValues(g + number_of_positive_gauss_points,data.PositiveInterfaceWeights[g],row(data.PositiveInterfaceN, g),data.PositiveInterfaceDNDX[g]);

                // Get the interface Gauss pt. unit noromal
                const auto &aux_unit_normal = data.PositiveInterfaceUnitNormals[g];

                // Compute Gauss pt. pressure
                const double p_gauss = inner_prod(data.N, data.Pressure);

                // Call the constitutive law to compute the shear contribution
                this->CalculateMaterialResponse(data);

                // Get the normal projection matrix in Voigt notation
                BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
                FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

                // Add the shear and pressure drag contributions
                const array_1d<double, Dim> shear_proj = data.Weight * prod(voigt_normal_proj_matrix, data.ShearStress);
                for (unsigned int i = 0; i < Dim ; ++i){
                    rOutput(i) -= shear_proj(i);
                }
                rOutput += data.Weight * p_gauss * aux_unit_normal;
            }

            // Integrate negative interface side drag
            const unsigned int n_int_neg_gauss = data.NegativeInterfaceWeights.size();
            for (unsigned int g = 0; g < n_int_neg_gauss; ++g) {

                // Update the Gauss pt. data
                data.UpdateGeometryValues(g + number_of_positive_gauss_points,data.NegativeInterfaceWeights[g],row(data.NegativeInterfaceN, g),data.NegativeInterfaceDNDX[g]);

                // Get the interface Gauss pt. unit noromal
                const auto &aux_unit_normal = data.NegativeInterfaceUnitNormals[g];

                // Compute Gauss pt. pressure
                const double p_gauss = inner_prod(data.N, data.Pressure);

                // Call the constitutive law to compute the shear contribution
                this->CalculateMaterialResponse(data);

                // Get the normal projection matrix in Voigt notation
                BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
                FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

                // Add the shear and pressure drag contributions
                const array_1d<double, Dim> shear_proj = data.Weight * prod(voigt_normal_proj_matrix, data.ShearStress);
                for (unsigned int i = 0; i < Dim ; ++i){
                    rOutput(i) -= shear_proj(i);
                }
                rOutput += data.Weight * p_gauss * aux_unit_normal;
            }
        }

    } else {
        TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::Calculate(
    const Variable<Vector> &rVariable,
    Vector& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::Calculate(
    const Variable<Matrix> &rVariable,
    Matrix& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <class TBaseElement>
int EmbeddedFluidElementDiscontinuous<TBaseElement>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    int out = EmbeddedDiscontinuousElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    return TBaseElement::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TBaseElement>
std::string EmbeddedFluidElementDiscontinuous<TBaseElement>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedFluidElementDiscontinuous #" << this->Id();
    return buffer.str();
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedFluidElementDiscontinuous" << Dim << "D" << NumNodes << "N"
             << std::endl
             << "on top of ";
    TBaseElement::PrintInfo(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// Operations

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::InitializeGeometryData(EmbeddedDiscontinuousElementData& rData) const
{
    rData.PositiveIndices.clear();
    rData.NegativeIndices.clear();

    // Number of positive and negative distance function values
    for (size_t i = 0; i < EmbeddedDiscontinuousElementData::NumNodes; ++i){
        if (rData.ElementalDistances[i] > 0.0){
            rData.NumPositiveNodes++;
            rData.PositiveIndices.push_back(i);
        } else {
            rData.NumNegativeNodes++;
            rData.NegativeIndices.push_back(i);
        }
    }

    if (rData.IsCut()){
        this->DefineCutGeometryData(rData);
    } else {
        this->DefineStandardGeometryData(rData);
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::DefineStandardGeometryData(EmbeddedDiscontinuousElementData& rData) const
{
    rData.NumNegativeNodes = 0;
    rData.NumPositiveNodes = NumNodes;
    this->CalculateGeometryData(rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::DefineCutGeometryData(EmbeddedDiscontinuousElementData& rData) const
{
    // Auxiliary distance vector for the element subdivision utility
    Vector elemental_distances = rData.ElementalDistances;

    ModifiedShapeFunctions::Pointer p_calculator =
        EmbeddedDiscontinuousInternals::GetShapeFunctionCalculator<EmbeddedDiscontinuousElementData::Dim, EmbeddedDiscontinuousElementData::NumNodes>(
            *this, 
            elemental_distances);

    // Positive side volume
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN, 
        rData.PositiveSideDNDX, 
        rData.PositiveSideWeights,
        GeometryData::GI_GAUSS_2);

    // Negative side volume
    p_calculator->ComputeNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeSideN, 
        rData.NegativeSideDNDX, 
        rData.NegativeSideWeights,
        GeometryData::GI_GAUSS_2);

    // Positive side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN,
        rData.PositiveInterfaceDNDX,
        rData.PositiveInterfaceWeights,
        GeometryData::GI_GAUSS_2);

    // Negative side interface
    p_calculator->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
        rData.NegativeInterfaceN,
        rData.NegativeInterfaceDNDX,
        rData.NegativeInterfaceWeights,
        GeometryData::GI_GAUSS_2);

    // Positive side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals,
        GeometryData::GI_GAUSS_2);

    // Negative side interface normals
    p_calculator->ComputeNegativeSideInterfaceAreaNormals(
        rData.NegativeInterfaceUnitNormals,
        GeometryData::GI_GAUSS_2);

    // Normalize the normals
    // Note: we calculate h here (and we don't use the value in rData.ElementSize)
    // because rData.ElementSize might still be uninitialized: some data classes define it at the Gauss point.
    double h = ElementSizeCalculator<Dim,NumNodes>::MinimumElementSize(this->GetGeometry());
    const double tolerance = std::pow(1e-3 * h, Dim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
    this->NormalizeInterfaceNormals(rData.NegativeInterfaceUnitNormals, tolerance);

    // Avoid numerical 0 in constant shape function gradients
    for (unsigned int i_gauss = 0; i_gauss < rData.PositiveSideDNDX.size(); ++i_gauss) {
        if (norm_frobenius(rData.PositiveSideDNDX[i_gauss]) < 1e-15) {
            rData.PositiveSideDNDX[i_gauss] = ZeroMatrix(rData.PositiveSideDNDX[i_gauss].size1(), rData.PositiveSideDNDX[i_gauss].size2());
        }
    }
    for (unsigned int i_gauss = 0; i_gauss < rData.NegativeSideDNDX.size(); ++i_gauss) {
        if (norm_frobenius(rData.NegativeSideDNDX[i_gauss]) < 1e-15) {
            rData.NegativeSideDNDX[i_gauss] = ZeroMatrix(rData.NegativeSideDNDX[i_gauss].size1(), rData.NegativeSideDNDX[i_gauss].size2());
        }
    }
    for (unsigned int i_gauss = 0; i_gauss < rData.PositiveInterfaceDNDX.size(); ++i_gauss) {
        if (norm_frobenius(rData.PositiveInterfaceDNDX[i_gauss]) < 1e-15) {
            rData.PositiveInterfaceDNDX[i_gauss] = ZeroMatrix(rData.PositiveInterfaceDNDX[i_gauss].size1(), rData.PositiveInterfaceDNDX[i_gauss].size2());
        }
    }
    for (unsigned int i_gauss = 0; i_gauss < rData.NegativeInterfaceDNDX.size(); ++i_gauss) {
        if (norm_frobenius(rData.NegativeInterfaceDNDX[i_gauss]) < 1e-15) {
            rData.NegativeInterfaceDNDX[i_gauss] = ZeroMatrix(rData.NegativeInterfaceDNDX[i_gauss].size1(), rData.NegativeInterfaceDNDX[i_gauss].size2());
        }
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::NormalizeInterfaceNormals(
    typename EmbeddedDiscontinuousElementData::InterfaceNormalsType& rNormals,
    double Tolerance) const
{
    for (unsigned int i = 0; i < rNormals.size(); ++i) {
        double norm = norm_2(rNormals[i]);
        rNormals[i] /= std::max(norm,Tolerance);
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::AddNormalPenaltyContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    EmbeddedDiscontinuousElementData& rData)
{
    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);

    // Compute the positive side LHS and RHS contributions
    const unsigned int number_of_positive_interface_integration_points = rData.PositiveInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_positive_interface_integration_points; ++g) {
        // Get the Gauss pt. data
        const double weight = rData.PositiveInterfaceWeights[g];
        rData.N = row(rData.PositiveInterfaceN, g);
        rData.DN_DX = rData.PositiveInterfaceDNDX[g];
        rData.Normal = rData.PositiveInterfaceUnitNormals[g];

        // Compute the Gauss pt. LHS contribution
        MatrixType aux_lhs = NitscheTermsLHS3D(rData);
        VectorType aux_rhs = NitscheTermsRHS3D(rData);
        rLHS += weight * aux_lhs;
        rRHS += weight * aux_rhs;
    }

    // Compute the negative side LHS and RHS contributions
    const unsigned int number_of_negative_interface_integration_points = rData.NegativeInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_negative_interface_integration_points; ++g) {
        // Get the Gauss pt. data
        const double weight = rData.NegativeInterfaceWeights[g];
        rData.N = row(rData.NegativeInterfaceN, g);
        rData.DN_DX = rData.NegativeInterfaceDNDX[g];
        rData.Normal = rData.NegativeInterfaceUnitNormals[g];

        // Compute the Gauss pt. LHS contribution
        MatrixType aux_lhs = NitscheTermsLHS3D(rData);
        VectorType aux_rhs = NitscheTermsRHS3D(rData);
        rLHS += weight * aux_lhs;
        rRHS += weight * aux_rhs;
    }
}

template <class TBaseElement>
typename EmbeddedFluidElementDiscontinuous<TBaseElement>::MatrixType EmbeddedFluidElementDiscontinuous<TBaseElement>::NitscheTermsLHS2D(
    EmbeddedDiscontinuousElementData& rData)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // If there is embedded velocity, substract it to the previous iteration solution
    array_1d<double,3> g = ZeroVector(3);
    if (this->Has(EMBEDDED_VELOCITY)){
        g = this->GetValue(EMBEDDED_VELOCITY);
    }

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;
    const auto& normal = rData.Normal;

    // Get Nitsche imposition values
    const double eps = rData.SlipLength;;
    const double gamma = rData.PenaltyCoefficient;
    const double adjoint = 1.0;

    MatrixType lhs = ZeroMatrix(9,9);

    const double clhs0 =             pow(normal[0], 2);
const double clhs1 =             pow(N[0], 2);
const double clhs2 =             1.0/gamma;
const double clhs3 =             1.0/h;
const double clhs4 =             clhs1*clhs2*clhs3*mu;
const double clhs5 =             clhs0 - 1;
const double clhs6 =             1.0/(eps + gamma*h);
const double clhs7 =             clhs1*clhs6*mu;
const double clhs8 =             N[0]*clhs6*mu;
const double clhs9 =             1.0*DN(0,1)*adjoint*clhs0*gamma*h*normal[1];
const double clhs10 =             N[0]*adjoint;
const double clhs11 =             normal[0]*normal[1];
const double clhs12 =             C(0,2)*DN(0,0);
const double clhs13 =             C(2,2)*DN(0,1) + clhs12;
const double clhs14 =             C(1,2)*DN(0,1);
const double clhs15 =             clhs13*normal[0] + normal[1]*(C(0,1)*DN(0,0) + clhs14);
const double clhs16 =             clhs11*clhs15;
const double clhs17 =             N[0]*adjoint*clhs0;
const double clhs18 =             clhs13*normal[1] + normal[0]*(C(0,0)*DN(0,0) + C(0,2)*DN(0,1));
const double clhs19 =             DN(0,0)*normal[0];
const double clhs20 =             DN(0,1)*normal[1];
const double clhs21 =             2.0*clhs19 + 1.0*clhs20;
const double clhs22 =             N[0]*adjoint*clhs5*clhs6*gamma*h*mu;
const double clhs23 =             h*rho*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2)) + mu + pow(h, 2)*rho/dt;
const double clhs24 =             clhs1*clhs2*clhs23*clhs3;
const double clhs25 =             N[0]*clhs6*eps;
const double clhs26 =             clhs16 + clhs18*clhs5;
const double clhs27 =             1.0*DN(0,1)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs28 =             clhs11*clhs18;
const double clhs29 =             pow(normal[1], 2);
const double clhs30 =             clhs29 - 1;
const double clhs31 =             clhs15*clhs30 + clhs28;
const double clhs32 =             adjoint*clhs21*clhs6*eps*gamma*h;
const double clhs33 =             clhs2*clhs3*mu*normal[0]*normal[1];
const double clhs34 =             clhs6*mu*normal[0]*normal[1];
const double clhs35 =             clhs2*clhs23*clhs3*normal[0]*normal[1];
const double clhs36 =             clhs1*clhs33 - clhs1*clhs34 + clhs1*clhs35;
const double clhs37 =             1.0*DN(0,1)*adjoint*clhs30*gamma*h*normal[0];
const double clhs38 =             N[0]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[1];
const double clhs39 =             N[0]*adjoint*clhs29;
const double clhs40 =             C(2,2)*DN(0,0) + clhs14;
const double clhs41 =             clhs40*normal[0] + normal[1]*(C(1,1)*DN(0,1) + C(1,2)*DN(0,0));
const double clhs42 =             clhs11*clhs41;
const double clhs43 =             -clhs0 + 1;
const double clhs44 =             clhs40*normal[1] + normal[0]*(C(0,1)*DN(0,1) + clhs12);
const double clhs45 =             clhs42 - clhs43*clhs44;
const double clhs46 =             clhs11*clhs44;
const double clhs47 =             -clhs29 + 1;
const double clhs48 =             -clhs41*clhs47 + clhs46;
const double clhs49 =             N[0]*N[1]*clhs2*clhs3*mu;
const double clhs50 =             N[0]*N[1]*clhs6*mu;
const double clhs51 =             N[0]*N[1]*clhs2*clhs23*clhs3;
const double clhs52 =             clhs0*clhs49 + clhs0*clhs51 - clhs5*clhs50;
const double clhs53 =             N[1]*clhs6*mu;
const double clhs54 =             N[1]*adjoint;
const double clhs55 =             N[1]*adjoint*clhs0;
const double clhs56 =             N[1]*adjoint*clhs5*clhs6*gamma*h*mu;
const double clhs57 =             C(0,2)*DN(1,0);
const double clhs58 =             C(2,2)*DN(1,1) + clhs57;
const double clhs59 =             C(1,2)*DN(1,1);
const double clhs60 =             clhs58*normal[0] + normal[1]*(C(0,1)*DN(1,0) + clhs59);
const double clhs61 =             clhs11*clhs60;
const double clhs62 =             clhs58*normal[1] + normal[0]*(C(0,0)*DN(1,0) + C(0,2)*DN(1,1));
const double clhs63 =             clhs5*clhs62 + clhs61;
const double clhs64 =             clhs11*clhs62;
const double clhs65 =             clhs30*clhs60 + clhs64;
const double clhs66 =             N[0]*clhs2*clhs3*mu*normal[0]*normal[1];
const double clhs67 =             N[1]*clhs66 - clhs11*clhs50 + clhs11*clhs51;
const double clhs68 =             N[1]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[1];
const double clhs69 =             N[1]*adjoint*clhs29;
const double clhs70 =             C(2,2)*DN(1,0) + clhs59;
const double clhs71 =             clhs70*normal[0] + normal[1]*(C(1,1)*DN(1,1) + C(1,2)*DN(1,0));
const double clhs72 =             clhs11*clhs71;
const double clhs73 =             clhs70*normal[1] + normal[0]*(C(0,1)*DN(1,1) + clhs57);
const double clhs74 =             -clhs43*clhs73 + clhs72;
const double clhs75 =             clhs11*clhs73;
const double clhs76 =             -clhs47*clhs71 + clhs75;
const double clhs77 =             N[0]*N[2]*clhs2*clhs3*mu;
const double clhs78 =             N[0]*N[2]*clhs6*mu;
const double clhs79 =             N[0]*N[2]*clhs2*clhs23*clhs3;
const double clhs80 =             clhs0*clhs77 + clhs0*clhs79 - clhs5*clhs78;
const double clhs81 =             N[2]*clhs6*mu;
const double clhs82 =             N[2]*adjoint*normal[0]*normal[1];
const double clhs83 =             N[2]*adjoint*clhs0;
const double clhs84 =             N[2]*adjoint*clhs5*clhs6*gamma*h*mu;
const double clhs85 =             C(0,2)*DN(2,0);
const double clhs86 =             C(2,2)*DN(2,1) + clhs85;
const double clhs87 =             C(1,2)*DN(2,1);
const double clhs88 =             clhs86*normal[0] + normal[1]*(C(0,1)*DN(2,0) + clhs87);
const double clhs89 =             clhs11*clhs88;
const double clhs90 =             clhs86*normal[1] + normal[0]*(C(0,0)*DN(2,0) + C(0,2)*DN(2,1));
const double clhs91 =             clhs5*clhs90 + clhs89;
const double clhs92 =             clhs11*clhs90;
const double clhs93 =             clhs30*clhs88 + clhs92;
const double clhs94 =             N[2]*normal[0]*normal[1];
const double clhs95 =             N[0]*clhs2*clhs23*clhs3*clhs94 + N[2]*clhs66 - clhs8*clhs94;
const double clhs96 =             N[2]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[1];
const double clhs97 =             N[2]*adjoint*clhs29;
const double clhs98 =             C(2,2)*DN(2,0) + clhs87;
const double clhs99 =             clhs98*normal[0] + normal[1]*(C(1,1)*DN(2,1) + C(1,2)*DN(2,0));
const double clhs100 =             clhs11*clhs99;
const double clhs101 =             clhs98*normal[1] + normal[0]*(C(0,1)*DN(2,1) + clhs85);
const double clhs102 =             clhs100 - clhs101*clhs43;
const double clhs103 =             clhs101*clhs11;
const double clhs104 =             clhs103 - clhs47*clhs99;
const double clhs105 =             1.0*DN(0,0)*adjoint*clhs5*gamma*h*normal[1];
const double clhs106 =             1.0*clhs19 + 2.0*clhs20;
const double clhs107 =             -clhs15*clhs47 + clhs28;
const double clhs108 =             1.0*DN(0,0)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs109 =             clhs16 - clhs18*clhs43;
const double clhs110 =             adjoint*clhs106*clhs6*eps*gamma*h;
const double clhs111 =             1.0*DN(0,0)*adjoint*clhs29*gamma*h*normal[0];
const double clhs112 =             N[0]*adjoint*clhs30*clhs6*gamma*h*mu;
const double clhs113 =             clhs30*clhs41 + clhs46;
const double clhs114 =             clhs42 + clhs44*clhs5;
const double clhs115 =             -clhs47*clhs60 + clhs64;
const double clhs116 =             -clhs43*clhs62 + clhs61;
const double clhs117 =             clhs29*clhs49 + clhs29*clhs51 - clhs30*clhs50;
const double clhs118 =             N[1]*adjoint*clhs30*clhs6*gamma*h*mu;
const double clhs119 =             clhs30*clhs71 + clhs75;
const double clhs120 =             clhs5*clhs73 + clhs72;
const double clhs121 =             -clhs47*clhs88 + clhs92;
const double clhs122 =             -clhs43*clhs90 + clhs89;
const double clhs123 =             clhs29*clhs77 + clhs29*clhs79 - clhs30*clhs78;
const double clhs124 =             N[2]*adjoint*clhs30*clhs6*gamma*h*mu;
const double clhs125 =             clhs103 + clhs30*clhs99;
const double clhs126 =             clhs100 + clhs101*clhs5;
const double clhs127 =             clhs0 + clhs29;
const double clhs128 =             clhs1*clhs127;
const double clhs129 =             N[0]*N[1]*clhs127;
const double clhs130 =             -clhs129*normal[0];
const double clhs131 =             -clhs129*normal[1];
const double clhs132 =             N[0]*N[2]*clhs127;
const double clhs133 =             -clhs132*normal[0];
const double clhs134 =             -clhs132*normal[1];
const double clhs135 =             1.0*DN(1,1)*adjoint*clhs0*gamma*h*normal[1];
const double clhs136 =             DN(1,0)*normal[0];
const double clhs137 =             DN(1,1)*normal[1];
const double clhs138 =             2.0*clhs136 + 1.0*clhs137;
const double clhs139 =             N[1]*clhs6*eps;
const double clhs140 =             1.0*DN(1,1)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs141 =             adjoint*clhs138*clhs6*eps*gamma*h;
const double clhs142 =             1.0*DN(1,1)*adjoint*clhs30*gamma*h*normal[0];
const double clhs143 =             pow(N[1], 2);
const double clhs144 =             clhs143*clhs2*clhs3*mu;
const double clhs145 =             clhs143*clhs6*mu;
const double clhs146 =             clhs143*clhs2*clhs23*clhs3;
const double clhs147 =             clhs143*clhs33 - clhs143*clhs34 + clhs143*clhs35;
const double clhs148 =             N[1]*N[2]*clhs2*clhs3*mu;
const double clhs149 =             N[1]*N[2]*clhs6*mu;
const double clhs150 =             N[1]*N[2]*clhs2*clhs23*clhs3;
const double clhs151 =             clhs0*clhs148 + clhs0*clhs150 - clhs149*clhs5;
const double clhs152 =             N[1]*clhs2*clhs23*clhs3*clhs94 + N[1]*clhs2*clhs3*clhs94*mu - clhs53*clhs94;
const double clhs153 =             1.0*DN(1,0)*adjoint*clhs5*gamma*h*normal[1];
const double clhs154 =             1.0*clhs136 + 2.0*clhs137;
const double clhs155 =             1.0*DN(1,0)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs156 =             adjoint*clhs154*clhs6*eps*gamma*h;
const double clhs157 =             1.0*DN(1,0)*adjoint*clhs29*gamma*h*normal[0];
const double clhs158 =             clhs148*clhs29 - clhs149*clhs30 + clhs150*clhs29;
const double clhs159 =             clhs127*clhs143;
const double clhs160 =             N[1]*N[2]*clhs127;
const double clhs161 =             -clhs160*normal[0];
const double clhs162 =             -clhs160*normal[1];
const double clhs163 =             1.0*DN(2,1)*adjoint*clhs0*gamma*h*normal[1];
const double clhs164 =             DN(2,0)*normal[0];
const double clhs165 =             DN(2,1)*normal[1];
const double clhs166 =             2.0*clhs164 + 1.0*clhs165;
const double clhs167 =             N[2]*clhs6*eps;
const double clhs168 =             1.0*DN(2,1)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs169 =             adjoint*clhs166*clhs6*eps*gamma*h;
const double clhs170 =             1.0*DN(2,1)*adjoint*clhs30*gamma*h*normal[0];
const double clhs171 =             pow(N[2], 2);
const double clhs172 =             clhs171*clhs2*clhs3*mu;
const double clhs173 =             clhs171*clhs6*mu;
const double clhs174 =             clhs171*clhs2*clhs23*clhs3;
const double clhs175 =             clhs171*clhs33 - clhs171*clhs34 + clhs171*clhs35;
const double clhs176 =             1.0*DN(2,0)*adjoint*clhs5*gamma*h*normal[1];
const double clhs177 =             1.0*clhs164 + 2.0*clhs165;
const double clhs178 =             1.0*DN(2,0)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs179 =             adjoint*clhs177*clhs6*eps*gamma*h;
const double clhs180 =             1.0*DN(2,0)*adjoint*clhs29*gamma*h*normal[0];
const double clhs181 =             clhs127*clhs171;
            lhs(0,0)=clhs0*clhs24 + clhs0*clhs4 - clhs10*clhs16 - clhs17*clhs18 + clhs21*clhs22 - clhs25*clhs26 + clhs26*clhs32 + clhs27*clhs31 - clhs5*clhs7 + clhs8*clhs9;
            lhs(0,1)=-clhs10*clhs28 - clhs15*clhs39 + clhs21*clhs38 - clhs25*clhs45 + clhs27*clhs48 + clhs32*clhs45 + clhs36 + clhs37*clhs8;
            lhs(0,2)=0;
            lhs(0,3)=-clhs16*clhs54 - clhs18*clhs55 + clhs21*clhs56 - clhs25*clhs63 + clhs27*clhs65 + clhs32*clhs63 + clhs52 + clhs53*clhs9;
            lhs(0,4)=-clhs15*clhs69 + clhs21*clhs68 - clhs25*clhs74 + clhs27*clhs76 - clhs28*clhs54 + clhs32*clhs74 + clhs37*clhs53 + clhs67;
            lhs(0,5)=0;
            lhs(0,6)=-clhs15*clhs82 - clhs18*clhs83 + clhs21*clhs84 - clhs25*clhs91 + clhs27*clhs93 + clhs32*clhs91 + clhs80 + clhs81*clhs9;
            lhs(0,7)=-clhs102*clhs25 + clhs102*clhs32 + clhs104*clhs27 - clhs15*clhs97 - clhs18*clhs82 + clhs21*clhs96 + clhs37*clhs81 + clhs95;
            lhs(0,8)=0;
            lhs(1,0)=-clhs10*clhs42 + clhs105*clhs8 + clhs106*clhs38 + clhs107*clhs110 - clhs107*clhs25 + clhs108*clhs109 - clhs17*clhs44 + clhs36;
            lhs(1,1)=-clhs10*clhs46 + clhs106*clhs112 + clhs108*clhs114 + clhs110*clhs113 + clhs111*clhs8 - clhs113*clhs25 + clhs24*clhs29 + clhs29*clhs4 - clhs30*clhs7 - clhs39*clhs41;
            lhs(1,2)=0;
            lhs(1,3)=clhs105*clhs53 + clhs106*clhs68 + clhs108*clhs116 + clhs110*clhs115 - clhs115*clhs25 - clhs42*clhs54 - clhs44*clhs55 + clhs67;
            lhs(1,4)=clhs106*clhs118 + clhs108*clhs120 + clhs110*clhs119 + clhs111*clhs53 + clhs117 - clhs119*clhs25 - clhs41*clhs69 - clhs46*clhs54;
            lhs(1,5)=0;
            lhs(1,6)=clhs105*clhs81 + clhs106*clhs96 + clhs108*clhs122 + clhs110*clhs121 - clhs121*clhs25 - clhs41*clhs82 - clhs44*clhs83 + clhs95;
            lhs(1,7)=clhs106*clhs124 + clhs108*clhs126 + clhs110*clhs125 + clhs111*clhs81 + clhs123 - clhs125*clhs25 - clhs41*clhs97 - clhs44*clhs82;
            lhs(1,8)=0;
            lhs(2,0)=-clhs128*normal[0];
            lhs(2,1)=-clhs128*normal[1];
            lhs(2,2)=0;
            lhs(2,3)=clhs130;
            lhs(2,4)=clhs131;
            lhs(2,5)=0;
            lhs(2,6)=clhs133;
            lhs(2,7)=clhs134;
            lhs(2,8)=0;
            lhs(3,0)=-clhs10*clhs61 + clhs135*clhs8 + clhs138*clhs22 - clhs139*clhs26 + clhs140*clhs31 + clhs141*clhs26 - clhs17*clhs62 + clhs52;
            lhs(3,1)=-clhs10*clhs64 + clhs138*clhs38 - clhs139*clhs45 + clhs140*clhs48 + clhs141*clhs45 + clhs142*clhs8 - clhs39*clhs60 + clhs67;
            lhs(3,2)=0;
            lhs(3,3)=clhs0*clhs144 + clhs0*clhs146 + clhs135*clhs53 + clhs138*clhs56 - clhs139*clhs63 + clhs140*clhs65 + clhs141*clhs63 - clhs145*clhs5 - clhs54*clhs61 - clhs55*clhs62;
            lhs(3,4)=clhs138*clhs68 - clhs139*clhs74 + clhs140*clhs76 + clhs141*clhs74 + clhs142*clhs53 + clhs147 - clhs54*clhs64 - clhs60*clhs69;
            lhs(3,5)=0;
            lhs(3,6)=clhs135*clhs81 + clhs138*clhs84 - clhs139*clhs91 + clhs140*clhs93 + clhs141*clhs91 + clhs151 - clhs60*clhs82 - clhs62*clhs83;
            lhs(3,7)=-clhs102*clhs139 + clhs102*clhs141 + clhs104*clhs140 + clhs138*clhs96 + clhs142*clhs81 + clhs152 - clhs60*clhs97 - clhs62*clhs82;
            lhs(3,8)=0;
            lhs(4,0)=-clhs10*clhs72 - clhs107*clhs139 + clhs107*clhs156 + clhs109*clhs155 + clhs153*clhs8 + clhs154*clhs38 - clhs17*clhs73 + clhs67;
            lhs(4,1)=-clhs10*clhs75 + clhs112*clhs154 - clhs113*clhs139 + clhs113*clhs156 + clhs114*clhs155 + clhs117 + clhs157*clhs8 - clhs39*clhs71;
            lhs(4,2)=0;
            lhs(4,3)=-clhs115*clhs139 + clhs115*clhs156 + clhs116*clhs155 + clhs147 + clhs153*clhs53 + clhs154*clhs68 - clhs54*clhs72 - clhs55*clhs73;
            lhs(4,4)=clhs118*clhs154 - clhs119*clhs139 + clhs119*clhs156 + clhs120*clhs155 + clhs144*clhs29 - clhs145*clhs30 + clhs146*clhs29 + clhs157*clhs53 - clhs54*clhs75 - clhs69*clhs71;
            lhs(4,5)=0;
            lhs(4,6)=-clhs121*clhs139 + clhs121*clhs156 + clhs122*clhs155 + clhs152 + clhs153*clhs81 + clhs154*clhs96 - clhs71*clhs82 - clhs73*clhs83;
            lhs(4,7)=clhs124*clhs154 - clhs125*clhs139 + clhs125*clhs156 + clhs126*clhs155 + clhs157*clhs81 + clhs158 - clhs71*clhs97 - clhs73*clhs82;
            lhs(4,8)=0;
            lhs(5,0)=clhs130;
            lhs(5,1)=clhs131;
            lhs(5,2)=0;
            lhs(5,3)=-clhs159*normal[0];
            lhs(5,4)=-clhs159*normal[1];
            lhs(5,5)=0;
            lhs(5,6)=clhs161;
            lhs(5,7)=clhs162;
            lhs(5,8)=0;
            lhs(6,0)=-clhs10*clhs89 + clhs163*clhs8 + clhs166*clhs22 - clhs167*clhs26 + clhs168*clhs31 + clhs169*clhs26 - clhs17*clhs90 + clhs80;
            lhs(6,1)=-clhs10*clhs92 + clhs166*clhs38 - clhs167*clhs45 + clhs168*clhs48 + clhs169*clhs45 + clhs170*clhs8 - clhs39*clhs88 + clhs95;
            lhs(6,2)=0;
            lhs(6,3)=clhs151 + clhs163*clhs53 + clhs166*clhs56 - clhs167*clhs63 + clhs168*clhs65 + clhs169*clhs63 - clhs54*clhs89 - clhs55*clhs90;
            lhs(6,4)=clhs152 + clhs166*clhs68 - clhs167*clhs74 + clhs168*clhs76 + clhs169*clhs74 + clhs170*clhs53 - clhs54*clhs92 - clhs69*clhs88;
            lhs(6,5)=0;
            lhs(6,6)=clhs0*clhs172 + clhs0*clhs174 + clhs163*clhs81 + clhs166*clhs84 - clhs167*clhs91 + clhs168*clhs93 + clhs169*clhs91 - clhs173*clhs5 - clhs82*clhs88 - clhs83*clhs90;
            lhs(6,7)=-clhs102*clhs167 + clhs102*clhs169 + clhs104*clhs168 + clhs166*clhs96 + clhs170*clhs81 + clhs175 - clhs82*clhs90 - clhs88*clhs97;
            lhs(6,8)=0;
            lhs(7,0)=-clhs10*clhs100 - clhs101*clhs17 - clhs107*clhs167 + clhs107*clhs179 + clhs109*clhs178 + clhs176*clhs8 + clhs177*clhs38 + clhs95;
            lhs(7,1)=-clhs10*clhs103 + clhs112*clhs177 - clhs113*clhs167 + clhs113*clhs179 + clhs114*clhs178 + clhs123 + clhs180*clhs8 - clhs39*clhs99;
            lhs(7,2)=0;
            lhs(7,3)=-clhs100*clhs54 - clhs101*clhs55 - clhs115*clhs167 + clhs115*clhs179 + clhs116*clhs178 + clhs152 + clhs176*clhs53 + clhs177*clhs68;
            lhs(7,4)=-clhs103*clhs54 + clhs118*clhs177 - clhs119*clhs167 + clhs119*clhs179 + clhs120*clhs178 + clhs158 + clhs180*clhs53 - clhs69*clhs99;
            lhs(7,5)=0;
            lhs(7,6)=-clhs101*clhs83 - clhs121*clhs167 + clhs121*clhs179 + clhs122*clhs178 + clhs175 + clhs176*clhs81 + clhs177*clhs96 - clhs82*clhs99;
            lhs(7,7)=-clhs101*clhs82 + clhs124*clhs177 - clhs125*clhs167 + clhs125*clhs179 + clhs126*clhs178 + clhs172*clhs29 - clhs173*clhs30 + clhs174*clhs29 + clhs180*clhs81 - clhs97*clhs99;
            lhs(7,8)=0;
            lhs(8,0)=clhs133;
            lhs(8,1)=clhs134;
            lhs(8,2)=0;
            lhs(8,3)=clhs161;
            lhs(8,4)=clhs162;
            lhs(8,5)=0;
            lhs(8,6)=-clhs181*normal[0];
            lhs(8,7)=-clhs181*normal[1];
            lhs(8,8)=0;


    return lhs;
}

template <class TBaseElement>
typename EmbeddedFluidElementDiscontinuous<TBaseElement>::MatrixType EmbeddedFluidElementDiscontinuous<TBaseElement>::NitscheTermsLHS3D(
    EmbeddedDiscontinuousElementData& rData)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;

    const double dyn_tau = rData.DynamicTau;

    const auto vconv = rData.Velocity - rData.MeshVelocity;

    // If there is embedded velocity, substract it to the previous iteration solution
    array_1d<double,3> g = ZeroVector(3);
    if (this->Has(EMBEDDED_VELOCITY)){
        g = this->GetValue(EMBEDDED_VELOCITY);
    }

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;
    const auto& normal = rData.Normal;

    // Get Nitsche imposition values
    const double eps = rData.SlipLength;;
    const double gamma = rData.PenaltyCoefficient;
    const double adjoint = 1.0;

    MatrixType lhs = ZeroMatrix(16,16);

    const double clhs0 =             pow(normal[0], 2);
const double clhs1 =             pow(N[0], 2);
const double clhs2 =             1.0/gamma;
const double clhs3 =             1.0/h;
const double clhs4 =             clhs1*clhs2*clhs3*mu;
const double clhs5 =             clhs0 - 1;
const double clhs6 =             1.0/(eps + gamma*h);
const double clhs7 =             clhs1*clhs6*mu;
const double clhs8 =             1.0*DN(0,1)*N[0]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs9 =             1.0*DN(0,2)*N[0]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs10 =             DN(0,0)*normal[0];
const double clhs11 =             DN(0,1)*normal[1];
const double clhs12 =             1.0*clhs11;
const double clhs13 =             DN(0,2)*normal[2];
const double clhs14 =             1.0*clhs13;
const double clhs15 =             2.0*clhs10 + clhs12 + clhs14;
const double clhs16 =             N[0]*adjoint*clhs5*clhs6*gamma*h*mu;
const double clhs17 =             N[0]*adjoint;
const double clhs18 =             normal[0]*normal[1];
const double clhs19 =             C(0,3)*DN(0,0);
const double clhs20 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs19;
const double clhs21 =             C(1,3)*DN(0,1);
const double clhs22 =             C(3,4)*DN(0,1);
const double clhs23 =             C(4,5)*DN(0,2);
const double clhs24 =             C(0,4)*DN(0,0) + clhs22 + clhs23;
const double clhs25 =             clhs20*normal[0] + clhs24*normal[2] + normal[1]*(C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs21);
const double clhs26 =             clhs18*clhs25;
const double clhs27 =             normal[0]*normal[2];
const double clhs28 =             C(0,5)*DN(0,0);
const double clhs29 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs28;
const double clhs30 =             C(2,5)*DN(0,2);
const double clhs31 =             clhs24*normal[1] + clhs29*normal[0] + normal[2]*(C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs30);
const double clhs32 =             clhs27*clhs31;
const double clhs33 =             N[0]*adjoint*clhs0;
const double clhs34 =             clhs20*normal[1] + clhs29*normal[2] + normal[0]*(C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2));
const double clhs35 =             h*rho*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2)) + mu + pow(h, 2)*rho/dt;
const double clhs36 =             clhs1*clhs2*clhs3*clhs35;
const double clhs37 =             N[0]*clhs6*eps;
const double clhs38 =             -clhs0 + 1;
const double clhs39 =             clhs26 + clhs32 - clhs34*clhs38;
const double clhs40 =             1.0*DN(0,1)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs41 =             clhs18*clhs34;
const double clhs42 =             normal[1]*normal[2];
const double clhs43 =             clhs31*clhs42;
const double clhs44 =             pow(normal[1], 2);
const double clhs45 =             -clhs44 + 1;
const double clhs46 =             -clhs25*clhs45 + clhs41 + clhs43;
const double clhs47 =             1.0*DN(0,2)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs48 =             clhs27*clhs34;
const double clhs49 =             clhs25*clhs42;
const double clhs50 =             pow(normal[2], 2);
const double clhs51 =             -clhs50 + 1;
const double clhs52 =             -clhs31*clhs51 + clhs48 + clhs49;
const double clhs53 =             adjoint*clhs15*clhs6*eps*gamma*h;
const double clhs54 =             clhs1*normal[0];
const double clhs55 =             clhs2*clhs3*mu*normal[1];
const double clhs56 =             clhs6*mu*normal[1];
const double clhs57 =             1.0*DN(0,2)*adjoint*gamma*h*normal[0]*normal[2];
const double clhs58 =             N[0]*clhs6*mu*normal[1];
const double clhs59 =             clhs2*clhs3*clhs35*normal[1];
const double clhs60 =             clhs54*clhs55 - clhs54*clhs56 + clhs54*clhs59 + clhs57*clhs58;
const double clhs61 =             clhs44 - 1;
const double clhs62 =             1.0*DN(0,1)*adjoint*clhs6*clhs61*gamma*h*mu*normal[0];
const double clhs63 =             N[0]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[1];
const double clhs64 =             N[0]*adjoint*clhs44;
const double clhs65 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs21;
const double clhs66 =             C(1,4)*DN(0,1);
const double clhs67 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs66;
const double clhs68 =             clhs65*normal[0] + clhs67*normal[2] + normal[1]*(C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2));
const double clhs69 =             clhs18*clhs68;
const double clhs70 =             C(3,5)*DN(0,0);
const double clhs71 =             C(1,5)*DN(0,1) + clhs23 + clhs70;
const double clhs72 =             C(2,4)*DN(0,2);
const double clhs73 =             clhs67*normal[1] + clhs71*normal[0] + normal[2]*(C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs72);
const double clhs74 =             clhs27*clhs73;
const double clhs75 =             clhs65*normal[1] + clhs71*normal[2] + normal[0]*(C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs19);
const double clhs76 =             -clhs38*clhs75 + clhs69 + clhs74;
const double clhs77 =             clhs18*clhs75;
const double clhs78 =             clhs42*clhs73;
const double clhs79 =             -clhs45*clhs68 + clhs77 + clhs78;
const double clhs80 =             clhs27*clhs75;
const double clhs81 =             clhs42*clhs68;
const double clhs82 =             -clhs51*clhs73 + clhs80 + clhs81;
const double clhs83 =             clhs2*clhs3*mu*normal[2];
const double clhs84 =             clhs6*mu*normal[2];
const double clhs85 =             1.0*DN(0,1)*adjoint*gamma*h*normal[0]*normal[1];
const double clhs86 =             N[0]*clhs6*mu*normal[2];
const double clhs87 =             clhs2*clhs3*clhs35*normal[2];
const double clhs88 =             clhs54*clhs83 - clhs54*clhs84 + clhs54*clhs87 + clhs85*clhs86;
const double clhs89 =             clhs50 - 1;
const double clhs90 =             1.0*DN(0,2)*adjoint*clhs6*clhs89*gamma*h*mu*normal[0];
const double clhs91 =             N[0]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[2];
const double clhs92 =             N[0]*adjoint*clhs50;
const double clhs93 =             C(2,3)*DN(0,2) + clhs22 + clhs70;
const double clhs94 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs72;
const double clhs95 =             clhs93*normal[0] + clhs94*normal[2] + normal[1]*(C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs66);
const double clhs96 =             clhs18*clhs95;
const double clhs97 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs30;
const double clhs98 =             clhs94*normal[1] + clhs97*normal[0] + normal[2]*(C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0));
const double clhs99 =             clhs27*clhs98;
const double clhs100 =             clhs93*normal[1] + clhs97*normal[2] + normal[0]*(C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs28);
const double clhs101 =             -clhs100*clhs38 + clhs96 + clhs99;
const double clhs102 =             clhs100*clhs18;
const double clhs103 =             clhs42*clhs98;
const double clhs104 =             clhs102 + clhs103 - clhs45*clhs95;
const double clhs105 =             clhs100*clhs27;
const double clhs106 =             clhs42*clhs95;
const double clhs107 =             clhs105 + clhs106 - clhs51*clhs98;
const double clhs108 =             N[0]*N[1]*clhs2*clhs3*mu;
const double clhs109 =             N[0]*N[1]*clhs6*mu;
const double clhs110 =             N[0]*N[1]*clhs2*clhs3*clhs35;
const double clhs111 =             clhs0*clhs108 + clhs0*clhs110 - clhs109*clhs5;
const double clhs112 =             1.0*DN(0,1)*N[1]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs113 =             1.0*DN(0,2)*N[1]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs114 =             N[1]*adjoint*clhs5*clhs6*gamma*h*mu;
const double clhs115 =             N[1]*adjoint;
const double clhs116 =             N[1]*adjoint*clhs0;
const double clhs117 =             C(0,3)*DN(1,0);
const double clhs118 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs117;
const double clhs119 =             C(1,3)*DN(1,1);
const double clhs120 =             C(3,4)*DN(1,1);
const double clhs121 =             C(4,5)*DN(1,2);
const double clhs122 =             C(0,4)*DN(1,0) + clhs120 + clhs121;
const double clhs123 =             clhs118*normal[0] + clhs122*normal[2] + normal[1]*(C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs119);
const double clhs124 =             clhs123*clhs18;
const double clhs125 =             C(0,5)*DN(1,0);
const double clhs126 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs125;
const double clhs127 =             C(2,5)*DN(1,2);
const double clhs128 =             clhs122*normal[1] + clhs126*normal[0] + normal[2]*(C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs127);
const double clhs129 =             clhs128*clhs27;
const double clhs130 =             clhs118*normal[1] + clhs126*normal[2] + normal[0]*(C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2));
const double clhs131 =             clhs124 + clhs129 - clhs130*clhs38;
const double clhs132 =             clhs130*clhs18;
const double clhs133 =             clhs128*clhs42;
const double clhs134 =             -clhs123*clhs45 + clhs132 + clhs133;
const double clhs135 =             clhs130*clhs27;
const double clhs136 =             clhs123*clhs42;
const double clhs137 =             -clhs128*clhs51 + clhs135 + clhs136;
const double clhs138 =             N[0]*N[1]*clhs2*clhs3*mu*normal[0];
const double clhs139 =             clhs138*normal[1];
const double clhs140 =             N[0]*N[1]*clhs6*mu*normal[0];
const double clhs141 =             -clhs140*normal[1];
const double clhs142 =             N[1]*clhs6*mu*normal[1];
const double clhs143 =             clhs110*clhs18;
const double clhs144 =             clhs139 + clhs141 + clhs142*clhs57 + clhs143;
const double clhs145 =             N[1]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[1];
const double clhs146 =             N[1]*adjoint*clhs44;
const double clhs147 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs119;
const double clhs148 =             C(1,4)*DN(1,1);
const double clhs149 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs148;
const double clhs150 =             clhs147*normal[0] + clhs149*normal[2] + normal[1]*(C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2));
const double clhs151 =             clhs150*clhs18;
const double clhs152 =             C(3,5)*DN(1,0);
const double clhs153 =             C(1,5)*DN(1,1) + clhs121 + clhs152;
const double clhs154 =             C(2,4)*DN(1,2);
const double clhs155 =             clhs149*normal[1] + clhs153*normal[0] + normal[2]*(C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs154);
const double clhs156 =             clhs155*clhs27;
const double clhs157 =             clhs147*normal[1] + clhs153*normal[2] + normal[0]*(C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs117);
const double clhs158 =             clhs151 + clhs156 - clhs157*clhs38;
const double clhs159 =             clhs157*clhs18;
const double clhs160 =             clhs155*clhs42;
const double clhs161 =             -clhs150*clhs45 + clhs159 + clhs160;
const double clhs162 =             clhs157*clhs27;
const double clhs163 =             clhs150*clhs42;
const double clhs164 =             -clhs155*clhs51 + clhs162 + clhs163;
const double clhs165 =             clhs138*normal[2];
const double clhs166 =             -clhs140*normal[2];
const double clhs167 =             N[1]*clhs6*mu*normal[2];
const double clhs168 =             clhs110*clhs27;
const double clhs169 =             clhs165 + clhs166 + clhs167*clhs85 + clhs168;
const double clhs170 =             N[1]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[2];
const double clhs171 =             N[1]*adjoint*clhs50;
const double clhs172 =             C(2,3)*DN(1,2) + clhs120 + clhs152;
const double clhs173 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs154;
const double clhs174 =             clhs172*normal[0] + clhs173*normal[2] + normal[1]*(C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs148);
const double clhs175 =             clhs174*clhs18;
const double clhs176 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs127;
const double clhs177 =             clhs173*normal[1] + clhs176*normal[0] + normal[2]*(C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0));
const double clhs178 =             clhs177*clhs27;
const double clhs179 =             clhs172*normal[1] + clhs176*normal[2] + normal[0]*(C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs125);
const double clhs180 =             clhs175 + clhs178 - clhs179*clhs38;
const double clhs181 =             clhs179*clhs18;
const double clhs182 =             clhs177*clhs42;
const double clhs183 =             -clhs174*clhs45 + clhs181 + clhs182;
const double clhs184 =             clhs179*clhs27;
const double clhs185 =             clhs174*clhs42;
const double clhs186 =             -clhs177*clhs51 + clhs184 + clhs185;
const double clhs187 =             N[0]*N[2]*clhs2*clhs3*mu;
const double clhs188 =             N[0]*N[2]*clhs6*mu;
const double clhs189 =             N[0]*N[2]*clhs2*clhs3*clhs35;
const double clhs190 =             clhs0*clhs187 + clhs0*clhs189 - clhs188*clhs5;
const double clhs191 =             1.0*DN(0,1)*N[2]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs192 =             1.0*DN(0,2)*N[2]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs193 =             N[2]*adjoint*clhs5*clhs6*gamma*h*mu;
const double clhs194 =             N[2]*adjoint*normal[0]*normal[1];
const double clhs195 =             N[2]*adjoint*normal[0]*normal[2];
const double clhs196 =             N[2]*adjoint*clhs0;
const double clhs197 =             C(0,3)*DN(2,0);
const double clhs198 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs197;
const double clhs199 =             C(1,3)*DN(2,1);
const double clhs200 =             C(3,4)*DN(2,1);
const double clhs201 =             C(4,5)*DN(2,2);
const double clhs202 =             C(0,4)*DN(2,0) + clhs200 + clhs201;
const double clhs203 =             clhs198*normal[0] + clhs202*normal[2] + normal[1]*(C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs199);
const double clhs204 =             clhs18*clhs203;
const double clhs205 =             C(0,5)*DN(2,0);
const double clhs206 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs205;
const double clhs207 =             C(2,5)*DN(2,2);
const double clhs208 =             clhs202*normal[1] + clhs206*normal[0] + normal[2]*(C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs207);
const double clhs209 =             clhs208*clhs27;
const double clhs210 =             clhs198*normal[1] + clhs206*normal[2] + normal[0]*(C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2));
const double clhs211 =             clhs204 + clhs209 - clhs210*clhs38;
const double clhs212 =             clhs18*clhs210;
const double clhs213 =             clhs208*clhs42;
const double clhs214 =             -clhs203*clhs45 + clhs212 + clhs213;
const double clhs215 =             clhs210*clhs27;
const double clhs216 =             clhs203*clhs42;
const double clhs217 =             -clhs208*clhs51 + clhs215 + clhs216;
const double clhs218 =             N[2]*normal[0];
const double clhs219 =             N[0]*clhs2*clhs3*mu*normal[1];
const double clhs220 =             clhs218*clhs219;
const double clhs221 =             -clhs218*clhs58;
const double clhs222 =             1.0*DN(0,2)*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs223 =             N[0]*N[2]*clhs2*clhs3*clhs35*normal[0];
const double clhs224 =             clhs223*normal[1];
const double clhs225 =             clhs218*clhs222 + clhs220 + clhs221 + clhs224;
const double clhs226 =             1.0*DN(0,1)*adjoint*clhs6*clhs61*gamma*h*mu;
const double clhs227 =             N[2]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[1];
const double clhs228 =             N[2]*adjoint*normal[1]*normal[2];
const double clhs229 =             N[2]*adjoint*clhs44;
const double clhs230 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs199;
const double clhs231 =             C(1,4)*DN(2,1);
const double clhs232 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs231;
const double clhs233 =             clhs230*normal[0] + clhs232*normal[2] + normal[1]*(C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2));
const double clhs234 =             clhs18*clhs233;
const double clhs235 =             C(3,5)*DN(2,0);
const double clhs236 =             C(1,5)*DN(2,1) + clhs201 + clhs235;
const double clhs237 =             C(2,4)*DN(2,2);
const double clhs238 =             clhs232*normal[1] + clhs236*normal[0] + normal[2]*(C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs237);
const double clhs239 =             clhs238*clhs27;
const double clhs240 =             clhs230*normal[1] + clhs236*normal[2] + normal[0]*(C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs197);
const double clhs241 =             clhs234 + clhs239 - clhs240*clhs38;
const double clhs242 =             clhs18*clhs240;
const double clhs243 =             clhs238*clhs42;
const double clhs244 =             -clhs233*clhs45 + clhs242 + clhs243;
const double clhs245 =             clhs240*clhs27;
const double clhs246 =             clhs233*clhs42;
const double clhs247 =             -clhs238*clhs51 + clhs245 + clhs246;
const double clhs248 =             N[0]*clhs2*clhs3*mu*normal[2];
const double clhs249 =             clhs218*clhs248;
const double clhs250 =             -clhs218*clhs86;
const double clhs251 =             1.0*DN(0,1)*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs252 =             clhs223*normal[2];
const double clhs253 =             clhs218*clhs251 + clhs249 + clhs250 + clhs252;
const double clhs254 =             1.0*DN(0,2)*adjoint*clhs6*clhs89*gamma*h*mu;
const double clhs255 =             N[2]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[2];
const double clhs256 =             N[2]*adjoint*clhs50;
const double clhs257 =             C(2,3)*DN(2,2) + clhs200 + clhs235;
const double clhs258 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs237;
const double clhs259 =             clhs257*normal[0] + clhs258*normal[2] + normal[1]*(C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs231);
const double clhs260 =             clhs18*clhs259;
const double clhs261 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs207;
const double clhs262 =             clhs258*normal[1] + clhs261*normal[0] + normal[2]*(C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0));
const double clhs263 =             clhs262*clhs27;
const double clhs264 =             clhs257*normal[1] + clhs261*normal[2] + normal[0]*(C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs205);
const double clhs265 =             clhs260 + clhs263 - clhs264*clhs38;
const double clhs266 =             clhs18*clhs264;
const double clhs267 =             clhs262*clhs42;
const double clhs268 =             -clhs259*clhs45 + clhs266 + clhs267;
const double clhs269 =             clhs264*clhs27;
const double clhs270 =             clhs259*clhs42;
const double clhs271 =             -clhs262*clhs51 + clhs269 + clhs270;
const double clhs272 =             N[0]*N[3]*clhs2*clhs3*mu;
const double clhs273 =             N[0]*N[3]*clhs6*mu;
const double clhs274 =             N[0]*N[3]*clhs2*clhs3*clhs35;
const double clhs275 =             clhs0*clhs272 + clhs0*clhs274 - clhs273*clhs5;
const double clhs276 =             1.0*DN(0,1)*N[3]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs277 =             1.0*DN(0,2)*N[3]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs278 =             N[3]*adjoint*clhs5*clhs6*gamma*h*mu;
const double clhs279 =             N[3]*adjoint*normal[0]*normal[1];
const double clhs280 =             N[3]*adjoint*normal[0]*normal[2];
const double clhs281 =             N[3]*adjoint*clhs0;
const double clhs282 =             C(0,3)*DN(3,0);
const double clhs283 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs282;
const double clhs284 =             C(1,3)*DN(3,1);
const double clhs285 =             C(3,4)*DN(3,1);
const double clhs286 =             C(4,5)*DN(3,2);
const double clhs287 =             C(0,4)*DN(3,0) + clhs285 + clhs286;
const double clhs288 =             clhs283*normal[0] + clhs287*normal[2] + normal[1]*(C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs284);
const double clhs289 =             clhs18*clhs288;
const double clhs290 =             C(0,5)*DN(3,0);
const double clhs291 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs290;
const double clhs292 =             C(2,5)*DN(3,2);
const double clhs293 =             clhs287*normal[1] + clhs291*normal[0] + normal[2]*(C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs292);
const double clhs294 =             clhs27*clhs293;
const double clhs295 =             clhs283*normal[1] + clhs291*normal[2] + normal[0]*(C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2));
const double clhs296 =             clhs289 + clhs294 - clhs295*clhs38;
const double clhs297 =             clhs18*clhs295;
const double clhs298 =             clhs293*clhs42;
const double clhs299 =             -clhs288*clhs45 + clhs297 + clhs298;
const double clhs300 =             clhs27*clhs295;
const double clhs301 =             clhs288*clhs42;
const double clhs302 =             -clhs293*clhs51 + clhs300 + clhs301;
const double clhs303 =             N[3]*normal[0];
const double clhs304 =             clhs219*clhs303;
const double clhs305 =             -clhs303*clhs58;
const double clhs306 =             N[0]*N[3]*clhs2*clhs3*clhs35*normal[0]*normal[1];
const double clhs307 =             clhs222*clhs303 + clhs304 + clhs305 + clhs306;
const double clhs308 =             N[3]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[1];
const double clhs309 =             N[3]*adjoint*normal[1]*normal[2];
const double clhs310 =             N[3]*adjoint*clhs44;
const double clhs311 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs284;
const double clhs312 =             C(1,4)*DN(3,1);
const double clhs313 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs312;
const double clhs314 =             clhs311*normal[0] + clhs313*normal[2] + normal[1]*(C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2));
const double clhs315 =             clhs18*clhs314;
const double clhs316 =             C(3,5)*DN(3,0);
const double clhs317 =             C(1,5)*DN(3,1) + clhs286 + clhs316;
const double clhs318 =             C(2,4)*DN(3,2);
const double clhs319 =             clhs313*normal[1] + clhs317*normal[0] + normal[2]*(C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs318);
const double clhs320 =             clhs27*clhs319;
const double clhs321 =             clhs311*normal[1] + clhs317*normal[2] + normal[0]*(C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs282);
const double clhs322 =             clhs315 + clhs320 - clhs321*clhs38;
const double clhs323 =             clhs18*clhs321;
const double clhs324 =             clhs319*clhs42;
const double clhs325 =             -clhs314*clhs45 + clhs323 + clhs324;
const double clhs326 =             clhs27*clhs321;
const double clhs327 =             clhs314*clhs42;
const double clhs328 =             -clhs319*clhs51 + clhs326 + clhs327;
const double clhs329 =             clhs248*clhs303;
const double clhs330 =             -clhs303*clhs86;
const double clhs331 =             N[0]*clhs2*clhs3*clhs35*normal[2];
const double clhs332 =             clhs303*clhs331;
const double clhs333 =             clhs251*clhs303 + clhs329 + clhs330 + clhs332;
const double clhs334 =             N[3]*adjoint*clhs6*gamma*h*mu*normal[0]*normal[2];
const double clhs335 =             N[3]*adjoint*clhs50;
const double clhs336 =             C(2,3)*DN(3,2) + clhs285 + clhs316;
const double clhs337 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs318;
const double clhs338 =             clhs336*normal[0] + clhs337*normal[2] + normal[1]*(C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs312);
const double clhs339 =             clhs18*clhs338;
const double clhs340 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs292;
const double clhs341 =             clhs337*normal[1] + clhs340*normal[0] + normal[2]*(C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0));
const double clhs342 =             clhs27*clhs341;
const double clhs343 =             clhs336*normal[1] + clhs340*normal[2] + normal[0]*(C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs290);
const double clhs344 =             clhs339 + clhs342 - clhs343*clhs38;
const double clhs345 =             clhs18*clhs343;
const double clhs346 =             clhs341*clhs42;
const double clhs347 =             -clhs338*clhs45 + clhs345 + clhs346;
const double clhs348 =             clhs27*clhs343;
const double clhs349 =             clhs338*clhs42;
const double clhs350 =             -clhs341*clhs51 + clhs348 + clhs349;
const double clhs351 =             1.0*DN(0,0)*adjoint*clhs5*gamma*h;
const double clhs352 =             1.0*clhs10;
const double clhs353 =             2.0*clhs11 + clhs14 + clhs352;
const double clhs354 =             1.0*DN(0,0)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs355 =             1.0*DN(0,2)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs356 =             adjoint*clhs353*clhs6*eps*gamma*h;
const double clhs357 =             1.0*DN(0,0)*N[0]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs358 =             N[0]*adjoint*clhs6*clhs61*gamma*h*mu;
const double clhs359 =             clhs2*clhs3*mu*normal[1]*normal[2];
const double clhs360 =             clhs6*mu*normal[1]*normal[2];
const double clhs361 =             1.0*DN(0,0)*adjoint*gamma*h*normal[0];
const double clhs362 =             N[0]*clhs6*mu*normal[1]*normal[2];
const double clhs363 =             clhs2*clhs3*clhs35*normal[1]*normal[2];
const double clhs364 =             clhs1*clhs359 - clhs1*clhs360 + clhs1*clhs363 + clhs361*clhs362;
const double clhs365 =             1.0*DN(0,2)*adjoint*clhs89*gamma*h;
const double clhs366 =             N[0]*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs367 =             clhs108*clhs44 - clhs109*clhs61 + clhs110*clhs44;
const double clhs368 =             1.0*DN(0,0)*N[1]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs369 =             N[1]*adjoint*clhs6*clhs61*gamma*h*mu;
const double clhs370 =             N[0]*clhs2*clhs3*mu*normal[1]*normal[2];
const double clhs371 =             N[1]*clhs370;
const double clhs372 =             -N[1]*clhs362;
const double clhs373 =             N[1]*clhs6*mu*normal[1]*normal[2];
const double clhs374 =             clhs110*clhs42;
const double clhs375 =             clhs361*clhs373 + clhs371 + clhs372 + clhs374;
const double clhs376 =             N[1]*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs377 =             N[2]*clhs6*mu*normal[1];
const double clhs378 =             clhs187*clhs44 - clhs188*clhs61 + clhs189*clhs44;
const double clhs379 =             1.0*DN(0,0)*N[2]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs380 =             N[2]*adjoint*clhs6*clhs61*gamma*h*mu;
const double clhs381 =             N[2]*clhs370;
const double clhs382 =             N[2]*normal[1];
const double clhs383 =             -clhs382*clhs86;
const double clhs384 =             1.0*DN(0,0)*adjoint*clhs6*gamma*h*mu*normal[0]*normal[2];
const double clhs385 =             clhs331*clhs382;
const double clhs386 =             clhs381 + clhs382*clhs384 + clhs383 + clhs385;
const double clhs387 =             N[2]*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs388 =             N[3]*clhs6*mu*normal[1];
const double clhs389 =             clhs272*clhs44 - clhs273*clhs61 + clhs274*clhs44;
const double clhs390 =             1.0*DN(0,0)*N[3]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs391 =             N[3]*adjoint*clhs6*clhs61*gamma*h*mu;
const double clhs392 =             N[3]*clhs370;
const double clhs393 =             -N[3]*clhs362;
const double clhs394 =             N[3]*normal[1];
const double clhs395 =             clhs331*clhs394;
const double clhs396 =             clhs384*clhs394 + clhs392 + clhs393 + clhs395;
const double clhs397 =             N[3]*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs398 =             clhs12 + 2.0*clhs13 + clhs352;
const double clhs399 =             1.0*DN(0,0)*adjoint*clhs6*eps*gamma*h*normal[2];
const double clhs400 =             1.0*DN(0,1)*adjoint*clhs6*eps*gamma*h*normal[2];
const double clhs401 =             adjoint*clhs398*clhs6*eps*gamma*h;
const double clhs402 =             1.0*DN(0,1)*adjoint*clhs61*gamma*h;
const double clhs403 =             N[0]*adjoint*clhs6*clhs89*gamma*h*mu;
const double clhs404 =             clhs108*clhs50 - clhs109*clhs89 + clhs110*clhs50;
const double clhs405 =             N[1]*adjoint*clhs6*clhs89*gamma*h*mu;
const double clhs406 =             1.0*DN(0,0)*adjoint*clhs5*clhs6*gamma*h*mu*normal[2];
const double clhs407 =             1.0*DN(0,1)*adjoint*clhs6*clhs61*gamma*h*mu*normal[2];
const double clhs408 =             clhs187*clhs50 - clhs188*clhs89 + clhs189*clhs50;
const double clhs409 =             N[2]*adjoint*clhs6*clhs89*gamma*h*mu;
const double clhs410 =             clhs272*clhs50 - clhs273*clhs89 + clhs274*clhs50;
const double clhs411 =             N[3]*adjoint*clhs6*clhs89*gamma*h*mu;
const double clhs412 =             clhs0 + clhs44 + clhs50;
const double clhs413 =             clhs412*normal[1];
const double clhs414 =             clhs412*normal[2];
const double clhs415 =             N[0]*N[1]*clhs412;
const double clhs416 =             -clhs415*normal[0];
const double clhs417 =             -clhs415*normal[1];
const double clhs418 =             -clhs415*normal[2];
const double clhs419 =             N[0]*clhs412;
const double clhs420 =             -clhs218*clhs419;
const double clhs421 =             -clhs382*clhs419;
const double clhs422 =             N[0]*clhs412*normal[2];
const double clhs423 =             -N[2]*clhs422;
const double clhs424 =             -clhs303*clhs419;
const double clhs425 =             -clhs394*clhs419;
const double clhs426 =             -N[3]*clhs422;
const double clhs427 =             1.0*DN(1,1)*N[0]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs428 =             1.0*DN(1,2)*N[0]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs429 =             DN(1,0)*normal[0];
const double clhs430 =             DN(1,1)*normal[1];
const double clhs431 =             1.0*clhs430;
const double clhs432 =             DN(1,2)*normal[2];
const double clhs433 =             1.0*clhs432;
const double clhs434 =             2.0*clhs429 + clhs431 + clhs433;
const double clhs435 =             N[1]*clhs6*eps;
const double clhs436 =             1.0*DN(1,1)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs437 =             1.0*DN(1,2)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs438 =             adjoint*clhs434*clhs6*eps*gamma*h;
const double clhs439 =             1.0*DN(1,2)*adjoint*gamma*h*normal[0]*normal[2];
const double clhs440 =             clhs139 + clhs141 + clhs143 + clhs439*clhs58;
const double clhs441 =             1.0*DN(1,1)*adjoint*clhs6*clhs61*gamma*h*mu*normal[0];
const double clhs442 =             1.0*DN(1,1)*adjoint*gamma*h*normal[0]*normal[1];
const double clhs443 =             clhs165 + clhs166 + clhs168 + clhs442*clhs86;
const double clhs444 =             1.0*DN(1,2)*adjoint*clhs6*clhs89*gamma*h*mu*normal[0];
const double clhs445 =             pow(N[1], 2);
const double clhs446 =             clhs2*clhs3*clhs445*mu;
const double clhs447 =             clhs445*clhs6*mu;
const double clhs448 =             1.0*DN(1,1)*N[1]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs449 =             1.0*DN(1,2)*N[1]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs450 =             clhs2*clhs3*clhs35*clhs445;
const double clhs451 =             clhs445*normal[0];
const double clhs452 =             clhs142*clhs439 + clhs451*clhs55 - clhs451*clhs56 + clhs451*clhs59;
const double clhs453 =             clhs167*clhs442 + clhs451*clhs83 - clhs451*clhs84 + clhs451*clhs87;
const double clhs454 =             N[1]*N[2]*clhs2*clhs3*mu;
const double clhs455 =             N[1]*N[2]*clhs6*mu;
const double clhs456 =             N[1]*N[2]*clhs2*clhs3*clhs35;
const double clhs457 =             clhs0*clhs454 + clhs0*clhs456 - clhs455*clhs5;
const double clhs458 =             1.0*DN(1,1)*N[2]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs459 =             1.0*DN(1,2)*N[2]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs460 =             N[1]*N[2]*clhs2*clhs3*mu*normal[0];
const double clhs461 =             clhs460*normal[1];
const double clhs462 =             N[1]*N[2]*clhs6*mu*normal[0];
const double clhs463 =             -clhs462*normal[1];
const double clhs464 =             1.0*DN(1,2)*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs465 =             N[1]*N[2]*clhs2*clhs3*clhs35*normal[0];
const double clhs466 =             clhs465*normal[1];
const double clhs467 =             clhs218*clhs464 + clhs461 + clhs463 + clhs466;
const double clhs468 =             1.0*DN(1,1)*adjoint*clhs6*clhs61*gamma*h*mu;
const double clhs469 =             clhs460*normal[2];
const double clhs470 =             -clhs462*normal[2];
const double clhs471 =             1.0*DN(1,1)*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs472 =             clhs465*normal[2];
const double clhs473 =             clhs218*clhs471 + clhs469 + clhs470 + clhs472;
const double clhs474 =             1.0*DN(1,2)*adjoint*clhs6*clhs89*gamma*h*mu;
const double clhs475 =             N[1]*N[3]*clhs2*clhs3*mu;
const double clhs476 =             N[1]*N[3]*clhs6*mu;
const double clhs477 =             N[1]*N[3]*clhs2*clhs3*clhs35;
const double clhs478 =             clhs0*clhs475 + clhs0*clhs477 - clhs476*clhs5;
const double clhs479 =             1.0*DN(1,1)*N[3]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs480 =             1.0*DN(1,2)*N[3]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs481 =             N[1]*clhs2*clhs3*clhs303*mu*normal[1];
const double clhs482 =             N[1]*N[3]*normal[0];
const double clhs483 =             -clhs482*clhs56;
const double clhs484 =             clhs482*clhs59;
const double clhs485 =             clhs303*clhs464 + clhs481 + clhs483 + clhs484;
const double clhs486 =             clhs482*clhs83;
const double clhs487 =             -clhs482*clhs84;
const double clhs488 =             clhs482*clhs87;
const double clhs489 =             clhs303*clhs471 + clhs486 + clhs487 + clhs488;
const double clhs490 =             1.0*DN(1,0)*adjoint*clhs5*gamma*h;
const double clhs491 =             1.0*clhs429;
const double clhs492 =             2.0*clhs430 + clhs433 + clhs491;
const double clhs493 =             1.0*DN(1,0)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs494 =             1.0*DN(1,2)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs495 =             adjoint*clhs492*clhs6*eps*gamma*h;
const double clhs496 =             1.0*DN(1,0)*N[0]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs497 =             1.0*DN(1,0)*adjoint*gamma*h*normal[0];
const double clhs498 =             clhs362*clhs497 + clhs371 + clhs372 + clhs374;
const double clhs499 =             1.0*DN(1,2)*adjoint*clhs89*gamma*h;
const double clhs500 =             1.0*DN(1,0)*N[1]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs501 =             clhs359*clhs445 - clhs360*clhs445 + clhs363*clhs445 + clhs373*clhs497;
const double clhs502 =             clhs44*clhs454 + clhs44*clhs456 - clhs455*clhs61;
const double clhs503 =             1.0*DN(1,0)*N[2]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs504 =             N[1]*clhs2*clhs3*mu*normal[1]*normal[2];
const double clhs505 =             N[2]*clhs504;
const double clhs506 =             -clhs167*clhs382;
const double clhs507 =             1.0*DN(1,0)*adjoint*clhs6*gamma*h*mu*normal[0]*normal[2];
const double clhs508 =             N[1]*clhs2*clhs3*clhs35*normal[2];
const double clhs509 =             clhs382*clhs508;
const double clhs510 =             clhs382*clhs507 + clhs505 + clhs506 + clhs509;
const double clhs511 =             clhs44*clhs475 + clhs44*clhs477 - clhs476*clhs61;
const double clhs512 =             1.0*DN(1,0)*N[3]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs513 =             N[3]*clhs504;
const double clhs514 =             -clhs167*clhs394;
const double clhs515 =             clhs394*clhs508;
const double clhs516 =             clhs394*clhs507 + clhs513 + clhs514 + clhs515;
const double clhs517 =             clhs431 + 2.0*clhs432 + clhs491;
const double clhs518 =             1.0*DN(1,0)*adjoint*clhs6*eps*gamma*h*normal[2];
const double clhs519 =             1.0*DN(1,1)*adjoint*clhs6*eps*gamma*h*normal[2];
const double clhs520 =             adjoint*clhs517*clhs6*eps*gamma*h;
const double clhs521 =             1.0*DN(1,1)*adjoint*clhs61*gamma*h;
const double clhs522 =             1.0*DN(1,0)*adjoint*clhs5*clhs6*gamma*h*mu*normal[2];
const double clhs523 =             1.0*DN(1,1)*adjoint*clhs6*clhs61*gamma*h*mu*normal[2];
const double clhs524 =             clhs454*clhs50 - clhs455*clhs89 + clhs456*clhs50;
const double clhs525 =             clhs475*clhs50 - clhs476*clhs89 + clhs477*clhs50;
const double clhs526 =             N[1]*clhs412;
const double clhs527 =             -clhs218*clhs526;
const double clhs528 =             -clhs382*clhs526;
const double clhs529 =             N[1]*clhs412*normal[2];
const double clhs530 =             -N[2]*clhs529;
const double clhs531 =             -clhs412*clhs482;
const double clhs532 =             -clhs394*clhs526;
const double clhs533 =             -N[3]*clhs529;
const double clhs534 =             1.0*DN(2,1)*N[0]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs535 =             1.0*DN(2,2)*N[0]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs536 =             DN(2,0)*normal[0];
const double clhs537 =             DN(2,1)*normal[1];
const double clhs538 =             1.0*clhs537;
const double clhs539 =             DN(2,2)*normal[2];
const double clhs540 =             1.0*clhs539;
const double clhs541 =             2.0*clhs536 + clhs538 + clhs540;
const double clhs542 =             N[2]*clhs6*eps;
const double clhs543 =             1.0*DN(2,1)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs544 =             1.0*DN(2,2)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs545 =             adjoint*clhs541*clhs6*eps*gamma*h;
const double clhs546 =             1.0*DN(2,2)*adjoint*gamma*h*normal[0]*normal[2];
const double clhs547 =             clhs220 + clhs221 + clhs224 + clhs546*clhs58;
const double clhs548 =             1.0*DN(2,1)*adjoint*clhs6*clhs61*gamma*h*mu*normal[0];
const double clhs549 =             1.0*DN(2,1)*adjoint*gamma*h*normal[0]*normal[1];
const double clhs550 =             clhs249 + clhs250 + clhs252 + clhs549*clhs86;
const double clhs551 =             1.0*DN(2,2)*adjoint*clhs6*clhs89*gamma*h*mu*normal[0];
const double clhs552 =             1.0*DN(2,1)*N[1]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs553 =             1.0*DN(2,2)*N[1]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs554 =             clhs142*clhs546 + clhs461 + clhs463 + clhs466;
const double clhs555 =             clhs167*clhs549 + clhs469 + clhs470 + clhs472;
const double clhs556 =             pow(N[2], 2);
const double clhs557 =             clhs2*clhs3*clhs556*mu;
const double clhs558 =             clhs556*clhs6*mu;
const double clhs559 =             1.0*DN(2,1)*N[2]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs560 =             1.0*DN(2,2)*N[2]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs561 =             clhs2*clhs3*clhs35*clhs556;
const double clhs562 =             clhs556*normal[0];
const double clhs563 =             1.0*DN(2,2)*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs564 =             clhs218*clhs563 + clhs55*clhs562 - clhs56*clhs562 + clhs562*clhs59;
const double clhs565 =             1.0*DN(2,1)*adjoint*clhs6*clhs61*gamma*h*mu;
const double clhs566 =             1.0*DN(2,1)*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs567 =             clhs218*clhs566 + clhs562*clhs83 - clhs562*clhs84 + clhs562*clhs87;
const double clhs568 =             1.0*DN(2,2)*adjoint*clhs6*clhs89*gamma*h*mu;
const double clhs569 =             N[2]*N[3]*clhs2*clhs3*mu;
const double clhs570 =             N[2]*N[3]*clhs6*mu;
const double clhs571 =             N[2]*N[3]*clhs2*clhs3*clhs35;
const double clhs572 =             clhs0*clhs569 + clhs0*clhs571 - clhs5*clhs570;
const double clhs573 =             1.0*DN(2,1)*N[3]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs574 =             1.0*DN(2,2)*N[3]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs575 =             N[2]*N[3]*normal[0];
const double clhs576 =             clhs55*clhs575;
const double clhs577 =             -clhs56*clhs575;
const double clhs578 =             clhs575*clhs59;
const double clhs579 =             clhs303*clhs563 + clhs576 + clhs577 + clhs578;
const double clhs580 =             clhs575*clhs83;
const double clhs581 =             -clhs575*clhs84;
const double clhs582 =             clhs575*clhs87;
const double clhs583 =             clhs303*clhs566 + clhs580 + clhs581 + clhs582;
const double clhs584 =             1.0*DN(2,0)*adjoint*clhs5*gamma*h;
const double clhs585 =             1.0*clhs536;
const double clhs586 =             2.0*clhs537 + clhs540 + clhs585;
const double clhs587 =             1.0*DN(2,0)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs588 =             1.0*DN(2,2)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs589 =             adjoint*clhs586*clhs6*eps*gamma*h;
const double clhs590 =             1.0*DN(2,0)*N[0]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs591 =             1.0*DN(2,0)*adjoint*gamma*h*normal[0];
const double clhs592 =             clhs362*clhs591 + clhs381 + clhs383 + clhs385;
const double clhs593 =             1.0*DN(2,2)*adjoint*clhs89*gamma*h;
const double clhs594 =             1.0*DN(2,0)*N[1]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs595 =             clhs373*clhs591 + clhs505 + clhs506 + clhs509;
const double clhs596 =             1.0*DN(2,0)*N[2]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs597 =             1.0*DN(2,0)*adjoint*clhs6*gamma*h*mu*normal[0]*normal[2];
const double clhs598 =             clhs359*clhs556 - clhs360*clhs556 + clhs363*clhs556 + clhs382*clhs597;
const double clhs599 =             clhs44*clhs569 + clhs44*clhs571 - clhs570*clhs61;
const double clhs600 =             1.0*DN(2,0)*N[3]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs601 =             N[2]*N[3]*normal[1];
const double clhs602 =             clhs601*clhs83;
const double clhs603 =             -clhs601*clhs84;
const double clhs604 =             clhs601*clhs87;
const double clhs605 =             clhs394*clhs597 + clhs602 + clhs603 + clhs604;
const double clhs606 =             clhs538 + 2.0*clhs539 + clhs585;
const double clhs607 =             1.0*DN(2,0)*adjoint*clhs6*eps*gamma*h*normal[2];
const double clhs608 =             1.0*DN(2,1)*adjoint*clhs6*eps*gamma*h*normal[2];
const double clhs609 =             adjoint*clhs6*clhs606*eps*gamma*h;
const double clhs610 =             1.0*DN(2,1)*adjoint*clhs61*gamma*h;
const double clhs611 =             1.0*DN(2,0)*adjoint*clhs5*clhs6*gamma*h*mu*normal[2];
const double clhs612 =             1.0*DN(2,1)*adjoint*clhs6*clhs61*gamma*h*mu*normal[2];
const double clhs613 =             clhs50*clhs569 + clhs50*clhs571 - clhs570*clhs89;
const double clhs614 =             -clhs412*clhs575;
const double clhs615 =             -clhs412*clhs601;
const double clhs616 =             -N[2]*N[3]*clhs414;
const double clhs617 =             1.0*DN(3,1)*N[0]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs618 =             1.0*DN(3,2)*N[0]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs619 =             DN(3,0)*normal[0];
const double clhs620 =             DN(3,1)*normal[1];
const double clhs621 =             1.0*clhs620;
const double clhs622 =             DN(3,2)*normal[2];
const double clhs623 =             1.0*clhs622;
const double clhs624 =             2.0*clhs619 + clhs621 + clhs623;
const double clhs625 =             N[3]*clhs6*eps;
const double clhs626 =             1.0*DN(3,1)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs627 =             1.0*DN(3,2)*adjoint*clhs6*eps*gamma*h*normal[0];
const double clhs628 =             adjoint*clhs6*clhs624*eps*gamma*h;
const double clhs629 =             1.0*DN(3,2)*adjoint*gamma*h*normal[0]*normal[2];
const double clhs630 =             clhs304 + clhs305 + clhs306 + clhs58*clhs629;
const double clhs631 =             1.0*DN(3,1)*adjoint*clhs6*clhs61*gamma*h*mu*normal[0];
const double clhs632 =             1.0*DN(3,1)*adjoint*gamma*h*normal[0]*normal[1];
const double clhs633 =             clhs329 + clhs330 + clhs332 + clhs632*clhs86;
const double clhs634 =             1.0*DN(3,2)*adjoint*clhs6*clhs89*gamma*h*mu*normal[0];
const double clhs635 =             1.0*DN(3,1)*N[1]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs636 =             1.0*DN(3,2)*N[1]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs637 =             clhs142*clhs629 + clhs481 + clhs483 + clhs484;
const double clhs638 =             clhs167*clhs632 + clhs486 + clhs487 + clhs488;
const double clhs639 =             1.0*DN(3,1)*N[2]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs640 =             1.0*DN(3,2)*N[2]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs641 =             1.0*DN(3,2)*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs642 =             clhs218*clhs641 + clhs576 + clhs577 + clhs578;
const double clhs643 =             1.0*DN(3,1)*adjoint*clhs6*clhs61*gamma*h*mu;
const double clhs644 =             1.0*DN(3,1)*adjoint*clhs6*gamma*h*mu*normal[1]*normal[2];
const double clhs645 =             clhs218*clhs644 + clhs580 + clhs581 + clhs582;
const double clhs646 =             1.0*DN(3,2)*adjoint*clhs6*clhs89*gamma*h*mu;
const double clhs647 =             pow(N[3], 2);
const double clhs648 =             clhs2*clhs3*clhs647*mu;
const double clhs649 =             clhs6*clhs647*mu;
const double clhs650 =             1.0*DN(3,1)*N[3]*adjoint*clhs6*gamma*h*mu*normal[1];
const double clhs651 =             1.0*DN(3,2)*N[3]*adjoint*clhs6*gamma*h*mu*normal[2];
const double clhs652 =             clhs2*clhs3*clhs35*clhs647;
const double clhs653 =             clhs647*normal[0];
const double clhs654 =             clhs303*clhs641 + clhs55*clhs653 - clhs56*clhs653 + clhs59*clhs653;
const double clhs655 =             clhs303*clhs644 + clhs653*clhs83 - clhs653*clhs84 + clhs653*clhs87;
const double clhs656 =             1.0*DN(3,0)*adjoint*clhs5*gamma*h;
const double clhs657 =             1.0*clhs619;
const double clhs658 =             2.0*clhs620 + clhs623 + clhs657;
const double clhs659 =             1.0*DN(3,0)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs660 =             1.0*DN(3,2)*adjoint*clhs6*eps*gamma*h*normal[1];
const double clhs661 =             adjoint*clhs6*clhs658*eps*gamma*h;
const double clhs662 =             1.0*DN(3,0)*N[0]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs663 =             1.0*DN(3,0)*adjoint*gamma*h*normal[0];
const double clhs664 =             clhs362*clhs663 + clhs392 + clhs393 + clhs395;
const double clhs665 =             1.0*DN(3,2)*adjoint*clhs89*gamma*h;
const double clhs666 =             1.0*DN(3,0)*N[1]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs667 =             clhs373*clhs663 + clhs513 + clhs514 + clhs515;
const double clhs668 =             1.0*DN(3,0)*N[2]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs669 =             1.0*DN(3,0)*adjoint*clhs6*gamma*h*mu*normal[0]*normal[2];
const double clhs670 =             clhs382*clhs669 + clhs602 + clhs603 + clhs604;
const double clhs671 =             1.0*DN(3,0)*N[3]*adjoint*clhs6*gamma*h*mu*normal[0];
const double clhs672 =             clhs359*clhs647 - clhs360*clhs647 + clhs363*clhs647 + clhs394*clhs669;
const double clhs673 =             clhs621 + 2.0*clhs622 + clhs657;
const double clhs674 =             1.0*DN(3,0)*adjoint*clhs6*eps*gamma*h*normal[2];
const double clhs675 =             1.0*DN(3,1)*adjoint*clhs6*eps*gamma*h*normal[2];
const double clhs676 =             adjoint*clhs6*clhs673*eps*gamma*h;
const double clhs677 =             1.0*DN(3,1)*adjoint*clhs61*gamma*h;
const double clhs678 =             1.0*DN(3,0)*adjoint*clhs5*clhs6*gamma*h*mu*normal[2];
const double clhs679 =             1.0*DN(3,1)*adjoint*clhs6*clhs61*gamma*h*mu*normal[2];
            lhs(0,0)=clhs0*clhs36 + clhs0*clhs4 + clhs0*clhs8 + clhs0*clhs9 + clhs15*clhs16 - clhs17*clhs26 - clhs17*clhs32 - clhs33*clhs34 - clhs37*clhs39 + clhs39*clhs53 + clhs40*clhs46 + clhs47*clhs52 - clhs5*clhs7;
            lhs(0,1)=N[0]*clhs62 + clhs15*clhs63 - clhs17*clhs41 - clhs17*clhs43 - clhs25*clhs64 - clhs37*clhs76 + clhs40*clhs79 + clhs47*clhs82 + clhs53*clhs76 + clhs60;
            lhs(0,2)=N[0]*clhs90 - clhs101*clhs37 + clhs101*clhs53 + clhs104*clhs40 + clhs107*clhs47 + clhs15*clhs91 - clhs17*clhs48 - clhs17*clhs49 - clhs31*clhs92 + clhs88;
            lhs(0,3)=0;
            lhs(0,4)=clhs0*clhs112 + clhs0*clhs113 + clhs111 + clhs114*clhs15 - clhs115*clhs26 - clhs115*clhs32 - clhs116*clhs34 - clhs131*clhs37 + clhs131*clhs53 + clhs134*clhs40 + clhs137*clhs47;
            lhs(0,5)=N[1]*clhs62 - clhs115*clhs41 - clhs115*clhs43 + clhs144 + clhs145*clhs15 - clhs146*clhs25 - clhs158*clhs37 + clhs158*clhs53 + clhs161*clhs40 + clhs164*clhs47;
            lhs(0,6)=N[1]*clhs90 - clhs115*clhs48 - clhs115*clhs49 + clhs15*clhs170 + clhs169 - clhs171*clhs31 - clhs180*clhs37 + clhs180*clhs53 + clhs183*clhs40 + clhs186*clhs47;
            lhs(0,7)=0;
            lhs(0,8)=clhs0*clhs191 + clhs0*clhs192 + clhs15*clhs193 + clhs190 - clhs194*clhs25 - clhs195*clhs31 - clhs196*clhs34 - clhs211*clhs37 + clhs211*clhs53 + clhs214*clhs40 + clhs217*clhs47;
            lhs(0,9)=clhs15*clhs227 - clhs194*clhs34 + clhs218*clhs226 + clhs225 - clhs228*clhs31 - clhs229*clhs25 - clhs241*clhs37 + clhs241*clhs53 + clhs244*clhs40 + clhs247*clhs47;
            lhs(0,10)=clhs15*clhs255 - clhs195*clhs34 + clhs218*clhs254 - clhs228*clhs25 + clhs253 - clhs256*clhs31 - clhs265*clhs37 + clhs265*clhs53 + clhs268*clhs40 + clhs271*clhs47;
            lhs(0,11)=0;
            lhs(0,12)=clhs0*clhs276 + clhs0*clhs277 + clhs15*clhs278 - clhs25*clhs279 + clhs275 - clhs280*clhs31 - clhs281*clhs34 - clhs296*clhs37 + clhs296*clhs53 + clhs299*clhs40 + clhs302*clhs47;
            lhs(0,13)=clhs15*clhs308 + clhs226*clhs303 - clhs25*clhs310 - clhs279*clhs34 + clhs307 - clhs309*clhs31 - clhs322*clhs37 + clhs322*clhs53 + clhs325*clhs40 + clhs328*clhs47;
            lhs(0,14)=clhs15*clhs334 - clhs25*clhs309 + clhs254*clhs303 - clhs280*clhs34 - clhs31*clhs335 + clhs333 - clhs344*clhs37 + clhs344*clhs53 + clhs347*clhs40 + clhs350*clhs47;
            lhs(0,15)=0;
            lhs(1,0)=-clhs17*clhs69 - clhs17*clhs74 - clhs33*clhs75 + clhs351*clhs58 + clhs353*clhs63 + clhs354*clhs39 + clhs355*clhs52 + clhs356*clhs46 - clhs37*clhs46 + clhs60;
            lhs(1,1)=-clhs17*clhs77 - clhs17*clhs78 + clhs353*clhs358 + clhs354*clhs76 + clhs355*clhs82 + clhs356*clhs79 + clhs357*clhs44 + clhs36*clhs44 - clhs37*clhs79 + clhs4*clhs44 + clhs44*clhs9 - clhs61*clhs7 - clhs64*clhs68;
            lhs(1,2)=clhs101*clhs354 + clhs104*clhs356 - clhs104*clhs37 + clhs107*clhs355 - clhs17*clhs80 - clhs17*clhs81 + clhs353*clhs366 + clhs364 + clhs365*clhs58 - clhs73*clhs92;
            lhs(1,3)=0;
            lhs(1,4)=-clhs115*clhs69 - clhs115*clhs74 - clhs116*clhs75 + clhs131*clhs354 + clhs134*clhs356 - clhs134*clhs37 + clhs137*clhs355 + clhs142*clhs351 + clhs144 + clhs145*clhs353;
            lhs(1,5)=clhs113*clhs44 - clhs115*clhs77 - clhs115*clhs78 - clhs146*clhs68 + clhs158*clhs354 + clhs161*clhs356 - clhs161*clhs37 + clhs164*clhs355 + clhs353*clhs369 + clhs367 + clhs368*clhs44;
            lhs(1,6)=-clhs115*clhs80 - clhs115*clhs81 + clhs142*clhs365 - clhs171*clhs73 + clhs180*clhs354 + clhs183*clhs356 - clhs183*clhs37 + clhs186*clhs355 + clhs353*clhs376 + clhs375;
            lhs(1,7)=0;
            lhs(1,8)=-clhs194*clhs68 - clhs195*clhs73 - clhs196*clhs75 + clhs211*clhs354 + clhs214*clhs356 - clhs214*clhs37 + clhs217*clhs355 + clhs225 + clhs227*clhs353 + clhs351*clhs377;
            lhs(1,9)=clhs192*clhs44 - clhs194*clhs75 - clhs228*clhs73 - clhs229*clhs68 + clhs241*clhs354 + clhs244*clhs356 - clhs244*clhs37 + clhs247*clhs355 + clhs353*clhs380 + clhs378 + clhs379*clhs44;
            lhs(1,10)=-clhs195*clhs75 - clhs228*clhs68 - clhs256*clhs73 + clhs265*clhs354 + clhs268*clhs356 - clhs268*clhs37 + clhs271*clhs355 + clhs353*clhs387 + clhs365*clhs377 + clhs386;
            lhs(1,11)=0;
            lhs(1,12)=-clhs279*clhs68 - clhs280*clhs73 - clhs281*clhs75 + clhs296*clhs354 + clhs299*clhs356 - clhs299*clhs37 + clhs302*clhs355 + clhs307 + clhs308*clhs353 + clhs351*clhs388;
            lhs(1,13)=clhs277*clhs44 - clhs279*clhs75 - clhs309*clhs73 - clhs310*clhs68 + clhs322*clhs354 + clhs325*clhs356 - clhs325*clhs37 + clhs328*clhs355 + clhs353*clhs391 + clhs389 + clhs390*clhs44;
            lhs(1,14)=-clhs280*clhs75 - clhs309*clhs68 - clhs335*clhs73 + clhs344*clhs354 + clhs347*clhs356 - clhs347*clhs37 + clhs350*clhs355 + clhs353*clhs397 + clhs365*clhs388 + clhs396;
            lhs(1,15)=0;
            lhs(2,0)=-clhs100*clhs33 - clhs17*clhs96 - clhs17*clhs99 + clhs351*clhs86 - clhs37*clhs52 + clhs39*clhs399 + clhs398*clhs91 + clhs400*clhs46 + clhs401*clhs52 + clhs88;
            lhs(2,1)=-clhs102*clhs17 - clhs103*clhs17 + clhs364 + clhs366*clhs398 - clhs37*clhs82 + clhs399*clhs76 + clhs400*clhs79 + clhs401*clhs82 + clhs402*clhs86 - clhs64*clhs95;
            lhs(2,2)=clhs101*clhs399 + clhs104*clhs400 - clhs105*clhs17 - clhs106*clhs17 - clhs107*clhs37 + clhs107*clhs401 + clhs357*clhs50 + clhs36*clhs50 + clhs398*clhs403 + clhs4*clhs50 + clhs50*clhs8 - clhs7*clhs89 - clhs92*clhs98;
            lhs(2,3)=0;
            lhs(2,4)=-clhs100*clhs116 - clhs115*clhs96 - clhs115*clhs99 + clhs131*clhs399 + clhs134*clhs400 - clhs137*clhs37 + clhs137*clhs401 + clhs167*clhs351 + clhs169 + clhs170*clhs398;
            lhs(2,5)=-clhs102*clhs115 - clhs103*clhs115 - clhs146*clhs95 + clhs158*clhs399 + clhs161*clhs400 - clhs164*clhs37 + clhs164*clhs401 + clhs167*clhs402 + clhs375 + clhs376*clhs398;
            lhs(2,6)=-clhs105*clhs115 - clhs106*clhs115 + clhs112*clhs50 - clhs171*clhs98 + clhs180*clhs399 + clhs183*clhs400 - clhs186*clhs37 + clhs186*clhs401 + clhs368*clhs50 + clhs398*clhs405 + clhs404;
            lhs(2,7)=0;
            lhs(2,8)=N[2]*clhs406 - clhs100*clhs196 - clhs194*clhs95 - clhs195*clhs98 + clhs211*clhs399 + clhs214*clhs400 - clhs217*clhs37 + clhs217*clhs401 + clhs253 + clhs255*clhs398;
            lhs(2,9)=N[2]*clhs407 - clhs100*clhs194 - clhs228*clhs98 - clhs229*clhs95 + clhs241*clhs399 + clhs244*clhs400 - clhs247*clhs37 + clhs247*clhs401 + clhs386 + clhs387*clhs398;
            lhs(2,10)=-clhs100*clhs195 + clhs191*clhs50 - clhs228*clhs95 - clhs256*clhs98 + clhs265*clhs399 + clhs268*clhs400 - clhs271*clhs37 + clhs271*clhs401 + clhs379*clhs50 + clhs398*clhs409 + clhs408;
            lhs(2,11)=0;
            lhs(2,12)=N[3]*clhs406 - clhs100*clhs281 - clhs279*clhs95 - clhs280*clhs98 + clhs296*clhs399 + clhs299*clhs400 - clhs302*clhs37 + clhs302*clhs401 + clhs333 + clhs334*clhs398;
            lhs(2,13)=N[3]*clhs407 - clhs100*clhs279 - clhs309*clhs98 - clhs310*clhs95 + clhs322*clhs399 + clhs325*clhs400 - clhs328*clhs37 + clhs328*clhs401 + clhs396 + clhs397*clhs398;
            lhs(2,14)=-clhs100*clhs280 + clhs276*clhs50 - clhs309*clhs95 - clhs335*clhs98 + clhs344*clhs399 + clhs347*clhs400 - clhs350*clhs37 + clhs350*clhs401 + clhs390*clhs50 + clhs398*clhs411 + clhs410;
            lhs(2,15)=0;
            lhs(3,0)=-clhs412*clhs54;
            lhs(3,1)=-clhs1*clhs413;
            lhs(3,2)=-clhs1*clhs414;
            lhs(3,3)=0;
            lhs(3,4)=clhs416;
            lhs(3,5)=clhs417;
            lhs(3,6)=clhs418;
            lhs(3,7)=0;
            lhs(3,8)=clhs420;
            lhs(3,9)=clhs421;
            lhs(3,10)=clhs423;
            lhs(3,11)=0;
            lhs(3,12)=clhs424;
            lhs(3,13)=clhs425;
            lhs(3,14)=clhs426;
            lhs(3,15)=0;
            lhs(4,0)=clhs0*clhs427 + clhs0*clhs428 + clhs111 - clhs124*clhs17 - clhs129*clhs17 - clhs130*clhs33 + clhs16*clhs434 - clhs39*clhs435 + clhs39*clhs438 + clhs436*clhs46 + clhs437*clhs52;
            lhs(4,1)=N[0]*clhs441 - clhs123*clhs64 - clhs132*clhs17 - clhs133*clhs17 + clhs434*clhs63 - clhs435*clhs76 + clhs436*clhs79 + clhs437*clhs82 + clhs438*clhs76 + clhs440;
            lhs(4,2)=N[0]*clhs444 - clhs101*clhs435 + clhs101*clhs438 + clhs104*clhs436 + clhs107*clhs437 - clhs128*clhs92 - clhs135*clhs17 - clhs136*clhs17 + clhs434*clhs91 + clhs443;
            lhs(4,3)=0;
            lhs(4,4)=clhs0*clhs446 + clhs0*clhs448 + clhs0*clhs449 + clhs0*clhs450 + clhs114*clhs434 - clhs115*clhs124 - clhs115*clhs129 - clhs116*clhs130 - clhs131*clhs435 + clhs131*clhs438 + clhs134*clhs436 + clhs137*clhs437 - clhs447*clhs5;
            lhs(4,5)=N[1]*clhs441 - clhs115*clhs132 - clhs115*clhs133 - clhs123*clhs146 + clhs145*clhs434 - clhs158*clhs435 + clhs158*clhs438 + clhs161*clhs436 + clhs164*clhs437 + clhs452;
            lhs(4,6)=N[1]*clhs444 - clhs115*clhs135 - clhs115*clhs136 - clhs128*clhs171 + clhs170*clhs434 - clhs180*clhs435 + clhs180*clhs438 + clhs183*clhs436 + clhs186*clhs437 + clhs453;
            lhs(4,7)=0;
            lhs(4,8)=clhs0*clhs458 + clhs0*clhs459 - clhs123*clhs194 - clhs128*clhs195 - clhs130*clhs196 + clhs193*clhs434 - clhs211*clhs435 + clhs211*clhs438 + clhs214*clhs436 + clhs217*clhs437 + clhs457;
            lhs(4,9)=-clhs123*clhs229 - clhs128*clhs228 - clhs130*clhs194 + clhs218*clhs468 + clhs227*clhs434 - clhs241*clhs435 + clhs241*clhs438 + clhs244*clhs436 + clhs247*clhs437 + clhs467;
            lhs(4,10)=-clhs123*clhs228 - clhs128*clhs256 - clhs130*clhs195 + clhs218*clhs474 + clhs255*clhs434 - clhs265*clhs435 + clhs265*clhs438 + clhs268*clhs436 + clhs271*clhs437 + clhs473;
            lhs(4,11)=0;
            lhs(4,12)=clhs0*clhs479 + clhs0*clhs480 - clhs123*clhs279 - clhs128*clhs280 - clhs130*clhs281 + clhs278*clhs434 - clhs296*clhs435 + clhs296*clhs438 + clhs299*clhs436 + clhs302*clhs437 + clhs478;
            lhs(4,13)=-clhs123*clhs310 - clhs128*clhs309 - clhs130*clhs279 + clhs303*clhs468 + clhs308*clhs434 - clhs322*clhs435 + clhs322*clhs438 + clhs325*clhs436 + clhs328*clhs437 + clhs485;
            lhs(4,14)=-clhs123*clhs309 - clhs128*clhs335 - clhs130*clhs280 + clhs303*clhs474 + clhs334*clhs434 - clhs344*clhs435 + clhs344*clhs438 + clhs347*clhs436 + clhs350*clhs437 + clhs489;
            lhs(4,15)=0;
            lhs(5,0)=-clhs151*clhs17 - clhs156*clhs17 - clhs157*clhs33 + clhs39*clhs493 - clhs435*clhs46 + clhs440 + clhs46*clhs495 + clhs490*clhs58 + clhs492*clhs63 + clhs494*clhs52;
            lhs(5,1)=-clhs150*clhs64 - clhs159*clhs17 - clhs160*clhs17 + clhs358*clhs492 + clhs367 + clhs428*clhs44 - clhs435*clhs79 + clhs44*clhs496 + clhs493*clhs76 + clhs494*clhs82 + clhs495*clhs79;
            lhs(5,2)=clhs101*clhs493 - clhs104*clhs435 + clhs104*clhs495 + clhs107*clhs494 - clhs155*clhs92 - clhs162*clhs17 - clhs163*clhs17 + clhs366*clhs492 + clhs498 + clhs499*clhs58;
            lhs(5,3)=0;
            lhs(5,4)=-clhs115*clhs151 - clhs115*clhs156 - clhs116*clhs157 + clhs131*clhs493 - clhs134*clhs435 + clhs134*clhs495 + clhs137*clhs494 + clhs142*clhs490 + clhs145*clhs492 + clhs452;
            lhs(5,5)=-clhs115*clhs159 - clhs115*clhs160 - clhs146*clhs150 + clhs158*clhs493 - clhs161*clhs435 + clhs161*clhs495 + clhs164*clhs494 + clhs369*clhs492 + clhs44*clhs446 + clhs44*clhs449 + clhs44*clhs450 + clhs44*clhs500 - clhs447*clhs61;
            lhs(5,6)=-clhs115*clhs162 - clhs115*clhs163 + clhs142*clhs499 - clhs155*clhs171 + clhs180*clhs493 - clhs183*clhs435 + clhs183*clhs495 + clhs186*clhs494 + clhs376*clhs492 + clhs501;
            lhs(5,7)=0;
            lhs(5,8)=-clhs150*clhs194 - clhs155*clhs195 - clhs157*clhs196 + clhs211*clhs493 - clhs214*clhs435 + clhs214*clhs495 + clhs217*clhs494 + clhs227*clhs492 + clhs377*clhs490 + clhs467;
            lhs(5,9)=-clhs150*clhs229 - clhs155*clhs228 - clhs157*clhs194 + clhs241*clhs493 - clhs244*clhs435 + clhs244*clhs495 + clhs247*clhs494 + clhs380*clhs492 + clhs44*clhs459 + clhs44*clhs503 + clhs502;
            lhs(5,10)=-clhs150*clhs228 - clhs155*clhs256 - clhs157*clhs195 + clhs265*clhs493 - clhs268*clhs435 + clhs268*clhs495 + clhs271*clhs494 + clhs377*clhs499 + clhs387*clhs492 + clhs510;
            lhs(5,11)=0;
            lhs(5,12)=-clhs150*clhs279 - clhs155*clhs280 - clhs157*clhs281 + clhs296*clhs493 - clhs299*clhs435 + clhs299*clhs495 + clhs302*clhs494 + clhs308*clhs492 + clhs388*clhs490 + clhs485;
            lhs(5,13)=-clhs150*clhs310 - clhs155*clhs309 - clhs157*clhs279 + clhs322*clhs493 - clhs325*clhs435 + clhs325*clhs495 + clhs328*clhs494 + clhs391*clhs492 + clhs44*clhs480 + clhs44*clhs512 + clhs511;
            lhs(5,14)=-clhs150*clhs309 - clhs155*clhs335 - clhs157*clhs280 + clhs344*clhs493 - clhs347*clhs435 + clhs347*clhs495 + clhs350*clhs494 + clhs388*clhs499 + clhs397*clhs492 + clhs516;
            lhs(5,15)=0;
            lhs(6,0)=-clhs17*clhs175 - clhs17*clhs178 - clhs179*clhs33 + clhs39*clhs518 - clhs435*clhs52 + clhs443 + clhs46*clhs519 + clhs490*clhs86 + clhs517*clhs91 + clhs52*clhs520;
            lhs(6,1)=-clhs17*clhs181 - clhs17*clhs182 - clhs174*clhs64 + clhs366*clhs517 - clhs435*clhs82 + clhs498 + clhs518*clhs76 + clhs519*clhs79 + clhs520*clhs82 + clhs521*clhs86;
            lhs(6,2)=clhs101*clhs518 + clhs104*clhs519 - clhs107*clhs435 + clhs107*clhs520 - clhs17*clhs184 - clhs17*clhs185 - clhs177*clhs92 + clhs403*clhs517 + clhs404 + clhs427*clhs50 + clhs496*clhs50;
            lhs(6,3)=0;
            lhs(6,4)=-clhs115*clhs175 - clhs115*clhs178 - clhs116*clhs179 + clhs131*clhs518 + clhs134*clhs519 - clhs137*clhs435 + clhs137*clhs520 + clhs167*clhs490 + clhs170*clhs517 + clhs453;
            lhs(6,5)=-clhs115*clhs181 - clhs115*clhs182 - clhs146*clhs174 + clhs158*clhs518 + clhs161*clhs519 - clhs164*clhs435 + clhs164*clhs520 + clhs167*clhs521 + clhs376*clhs517 + clhs501;
            lhs(6,6)=-clhs115*clhs184 - clhs115*clhs185 - clhs171*clhs177 + clhs180*clhs518 + clhs183*clhs519 - clhs186*clhs435 + clhs186*clhs520 + clhs405*clhs517 + clhs446*clhs50 - clhs447*clhs89 + clhs448*clhs50 + clhs450*clhs50 + clhs50*clhs500;
            lhs(6,7)=0;
            lhs(6,8)=N[2]*clhs522 - clhs174*clhs194 - clhs177*clhs195 - clhs179*clhs196 + clhs211*clhs518 + clhs214*clhs519 - clhs217*clhs435 + clhs217*clhs520 + clhs255*clhs517 + clhs473;
            lhs(6,9)=N[2]*clhs523 - clhs174*clhs229 - clhs177*clhs228 - clhs179*clhs194 + clhs241*clhs518 + clhs244*clhs519 - clhs247*clhs435 + clhs247*clhs520 + clhs387*clhs517 + clhs510;
            lhs(6,10)=-clhs174*clhs228 - clhs177*clhs256 - clhs179*clhs195 + clhs265*clhs518 + clhs268*clhs519 - clhs271*clhs435 + clhs271*clhs520 + clhs409*clhs517 + clhs458*clhs50 + clhs50*clhs503 + clhs524;
            lhs(6,11)=0;
            lhs(6,12)=N[3]*clhs522 - clhs174*clhs279 - clhs177*clhs280 - clhs179*clhs281 + clhs296*clhs518 + clhs299*clhs519 - clhs302*clhs435 + clhs302*clhs520 + clhs334*clhs517 + clhs489;
            lhs(6,13)=N[3]*clhs523 - clhs174*clhs310 - clhs177*clhs309 - clhs179*clhs279 + clhs322*clhs518 + clhs325*clhs519 - clhs328*clhs435 + clhs328*clhs520 + clhs397*clhs517 + clhs516;
            lhs(6,14)=-clhs174*clhs309 - clhs177*clhs335 - clhs179*clhs280 + clhs344*clhs518 + clhs347*clhs519 - clhs350*clhs435 + clhs350*clhs520 + clhs411*clhs517 + clhs479*clhs50 + clhs50*clhs512 + clhs525;
            lhs(6,15)=0;
            lhs(7,0)=clhs416;
            lhs(7,1)=clhs417;
            lhs(7,2)=clhs418;
            lhs(7,3)=0;
            lhs(7,4)=-clhs412*clhs451;
            lhs(7,5)=-clhs413*clhs445;
            lhs(7,6)=-clhs414*clhs445;
            lhs(7,7)=0;
            lhs(7,8)=clhs527;
            lhs(7,9)=clhs528;
            lhs(7,10)=clhs530;
            lhs(7,11)=0;
            lhs(7,12)=clhs531;
            lhs(7,13)=clhs532;
            lhs(7,14)=clhs533;
            lhs(7,15)=0;
            lhs(8,0)=clhs0*clhs534 + clhs0*clhs535 + clhs16*clhs541 - clhs17*clhs204 - clhs17*clhs209 + clhs190 - clhs210*clhs33 - clhs39*clhs542 + clhs39*clhs545 + clhs46*clhs543 + clhs52*clhs544;
            lhs(8,1)=N[0]*clhs548 - clhs17*clhs212 - clhs17*clhs213 - clhs203*clhs64 + clhs541*clhs63 - clhs542*clhs76 + clhs543*clhs79 + clhs544*clhs82 + clhs545*clhs76 + clhs547;
            lhs(8,2)=N[0]*clhs551 - clhs101*clhs542 + clhs101*clhs545 + clhs104*clhs543 + clhs107*clhs544 - clhs17*clhs215 - clhs17*clhs216 - clhs208*clhs92 + clhs541*clhs91 + clhs550;
            lhs(8,3)=0;
            lhs(8,4)=clhs0*clhs552 + clhs0*clhs553 + clhs114*clhs541 - clhs115*clhs204 - clhs115*clhs209 - clhs116*clhs210 - clhs131*clhs542 + clhs131*clhs545 + clhs134*clhs543 + clhs137*clhs544 + clhs457;
            lhs(8,5)=N[1]*clhs548 - clhs115*clhs212 - clhs115*clhs213 + clhs145*clhs541 - clhs146*clhs203 - clhs158*clhs542 + clhs158*clhs545 + clhs161*clhs543 + clhs164*clhs544 + clhs554;
            lhs(8,6)=N[1]*clhs551 - clhs115*clhs215 - clhs115*clhs216 + clhs170*clhs541 - clhs171*clhs208 - clhs180*clhs542 + clhs180*clhs545 + clhs183*clhs543 + clhs186*clhs544 + clhs555;
            lhs(8,7)=0;
            lhs(8,8)=clhs0*clhs557 + clhs0*clhs559 + clhs0*clhs560 + clhs0*clhs561 + clhs193*clhs541 - clhs194*clhs203 - clhs195*clhs208 - clhs196*clhs210 - clhs211*clhs542 + clhs211*clhs545 + clhs214*clhs543 + clhs217*clhs544 - clhs5*clhs558;
            lhs(8,9)=-clhs194*clhs210 - clhs203*clhs229 - clhs208*clhs228 + clhs218*clhs565 + clhs227*clhs541 - clhs241*clhs542 + clhs241*clhs545 + clhs244*clhs543 + clhs247*clhs544 + clhs564;
            lhs(8,10)=-clhs195*clhs210 - clhs203*clhs228 - clhs208*clhs256 + clhs218*clhs568 + clhs255*clhs541 - clhs265*clhs542 + clhs265*clhs545 + clhs268*clhs543 + clhs271*clhs544 + clhs567;
            lhs(8,11)=0;
            lhs(8,12)=clhs0*clhs573 + clhs0*clhs574 - clhs203*clhs279 - clhs208*clhs280 - clhs210*clhs281 + clhs278*clhs541 - clhs296*clhs542 + clhs296*clhs545 + clhs299*clhs543 + clhs302*clhs544 + clhs572;
            lhs(8,13)=-clhs203*clhs310 - clhs208*clhs309 - clhs210*clhs279 + clhs303*clhs565 + clhs308*clhs541 - clhs322*clhs542 + clhs322*clhs545 + clhs325*clhs543 + clhs328*clhs544 + clhs579;
            lhs(8,14)=-clhs203*clhs309 - clhs208*clhs335 - clhs210*clhs280 + clhs303*clhs568 + clhs334*clhs541 - clhs344*clhs542 + clhs344*clhs545 + clhs347*clhs543 + clhs350*clhs544 + clhs583;
            lhs(8,15)=0;
            lhs(9,0)=-clhs17*clhs234 - clhs17*clhs239 - clhs240*clhs33 + clhs39*clhs587 - clhs46*clhs542 + clhs46*clhs589 + clhs52*clhs588 + clhs547 + clhs58*clhs584 + clhs586*clhs63;
            lhs(9,1)=-clhs17*clhs242 - clhs17*clhs243 - clhs233*clhs64 + clhs358*clhs586 + clhs378 + clhs44*clhs535 + clhs44*clhs590 - clhs542*clhs79 + clhs587*clhs76 + clhs588*clhs82 + clhs589*clhs79;
            lhs(9,2)=clhs101*clhs587 - clhs104*clhs542 + clhs104*clhs589 + clhs107*clhs588 - clhs17*clhs245 - clhs17*clhs246 - clhs238*clhs92 + clhs366*clhs586 + clhs58*clhs593 + clhs592;
            lhs(9,3)=0;
            lhs(9,4)=-clhs115*clhs234 - clhs115*clhs239 - clhs116*clhs240 + clhs131*clhs587 - clhs134*clhs542 + clhs134*clhs589 + clhs137*clhs588 + clhs142*clhs584 + clhs145*clhs586 + clhs554;
            lhs(9,5)=-clhs115*clhs242 - clhs115*clhs243 - clhs146*clhs233 + clhs158*clhs587 - clhs161*clhs542 + clhs161*clhs589 + clhs164*clhs588 + clhs369*clhs586 + clhs44*clhs553 + clhs44*clhs594 + clhs502;
            lhs(9,6)=-clhs115*clhs245 - clhs115*clhs246 + clhs142*clhs593 - clhs171*clhs238 + clhs180*clhs587 - clhs183*clhs542 + clhs183*clhs589 + clhs186*clhs588 + clhs376*clhs586 + clhs595;
            lhs(9,7)=0;
            lhs(9,8)=-clhs194*clhs233 - clhs195*clhs238 - clhs196*clhs240 + clhs211*clhs587 - clhs214*clhs542 + clhs214*clhs589 + clhs217*clhs588 + clhs227*clhs586 + clhs377*clhs584 + clhs564;
            lhs(9,9)=-clhs194*clhs240 - clhs228*clhs238 - clhs229*clhs233 + clhs241*clhs587 - clhs244*clhs542 + clhs244*clhs589 + clhs247*clhs588 + clhs380*clhs586 + clhs44*clhs557 + clhs44*clhs560 + clhs44*clhs561 + clhs44*clhs596 - clhs558*clhs61;
            lhs(9,10)=-clhs195*clhs240 - clhs228*clhs233 - clhs238*clhs256 + clhs265*clhs587 - clhs268*clhs542 + clhs268*clhs589 + clhs271*clhs588 + clhs377*clhs593 + clhs387*clhs586 + clhs598;
            lhs(9,11)=0;
            lhs(9,12)=-clhs233*clhs279 - clhs238*clhs280 - clhs240*clhs281 + clhs296*clhs587 - clhs299*clhs542 + clhs299*clhs589 + clhs302*clhs588 + clhs308*clhs586 + clhs388*clhs584 + clhs579;
            lhs(9,13)=-clhs233*clhs310 - clhs238*clhs309 - clhs240*clhs279 + clhs322*clhs587 - clhs325*clhs542 + clhs325*clhs589 + clhs328*clhs588 + clhs391*clhs586 + clhs44*clhs574 + clhs44*clhs600 + clhs599;
            lhs(9,14)=-clhs233*clhs309 - clhs238*clhs335 - clhs240*clhs280 + clhs344*clhs587 - clhs347*clhs542 + clhs347*clhs589 + clhs350*clhs588 + clhs388*clhs593 + clhs397*clhs586 + clhs605;
            lhs(9,15)=0;
            lhs(10,0)=-clhs17*clhs260 - clhs17*clhs263 - clhs264*clhs33 + clhs39*clhs607 + clhs46*clhs608 - clhs52*clhs542 + clhs52*clhs609 + clhs550 + clhs584*clhs86 + clhs606*clhs91;
            lhs(10,1)=-clhs17*clhs266 - clhs17*clhs267 - clhs259*clhs64 + clhs366*clhs606 - clhs542*clhs82 + clhs592 + clhs607*clhs76 + clhs608*clhs79 + clhs609*clhs82 + clhs610*clhs86;
            lhs(10,2)=clhs101*clhs607 + clhs104*clhs608 - clhs107*clhs542 + clhs107*clhs609 - clhs17*clhs269 - clhs17*clhs270 - clhs262*clhs92 + clhs403*clhs606 + clhs408 + clhs50*clhs534 + clhs50*clhs590;
            lhs(10,3)=0;
            lhs(10,4)=-clhs115*clhs260 - clhs115*clhs263 - clhs116*clhs264 + clhs131*clhs607 + clhs134*clhs608 - clhs137*clhs542 + clhs137*clhs609 + clhs167*clhs584 + clhs170*clhs606 + clhs555;
            lhs(10,5)=-clhs115*clhs266 - clhs115*clhs267 - clhs146*clhs259 + clhs158*clhs607 + clhs161*clhs608 - clhs164*clhs542 + clhs164*clhs609 + clhs167*clhs610 + clhs376*clhs606 + clhs595;
            lhs(10,6)=-clhs115*clhs269 - clhs115*clhs270 - clhs171*clhs262 + clhs180*clhs607 + clhs183*clhs608 - clhs186*clhs542 + clhs186*clhs609 + clhs405*clhs606 + clhs50*clhs552 + clhs50*clhs594 + clhs524;
            lhs(10,7)=0;
            lhs(10,8)=N[2]*clhs611 - clhs194*clhs259 - clhs195*clhs262 - clhs196*clhs264 + clhs211*clhs607 + clhs214*clhs608 - clhs217*clhs542 + clhs217*clhs609 + clhs255*clhs606 + clhs567;
            lhs(10,9)=N[2]*clhs612 - clhs194*clhs264 - clhs228*clhs262 - clhs229*clhs259 + clhs241*clhs607 + clhs244*clhs608 - clhs247*clhs542 + clhs247*clhs609 + clhs387*clhs606 + clhs598;
            lhs(10,10)=-clhs195*clhs264 - clhs228*clhs259 - clhs256*clhs262 + clhs265*clhs607 + clhs268*clhs608 - clhs271*clhs542 + clhs271*clhs609 + clhs409*clhs606 + clhs50*clhs557 + clhs50*clhs559 + clhs50*clhs561 + clhs50*clhs596 - clhs558*clhs89;
            lhs(10,11)=0;
            lhs(10,12)=N[3]*clhs611 - clhs259*clhs279 - clhs262*clhs280 - clhs264*clhs281 + clhs296*clhs607 + clhs299*clhs608 - clhs302*clhs542 + clhs302*clhs609 + clhs334*clhs606 + clhs583;
            lhs(10,13)=N[3]*clhs612 - clhs259*clhs310 - clhs262*clhs309 - clhs264*clhs279 + clhs322*clhs607 + clhs325*clhs608 - clhs328*clhs542 + clhs328*clhs609 + clhs397*clhs606 + clhs605;
            lhs(10,14)=-clhs259*clhs309 - clhs262*clhs335 - clhs264*clhs280 + clhs344*clhs607 + clhs347*clhs608 - clhs350*clhs542 + clhs350*clhs609 + clhs411*clhs606 + clhs50*clhs573 + clhs50*clhs600 + clhs613;
            lhs(10,15)=0;
            lhs(11,0)=clhs420;
            lhs(11,1)=clhs421;
            lhs(11,2)=clhs423;
            lhs(11,3)=0;
            lhs(11,4)=clhs527;
            lhs(11,5)=clhs528;
            lhs(11,6)=clhs530;
            lhs(11,7)=0;
            lhs(11,8)=-clhs412*clhs562;
            lhs(11,9)=-clhs413*clhs556;
            lhs(11,10)=-clhs414*clhs556;
            lhs(11,11)=0;
            lhs(11,12)=clhs614;
            lhs(11,13)=clhs615;
            lhs(11,14)=clhs616;
            lhs(11,15)=0;
            lhs(12,0)=clhs0*clhs617 + clhs0*clhs618 + clhs16*clhs624 - clhs17*clhs289 - clhs17*clhs294 + clhs275 - clhs295*clhs33 - clhs39*clhs625 + clhs39*clhs628 + clhs46*clhs626 + clhs52*clhs627;
            lhs(12,1)=N[0]*clhs631 - clhs17*clhs297 - clhs17*clhs298 - clhs288*clhs64 + clhs624*clhs63 - clhs625*clhs76 + clhs626*clhs79 + clhs627*clhs82 + clhs628*clhs76 + clhs630;
            lhs(12,2)=N[0]*clhs634 - clhs101*clhs625 + clhs101*clhs628 + clhs104*clhs626 + clhs107*clhs627 - clhs17*clhs300 - clhs17*clhs301 - clhs293*clhs92 + clhs624*clhs91 + clhs633;
            lhs(12,3)=0;
            lhs(12,4)=clhs0*clhs635 + clhs0*clhs636 + clhs114*clhs624 - clhs115*clhs289 - clhs115*clhs294 - clhs116*clhs295 - clhs131*clhs625 + clhs131*clhs628 + clhs134*clhs626 + clhs137*clhs627 + clhs478;
            lhs(12,5)=N[1]*clhs631 - clhs115*clhs297 - clhs115*clhs298 + clhs145*clhs624 - clhs146*clhs288 - clhs158*clhs625 + clhs158*clhs628 + clhs161*clhs626 + clhs164*clhs627 + clhs637;
            lhs(12,6)=N[1]*clhs634 - clhs115*clhs300 - clhs115*clhs301 + clhs170*clhs624 - clhs171*clhs293 - clhs180*clhs625 + clhs180*clhs628 + clhs183*clhs626 + clhs186*clhs627 + clhs638;
            lhs(12,7)=0;
            lhs(12,8)=clhs0*clhs639 + clhs0*clhs640 + clhs193*clhs624 - clhs194*clhs288 - clhs195*clhs293 - clhs196*clhs295 - clhs211*clhs625 + clhs211*clhs628 + clhs214*clhs626 + clhs217*clhs627 + clhs572;
            lhs(12,9)=-clhs194*clhs295 + clhs218*clhs643 + clhs227*clhs624 - clhs228*clhs293 - clhs229*clhs288 - clhs241*clhs625 + clhs241*clhs628 + clhs244*clhs626 + clhs247*clhs627 + clhs642;
            lhs(12,10)=-clhs195*clhs295 + clhs218*clhs646 - clhs228*clhs288 + clhs255*clhs624 - clhs256*clhs293 - clhs265*clhs625 + clhs265*clhs628 + clhs268*clhs626 + clhs271*clhs627 + clhs645;
            lhs(12,11)=0;
            lhs(12,12)=clhs0*clhs648 + clhs0*clhs650 + clhs0*clhs651 + clhs0*clhs652 + clhs278*clhs624 - clhs279*clhs288 - clhs280*clhs293 - clhs281*clhs295 - clhs296*clhs625 + clhs296*clhs628 + clhs299*clhs626 + clhs302*clhs627 - clhs5*clhs649;
            lhs(12,13)=-clhs279*clhs295 - clhs288*clhs310 - clhs293*clhs309 + clhs303*clhs643 + clhs308*clhs624 - clhs322*clhs625 + clhs322*clhs628 + clhs325*clhs626 + clhs328*clhs627 + clhs654;
            lhs(12,14)=-clhs280*clhs295 - clhs288*clhs309 - clhs293*clhs335 + clhs303*clhs646 + clhs334*clhs624 - clhs344*clhs625 + clhs344*clhs628 + clhs347*clhs626 + clhs350*clhs627 + clhs655;
            lhs(12,15)=0;
            lhs(13,0)=-clhs17*clhs315 - clhs17*clhs320 - clhs321*clhs33 + clhs39*clhs659 - clhs46*clhs625 + clhs46*clhs661 + clhs52*clhs660 + clhs58*clhs656 + clhs63*clhs658 + clhs630;
            lhs(13,1)=-clhs17*clhs323 - clhs17*clhs324 - clhs314*clhs64 + clhs358*clhs658 + clhs389 + clhs44*clhs618 + clhs44*clhs662 - clhs625*clhs79 + clhs659*clhs76 + clhs660*clhs82 + clhs661*clhs79;
            lhs(13,2)=clhs101*clhs659 - clhs104*clhs625 + clhs104*clhs661 + clhs107*clhs660 - clhs17*clhs326 - clhs17*clhs327 - clhs319*clhs92 + clhs366*clhs658 + clhs58*clhs665 + clhs664;
            lhs(13,3)=0;
            lhs(13,4)=-clhs115*clhs315 - clhs115*clhs320 - clhs116*clhs321 + clhs131*clhs659 - clhs134*clhs625 + clhs134*clhs661 + clhs137*clhs660 + clhs142*clhs656 + clhs145*clhs658 + clhs637;
            lhs(13,5)=-clhs115*clhs323 - clhs115*clhs324 - clhs146*clhs314 + clhs158*clhs659 - clhs161*clhs625 + clhs161*clhs661 + clhs164*clhs660 + clhs369*clhs658 + clhs44*clhs636 + clhs44*clhs666 + clhs511;
            lhs(13,6)=-clhs115*clhs326 - clhs115*clhs327 + clhs142*clhs665 - clhs171*clhs319 + clhs180*clhs659 - clhs183*clhs625 + clhs183*clhs661 + clhs186*clhs660 + clhs376*clhs658 + clhs667;
            lhs(13,7)=0;
            lhs(13,8)=-clhs194*clhs314 - clhs195*clhs319 - clhs196*clhs321 + clhs211*clhs659 - clhs214*clhs625 + clhs214*clhs661 + clhs217*clhs660 + clhs227*clhs658 + clhs377*clhs656 + clhs642;
            lhs(13,9)=-clhs194*clhs321 - clhs228*clhs319 - clhs229*clhs314 + clhs241*clhs659 - clhs244*clhs625 + clhs244*clhs661 + clhs247*clhs660 + clhs380*clhs658 + clhs44*clhs640 + clhs44*clhs668 + clhs599;
            lhs(13,10)=-clhs195*clhs321 - clhs228*clhs314 - clhs256*clhs319 + clhs265*clhs659 - clhs268*clhs625 + clhs268*clhs661 + clhs271*clhs660 + clhs377*clhs665 + clhs387*clhs658 + clhs670;
            lhs(13,11)=0;
            lhs(13,12)=-clhs279*clhs314 - clhs280*clhs319 - clhs281*clhs321 + clhs296*clhs659 - clhs299*clhs625 + clhs299*clhs661 + clhs302*clhs660 + clhs308*clhs658 + clhs388*clhs656 + clhs654;
            lhs(13,13)=-clhs279*clhs321 - clhs309*clhs319 - clhs310*clhs314 + clhs322*clhs659 - clhs325*clhs625 + clhs325*clhs661 + clhs328*clhs660 + clhs391*clhs658 + clhs44*clhs648 + clhs44*clhs651 + clhs44*clhs652 + clhs44*clhs671 - clhs61*clhs649;
            lhs(13,14)=-clhs280*clhs321 - clhs309*clhs314 - clhs319*clhs335 + clhs344*clhs659 - clhs347*clhs625 + clhs347*clhs661 + clhs350*clhs660 + clhs388*clhs665 + clhs397*clhs658 + clhs672;
            lhs(13,15)=0;
            lhs(14,0)=-clhs17*clhs339 - clhs17*clhs342 - clhs33*clhs343 + clhs39*clhs674 + clhs46*clhs675 - clhs52*clhs625 + clhs52*clhs676 + clhs633 + clhs656*clhs86 + clhs673*clhs91;
            lhs(14,1)=-clhs17*clhs345 - clhs17*clhs346 - clhs338*clhs64 + clhs366*clhs673 - clhs625*clhs82 + clhs664 + clhs674*clhs76 + clhs675*clhs79 + clhs676*clhs82 + clhs677*clhs86;
            lhs(14,2)=clhs101*clhs674 + clhs104*clhs675 - clhs107*clhs625 + clhs107*clhs676 - clhs17*clhs348 - clhs17*clhs349 - clhs341*clhs92 + clhs403*clhs673 + clhs410 + clhs50*clhs617 + clhs50*clhs662;
            lhs(14,3)=0;
            lhs(14,4)=-clhs115*clhs339 - clhs115*clhs342 - clhs116*clhs343 + clhs131*clhs674 + clhs134*clhs675 - clhs137*clhs625 + clhs137*clhs676 + clhs167*clhs656 + clhs170*clhs673 + clhs638;
            lhs(14,5)=-clhs115*clhs345 - clhs115*clhs346 - clhs146*clhs338 + clhs158*clhs674 + clhs161*clhs675 - clhs164*clhs625 + clhs164*clhs676 + clhs167*clhs677 + clhs376*clhs673 + clhs667;
            lhs(14,6)=-clhs115*clhs348 - clhs115*clhs349 - clhs171*clhs341 + clhs180*clhs674 + clhs183*clhs675 - clhs186*clhs625 + clhs186*clhs676 + clhs405*clhs673 + clhs50*clhs635 + clhs50*clhs666 + clhs525;
            lhs(14,7)=0;
            lhs(14,8)=N[2]*clhs678 - clhs194*clhs338 - clhs195*clhs341 - clhs196*clhs343 + clhs211*clhs674 + clhs214*clhs675 - clhs217*clhs625 + clhs217*clhs676 + clhs255*clhs673 + clhs645;
            lhs(14,9)=N[2]*clhs679 - clhs194*clhs343 - clhs228*clhs341 - clhs229*clhs338 + clhs241*clhs674 + clhs244*clhs675 - clhs247*clhs625 + clhs247*clhs676 + clhs387*clhs673 + clhs670;
            lhs(14,10)=-clhs195*clhs343 - clhs228*clhs338 - clhs256*clhs341 + clhs265*clhs674 + clhs268*clhs675 - clhs271*clhs625 + clhs271*clhs676 + clhs409*clhs673 + clhs50*clhs639 + clhs50*clhs668 + clhs613;
            lhs(14,11)=0;
            lhs(14,12)=N[3]*clhs678 - clhs279*clhs338 - clhs280*clhs341 - clhs281*clhs343 + clhs296*clhs674 + clhs299*clhs675 - clhs302*clhs625 + clhs302*clhs676 + clhs334*clhs673 + clhs655;
            lhs(14,13)=N[3]*clhs679 - clhs279*clhs343 - clhs309*clhs341 - clhs310*clhs338 + clhs322*clhs674 + clhs325*clhs675 - clhs328*clhs625 + clhs328*clhs676 + clhs397*clhs673 + clhs672;
            lhs(14,14)=-clhs280*clhs343 - clhs309*clhs338 - clhs335*clhs341 + clhs344*clhs674 + clhs347*clhs675 - clhs350*clhs625 + clhs350*clhs676 + clhs411*clhs673 + clhs50*clhs648 + clhs50*clhs650 + clhs50*clhs652 + clhs50*clhs671 - clhs649*clhs89;
            lhs(14,15)=0;
            lhs(15,0)=clhs424;
            lhs(15,1)=clhs425;
            lhs(15,2)=clhs426;
            lhs(15,3)=0;
            lhs(15,4)=clhs531;
            lhs(15,5)=clhs532;
            lhs(15,6)=clhs533;
            lhs(15,7)=0;
            lhs(15,8)=clhs614;
            lhs(15,9)=clhs615;
            lhs(15,10)=clhs616;
            lhs(15,11)=0;
            lhs(15,12)=-clhs412*clhs653;
            lhs(15,13)=-clhs413*clhs647;
            lhs(15,14)=-clhs414*clhs647;
            lhs(15,15)=0;


    return lhs;
}

template <class TBaseElement>
typename EmbeddedFluidElementDiscontinuous<TBaseElement>::VectorType EmbeddedFluidElementDiscontinuous<TBaseElement>::NitscheTermsRHS2D(
    EmbeddedDiscontinuousElementData& rData)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;

    const double dt = rData.DeltaTime;

    const auto& v = rData.Velocity;
    const auto& vmesh = rData.MeshVelocity;
    const auto& vconv = v - vmesh;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // If there is embedded velocity, substract it to the previous iteration solution
    array_1d<double,3> g = ZeroVector(3);
    if (this->Has(EMBEDDED_VELOCITY)){
        g = this->GetValue(EMBEDDED_VELOCITY);
    }

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;
    const auto& normal = rData.Normal;

    // Get Nitsche imposition values
    const double eps = rData.SlipLength;;
    const double gamma = rData.PenaltyCoefficient;
    const double adjoint = 1.0;

    VectorType rhs = ZeroVector(9);

    const double crhs0 =             1.0/gamma;
const double crhs1 =             1.0/h;
const double crhs2 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) - g[0];
const double crhs3 =             crhs2*normal[0];
const double crhs4 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) - g[1];
const double crhs5 =             crhs4*normal[1];
const double crhs6 =             crhs3 + crhs5;
const double crhs7 =             N[0]*crhs0*crhs1*crhs6*mu;
const double crhs8 =             1.0/(eps + gamma*h);
const double crhs9 =             N[0]*crhs8*mu;
const double crhs10 =             pow(normal[0], 2);
const double crhs11 =             -crhs10 + 1;
const double crhs12 =             -crhs11*crhs2 + crhs5*normal[0];
const double crhs13 =             C(0,2)*DN(0,0);
const double crhs14 =             C(2,2)*DN(0,1) + crhs13;
const double crhs15 =             adjoint*crhs6*normal[0];
const double crhs16 =             C(1,2)*DN(0,1);
const double crhs17 =             adjoint*crhs6*normal[1];
const double crhs18 =             pow(normal[1], 2);
const double crhs19 =             -crhs18 + 1;
const double crhs20 =             -crhs19*crhs4 + crhs3*normal[1];
const double crhs21 =             1.0*adjoint*crhs20*crhs8*gamma*h*mu*normal[0];
const double crhs22 =             DN(0,0)*normal[0];
const double crhs23 =             DN(0,1)*normal[1];
const double crhs24 =             2.0*crhs22 + 1.0*crhs23;
const double crhs25 =             adjoint*crhs12*crhs8*gamma*h*mu;
const double crhs26 =             h*rho*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1), 2)) + mu + pow(h, 2)*rho/dt;
const double crhs27 =             N[0]*crhs0*crhs1*crhs26*crhs6;
const double crhs28 =             N[0]*crhs8*eps;
const double crhs29 =             normal[0]*normal[1];
const double crhs30 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs31 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs32 =             DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0);
const double crhs33 =             C(0,2)*crhs30 + C(1,2)*crhs31 + C(2,2)*crhs32;
const double crhs34 =             crhs33*normal[0] + normal[1]*(C(0,1)*crhs30 + C(1,1)*crhs31 + C(1,2)*crhs32);
const double crhs35 =             crhs33*normal[1] + normal[0]*(C(0,0)*crhs30 + C(0,1)*crhs31 + C(0,2)*crhs32);
const double crhs36 =             -crhs11*crhs35 + crhs29*crhs34;
const double crhs37 =             -crhs19*crhs34 + crhs29*crhs35;
const double crhs38 =             1.0*adjoint*crhs37*crhs8*eps*gamma*h*normal[0];
const double crhs39 =             adjoint*crhs36*crhs8*eps*gamma*h;
const double crhs40 =             C(2,2)*DN(0,0) + crhs16;
const double crhs41 =             1.0*adjoint*crhs12*crhs8*gamma*h*mu*normal[1];
const double crhs42 =             1.0*crhs22 + 2.0*crhs23;
const double crhs43 =             adjoint*crhs20*crhs8*gamma*h*mu;
const double crhs44 =             1.0*adjoint*crhs36*crhs8*eps*gamma*h*normal[1];
const double crhs45 =             adjoint*crhs37*crhs8*eps*gamma*h;
const double crhs46 =             crhs6*(crhs10 + crhs18);
const double crhs47 =             N[1]*crhs0*crhs1*crhs6*mu;
const double crhs48 =             N[1]*crhs8*mu;
const double crhs49 =             C(0,2)*DN(1,0);
const double crhs50 =             C(2,2)*DN(1,1) + crhs49;
const double crhs51 =             C(1,2)*DN(1,1);
const double crhs52 =             DN(1,0)*normal[0];
const double crhs53 =             DN(1,1)*normal[1];
const double crhs54 =             2.0*crhs52 + 1.0*crhs53;
const double crhs55 =             N[1]*crhs0*crhs1*crhs26*crhs6;
const double crhs56 =             N[1]*crhs8*eps;
const double crhs57 =             C(2,2)*DN(1,0) + crhs51;
const double crhs58 =             1.0*crhs52 + 2.0*crhs53;
const double crhs59 =             N[2]*crhs0*crhs1*crhs6*mu;
const double crhs60 =             N[2]*crhs8*mu;
const double crhs61 =             C(0,2)*DN(2,0);
const double crhs62 =             C(2,2)*DN(2,1) + crhs61;
const double crhs63 =             C(1,2)*DN(2,1);
const double crhs64 =             DN(2,0)*normal[0];
const double crhs65 =             DN(2,1)*normal[1];
const double crhs66 =             2.0*crhs64 + 1.0*crhs65;
const double crhs67 =             N[2]*crhs0*crhs1*crhs26*crhs6;
const double crhs68 =             N[2]*crhs8*eps;
const double crhs69 =             C(2,2)*DN(2,0) + crhs63;
const double crhs70 =             1.0*crhs64 + 2.0*crhs65;
            rhs[0]=-DN(0,1)*crhs21 - DN(0,1)*crhs38 + crhs12*crhs9 + crhs15*(crhs14*normal[1] + normal[0]*(C(0,0)*DN(0,0) + C(0,2)*DN(0,1))) + crhs17*(crhs14*normal[0] + normal[1]*(C(0,1)*DN(0,0) + crhs16)) - crhs24*crhs25 - crhs24*crhs39 - crhs27*normal[0] + crhs28*crhs36 - crhs7*normal[0];
            rhs[1]=-DN(0,0)*crhs41 - DN(0,0)*crhs44 + crhs15*(crhs40*normal[1] + normal[0]*(C(0,1)*DN(0,1) + crhs13)) + crhs17*(crhs40*normal[0] + normal[1]*(C(1,1)*DN(0,1) + C(1,2)*DN(0,0))) + crhs20*crhs9 - crhs27*normal[1] + crhs28*crhs37 - crhs42*crhs43 - crhs42*crhs45 - crhs7*normal[1];
            rhs[2]=N[0]*crhs46;
            rhs[3]=-DN(1,1)*crhs21 - DN(1,1)*crhs38 + crhs12*crhs48 + crhs15*(crhs50*normal[1] + normal[0]*(C(0,0)*DN(1,0) + C(0,2)*DN(1,1))) + crhs17*(crhs50*normal[0] + normal[1]*(C(0,1)*DN(1,0) + crhs51)) - crhs25*crhs54 + crhs36*crhs56 - crhs39*crhs54 - crhs47*normal[0] - crhs55*normal[0];
            rhs[4]=-DN(1,0)*crhs41 - DN(1,0)*crhs44 + crhs15*(crhs57*normal[1] + normal[0]*(C(0,1)*DN(1,1) + crhs49)) + crhs17*(crhs57*normal[0] + normal[1]*(C(1,1)*DN(1,1) + C(1,2)*DN(1,0))) + crhs20*crhs48 + crhs37*crhs56 - crhs43*crhs58 - crhs45*crhs58 - crhs47*normal[1] - crhs55*normal[1];
            rhs[5]=N[1]*crhs46;
            rhs[6]=-DN(2,1)*crhs21 - DN(2,1)*crhs38 + crhs12*crhs60 + crhs15*(crhs62*normal[1] + normal[0]*(C(0,0)*DN(2,0) + C(0,2)*DN(2,1))) + crhs17*(crhs62*normal[0] + normal[1]*(C(0,1)*DN(2,0) + crhs63)) - crhs25*crhs66 + crhs36*crhs68 - crhs39*crhs66 - crhs59*normal[0] - crhs67*normal[0];
            rhs[7]=-DN(2,0)*crhs41 - DN(2,0)*crhs44 + crhs15*(crhs69*normal[1] + normal[0]*(C(0,1)*DN(2,1) + crhs61)) + crhs17*(crhs69*normal[0] + normal[1]*(C(1,1)*DN(2,1) + C(1,2)*DN(2,0))) + crhs20*crhs60 + crhs37*crhs68 - crhs43*crhs70 - crhs45*crhs70 - crhs59*normal[1] - crhs67*normal[1];
            rhs[8]=N[2]*crhs46;


    return rhs;
}

template <class TBaseElement>
typename EmbeddedFluidElementDiscontinuous<TBaseElement>::VectorType EmbeddedFluidElementDiscontinuous<TBaseElement>::NitscheTermsRHS3D(
    EmbeddedDiscontinuousElementData& rData)
{
    const double rho = rData.Density;
    const double mu = rData.EffectiveViscosity;

    const double h = rData.ElementSize;
    const double c = rData.SoundVelocity;

    const double dt = rData.DeltaTime;

    const auto& v = rData.Velocity;
    const auto& vmesh = rData.MeshVelocity;
    const auto& vconv = v - vmesh;
    const auto& p = rData.Pressure;
    const auto& stress = rData.ShearStress;

    // If there is embedded velocity, substract it to the previous iteration solution
    array_1d<double,3> g = ZeroVector(3);
    if (this->Has(EMBEDDED_VELOCITY)){
        g = this->GetValue(EMBEDDED_VELOCITY);
    }

    // Get constitutive matrix
    const Matrix& C = rData.C;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN_DX;
    const auto& normal = rData.Normal;

    // Get Nitsche imposition values
    const double eps = rData.SlipLength;;
    const double gamma = rData.PenaltyCoefficient;
    const double adjoint = 1.0;

    VectorType rhs = ZeroVector(16);

    const double crhs0 =             1.0/gamma;
const double crhs1 =             1.0/h;
const double crhs2 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0) - g[0];
const double crhs3 =             crhs2*normal[0];
const double crhs4 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1) - g[1];
const double crhs5 =             crhs4*normal[1];
const double crhs6 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2) - g[2];
const double crhs7 =             crhs6*normal[2];
const double crhs8 =             crhs3 + crhs5 + crhs7;
const double crhs9 =             N[0]*crhs0*crhs1*crhs8*mu;
const double crhs10 =             1.0/(eps + gamma*h);
const double crhs11 =             N[0]*crhs10*mu;
const double crhs12 =             pow(normal[0], 2);
const double crhs13 =             -crhs12 + 1;
const double crhs14 =             -crhs13*crhs2 + crhs5*normal[0] + crhs7*normal[0];
const double crhs15 =             pow(normal[1], 2);
const double crhs16 =             -crhs15 + 1;
const double crhs17 =             -crhs16*crhs4 + crhs3*normal[1] + crhs7*normal[1];
const double crhs18 =             1.0*DN(0,1)*adjoint*crhs10*crhs17*gamma*h*mu;
const double crhs19 =             pow(normal[2], 2);
const double crhs20 =             -crhs19 + 1;
const double crhs21 =             -crhs20*crhs6 + crhs3*normal[2] + crhs5*normal[2];
const double crhs22 =             1.0*DN(0,2)*adjoint*crhs10*crhs21*gamma*h*mu;
const double crhs23 =             DN(0,0)*normal[0];
const double crhs24 =             DN(0,1)*normal[1];
const double crhs25 =             1.0*crhs24;
const double crhs26 =             DN(0,2)*normal[2];
const double crhs27 =             1.0*crhs26;
const double crhs28 =             2.0*crhs23 + crhs25 + crhs27;
const double crhs29 =             adjoint*crhs10*crhs14*gamma*h*mu;
const double crhs30 =             C(0,3)*DN(0,0);
const double crhs31 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + crhs30;
const double crhs32 =             C(0,5)*DN(0,0);
const double crhs33 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + crhs32;
const double crhs34 =             adjoint*crhs8*normal[0];
const double crhs35 =             C(1,3)*DN(0,1);
const double crhs36 =             C(3,4)*DN(0,1);
const double crhs37 =             C(4,5)*DN(0,2);
const double crhs38 =             C(0,4)*DN(0,0) + crhs36 + crhs37;
const double crhs39 =             adjoint*crhs8*normal[1];
const double crhs40 =             C(2,5)*DN(0,2);
const double crhs41 =             adjoint*crhs8*normal[2];
const double crhs42 =             h*rho*sqrt(pow(N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0), 2) + pow(N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1), 2) + pow(N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2), 2)) + mu + pow(h, 2)*rho/dt;
const double crhs43 =             N[0]*crhs0*crhs1*crhs42*crhs8;
const double crhs44 =             N[0]*crhs10*eps;
const double crhs45 =             normal[0]*normal[1];
const double crhs46 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs47 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs48 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs49 =             DN(0,0)*v(0,1) + DN(0,1)*v(0,0) + DN(1,0)*v(1,1) + DN(1,1)*v(1,0) + DN(2,0)*v(2,1) + DN(2,1)*v(2,0) + DN(3,0)*v(3,1) + DN(3,1)*v(3,0);
const double crhs50 =             DN(0,1)*v(0,2) + DN(0,2)*v(0,1) + DN(1,1)*v(1,2) + DN(1,2)*v(1,1) + DN(2,1)*v(2,2) + DN(2,2)*v(2,1) + DN(3,1)*v(3,2) + DN(3,2)*v(3,1);
const double crhs51 =             DN(0,0)*v(0,2) + DN(0,2)*v(0,0) + DN(1,0)*v(1,2) + DN(1,2)*v(1,0) + DN(2,0)*v(2,2) + DN(2,2)*v(2,0) + DN(3,0)*v(3,2) + DN(3,2)*v(3,0);
const double crhs52 =             C(0,3)*crhs46 + C(1,3)*crhs47 + C(2,3)*crhs48 + C(3,3)*crhs49 + C(3,4)*crhs50 + C(3,5)*crhs51;
const double crhs53 =             C(0,4)*crhs46 + C(1,4)*crhs47 + C(2,4)*crhs48 + C(3,4)*crhs49 + C(4,4)*crhs50 + C(4,5)*crhs51;
const double crhs54 =             crhs52*normal[0] + crhs53*normal[2] + normal[1]*(C(0,1)*crhs46 + C(1,1)*crhs47 + C(1,2)*crhs48 + C(1,3)*crhs49 + C(1,4)*crhs50 + C(1,5)*crhs51);
const double crhs55 =             normal[0]*normal[2];
const double crhs56 =             C(0,5)*crhs46 + C(1,5)*crhs47 + C(2,5)*crhs48 + C(3,5)*crhs49 + C(4,5)*crhs50 + C(5,5)*crhs51;
const double crhs57 =             crhs53*normal[1] + crhs56*normal[0] + normal[2]*(C(0,2)*crhs46 + C(1,2)*crhs47 + C(2,2)*crhs48 + C(2,3)*crhs49 + C(2,4)*crhs50 + C(2,5)*crhs51);
const double crhs58 =             crhs52*normal[1] + crhs56*normal[2] + normal[0]*(C(0,0)*crhs46 + C(0,1)*crhs47 + C(0,2)*crhs48 + C(0,3)*crhs49 + C(0,4)*crhs50 + C(0,5)*crhs51);
const double crhs59 =             -crhs13*crhs58 + crhs45*crhs54 + crhs55*crhs57;
const double crhs60 =             normal[1]*normal[2];
const double crhs61 =             -crhs16*crhs54 + crhs45*crhs58 + crhs57*crhs60;
const double crhs62 =             1.0*DN(0,1)*adjoint*crhs10*crhs61*eps*gamma*h;
const double crhs63 =             -crhs20*crhs57 + crhs54*crhs60 + crhs55*crhs58;
const double crhs64 =             1.0*DN(0,2)*adjoint*crhs10*crhs63*eps*gamma*h;
const double crhs65 =             adjoint*crhs10*crhs59*eps*gamma*h;
const double crhs66 =             1.0*DN(0,0)*adjoint*crhs10*crhs14*gamma*h*mu;
const double crhs67 =             1.0*crhs23;
const double crhs68 =             2.0*crhs24 + crhs27 + crhs67;
const double crhs69 =             adjoint*crhs10*crhs17*gamma*h*mu;
const double crhs70 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + crhs35;
const double crhs71 =             C(3,5)*DN(0,0);
const double crhs72 =             C(1,5)*DN(0,1) + crhs37 + crhs71;
const double crhs73 =             C(1,4)*DN(0,1);
const double crhs74 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + crhs73;
const double crhs75 =             C(2,4)*DN(0,2);
const double crhs76 =             1.0*DN(0,0)*adjoint*crhs10*crhs59*eps*gamma*h;
const double crhs77 =             adjoint*crhs10*crhs61*eps*gamma*h;
const double crhs78 =             crhs25 + 2.0*crhs26 + crhs67;
const double crhs79 =             adjoint*crhs10*crhs21*gamma*h*mu;
const double crhs80 =             C(2,3)*DN(0,2) + crhs36 + crhs71;
const double crhs81 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + crhs40;
const double crhs82 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + crhs75;
const double crhs83 =             adjoint*crhs10*crhs63*eps*gamma*h;
const double crhs84 =             crhs8*(crhs12 + crhs15 + crhs19);
const double crhs85 =             N[1]*crhs0*crhs1*crhs8*mu;
const double crhs86 =             N[1]*crhs10*mu;
const double crhs87 =             1.0*DN(1,1)*adjoint*crhs10*crhs17*gamma*h*mu;
const double crhs88 =             1.0*DN(1,2)*adjoint*crhs10*crhs21*gamma*h*mu;
const double crhs89 =             DN(1,0)*normal[0];
const double crhs90 =             DN(1,1)*normal[1];
const double crhs91 =             1.0*crhs90;
const double crhs92 =             DN(1,2)*normal[2];
const double crhs93 =             1.0*crhs92;
const double crhs94 =             2.0*crhs89 + crhs91 + crhs93;
const double crhs95 =             C(0,3)*DN(1,0);
const double crhs96 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + crhs95;
const double crhs97 =             C(0,5)*DN(1,0);
const double crhs98 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + crhs97;
const double crhs99 =             C(1,3)*DN(1,1);
const double crhs100 =             C(3,4)*DN(1,1);
const double crhs101 =             C(4,5)*DN(1,2);
const double crhs102 =             C(0,4)*DN(1,0) + crhs100 + crhs101;
const double crhs103 =             C(2,5)*DN(1,2);
const double crhs104 =             N[1]*crhs0*crhs1*crhs42*crhs8;
const double crhs105 =             N[1]*crhs10*eps;
const double crhs106 =             1.0*DN(1,1)*adjoint*crhs10*crhs61*eps*gamma*h;
const double crhs107 =             1.0*DN(1,2)*adjoint*crhs10*crhs63*eps*gamma*h;
const double crhs108 =             1.0*DN(1,0)*adjoint*crhs10*crhs14*gamma*h*mu;
const double crhs109 =             1.0*crhs89;
const double crhs110 =             crhs109 + 2.0*crhs90 + crhs93;
const double crhs111 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + crhs99;
const double crhs112 =             C(3,5)*DN(1,0);
const double crhs113 =             C(1,5)*DN(1,1) + crhs101 + crhs112;
const double crhs114 =             C(1,4)*DN(1,1);
const double crhs115 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + crhs114;
const double crhs116 =             C(2,4)*DN(1,2);
const double crhs117 =             1.0*DN(1,0)*adjoint*crhs10*crhs59*eps*gamma*h;
const double crhs118 =             crhs109 + crhs91 + 2.0*crhs92;
const double crhs119 =             C(2,3)*DN(1,2) + crhs100 + crhs112;
const double crhs120 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + crhs103;
const double crhs121 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + crhs116;
const double crhs122 =             N[2]*crhs0*crhs1*crhs8*mu;
const double crhs123 =             N[2]*crhs10*mu;
const double crhs124 =             1.0*DN(2,1)*adjoint*crhs10*crhs17*gamma*h*mu;
const double crhs125 =             1.0*DN(2,2)*adjoint*crhs10*crhs21*gamma*h*mu;
const double crhs126 =             DN(2,0)*normal[0];
const double crhs127 =             DN(2,1)*normal[1];
const double crhs128 =             1.0*crhs127;
const double crhs129 =             DN(2,2)*normal[2];
const double crhs130 =             1.0*crhs129;
const double crhs131 =             2.0*crhs126 + crhs128 + crhs130;
const double crhs132 =             C(0,3)*DN(2,0);
const double crhs133 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + crhs132;
const double crhs134 =             C(0,5)*DN(2,0);
const double crhs135 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + crhs134;
const double crhs136 =             C(1,3)*DN(2,1);
const double crhs137 =             C(3,4)*DN(2,1);
const double crhs138 =             C(4,5)*DN(2,2);
const double crhs139 =             C(0,4)*DN(2,0) + crhs137 + crhs138;
const double crhs140 =             C(2,5)*DN(2,2);
const double crhs141 =             N[2]*crhs0*crhs1*crhs42*crhs8;
const double crhs142 =             N[2]*crhs10*eps;
const double crhs143 =             1.0*DN(2,1)*adjoint*crhs10*crhs61*eps*gamma*h;
const double crhs144 =             1.0*DN(2,2)*adjoint*crhs10*crhs63*eps*gamma*h;
const double crhs145 =             1.0*DN(2,0)*adjoint*crhs10*crhs14*gamma*h*mu;
const double crhs146 =             1.0*crhs126;
const double crhs147 =             2.0*crhs127 + crhs130 + crhs146;
const double crhs148 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + crhs136;
const double crhs149 =             C(3,5)*DN(2,0);
const double crhs150 =             C(1,5)*DN(2,1) + crhs138 + crhs149;
const double crhs151 =             C(1,4)*DN(2,1);
const double crhs152 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + crhs151;
const double crhs153 =             C(2,4)*DN(2,2);
const double crhs154 =             1.0*DN(2,0)*adjoint*crhs10*crhs59*eps*gamma*h;
const double crhs155 =             crhs128 + 2.0*crhs129 + crhs146;
const double crhs156 =             C(2,3)*DN(2,2) + crhs137 + crhs149;
const double crhs157 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + crhs140;
const double crhs158 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + crhs153;
const double crhs159 =             N[3]*crhs0*crhs1*crhs8*mu;
const double crhs160 =             N[3]*crhs10*mu;
const double crhs161 =             1.0*DN(3,1)*adjoint*crhs10*crhs17*gamma*h*mu;
const double crhs162 =             1.0*DN(3,2)*adjoint*crhs10*crhs21*gamma*h*mu;
const double crhs163 =             DN(3,0)*normal[0];
const double crhs164 =             DN(3,1)*normal[1];
const double crhs165 =             1.0*crhs164;
const double crhs166 =             DN(3,2)*normal[2];
const double crhs167 =             1.0*crhs166;
const double crhs168 =             2.0*crhs163 + crhs165 + crhs167;
const double crhs169 =             C(0,3)*DN(3,0);
const double crhs170 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + crhs169;
const double crhs171 =             C(0,5)*DN(3,0);
const double crhs172 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + crhs171;
const double crhs173 =             C(1,3)*DN(3,1);
const double crhs174 =             C(3,4)*DN(3,1);
const double crhs175 =             C(4,5)*DN(3,2);
const double crhs176 =             C(0,4)*DN(3,0) + crhs174 + crhs175;
const double crhs177 =             C(2,5)*DN(3,2);
const double crhs178 =             N[3]*crhs0*crhs1*crhs42*crhs8;
const double crhs179 =             N[3]*crhs10*eps;
const double crhs180 =             1.0*DN(3,1)*adjoint*crhs10*crhs61*eps*gamma*h;
const double crhs181 =             1.0*DN(3,2)*adjoint*crhs10*crhs63*eps*gamma*h;
const double crhs182 =             1.0*DN(3,0)*adjoint*crhs10*crhs14*gamma*h*mu;
const double crhs183 =             1.0*crhs163;
const double crhs184 =             2.0*crhs164 + crhs167 + crhs183;
const double crhs185 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + crhs173;
const double crhs186 =             C(3,5)*DN(3,0);
const double crhs187 =             C(1,5)*DN(3,1) + crhs175 + crhs186;
const double crhs188 =             C(1,4)*DN(3,1);
const double crhs189 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + crhs188;
const double crhs190 =             C(2,4)*DN(3,2);
const double crhs191 =             1.0*DN(3,0)*adjoint*crhs10*crhs59*eps*gamma*h;
const double crhs192 =             crhs165 + 2.0*crhs166 + crhs183;
const double crhs193 =             C(2,3)*DN(3,2) + crhs174 + crhs186;
const double crhs194 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + crhs177;
const double crhs195 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + crhs190;
            rhs[0]=crhs11*crhs14 - crhs18*normal[0] - crhs22*normal[0] - crhs28*crhs29 - crhs28*crhs65 + crhs34*(crhs31*normal[1] + crhs33*normal[2] + normal[0]*(C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2))) + crhs39*(crhs31*normal[0] + crhs38*normal[2] + normal[1]*(C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + crhs35)) + crhs41*(crhs33*normal[0] + crhs38*normal[1] + normal[2]*(C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + crhs40)) - crhs43*normal[0] + crhs44*crhs59 - crhs62*normal[0] - crhs64*normal[0] - crhs9*normal[0];
            rhs[1]=crhs11*crhs17 - crhs22*normal[1] + crhs34*(crhs70*normal[1] + crhs72*normal[2] + normal[0]*(C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + crhs30)) + crhs39*(crhs70*normal[0] + crhs74*normal[2] + normal[1]*(C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2))) + crhs41*(crhs72*normal[0] + crhs74*normal[1] + normal[2]*(C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + crhs75)) - crhs43*normal[1] + crhs44*crhs61 - crhs64*normal[1] - crhs66*normal[1] - crhs68*crhs69 - crhs68*crhs77 - crhs76*normal[1] - crhs9*normal[1];
            rhs[2]=crhs11*crhs21 - crhs18*normal[2] + crhs34*(crhs80*normal[1] + crhs81*normal[2] + normal[0]*(C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + crhs32)) + crhs39*(crhs80*normal[0] + crhs82*normal[2] + normal[1]*(C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + crhs73)) + crhs41*(crhs81*normal[0] + crhs82*normal[1] + normal[2]*(C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0))) - crhs43*normal[2] + crhs44*crhs63 - crhs62*normal[2] - crhs66*normal[2] - crhs76*normal[2] - crhs78*crhs79 - crhs78*crhs83 - crhs9*normal[2];
            rhs[3]=N[0]*crhs84;
            rhs[4]=-crhs104*normal[0] + crhs105*crhs59 - crhs106*normal[0] - crhs107*normal[0] + crhs14*crhs86 - crhs29*crhs94 + crhs34*(crhs96*normal[1] + crhs98*normal[2] + normal[0]*(C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2))) + crhs39*(crhs102*normal[2] + crhs96*normal[0] + normal[1]*(C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + crhs99)) + crhs41*(crhs102*normal[1] + crhs98*normal[0] + normal[2]*(C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + crhs103)) - crhs65*crhs94 - crhs85*normal[0] - crhs87*normal[0] - crhs88*normal[0];
            rhs[5]=-crhs104*normal[1] + crhs105*crhs61 - crhs107*normal[1] - crhs108*normal[1] - crhs110*crhs69 - crhs110*crhs77 - crhs117*normal[1] + crhs17*crhs86 + crhs34*(crhs111*normal[1] + crhs113*normal[2] + normal[0]*(C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + crhs95)) + crhs39*(crhs111*normal[0] + crhs115*normal[2] + normal[1]*(C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2))) + crhs41*(crhs113*normal[0] + crhs115*normal[1] + normal[2]*(C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + crhs116)) - crhs85*normal[1] - crhs88*normal[1];
            rhs[6]=-crhs104*normal[2] + crhs105*crhs63 - crhs106*normal[2] - crhs108*normal[2] - crhs117*normal[2] - crhs118*crhs79 - crhs118*crhs83 + crhs21*crhs86 + crhs34*(crhs119*normal[1] + crhs120*normal[2] + normal[0]*(C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + crhs97)) + crhs39*(crhs119*normal[0] + crhs121*normal[2] + normal[1]*(C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + crhs114)) + crhs41*(crhs120*normal[0] + crhs121*normal[1] + normal[2]*(C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0))) - crhs85*normal[2] - crhs87*normal[2];
            rhs[7]=N[1]*crhs84;
            rhs[8]=-crhs122*normal[0] + crhs123*crhs14 - crhs124*normal[0] - crhs125*normal[0] - crhs131*crhs29 - crhs131*crhs65 - crhs141*normal[0] + crhs142*crhs59 - crhs143*normal[0] - crhs144*normal[0] + crhs34*(crhs133*normal[1] + crhs135*normal[2] + normal[0]*(C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2))) + crhs39*(crhs133*normal[0] + crhs139*normal[2] + normal[1]*(C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + crhs136)) + crhs41*(crhs135*normal[0] + crhs139*normal[1] + normal[2]*(C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + crhs140));
            rhs[9]=-crhs122*normal[1] + crhs123*crhs17 - crhs125*normal[1] - crhs141*normal[1] + crhs142*crhs61 - crhs144*normal[1] - crhs145*normal[1] - crhs147*crhs69 - crhs147*crhs77 - crhs154*normal[1] + crhs34*(crhs148*normal[1] + crhs150*normal[2] + normal[0]*(C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + crhs132)) + crhs39*(crhs148*normal[0] + crhs152*normal[2] + normal[1]*(C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2))) + crhs41*(crhs150*normal[0] + crhs152*normal[1] + normal[2]*(C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + crhs153));
            rhs[10]=-crhs122*normal[2] + crhs123*crhs21 - crhs124*normal[2] - crhs141*normal[2] + crhs142*crhs63 - crhs143*normal[2] - crhs145*normal[2] - crhs154*normal[2] - crhs155*crhs79 - crhs155*crhs83 + crhs34*(crhs156*normal[1] + crhs157*normal[2] + normal[0]*(C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + crhs134)) + crhs39*(crhs156*normal[0] + crhs158*normal[2] + normal[1]*(C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + crhs151)) + crhs41*(crhs157*normal[0] + crhs158*normal[1] + normal[2]*(C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0)));
            rhs[11]=N[2]*crhs84;
            rhs[12]=crhs14*crhs160 - crhs159*normal[0] - crhs161*normal[0] - crhs162*normal[0] - crhs168*crhs29 - crhs168*crhs65 - crhs178*normal[0] + crhs179*crhs59 - crhs180*normal[0] - crhs181*normal[0] + crhs34*(crhs170*normal[1] + crhs172*normal[2] + normal[0]*(C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2))) + crhs39*(crhs170*normal[0] + crhs176*normal[2] + normal[1]*(C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + crhs173)) + crhs41*(crhs172*normal[0] + crhs176*normal[1] + normal[2]*(C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + crhs177));
            rhs[13]=-crhs159*normal[1] + crhs160*crhs17 - crhs162*normal[1] - crhs178*normal[1] + crhs179*crhs61 - crhs181*normal[1] - crhs182*normal[1] - crhs184*crhs69 - crhs184*crhs77 - crhs191*normal[1] + crhs34*(crhs185*normal[1] + crhs187*normal[2] + normal[0]*(C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + crhs169)) + crhs39*(crhs185*normal[0] + crhs189*normal[2] + normal[1]*(C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2))) + crhs41*(crhs187*normal[0] + crhs189*normal[1] + normal[2]*(C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + crhs190));
            rhs[14]=-crhs159*normal[2] + crhs160*crhs21 - crhs161*normal[2] - crhs178*normal[2] + crhs179*crhs63 - crhs180*normal[2] - crhs182*normal[2] - crhs191*normal[2] - crhs192*crhs79 - crhs192*crhs83 + crhs34*(crhs193*normal[1] + crhs194*normal[2] + normal[0]*(C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + crhs171)) + crhs39*(crhs193*normal[0] + crhs195*normal[2] + normal[1]*(C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + crhs188)) + crhs41*(crhs194*normal[0] + crhs195*normal[1] + normal[2]*(C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0)));
            rhs[15]=N[3]*crhs84;


    return rhs;
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::AddNormalSymmetricCounterpartContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedDiscontinuousElementData& rData) const
{
    // // Obtain the previous iteration velocity solution
    // array_1d<double,LocalSize> values;
    // this->GetCurrentValuesVector(rData,values);

    // // If there is embedded velocity, substract it to the previous iteration solution
    // if (this->Has(EMBEDDED_VELOCITY)) {
    //     const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
    //     array_1d<double, LocalSize> embedded_vel_exp = ZeroVector(LocalSize);

    //     for (unsigned int i = 0; i < NumNodes; ++i) {
    //         for (unsigned int comp = 0; comp < Dim; ++comp) {
    //             embedded_vel_exp(i*BlockSize + comp) = embedded_vel(comp);
    //         }
    //     }

    //     noalias(values) -= embedded_vel_exp;
    // }

    // // Set if the shear stress term is adjoint consistent (1.0) or not (-1.0)
    // const double adjoint_consistency = -1.0;

    // // Set an auxiliar array to compute the LHS contribution
    // BoundedMatrix<double, LocalSize, LocalSize> aux_LHS = ZeroMatrix(LocalSize, LocalSize);

    // // Compute positive side LHS contribution
    // const unsigned int number_of_positive_interface_integration_points = rData.PositiveInterfaceWeights.size();
    // for (unsigned int g = 0; g < number_of_positive_interface_integration_points; ++g){
    //     // Get the Gauss pt. data
    //     const double weight = rData.PositiveInterfaceWeights[g];
    //     const auto aux_N = row(rData.PositiveInterfaceN, g);
    //     const BoundedMatrix<double, NumNodes, Dim> &aux_DN_DX = rData.PositiveInterfaceDNDX[g];
    //     const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

    //     // Fill the pressure to Voigt notation operator normal projected matrix
    //     BoundedMatrix<double, LocalSize, Dim> trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(LocalSize, Dim);
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         for (unsigned int comp = 0; comp < Dim; ++comp){
    //             trans_pres_to_voigt_matrix_normal_op(i*BlockSize + Dim, comp) = aux_N(i)*aux_unit_normal(comp);
    //         }
    //     }

    //     // Set the shape functions auxiliar matrix
    //     BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         for (unsigned int comp = 0; comp < Dim; ++comp){
    //             N_mat(comp, i*BlockSize + comp) = aux_N(i);
    //         }
    //     }

    //     // Set the current Gauss pt. strain matrix
    //     BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    //     FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

    //     // Set the normal projection matrix (n x n)
    //     BoundedMatrix<double, Dim, Dim> normal_proj_matrix;
    //     FluidElementUtilities<NumNodes>::SetNormalProjectionMatrix(aux_unit_normal, normal_proj_matrix);

    //     // Get the normal projection matrix in Voigt notation
    //     BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
    //     FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

    //     // Compute some Gauss pt. auxiliar matrices
    //     const BoundedMatrix<double, LocalSize, StrainSize> aux_matrix_BC = prod(trans(B_matrix), trans(rData.C));
    //     const BoundedMatrix<double, StrainSize, Dim> aux_matrix_APnorm = prod(trans(voigt_normal_proj_matrix), normal_proj_matrix);
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BCAPnorm = prod(aux_matrix_BC, aux_matrix_APnorm);

    //     // Contribution coming fron the shear stress operator
    //     noalias(aux_LHS) -= adjoint_consistency*weight*prod(aux_matrix_BCAPnorm, N_mat);

    //     // Contribution coming from the pressure terms
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_VPnorm = prod(trans_pres_to_voigt_matrix_normal_op, normal_proj_matrix);
    //     noalias(aux_LHS) -= weight*prod(aux_matrix_VPnorm, N_mat);
    // }

    // // Compute negative side LHS contribution
    // const unsigned int number_of_negative_interface_integration_points = rData.NegativeInterfaceWeights.size();
    // for (unsigned int g = 0; g < number_of_negative_interface_integration_points; ++g){
    //     // Get the Gauss pt. data
    //     const double weight = rData.NegativeInterfaceWeights[g];
    //     const auto aux_N = row(rData.NegativeInterfaceN, g);
    //     const BoundedMatrix<double, NumNodes, Dim> &aux_DN_DX = rData.NegativeInterfaceDNDX[g];
    //     const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

    //     // Fill the pressure to Voigt notation operator normal projected matrix
    //     BoundedMatrix<double, LocalSize, Dim> trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(LocalSize, Dim);
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         for (unsigned int comp = 0; comp < Dim; ++comp){
    //             trans_pres_to_voigt_matrix_normal_op(i*BlockSize + Dim, comp) = aux_N(i)*aux_unit_normal(comp);
    //         }
    //     }

    //     // Set the shape functions auxiliar matrix
    //     BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         for (unsigned int comp = 0; comp < Dim; ++comp){
    //             N_mat(comp, i*BlockSize + comp) = aux_N(i);
    //         }
    //     }

    //     // Set the current Gauss pt. strain matrix
    //     BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    //     FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

    //     // Set the normal projection matrix (n x n)
    //     BoundedMatrix<double, Dim, Dim> normal_proj_matrix;
    //     FluidElementUtilities<NumNodes>::SetNormalProjectionMatrix(aux_unit_normal, normal_proj_matrix);

    //     // Get the normal projection matrix in Voigt notation
    //     BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
    //     FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

    //     // Compute some Gauss pt. auxiliar matrices
    //     const BoundedMatrix<double, LocalSize, StrainSize> aux_matrix_BC = prod(trans(B_matrix), trans(rData.C));
    //     const BoundedMatrix<double, StrainSize, Dim> aux_matrix_APnorm = prod(trans(voigt_normal_proj_matrix), normal_proj_matrix);
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BCAPnorm = prod(aux_matrix_BC, aux_matrix_APnorm);

    //     // Contribution coming fron the shear stress operator
    //     noalias(aux_LHS) -= adjoint_consistency*weight*prod(aux_matrix_BCAPnorm, N_mat);

    //     // Contribution coming from the pressure terms
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_VPnorm = prod(trans_pres_to_voigt_matrix_normal_op, normal_proj_matrix);
    //     noalias(aux_LHS) -= weight*prod(aux_matrix_VPnorm, N_mat);
    // }

    // // LHS outside Nitsche contribution assembly
    // noalias(rLHS) += aux_LHS;

    // // RHS outside Nitsche contribution assembly
    // // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
    // noalias(rRHS) -= prod(aux_LHS, values);
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::AddTangentialPenaltyContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedDiscontinuousElementData& rData) const
{
    // // Obtain the previous iteration velocity solution
    // array_1d<double,LocalSize> values;
    // this->GetCurrentValuesVector(rData, values);

    // // Compute the Nitsche tangential imposition penalty coefficients
    // std::pair<const double, const double> pen_coefs = this->ComputeTangentialPenaltyCoefficients(rData);

    // // Declare auxiliar arrays
    // BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_1 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution coming from the tangential component of the Cauchy stress vector
    // BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_2 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution generated by the viscous shear force generated by the velocity

    // // Compute positive side LHS contribution 
    // const unsigned int number_of_positive_interface_integration_points = rData.PositiveInterfaceWeights.size();
    // for (unsigned int g = 0; g < number_of_positive_interface_integration_points; ++g){
    //     // Get the Gauss pt. data
    //     const double weight = rData.PositiveInterfaceWeights[g];
    //     const auto aux_N = row(rData.PositiveInterfaceN, g);
    //     const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.PositiveInterfaceDNDX[g];
    //     const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

    //     // Set the shape functions auxiliar matrices
    //     BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         for (unsigned int comp = 0; comp < Dim; ++comp){
    //             N_mat(comp, i*BlockSize + comp) = aux_N(i);
    //         }
    //     }
    //     BoundedMatrix<double, LocalSize, Dim> N_mat_trans = trans(N_mat);

    //     // Set the tangential projection matrix (I - n x n)
    //     BoundedMatrix<double, Dim, Dim> tang_proj_matrix;
    //     FluidElementUtilities<NumNodes>::SetTangentialProjectionMatrix(aux_unit_normal, tang_proj_matrix);

    //     // Set the current Gauss pt. strain matrix
    //     BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    //     FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

    //     // Get the normal projection matrix in Voigt notation
    //     BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
    //     FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

    //     // Compute some Gauss pt. auxiliar matrices
    //     const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_CB = prod(rData.C, B_matrix);
    //     const BoundedMatrix<double, StrainSize, Dim> aux_matrix_PtangA = prod(tang_proj_matrix, voigt_normal_proj_matrix);
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_PtangACB = prod(aux_matrix_PtangA, aux_matrix_CB);

    //     // Contribution coming from the traction vector tangencial component
    //     noalias(aux_LHS_1) += pen_coefs.first*weight*prod(N_mat_trans, aux_matrix_PtangACB);

    //     // Contribution coming from the shear force generated by the velocity jump
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_N_trans_tang = prod(N_mat_trans, tang_proj_matrix);
    //     noalias(aux_LHS_2) += pen_coefs.second*weight*prod(aux_matrix_N_trans_tang, N_mat);
    // }

    // // Compute negative side LHS contribution 
    // const unsigned int number_of_negative_interface_integration_points = rData.NegativeInterfaceWeights.size();
    // for (unsigned int g = 0; g < number_of_negative_interface_integration_points; ++g){
    //     // Get the Gauss pt. data
    //     const double weight = rData.NegativeInterfaceWeights[g];
    //     const auto aux_N = row(rData.NegativeInterfaceN, g);
    //     const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.NegativeInterfaceDNDX[g];
    //     const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

    //     // Set the shape functions auxiliar matrices
    //     BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         for (unsigned int comp = 0; comp < Dim; ++comp){
    //             N_mat(comp, i*BlockSize + comp) = aux_N(i);
    //         }
    //     }
    //     BoundedMatrix<double, LocalSize, Dim> N_mat_trans = trans(N_mat);

    //     // Set the tangential projection matrix (I - n x n)
    //     BoundedMatrix<double, Dim, Dim> tang_proj_matrix;
    //     FluidElementUtilities<NumNodes>::SetTangentialProjectionMatrix(aux_unit_normal, tang_proj_matrix);

    //     // Set the current Gauss pt. strain matrix
    //     BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    //     FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

    //     // Get the normal projection matrix in Voigt notation
    //     BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
    //     FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

    //     // Compute some Gauss pt. auxiliar matrices
    //     const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_CB = prod(rData.C, B_matrix);
    //     const BoundedMatrix<double, StrainSize, Dim> aux_matrix_PtangA = prod(tang_proj_matrix, voigt_normal_proj_matrix);
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_PtangACB = prod(aux_matrix_PtangA, aux_matrix_CB);

    //     // Contribution coming from the traction vector tangencial component
    //     noalias(aux_LHS_1) += pen_coefs.first*weight*prod(N_mat_trans, aux_matrix_PtangACB);

    //     // Contribution coming from the shear force generated by the velocity jump
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_N_trans_tang = prod(N_mat_trans, tang_proj_matrix);
    //     noalias(aux_LHS_2) += pen_coefs.second*weight*prod(aux_matrix_N_trans_tang, N_mat);
    // }

    // // LHS outside Nitsche contribution assembly
    // noalias(rLHS) += aux_LHS_1;
    // noalias(rLHS) += aux_LHS_2;

    // // RHS outside Nitsche contribution assembly
    // // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
    // noalias(rRHS) -= prod(aux_LHS_1, values);
    // noalias(rRHS) -= prod(aux_LHS_2, values);

    // // If level set velocity is not 0, add its contribution to the RHS
    // if (this->Has(EMBEDDED_VELOCITY)) {
    //     const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
    //     array_1d<double, LocalSize> embedded_vel_exp = ZeroVector(LocalSize);

    //     for (unsigned int i = 0; i < NumNodes; ++i) {
    //         for (unsigned int comp = 0; comp < Dim; ++comp) {
    //             embedded_vel_exp(i*BlockSize + comp) = embedded_vel(comp);
    //         }
    //     }
    //     noalias(rRHS) += prod(aux_LHS_2, embedded_vel_exp);
    // }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::AddTangentialSymmetricCounterpartContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedDiscontinuousElementData& rData) const
{
    // // Obtain the previous iteration velocity solution
    // array_1d<double,LocalSize> values;
    // this->GetCurrentValuesVector(rData, values);

    // // Set if the shear stress term is adjoint consistent (1.0) or not (-1.0)
    // const double adjoint_consistency = -1.0;

    // // Compute the coefficients
    // std::pair<const double, const double> nitsche_coefs = this->ComputeTangentialNitscheCoefficients(rData);

    // // Declare auxiliar arrays
    // BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_1 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution coming from the tangential component of the Cauchy stress vector
    // BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_2 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution generated by the viscous shear force generated by the velocity

    // // Compute positive side LHS contribution
    // const unsigned int number_of_positive_interface_integration_points = rData.PositiveInterfaceWeights.size();
    // for (unsigned int g = 0; g < number_of_positive_interface_integration_points; ++g){
    //     // Get the Gauss pt. data
    //     const double weight = rData.PositiveInterfaceWeights[g];
    //     const auto aux_N = row(rData.PositiveInterfaceN, g);
    //     const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.PositiveInterfaceDNDX[g];
    //     const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

    //     // Set the shape functions auxiliar matrices
    //     BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         for (unsigned int comp = 0; comp < Dim; ++comp){
    //             N_mat(comp, i*BlockSize + comp) = aux_N(i);
    //         }
    //     }

    //     // Set the current Gauss pt. strain matrix
    //     BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    //     FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

    //     // Set the tangential projection matrix (I - n x n)
    //     BoundedMatrix<double, Dim, Dim> tang_proj_matrix;
    //     FluidElementUtilities<NumNodes>::SetTangentialProjectionMatrix(aux_unit_normal, tang_proj_matrix);

    //     // Get the normal projection matrix in Voigt notation
    //     BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
    //     FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

    //     // Compute some Gauss pt. auxiliar matrices
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BtransAtrans = prod(trans(B_matrix), trans(voigt_normal_proj_matrix));
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BtransAtransPtan = prod(aux_matrix_BtransAtrans, tang_proj_matrix);
    //     const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_CB = prod(rData.C, B_matrix);
    //     const BoundedMatrix<double, Dim, LocalSize> aux_matrix_ACB = prod(voigt_normal_proj_matrix, aux_matrix_CB);
    //     const BoundedMatrix<double, LocalSize, LocalSize> aux_matrix_BtransAtransPtanACB = prod(aux_matrix_BtransAtransPtan, aux_matrix_ACB);

    //     // Contribution coming from the traction vector tangencial component
    //     noalias(aux_LHS_1) -= adjoint_consistency*nitsche_coefs.first*weight*aux_matrix_BtransAtransPtanACB;

    //     // Contribution coming from the shear force generated by the velocity jump
    //     noalias(aux_LHS_2) -= adjoint_consistency*nitsche_coefs.second*weight*prod(aux_matrix_BtransAtransPtan, N_mat);
    // }

    // // Compute negative side LHS contribution
    // const unsigned int number_of_negative_interface_integration_points = rData.NegativeInterfaceWeights.size();
    // for (unsigned int g = 0; g < number_of_negative_interface_integration_points; ++g){
    //     // Get the Gauss pt. data
    //     const double weight = rData.NegativeInterfaceWeights[g];
    //     const auto aux_N = row(rData.NegativeInterfaceN, g);
    //     const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.NegativeInterfaceDNDX[g];
    //     const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

    //     // Set the shape functions auxiliar matrices
    //     BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
    //     for (unsigned int i = 0; i < NumNodes; ++i){
    //         for (unsigned int comp = 0; comp < Dim; ++comp){
    //             N_mat(comp, i*BlockSize + comp) = aux_N(i);
    //         }
    //     }

    //     // Set the current Gauss pt. strain matrix
    //     BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
    //     FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

    //     // Set the tangential projection matrix (I - n x n)
    //     BoundedMatrix<double, Dim, Dim> tang_proj_matrix;
    //     FluidElementUtilities<NumNodes>::SetTangentialProjectionMatrix(aux_unit_normal, tang_proj_matrix);

    //     // Get the normal projection matrix in Voigt notation
    //     BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
    //     FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

    //     // Compute some Gauss pt. auxiliar matrices
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BtransAtrans = prod(trans(B_matrix), trans(voigt_normal_proj_matrix));
    //     const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BtransAtransPtan = prod(aux_matrix_BtransAtrans, tang_proj_matrix);
    //     const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_CB = prod(rData.C, B_matrix);
    //     const BoundedMatrix<double, Dim, LocalSize> aux_matrix_ACB = prod(voigt_normal_proj_matrix, aux_matrix_CB);
    //     const BoundedMatrix<double, LocalSize, LocalSize> aux_matrix_BtransAtransPtanACB = prod(aux_matrix_BtransAtransPtan, aux_matrix_ACB);

    //     // Contribution coming from the traction vector tangencial component
    //     noalias(aux_LHS_1) -= adjoint_consistency*nitsche_coefs.first*weight*aux_matrix_BtransAtransPtanACB;

    //     // Contribution coming from the shear force generated by the velocity jump
    //     noalias(aux_LHS_2) -= adjoint_consistency*nitsche_coefs.second*weight*prod(aux_matrix_BtransAtransPtan, N_mat);
    // }

    // // LHS outside Nitsche contribution assembly
    // noalias(rLHS) += aux_LHS_1;
    // noalias(rLHS) += aux_LHS_2;

    // // RHS outside Nitsche contribution assembly
    // // If level set velocity is not 0, add its contribution to the RHS
    // if (this->Has(EMBEDDED_VELOCITY)) {
    //     const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
    //     array_1d<double, LocalSize> embedded_vel_exp = ZeroVector(LocalSize);

    //     for (unsigned int i = 0; i < NumNodes; ++i) {
    //         for (unsigned int comp = 0; comp < Dim; ++comp) {
    //             embedded_vel_exp(i*BlockSize + comp) = embedded_vel(comp);
    //         }
    //     }

    //     noalias(rRHS) += prod(aux_LHS_2, embedded_vel_exp);
    // }

    // // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
    // noalias(rRHS) -= prod(aux_LHS_1, values);
    // noalias(rRHS) -= prod(aux_LHS_2, values);

}

template <class TBaseElement>
double EmbeddedFluidElementDiscontinuous<TBaseElement>::ComputeNormalPenaltyCoefficient(const EmbeddedDiscontinuousElementData& rData) const
{
    // Compute the element average velocity norm
    double v_norm = 0.0;
    for (unsigned int comp = 0; comp < Dim; ++comp){
        double aux_vel = 0.0;
        for (unsigned int j = 0; j < NumNodes; ++j){
            aux_vel += rData.Velocity(j,comp);
        }
        aux_vel /= NumNodes;
        v_norm += aux_vel*aux_vel;
    }
    v_norm = std::sqrt(v_norm);

    // Compute the Nitsche coefficient (including the Winter stabilization term)
    const double h = rData.ElementSize;
    const double avg_rho = rData.Density;
    const double eff_mu = rData.EffectiveViscosity;
    const double penalty = rData.PenaltyCoefficient;
    const double cons_coef = (eff_mu + eff_mu + avg_rho*v_norm*h + avg_rho*h*h/rData.DeltaTime)/(h*penalty);

    return cons_coef;
}

template <class TBaseElement>
std::pair<const double, const double> EmbeddedFluidElementDiscontinuous<TBaseElement>::ComputeTangentialPenaltyCoefficients(const EmbeddedDiscontinuousElementData& rData) const
{
    const double slip_length = rData.SlipLength;;
    const double penalty = rData.PenaltyCoefficient;

    const double h = rData.ElementSize;
    const double eff_mu = rData.EffectiveViscosity;
    const double coeff_1 = slip_length / (slip_length + penalty*h);
    const double coeff_2 = eff_mu / (slip_length + penalty*h);

    std::pair<const double, const double> pen_coeffs(coeff_1, coeff_2);

    return pen_coeffs;
}

template <class TBaseElement>
std::pair<const double, const double> EmbeddedFluidElementDiscontinuous<TBaseElement>::ComputeTangentialNitscheCoefficients(const EmbeddedDiscontinuousElementData& rData) const
{
    const double slip_length = rData.SlipLength;;
    const double penalty = rData.PenaltyCoefficient;

    const double h = rData.ElementSize;
    const double eff_mu = rData.EffectiveViscosity;
    const double coeff_1 = slip_length*penalty*h / (slip_length + penalty*h);
    const double coeff_2 = eff_mu*penalty*h / (slip_length + penalty*h);

    std::pair<const double, const double> pen_coeffs(coeff_1, coeff_2);

    return pen_coeffs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement);
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace EmbeddedDiscontinuousInternals {

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<2, 3>(const Element& rElement, const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Triangle2D3AusasModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<3, 4>(const Element& rElement, const Vector& rElementalDistances)
{
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4AusasModifiedShapeFunctions(rElement.pGetGeometry(), rElementalDistances));
}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

// template class EmbeddedFluidElementDiscontinuous< QSVMS< TimeIntegratedQSVMSData<2,3> > >;
// template class EmbeddedFluidElementDiscontinuous< QSVMS< TimeIntegratedQSVMSData<3,4> > >;

template class EmbeddedFluidElementDiscontinuous< SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > >;
template class EmbeddedFluidElementDiscontinuous< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}