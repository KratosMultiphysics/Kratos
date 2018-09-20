#include "includes/kratos_flags.h"

#include "custom_elements/embedded_fluid_element.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/symbolic_navier_stokes.h"

#include "custom_utilities/embedded_data.h"
#include "custom_utilities/element_size_calculator.h"
#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::EmbeddedFluidElement(IndexType NewId):
    TBaseElement(NewId)
{}

template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::EmbeddedFluidElement(IndexType NewId, const NodesArrayType& ThisNodes):
    TBaseElement(NewId,ThisNodes)
{}


template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::EmbeddedFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry):
    TBaseElement(NewId,pGeometry)
{}

template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::EmbeddedFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties):
    TBaseElement(NewId,pGeometry,pProperties)
{}


template< class TBaseElement >
EmbeddedFluidElement<TBaseElement>::~EmbeddedFluidElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TBaseElement >
Element::Pointer EmbeddedFluidElement<TBaseElement>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<EmbeddedFluidElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TBaseElement >
Element::Pointer EmbeddedFluidElement<TBaseElement>::Create(IndexType NewId,Geometry<NodeType>::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_shared<EmbeddedFluidElement>(NewId, pGeom, pProperties);
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {

    // Resize and intialize output
    if (rLeftHandSideMatrix.size1() != LocalSize)
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

    if (rRightHandSideVector.size() != LocalSize)
        rRightHandSideVector.resize(LocalSize, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    EmbeddedElementData data;
    data.Initialize(*this, rCurrentProcessInfo);
    this->InitializeGeometryData(data);

    // Iterate over integration points on the volume
    const unsigned int number_of_positive_gauss_points =
        data.PositiveSideWeights.size();
    for (unsigned int g = 0; g < number_of_positive_gauss_points; g++) {
        data.UpdateGeometryValues(g, data.PositiveSideWeights[g],
            row(data.PositiveSideN, g), data.PositiveSideDNDX[g]);

        this->CalculateMaterialResponse(data);

        this->AddTimeIntegratedSystem(
            data, rLeftHandSideMatrix, rRightHandSideVector);
    }

    if ( data.IsCut() ) {
        // Iterate over integration points on the boundary
        const unsigned int number_of_interface_gauss_points =
            data.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < number_of_interface_gauss_points; g++) {
            data.UpdateGeometryValues(g + number_of_positive_gauss_points, data.PositiveInterfaceWeights[g],
                row(data.PositiveInterfaceN, g), data.PositiveInterfaceDNDX[g]);

            this->CalculateMaterialResponse(data);

            this->AddBoundaryIntegral(data, data.PositiveInterfaceUnitNormals[g],
                rLeftHandSideMatrix, rRightHandSideVector);
        }

        if (this->Is(SLIP)){
            // Nitsche Navier-Slip boundary condition implementation (Winter, 2018)
            AddSlipNormalPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo);
            AddSlipNormalSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, data); // NOTE: IMPLEMENT THE SKEW-SYMMETRIC ADJOINT IF IT IS NEEDED IN THE FUTURE. CREATE A IS_SKEW_SYMMETRIC ELEMENTAL FLAG.
            AddSlipTangentialPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo);
            AddSlipTangentialSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo); // NOTE: IMPLEMENT THE SKEW-SYMMETRIC ADJOINT IF IT IS NEEDED IN THE FUTURE. CREATE A IS_SKEW_SYMMETRIC ELEMENTAL FLAG.
        } else {
            // First, compute and assemble the penalty level set BC imposition contribution
            // Secondly, compute and assemble the modified Nitsche method level set BC imposition contribution (Codina and Baiges, 2009)
            // Note that the Nistche contribution has to be computed the last since it drops the outer nodes rows previous constributions
            AddBoundaryConditionPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data, rCurrentProcessInfo);
            DropOuterNodesVelocityContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
            AddBoundaryConditionModifiedNitscheContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
        }
    }
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::Calculate(
    const Variable<double> &rVariable,
    double& rOutput,
    const ProcessInfo &rCurrentProcessInfo) {

    rOutput = 0.0;

    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo) {

    rOutput = ZeroVector(3);

    // If the element is split, integrate sigma.n over the interface
    // Note that in the ausas formulation, both interface sides need to be integrated
    if (rVariable == DRAG_FORCE) {

        EmbeddedElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);
        const unsigned int number_of_positive_gauss_points = data.PositiveSideWeights.size();

        if ( data.IsCut() ){
            // Integrate positive interface side drag
            const unsigned int n_int_pos_gauss = data.PositiveInterfaceWeights.size();
            for (unsigned int g = 0; g < n_int_pos_gauss; ++g) {

                // Update the Gauss pt. data
                data.UpdateGeometryValues(g + number_of_positive_gauss_points,
                   data.PositiveInterfaceWeights[g],row(data.PositiveInterfaceN, g),data.PositiveInterfaceDNDX[g]);

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
        }

    } else {
        TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::Calculate(
    const Variable<Vector> &rVariable,
    Vector& rOutput,
    const ProcessInfo &rCurrentProcessInfo) {

    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::Calculate(
    const Variable<Matrix> &rVariable,
    Matrix& rOutput,
    const ProcessInfo &rCurrentProcessInfo) {

    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <class TBaseElement>
int EmbeddedFluidElement<TBaseElement>::Check(
    const ProcessInfo& rCurrentProcessInfo) {

    int out = EmbeddedElementData::Check(*this,rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    return TBaseElement::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TBaseElement>
std::string EmbeddedFluidElement<TBaseElement>::Info() const {
    std::stringstream buffer;
    buffer << "EmbeddedFluidElement #" << this->Id();
    return buffer.str();
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::PrintInfo(
    std::ostream& rOStream) const {
    rOStream << "EmbeddedFluidElement" << Dim << "D" << NumNodes << "N"
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
void EmbeddedFluidElement<TBaseElement>::InitializeGeometryData(
    EmbeddedElementData& rData) const {

    rData.PositiveIndices.clear();
    rData.NegativeIndices.clear();

    // Number of positive and negative distance function values
    for (size_t i = 0; i < EmbeddedElementData::NumNodes; ++i) {

        if (rData.Distance[i] > 0.0) {
            rData.NumPositiveNodes++;
            rData.PositiveIndices.push_back(i);
        }
        else {
            rData.NumNegativeNodes++;
            rData.NegativeIndices.push_back(i);
        }
    }

    if ( rData.IsCut() ) {
        this->DefineCutGeometryData(rData);
    }
    else {
        this->DefineStandardGeometryData(rData);
    }
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::DefineStandardGeometryData(
    EmbeddedElementData& rData) const {

    this->CalculateGeometryData(
        rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);
    rData.NumPositiveNodes = NumNodes;
    rData.NumNegativeNodes = 0;
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::DefineCutGeometryData(
    EmbeddedElementData& rData) const {

    // Auxiliary distance vector for the element subdivision utility
    Vector distances = rData.Distance;

    ModifiedShapeFunctions::Pointer p_calculator =
        Internals::GetShapeFunctionCalculator<EmbeddedElementData::Dim,
            EmbeddedElementData::NumNodes>(*this, distances);

    // Fluid side
    p_calculator->ComputePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveSideN, rData.PositiveSideDNDX, rData.PositiveSideWeights,
        GeometryData::GI_GAUSS_2);

    // Fluid side interface
    p_calculator->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(
        rData.PositiveInterfaceN, rData.PositiveInterfaceDNDX,
        rData.PositiveInterfaceWeights, GeometryData::GI_GAUSS_2);

    // Fluid side interface normals
    p_calculator->ComputePositiveSideInterfaceAreaNormals(
        rData.PositiveInterfaceUnitNormals, GeometryData::GI_GAUSS_2);

    // Normalize the normals
    // Note: we calculate h here (and we don't use the value in rData.ElementSize)
    // because rData.ElementSize might still be uninitialized: some data classes define it at the Gauss point.
    double h = ElementSizeCalculator<Dim,NumNodes>::MinimumElementSize(this->GetGeometry());
    const double tolerance = std::pow(1e-3 * h,Dim-1);
    this->NormalizeInterfaceNormals(rData.PositiveInterfaceUnitNormals, tolerance);
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::NormalizeInterfaceNormals(
    typename EmbeddedElementData::InterfaceNormalsType& rNormals,
    double Tolerance) const {
    for (unsigned int i = 0; i < rNormals.size(); ++i) {
        double norm = norm_2(rNormals[i]);
        rNormals[i] /= std::max(norm,Tolerance);
    }
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::AddSlipNormalPenaltyContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedElementData& rData,
    const ProcessInfo &rCurrentProcessInfo) const
{
    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);

    // If there is embedded velocity, substract it to the previous iteration solution
    if (this->Has(EMBEDDED_VELOCITY)) {
        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        array_1d<double, LocalSize> embedded_vel_exp = ZeroVector(LocalSize);

        for (unsigned int i = 0; i < NumNodes; ++i) {
            for (unsigned int comp = 0; comp < Dim; ++comp) {
                embedded_vel_exp(i*BlockSize + comp) = embedded_vel(comp);
            }
        }

        noalias(values) -= embedded_vel_exp;
    }

    // Compute the Nitsche normal imposition penalty coefficient
    const double pen_coef = this->ComputeSlipNormalPenaltyCoefficient(rData, rCurrentProcessInfo);

    // Compute LHS contribution
    // BoundedMatrix<double, LocalSize, LocalSize> aux_LHS = ZeroMatrix(LocalSize, LocalSize);
    const unsigned int number_of_integration_points = rData.PositiveInterfaceWeights.size();

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        // Get the Gauss pt. data
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto aux_N = row(rData.PositiveInterfaceN, g);
        const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

        // Compute the Gauss pt. LHS contribution
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int j = 0; j < NumNodes; ++j){
                for (unsigned int m = 0; m < Dim; ++m){
                    const unsigned int row = i * BlockSize + m;
                    for (unsigned int n = 0; n < Dim; ++n){
                        const unsigned int col = j * BlockSize + n;
                        #ifdef KRATOS_USE_AMATRIX
                        double lhs_ij = pen_coef*weight*aux_N[i]*aux_unit_normal(m)*aux_unit_normal(n)*aux_N[j];
                        #else
                        double lhs_ij = pen_coef*weight*aux_N(i)*aux_unit_normal(m)*aux_unit_normal(n)*aux_N(j);
                        #endif
                        rLHS(row, col) += lhs_ij;
                        rRHS(row) -= lhs_ij*values(col);
                    }
                }
            }
        }
    }
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::AddSlipNormalSymmetricCounterpartContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedElementData& rData) const {

    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);

    // If there is embedded velocity, substract it to the previous iteration solution
    if (this->Has(EMBEDDED_VELOCITY)) {
        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        array_1d<double, LocalSize> embedded_vel_exp = ZeroVector(LocalSize);

        for (unsigned int i = 0; i < NumNodes; ++i) {
            for (unsigned int comp = 0; comp < Dim; ++comp) {
                embedded_vel_exp(i*BlockSize + comp) = embedded_vel(comp);
            }
        }

        noalias(values) -= embedded_vel_exp;
    }

    // Set if the shear stress term is adjoint consistent (1.0) or not (-1.0)
    const double adjoint_consistency = -1.0;

    // Compute LHS contribution
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS = ZeroMatrix(LocalSize, LocalSize);
    const unsigned int number_of_integration_points = rData.PositiveInterfaceWeights.size();

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        // Get the Gauss pt. data
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto aux_N = row(rData.PositiveInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> &aux_DN_DX = rData.PositiveInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

        // Fill the pressure to Voigt notation operator normal projected matrix
        BoundedMatrix<double, LocalSize, Dim> trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(LocalSize, Dim);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                #ifdef KRATOS_USE_AMATRIX
                trans_pres_to_voigt_matrix_normal_op(i*BlockSize + Dim, comp) = aux_N[i]*aux_unit_normal(comp);
                #else
                trans_pres_to_voigt_matrix_normal_op(i*BlockSize + Dim, comp) = aux_N(i)*aux_unit_normal(comp);
                #endif
            }
        }

        // Set the shape functions auxiliar matrix
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                #ifdef KRATOS_USE_AMATRIX
                N_mat(comp, i*BlockSize + comp) = aux_N[i];
                #else
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
                #endif
            }
        }

        // Set the current Gauss pt. strain matrix
        BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
        FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

        // Set the normal projection matrix (n x n)
        BoundedMatrix<double, Dim, Dim> normal_proj_matrix;
        FluidElementUtilities<NumNodes>::SetNormalProjectionMatrix(aux_unit_normal, normal_proj_matrix);

        // Get the normal projection matrix in Voigt notation
        BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
        FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

        // Compute some Gauss pt. auxiliar matrices
        const BoundedMatrix<double, LocalSize, StrainSize> aux_matrix_BC = prod(trans(B_matrix), trans(rData.C));
        const BoundedMatrix<double, StrainSize, Dim> aux_matrix_APnorm = prod(trans(voigt_normal_proj_matrix), normal_proj_matrix);
        const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BCAPnorm = prod(aux_matrix_BC, aux_matrix_APnorm);

        // Contribution coming fron the shear stress operator
        noalias(aux_LHS) -= adjoint_consistency*weight*prod(aux_matrix_BCAPnorm, N_mat);

        // Contribution coming from the pressure terms
        const BoundedMatrix<double, LocalSize, Dim> aux_matrix_VPnorm = prod(trans_pres_to_voigt_matrix_normal_op, normal_proj_matrix);
        noalias(aux_LHS) -= weight*prod(aux_matrix_VPnorm, N_mat);
    }

    // LHS outside Nitsche contribution assembly
    noalias(rLHS) += aux_LHS;

    // RHS outside Nitsche contribution assembly
    // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
    noalias(rRHS) -= prod(aux_LHS, values);
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::AddSlipTangentialPenaltyContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedElementData& rData,
    const ProcessInfo &rCurrentProcessInfo) const
{
    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);

    // Compute the Nitsche tangential imposition penalty coefficients
    std::pair<const double, const double> pen_coefs = this->ComputeSlipTangentialPenaltyCoefficients(rData, rCurrentProcessInfo);

    // Declare auxiliar arrays
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_1 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution coming from the tangential component of the Cauchy stress vector
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_2 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution generated by the viscous shear force generated by the velocity
    const unsigned int number_of_integration_points = rData.PositiveInterfaceWeights.size();

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        // Get the Gauss pt. data
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto aux_N = row(rData.PositiveInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.PositiveInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

        // Set the shape functions auxiliar matrices
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                #ifdef KRATOS_USE_AMATRIX
                N_mat(comp, i*BlockSize + comp) = aux_N[i];
                #else
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
                #endif
            }
        }
        BoundedMatrix<double, LocalSize, Dim> N_mat_trans = trans(N_mat);

        // Set the tangential projection matrix (I - n x n)
        BoundedMatrix<double, Dim, Dim> tang_proj_matrix;
        FluidElementUtilities<NumNodes>::SetTangentialProjectionMatrix(aux_unit_normal, tang_proj_matrix);

        // Set the current Gauss pt. strain matrix
        BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
        FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

        // Get the normal projection matrix in Voigt notation
        BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
        FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

        // Compute some Gauss pt. auxiliar matrices
        const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_CB = prod(rData.C, B_matrix);
        const BoundedMatrix<double, StrainSize, Dim> aux_matrix_PtangA = prod(tang_proj_matrix, voigt_normal_proj_matrix);
        const BoundedMatrix<double, LocalSize, Dim> aux_matrix_PtangACB = prod(aux_matrix_PtangA, aux_matrix_CB);

        // Contribution coming from the traction vector tangencial component
        noalias(aux_LHS_1) += pen_coefs.first*weight*prod(N_mat_trans, aux_matrix_PtangACB);

        // Contribution coming from the shear force generated by the velocity jump
        const BoundedMatrix<double, LocalSize, Dim> aux_matrix_N_trans_tang = prod(N_mat_trans, tang_proj_matrix);
        noalias(aux_LHS_2) += pen_coefs.second*weight*prod(aux_matrix_N_trans_tang, N_mat);
    }

    // LHS outside Nitsche contribution assembly
    noalias(rLHS) += aux_LHS_1;
    noalias(rLHS) += aux_LHS_2;

    // RHS outside Nitsche contribution assembly
    // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
    noalias(rRHS) -= prod(aux_LHS_1, values);
    noalias(rRHS) -= prod(aux_LHS_2, values);

    // If level set velocity is not 0, add its contribution to the RHS
    if (this->Has(EMBEDDED_VELOCITY)) {
        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        array_1d<double, LocalSize> embedded_vel_exp = ZeroVector(LocalSize);

        for (unsigned int i = 0; i < NumNodes; ++i) {
            for (unsigned int comp = 0; comp < Dim; ++comp) {
                embedded_vel_exp(i*BlockSize + comp) = embedded_vel(comp);
            }
        }
        noalias(rRHS) += prod(aux_LHS_2, embedded_vel_exp);
    }
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::AddSlipTangentialSymmetricCounterpartContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedElementData& rData,
    const ProcessInfo &rCurrentProcessInfo) const
{
    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);

    // Set if the shear stress term is adjoint consistent (1.0) or not (-1.0)
    const double adjoint_consistency = -1.0;

    // Compute the coefficients
    std::pair<const double, const double> nitsche_coefs = this->ComputeSlipTangentialNitscheCoefficients(rData, rCurrentProcessInfo);

    // Declare auxiliar arrays
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_1 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution coming from the tangential component of the Cauchy stress vector
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_2 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution generated by the viscous shear force generated by the velocity

    const unsigned int number_of_integration_points = rData.PositiveInterfaceWeights.size();

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        // Get the Gauss pt. data
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto aux_N = row(rData.PositiveInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.PositiveInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

        // Set the shape functions auxiliar matrices
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                #ifdef KRATOS_USE_AMATRIX
                N_mat(comp, i*BlockSize + comp) = aux_N[i];
                #else
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
                #endif
            }
        }

        // Set the current Gauss pt. strain matrix
        BoundedMatrix<double, StrainSize, LocalSize> B_matrix = ZeroMatrix(StrainSize, LocalSize);
        FluidElementUtilities<NumNodes>::GetStrainMatrix(aux_DN_DX, B_matrix);

        // Set the tangential projection matrix (I - n x n)
        BoundedMatrix<double, Dim, Dim> tang_proj_matrix;
        FluidElementUtilities<NumNodes>::SetTangentialProjectionMatrix(aux_unit_normal, tang_proj_matrix);

        // Get the normal projection matrix in Voigt notation
        BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix = ZeroMatrix(Dim, StrainSize);
        FluidElementUtilities<NumNodes>::VoigtTransformForProduct(aux_unit_normal, voigt_normal_proj_matrix);

        // Compute some Gauss pt. auxiliar matrices
        const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BtransAtrans = prod(trans(B_matrix), trans(voigt_normal_proj_matrix));
        const BoundedMatrix<double, LocalSize, Dim> aux_matrix_BtransAtransPtan = prod(aux_matrix_BtransAtrans, tang_proj_matrix);
        const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_CB = prod(rData.C, B_matrix);
        const BoundedMatrix<double, Dim, LocalSize> aux_matrix_ACB = prod(voigt_normal_proj_matrix, aux_matrix_CB);
        const BoundedMatrix<double, LocalSize, LocalSize> aux_matrix_BtransAtransPtanACB = prod(aux_matrix_BtransAtransPtan, aux_matrix_ACB);

        // Contribution coming from the traction vector tangencial component
        noalias(aux_LHS_1) -= adjoint_consistency*nitsche_coefs.first*weight*aux_matrix_BtransAtransPtanACB;

        // Contribution coming from the shear force generated by the velocity jump
        noalias(aux_LHS_2) -= adjoint_consistency*nitsche_coefs.second*weight*prod(aux_matrix_BtransAtransPtan, N_mat);

    }

    // LHS outside Nitsche contribution assembly
    noalias(rLHS) += aux_LHS_1;
    noalias(rLHS) += aux_LHS_2;

    // RHS outside Nitsche contribution assembly
    // If level set velocity is not 0, add its contribution to the RHS
    if (this->Has(EMBEDDED_VELOCITY)) {
        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        array_1d<double, LocalSize> embedded_vel_exp = ZeroVector(LocalSize);

        for (unsigned int i = 0; i < NumNodes; ++i) {
            for (unsigned int comp = 0; comp < Dim; ++comp) {
                embedded_vel_exp(i*BlockSize + comp) = embedded_vel(comp);
            }
        }

        noalias(rRHS) += prod(aux_LHS_2, embedded_vel_exp);
    }

    // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
    noalias(rRHS) -= prod(aux_LHS_1, values);
    noalias(rRHS) -= prod(aux_LHS_2, values);

}

template <class TBaseElement>
double EmbeddedFluidElement<TBaseElement>::ComputeSlipNormalPenaltyCoefficient(
    const EmbeddedElementData& rData,
    const ProcessInfo &rCurrentProcessInfo) const
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
    const double avg_rho = rData.Density;
    const double eff_mu = rData.EffectiveViscosity;
    const double h = rData.ElementSize;
    const double penalty = 1.0/rCurrentProcessInfo[PENALTY_COEFFICIENT];
    const double cons_coef = (eff_mu + eff_mu + avg_rho*v_norm*h + avg_rho*h*h/rData.DeltaTime)/(h*penalty);

    return cons_coef;
}

template <class TBaseElement>
std::pair<const double, const double> EmbeddedFluidElement<TBaseElement>::ComputeSlipTangentialPenaltyCoefficients(
    const EmbeddedElementData& rData,
    const ProcessInfo &rCurrentProcessInfo) const
{
    const double slip_length = rCurrentProcessInfo[SLIP_LENGTH];
    const double penalty = 1.0/rCurrentProcessInfo[PENALTY_COEFFICIENT];

    const double eff_mu = rData.EffectiveViscosity;
    const double h = rData.ElementSize;
    const double coeff_1 = slip_length / (slip_length + penalty*h);
    const double coeff_2 = eff_mu / (slip_length + penalty*h);

    std::pair<const double, const double> pen_coeffs(coeff_1, coeff_2);

    return pen_coeffs;
}

template <class TBaseElement>
std::pair<const double, const double> EmbeddedFluidElement<TBaseElement>::ComputeSlipTangentialNitscheCoefficients(
    const EmbeddedElementData& rData,
    const ProcessInfo &rCurrentProcessInfo) const 
{
    const double slip_length = rCurrentProcessInfo[SLIP_LENGTH];
    const double penalty = 1.0/rCurrentProcessInfo[PENALTY_COEFFICIENT];

    const double eff_mu = rData.EffectiveViscosity;
    const double h = rData.ElementSize;
    const double coeff_1 = slip_length*penalty*h / (slip_length + penalty*h);
    const double coeff_2 = eff_mu*penalty*h / (slip_length + penalty*h);

    std::pair<const double, const double> pen_coeffs(coeff_1, coeff_2);

    return pen_coeffs;
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::AddBoundaryConditionPenaltyContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedElementData& rData,
    const ProcessInfo &rCurrentProcessInfo) const
{
    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);

    // Set the penalty matrix
    BoundedMatrix<double,NumNodes,NumNodes> p_gamma = ZeroMatrix(NumNodes, NumNodes);

    const unsigned int number_of_interface_gauss_points = rData.PositiveInterfaceWeights.size();

    for (unsigned int g = 0; g < number_of_interface_gauss_points; ++g) {
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto shape_functions = row(rData.PositiveInterfaceN,g);
        p_gamma += weight*outer_prod(shape_functions,shape_functions);
    }

    // Multiply the penalty matrix by the penalty coefficient
    double penalty_coefficient = this->ComputePenaltyCoefficient(rData, rCurrentProcessInfo);
    p_gamma *= penalty_coefficient;

    MatrixType penalty_lhs = ZeroMatrix(LocalSize, LocalSize);

    // LHS penalty contribution assembly (symmetric mass matrix)
    for (unsigned int i = 0; i<NumNodes; i++) {
        // Diagonal terms
        for (unsigned int comp = 0; comp<Dim; comp++) {
            penalty_lhs(i*BlockSize+comp, i*BlockSize+comp) = p_gamma(i,i);
        }
        // Off-diagonal terms
        for (unsigned int j = i+1; j<NumNodes; j++) {
            for (unsigned int comp = 0; comp<Dim; comp++) {
                penalty_lhs(i*BlockSize+comp, j*BlockSize+comp) = p_gamma(i,j);
                penalty_lhs(j*BlockSize+comp, i*BlockSize+comp) = p_gamma(i,j);
            }
        }
    }

    noalias(rLHS) += penalty_lhs;

    // RHS penalty contribution assembly
    if (this->Has(EMBEDDED_VELOCITY)) {
        VectorType penalty_rhs = ZeroVector(LocalSize);
        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        array_1d<double, LocalSize> aux_embedded_vel = ZeroVector(LocalSize);

        for (unsigned int i=0; i<NumNodes; i++) {
            for (unsigned int comp=0; comp<Dim; comp++) {
                aux_embedded_vel(i*BlockSize+comp) = embedded_vel(comp);
            }
        }

        noalias(rRHS) += prod(penalty_lhs, aux_embedded_vel);
    }

    noalias(rRHS) -= prod(penalty_lhs, values); // Residual contribution assembly
}

template <class TBaseElement>
double EmbeddedFluidElement<TBaseElement>::ComputePenaltyCoefficient(
    const EmbeddedElementData& rData,
    const ProcessInfo &rCurrentProcessInfo) const
{
    // Compute the intersection area using the Gauss pts. weights
    double intersection_area = 0.0;
    for (unsigned int g = 0; g < rData.PositiveInterfaceWeights.size(); ++g) {
        intersection_area += rData.PositiveInterfaceWeights[g];
    }

    // Compute the element average velocity value
    array_1d<double, Dim> avg_vel = ZeroVector(Dim);

    for (unsigned int i = 0; i < NumNodes; ++i) {
        avg_vel += row(rData.Velocity, i);
    }

    constexpr double weight = 1./double(NumNodes);
    avg_vel *= weight;

    const double v_norm = norm_2(avg_vel);

    // Compute the penalty constant
    double h = rData.ElementSize;
    const double rho = rData.Density;
    const double eff_mu = rData.EffectiveViscosity;
    const double pen_cons = rho*std::pow(h, Dim)/rData.DeltaTime +
                                rho*eff_mu*std::pow(h,Dim-2) +
                                rho*v_norm*std::pow(h, Dim-1);

    // Return the penalty coefficient
    const double K = rCurrentProcessInfo[PENALTY_COEFFICIENT];
    const double pen_coef = K * pen_cons / intersection_area;

    return pen_coef;
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::DropOuterNodesVelocityContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedElementData& rData) const {
    // Set the LHS and RHS u_out rows to zero (outside nodes used to impose the BC)
    for (unsigned int i=0; i<rData.NumNegativeNodes ; ++i) {
        const unsigned int out_node_row_id = rData.NegativeIndices[i];

        for (unsigned int j=0; j<Dim; ++j) {
            // LHS matrix u_out zero set (note that just the velocity rows are set to 0)
            for (unsigned int col = 0; col<LocalSize; col++) {
                rLHS(out_node_row_id*BlockSize+j, col) = 0.0;
            }

            // RHS vector u_out zero set (note that just the velocity rows are set to 0)
            rRHS(out_node_row_id*BlockSize+j) = 0.0;
        }
    }
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::AddBoundaryConditionModifiedNitscheContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedElementData& rData) const {

    // Obtain the previous iteration velocity solution
    array_1d<double, LocalSize> values;
    this->GetCurrentValuesVector(rData, values);

    // Compute the BCs imposition matrices
    MatrixType M_gamma = ZeroMatrix(rData.NumNegativeNodes, rData.NumNegativeNodes);  // Outside nodes matrix (Nitsche contribution)
    MatrixType N_gamma = ZeroMatrix(rData.NumNegativeNodes, rData.NumPositiveNodes);  // Interior nodes matrix (Nitsche contribution)
    MatrixType f_gamma = ZeroMatrix(rData.NumNegativeNodes, NumNodes);    // Matrix to compute the RHS (Nitsche contribution)

    VectorType aux_out(rData.NumNegativeNodes);
    VectorType aux_int(rData.NumPositiveNodes);

    const unsigned int number_of_integration_points = rData.PositiveInterfaceWeights.size();

    for (unsigned int g = 0; g < number_of_integration_points; g++) {
        const double weight = rData.PositiveInterfaceWeights[g];

        const auto aux_cut = row(rData.PositiveInterfaceN, g);

        for (unsigned int i_out = 0; i_out < rData.NumNegativeNodes; i_out++) {
            const unsigned int i_out_nodeid = rData.NegativeIndices[i_out];
            #ifdef KRATOS_USE_AMATRIX
            aux_out(i_out) = aux_cut[i_out_nodeid];
            #else
            aux_out(i_out) = aux_cut(i_out_nodeid);
            #endif
        }

        for (unsigned int i_int = 0; i_int < rData.NumPositiveNodes; ++i_int) {
            const unsigned int i_int_nodeid = rData.PositiveIndices[i_int];
            #ifdef KRATOS_USE_AMATRIX
            aux_int(i_int) = aux_cut[i_int_nodeid];
            #else
            aux_int(i_int) = aux_cut(i_int_nodeid);
            #endif
        }

        M_gamma += weight*outer_prod(aux_out,aux_out);
        N_gamma += weight*outer_prod(aux_out,aux_int);
        f_gamma += weight*outer_prod(aux_out,aux_cut);
    }

    MatrixType nitsche_lhs = ZeroMatrix(LocalSize, LocalSize);

    // LHS outside nodes contribution assembly
    // Outer nodes contribution assembly
    for (unsigned int i = 0; i<rData.NumNegativeNodes; i++) {
        unsigned int out_node_row_id = rData.NegativeIndices[i];

        for (unsigned int j = 0; j<rData.NumNegativeNodes; j++) {
            unsigned int out_node_col_id = rData.NegativeIndices[j];

            for (unsigned int comp = 0; comp<Dim; comp++) {
                nitsche_lhs(out_node_row_id*BlockSize+comp, out_node_col_id*BlockSize+comp) = M_gamma(i, j);
            }
        }
    }

    // Interior nodes contribution assembly
    for (unsigned int i = 0; i<rData.NumNegativeNodes; i++) {
        unsigned int out_node_row_id = rData.NegativeIndices[i];

        for (unsigned int j = 0; j<rData.NumPositiveNodes; j++) {
            unsigned int int_node_col_id = rData.PositiveIndices[j];

            for (unsigned int comp = 0; comp<Dim; comp++) {
                nitsche_lhs(out_node_row_id*BlockSize+comp, int_node_col_id*BlockSize+comp) = N_gamma(i, j);
            }
        }
    }

    // LHS outside Nitsche contribution assembly
    noalias(rLHS) += nitsche_lhs;

    // RHS outside Nitsche contribution assembly
    // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
    noalias(rRHS) -= prod(nitsche_lhs, values);

    // Compute f_gamma if level set velocity is not 0
    if (this->Has(EMBEDDED_VELOCITY)) {
        nitsche_lhs.clear();

        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        array_1d<double, LocalSize> aux_embedded_vel = ZeroVector(LocalSize);

        for (unsigned int i=0; i<NumNodes; i++) {
            aux_embedded_vel(i*BlockSize) = embedded_vel(0);
            aux_embedded_vel(i*BlockSize+1) = embedded_vel(1);
            aux_embedded_vel(i*BlockSize+2) = embedded_vel(2);
        }

        // Asemble the RHS f_gamma contribution
        for (unsigned int i=0; i<rData.NumNegativeNodes; i++) {
            unsigned int out_node_row_id = rData.NegativeIndices[i];

            for (unsigned int j=0; j<NumNodes; j++) {
                for (unsigned int comp = 0; comp<Dim; comp++) {
                    nitsche_lhs(out_node_row_id*BlockSize+comp, j*BlockSize+comp) = f_gamma(i,j);
                }
            }
        }

        noalias(rRHS) += prod(nitsche_lhs, aux_embedded_vel);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::save(Serializer& rSerializer) const {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement);
}

template <class TBaseElement>
void EmbeddedFluidElement<TBaseElement>::load(Serializer& rSerializer) {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for template specialization
///////////////////////////////////////////////////////////////////////////////////////////////////

namespace Internals {

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<2, 3>(
    const Element& rElement, const Vector& rDistance) {
    return ModifiedShapeFunctions::Pointer(new Triangle2D3ModifiedShapeFunctions(rElement.pGetGeometry(),rDistance));
}

template <>
ModifiedShapeFunctions::Pointer GetShapeFunctionCalculator<3, 4>(
    const Element& rElement, const Vector& rDistance) {
    return ModifiedShapeFunctions::Pointer(new Tetrahedra3D4ModifiedShapeFunctions(rElement.pGetGeometry(),rDistance));
}

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class EmbeddedFluidElement< QSVMS< TimeIntegratedQSVMSData<2,3> > >;
template class EmbeddedFluidElement< QSVMS< TimeIntegratedQSVMSData<3,4> > >;

template class EmbeddedFluidElement< SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > >;
template class EmbeddedFluidElement< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}