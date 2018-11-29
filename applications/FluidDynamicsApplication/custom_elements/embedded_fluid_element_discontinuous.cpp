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
        AddNormalSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, data); // NOTE: IMPLEMENT THE SKEW-SYMMETRIC ADJOINT IF IT IS NEEDED IN THE FUTURE. CREATE A IS_SKEW_SYMMETRIC ELEMENTAL FLAG.
        AddTangentialPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
        AddTangentialSymmetricCounterpartContribution(rLeftHandSideMatrix, rRightHandSideVector, data); // NOTE: IMPLEMENT THE SKEW-SYMMETRIC ADJOINT IF IT IS NEEDED IN THE FUTURE. CREATE A IS_SKEW_SYMMETRIC ELEMENTAL FLAG.
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
    const EmbeddedDiscontinuousElementData& rData) const
{
    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);

    // If there is embedded velocity, substract it to the previous iteration solution
    if (this->Has(EMBEDDED_VELOCITY)){
        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        array_1d<double, LocalSize> embedded_vel_exp(LocalSize, 0.0);

        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                embedded_vel_exp(i*BlockSize + comp) = embedded_vel(comp);
            }
        }

        noalias(values) -= embedded_vel_exp;
    }

    // Compute the Nitsche normal imposition penalty coefficient
    const double pen_coef = this->ComputeNormalPenaltyCoefficient(rData);

    // Compute the positive side LHS and RHS contributions
    const unsigned int number_of_positive_interface_integration_points = rData.PositiveInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_positive_interface_integration_points; ++g) {
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
                        const double aux = pen_coef*weight*aux_N(i)*aux_unit_normal(m)*aux_unit_normal(n)*aux_N(j); 
                        rLHS(row, col) += aux;
                        rRHS(row) -= aux*values(col);
                    }
                }
            }
        }
    }

    // Compute the negative side LHS and RHS contributions
    const unsigned int number_of_negative_interface_integration_points = rData.NegativeInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_negative_interface_integration_points; ++g) {
        // Get the Gauss pt. data
        const double weight = rData.NegativeInterfaceWeights[g];
        const auto aux_N = row(rData.NegativeInterfaceN, g);
        const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

        // Compute the Gauss pt. LHS contribution
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int j = 0; j < NumNodes; ++j){
                for (unsigned int m = 0; m < Dim; ++m){
                    const unsigned int row = i * BlockSize + m;
                    for (unsigned int n = 0; n < Dim; ++n){
                        const unsigned int col = j * BlockSize + n;
                        const double aux = pen_coef*weight*aux_N(i)*aux_unit_normal(m)*aux_unit_normal(n)*aux_N(j);
                        rLHS(row, col) += aux;
                        rRHS(row) -= aux*values(col);
                    }
                }
            }
        }
    }
}

template <class TBaseElement>
void EmbeddedFluidElementDiscontinuous<TBaseElement>::AddNormalSymmetricCounterpartContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedDiscontinuousElementData& rData) const
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

    // Set if the shear stress term is adjoint consistent (1.0) or not (-1.0)
    const double adjoint_consistency = -1.0;

    // Set an auxiliar array to compute the LHS contribution
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS = ZeroMatrix(LocalSize, LocalSize);

    // Compute positive side LHS contribution
    const unsigned int number_of_positive_interface_integration_points = rData.PositiveInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_positive_interface_integration_points; ++g){
        // Get the Gauss pt. data
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto aux_N = row(rData.PositiveInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> &aux_DN_DX = rData.PositiveInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

        // Fill the pressure to Voigt notation operator normal projected matrix
        BoundedMatrix<double, LocalSize, Dim> trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(LocalSize, Dim);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                trans_pres_to_voigt_matrix_normal_op(i*BlockSize + Dim, comp) = aux_N(i)*aux_unit_normal(comp);
            }
        }

        // Set the shape functions auxiliar matrix
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
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

    // Compute negative side LHS contribution
    const unsigned int number_of_negative_interface_integration_points = rData.NegativeInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_negative_interface_integration_points; ++g){
        // Get the Gauss pt. data
        const double weight = rData.NegativeInterfaceWeights[g];
        const auto aux_N = row(rData.NegativeInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> &aux_DN_DX = rData.NegativeInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

        // Fill the pressure to Voigt notation operator normal projected matrix
        BoundedMatrix<double, LocalSize, Dim> trans_pres_to_voigt_matrix_normal_op = ZeroMatrix(LocalSize, Dim);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                trans_pres_to_voigt_matrix_normal_op(i*BlockSize + Dim, comp) = aux_N(i)*aux_unit_normal(comp);
            }
        }

        // Set the shape functions auxiliar matrix
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
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
void EmbeddedFluidElementDiscontinuous<TBaseElement>::AddTangentialPenaltyContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedDiscontinuousElementData& rData) const
{
    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);

    // Compute the Nitsche tangential imposition penalty coefficients
    std::pair<const double, const double> pen_coefs = this->ComputeTangentialPenaltyCoefficients(rData);

    // Declare auxiliar arrays
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_1 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution coming from the tangential component of the Cauchy stress vector
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_2 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution generated by the viscous shear force generated by the velocity

    // Compute positive side LHS contribution 
    const unsigned int number_of_positive_interface_integration_points = rData.PositiveInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_positive_interface_integration_points; ++g){
        // Get the Gauss pt. data
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto aux_N = row(rData.PositiveInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.PositiveInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

        // Set the shape functions auxiliar matrices
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
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

    // Compute negative side LHS contribution 
    const unsigned int number_of_negative_interface_integration_points = rData.NegativeInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_negative_interface_integration_points; ++g){
        // Get the Gauss pt. data
        const double weight = rData.NegativeInterfaceWeights[g];
        const auto aux_N = row(rData.NegativeInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.NegativeInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

        // Set the shape functions auxiliar matrices
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
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
void EmbeddedFluidElementDiscontinuous<TBaseElement>::AddTangentialSymmetricCounterpartContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedDiscontinuousElementData& rData) const
{
    // Obtain the previous iteration velocity solution
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData, values);

    // Set if the shear stress term is adjoint consistent (1.0) or not (-1.0)
    const double adjoint_consistency = -1.0;

    // Compute the coefficients
    std::pair<const double, const double> nitsche_coefs = this->ComputeTangentialNitscheCoefficients(rData);

    // Declare auxiliar arrays
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_1 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution coming from the tangential component of the Cauchy stress vector
    BoundedMatrix<double, LocalSize, LocalSize> aux_LHS_2 = ZeroMatrix(LocalSize, LocalSize); // Adds the contribution generated by the viscous shear force generated by the velocity

    // Compute positive side LHS contribution
    const unsigned int number_of_positive_interface_integration_points = rData.PositiveInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_positive_interface_integration_points; ++g){
        // Get the Gauss pt. data
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto aux_N = row(rData.PositiveInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.PositiveInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.PositiveInterfaceUnitNormals[g];

        // Set the shape functions auxiliar matrices
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
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

    // Compute negative side LHS contribution
    const unsigned int number_of_negative_interface_integration_points = rData.NegativeInterfaceWeights.size();
    for (unsigned int g = 0; g < number_of_negative_interface_integration_points; ++g){
        // Get the Gauss pt. data
        const double weight = rData.NegativeInterfaceWeights[g];
        const auto aux_N = row(rData.NegativeInterfaceN, g);
        const BoundedMatrix<double, NumNodes, Dim> aux_DN_DX = rData.NegativeInterfaceDNDX[g];
        const auto &aux_unit_normal = rData.NegativeInterfaceUnitNormals[g];

        // Set the shape functions auxiliar matrices
        BoundedMatrix<double, Dim, LocalSize> N_mat = ZeroMatrix(Dim, LocalSize);
        for (unsigned int i = 0; i < NumNodes; ++i){
            for (unsigned int comp = 0; comp < Dim; ++comp){
                N_mat(comp, i*BlockSize + comp) = aux_N(i);
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

template class EmbeddedFluidElementDiscontinuous< QSVMS< TimeIntegratedQSVMSData<2,3> > >;
template class EmbeddedFluidElementDiscontinuous< QSVMS< TimeIntegratedQSVMSData<3,4> > >;

template class EmbeddedFluidElementDiscontinuous< SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > >;
template class EmbeddedFluidElementDiscontinuous< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}