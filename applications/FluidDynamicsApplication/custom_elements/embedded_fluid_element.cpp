#include "custom_elements/embedded_fluid_element.h"
//#include "custom_elements/qs_vms.h"
#include "custom_elements/symbolic_navier_stokes.h"

#include "custom_utilities/embedded_data.h"
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
        data.UpdateGeometryValues(data.PositiveSideWeights[g],
            row(data.PositiveSideN, g), data.PositiveSideDNDX[g]);

        this->AddTimeIntegratedSystem(
            data, rLeftHandSideMatrix, rRightHandSideVector);
    }

    if ( data.IsCut() ) {
        // Iterate over integration points on the boundary
        const unsigned int number_of_interface_gauss_points =
            data.PositiveInterfaceWeights.size();
        for (unsigned int g = 0; g < number_of_interface_gauss_points; g++) {
            data.UpdateGeometryValues(data.PositiveInterfaceWeights[g],
                row(data.PositiveInterfaceN, g), data.PositiveInterfaceDNDX[g]);
            this->AddBoundaryIntegral(data, data.PositiveInterfaceUnitNormals[g],
                rLeftHandSideMatrix, rRightHandSideVector);
        }

        // First, compute and assemble the penalty level set BC imposition contribution
        // Secondly, compute and assemble the modified Nitche method level set BC imposition contribution (Codina and Baiges, 2009)
        // Note that the Nistche contribution has to be computed the last since it drops the outer nodes rows previous constributions
        AddBoundaryConditionPenaltyContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
        DropOuterNodesVelocityContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
        AddBoundaryConditionModifiedNitscheContribution(rLeftHandSideMatrix, rRightHandSideVector, data);
    }
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
    const double tolerance = std::pow(1e-3 * this->ElementSize(),Dim-1);
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
void EmbeddedFluidElement<TBaseElement>::AddBoundaryConditionPenaltyContribution(
    MatrixType& rLHS,
    VectorType& rRHS,
    const EmbeddedElementData& rData) const {
    
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);
    // Obtain the previous iteration velocity solution

    // Set the penalty matrix
    bounded_matrix<double,NumNodes,NumNodes> p_gamma = ZeroMatrix(NumNodes, NumNodes);

    const unsigned int number_of_interface_gauss_points = rData.PositiveInterfaceWeights.size();

    for (unsigned int g = 0; g < number_of_interface_gauss_points; ++g) {
        const double weight = rData.PositiveInterfaceWeights[g];
        const auto shape_functions = row(rData.PositiveInterfaceN,g);
        p_gamma += weight*outer_prod(shape_functions,shape_functions);
    }

    // Multiply the penalty matrix by the penalty coefficient
    double penalty_coefficient = this->ComputePenaltyCoefficient(rData);
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
    const EmbeddedElementData& rData) const {

    // Compute the intersection area using the Gauss pts. weights
    double intersection_area = 0.0;
    for (unsigned int g = 0; g < rData.PositiveInterfaceWeights.size(); ++g) {
        intersection_area += rData.PositiveInterfaceWeights[g];
    }

    // Compute the element average values
    double avg_rho = 0.0;
    double avg_visc = 0.0;
    array_1d<double, Dim> avg_vel(Dim,0.0);

    for (unsigned int i = 0; i < NumNodes; ++i) {
        avg_rho += rData.Density[i];
        avg_visc += rData.DynamicViscosity[i];
        avg_vel += row(rData.Velocity, i);
    }

    constexpr double weight = 1./double(NumNodes);
    avg_rho *= weight;
    avg_visc *= weight;
    avg_vel *= weight;

    const double v_norm = norm_2(avg_vel);

    // Compute the penalty constant
    double h = this->ElementSize();
    const double pen_cons = avg_rho*std::pow(h, Dim)/rData.DeltaTime +
                                avg_rho*avg_visc*std::pow(h,Dim-2) +
                                avg_rho*v_norm*std::pow(h, Dim-1);

    // Return the penalty coefficient
    constexpr double K = 10.0;
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
            aux_out(i_out) = aux_cut(i_out_nodeid);
        }

        for (unsigned int i_int = 0; i_int < rData.NumPositiveNodes; ++i_int) {
            const unsigned int i_int_nodeid = rData.PositiveIndices[i_int];
            aux_int(i_int) = aux_cut(i_int_nodeid);
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

    // LHS outside Nitche contribution assembly
    noalias(rLHS) += nitsche_lhs;

    // RHS outside Nitche contribution assembly
    // Note that since we work with a residualbased formulation, the RHS is f_gamma - LHS*prev_sol
    noalias(rRHS) -= prod(nitsche_lhs, values);

    // Compute f_gamma if level set velocity is not 0
    if (this->Has(EMBEDDED_VELOCITY)) {
        nitsche_lhs.clear();

        const array_1d<double, 3 >& embedded_vel = this->GetValue(EMBEDDED_VELOCITY);
        array_1d<double, LocalSize> aux_embedded_vel(LocalSize,0.0);

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

//template class EmbeddedFluidElement< QSVMS< TimeIntegratedQSVMSData<3,4> > >;
template class EmbeddedFluidElement< SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > >;
template class EmbeddedFluidElement< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > >;

}