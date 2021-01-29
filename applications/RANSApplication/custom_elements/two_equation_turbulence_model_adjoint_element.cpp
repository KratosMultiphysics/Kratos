//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "utilities/adjoint_extensions.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/indirect_scalar.h"

// Application includes
#include "custom_constitutive/fluid_adjoint_constitutive_law.h"
#include "fluid_dynamics_application_variables.h"

#include "custom_elements/data_containers/k_epsilon/qsvms_rfc_adjoint_element_data.h"

// Include base h
#include "two_equation_turbulence_model_adjoint_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::ThisExtensions(Element* pElement)
    : mpElement{pElement}
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    // auto& r_node = mpElement->GetGeometry()[NodeId];
    // rVector.resize(TDim + 1);

    // const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    // for (unsigned int i = 0; i < TDim; ++i) {
    //     rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative(), Step);
    // }

    // rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    // auto& r_node = mpElement->GetGeometry()[NodeId];
    // rVector.resize(TDim + 1);

    // const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    // for (unsigned int i = 0; i < TDim; ++i) {
    //     rVector[i] = MakeIndirectScalar(
    //         r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative(), Step);
    // }

    // rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    // auto& r_node = mpElement->GetGeometry()[NodeId];
    // rVector.resize(TDim + 1);

    // const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    // for (unsigned int i = 0; i < TDim; ++i) {
    //     rVector[i] = MakeIndirectScalar(
    //         r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(),
    //         Step);
    // }

    // rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    // rVariables.resize(1);
    // rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    // rVariables.resize(1);
    // rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    // rVariables.resize(1);
    // rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::TwoEquationTurbulenceModelAdjointElement(IndexType NewId)
    : BaseType(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::TwoEquationTurbulenceModelAdjointElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::TwoEquationTurbulenceModelAdjointElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::~TwoEquationTurbulenceModelAdjointElement()
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::Clone(
    IndexType NewId,
    NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
int TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    // return TAdjointElementData::Check(this->GetGeometry(), rCurrentProcessInfo);
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::EquationIdVector(
    EquationIdVectorType& rElementalEquationIdList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // if (rElementalEquationIdList.size() != TElementLocalSize) {
    //     rElementalEquationIdList.resize(TElementLocalSize, false);
    // }

    // const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

    // IndexType local_index = 0;
    // for (IndexType i = 0; i < TNumNodes; ++i) {
    //     const auto& r_node = this->GetGeometry()[i];
    //     for (const auto p_variable : r_variables_list) {
    //         rElementalEquationIdList[local_index++] =
    //             r_node.GetDof(*p_variable).EquationId();
    //     }
    // }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    // if (rElementalDofList.size() != TElementLocalSize) {
    //     rElementalDofList.resize(TElementLocalSize);
    // }

    // const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

    // IndexType local_index = 0;
    // for (IndexType i = 0; i < TNumNodes; ++i) {
    //     const auto& r_node = this->GetGeometry()[i];
    //     for (const auto p_variable : r_variables_list) {
    //         rElementalDofList[local_index++] = r_node.pGetDof(*p_variable);
    //     }
    // }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetValuesVector(
    VectorType& rValues,
    int Step) const
{
    // if (rValues.size() != TElementLocalSize) {
    //     rValues.resize(TElementLocalSize, false);
    // }

    // const auto& r_geometry = this->GetGeometry();
    // IndexType local_index = 0;
    // for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
    //     const auto& r_node = r_geometry[i_node];
    //     const auto& r_velocity =
    //         r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, Step);
    //     for (IndexType d = 0; d < TDim; ++d) {
    //         rValues[local_index++] = r_velocity[d];
    //     }
    //     rValues[local_index++] =
    //         r_node.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
    // }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetFirstDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    // if (rValues.size() != TElementLocalSize) {
    //     rValues.resize(TElementLocalSize, false);
    // }
    // rValues.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetSecondDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    // if (rValues.size() != TElementLocalSize) {
    //     rValues.resize(TElementLocalSize);
    // }

    // const auto& r_geometry = this->GetGeometry();
    // IndexType local_index = 0;
    // for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
    //     const auto& r_acceleration =
    //         r_geometry[i_node].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, Step);
    //     for (IndexType d = 0; d < TDim; ++d) {
    //         rValues[local_index++] = r_acceleration[d];
    //     }
    //     rValues[local_index++] = 0.0; // pressure dof
    // }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // If we are restarting, the constitutive law will be already defined
    if (mpFluidConstitutiveLaw == nullptr) {
        const Properties& r_properties = this->GetProperties();

        KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
            << "In initialization of Element " << this->Info()
            << ": No CONSTITUTIVE_LAW defined for property "
            << r_properties.Id() << "." << std::endl;

        // Here we can do down casting because, it should be always a FluidConstitutiveLaw
        mpFluidConstitutiveLaw = std::static_pointer_cast<FluidConstitutiveLaw>(
            r_properties[CONSTITUTIVE_LAW]->Clone());

        const GeometryType& r_geometry = this->GetGeometry();
        const auto& r_shape_functions =
            r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

        mpFluidConstitutiveLaw->InitializeMaterial(r_properties, r_geometry,
                                                   row(r_shape_functions, 0));
    }

    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize and initialize output
    if (rLeftHandSideMatrix.size1() != TElementLocalSize)
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

    if (rRightHandSideVector.size() != TElementLocalSize)
        rRightHandSideVector.resize(TElementLocalSize, false);

    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "TwoEquationTurbulenceModelAdjointElement::"
                    "CalculateRightHandSide method is not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLocalVelocityContribution(
    MatrixType &rDampMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    // KRATOS_TRY

    // Vector gauss_weights;
    // Matrix shape_functions;
    // ShapeFunctionDerivativesArrayType shape_derivatives;
    // this->CalculateGeometryData(gauss_weights, shape_functions,
    //                             shape_derivatives, TAdjointElementData::GetIntegrationMethod());

    // typename TAdjointElementData::Primal::Data element_data(*this, *mpFluidConstitutiveLaw);
    // typename TAdjointElementData::Primal::ResidualContributions residual_contributions(element_data);

    // residual_contributions.Initialize(rDampMatrix, rCurrentProcessInfo);
    // element_data.Initialize(rCurrentProcessInfo);

    // BoundedVector<double, TElementLocalSize> residual;
    // residual.clear();

    // for (IndexType g = 0; g < gauss_weights.size(); ++g) {
    //     const Vector& N = row(shape_functions, g);
    //     const Matrix& dNdX = shape_derivatives[g];
    //     const double weight = gauss_weights[g];

    //     element_data.CalculateGaussPointData(weight, N, dNdX);
    //     residual_contributions.AddResidualContribution(residual, weight, N, dNdX);
    // }

    // noalias(rRightHandSideVector) = residual;

    // KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
    //     rLeftHandSideMatrix.size2() != TElementLocalSize) {
    //     rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    // }

    // rLeftHandSideMatrix.clear();

    // AddFluidFirstDerivatives<typename TAdjointElementData::StateDerivatives::FirstDerivatives>(
    //     rLeftHandSideMatrix, TAdjointElementData::GetIntegrationMethod(), rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
    //     rLeftHandSideMatrix.size2() != TElementLocalSize) {
    //     rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    // }

    // rLeftHandSideMatrix.clear();

    // AddFluidSecondDerivatives<typename TAdjointElementData::StateDerivatives::SecondDerivatives>(
    //     rLeftHandSideMatrix, TAdjointElementData::GetIntegrationMethod(), rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    rMassMatrix.resize(0, 0);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "TwoEquationTurbulenceModelAdjointElement::"
                    "CalculateDampingMatrix method is not implemented.";

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    // KRATOS_TRY

    // if (rSensitivityVariable == SHAPE_SENSITIVITY) {
    //     using derivatives_type = typename TAdjointElementData::SensitivityDerivatives;
    //     constexpr IndexType shape_derivatives_size =
    //         TNumNodes * derivatives_type::Shape::TDerivativeDimension;

    //     if (rOutput.size1() != shape_derivatives_size || rOutput.size2() != TElementLocalSize) {
    //         rOutput.resize(shape_derivatives_size, TElementLocalSize, false);
    //     }

    //     rOutput.clear();
    //     AddFluidShapeDerivatives<derivatives_type>(
    //         rOutput, TAdjointElementData::GetIntegrationMethod(), rCurrentProcessInfo);
    // } else {
    //     KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
    //                  << " not supported." << std::endl;
    // }

    // KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
std::string TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "TwoEquationTurbulenceModelAdjointElement #" << Element::Id();
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "TwoEquationTurbulenceModelAdjointElement #" << Element::Id();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::PrintData(std::ostream& rOStream) const
{
    Element::pGetGeometry()->PrintData(rOStream);
}

///@}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidFirstDerivatives(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // KRATOS_TRY

    // Vector gauss_weights;
    // Matrix shape_functions;
    // ShapeFunctionDerivativesArrayType shape_derivatives;
    // this->CalculateGeometryData(gauss_weights, shape_functions,
    //                             shape_derivatives, rIntegrationMethod);

    // typename TDerivativesType::Data element_data(*this, *mpFluidConstitutiveLaw);
    // typename TDerivativesType::U u_derivatives(element_data);
    // typename TDerivativesType::P p_derivatives(element_data);

    // u_derivatives.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
    // p_derivatives.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
    // element_data.Initialize(rCurrentProcessInfo);

    // BoundedVector<double, TElementLocalSize> residual;
    // BoundedMatrix<double, TNumNodes, TDim> dNdXDerivative = ZeroMatrix(TNumNodes, TDim);

    // for (IndexType g = 0; g < gauss_weights.size(); ++g) {
    //     const Vector& N = row(shape_functions, g);
    //     const Matrix& dNdX = shape_derivatives[g];
    //     const double weight = gauss_weights[g];

    //     element_data.CalculateGaussPointData(weight, N, dNdX);
    //     for (IndexType c = 0; c < TNumNodes; ++c) {
    //         const IndexType derivative_block = c * TBlockSize;

    //         for (IndexType k = 0; k < TDim; ++k) {
    //             u_derivatives.CalculateResidualDerivative(
    //                 residual, c, k, weight, N, dNdX, 0.0, 0.0, dNdXDerivative);
    //             row(rLeftHandSideMatrix, derivative_block + k) += residual;
    //         }

    //         p_derivatives.CalculateResidualDerivative(
    //             residual, c, 0, weight, N, dNdX, 0.0, 0.0, dNdXDerivative);
    //         row(rLeftHandSideMatrix, derivative_block + TDim) += residual;
    //     }
    // }

    // KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidSecondDerivatives(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    // KRATOS_TRY

    // Vector gauss_weights;
    // Matrix shape_functions;
    // ShapeFunctionDerivativesArrayType shape_derivatives;
    // this->CalculateGeometryData(gauss_weights, shape_functions,
    //                             shape_derivatives, rIntegrationMethod);

    // TDerivativesType element_data(*this, *mpFluidConstitutiveLaw);

    // element_data.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);

    // for (IndexType g = 0; g < gauss_weights.size(); ++g) {
    //     const Vector& N = row(shape_functions, g);
    //     const Matrix& dNdX = shape_derivatives[g];
    //     const double weight = gauss_weights[g];

    //     element_data.AddResidualDerivativeContributions(rLeftHandSideMatrix,
    //                                                     weight, N, dNdX);
    // }

    // KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidShapeDerivatives(
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    // KRATOS_TRY

    // Vector gauss_weights;
    // Matrix shape_functions;
    // ShapeFunctionDerivativesArrayType shape_derivatives;
    // this->CalculateGeometryData(gauss_weights, shape_functions,
    //                             shape_derivatives, rIntegrationMethod);

    // typename TDerivativesType::Data element_data(*this, *mpFluidConstitutiveLaw);
    // typename TDerivativesType::Shape derivatives(element_data);

    // derivatives.Initialize(rOutput, rCurrentProcessInfo);
    // element_data.Initialize(rCurrentProcessInfo);

    // BoundedVector<double, TElementLocalSize> residual;
    // BoundedMatrix<double, TNumNodes, TDim> dNdXDerivative = ZeroMatrix(TNumNodes, TDim);

    // for (IndexType g = 0; g < gauss_weights.size(); ++g) {
    //     const Vector& N = row(shape_functions, g);
    //     const Matrix& dNdX = shape_derivatives[g];
    //     const double weight = gauss_weights[g];

    //     element_data.CalculateGaussPointData(weight, N, dNdX);

    //     Geometry<Point>::JacobiansType J;
    //     this->GetGeometry().Jacobian(J, rIntegrationMethod);
    //     const auto& DN_De =
    //         this->GetGeometry().ShapeFunctionsLocalGradients(rIntegrationMethod);

    //     GeometricalSensitivityUtility::ShapeFunctionsGradientType dNdX_deriv;
    //     const Matrix& rJ = J[g];
    //     const Matrix& rDN_De = DN_De[g];
    //     const double inv_detJ = 1.0 / MathUtils<double>::DetMat(rJ);
    //     GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

    //     ShapeParameter deriv;
    //     for (deriv.NodeIndex = 0; deriv.NodeIndex < TNumNodes; ++deriv.NodeIndex) {
    //         const IndexType derivative_row =
    //             deriv.NodeIndex * TDerivativesType::Shape::TDerivativeDimension;
    //         for (deriv.Direction = 0;
    //              deriv.Direction < TDerivativesType::Shape::TDerivativeDimension;
    //              ++deriv.Direction) {
    //             double detJ_deriv;
    //             geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, dNdX_deriv);
    //             const double weight_deriv = detJ_deriv * inv_detJ * weight;

    //             derivatives.CalculateResidualDerivative(
    //                 residual, deriv.NodeIndex, deriv.Direction, weight, N, dNdX,
    //                 weight_deriv, detJ_deriv, dNdX_deriv);
    //             row(rOutput, derivative_row + deriv.Direction) += residual;
    //         }
    //     }
    // }

    // KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddTurbulenceFirstDerivatives(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddTurbulenceSecondDerivatives(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddTurbulenceShapeDerivatives(
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void TwoEquationTurbulenceModelAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateGeometryData(
    Vector& rGaussWeights,
    Matrix& rNContainer,
    ShapeFunctionDerivativesArrayType& rDN_DX,
    const GeometryData::IntegrationMethod& rIntegrationMethod) const
{
    const auto& r_geometry = this->GetGeometry();
    const IndexType number_of_gauss_points =
        r_geometry.IntegrationPointsNumber(rIntegrationMethod);

    Vector DetJ;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, rIntegrationMethod);

    if (rNContainer.size1() != number_of_gauss_points || rNContainer.size2() != TNumNodes) {
        rNContainer.resize(number_of_gauss_points, TNumNodes, false);
    }
    rNContainer = r_geometry.ShapeFunctionsValues(rIntegrationMethod);

    const auto& IntegrationPoints = r_geometry.IntegrationPoints(rIntegrationMethod);

    if (rGaussWeights.size() != number_of_gauss_points) {
        rGaussWeights.resize(number_of_gauss_points, false);
    }

    for (IndexType g = 0; g < number_of_gauss_points; ++g) {
        rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
    }
}

// template instantiations
template class TwoEquationTurbulenceModelAdjointElement<2, 3, KEpsilonElementData::QSVMSRFCAdjointElementData<2, 3>>;
template class TwoEquationTurbulenceModelAdjointElement<3, 4, KEpsilonElementData::QSVMSRFCAdjointElementData<3, 4>>;

} // namespace Kratos