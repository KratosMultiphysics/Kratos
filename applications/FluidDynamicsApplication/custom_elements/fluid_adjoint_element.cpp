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
#include "custom_elements/data_containers/fluid_adjoint_derivatives.h"
#include "fluid_dynamics_application_variables.h"

#include "custom_elements/data_containers/qs_vms/qs_vms_adjoint_element_data.h"

// Include base h
#include "fluid_adjoint_element.h"

namespace Kratos
{
namespace FluidAdjointElementHelperUtilities
{
template<class TDerivativeDataContainerType, class TCombinedElementDataContainerType>
void AddFluidShapeDerivativesImpl(
    TDerivativeDataContainerType& rDerivativeDataHolderType,
    TCombinedElementDataContainerType& rCombinedElementDataContainerType,
    GeometricalSensitivityUtility& rGeometricalSensitivityUtility,
    ShapeParameter& rShapeParameter,
    GeometricalSensitivityUtility::ShapeFunctionsGradientType& rdNdXDerivative,
    const double InvDetJ,
    const double W,
    const Vector& N,
    const Matrix& dNdX)
{
    double detJ_derivative;
    rGeometricalSensitivityUtility.CalculateSensitivity(rShapeParameter, detJ_derivative, rdNdXDerivative);
    const double W_derivative = detJ_derivative * InvDetJ * W;

    rDerivativeDataHolderType.CalculateGaussPointResidualsDerivativeContributions(
        rDerivativeDataHolderType.GetSubVector(),
        rDerivativeDataHolderType.GetElementDataContainer(rCombinedElementDataContainerType),
        rShapeParameter.NodeIndex, W, N, dNdX, W_derivative, detJ_derivative, rdNdXDerivative);

    rShapeParameter.Direction++;
}
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::ThisExtensions(Element* pElement)
    : mpElement{pElement}
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(TDim + 1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative(), Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(TDim + 1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(
            r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative(), Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(TDim + 1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    for (unsigned int i = 0; i < TDim; ++i) {
        rVector[i] = MakeIndirectScalar(
            r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(),
            Step);
    }

    rVector[TDim] = IndirectScalar<double>{}; // pressure
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::FluidAdjointElement(IndexType NewId)
    : BaseType(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::FluidAdjointElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::FluidAdjointElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::~FluidAdjointElement()
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FluidAdjointElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FluidAdjointElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Clone(
    IndexType NewId,
    NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<FluidAdjointElement>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
int FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    TAdjointElementData::Check(*this, rCurrentProcessInfo);
    return BaseType::Check(rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::EquationIdVector(
    EquationIdVectorType& rElementalEquationIdList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rElementalEquationIdList.size() != TElementLocalSize) {
        rElementalEquationIdList.resize(TElementLocalSize, false);
    }

    const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        for (const auto p_variable : r_variables_list) {
            rElementalEquationIdList[local_index++] =
                r_node.GetDof(*p_variable).EquationId();
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != TElementLocalSize) {
        rElementalDofList.resize(TElementLocalSize);
    }

    const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        for (const auto p_variable : r_variables_list) {
            rElementalDofList[local_index++] = r_node.pGetDof(*p_variable);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetValuesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize, false);
    }

    const auto& r_geometry = this->GetGeometry();
    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        const auto& r_velocity =
            r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, Step);
        for (IndexType d = 0; d < TDim; ++d) {
            rValues[local_index++] = r_velocity[d];
        }
        rValues[local_index++] =
            r_node.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetFirstDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize, false);
    }
    rValues.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetSecondDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize);
    }

    const auto& r_geometry = this->GetGeometry();
    IndexType local_index = 0;
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const auto& r_acceleration =
            r_geometry[i_node].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, Step);
        for (IndexType d = 0; d < TDim; ++d) {
            rValues[local_index++] = r_acceleration[d];
        }
        rValues[local_index++] = 0.0; // pressure dof
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // If we are restarting, the constitutive law will be already defined
    if (mpConstitutiveLaw == nullptr) {
        const Properties& r_properties = this->GetProperties();

        KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
            << "In initialization of Element " << this->Info()
            << ": No CONSTITUTIVE_LAW defined for property "
            << r_properties.Id() << "." << std::endl;

        mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();

        const GeometryType& r_geometry = this->GetGeometry();
        const auto& r_shape_functions =
            r_geometry.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);

        mpConstitutiveLaw->InitializeMaterial(r_properties, r_geometry,
                                              row(r_shape_functions, 0));
    }

    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize and initialize output
    if (rLeftHandSideMatrix.size1() != TElementLocalSize || rLeftHandSideMatrix.size2() != TElementLocalSize)
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

    rLeftHandSideMatrix.clear();

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize || rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "FluidAdjointElement::"
                    "CalculateRightHandSide method is not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLocalVelocityContribution(
    MatrixType &rDampMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rDampMatrix.size1() != TElementLocalSize || rDampMatrix.size2() != TElementLocalSize) {
        rDampMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rDampMatrix.clear();

    if (rRightHandSideVector.size() != TElementLocalSize)
        rRightHandSideVector.resize(TElementLocalSize, false);

    rRightHandSideVector.clear();

    AddFluidResidualsContributions(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    AddFluidFirstDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    AddFluidSecondDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    rMassMatrix.resize(0, 0);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "FluidAdjointElement::"
                    "CalculateDampingMatrix method is not implemented.";

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY) {
        if (rOutput.size1() != TCoordsLocalSize || rOutput.size2() != TElementLocalSize) {
            rOutput.resize(TCoordsLocalSize, TElementLocalSize, false);
        }

        rOutput.clear();
        AddFluidShapeDerivatives(rOutput, rCurrentProcessInfo);
    } else {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Calculate(
    const Variable<Vector>& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == PRIMAL_RELAXED_SECOND_DERIVATIVE_VALUES) {
        if (rOutput.size() != TElementLocalSize) {
            rOutput.resize(TElementLocalSize);
        }

        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const auto& r_node = this->GetGeometry()[i_node];
            rOutput[local_index++] = r_node.GetValue(RELAXED_ACCELERATION_X);
            rOutput[local_index++] = r_node.GetValue(RELAXED_ACCELERATION_Y);

            if (TDim == 3)
                rOutput[local_index++] = r_node.GetValue(RELAXED_ACCELERATION_Z);

            rOutput[local_index++] = 0.0;
        }
    } else {
        KRATOS_ERROR << "Unsupported variable requested Calculate method. [ rVariable.Name() = " << rVariable.Name() << " ].\n";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Calculate(
    const Variable<Matrix>& rVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == PRIMAL_STEADY_RESIDUAL_FIRST_DERIVATIVES) {

        if (rOutput.size1() != TElementLocalSize || rOutput.size2() != TElementLocalSize) {
            rOutput.resize(TElementLocalSize, TElementLocalSize, false);
        }

        rOutput.clear();

        AddFluidFirstDerivatives(rOutput, rCurrentProcessInfo, 0.0);
    } else {
        KRATOS_ERROR << "Unsupported variable requested for Calculate method. "
                        "[ rVariable.Name() = "
                     << rVariable.Name() << " ].\n";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
std::string FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "FluidAdjointElement #" << Element::Id();
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "FluidAdjointElement #" << Element::Id();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::PrintData(std::ostream& rOStream) const
{
    Element::pGetGeometry()->PrintData(rOStream);
}

///@}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidResidualsContributions(
    VectorType& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = TAdjointElementData::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using residual_data_container_type = typename TAdjointElementData::Residual;

    // create data holders tuple
    typename residual_data_container_type::CombinedElementDataContainerType combined_element_data;
    // create derivatives tuple
    typename residual_data_container_type::CombinedCalculationContainersType combined_residuals;

    // initialize element data holders
    std::apply([&](auto&... args) {((args.Initialize(*this, *mpConstitutiveLaw, rCurrentProcessInfo)), ...);}, combined_element_data);

    // initialize residual data holders
    std::apply([](auto&... args) {((args.GetSubVector().clear()), ...);}, combined_residuals);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];
        const double W = Ws[g];

        // calculate gauss point data in each data holder
        std::apply([W, &N, &dNdX](auto&&... args) {((args.CalculateGaussPointData(W, N, dNdX)), ...);}, combined_element_data);

        // calculate gauss point residual contributions
        std::apply([W, &N, &dNdX, &combined_element_data](auto&&... args) {((
            args.AddGaussPointResidualsContributions(
                args.GetSubVector(),
                args.GetElementDataContainer(combined_element_data),
                W, N, dNdX)), ...);}, combined_residuals);
    }

    // assemble vector to matrix
    std::apply([&rOutput](auto&&... args) {
        ((args.template AssembleSubVectorToVector<TBlockSize>(rOutput)), ...);}, combined_residuals);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidFirstDerivatives(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo,
    const double MassTermsDerivativesWeight)
{
    KRATOS_TRY

    const auto& integration_method = TAdjointElementData::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using derivative_data_container_type = typename TAdjointElementData::ResidualStateVariableFirstDerivatives;

    // create data holders tuple
    typename derivative_data_container_type::CombinedElementDataContainerType combined_element_data;
    // create derivatives tuple
    typename derivative_data_container_type::CombinedCalculationContainersType combined_derivatives;

    BoundedMatrix<double, TNumNodes, TDim> dNdXDerivative = ZeroMatrix(TNumNodes, TDim);

    // initialize element data holders
    std::apply([&](auto&... args) {((args.Initialize(*this, *mpConstitutiveLaw, rCurrentProcessInfo)), ...);}, combined_element_data);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        // calculate gauss point data in each data holder
        std::apply([W, &N, &dNdX](auto&&... args) {((args.CalculateGaussPointData(W, N, dNdX)), ...);}, combined_element_data);

        for (IndexType c = 0; c < TNumNodes; ++c) {
            // calculate derivatives
            std::apply([c, W, &N, &dNdX, &dNdXDerivative, MassTermsDerivativesWeight, &combined_element_data](auto&&... args) {
                ((args.CalculateGaussPointResidualsDerivativeContributions(
                    args.GetSubVector(),
                    args.GetElementDataContainer(combined_element_data),
                    c, W, N, dNdX, 0.0, 0.0, dNdXDerivative, MassTermsDerivativesWeight)), ...);}, combined_derivatives);

            // assemble vector to matrix
            std::apply([&rLeftHandSideMatrix, c](auto&&... args) {
                ((args.template AssembleSubVectorToMatrix<TBlockSize>(rLeftHandSideMatrix, c)), ...);}, combined_derivatives);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidSecondDerivatives(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = TAdjointElementData::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using derivative_data_container_type = typename TAdjointElementData::ResidualStateVariableSecondDerivatives;

    // create data holders tuple
    typename derivative_data_container_type::CombinedElementDataContainerType combined_element_data;
    // create derivatives tuple
    typename derivative_data_container_type::CombinedCalculationContainersType combined_derivatives;

    std::apply([&](auto&... args) {((args.Initialize(*this, *mpConstitutiveLaw, rCurrentProcessInfo)), ...);}, combined_element_data);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        // calculate gauss point data in each data holder
        std::apply([W, &N, &dNdX](auto&&... args) {((args.CalculateGaussPointData(W, N, dNdX)), ...);}, combined_element_data);

        for (IndexType c = 0; c < TNumNodes; ++c) {
            // calculate derivatives
            std::apply([c, W, &N, &dNdX, &combined_element_data](auto&&... args) {
                ((args.CalculateGaussPointResidualsDerivativeContributions(
                    args.GetSubVector(),
                    args.GetElementDataContainer(combined_element_data),
                    c, W, N, dNdX)), ...);}, combined_derivatives);

            // assemble vector to matrix
            std::apply([&rLeftHandSideMatrix, c](auto&&... args) {
                ((args.template AssembleSubVectorToMatrix<TBlockSize>(rLeftHandSideMatrix, c)), ...);}, combined_derivatives);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddFluidShapeDerivatives(
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = TAdjointElementData::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using derivative_data_container_type = typename TAdjointElementData::ResidualShapeDerivatives;

    // create data holders tuple
    typename derivative_data_container_type::CombinedElementDataContainerType combined_element_data;
    // create derivatives tuple
    typename derivative_data_container_type::CombinedCalculationContainersType combined_derivatives;

    std::apply([&](auto&... args) {((args.Initialize(*this, *mpConstitutiveLaw, rCurrentProcessInfo)), ...);}, combined_element_data);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];
        const double W = Ws[g];

        // calculate gauss point data in each data holder
        std::apply([W, &N, &dNdX](auto&&... args) {((args.CalculateGaussPointData(W, N, dNdX)), ...);}, combined_element_data);

        Geometry<Point>::JacobiansType J;
        this->GetGeometry().Jacobian(J, integration_method);
        const auto& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(integration_method);

        GeometricalSensitivityUtility::ShapeFunctionsGradientType dNdX_derivative;
        const Matrix& rJ = J[g];
        const Matrix& rDN_De = DN_De[g];
        const double inv_detJ = 1.0 / MathUtils<double>::Det(rJ);
        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        ShapeParameter deriv;
        for (deriv.NodeIndex = 0; deriv.NodeIndex < TNumNodes; ++deriv.NodeIndex) {
            deriv.Direction = 0;

            // calculate derivatives
            std::apply([&combined_element_data, &geom_sensitivity, &deriv, &dNdX_derivative, inv_detJ, W, &N, &dNdX](auto&&... args) {
                ((FluidAdjointElementHelperUtilities::AddFluidShapeDerivativesImpl(
                    args,
                    combined_element_data,
                    geom_sensitivity,
                    deriv,
                    dNdX_derivative,
                    inv_detJ, W, N, dNdX)), ...);}, combined_derivatives);

            // assemble vector to matrix
            std::apply([&rOutput, &deriv](auto&&... args) {
                ((args.template AssembleSubVectorToMatrix<TDim, TBlockSize>(rOutput, deriv.NodeIndex)), ...);}, combined_derivatives);
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void FluidAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateGeometryData(
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
template class FluidAdjointElement<2, 3, QSVMSAdjointElementData<2, 3>>;
template class FluidAdjointElement<2, 4, QSVMSAdjointElementData<2, 4>>;
template class FluidAdjointElement<3, 4, QSVMSAdjointElementData<3, 4>>;
template class FluidAdjointElement<3, 8, QSVMSAdjointElementData<3, 8>>;

} // namespace Kratos