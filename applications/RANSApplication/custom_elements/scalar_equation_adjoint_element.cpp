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
#include <sstream>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "utilities/adjoint_extensions.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/indirect_scalar.h"

// Application includes
#include "fluid_dynamics_application_variables.h"
#include "rans_application_variables.h"

// data containers
#include "custom_elements/data_containers/stabilization_validation/circular_convection_rfc_adjoint_element_data.h"
#include "custom_elements/data_containers/stabilization_validation/diffusion_rfc_adjoint_element_data.h"

// Include base h
#include "scalar_equation_adjoint_element.h"

namespace Kratos
{
template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::ThisExtensions(Element* pElement)
    : mpElement{pElement}
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetFirstDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    rVector[0] = MakeIndirectScalar(r_node, (*dofs_list[0]).GetTimeDerivative(), Step);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetSecondDerivativesVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    rVector[0] = MakeIndirectScalar(r_node, (*dofs_list[0]).GetTimeDerivative().GetTimeDerivative(), Step);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetAuxiliaryVector(
    std::size_t NodeId,
    std::vector<IndirectScalar<double>>& rVector,
    std::size_t Step)
{
    auto& r_node = mpElement->GetGeometry()[NodeId];
    rVector.resize(1);

    const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

    rVector[0] = MakeIndirectScalar(r_node, (*dofs_list[0]).GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(), Step);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetFirstDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &RANS_SCALAR_1_ADJOINT_2;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetSecondDerivativesVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &RANS_SCALAR_1_ADJOINT_3;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ThisExtensions::GetAuxiliaryVariables(
    std::vector<VariableData const*>& rVariables) const
{
    rVariables.resize(1);
    rVariables[0] = &RANS_AUX_ADJOINT_SCALAR_1;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ScalarEquationAdjointElement(IndexType NewId)
    : BaseType(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ScalarEquationAdjointElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::ScalarEquationAdjointElement(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::~ScalarEquationAdjointElement()
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>>(
        NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
Element::Pointer ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::Clone(
    IndexType NewId,
    NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>>(
        NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
int ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    TAdjointElementData::Check(*this, rCurrentProcessInfo);
    return 0.0;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::EquationIdVector(
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
            rElementalEquationIdList[local_index++] = r_node.GetDof(*p_variable).EquationId();
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetDofList(
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
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetValuesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize, false);
    }

    const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

    IndexType local_index = 0;
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        for (const auto p_variable : r_variables_list) {
            rValues[local_index++] = r_node.FastGetSolutionStepValue(*p_variable, Step);
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetFirstDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize, false);
    }
    rValues.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::GetSecondDerivativesVector(
    VectorType& rValues,
    int Step) const
{
    if (rValues.size() != TElementLocalSize) {
        rValues.resize(TElementLocalSize, false);
    }

    const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = this->GetGeometry()[i];
        rValues[i] = r_node.FastGetSolutionStepValue(r_variables_list[0]->GetTimeDerivative().GetTimeDerivative(), Step);
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    // If we are restarting, the constitutive law will be already defined
    if (mpConstitutiveLaw == nullptr) {
        const Properties& r_properties = this->GetProperties();

        KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
            << "In initialization of Element " << this->Info()
            << ": No CONSTITUTIVE_LAW defined for property "
            << r_properties.Id() << "." << std::endl;

        // Here we can do down casting because, it should be always a FluidConstitutiveLaw
        mpConstitutiveLaw = r_properties[CONSTITUTIVE_LAW]->Clone();

        const GeometryType& r_geometry = this->GetGeometry();
        const auto& r_shape_functions =
            row(r_geometry.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1), 0);

        // This constitutive law should return nu + nu_t
        mpConstitutiveLaw->InitializeMaterial(r_properties, r_geometry, r_shape_functions);

        // Now we set the constitutive law to be used by the RANS equations
        // because RANS equations only need nu from constitutive law

        const auto rans_cl_name = mpConstitutiveLaw->Info();

        KRATOS_ERROR_IF(rans_cl_name.substr(0, 4) != "Rans")
            << "Incompatible constitutive law is used. Please use constitutive "
            "laws which starts with \"Rans*\" [ Constitutive law "
            "name = "
            << rans_cl_name << " ].\n";

    }

    this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize and initialize output
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLeftHandSide(
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
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "ScalarEquationAdjointElement::"
                    "CalculateRightHandSide method is not implemented.";

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateLocalVelocityContribution(
    MatrixType &rDampMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo &rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != TElementLocalSize)
        rRightHandSideVector.resize(TElementLocalSize, false);

    rRightHandSideVector.clear();

    AddScalarEquationResidualsContributions(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateFirstDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    AddScalarEquationFirstDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
        rLeftHandSideMatrix.size2() != TElementLocalSize) {
        rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
    }

    rLeftHandSideMatrix.clear();

    AddScalarEquationSecondDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    rMassMatrix.resize(0, 0);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "ScalarEquationAdjointElement::"
                    "CalculateDampingMatrix method is not implemented.";

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateSensitivityMatrix(
    const Variable<array_1d<double, 3>>& rSensitivityVariable,
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rSensitivityVariable == SHAPE_SENSITIVITY) {
        if (rOutput.size1() != TCoordLocalSize || rOutput.size2() != TElementLocalSize) {
            rOutput.resize(TCoordLocalSize, TElementLocalSize, false);
        }

        rOutput.clear();
        AddScalarEquationShapeDerivatives(rOutput, rCurrentProcessInfo);
    } else {
        KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                     << " not supported." << std::endl;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::Calculate(
    const Variable<Vector>& rVariable,
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rVariable == PRIMAL_RELAXED_SECOND_DERIVATIVE_VALUES) {

        if (rOutput.size() != TElementLocalSize) {
            rOutput.resize(TElementLocalSize, false);
        }

        const auto& r_variables_list = TAdjointElementData::GetPrimalSecondDerivativeVariablesList();

        for (IndexType i = 0; i < TNumNodes; ++i) {
            const auto& r_node = this->GetGeometry()[i];
            rOutput[i] = r_node.FastGetSolutionStepValue(*r_variables_list[0]);
        }
    } else {
        KRATOS_ERROR << "Unsupported variable requested for Calculate method. "
                        "[ rVariable.Name() = "
                     << rVariable.Name() << " ].\n";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::Calculate(
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

        AddScalarEquationFirstDerivatives(rOutput, rCurrentProcessInfo, 0.0);

    } else {
        KRATOS_ERROR << "Unsupported variable requested for Calculate method. "
                        "[ rVariable.Name() = "
                     << rVariable.Name() << " ].\n";
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
std::string ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "ScalarEquationAdjointElement #" << Element::Id();
    return buffer.str();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ScalarEquationAdjointElement #" << Element::Id();
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::PrintData(std::ostream& rOStream) const
{
    Element::pGetGeometry()->PrintData(rOStream);
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddScalarEquationResidualsContributions(
    Vector& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    const auto& integration_method = ScalarEquationData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using ScalarEquationPrimal= typename ScalarEquationData::Primal;

    // create data holders for turbulence equations
    typename ScalarEquationPrimal::Data eq_data(*this, rCurrentProcessInfo);

    // create equation residual data holders
    typename ScalarEquationPrimal::ResidualsContributions eq_residuals(eq_data);

    VectorN residuals = ZeroVector(TNumNodes);

    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        eq_data.CalculateGaussPointData(W, N, dNdX);

        eq_residuals.AddGaussPointResidualsContributions(residuals, W, N, dNdX);
    }

    eq_data.CalculateDataAfterGaussPointPointLoop();

    eq_residuals.AddResidualsContributionsAfterGaussPointLoop(residuals);

    noalias(rOutput) += residuals;
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddScalarEquationFirstDerivatives(
    MatrixType& rOutput,
    const ProcessInfo& rCurrentProcessInfo,
    const double MassTermsDerivativesWeight)
{
    KRATOS_TRY

    const auto& integration_method = ScalarEquationData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using ScalarEquationDerivatives = typename ScalarEquationData::StateDerivatives::FirstDerivatives;

    // create data holders for turbulence equations
    typename ScalarEquationDerivatives::Data eq_data(*this, rCurrentProcessInfo);

    // create equation derivative data holders
    typename ScalarEquationDerivatives::ScalarVariable eq_derivative(eq_data);

    MatrixND dNdX_derivative = ZeroMatrix(TNumNodes, TDim);
    VectorN residual_derivatives;

    IndexType row_index = 0;
    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        eq_data.CalculateGaussPointData(W, N, dNdX);

        row_index = 0;
        for (IndexType c = 0; c < TNumNodes; ++c) {
            // add derivatives w.r.t. scalar variable
            eq_derivative.CalculateGaussPointResidualsDerivativeContributions(residual_derivatives, c, 0, W, N, dNdX, 0, 0, dNdX_derivative, MassTermsDerivativesWeight);
            AssembleSubVectorToMatrix(rOutput, row_index, 0, residual_derivatives);
            ++row_index;
        }
    }

    eq_data.CalculateDataAfterGaussPointPointLoop();

    // finalize derivative data holders
    row_index = 0;
    for (IndexType c = 0; c < TNumNodes; ++c) {
        // initializing derivative row w.r.t. turbulence variable 1
        eq_derivative.CalculateResidualsDerivativeContributionsAfterGaussPointPointLoop(residual_derivatives, c, 0);
        AssembleSubVectorToMatrix(rOutput, row_index, 0, residual_derivatives);
        ++row_index;
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddScalarEquationSecondDerivatives(
    MatrixType& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = ScalarEquationData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using ScalarEquationDerivatives = typename ScalarEquationData::StateDerivatives::SecondDerivatives;

    // create data holders for turbulence equations
    typename ScalarEquationDerivatives::Data eq_data(*this, rCurrentProcessInfo);

    // create equation derivative data holders
    typename ScalarEquationDerivatives::ScalarVariableRate eq_derivative(eq_data);

    MatrixND dNdX_derivative = ZeroMatrix(TNumNodes, TDim);
    VectorN residual_derivatives;

    IndexType row_index = 0;
    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        eq_data.CalculateGaussPointData(W, N, dNdX);

        row_index = 0;
        for (IndexType c = 0; c < TNumNodes; ++c) {
            // add derivatives w.r.t. turbulence variable 1
            eq_derivative.CalculateGaussPointResidualsDerivativeContributions(residual_derivatives, c, W, N, dNdX);
            AssembleSubVectorToMatrix(rOutput, row_index++, 0, residual_derivatives);
        }
    }

    eq_data.CalculateDataAfterGaussPointPointLoop();

    // finalize derivative data holders
    row_index = 0;
    for (IndexType c = 0; c < TNumNodes; ++c) {
        // initializing derivative row w.r.t. turbulence variable 1
        eq_derivative.CalculateResidualsDerivativeContributionsAfterGaussPointPointLoop(residual_derivatives, c);
        AssembleSubVectorToMatrix(rOutput, row_index++, 0, residual_derivatives);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::AddScalarEquationShapeDerivatives(
    Matrix& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& integration_method = ScalarEquationData::TResidualsDerivatives::GetIntegrationMethod();

    Vector Ws;
    Matrix Ns;
    ShapeFunctionDerivativesArrayType dNdXs;
    this->CalculateGeometryData(Ws, Ns, dNdXs, integration_method);

    using ScalarEquationDerivatives = typename ScalarEquationData::SensitivityDerivatives;

    // create data holders for turbulence equations
    typename ScalarEquationDerivatives::Data eq_data(*this, rCurrentProcessInfo);

    // create equation derivative data holders
    typename ScalarEquationDerivatives::Shape eq_derivative(eq_data);

    MatrixND dNdX_derivative = ZeroMatrix(TNumNodes, TDim);
    VectorN residual_derivatives;

    IndexType row_index = 0;
    for (IndexType g = 0; g < Ws.size(); ++g) {
        const double W = Ws[g];
        const Vector& N = row(Ns, g);
        const Matrix& dNdX = dNdXs[g];

        eq_data.CalculateGaussPointData(W, N, dNdX);

        Geometry<Point>::JacobiansType J;
        this->GetGeometry().Jacobian(J, integration_method);
        const auto& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(integration_method);

        GeometricalSensitivityUtility::ShapeFunctionsGradientType dNdX_deriv;
        const Matrix& rJ = J[g];
        const Matrix& rDN_De = DN_De[g];
        const double inv_detJ = 1.0 / MathUtils<double>::Det(rJ);
        GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

        row_index = 0;
        ShapeParameter deriv;
        for (deriv.NodeIndex = 0; deriv.NodeIndex < TNumNodes; ++deriv.NodeIndex) {
            for (deriv.Direction = 0; deriv.Direction < TDim; ++deriv.Direction) {
                double detJ_deriv;
                geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, dNdX_deriv);
                const double weight_deriv = detJ_deriv * inv_detJ * W;

                eq_derivative.CalculateGaussPointResidualsDerivativeContributions(residual_derivatives, deriv.NodeIndex, deriv.Direction, W, N, dNdX, weight_deriv, detJ_deriv, dNdX_deriv);
                AssembleSubVectorToMatrix(rOutput, row_index, 0, residual_derivatives);
                ++row_index;
            }
        }
    }

    eq_data.CalculateDataAfterGaussPointPointLoop();

    // finalize derivative data holders
    row_index = 0;
    for (IndexType c = 0; c < TNumNodes; ++c) {
        for (IndexType k = 0; k < TDim; ++k) {
            eq_derivative.CalculateResidualsDerivativeContributionsAfterGaussPointPointLoop(residual_derivatives, c, k);
            AssembleSubVectorToMatrix(rOutput, row_index, 0, residual_derivatives);
            ++row_index;
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
void ScalarEquationAdjointElement<TDim, TNumNodes, TAdjointElementData>::CalculateGeometryData(
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
// k-epsilon
template class ScalarEquationAdjointElement<2, 3, StabilizationValidationElementData::CircularConvectionRFCAdjointElementData>;
template class ScalarEquationAdjointElement<2, 3, StabilizationValidationElementData::DiffusionRFCAdjointElementData>;

} // namespace Kratos