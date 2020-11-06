//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_EVM_MONOLITHIC_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED)
#define KRATOS_EVM_MONOLITHIC_K_EPSILON_VMS_ADJOINT_ELEMENT_H_INCLUDED

// System includes
#include <algorithm>
#include <iterator>

// External includes

// Project includes
// #include "custom_elements/rans_evm_vms_adjoint_element.h"
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
// #include "rans_application_variables.h"
#include "utilities/adjoint_extensions.h"

// #include "custom_elements/evm_k_epsilon/evm_epsilon_adjoint_element.h"
// #include "custom_elements/evm_k_epsilon/evm_k_adjoint_element.h"
// #include "custom_elements/evm_k_epsilon/evm_k_epsilon_vms_adjoint_element.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointTurbulenceModelDataType>
class TwoEquationTurbulenceModelAdjointElement : public Element
{
    class ThisExtensions : public AdjointExtensions
    {
        Element* mpElement;

    public:
        explicit ThisExtensions(Element* pElement) : mpElement{pElement}
        {
        }

        void GetFirstDerivativesVector(std::size_t NodeId,
                                       std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step) override
        {
            // auto& r_node = mpElement->GetGeometry()[NodeId];
            // rVector.resize(TDim + 3);
            // std::size_t index = 0;
            // rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_X, Step);
            // rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Y, Step);
            // if (TDim == 3)
            // {
            //     rVector[index++] =
            //         MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_2_Z, Step);
            // }
            // rVector[index++] = IndirectScalar<double>{}; // pressure
            // rVector[index++] = MakeIndirectScalar(r_node, RANS_SCALAR_1_ADJOINT_2, Step);
            // rVector[index] = MakeIndirectScalar(r_node, RANS_SCALAR_2_ADJOINT_2, Step);
        }

        void GetSecondDerivativesVector(std::size_t NodeId,
                                        std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step) override
        {
            // auto& r_node = mpElement->GetGeometry()[NodeId];
            // rVector.resize(TDim + 3);
            // std::size_t index = 0;
            // rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_X, Step);
            // rVector[index++] = MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Y, Step);
            // if (TDim == 3)
            // {
            //     rVector[index++] =
            //         MakeIndirectScalar(r_node, ADJOINT_FLUID_VECTOR_3_Z, Step);
            // }
            // rVector[index++] = IndirectScalar<double>{}; // pressure
            // rVector[index++] = MakeIndirectScalar(r_node, RANS_SCALAR_1_ADJOINT_3, Step);
            // rVector[index] = MakeIndirectScalar(r_node, RANS_SCALAR_2_ADJOINT_3, Step);
        }

        void GetAuxiliaryVector(std::size_t NodeId,
                                std::vector<IndirectScalar<double>>& rVector,
                                std::size_t Step) override
        {
            // auto& r_node = mpElement->GetGeometry()[NodeId];
            // rVector.resize(TDim + 3);
            // std::size_t index = 0;
            // rVector[index++] =
            //     MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_X, Step);
            // rVector[index++] =
            //     MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Y, Step);
            // if (TDim == 3)
            // {
            //     rVector[index++] =
            //         MakeIndirectScalar(r_node, AUX_ADJOINT_FLUID_VECTOR_1_Z, Step);
            // }
            // rVector[index++] = IndirectScalar<double>{}; // pressure
            // rVector[index++] = MakeIndirectScalar(r_node, RANS_AUX_ADJOINT_SCALAR_1, Step);
            // rVector[index] = MakeIndirectScalar(r_node, RANS_AUX_ADJOINT_SCALAR_2, Step);
        }

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            // rVariables.resize(3);
            // rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
            // rVariables[1] = &RANS_SCALAR_1_ADJOINT_2;
            // rVariables[2] = &RANS_SCALAR_2_ADJOINT_2;
        }

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            // rVariables.resize(3);
            // rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
            // rVariables[1] = &RANS_SCALAR_1_ADJOINT_3;
            // rVariables[2] = &RANS_SCALAR_2_ADJOINT_3;
        }

        void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override
        {
            // rVariables.resize(3);
            // rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
            // rVariables[1] = &RANS_AUX_ADJOINT_SCALAR_1;
            // rVariables[2] = &RANS_AUX_ADJOINT_SCALAR_2;
        }
    };

public:
    ///@name Type Definitions
    ///@{

    // defining the base type
    typedef Element BaseType;
    // // defining the base adjoint base fluid element type
    // typedef EvmKEpsilonVMSAdjointElement<TDim, TNumNodes> AdjointFluidElement;
    // // defining the k element type
    // typedef EvmKAdjointElement<TDim, TNumNodes> AdjointKElement;
    // // defining the epsilon element type
    // typedef EvmEpsilonAdjointElement<TDim, TNumNodes> AdjointEpsilonElement;

    // constexpr static unsigned int TFluidBlockSize = (TDim + 1);

    // constexpr static unsigned int TFluidLocalSize = TFluidBlockSize * TNumNodes;

    // constexpr static unsigned int TKBlockSize = 1;

    // constexpr static unsigned int TKLocalSize = TKBlockSize * TNumNodes;

    // constexpr static unsigned int TEpsilonBlockSize = 1;

    // constexpr static unsigned int TEpsilonLocalSize = TEpsilonBlockSize * TNumNodes;

    // constexpr static unsigned int TElementBlockSize =
    //     (TFluidBlockSize + TKBlockSize + TEpsilonBlockSize);

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    // variable definitions
    typedef std::size_t IndexType;

    typedef Element::NodeType NodeType;

    typedef Element::NodesArrayType NodesArrayType;

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::VectorType VectorType;

    typedef Element::MatrixType MatrixType;

    constexpr static unsigned int TPrimalFluidEqOffset = 0;
    constexpr static unsigned int TPrimalTurbulenceEquation1EqOffset = TDim + 1;
    constexpr static unsigned int TPrimalTurbulenceEquation2EqOffset = TDim + 2;
    constexpr static unsigned int TElementBlockSize = TDim + 3;
    constexpr static unsigned int TElementLocalSize = TElementBlockSize * TNumNodes;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using TurbulenceEquation1Data = typename TAdjointTurbulenceModelDataType::TurbulenceEquation1;

    using TurbulenceEquation2Data = typename TAdjointTurbulenceModelDataType::TurbulenceEquation2;


    ///@name Pointer Definitions
    /// Pointer definition of TwoEquationTurbulenceModelAdjointElement
    KRATOS_CLASS_POINTER_DEFINITION(TwoEquationTurbulenceModelAdjointElement);

    ///@}

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    TwoEquationTurbulenceModelAdjointElement(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using Geometry
     */
    TwoEquationTurbulenceModelAdjointElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    TwoEquationTurbulenceModelAdjointElement(IndexType NewId,
                                           GeometryType::Pointer pGeometry,
                                           PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Destructor
     */
    ~TwoEquationTurbulenceModelAdjointElement() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void Initialize() override
    {
        this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointElement>(
            NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointElement>(
            NewId, pGeom, pProperties);
        KRATOS_CATCH("");
    }

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<TwoEquationTurbulenceModelAdjointElement>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        const auto& geometry = this->GetGeometry();
        // TAdjointTurbulenceModelDataType::FluidEquation::Check(geometry, rCurrentProcessInfo);
        TurbulenceEquation1Data::StateDerivatives::Check(geometry, rCurrentProcessInfo);
        TurbulenceEquation2Data::StateDerivatives::Check(geometry, rCurrentProcessInfo);
        return 0;
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(
        EquationIdVectorType& rElementalEquationIdList,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        if (rElementalEquationIdList.size() != TElementLocalSize) {
            rElementalEquationIdList.resize(TElementLocalSize, false);
        }

        IndexType local_index = 0;
        const auto& set_equation_ids =
            [&](const NodeType& rNode, const std::vector<const Variable<double>*>& rVariables) {
                for (const auto p_variable : rVariables) {
                    rElementalEquationIdList[local_index++] = rNode.GetDof(*p_variable).EquationId();
                }
            };

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            const auto& r_node = this->GetGeometry()[i];
            // set_equation_ids(r_node, TAdjointTurbulenceModelDataType::FluidEquation::GetDofVariablesList());
            set_equation_ids(r_node, TurbulenceEquation1Data::StateDerivatives::GetDofVariablesList());
            set_equation_ids(r_node, TurbulenceEquation2Data::StateDerivatives::GetDofVariablesList());
        }
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo) const override
    {
        if (rElementalDofList.size() != TElementLocalSize) {
            rElementalDofList.resize(TElementLocalSize);
        }

        IndexType local_index = 0;
        const auto& set_dof_variables =
            [&](const NodeType& rNode, const std::vector<const Variable<double>*>& rVariables) {
                for (const auto p_variable : rVariables) {
                    rElementalDofList[local_index++] = rNode.pGetDof(*p_variable);
                }
            };

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            const auto& r_node = this->GetGeometry()[i];
            // set_dof_variables(r_node, TAdjointTurbulenceModelDataType::FluidEquation::GetDofVariablesList());
            set_dof_variables(r_node, TurbulenceEquation1Data::StateDerivatives::GetDofVariablesList());
            set_dof_variables(r_node, TurbulenceEquation2Data::StateDerivatives::GetDofVariablesList());
        }
    }

    /// Returns the adjoint values stored in this element's nodes.
    void GetValuesVector(
        VectorType& rValues,
        int Step = 0) const override
    {
        GetFirstDerivativesVector(rValues, Step);
    }

    /// Returns the adjoint velocity values stored in this element's nodes.
    void GetFirstDerivativesVector(
        VectorType& rValues,
        int Step = 0) const override
    {
        if (rValues.size() != TElementLocalSize) {
            rValues.resize(TElementLocalSize);
        }

        IndexType local_index = 0;
        const auto& set_first_derivative_values =
            [&](const NodeType& rNode, const std::vector<const Variable<double>*>& rVariables) {
                for (const auto p_variable : rVariables) {
                    rValues[local_index++] = rNode.FastGetSolutionStepValue(p_variable->GetTimeDerivative());
                }
            };

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            const auto& r_node = this->GetGeometry()[i];
            // set_first_derivative_values(r_node, TAdjointTurbulenceModelDataType::FluidEquation::GetDofVariablesList());
            set_first_derivative_values(r_node, TurbulenceEquation1Data::StateDerivatives::GetDofVariablesList());
            set_first_derivative_values(r_node, TurbulenceEquation2Data::StateDerivatives::GetDofVariablesList());
        }
    }

    void GetSecondDerivativesVector(
        VectorType& rValues,
        int Step) const override
    {
        if (rValues.size() != TElementLocalSize) {
            rValues.resize(TElementLocalSize);
        }

        IndexType local_index = 0;
        const auto& set_second_derivative_values =
            [&](const NodeType& rNode, const std::vector<const Variable<double>*>& rVariables) {
                for (const auto p_variable : rVariables) {
                    rValues[local_index++] = rNode.FastGetSolutionStepValue(p_variable->GetTimeDerivative());
                }
            };

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            const auto& r_node = this->GetGeometry()[i];
            // set_second_derivative_values(r_node, TAdjointTurbulenceModelDataType::FluidEquation::GetDofVariablesList());
            set_second_derivative_values(r_node, TurbulenceEquation1Data::StateDerivatives::GetDofVariablesList());
            set_second_derivative_values(r_node, TurbulenceEquation2Data::StateDerivatives::GetDofVariablesList());
        }
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "TwoEquationTurbulenceModelAdjointElement::"
                        "CalculateLocalSystem method is not implemented.";

        KRATOS_CATCH("");
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize)
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

        rLeftHandSideMatrix.clear();
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "TwoEquationTurbulenceModelAdjointElement::"
                        "CalculateRightHandSide method is not implemented.";

        KRATOS_CATCH("");
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize)
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

        rLeftHandSideMatrix.clear();

        using turbulence_equation_1_derivatives_data = typename TurbulenceEquation1Data::StateDerivatives::FirstDerivatives;
        using turbulence_equation_2_derivatives_data = typename TurbulenceEquation2Data::StateDerivatives::FirstDerivatives;

        typename turbulence_equation_1_derivatives_data::Data equation_1_data(this->GetGeometry());
        typename turbulence_equation_1_derivatives_data::VelocityPressure    turbulence_equation_1_derivative_variable_0(equation_1_data);
        typename turbulence_equation_1_derivatives_data::TurbulenceVariable1 turbulence_equation_1_derivative_variable_1(equation_1_data);
        typename turbulence_equation_1_derivatives_data::TurbulenceVariable2 turbulence_equation_1_derivative_variable_2(equation_1_data);

        typename turbulence_equation_2_derivatives_data::Data equation_2_data(this->GetGeometry());
        typename turbulence_equation_2_derivatives_data::VelocityPressure    turbulence_equation_2_derivative_variable_0(equation_2_data);
        typename turbulence_equation_2_derivatives_data::TurbulenceVariable1 turbulence_equation_2_derivative_variable_1(equation_2_data);
        typename turbulence_equation_2_derivatives_data::TurbulenceVariable2 turbulence_equation_2_derivative_variable_2(equation_2_data);

        turbulence_equation_1_derivative_variable_0.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_1_derivative_variable_1.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_1_derivative_variable_2.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);

        turbulence_equation_2_derivative_variable_0.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_2_derivative_variable_1.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_2_derivative_variable_2.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);

        // fluid_data.Initialize(rCurrentProcessInfo);
        equation_1_data.Initialize(rCurrentProcessInfo);
        equation_2_data.Initialize(rCurrentProcessInfo);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        RansCalculationUtilities::CalculateGeometryData(this->GetGeometry(), TurbulenceEquation1Data::GetIntegrationMethod(), gauss_weights, shape_functions, shape_derivatives);

        for (IndexType g = 0; g < gauss_weights.size(); ++g) {
            const Vector& N = row(shape_functions, g);
            const Matrix& dNdX = shape_derivatives[g];
            const double weight = gauss_weights[g];

            // turbulence equation 1 derivatives
            equation_1_data.CalculateGaussPointData(weight, N, dNdX);
            turbulence_equation_1_derivative_variable_0.CalculateResidualDerivatives(rLeftHandSideMatrix, weight, N, dNdX);
            turbulence_equation_1_derivative_variable_1.CalculateResidualDerivatives(rLeftHandSideMatrix, weight, N, dNdX);
            turbulence_equation_1_derivative_variable_2.CalculateResidualDerivatives(rLeftHandSideMatrix, weight, N, dNdX);

            // turbulence equation 2 derivatives
            equation_2_data.CalculateGaussPointData(weight, N, dNdX);
            turbulence_equation_2_derivative_variable_0.CalculateResidualDerivatives(rLeftHandSideMatrix, weight, N, dNdX);
            turbulence_equation_2_derivative_variable_1.CalculateResidualDerivatives(rLeftHandSideMatrix, weight, N, dNdX);
            turbulence_equation_2_derivative_variable_2.CalculateResidualDerivatives(rLeftHandSideMatrix, weight, N, dNdX);
        }

        // fluid_data.Finalize(rCurrentProcessInfo);
        equation_1_data.Finalize(rCurrentProcessInfo);
        equation_2_data.Finalize(rCurrentProcessInfo);

        turbulence_equation_1_derivative_variable_0.Finalize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_1_derivative_variable_1.Finalize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_1_derivative_variable_2.Finalize(rLeftHandSideMatrix, rCurrentProcessInfo);

        turbulence_equation_2_derivative_variable_0.Finalize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_2_derivative_variable_1.Finalize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_2_derivative_variable_2.Finalize(rLeftHandSideMatrix, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize)
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);

        rLeftHandSideMatrix.clear();

        using turbulence_equation_1_derivatives_data = typename TurbulenceEquation1Data::StateDerivatives::SecondDerivatives;
        using turbulence_equation_2_derivatives_data = typename TurbulenceEquation2Data::StateDerivatives::SecondDerivatives;

        turbulence_equation_1_derivatives_data turbulence_equation_1(this->GetGeometry());
        turbulence_equation_2_derivatives_data turbulence_equation_2(this->GetGeometry());

        turbulence_equation_1.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_2.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        RansCalculationUtilities::CalculateGeometryData(this->GetGeometry(), TurbulenceEquation1Data::GetIntegrationMethod(), gauss_weights, shape_functions, shape_derivatives);

        for (IndexType g = 0; g < gauss_weights.size(); ++g) {
            const Vector& N = row(shape_functions, g);
            const Matrix& dNdX = shape_derivatives[g];
            const double weight = gauss_weights[g];

            turbulence_equation_1.CalculateResidualDerivatives(rLeftHandSideMatrix, weight, N, dNdX);
            turbulence_equation_2.CalculateResidualDerivatives(rLeftHandSideMatrix, weight, N, dNdX);
        }

        turbulence_equation_1.Finalize(rLeftHandSideMatrix, rCurrentProcessInfo);
        turbulence_equation_2.Finalize(rLeftHandSideMatrix, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "TwoEquationTurbulenceModelAdjointElement::"
                        "CalculateMassMatrix method is not implemented.";

        KRATOS_CATCH("")
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                const ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "TwoEquationTurbulenceModelAdjointElement::"
                        "CalculateDampingMatrix method is not implemented.";

        KRATOS_CATCH("")
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY) {
        } else {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "TwoEquationTurbulenceModelAdjointElement #" << Element::Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TwoEquationTurbulenceModelAdjointElement #" << Element::Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        Element::pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected Operations
    ///@{
    ///@}
private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}
};

///@}

} // namespace Kratos

#endif