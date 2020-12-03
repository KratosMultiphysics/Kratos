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

#if !defined(KRATOS_FLUID_ADJOINT_ELEMENT_H)
#define KRATOS_FLUID_ADJOINT_ELEMENT_H

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "utilities/adjoint_extensions.h"
#include "utilities/geometrical_sensitivity_utility.h"

// Application includes
#include "custom_constitutive/fluid_adjoint_constitutive_law.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TAdjointElementData>
class FluidAdjointElement : public Element
{
    class ThisExtensions : public AdjointExtensions
    {
        Element* mpElement;

    public:
        explicit ThisExtensions(Element* pElement)
            : mpElement{pElement}
        {
        }

        void GetFirstDerivativesVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(TDim + 1);

            const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

            for (unsigned int i = 0; i < TDim; ++i) {
                rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative(), Step);
            }

            rVector[TDim] = IndirectScalar<double>{}; // pressure
        }

        void GetSecondDerivativesVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(TDim + 1);

            const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

            for (unsigned int i = 0; i < TDim; ++i) {
                rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative(), Step);
            }

            rVector[TDim] = IndirectScalar<double>{}; // pressure
        }

        void GetAuxiliaryVector(
            std::size_t NodeId,
            std::vector<IndirectScalar<double>>& rVector,
            std::size_t Step) override
        {
            auto& r_node = mpElement->GetGeometry()[NodeId];
            rVector.resize(TDim + 1);

            const auto& dofs_list = TAdjointElementData::GetDofVariablesList();

            for (unsigned int i = 0; i < TDim; ++i) {
                rVector[i] = MakeIndirectScalar(r_node, (*dofs_list[i]).GetTimeDerivative().GetTimeDerivative().GetTimeDerivative(), Step);
            }

            rVector[TDim] = IndirectScalar<double>{}; // pressure
        }

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(1);
            rVariables[0] = &ADJOINT_FLUID_VECTOR_2;
        }

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(1);
            rVariables[0] = &ADJOINT_FLUID_VECTOR_3;
        }

        void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override
        {
            rVariables.resize(1);
            rVariables[0] = &AUX_ADJOINT_FLUID_VECTOR_1;
        }
    };

public:
    ///@name Type Definitions
    ///@{

    using BaseType = Element;

    using NodesArrayType = typename BaseType::NodesArrayType;

    using PropertiesType = typename BaseType::PropertiesType;

    using GeometryType = typename BaseType::GeometryType;

    using VectorType = typename BaseType::VectorType;

    using MatrixType = typename BaseType::MatrixType;

    using IndexType = std::size_t;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    constexpr static IndexType TBlockSize = TAdjointElementData::TBlockSize;

    constexpr static IndexType TElementLocalSize = TBlockSize * TNumNodes;

    KRATOS_CLASS_POINTER_DEFINITION(FluidAdjointElement);

    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    FluidAdjointElement(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using Geometry
     */
    FluidAdjointElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    FluidAdjointElement(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Destructor
     */
    ~FluidAdjointElement() override
    {
    }

    ///@}
    ///@name Operations
    ///@{

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
    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<FluidAdjointElement>(
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
    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<FluidAdjointElement>(
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
    Element::Pointer Clone(
        IndexType NewId,
        NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY
        return Kratos::make_intrusive<FluidAdjointElement>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
        KRATOS_CATCH("");
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        return TAdjointElementData::Check(this->GetGeometry(), rCurrentProcessInfo);
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

        const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

        IndexType local_index = 0;
        for (IndexType i = 0; i < TNumNodes; ++i) {
            const auto& r_node = this->GetGeometry()[i];
            for (const auto p_variable : r_variables_list) {
                rElementalEquationIdList[local_index++] = r_node.GetDof(*p_variable).EquationId();
            }
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

        const auto& r_variables_list = TAdjointElementData::GetDofVariablesList();

        IndexType local_index = 0;
        for (IndexType i = 0; i < TNumNodes; ++i) {
            const auto& r_node = this->GetGeometry()[i];
            for (const auto p_variable : r_variables_list) {
                rElementalDofList[local_index++] = r_node.pGetDof(*p_variable);
            }
        }
    }

    /// Returns the adjoint values stored in this element's nodes.
    void GetValuesVector(
        VectorType& rValues,
        int Step = 0) const override
    {
        if (rValues.size() != TElementLocalSize) {
            rValues.resize(TElementLocalSize, false);
        }

        const auto& r_geometry = this->GetGeometry();
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const auto& r_node = r_geometry[i_node];
            const auto& r_velocity = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, Step);
            for (IndexType d = 0; d < TDim; ++d) {
                rValues[local_index++] = r_velocity[d];
            }
            rValues[local_index++] = r_node.FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, Step);
        }
    }

    /// Returns the adjoint velocity values stored in this element's nodes.
    void GetFirstDerivativesVector(
        VectorType& rValues,
        int Step = 0) const override
    {
        if (rValues.size() != TElementLocalSize) {
            rValues.resize(TElementLocalSize, false);
        }
        rValues.clear();
    }

    void GetSecondDerivativesVector(
        VectorType& rValues,
        int Step) const override
    {
        if (rValues.size() != TElementLocalSize) {
            rValues.resize(TElementLocalSize);
        }

        const auto& r_geometry = this->GetGeometry();
        IndexType local_index = 0;
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
            const auto& r_acceleration = r_geometry[i_node].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_3, Step);
            for (IndexType d = 0; d < TDim; ++d) {
                rValues[local_index++] = r_acceleration[d];
            }
            rValues[local_index++] = 0.0; // pressure dof
        }
    }

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override {
        KRATOS_TRY;

        // If we are restarting, the constitutive law will be already defined
        if (mpFluidConstitutiveLaw == nullptr) {
            const Properties& r_properties = this->GetProperties();

            KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
                << "In initialization of Element " << this->Info()
                << ": No CONSTITUTIVE_LAW defined for property "
                << r_properties.Id() << "." << std::endl;

            // Here we can do down casting because, it should be always a FluidConstitutiveLaw
            mpFluidConstitutiveLaw = std::static_pointer_cast<FluidConstitutiveLaw>(r_properties[CONSTITUTIVE_LAW]->Clone());

            const GeometryType& r_geometry = this->GetGeometry();
            const auto& r_shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);

            mpFluidConstitutiveLaw->InitializeMaterial(r_properties,r_geometry,row(r_shape_functions,0));
        }

        this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));

        KRATOS_CATCH("");
    }


    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "FluidAdjointElement::"
                        "CalculateLocalSystem method is not implemented.";

        KRATOS_CATCH("");
    }

    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize) {
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
        }

        rLeftHandSideMatrix.clear();
    }

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "FluidAdjointElement::"
                        "CalculateRightHandSide method is not implemented.";

        KRATOS_CATCH("");
    }

    void CalculateFirstDerivativesLHS(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize) {
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
        }

        rLeftHandSideMatrix.clear();

        AddVelocityPressureDerivatives<typename TAdjointElementData::StateDerivatives::FirstDerivatives>(
            rLeftHandSideMatrix,
            TAdjointElementData::GetIntegrationMethod(),
            rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesLHS(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TElementLocalSize ||
            rLeftHandSideMatrix.size2() != TElementLocalSize) {
            rLeftHandSideMatrix.resize(TElementLocalSize, TElementLocalSize, false);
        }

        rLeftHandSideMatrix.clear();

        AddVelocityPressureSecondDerivatives<typename TAdjointElementData::StateDerivatives::SecondDerivatives>(
            rLeftHandSideMatrix,
            TAdjointElementData::GetIntegrationMethod(),
            rCurrentProcessInfo);
    }

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        const ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "FluidAdjointElement::"
                        "CalculateMassMatrix method is not implemented.";

        KRATOS_CATCH("")
    }

    void CalculateDampingMatrix(
        MatrixType& rDampingMatrix,
        const ProcessInfo& /*rCurrentProcessInfo*/) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "FluidAdjointElement::"
                        "CalculateDampingMatrix method is not implemented.";

        KRATOS_CATCH("")
    }

    void CalculateSensitivityMatrix(
        const Variable<array_1d<double, 3>>& rSensitivityVariable,
        Matrix& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY) {
            using derivatives_type = typename TAdjointElementData::SensitivityDerivatives;
            constexpr IndexType shape_derivatives_size = TNumNodes * derivatives_type::ShapeSensitivities::TDerivativeDimension;

            if (rOutput.size1() != shape_derivatives_size || rOutput.size2() != TElementLocalSize) {
                rOutput.resize(shape_derivatives_size, TElementLocalSize, false);
            }

            rOutput.clear();
            AddVelocityPressureSensitivityDerivatives<derivatives_type>(
                rOutput, TAdjointElementData::GetIntegrationMethod(), rCurrentProcessInfo);
        } else {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "FluidAdjointElement #" << Element::Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "FluidAdjointElement #" << Element::Id();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        Element::pGetGeometry()->PrintData(rOStream);
    }

    ///@}

protected:
    ///@name Protected Members
    ///@{

    FluidConstitutiveLaw::Pointer mpFluidConstitutiveLaw = nullptr;

    ///@}
    ///@name Protected Operations
    ///@{

    template<class TDerivativesType>
    void AddVelocityPressureDerivatives(
        MatrixType& rLeftHandSideMatrix,
        const GeometryData::IntegrationMethod& rIntegrationMethod,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives, rIntegrationMethod);

        typename TDerivativesType::Data element_data(*this, *mpFluidConstitutiveLaw);
        typename TDerivativesType::Velocity velocity_derivatives(element_data);
        typename TDerivativesType::Pressure pressure_derivatives(element_data);

        velocity_derivatives.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
        pressure_derivatives.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);
        element_data.Initialize(rCurrentProcessInfo);

        BoundedVector<double, TElementLocalSize> residual;
        BoundedMatrix<double, TNumNodes, TDim> dNdXDerivative = ZeroMatrix(TNumNodes, TDim);

        for (IndexType g = 0; g < gauss_weights.size(); ++g) {
            const Vector& N = row(shape_functions, g);
            const Matrix& dNdX = shape_derivatives[g];
            const double weight = gauss_weights[g];

            element_data.CalculateGaussPointData(weight, N, dNdX);
            for (IndexType c = 0; c < TNumNodes; ++c) {
                const IndexType derivative_block = c * TBlockSize;

                for (IndexType k = 0; k < TDim; ++k) {
                    velocity_derivatives.CalculateResidualDerivative(residual, c, k, weight, N, dNdX, 0.0, 0.0, dNdXDerivative);
                    row(rLeftHandSideMatrix, derivative_block + k) += residual;
                }

                pressure_derivatives.CalculateResidualDerivative(residual, c, 0, weight, N, dNdX, 0.0, 0.0, dNdXDerivative);
                row(rLeftHandSideMatrix, derivative_block + TDim) += residual;
            }
        }

        KRATOS_CATCH("");
    }

    template<class TDerivativesType>
    void AddVelocityPressureSecondDerivatives(
        MatrixType& rLeftHandSideMatrix,
        const GeometryData::IntegrationMethod& rIntegrationMethod,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives, rIntegrationMethod);

        TDerivativesType element_data(*this, *mpFluidConstitutiveLaw);

        element_data.Initialize(rLeftHandSideMatrix, rCurrentProcessInfo);

        for (IndexType g = 0; g < gauss_weights.size(); ++g) {
            const Vector& N = row(shape_functions, g);
            const Matrix& dNdX = shape_derivatives[g];
            const double weight = gauss_weights[g];

            element_data.AddResidualDerivativeContributions(rLeftHandSideMatrix, weight, N, dNdX);
        }

        KRATOS_CATCH("");
    }

    template<class TDerivativesType>
    void AddVelocityPressureSensitivityDerivatives(
        Matrix& rOutput,
        const GeometryData::IntegrationMethod& rIntegrationMethod,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives, rIntegrationMethod);

        typename TDerivativesType::Data element_data(*this, *mpFluidConstitutiveLaw);
        typename TDerivativesType::ShapeSensitivities derivatives(element_data);

        derivatives.Initialize(rOutput, rCurrentProcessInfo);
        element_data.Initialize(rCurrentProcessInfo);

        BoundedVector<double, TElementLocalSize> residual;
        BoundedMatrix<double, TNumNodes, TDim> dNdXDerivative = ZeroMatrix(TNumNodes, TDim);

        for (IndexType g = 0; g < gauss_weights.size(); ++g) {
            const Vector& N = row(shape_functions, g);
            const Matrix& dNdX = shape_derivatives[g];
            const double weight = gauss_weights[g];

            element_data.CalculateGaussPointData(weight, N, dNdX);

            Geometry<Point>::JacobiansType J;
            this->GetGeometry().Jacobian(J, rIntegrationMethod);
            const auto& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(rIntegrationMethod);

            GeometricalSensitivityUtility::ShapeFunctionsGradientType dNdX_deriv;
            const Matrix& rJ = J[g];
            const Matrix& rDN_De = DN_De[g];
            const double inv_detJ = 1.0 / MathUtils<double>::DetMat(rJ);
            GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

            ShapeParameter deriv;
            for (deriv.NodeIndex = 0; deriv.NodeIndex < TNumNodes; ++deriv.NodeIndex) {
                const IndexType derivative_row = deriv.NodeIndex * TDerivativesType::ShapeSensitivities::TDerivativeDimension;
                for (deriv.Direction = 0; deriv.Direction < TDerivativesType::ShapeSensitivities::TDerivativeDimension; ++deriv.Direction) {

                    double detJ_deriv;
                    geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, dNdX_deriv);
                    const double weight_deriv = detJ_deriv * inv_detJ * weight;

                    derivatives.CalculateResidualDerivative(residual, deriv.NodeIndex, deriv.Direction, weight, N, dNdX, weight_deriv, detJ_deriv, dNdX_deriv);
                    row(rOutput, derivative_row + deriv.Direction) += residual;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void CalculateGeometryData(
        Vector& rGaussWeights,
        Matrix& rNContainer,
        ShapeFunctionDerivativesArrayType& rDN_DX,
        const GeometryData::IntegrationMethod& rIntegrationMethod) const
    {
        const auto& r_geometry = this->GetGeometry();
        const IndexType number_of_gauss_points = r_geometry.IntegrationPointsNumber(rIntegrationMethod);

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


    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_FLUID_ADJOINT_ELEMENT_H