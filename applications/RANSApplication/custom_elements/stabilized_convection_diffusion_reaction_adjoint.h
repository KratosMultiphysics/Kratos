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

#if !defined(KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT)
#define KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT

// System includes
#include <array>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/element.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/time_discretization.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"
#include "rans_application_variables.h"
#include "stabilized_convection_diffusion_reaction.h"
#include "stabilized_convection_diffusion_reaction_adjoint_utilities.h"
#include "stabilized_convection_diffusion_reaction_utilities.h"

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

/**
 * @brief Discrete stabilized Convection-Diffusion-Reaction adjoint element
 *
 * This class provides derivatives with respect to velocity, pressure,
 * local coordinates and the scalar variable which is being used in the
 * element as the transport variable
 *
 * The sensitivities are calculated on the discrete version of the stabilized
 * equation given below
 *
 * \[
 *	\frac{\partial\phi}{\partial t} + \underline{u} \cdot \frac{\partial \phi}{\partial\underline{x}} - \nu_\phi \frac{\partial^2\phi}{\partial\underline{x}^2} + s_\phi \phi = f_\phi
 * \]
 *
 * Where, $\phi$ is the scalar variable, and $\underline{u}$, $\nu_\phi$,
 * $s_\phi$, $f_\phi$ are velocity, effective kinematic viscosity,
 * reaction coefficient, source term respectively.
 *
 * @tparam TDim                            Dimensionality of the element
 * @tparam TNumNodes                       Number of nodes in the element
 * @tparam TElementData                    Data container used to calculate derivatives
 */
template <unsigned int TDim, unsigned int TNumNodes, class TElementData>
class StabilizedConvectionDiffusionReactionAdjointElement
    : public StabilizedConvectionDiffusionReaction<TDim, TNumNodes, TElementData>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Element
    KRATOS_CLASS_POINTER_DEFINITION(StabilizedConvectionDiffusionReactionAdjointElement);

    /// base type: an GeometricalObject that automatically has a unique number
    typedef StabilizedConvectionDiffusionReaction<TDim, TNumNodes, TElementData> BaseType;

    /// definition of node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    typedef Properties PropertiesType;

    /// definition of the geometry type with given NodeType
    typedef Geometry<NodeType> GeometryType;

    /// definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Type definition for integration methods
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    typedef GeometryData GeometryDataType;

    typedef BoundedMatrix<double, TDim, TDim> BoundedMatrixDD;

    typedef BoundedMatrix<double, TNumNodes, TDim> BoundedMatrixND;

    typedef BoundedVector<double, TNumNodes> BoundedVectorN;

    using ConvectionDiffusionReactionAdjointDataType = TElementData;

    ///@}

    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        GeometryType::Pointer pGeometry,
                                                        PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    StabilizedConvectionDiffusionReactionAdjointElement(
        StabilizedConvectionDiffusionReactionAdjointElement const& rOther)
        : BaseType(rOther)
    {
    }

    ~StabilizedConvectionDiffusionReactionAdjointElement() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.

    StabilizedConvectionDiffusionReactionAdjointElement& operator=(
        StabilizedConvectionDiffusionReactionAdjointElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        // mpProperties = rOther.mpProperties;
        return *this;
    }

    ///@}
    ///@name Informations
    ///@{

    ///@}
    ///@name Operations
    ///@{
    /**
     * @brief It creates a new element pointer
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to Create base "
               "StabilizedConvectionDiffusionReactionAdjointElement instances."
            << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * @brief It creates a new element pointer
     * @param NewId the ID of the new element
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to Create base "
               "StabilizedConvectionDiffusionReactionAdjointElement instances."
            << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * @brief It creates a new element pointer and clones the previous element data
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to Clone base "
               "StabilizedConvectionDiffusionReactionAdjointElement instances."
            << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult the elemental equation ID vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override
    {
        std::array<std::size_t, TNumNodes> ids;
        this->EquationIdArray(ids, rCurrentProcessInfo);

        if (rResult.size() != TNumNodes)
            rResult.resize(TNumNodes, false);
        std::copy(ids.begin(), ids.end(), rResult.begin());
    }

    void EquationIdArray(std::array<std::size_t, TNumNodes>& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        const Variable<double>& r_adjoint_Variable = this->GetAdjointVariable();
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            rResult[i] = this->GetGeometry()[i].GetDof(r_adjoint_Variable).EquationId();
        }
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        std::array<Dof<double>::Pointer, TNumNodes> dofs;
        this->GetDofArray(dofs, rCurrentProcessInfo);

        if (rElementalDofList.size() != TNumNodes)
            rElementalDofList.resize(TNumNodes);

        std::copy(dofs.begin(), dofs.end(), rElementalDofList.begin());
    }

    void GetDofArray(std::array<Dof<double>::Pointer, TNumNodes>& rElementalDofList,
                     ProcessInfo& rCurrentProcessInfo)
    {
        const Variable<double>& r_adjoint_Variable = this->GetAdjointVariable();
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            rElementalDofList[i] = this->GetGeometry()[i].pGetDof(r_adjoint_Variable);
        }
    }

    void GetValuesVector(VectorType& rValues, int Step = 0) override
    {
        std::array<double, TNumNodes> values;
        this->GetValuesArray(values, Step);

        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);
        std::copy(values.begin(), values.end(), rValues.begin());
    }

    void GetValuesArray(std::array<double, TNumNodes>& rValues, int Step = 0)
    {
        const Variable<double>& r_adjoint_variable = this->GetAdjointVariable();
        const GeometryType& r_geometry = this->GetGeometry();
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            rValues[i] = r_geometry[i].FastGetSolutionStepValue(r_adjoint_variable, Step);
        }
    }

    void GetFirstDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        std::array<double, TNumNodes> values;
        this->GetFirstDerivativesArray(values, Step);
        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);
        std::copy(values.begin(), values.end(), rValues.begin());
    }

    void GetFirstDerivativesArray(std::array<double, TNumNodes>& rValues, int Step = 0)
    {
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            rValues[i] = 0.0;
        }
    }

    void GetSecondDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        std::array<double, TNumNodes> values;
        this->GetSecondDerivativesArray(values, Step);
        if (rValues.size() != TNumNodes)
            rValues.resize(TNumNodes, false);
        std::copy(values.begin(), values.end(), rValues.begin());
    }

    void GetSecondDerivativesArray(std::array<double, TNumNodes>& rValues, int Step = 0)
    {
        const Variable<double>& r_adjoint_second_variable =
            this->GetAdjointSecondVariable();
        const GeometryType& r_geometry = this->GetGeometry();
        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            rValues[i] = r_geometry[i].FastGetSolutionStepValue(
                r_adjoint_second_variable, Step);
        }
    }

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rRightHandSideVector the elemental right hand side
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "StabilizeConvectionDiffusionReactionAdjointElement::"
                        "CalculateLocalSystem method not implemented.";

        KRATOS_CATCH("");
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix the elemental left hand side matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        rLeftHandSideMatrix.clear();
    }

    void CalculateLeftHandSide(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo)
    {
        rLeftHandSideMatrix.clear();
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector the elemental right hand side vector
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "StabilizeConvectionDiffusionReactionAdjointElement::"
                        "CalculateRightHandSide method not implemented.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates the adjoint matrix for scalar variable
     *
     * This function returns the gradient of the elemental residual w.r.t.
     * scalar variable transposed:
     *
     * \f[
     *    \partial_{\mathbf{w}^n}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\mathbf{w}^n}(\mathbf{M}^n \dot{\mathbf{w}}^n)^T
     * \f]
     *
     * where \f$\mathbf{w}^n\f$ is the vector of nodal scalar
     * stored at the current step. For steady problems, the scalar rate variable
     * (\f$\dot{\mathbf{w}}^n\f$) must be set to zero on the nodes. For
     * the Bossak method, \f$\dot{\mathbf{w}}^{n-\alpha}\f$ must be stored in
     * the variable given by @GetPrimalRelaxedRateVariable().
     *
     * @param rLeftHandSideMatrix
     * @param rCurrentProcessInfo
     */
    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculateFirstDerivativesLHS(local_matrix, rCurrentProcessInfo);
        if (rLeftHandSideMatrix.size1() != local_matrix.size1() ||
            rLeftHandSideMatrix.size2() != local_matrix.size2())
            rLeftHandSideMatrix.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rLeftHandSideMatrix) = local_matrix;
    }

    void CalculateFirstDerivativesLHS(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo)
    {
        const Variable<double>& r_derivative_variable = this->GetPrimalVariable();
        CalculateResidualScalarDerivatives(
            r_derivative_variable, rLeftHandSideMatrix, rCurrentProcessInfo);
        AddPrimalDampingMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    }

    /**
     * @brief Calculate the adjoint matrix for scalar rate
     *
     * This function returns the gradient of the elemental residual w.r.t.
     * scalar rate variable:
     *
     * \f[
     *    \partial_{\dot{\mathbf{w}}^n}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\dot{\mathbf{w}}^n}(\mathbf{M}^n \dot{\mathbf{w}}^n)^T
     * \f]
     *
     * @param rLeftHandSideMatrix
     * @param rCurrentProcessInfo
     */
    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculateSecondDerivativesLHS(local_matrix, rCurrentProcessInfo);

        if (rLeftHandSideMatrix.size1() != local_matrix.size1() ||
            rLeftHandSideMatrix.size2() != local_matrix.size2())
            rLeftHandSideMatrix.resize(local_matrix.size1(), local_matrix.size2(), false);
        rLeftHandSideMatrix.clear();

        noalias(rLeftHandSideMatrix) = local_matrix;
    }

    void CalculateSecondDerivativesLHS(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo)
    {
        rLeftHandSideMatrix.clear();
        AddPrimalMassMatrix(rLeftHandSideMatrix, -1.0, rCurrentProcessInfo);
        AddPrimalSteadyTermScalarRateDerivatives(rLeftHandSideMatrix, rCurrentProcessInfo);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix the elemental mass matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "StabilizeConvectionDiffusionReactionAdjointElement::"
                        "CalculateMassMatrix method not implemented.";

        KRATOS_CATCH("");
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental damping matrix
     * @param rDampingMatrix the elemental damping matrix
     * @param rCurrentProcessInfo the current process info instance
     */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR << "StabilizeConvectionDiffusionReactionAdjointElement::"
                        "CalculateDampingMatrix method not implemented.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates the adjoint matrix for scalar variable
     *
     * This function returns the gradient of the elemental residual w.r.t.
     * velocity and pressure transposed:
     *
     * \f[
     *    \partial_{\mathbf{w}^n}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\mathbf{w}^n}(\mathbf{M}^n \dot{\mathbf{w}}^n)^T
     * \f]
     *
     * where \f$\mathbf{w}^n\f$ is the vector of nodal velocity and pressure
     * stored at the current step.
     *
     * @param rVariable
     * @param Output
     * @param rCurrentProcessInfo
     */

    void Calculate(const Variable<Matrix>& rVariable,
                   Matrix& rOutput,
                   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rVariable == RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE)
        {
            BoundedMatrix<double, TNumNodes * TDim, TNumNodes> local_matrix;
            CalculateElementTotalResidualVelocityDerivatives(local_matrix, rCurrentProcessInfo);

            if (rOutput.size1() != local_matrix.size1() ||
                rOutput.size2() != local_matrix.size2())
                rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

            noalias(rOutput) = local_matrix;
        }
        else
        {
            KRATOS_ERROR << "Unsupported variable "
                         << rVariable.Name() << " requested at StabilizedConvectionDiffusionReactionAdjoint::Calculate.";
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates the sensitivity matrix.
     *
     * \f[
     *    \partial_{\mathbf{s}}\mathbf{f}(\mathbf{w}^n)^T
     *  - \partial_{\mathbf{s}}(\mathbf{M}^n \dot{\mathbf{w}}^{n-\alpha})^T
     * \f]
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            BoundedMatrix<double, TNumNodes * TDim, TNumNodes> local_matrix;
            this->CalculateResidualShapeSensitivity(local_matrix, rCurrentProcessInfo);
            if (rOutput.size1() != local_matrix.size1() ||
                rOutput.size2() != local_matrix.size2())
                rOutput.resize(local_matrix.size1(), local_matrix.size2(), false);

            noalias(rOutput) = local_matrix;
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            this->CalculateResidualShapeSensitivity(rOutput, rCurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
        }

        KRATOS_CATCH("")
    }

    void CalculateElementTotalResidualVelocityDerivatives(
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rResidualDerivatives,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rResidualDerivatives.clear();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedVector<double, TNumNodes> primal_variable_gradient_convective_terms,
            velocity_convective_terms;

        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;

        BoundedMatrix<double, TNumNodes, TDim> effective_kinematic_viscosity_derivatives,
            reaction_derivatives, source_derivatives, element_length_derivatives,
            tau_derivatives, psi_one_derivatives, psi_two_derivatives, chi_derivatives,
            residual_derivatives, positivity_preservation_coeff_derivatives,
            stream_line_diffusion_coeff_derivatives,
            cross_wind_diffusion_coeff_derivatives,
            velocity_magnitude_derivatives, absolute_residual_derivatives,
            absolute_reaction_tilde_derivatives, s_derivatives;

        array_1d<double, 3> primal_variable_gradient;

        TElementData current_data;

        for (unsigned int g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            this->CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            this->CalculateElementData(current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            const double primal_variable_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);
            this->CalculateGradient(primal_variable_gradient, primal_variable,
                                    r_shape_derivatives);
            this->GetConvectionOperator(primal_variable_gradient_convective_terms,
                                        primal_variable_gradient, r_shape_derivatives);

            const double effective_kinematic_viscosity = this->CalculateEffectiveKinematicViscosity(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityVelocityDerivatives(
                effective_kinematic_viscosity_derivatives, current_data,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction = this->CalculateReactionTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            this->CalculateReactionTermVelocityDerivatives(
                reaction_derivatives, current_data, gauss_shape_functions,
                r_shape_derivatives, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            const double velocity_magnitude = norm_2(velocity);
            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateVelocityMagnitudeVelocityDerivative(
                velocity_magnitude_derivatives, velocity_magnitude, velocity,
                gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateElementLengthH2VelocityDerivative(
                element_length_derivatives, velocity_magnitude, velocity,
                velocity_magnitude_derivatives, contravariant_metric_tensor,
                gauss_shape_functions);

            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStabilizationTauVelocityDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction,
                element_length, velocity, contravariant_metric_tensor,
                effective_kinematic_viscosity_derivatives, reaction_derivatives,
                element_length_derivatives, gauss_shape_functions);

            const double source = this->CalculateSourceTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            this->CalculateSourceTermVelocityDerivatives(
                source_derivatives, current_data, gauss_shape_functions,
                r_shape_derivatives, rCurrentProcessInfo);

            const double velocity_dot_primal_variable_gradient =
                inner_prod(velocity, primal_variable_gradient);

            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double primal_variable_gradient_norm = norm_2(primal_variable_gradient);

            double chi{0.0}, k1{0.0}, k2{0.0}, residual{0.0},
                positivity_preservation_coeff{0.0};

            const double primal_variable_relaxed_rate = this->EvaluateInPoint(
                this->GetPrimalRelaxedRateVariable(), gauss_shape_functions);

            if (velocity_magnitude_square > 0.0 && primal_variable_gradient_norm > 0.0)
            {
                residual = primal_variable_relaxed_rate;
                residual += velocity_dot_primal_variable_gradient;
                residual += reaction * primal_variable_value;
                residual -= source;

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);

                positivity_preservation_coeff =
                    std::abs(residual) * chi /
                    (velocity_magnitude_square * primal_variable_gradient_norm);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateChiVelocityDerivatives(
                    chi_derivatives, chi, element_length, bossak_alpha, bossak_gamma,
                    delta_time, reaction, dynamic_tau, reaction_derivatives,
                    velocity_magnitude_derivatives, element_length_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateResidualVelocityDerivative(
                    residual_derivatives, primal_variable_value, primal_variable_gradient,
                    reaction_derivatives, source_derivatives, gauss_shape_functions);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueVectorDerivatives(
                    absolute_residual_derivatives, residual, residual_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePositivityPreservationCoefficientVelocityDerivatives(
                    positivity_preservation_coeff_derivatives, std::abs(residual),
                    primal_variable_gradient_norm, velocity_magnitude, chi, chi_derivatives,
                    absolute_residual_derivatives, velocity_magnitude_derivatives);

                const double reaction_tilde =
                    reaction + dynamic_tau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueVectorDerivatives(
                    absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

                const double psi_one =
                    StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
                        velocity_magnitude, tau, reaction_tilde);
                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiOneVelocityDerivatives(
                    psi_one_derivatives, velocity_magnitude, reaction_tilde,
                    tau, tau_derivatives, absolute_reaction_tilde_derivatives,
                    velocity_magnitude_derivatives);

                const double psi_two =
                    StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
                        reaction_tilde, tau, element_length);
                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiTwoVelocityDerivatives(
                    psi_two_derivatives, reaction_tilde, tau, element_length,
                    tau_derivatives, reaction_derivatives,
                    absolute_reaction_tilde_derivatives, element_length_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStreamLineDiffusionCoeffVelocityDerivatives(
                    stream_line_diffusion_coeff_derivatives, k1, element_length,
                    tau, velocity_magnitude, reaction_tilde, psi_one, psi_two,
                    velocity_magnitude_derivatives, psi_one_derivatives,
                    psi_two_derivatives, tau_derivatives, reaction_derivatives,
                    effective_kinematic_viscosity_derivatives, element_length_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateCrossWindDiffusionCoeffVelocityDerivatives(
                    cross_wind_diffusion_coeff_derivatives, k2, psi_one,
                    element_length, psi_one_derivatives, psi_two_derivatives,
                    effective_kinematic_viscosity_derivatives, element_length_derivatives);
            }
            else
            {
                positivity_preservation_coeff_derivatives.clear();
                cross_wind_diffusion_coeff_derivatives.clear();
                stream_line_diffusion_coeff_derivatives.clear();
            }

            const double s = std::abs(reaction);
            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueVectorDerivatives(
                s_derivatives, reaction, reaction_derivatives);

            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    for (unsigned int k = 0; k < TDim; ++k)
                    {
                        // adding damping matrix derivatives
                        double value = 0.0;

                        // convection term derivatives
                        value += gauss_shape_functions[a] * gauss_shape_functions[c] *
                                 primal_variable_gradient[k];

                        // adding diffusion term derivatives
                        value += effective_kinematic_viscosity_derivatives(c, k) *
                                 primal_variable_gradient_convective_terms[a];

                        // adding reaction term derivatives
                        value += reaction_derivatives(c, k) *
                                 gauss_shape_functions[a] * primal_variable_value;

                        // adding SUPG term derivatives
                        value += tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 velocity_dot_primal_variable_gradient;
                        value += tau *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) +
                                  s_derivatives(c, k) * gauss_shape_functions[a]) *
                                 velocity_dot_primal_variable_gradient;
                        value +=
                            tau *
                            (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                            gauss_shape_functions[c] * primal_variable_gradient[k];

                        value += tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 reaction * primal_variable_value;
                        value += tau *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) +
                                  s_derivatives(c, k) * gauss_shape_functions[a]) *
                                 reaction * primal_variable_value;
                        value += tau *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 reaction_derivatives(c, k) * primal_variable_value;

                        // adding cross wind dissipation derivatives
                        value += positivity_preservation_coeff_derivatives(c, k) *
                                 k2 * velocity_magnitude_square *
                                 primal_variable_gradient_convective_terms[a];
                        value += positivity_preservation_coeff *
                                 cross_wind_diffusion_coeff_derivatives(c, k) *
                                 velocity_magnitude_square *
                                 primal_variable_gradient_convective_terms[a];
                        value += positivity_preservation_coeff * k2 * 2 * velocity_magnitude *
                                 velocity_magnitude_derivatives(c, k) *
                                 primal_variable_gradient_convective_terms[a];
                        value -= positivity_preservation_coeff_derivatives(c, k) *
                                 k2 * velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient;
                        value -= positivity_preservation_coeff *
                                 cross_wind_diffusion_coeff_derivatives(c, k) *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient;
                        value -= positivity_preservation_coeff * k2 *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) *
                                  velocity_dot_primal_variable_gradient);
                        value -= positivity_preservation_coeff * k2 *
                                 (gauss_shape_functions[c] * velocity_convective_terms[a] *
                                  primal_variable_gradient[k]);

                        // adding stream line dissipation derivatives
                        value += positivity_preservation_coeff_derivatives(c, k) *
                                 k1 * velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient;
                        value += positivity_preservation_coeff *
                                 stream_line_diffusion_coeff_derivatives(c, k) *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient;
                        value += positivity_preservation_coeff * k1 *
                                 gauss_shape_functions[c] * r_shape_derivatives(a, k) *
                                 velocity_dot_primal_variable_gradient;
                        value += positivity_preservation_coeff * k1 *
                                 velocity_convective_terms[a] *
                                 gauss_shape_functions[c] * primal_variable_gradient[k];

                        // putting transposed values
                        rResidualDerivatives(c * TDim + k, a) -=
                            value * gauss_weights[g];

                        // adding source term derivatives
                        value = 0.0;

                        value += gauss_shape_functions[a] * source_derivatives(c, k);
                        value += tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 source;
                        value += tau *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) +
                                  s_derivatives(c, k) * gauss_shape_functions[a]) *
                                 source;
                        value += tau *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 source_derivatives(c, k);

                        // putting transposed values
                        rResidualDerivatives(c * TDim + k, a) +=
                            value * gauss_weights[g];

                        // adding mass term derivatives
                        value = 0.0;

                        value += tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  s * gauss_shape_functions[a]) *
                                 primal_variable_relaxed_rate;
                        value += tau *
                                 (gauss_shape_functions[c] * r_shape_derivatives(a, k) +
                                  s_derivatives(c, k) * gauss_shape_functions[a]) *
                                 primal_variable_relaxed_rate;

                        rResidualDerivatives(c * TDim + k, a) -=
                            value * gauss_weights[g];
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     */

    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(this->Id() < 1) << "StabilizedConvectionDiffusionReacti"
                                           "onAdjointElement found with Id 0 "
                                           "or negative"
                                        << std::endl;

        KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0)
            << "On StabilizedConvectionDiffusionReactionAdjointElement -> "
            << this->Id() << "; Area cannot be less than or equal to 0" << std::endl;

        const Variable<double>& r_primal_variable = this->GetPrimalVariable();
        const Variable<double>& r_primal_relaxed_rate_variable =
            this->GetPrimalRelaxedRateVariable();
        const Variable<double>& r_adjoint_variable = this->GetAdjointVariable();

        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            NodeType& r_node = this->GetGeometry()[iNode];
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_primal_variable, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_primal_relaxed_rate_variable, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_adjoint_variable, r_node);

            KRATOS_CHECK_DOF_IN_NODE(r_adjoint_variable, r_node);
        }

        return 0;

        KRATOS_CATCH("")
    }

    /**
     * @brief Calculate element residual scalar derivatives
     *
     * This method calculates transposed element residuals scalar partial derivatives.
     * In the steady regime, the the variable given by
     * @GetPrimalRelaxedRateVariable() should have zero values.
     * This method is called for each gauss point.
     *
     * This method should be implemented by the derrived class.
     *
     * @param rResidualDerivatives Residual derivatives, rows contains derivative variables, columns contains residual equations.
     * @param rDerivativeVariable  Scalar derivative variables
     * @param rCurrentProcessInfo  Current process info
     */
    void CalculateResidualScalarDerivatives(const Variable<double>& rDerivativeVariable,
                                            BoundedMatrix<double, TNumNodes, TNumNodes>& rResidualDerivatives,
                                            const ProcessInfo& rCurrentProcessInfo)
    {
        rResidualDerivatives.clear();
        AddPrimalSteadyTermScalarDerivatives(
            rResidualDerivatives, rDerivativeVariable, rCurrentProcessInfo);
        AddMassTermScalarDerivatives(rResidualDerivatives, rDerivativeVariable,
                                     rCurrentProcessInfo);
    }

    /**
     * @brief Get the adjoint variable
     *
     * This method returns the adjoint variable.
     *
     * This method should be implemented by the derrived class.
     *
     * @return const Variable<double>&
     */
    virtual const Variable<double>& GetAdjointVariable() const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base GetAdjointVariable method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        return RANS_SCALAR_1_ADJOINT_1;

        KRATOS_CATCH("");
    }

    /**
     * @brief Get the second adjoint variable
     *
     * This method returns the second adjoint variable.
     *
     * This method should be implemented by the derrived class.
     *
     * @return const Variable<double>&
     */
    virtual const Variable<double>& GetAdjointSecondVariable() const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base GetAdjointSecondVariable method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        return RANS_SCALAR_1_ADJOINT_1;

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates scalar partial derivatives of effective kinematic viscosity
     *
     * This method calculates partial derivative of the effective kinematic viscosity
     * for a given scalar variable. This method is called for each gauss point.
     *
     * \[
     *   \frac{\partial\nu_\phi}{\partial w}
     * \]
     *
     * Where $\nu_\phi$ is the effective kinematic viscosity, and $w$ is the derivative
     * variable
     *
     * This method should be implemented by the derrived class.
     *
     * @param rOutput              Output vector containing partial derivatives for each node
     * @param rDerivativeVariable  Derivative variable (i.e. $w$)
     * @param rCurrentData         Data required to calculate partial derivatives
     * @param rCurrentProcessInfo  Current process info
     */
    virtual void CalculateEffectiveKinematicViscosityScalarDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base "
                        "CalculateEffectiveKinematicViscosityScalarDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates scalar partial derivatives of reaction coefficient
     *
     * This method calculates partial derivative of the reaction coefficient
     * for a given scalar variable. This method is called for each gauss point.
     *
     * \[
     *   \frac{\partial s_\phi}{\partial w}
     * \]
     *
     * Where $s_\phi$ is the reaction coefficient, and $w$ is the derivative
     * variable
     *
     * This method should be implemented by the derrived class.
     *
     * @param rOutput              Output vector containing partial derivatives for each node
     * @param rDerivativeVariable  Derivative variable (i.e. $w$)
     * @param rCurrentData         Data required to calculate partial derivatives
     * @param rCurrentProcessInfo  Current process info
     */
    virtual void CalculateReactionTermScalarDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateReactionTermScalarDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates scalar partial derivatives of source term
     *
     * This method calculates partial derivative of the source term
     * for a given scalar variable. This method is called for each gauss point.
     *
     * \[
     *   \frac{\partial f_\phi}{\partial w}
     * \]
     *
     * Where $s_\phi$ is the source term, and $w$ is the derivative
     * variable
     *
     * This method should be implemented by the derrived class.
     *
     * @param rOutput                      Output vector containing partial derivatives for each node
     * @param rDerivativeVariable          Derivative variable (i.e. $w$)
     * @param rCurrentData                 Data required to calculate partial derivatives
     * @param rShapeFunctions              Shape functions of the gauss point
     * @param rShapeFunctionDerivatives    Shape function derivatives of the gauss point
     * @param rCurrentProcessInfo  Current process info
     */
    virtual void CalculateSourceTermScalarDerivatives(
        BoundedVector<double, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateSourceTermScalarDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculate velocity derivatives of effective kinematic viscosity
     *
     * This method calculates partial velocity derivatives of effective kinematic viscosity
     * This method is called for each gauss point.
     *
     * \[
     *   \frac{\partial \nu_\phi}{\partial \underline{u}}
     * \]
     *
     * This method should be implemented by the derrived class.
     *
     * @param rOutput                      Output vector containing partial derivatives for each node in rows, and each dimension in columns
     * @param rCurrentData                 Data required to calculate partial derivatives
     * @param rShapeFunctions              Shape functions of the gauss point
     * @param rShapeFunctionDerivatives    Shape function derivatives of the gauss point
     * @param rCurrentProcessInfo          Current process info
     */
    virtual void CalculateEffectiveKinematicViscosityVelocityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR
            << "Calling base "
               "CalculateEffectiveKinematicViscosityVelocityDerivatives "
               "method in "
               "StabilizedConvectionDiffusionReactionAdjointElement "
               "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculate velocity derivatives of reaction coefficient
     *
     * This method calculates partial velocity derivatives of reaction coefficient
     * This method is called for each gauss point.
     *
     * \[
     *   \frac{\partial s_\phi}{\partial \underline{u}}
     * \]
     *
     * This method should be implemented by the derrived class.
     *
     * @param rOutput                      Output vector containing partial derivatives for each node in rows, and each dimension in columns
     * @param rCurrentData                 Data required to calculate partial derivatives
     * @param rShapeFunctions              Shape functions of the gauss point
     * @param rShapeFunctionDerivatives    Shape function derivatives of the gauss point
     * @param rCurrentProcessInfo          Current process info
     */
    virtual void CalculateReactionTermVelocityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateReactionTermVelocityDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculate velocity derivatives of source term
     *
     * This method calculates partial velocity derivatives of source term
     * This method is called for each gauss point.
     *
     * \[
     *   \frac{\partial f_\phi}{\partial \underline{u}}
     * \]
     *
     * This method should be implemented by the derrived class.
     *
     * @param rOutput                      Output vector containing partial derivatives for each node in rows, and each dimension in columns
     * @param rCurrentData                 Data required to calculate partial derivatives
     * @param rShapeFunctions              Shape functions of the gauss point
     * @param rShapeFunctionDerivatives    Shape function derivatives of the gauss point
     * @param rCurrentProcessInfo          Current process info
     */
    virtual void CalculateSourceTermVelocityDerivatives(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateSourceTermVelocityDerivatives "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates shape sensitivity of effective kinematic viscosity
     *
     * Calculates shape derivative of effective kinematic viscosity
     * w.r.t. nodal coordinates. This method is called for each nodal coordinate,
     * its different dimensions and for each gauss point as well.
     *
     * This method should be implemented by the derrived class.
     *
     * @param rCurrentData                 Data required to calculate partial derivatives
     * @param rShapeFunctions              Shape functions of the gauss point
     * @param rShapeFunctionDerivatives    Shape function derivatives of the gauss point
     * @param rShapeDerivative             Current derivative (i.e. $x^c_k$, where $c$ is the node, $k$ is the dimension)
     * @param detJ_deriv                   Derivative of determinant of jacobian (i.e. $|J|$) w.r.t. $x^c_k$
     * @param rDN_Dx_deriv                 Derivative of shape function gradients w.r.t. $x^c_k$
     * @param rCurrentProcessInfo          Current process info
     * @return double                      Effective kinematic viscosity shape derivative w.r.t. $x^c_k$
     */
    virtual double CalculateEffectiveKinematicViscosityShapeSensitivity(
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base "
                        "CalculateEffectiveKinematicViscosityShapeSensitivity "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates shape sensitivity of reaction coefficient
     *
     * Calculates shape derivative of reaction coefficient
     * w.r.t. nodal coordinates. This method is called for each nodal coordinate,
     * its different dimensions and for each gauss point as well.
     *
     * This method should be implemented by the derrived class.
     *
     * @param rCurrentData                 Data required to calculate partial derivatives
     * @param rShapeFunctions              Shape functions of the gauss point
     * @param rShapeFunctionDerivatives    Shape function derivatives of the gauss point*
     * @param rShapeDerivative             Current derivative (i.e. $x^c_k$, where $c$ is the node, $k$ is the dimension)
     * @param detJ_deriv                   Derivative of determinant of jacobian (i.e. $|J|$) w.r.t. $x^c_k$
     * @param rDN_Dx_deriv                 Derivative of shape function gradients w.r.t. $x^c_k$
     * @param rCurrentProcessInfo          Current process info
     * @return double                      Reaction coefficient shape derivative w.r.t. $x^c_k$
     */
    virtual double CalculateReactionTermShapeSensitivity(
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateReactionTermShapeSensitivity "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates shape sensitivity of source term
     *
     * Calculates shape derivative of source term
     * w.r.t. nodal coordinates. This method is called for each nodal coordinate,
     * its different dimensions and for each gauss point as well.
     *
     * This method should be implemented by the derrived class.
     *
     * @param rCurrentData                 Data required to calculate partial derivatives
     * @param rShapeFunctions              Shape functions of the gauss point
     * @param rShapeFunctionDerivatives    Shape function derivatives of the gauss point
     * @param rShapeDerivative             Current derivative (i.e. $x^c_k$, where $c$ is the node, $k$ is the dimension)
     * @param detJ_deriv                   Derivative of determinant of jacobian (i.e. $|J|$) w.r.t. $x^c_k$
     * @param rDN_Dx_deriv                 Derivative of shape function gradients w.r.t. $x^c_k$
     * @param rCurrentProcessInfo          Current process info
     * @return double                      Source term shape derivative w.r.t. $x^c_k$
     */
    virtual double CalculateSourceTermShapeSensitivity(
        const TElementData& rCurrentData,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ShapeParameter& rShapeDerivative,
        const double detJ_deriv,
        const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_Dx_deriv,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateSourceTermShapeSensitivity "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    double GetDeltaTime(const ProcessInfo& rProcessInfo) const override
    {
        return -1.0 * rProcessInfo[DELTA_TIME];
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "StabilizedConvectionDiffusionReactionAdjointElement #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StabilizedConvectionDiffusionReactionAdjointElement #"
                 << this->Id();
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        this->pGetGeometry()->PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{
    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{

    void AddPrimalSteadyTermScalarDerivatives(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                              const Variable<double>& rDerivativeVariable,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;

        BoundedVector<double, TNumNodes> effective_kinematic_viscosity_derivatives,
            reaction_derivatives, source_derivatives, tau_derivatives,
            psi_one_derivatives, psi_two_derivatives, chi_derivatives,
            residual_derivatives, positivity_preserving_coeff_derivatives,
            streamline_diffusion_coeff_derivatives,
            crosswind_diffusion_coeff_derivatives, velocity_convective_terms,
            absolute_reaction_tilde_derivatives, scalar_convective_terms, s_derivatives,
            scalar_gradient_norm_derivative, absolute_residual_derivatives;

        array_1d<double, 3> scalar_gradient;

        TElementData current_data;

        for (unsigned int g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            this->CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            this->CalculateElementData(current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            const double scalar_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);

            this->CalculateGradient(scalar_gradient, primal_variable, r_shape_derivatives);
            this->GetConvectionOperator(scalar_convective_terms,
                                        scalar_gradient, r_shape_derivatives);

            const double effective_kinematic_viscosity = this->CalculateEffectiveKinematicViscosity(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityScalarDerivatives(
                effective_kinematic_viscosity_derivatives, rDerivativeVariable, current_data,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction = this->CalculateReactionTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            this->CalculateReactionTermScalarDerivatives(
                reaction_derivatives, rDerivativeVariable, current_data,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double source = this->CalculateSourceTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            this->CalculateSourceTermScalarDerivatives(
                source_derivatives, rDerivativeVariable, current_data,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);
            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStabilizationTauScalarDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction, element_length,
                effective_kinematic_viscosity_derivatives, reaction_derivatives);

            const double s = std::abs(reaction);
            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
                s_derivatives, reaction, reaction_derivatives);

            const double velocity_dot_scalar_gradient =
                inner_prod(velocity, scalar_gradient);

            const double velocity_magnitude = norm_2(velocity);
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double scalar_gradient_norm = norm_2(scalar_gradient);

            const double relaxed_scalar_rate = this->EvaluateInPoint(
                this->GetPrimalRelaxedRateVariable(), gauss_shape_functions);

            double chi{0.0}, k1{0.0}, k2{0.0}, residual{0.0},
                positivity_preserving_coeff{0.0};
            if (scalar_gradient_norm > 0.0 && velocity_magnitude_square > 0.0)
            {
                residual = relaxed_scalar_rate;
                residual += velocity_dot_scalar_gradient;
                residual += reaction * scalar_value;
                residual -= source;
                residual = std::abs(residual);

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);

                positivity_preserving_coeff =
                    residual * chi / (velocity_magnitude_square * scalar_gradient_norm);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateChiScalarDerivatives(
                    chi_derivatives, chi, element_length, bossak_alpha, bossak_gamma,
                    delta_time, reaction, dynamic_tau, reaction_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarGradientScalarDerivative(
                    scalar_gradient_norm_derivative, scalar_gradient, r_shape_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateResidualScalarDerivative(
                    residual_derivatives, scalar_value, reaction, velocity,
                    reaction_derivatives, source_derivatives, gauss_shape_functions,
                    r_shape_derivatives, primal_variable, rDerivativeVariable);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
                    absolute_residual_derivatives, residual, residual_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePositivityPreservationCoefficientScalarDerivatives(
                    positivity_preserving_coeff_derivatives, chi, residual,
                    scalar_gradient_norm, velocity_magnitude_square, chi_derivatives,
                    absolute_residual_derivatives, scalar_gradient_norm_derivative,
                    primal_variable, rDerivativeVariable);

                const double reaction_tilde =
                    reaction + dynamic_tau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
                    absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

                const double psi_one =
                    StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
                        velocity_magnitude, tau, reaction_tilde);
                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiOneScalarDerivatives(
                    psi_one_derivatives, velocity_magnitude, reaction_tilde,
                    tau, tau_derivatives, absolute_reaction_tilde_derivatives);

                const double psi_two =
                    StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
                        reaction_tilde, tau, element_length);
                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiTwoScalarDerivatives(
                    psi_two_derivatives, element_length, tau, reaction_tilde, tau_derivatives,
                    reaction_derivatives, absolute_reaction_tilde_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStreamLineDiffusionCoeffScalarDerivatives(
                    streamline_diffusion_coeff_derivatives, k1, element_length,
                    tau, velocity_magnitude, reaction_tilde, psi_one, psi_two,
                    psi_one_derivatives, psi_two_derivatives, tau_derivatives,
                    reaction_derivatives, effective_kinematic_viscosity_derivatives);

                StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateCrossWindDiffusionCoeffScalarDerivatives(
                    crosswind_diffusion_coeff_derivatives, k2, psi_one, psi_two,
                    element_length, effective_kinematic_viscosity, psi_one_derivatives,
                    psi_two_derivatives, effective_kinematic_viscosity_derivatives);
            }
            else
            {
                crosswind_diffusion_coeff_derivatives.clear();
                streamline_diffusion_coeff_derivatives.clear();
                positivity_preserving_coeff_derivatives.clear();
            }

            // calculating primal damping matrix scalar derivatives
            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double dNa_dNc = 0.0;
                    for (unsigned int i = 0; i < TDim; ++i)
                        dNa_dNc += r_shape_derivatives(a, i) * r_shape_derivatives(c, i);

                    double value = 0.0;

                    // adding derivative of the diffusion term
                    value += scalar_convective_terms[a] *
                             effective_kinematic_viscosity_derivatives[c];

                    // adding reaction term derivatives
                    value += gauss_shape_functions[a] * reaction_derivatives[c] * scalar_value;

                    // adding SUPG stabilization derivatives
                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             velocity_dot_scalar_gradient;
                    value += tau * s_derivatives[c] * gauss_shape_functions[a] *
                             velocity_dot_scalar_gradient;

                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction * scalar_value;
                    value += tau * (s_derivatives[c] * gauss_shape_functions[a]) *
                             reaction * scalar_value;
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction_derivatives[c] * scalar_value;

                    // Adding cross wind dissipation derivatives
                    value += positivity_preserving_coeff_derivatives[c] * k2 *
                             scalar_convective_terms[a] * velocity_magnitude_square;
                    value += positivity_preserving_coeff *
                             crosswind_diffusion_coeff_derivatives[c] *
                             scalar_convective_terms[a] * velocity_magnitude_square;

                    value -= positivity_preserving_coeff_derivatives[c] * k2 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;
                    value -= positivity_preserving_coeff *
                             crosswind_diffusion_coeff_derivatives[c] *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // Adding stream line dissipation derivatives
                    value += positivity_preserving_coeff_derivatives[c] * k1 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;
                    value += positivity_preserving_coeff *
                             streamline_diffusion_coeff_derivatives[c] *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // putting it in the transposed matrix
                    rLeftHandSideMatrix(c, a) += -1.0 * gauss_weights[g] * value;
                }
            }

            // calculating right hand side scalar derivatives
            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double value = 0.0;

                    value += gauss_shape_functions[a] * source_derivatives[c];
                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             source;
                    value += tau * (s_derivatives[c] * gauss_shape_functions[a]) * source;
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             source_derivatives[c];

                    // putting it in the transposed matrix
                    rLeftHandSideMatrix(c, a) += gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AddPrimalSteadyTermScalarRateDerivatives(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                                  const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            scalar_convective_terms, positivity_preserving_coeff_derivatives;
        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;
        array_1d<double, 3> scalar_gradient;

        TElementData current_data;

        for (unsigned int g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            this->CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            this->CalculateElementData(current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            const double scalar_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);

            this->CalculateGradient(scalar_gradient, primal_variable, r_shape_derivatives);
            this->GetConvectionOperator(scalar_convective_terms,
                                        scalar_gradient, r_shape_derivatives);

            const double effective_kinematic_viscosity = this->CalculateEffectiveKinematicViscosity(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction = this->CalculateReactionTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double source = this->CalculateSourceTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            const double velocity_dot_scalar_gradient =
                inner_prod(velocity, scalar_gradient);

            const double velocity_magnitude = norm_2(velocity);
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double scalar_gradient_norm = norm_2(scalar_gradient);

            const double relaxed_scalar_rate = this->EvaluateInPoint(
                this->GetPrimalRelaxedRateVariable(), gauss_shape_functions);

            double chi{0.0}, k1{0.0}, k2{0.0}, residual{0.0};
            if (scalar_gradient_norm > 0.0 && velocity_magnitude_square > 0.0)
            {
                residual = relaxed_scalar_rate;
                residual += velocity_dot_scalar_gradient;
                residual += reaction * scalar_value;
                residual -= source;

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);

                if (residual > 0.0)
                {
                    noalias(positivity_preserving_coeff_derivatives) =
                        gauss_shape_functions * chi /
                        (scalar_gradient_norm * velocity_magnitude_square);
                }
                else
                {
                    noalias(positivity_preserving_coeff_derivatives) =
                        -1.0 * gauss_shape_functions * chi /
                        (scalar_gradient_norm * velocity_magnitude_square);
                }
            }
            else
            {
                positivity_preserving_coeff_derivatives.clear();
            }

            // calculating primal damping matrix scalar derivatives
            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double value = 0.0;

                    // Adding cross wind dissipation derivatives
                    value += positivity_preserving_coeff_derivatives[c] * k2 *
                             scalar_convective_terms[a] * velocity_magnitude_square;

                    value -= positivity_preserving_coeff_derivatives[c] * k2 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // Adding stream line dissipation derivatives
                    value += positivity_preserving_coeff_derivatives[c] * k1 *
                             velocity_convective_terms[a] * velocity_dot_scalar_gradient;

                    // putting it in the transposed matrix
                    rLeftHandSideMatrix(c, a) += -1.0 * gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AddPrimalDampingMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rPrimalDampingMatrix,
                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculatePrimalDampingMatrix(local_matrix, rCurrentProcessInfo);
        local_matrix = trans(local_matrix);
        noalias(rPrimalDampingMatrix) -= local_matrix;

        KRATOS_CATCH("");
    }

    void AddPrimalMassMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rPrimalMassMatrix,
                             const double ScalingFactor,
                             const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculatePrimalMassMatrix(local_matrix, rCurrentProcessInfo);
        local_matrix = trans(local_matrix);
        noalias(rPrimalMassMatrix) += local_matrix * ScalingFactor;

        KRATOS_CATCH("");
    }

    void AddMassTermScalarDerivatives(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                      const Variable<double>& rDerivativeVariable,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            tau_derivatives, effective_kinematic_viscosity_derivatives,
            reaction_derivatives, s_derivatives;

        TElementData current_data;

        for (unsigned int g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            this->CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            this->CalculateElementData(current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            const double effective_kinematic_viscosity = this->CalculateEffectiveKinematicViscosity(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityScalarDerivatives(
                effective_kinematic_viscosity_derivatives, rDerivativeVariable, current_data,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction = this->CalculateReactionTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            this->CalculateReactionTermScalarDerivatives(
                reaction_derivatives, rDerivativeVariable, current_data,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double relaxed_scalar_rate = this->EvaluateInPoint(
                this->GetPrimalRelaxedRateVariable(), gauss_shape_functions);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);
            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStabilizationTauScalarDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction, element_length,
                effective_kinematic_viscosity_derivatives, reaction_derivatives);

            const double s = std::abs(reaction);
            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarValueScalarDerivatives(
                s_derivatives, reaction, reaction_derivatives);

            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double value = 0.0;

                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             relaxed_scalar_rate;
                    value += tau * (s_derivatives[c] * gauss_shape_functions[a]) * relaxed_scalar_rate;

                    rLeftHandSideMatrix(c, a) += -1.0 * gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void CalculateResidualShapeSensitivity(BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
                                           const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rOutput.clear();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        Geometry<NodeType>& r_geometry = this->GetGeometry();

        ShapeParameter deriv;

        Geometry<Point>::JacobiansType J;
        r_geometry.Jacobian(J, this->GetIntegrationMethod());

        Geometry<Point>::ShapeFunctionsGradientsType DN_De;
        DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        const GeometryType::ShapeFunctionsGradientsType& r_dn_de =
            this->GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor,
            parameter_derivatives_shape_derivs, contravariant_metric_tensor_deriv;

        TElementData current_data;

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            primal_variable_gradient_convective_terms,
            convective_primal_variable_gradient_terms_deriv,
            convective_deriv_primal_variable_gradient_terms,
            velocity_convective_terms_deriv,
            primal_variable_gradient_convective_terms_deriv,
            primal_variable_gradient_deriv_convective_terms;

        Vector primal_variable_relaxed_rate_nodal_values(TNumNodes),
            primal_variable_nodal_values(TNumNodes);

        GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_deriv;

        array_1d<double, 3> primal_variable_gradient, primal_variable_gradient_deriv;

        for (unsigned int g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& gauss_r_dn_de = r_dn_de[g];
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);
            const double gauss_weight = gauss_weights[g];

            const Matrix& rJ = J[g];
            const Matrix& rDN_De = DN_De[g];
            const double inv_detJ = 1.0 / MathUtils<double>::DetMat(rJ);
            GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];

            this->CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives_g);

            this->CalculateElementData(current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            const double velocity_magnitude = norm_2(velocity);
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            const double primal_variable_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);
            this->CalculateGradient(primal_variable_gradient, primal_variable,
                                    r_shape_derivatives);
            const double primal_variable_gradient_norm = norm_2(primal_variable_gradient);
            const double velocity_dot_primal_variable_gradient =
                inner_prod(velocity, primal_variable_gradient);
            this->GetConvectionOperator(primal_variable_gradient_convective_terms,
                                        primal_variable_gradient, r_shape_derivatives);

            const double primal_variable_relaxed_rate = this->EvaluateInPoint(
                this->GetPrimalRelaxedRateVariable(), gauss_shape_functions);

            RansVariableUtilities::GetNodalArray(
                primal_variable_relaxed_rate_nodal_values, *this,
                this->GetPrimalRelaxedRateVariable());

            const double effective_kinematic_viscosity = this->CalculateEffectiveKinematicViscosity(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            const double reaction = this->CalculateReactionTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            const double source = this->CalculateSourceTerm(
                current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            double chi{0.0}, stream_line_diffusion_coeff{0.0},
                cross_wind_diffusion_coeff{0.0}, residual{0.0},
                positivity_preserving_coeff{0.0};

            if (primal_variable_gradient_norm > 0.0 && velocity_magnitude_square > 0.0)
            {
                residual = primal_variable_relaxed_rate;
                residual += velocity_dot_primal_variable_gradient;
                residual += reaction * primal_variable_value;
                residual -= source;
                residual = std::abs(residual);

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, stream_line_diffusion_coeff, cross_wind_diffusion_coeff,
                    velocity_magnitude, tau, effective_kinematic_viscosity, reaction,
                    bossak_alpha, bossak_gamma, delta_time, element_length, dynamic_tau);

                positivity_preserving_coeff =
                    chi * residual / (primal_variable_gradient_norm * velocity_magnitude_square);
            }

            RansVariableUtilities::GetNodalArray(primal_variable_nodal_values,
                                                 *this, primal_variable);
            const double psi_one = StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
                velocity_magnitude, tau,
                reaction + dynamic_tau * (1 - bossak_alpha) / (bossak_gamma * delta_time));

            const double psi_two = StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
                reaction + dynamic_tau * (1 - bossak_alpha) / (bossak_gamma * delta_time),
                tau, element_length);

            const double s = std::abs(reaction);

            for (unsigned int c = 0; c < TNumNodes; ++c)
            {
                const unsigned int block_size = c * TDim;
                for (unsigned int k = 0; k < TDim; ++k)
                {
                    double stream_line_diffusion_coeff_deriv{0.0},
                        cross_wind_diffusion_coeff_deriv{0.0},
                        positivity_preserving_coeff_deriv{0.0};

                    deriv.NodeIndex = c;
                    deriv.Direction = k;

                    double detJ_deriv;
                    geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);
                    const double gauss_weight_deriv =
                        detJ_deriv * inv_detJ * gauss_weights[g];

                    RansCalculationUtilities::CalculateGeometryParameterDerivativesShapeSensitivity(
                        parameter_derivatives_shape_derivs, deriv,
                        gauss_r_dn_de, r_parameter_derivatives_g);

                    noalias(contravariant_metric_tensor_deriv) =
                        prod(trans(parameter_derivatives_shape_derivs), r_parameter_derivatives_g) +
                        prod(trans(r_parameter_derivatives_g), parameter_derivatives_shape_derivs);

                    this->CalculateGradient(primal_variable_gradient_deriv,
                                            primal_variable, DN_DX_deriv);
                    const double velocity_dot_primal_variable_gradient_deriv =
                        inner_prod(velocity, primal_variable_gradient_deriv);

                    const double reaction_deriv = this->CalculateReactionTermShapeSensitivity(
                        current_data, gauss_shape_functions, r_shape_derivatives,
                        deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);
                    const double effective_kinematic_viscosity_deriv =
                        this->CalculateEffectiveKinematicViscosityShapeSensitivity(
                            current_data, gauss_shape_functions, r_shape_derivatives,
                            deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);

                    this->GetConvectionOperator(
                        convective_primal_variable_gradient_terms_deriv,
                        primal_variable_gradient_deriv, r_shape_derivatives);
                    this->GetConvectionOperator(convective_deriv_primal_variable_gradient_terms,
                                                primal_variable_gradient, DN_DX_deriv);

                    const double element_length_deriv =
                        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateElementLengthH2ShapeSensitivity(
                            velocity_magnitude, velocity, contravariant_metric_tensor,
                            contravariant_metric_tensor_deriv);

                    const double tau_deriv =
                        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStabilizationTauShapeSensitivity(
                            tau, velocity_magnitude, element_length,
                            element_length_deriv, effective_kinematic_viscosity,
                            effective_kinematic_viscosity_deriv, reaction, reaction_deriv);

                    const double chi_deriv =
                        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateChiShapeSensitivity(
                            chi, reaction, reaction_deriv, element_length, element_length_deriv,
                            bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

                    const double primal_variable_gradient_norm_deriv =
                        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateAbsoluteScalarGradientShapeSensitivity(
                            primal_variable_gradient, DN_DX_deriv, primal_variable_nodal_values);

                    const double source_deriv = this->CalculateSourceTermShapeSensitivity(
                        current_data, gauss_shape_functions, r_shape_derivatives,
                        deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);

                    const double residual_deriv =
                        StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateResidualShapeSensitivity(
                            residual, velocity, DN_DX_deriv, primal_variable_value,
                            primal_variable_nodal_values, reaction_deriv, source_deriv);

                    if (primal_variable_gradient_norm > 0.0 && velocity_magnitude_square > 0.0)
                    {
                        positivity_preserving_coeff_deriv =
                            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePositivityPreservationCoefficientShapeSensitivity(
                                chi, chi_deriv, residual, residual_deriv,
                                velocity_magnitude_square, primal_variable_gradient_norm,
                                primal_variable_gradient_norm_deriv);

                        const double psi_one_deriv =
                            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiOneShapeSensitivity(
                                tau, tau_deriv, velocity_magnitude, reaction, reaction_deriv,
                                bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

                        const double psi_two_deriv =
                            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculatePsiTwoShapeSensitivity(
                                psi_two, element_length, element_length_deriv,
                                reaction, reaction_deriv, tau, tau_deriv,
                                bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

                        stream_line_diffusion_coeff_deriv =
                            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateStreamLineDiffusionCoeffShapeSensitivity(
                                stream_line_diffusion_coeff, psi_one, psi_one_deriv,
                                tau, tau_deriv, velocity_magnitude, reaction,
                                reaction_deriv, element_length, element_length_deriv,
                                effective_kinematic_viscosity_deriv, psi_two_deriv,
                                bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

                        cross_wind_diffusion_coeff_deriv =
                            StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateCrossWindDiffusionCoeffShapeSensitivity(
                                cross_wind_diffusion_coeff, psi_one, psi_one_deriv,
                                element_length, element_length_deriv,
                                effective_kinematic_viscosity_deriv, psi_two_deriv);
                    }

                    double s_deriv = 0.0;
                    if (reaction > 0.0)
                    {
                        s_deriv = reaction_deriv;
                    }
                    else
                    {
                        s_deriv = -reaction_deriv;
                    }

                    this->GetConvectionOperator(velocity_convective_terms_deriv,
                                                velocity, DN_DX_deriv);
                    this->GetConvectionOperator(primal_variable_gradient_convective_terms_deriv,
                                                primal_variable_gradient, DN_DX_deriv);
                    this->GetConvectionOperator(
                        primal_variable_gradient_deriv_convective_terms,
                        primal_variable_gradient_deriv, r_shape_derivatives);

                    for (unsigned int a = 0; a < TNumNodes; ++a)
                    {
                        const double tau_operator = velocity_convective_terms[a] +
                                                    s * gauss_shape_functions[a];
                        const double tau_operator_deriv =
                            velocity_convective_terms_deriv[a] +
                            s_deriv * gauss_shape_functions[a];

                        double value = 0.0;

                        // adding convective term shape sensitivity
                        value += gauss_shape_functions[a] *
                                 velocity_dot_primal_variable_gradient_deriv * gauss_weight;
                        value += gauss_shape_functions[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight_deriv;

                        // adding reaction term shape sensitivity
                        value += gauss_shape_functions[a] * primal_variable_value *
                                 reaction_deriv * gauss_weight;
                        value += gauss_shape_functions[a] * primal_variable_value *
                                 reaction * gauss_weight_deriv;

                        // adding diffusion term shape sensitivity
                        value += effective_kinematic_viscosity_deriv *
                                 primal_variable_gradient_convective_terms[a] * gauss_weight;
                        value += effective_kinematic_viscosity *
                                 convective_deriv_primal_variable_gradient_terms[a] *
                                 gauss_weight;
                        value += effective_kinematic_viscosity *
                                 convective_primal_variable_gradient_terms_deriv[a] *
                                 gauss_weight;
                        value += effective_kinematic_viscosity *
                                 primal_variable_gradient_convective_terms[a] *
                                 gauss_weight_deriv;

                        // adding SUPG term derivatives
                        value += tau_deriv * tau_operator *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += tau * tau_operator_deriv *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += tau * tau_operator *
                                 velocity_dot_primal_variable_gradient_deriv * gauss_weight;
                        value += tau * tau_operator *
                                 velocity_dot_primal_variable_gradient * gauss_weight_deriv;

                        value += tau_deriv * tau_operator * reaction *
                                 primal_variable_value * gauss_weight;
                        value += tau * tau_operator_deriv * reaction *
                                 primal_variable_value * gauss_weight;
                        value += tau * tau_operator * reaction_deriv *
                                 primal_variable_value * gauss_weight;
                        value += tau * tau_operator * reaction *
                                 primal_variable_value * gauss_weight_deriv;

                        // adding cross wind dissipation term derivatives
                        value += velocity_magnitude_square * positivity_preserving_coeff_deriv *
                                 cross_wind_diffusion_coeff *
                                 primal_variable_gradient_convective_terms[a] * gauss_weight;
                        value += velocity_magnitude_square * positivity_preserving_coeff *
                                 cross_wind_diffusion_coeff_deriv *
                                 primal_variable_gradient_convective_terms[a] * gauss_weight;
                        value += velocity_magnitude_square * positivity_preserving_coeff *
                                 cross_wind_diffusion_coeff *
                                 primal_variable_gradient_convective_terms_deriv[a] *
                                 gauss_weight;
                        value += velocity_magnitude_square * positivity_preserving_coeff *
                                 cross_wind_diffusion_coeff *
                                 primal_variable_gradient_deriv_convective_terms[a] *
                                 gauss_weight;
                        value += velocity_magnitude_square * positivity_preserving_coeff *
                                 cross_wind_diffusion_coeff *
                                 primal_variable_gradient_convective_terms[a] *
                                 gauss_weight_deriv;

                        value -= positivity_preserving_coeff_deriv * cross_wind_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value -= positivity_preserving_coeff * cross_wind_diffusion_coeff_deriv *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value -= positivity_preserving_coeff * cross_wind_diffusion_coeff *
                                 velocity_convective_terms_deriv[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value -= positivity_preserving_coeff * cross_wind_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient_deriv * gauss_weight;
                        value -= positivity_preserving_coeff * cross_wind_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight_deriv;

                        // adding stream line diffusion term derivatives
                        value += positivity_preserving_coeff_deriv * stream_line_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += positivity_preserving_coeff * stream_line_diffusion_coeff_deriv *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += positivity_preserving_coeff * stream_line_diffusion_coeff *
                                 velocity_convective_terms_deriv[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight;
                        value += positivity_preserving_coeff * stream_line_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient_deriv * gauss_weight;
                        value += positivity_preserving_coeff * stream_line_diffusion_coeff *
                                 velocity_convective_terms[a] *
                                 velocity_dot_primal_variable_gradient * gauss_weight_deriv;

                        // adding right hand side terms shape gradient terms

                        value -= gauss_shape_functions[a] * source_deriv * gauss_weight;
                        value -= gauss_shape_functions[a] * source * gauss_weight_deriv;

                        value -= tau_deriv * tau_operator * source * gauss_weight;
                        value -= tau * tau_operator_deriv * source * gauss_weight;
                        value -= tau * tau_operator * source_deriv * gauss_weight;
                        value -= tau * tau_operator * source * gauss_weight_deriv;

                        // adding mass shape gradient terms
                        value += gauss_weight_deriv *
                                 primal_variable_relaxed_rate_nodal_values[a] / TNumNodes;
                        value += tau_deriv * tau_operator *
                                 primal_variable_relaxed_rate * gauss_weight;
                        value += tau * tau_operator_deriv *
                                 primal_variable_relaxed_rate * gauss_weight;
                        value += tau * tau_operator * primal_variable_relaxed_rate * gauss_weight_deriv;

                        rOutput(block_size + k, a) -= value;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_TRY

        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Private  Access
    ///@{
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class Element

///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes, class TElementData>
inline std::istream& operator>>(
    std::istream& rIStream,
    StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, TElementData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TElementData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

} // namespace Kratos.
#endif // KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT defined
