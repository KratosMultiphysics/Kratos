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

#if !defined(KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_H_INCLUDED)
#define KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/geometry_utilities.h"
#include "utilities/time_discretization.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
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

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
class StabilizedConvectionDiffusionReaction : public Element
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = Element;

    /// Node type (default is: Node<3>)
    using NodeType = Node<3>;

    /// Geometry type (using with given NodeType)
    using GeometryType = Geometry<NodeType>;

    /// Definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    /// Vector type for local contributions to the linear system
    using VectorType = Vector;

    /// Matrix type for local contributions to the linear system
    using MatrixType = Matrix;

    using IndexType = std::size_t;

    using EquationIdVectorType = std::vector<IndexType>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    /// Type for an array of shape function gradient matrices
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    using ConvectionDiffusionReactionDataType = TConvectionDiffusionReactionData;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of StabilizedConvectionDiffusionReaction
    KRATOS_CLASS_POINTER_DEFINITION(StabilizedConvectionDiffusionReaction);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit StabilizedConvectionDiffusionReaction(IndexType NewId = 0)
        : Element(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    StabilizedConvectionDiffusionReaction(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    StabilizedConvectionDiffusionReaction(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    StabilizedConvectionDiffusionReaction(IndexType NewId,
                                          GeometryType::Pointer pGeometry,
                                          PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    StabilizedConvectionDiffusionReaction(StabilizedConvectionDiffusionReaction const& rOther)
        : Element(rOther)
    {
    }

    /**
     * Destructor
     */
    ~StabilizedConvectionDiffusionReaction() override = default;

    ///@}
    ///@name Operators
    ///@{

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
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to Create base "
                        "StabilizedConvectionDiffusionReaction instances."
                     << std::endl;
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
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to Create base "
                        "StabilizedConvectionDiffusionReaction instances."
                     << std::endl;
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
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to Clone base "
                        "StabilizedConvectionDiffusionReaction instances."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReaction "
                        "EquationIdVector method. Please implement it in the "
                        "derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReaction GetDofList "
                        "method. Please implement it in the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide
     * methods they can be managed internally with a private method to do the
     * same calculations only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override
    {
        // Check sizes and initialize matrix
        if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes)
            rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

        // Calculate RHS
        this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rRightHandSideVector.size() != TNumNodes)
            rRightHandSideVector.resize(TNumNodes, false);

        noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const IndexType num_gauss_points = gauss_weights.size();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector gauss_shape_functions = row(shape_functions, g);

            this->CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            const array_1d<double, 3> velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);

            TConvectionDiffusionReactionData r_current_data;
            this->CalculateElementData(r_current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);
            const double effective_kinematic_viscosity = this->CalculateEffectiveKinematicViscosity(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction = this->CalculateReactionTerm(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            const double source =
                this->CalculateSourceTerm(r_current_data, gauss_shape_functions,
                                          r_shape_derivatives, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            const double s = std::abs(reaction);

            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                double value = 0.0;

                value += gauss_shape_functions[a] * source;

                // Add supg stabilization terms
                value += (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                         tau * source;

                rRightHandSideVector[a] += gauss_weights[g] * value;
            }
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief CalculateLocalVelocityContribution Calculate the local contribution in terms of velocity and pressure.
     * @param rDampMatrix Local finite element system matrix (output)
     * @param rRightHandSideVector Local finite element residual vector (output)
     * @param rCurrentProcessInfo Current ProcessInfo values (input)
     */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
                                            VectorType& rRightHandSideVector,
                                            ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);

        // Now calculate an additional contribution to the residual: r -= rDampingMatrix * (u,p)
        VectorType U;
        this->GetValuesVector(U);
        noalias(rRightHandSideVector) -= prod(rDampingMatrix, U);
    }

    /**
     * ELEMENTS inherited from this class must implement this methods
     * if they need to add dynamic element contributions
     * note: second derivatives means the accelerations if the displacements are the dof of the analysis
     * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
     * CalculateSecondDerivativesContributions,
     * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
     */

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculatePrimalMassMatrix(local_matrix, rCurrentProcessInfo);

        if (rMassMatrix.size1() != TNumNodes || rMassMatrix.size2() != TNumNodes)
            rMassMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rMassMatrix) = local_matrix;
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental damping matrix
     * @param rDampingMatrix: the elemental damping matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculatePrimalDampingMatrix(local_matrix, rCurrentProcessInfo);

        if (rDampingMatrix.size1() != TNumNodes || rDampingMatrix.size2() != TNumNodes)
            rDampingMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rDampingMatrix) = local_matrix;
    }

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int check = BaseType::Check(rCurrentProcessInfo);

        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, this->GetGeometry()[iNode]);
        }

        return check;

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates scalar value for given gauss point
     *
     * @param rVariable      Scalar variable
     * @param rShapeFunction Gauss point shape functions
     * @param Step           Step
     * @return double        Gauss point scalar value
     */
    double EvaluateInPoint(const Variable<double>& rVariable,
                           const Vector& rShapeFunction,
                           const int Step = 0) const
    {
        return RansCalculationUtilities::EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
    }

    /**
     * @brief Calculates vector value for given gauss point
     *
     * @param rVariable            Vector variable
     * @param rShapeFunction       Gauss point shape functions
     * @param Step                 Step
     * @return array_1d<double, 3> Gauss point vector value
     */
    array_1d<double, 3> EvaluateInPoint(const Variable<array_1d<double, 3>>& rVariable,
                                        const Vector& rShapeFunction,
                                        const int Step = 0) const
    {
        return RansCalculationUtilities::EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
    }

    /**
     * @brief Get the Divergence Operator object
     *
     * Calculates divergence of a vector at a gauss point
     *
     * @param rVariable          Vector variable
     * @param rShapeDerivatives  Shape derivatives at gauss point
     * @param Step               time step
     * @return double            Divergence of the variable
     */
    double GetDivergenceOperator(const Variable<array_1d<double, 3>>& rVariable,
                                 const Matrix& rShapeDerivatives,
                                 const int Step = 0) const
    {
        double value = 0.0;
        const GeometryType& r_geometry = this->GetGeometry();

        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            const array_1d<double, 3>& r_value =
                r_geometry[i].FastGetSolutionStepValue(rVariable, Step);
            for (IndexType j = 0; j < TDim; ++j)
            {
                value += r_value[j] * rShapeDerivatives(i, j);
            }
        }

        return value;
    }

    /**
     * @brief Calculates all the data required by the element.
     *
     * This method is used to calculate and store all the required
     * quantities under each gauss point. This method is called
     * for each gauss point, before calculating the derivatives
     *
     * This method should be implemented by the derrived class.
     *
     * @param rData                      Element data container
     * @param rShapeFunctions            Gauss point shape functions
     * @param rShapeFunctionDerivatives  Gauss point shape function derivatives
     * @param rCurrentProcessInfo        Current process info
     */
    virtual void CalculateElementData(TConvectionDiffusionReactionData& rData,
                                      const Vector& rShapeFunctions,
                                      const Matrix& rShapeFunctionDerivatives,
                                      const ProcessInfo& rCurrentProcessInfo,
                                      const int Step = 0) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReaction "
                        "CalculateElementData method. "
                        "Please implement it in the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * @brief Calculate effective kinematic viscosity
     *
     * Calculate effective kinematic viscosity (i.e. $\nu_\phi$)
     * This method is called for each gauss point.
     * This method should be implemented by the derrived class.
     *
     * @param rData
     * @param rShapeFunctions
     * @param rShapeFunctionDerivatives
     * @param rCurrentProcessInfo
     * @param Step
     * @return double
     */
    virtual double CalculateEffectiveKinematicViscosity(const TConvectionDiffusionReactionData& rData,
                                                        const Vector& rShapeFunctions,
                                                        const Matrix& rShapeFunctionDerivatives,
                                                        const ProcessInfo& rCurrentProcessInfo,
                                                        const int Step = 0) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReaction "
                        "CalculateEffectiveKinematicViscosity method. "
                        "Please implement it in the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * @brief Get the primal scalar variable
     *
     * This returns the scalar variable ($\phi$) used in stabilized
     * convection-diffusion-reaction transport equation.
     *
     * This method should be implemented by the derrived class.
     *
     * @return const Variable<double>&
     */
    virtual const Variable<double>& GetPrimalVariable() const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base GetPrimalVariable method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        return PRESSURE;

        KRATOS_CATCH("");
    }

    /**
     * @brief Get the primal relaxed rate variable
     *
     * This method returns the relaxed scalar rate variable ($\dot{\phi}_r$) calculated according
     * to following equation, where $n$ is the time step:
     *
     * \[
     *      \dot{\phi}_r = \left(1-\alpha_{bossak}\right)\dot{\phi}^{n} + \alpha_{bossak}\dot{\phi}^{n-1}
     * \]
     *
     * This method should be implemented by the derrived class.
     *
     * @return const Variable<double>&
     */
    virtual const Variable<double>& GetPrimalRelaxedRateVariable() const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base GetPrimalRelaxedRateVariable method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        return PRESSURE;

        KRATOS_CATCH("");
    }

    double GetScalarVariableGradientNorm(const Matrix& rShapeFunctionDerivatives,
                                         const int Step = 0) const
    {
        KRATOS_TRY;

        array_1d<double, 3> scalar_variable_gradient;
        this->CalculateGradient(scalar_variable_gradient, this->GetPrimalVariable(),
                                rShapeFunctionDerivatives, Step);
        return norm_2(scalar_variable_gradient);

        KRATOS_CATCH("");
    }

    double GetScalarVariableRelaxedAcceleration(const Vector& rShapeFunctions,
                                                const int Step = 0) const
    {
        KRATOS_TRY;

        return this->EvaluateInPoint(this->GetPrimalRelaxedRateVariable(),
                                     rShapeFunctions, Step);

        KRATOS_CATCH("");
    }

    /**
     * @brief Get the Geometry Parameter Derivatives object
     *
     * This method calculates partial derivatives of parametric coordinates(i.e. $\underline{\xi}$) of element
     * w.r.t. physical coordinates (i.e. $\underline{x}$)
     *
     * \[
     *      \frac{\partial \underline{\xi}}{\partial \underline{x}}
     * \]
     *
     * @return ShapeFunctionDerivativesArrayType
     */
    ShapeFunctionDerivativesArrayType GetGeometryParameterDerivatives() const
    {
        const GeometryType& r_geometry = this->GetGeometry();
        return RansCalculationUtilities::CalculateGeometryParameterDerivatives(
            r_geometry, this->GetIntegrationMethod());
    }

    /**
     * @brief Calculates reaction coefficient
     *
     * This method calculates reaction coefficient (i.e. $s_\phi$).
     * This method is called for each gauss point.
     * This method should be implemented by the derrived class.
     *
     * @param rData
     * @param rShapeFunctions
     * @param rShapeFunctionDerivatives
     * @param rCurrentProcessInfo
     * @param Step
     * @return double
     */
    virtual double CalculateReactionTerm(const TConvectionDiffusionReactionData& rData,
                                         const Vector& rShapeFunctions,
                                         const Matrix& rShapeFunctionDerivatives,
                                         const ProcessInfo& rCurrentProcessInfo,
                                         const int Step = 0) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReaction "
                        "CalculateReactionTerm method. Please implement it in "
                        "the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates source term
     *
     * This method calculates source term (i.e. $f_\phi$).
     * This method is called for each gauss point.
     * This method should be implemented by the derrived class.
     *
     * @param rData
     * @param rShapeFunctions
     * @param rShapeFunctionDerivatives
     * @param rCurrentProcessInfo
     * @param Step
     * @return double
     */
    virtual double CalculateSourceTerm(const TConvectionDiffusionReactionData& rData,
                                       const Vector& rShapeFunctions,
                                       const Matrix& rShapeFunctionDerivatives,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       const int Step = 0) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReaction "
                        "CalculateSourceTerm method. Please implement it in "
                        "the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    void CalculatePrimalDampingMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rDampingMatrix,
                                      const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        rDampingMatrix.clear();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const IndexType num_gauss_points = gauss_weights.size();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const double eps = std::numeric_limits<double>::epsilon();

        array_1d<double, 3> variable_gradient;
        const Variable<double>& primal_variable = this->GetPrimalVariable();

        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector gauss_shape_functions = row(shape_functions, g);

            this->CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);
            const double velocity_magnitude = norm_2(velocity);

            TConvectionDiffusionReactionData r_current_data;
            this->CalculateElementData(r_current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);
            const double effective_kinematic_viscosity = this->CalculateEffectiveKinematicViscosity(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);
            const double variable_gradient_norm =
                this->GetScalarVariableGradientNorm(r_shape_derivatives);
            const double relaxed_variable_acceleration =
                this->GetScalarVariableRelaxedAcceleration(gauss_shape_functions);
            this->CalculateGradient(variable_gradient, primal_variable, r_shape_derivatives);

            const double reaction = this->CalculateReactionTerm(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            // Calculate residual for cross wind dissipation coefficient
            double positivity_preserving_coefficient{0.0}, k1{0.0}, k2{0.0}, chi{0.0};
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

            const double velocity_dot_variable_gradient =
                inner_prod(velocity, variable_gradient);
            const double variable_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);

            if (variable_gradient_norm > eps && velocity_magnitude_square > eps)
            {
                const double source = this->CalculateSourceTerm(
                    r_current_data, gauss_shape_functions, r_shape_derivatives,
                    rCurrentProcessInfo);

                double residual = relaxed_variable_acceleration;
                residual += velocity_dot_variable_gradient;
                residual += reaction * variable_value;
                residual -= source;
                residual = std::abs(residual);

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);

                positivity_preserving_coefficient = residual * chi / (variable_gradient_norm * velocity_magnitude_square);
            }

            const double s = std::abs(reaction);

            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                for (IndexType b = 0; b < TNumNodes; ++b)
                {
                    double dNa_dNb = 0.0;
                    for (IndexType i = 0; i < TDim; ++i)
                        dNa_dNb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);

                    double value = 0.0;

                    value += gauss_shape_functions[a] * velocity_convective_terms[b];
                    value += gauss_shape_functions[a] * reaction *
                             gauss_shape_functions[b];
                    value += effective_kinematic_viscosity * dNa_dNb;

                    // Adding SUPG stabilization terms
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             velocity_convective_terms[b];
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction * gauss_shape_functions[b];

                    // Adding cross wind dissipation
                    value += positivity_preserving_coefficient * k2 * dNa_dNb *
                             velocity_magnitude_square;
                    value -= positivity_preserving_coefficient * k2 *
                             velocity_convective_terms[a] *
                             velocity_convective_terms[b];

                    // Adding stream line dissipation
                    value += positivity_preserving_coefficient * k1 *
                             velocity_convective_terms[a] *
                             velocity_convective_terms[b];

                    rDampingMatrix(a, b) += gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void CalculatePrimalMassMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix,
                                   const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rMassMatrix.clear();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const IndexType num_gauss_points = gauss_weights.size();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        BoundedMatrix<double, TDim, TDim> contravariant_metric_tensor;

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector gauss_shape_functions = row(shape_functions, g);

            this->CalculateContravariantMetricTensor(
                contravariant_metric_tensor, r_parameter_derivatives[g]);

            const double mass = gauss_weights[g] * (1.0 / TNumNodes);
            this->AddLumpedMassMatrix(rMassMatrix, mass);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            TConvectionDiffusionReactionData r_current_data;
            this->CalculateElementData(r_current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);
            const double effective_kinematic_viscosity = this->CalculateEffectiveKinematicViscosity(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction = this->CalculateReactionTerm(
                r_current_data, gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            const double s = std::abs(reaction);

            // Add mass stabilization terms
            for (IndexType i = 0; i < TNumNodes; ++i)
                for (IndexType j = 0; j < TNumNodes; ++j)
                    rMassMatrix(i, j) +=
                        gauss_weights[g] * tau *
                        (velocity_convective_terms[i] + s * gauss_shape_functions[i]) *
                        gauss_shape_functions[j];
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Get the Convection Operator object
     *
     * Calculates convection operator given by following equation
     *
     * \[
     *  w_i\frac{\partial N^a}{\partial x_i}
     * \]
     *
     * $w_i$ being the $i^{th}$ dimension of $\underline{w}$ vector, $N^a$ being the
     * shape function of $a^{th}$ node, $x_i$ being the $i^{th}$ dimension
     * of local coordinates
     *
     * @param rOutput           Vector of results
     * @param rVector           Input vector (i.e. $\underline{w}$)
     * @param rShapeDerivatives Shape function derivatives w.r.t. physical coordinates
     */
    void GetConvectionOperator(BoundedVector<double, TNumNodes>& rOutput,
                               const array_1d<double, 3>& rVector,
                               const Matrix& rShapeDerivatives) const
    {
        rOutput.clear();
        for (IndexType i = 0; i < TNumNodes; ++i)
            for (IndexType j = 0; j < TDim; ++j)
            {
                rOutput[i] += rVector[j] * rShapeDerivatives(i, j);
            }
    }

    /**
     * @brief Calculate gradient matrix for a vector
     *
     * Calculates the gradient matrix for a given vector variable.
     *
     * @param rOutput            Output matrix, rows contain the given vector indices, columns containt physical coordinate dimensions
     * @param rVariable          Vector variable
     * @param rShapeDerivatives  Shape function derivatives at the gauss point
     * @param Step               Time step
     */
    void CalculateGradient(BoundedMatrix<double, TDim, TDim>& rOutput,
                           const Variable<array_1d<double, 3>>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const
    {
        const GeometryType& r_geometry = this->GetGeometry();

        RansCalculationUtilities::CalculateGradient<TDim>(
            rOutput, r_geometry, rVariable, rShapeDerivatives, Step);
    }

    /**
     * @brief Calculate gradient vector for a scalar
     *
     * Calculates the gradient vector for a given scalar variable.
     *
     * @param rOutput            Output vector
     * @param rVariable          Scalar variable
     * @param rShapeDerivatives  Shape function derivatives at the gauss point
     * @param Step               Time step
     */
    void CalculateGradient(array_1d<double, 3>& rOutput,
                           const Variable<double>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const
    {
        const GeometryType& r_geometry = this->GetGeometry();
        RansCalculationUtilities::CalculateGradient(
            rOutput, r_geometry, rVariable, rShapeDerivatives, Step);
    }

    void CalculateContravariantMetricTensor(BoundedMatrix<double, TDim, TDim>& rOutput,
                                            const Matrix& rParameterDerivatives) const
    {
        noalias(rOutput) = prod(trans(rParameterDerivatives), rParameterDerivatives);
    }

    virtual double GetDeltaTime(const ProcessInfo& rProcessInfo) const
    {
        return rProcessInfo[DELTA_TIME];
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
        buffer << "StabilizedConvectionDiffusionReaction #" << Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StabilizedConvectionDiffusionReaction #" << Id();
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
    /**
     * @brief Calculates shape function data for this element
     *
     * @param rGaussWeights Gauss point weights list
     * @param rNContainer   Shape function values. Each row contains shape functions for respective gauss point
     * @param rDN_DX        List of matrices containing shape function derivatives for each gauss point
     */
    virtual void CalculateGeometryData(Vector& rGaussWeights,
                                       Matrix& rNContainer,
                                       ShapeFunctionDerivativesArrayType& rDN_DX) const
    {
        const GeometryType& r_geometry = this->GetGeometry();

        RansCalculationUtilities::CalculateGeometryData(
            r_geometry, this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
    }

    void AddLumpedMassMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rMassMatrix,
                             const double Mass) const
    {
        for (IndexType iNode = 0; iNode < TNumNodes; ++iNode)
            rMassMatrix(iNode, iNode) += Mass;
    }

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

}; // Class StabilizedConvectionDiffusionReaction

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::istream& operator>>(
    std::istream& rIStream,
    StabilizedConvectionDiffusionReaction<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const StabilizedConvectionDiffusionReaction<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_H_INCLUDED defined
