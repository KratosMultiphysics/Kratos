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

#if !defined(KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED)
#define KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "includes/checks.h"
#include "includes/element.h"
#include "rans_modelling_application_variables.h"
#include "stabilized_convection_diffusion_reaction_utilities.h"
#include "utilities/time_discretization.h"
#include "includes/cfd_variables.h"
#include "utilities/geometry_utilities.h"
#include "containers/array_1d.h"

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
class StabilizedConvectionDiffusionReactionElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    /// Node type (default is: Node<3>)
    typedef Node<3> NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Type for shape function values container
    typedef MatrixRow<Matrix> ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    typedef TConvectionDiffusionReactionData ConvectionDiffusionReactionDataType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of StabilizedConvectionDiffusionReactionElement
    KRATOS_CLASS_POINTER_DEFINITION(StabilizedConvectionDiffusionReactionElement);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    StabilizedConvectionDiffusionReactionElement(IndexType NewId = 0)
        : Element(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    StabilizedConvectionDiffusionReactionElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    StabilizedConvectionDiffusionReactionElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    StabilizedConvectionDiffusionReactionElement(IndexType NewId,
                                                 GeometryType::Pointer pGeometry,
                                                 PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    StabilizedConvectionDiffusionReactionElement(StabilizedConvectionDiffusionReactionElement const& rOther)
        : Element(rOther)
    {
    }

    /**
     * Destructor
     */
    ~StabilizedConvectionDiffusionReactionElement() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    StabilizedConvectionDiffusionReactionElement& operator=(
        StabilizedConvectionDiffusionReactionElement const& rOther)
    {
        BaseType::operator=(rOther);
        Flags::operator=(rOther);
        // mpProperties = rOther.mpProperties;
        return *this;
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
    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override
    {
        KRATOS_TRY;
        KRATOS_ERROR
            << "Attempting to Create base "
               "StabilizedConvectionDiffusionReactionElement instances."
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
        KRATOS_ERROR
            << "Attempting to Create base "
               "StabilizedConvectionDiffusionReactionElement instances."
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
        KRATOS_TRY
        return Kratos::make_intrusive<StabilizedConvectionDiffusionReactionElement>(
            NewId, GetGeometry().Create(ThisNodes), pGetProperties());
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
                        "StabilizedConvectionDiffusionReactionElement "
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
        KRATOS_ERROR
            << "Attempting to call base "
               "StabilizedConvectionDiffusionReactionElement GetDofList "
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
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override
    {
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
        const unsigned int num_gauss_points = gauss_weights.size();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                               r_parameter_derivatives_g.size2());
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const array_1d<double, 3> velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);

            TConvectionDiffusionReactionData r_current_data;
            double effective_kinematic_viscosity;
            this->CalculateConvectionDiffusionReactionData(
                r_current_data, effective_kinematic_viscosity,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(r_current_data, rCurrentProcessInfo);
            const double source =
                this->CalculateSourceTerm(r_current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            const double s = std::abs(reaction);

            for (unsigned int a = 0; a < TNumNodes; ++a)
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
     * this is called during the assembling process in order
     * to calculate the first derivatives contributions for the LHS and RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                VectorType& rRightHandSideVector,
                                                ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != 0)
            rLeftHandSideMatrix.resize(0, 0, false);
        if (rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0, false);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix for the first derivatives constributions
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != 0)
            rLeftHandSideMatrix.resize(0, 0, false);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector for the first derivatives constributions
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        if (rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0, false);
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
     * to calculate the second derivative contributions for the LHS and RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
                                                 VectorType& rRightHandSideVector,
                                                 ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != 0)
            rLeftHandSideMatrix.resize(0, 0, false);
        if (rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0, false);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix for the second derivatives constributions
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        if (rLeftHandSideMatrix.size1() != 0)
            rLeftHandSideMatrix.resize(0, 0, false);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector for the second derivatives constributions
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        if (rRightHandSideVector.size() != 0)
            rRightHandSideVector.resize(0, false);
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rMassMatrix.size1() != TNumNodes || rMassMatrix.size2() != TNumNodes)
            rMassMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rMassMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const unsigned int num_gauss_points = gauss_weights.size();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                               r_parameter_derivatives_g.size2());
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const double mass = gauss_weights[g] / TNumNodes;
            this->AddLumpedMassMatrix(rMassMatrix, mass);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            TConvectionDiffusionReactionData r_current_data;
            double effective_kinematic_viscosity;
            this->CalculateConvectionDiffusionReactionData(
                r_current_data, effective_kinematic_viscosity,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(r_current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor, reaction,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

            const double s = std::abs(reaction);

            // Add mass stabilization terms
            for (unsigned int i = 0; i < TNumNodes; ++i)
                for (unsigned int j = 0; j < TNumNodes; ++j)
                    rMassMatrix(i, j) +=
                        gauss_weights[g] * tau *
                        (velocity_convective_terms[i] + s * gauss_shape_functions[i]) *
                        gauss_shape_functions[j];
        }

        KRATOS_CATCH("");
    }

    /**
     * this is called during the assembling process in order
     * to calculate the elemental damping matrix
     * @param rDampingMatrix: the elemental damping matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rDampingMatrix.size1() != TNumNodes || rDampingMatrix.size2() != TNumNodes)
            rDampingMatrix.resize(TNumNodes, TNumNodes, false);

        noalias(rDampingMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const unsigned int num_gauss_points = gauss_weights.size();

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            Matrix contravariant_metric_tensor(r_parameter_derivatives_g.size1(),
                                               r_parameter_derivatives_g.size2());
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const array_1d<double, 3> velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);
            const double velocity_magnitude = norm_2(velocity);

            TConvectionDiffusionReactionData r_current_data;
            double effective_kinematic_viscosity, variable_gradient_norm,
                relaxed_variable_acceleration;
            this->CalculateConvectionDiffusionReactionData(
                r_current_data, effective_kinematic_viscosity,
                variable_gradient_norm, relaxed_variable_acceleration,
                gauss_shape_functions, r_shape_derivatives, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(r_current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            // Calculate residual for cross wind dissipation coefficient
            double cross_wind_diffusion{0.0}, stream_line_diffusion{0.0};
            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

            Vector nodal_variable;
            this->GetValuesVector(nodal_variable);

            if (variable_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                const double source =
                    this->CalculateSourceTerm(r_current_data, rCurrentProcessInfo);

                double residual = relaxed_variable_acceleration;
                residual += inner_prod(velocity_convective_terms, nodal_variable);
                residual += reaction * inner_prod(gauss_shape_functions, nodal_variable);
                residual -= source;
                residual = std::abs(residual);
                residual /= variable_gradient_norm;

                double chi, k1, k2;
                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau, effective_kinematic_viscosity,
                    reaction, bossak_alpha, bossak_gamma, delta_time, element_length, dynamic_tau);

                stream_line_diffusion = residual * chi * k1 / velocity_magnitude_square;
                cross_wind_diffusion = residual * chi * k2 / velocity_magnitude_square;
            }

            const double s = std::abs(reaction);

            for (unsigned int a = 0; a < TNumNodes; a++)
            {
                for (unsigned int b = 0; b < TNumNodes; b++)
                {
                    double dNa_dNb = 0.0;
                    for (unsigned int i = 0; i < TDim; i++)
                        dNa_dNb += r_shape_derivatives(a, i) * r_shape_derivatives(b, i);

                    double value = 0.0;

                    value += gauss_shape_functions[a] * velocity_convective_terms[b];
                    value += gauss_shape_functions[a] * reaction *
                             gauss_shape_functions[b]; // * positive_values_list[b];
                    value += effective_kinematic_viscosity * dNa_dNb;

                    // Adding SUPG stabilization terms
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             velocity_convective_terms[b];
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction * gauss_shape_functions[b]; // * positive_values_list[b];

                    // Adding cross wind dissipation
                    value += cross_wind_diffusion * dNa_dNb * velocity_magnitude_square;
                    value -= cross_wind_diffusion * velocity_convective_terms[a] *
                             velocity_convective_terms[b];

                    // Adding stream line dissipation
                    value += stream_line_diffusion * velocity_convective_terms[a] *
                             velocity_convective_terms[b];

                    rDampingMatrix(a, b) += gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    /// Implementation of Calculate to compute an error estimate.
    /**
     * If rVariable == ERROR_RATIO, this function will provide an a posteriori
     * estimate of the norm of the subscale velocity, calculated as TauOne*||MomentumResidual||.
     * Note that the residual of the momentum equation is evaluated at the element center
     * and that the result has units of velocity (L/T).
     * The error estimate both saved as the elemental ERROR_RATIO variable and returned as rOutput.
     * If rVARIABLE == NODAL_AREA, the element's contribution to nodal area is added to its nodes.
     * @param rVariable Use ERROR_RATIO or NODAL_AREA
     * @param rOutput Returns the error estimate for ERROR_RATIO, unused for NODAL_AREA
     * @param rCurrentProcessInfo Process info instance (will be checked for OSS_SWITCH)
     * @see MarkForRefinement for a use of the error ratio
     */
    void Calculate(const Variable<double>& rVariable,
                   double& rOutput,
                   const ProcessInfo& rCurrentProcessInfo) override
    {
        // if (rVariable == NODAL_AREA)
        // {
        //     // Get the element's geometric parameters
        //     double Area;
        //     array_1d<double, TNumNodes> N;
        //     BoundedMatrix<double, TNumNodes, TDim> DN_DX;
        //     GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        //     // Carefully write results to nodal variables, to avoid parallelism problems
        //     for (unsigned int i = 0; i < TNumNodes; ++i)
        //     {
        //         this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
        //         this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];
        //         this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
        //     }
        // }
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

        KRATOS_ERROR_IF(this->Id() < 1)
            << "StabilizedConvectionDiffusionReactionElement found with Id 0 "
               "or negative"
            << std::endl;

        KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0)
            << "On StabilizedConvectionDiffusionReactionElement -> " << this->Id()
            << "; Area cannot be less than or equal to 0" << std::endl;

        KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME);
        KRATOS_CHECK_VARIABLE_KEY(BOSSAK_ALPHA);
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);

        for (IndexType iNode = 0; iNode < this->GetGeometry().size(); ++iNode)
        {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, this->GetGeometry()[iNode]);
        }

        return 0;

        KRATOS_CATCH("");
    }

    double EvaluateInPoint(const Variable<double>& rVariable,
                           const Vector& rShapeFunction,
                           const int Step = 0) const
    {
        return RansCalculationUtilities().EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
    }

    array_1d<double, 3> EvaluateInPoint(const Variable<array_1d<double, 3>>& rVariable,
                                        const Vector& rShapeFunction,
                                        const int Step = 0) const
    {
        return RansCalculationUtilities().EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
    }

    double GetDivergenceOperator(const Variable<array_1d<double, 3>>& rVariable,
                                 const Matrix& rShapeDerivatives,
                                 const int Step = 0) const
    {
        double value = 0.0;
        const GeometryType& r_geometry = this->GetGeometry();

        for (unsigned int i = 0; i < TNumNodes; ++i)
        {
            const array_1d<double, 3>& r_value =
                r_geometry[i].FastGetSolutionStepValue(rVariable, Step);
            for (unsigned int j = 0; j < TDim; ++j)
            {
                value += r_value[j] * rShapeDerivatives(i, j);
            }
        }

        return value;
    }

    virtual void CalculateConvectionDiffusionReactionData(
        TConvectionDiffusionReactionData& rData,
        double& rEffectiveKinematicViscosity,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo,
        const int Step = 0) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReactionElement "
                        "CalculateConvectionDiffusionReactionData method. "
                        "Please implement it in the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    virtual void CalculateConvectionDiffusionReactionData(
        TConvectionDiffusionReactionData& rData,
        double& rEffectiveKinematicViscosity,
        double& rVariableGradientNorm,
        double& rVariableRelaxedAcceleration,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives,
        const ProcessInfo& rCurrentProcessInfo,
        const int Step = 0) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReactionElement "
                        "CalculateConvectionDiffusionReactionData method. "
                        "Please implement it in the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    ShapeFunctionDerivativesArrayType GetGeometryParameterDerivatives() const
    {
        const GeometryType& r_geometry = this->GetGeometry();
        return RansCalculationUtilities().CalculateGeometryParameterDerivatives(
            r_geometry, this->GetIntegrationMethod());
    }

    virtual double CalculateReactionTerm(const TConvectionDiffusionReactionData& rData,
                                         const ProcessInfo& rCurrentProcessInfo,
                                         const int Step = 0) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReactionElement "
                        "CalculateReactionTerm method. Please implement it in "
                        "the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
    }

    virtual double CalculateSourceTerm(const TConvectionDiffusionReactionData& rData,
                                       const ProcessInfo& rCurrentProcessInfo,
                                       const int Step = 0) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReactionElement "
                        "CalculateSourceTerm method. Please implement it in "
                        "the derrived class."
                     << std::endl;
        KRATOS_CATCH("");
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
        buffer << "StabilizedConvectionDiffusionReactionElement #" << Id();
        return buffer.str();
    }
    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StabilizedConvectionDiffusionReactionElement #" << Id();
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

    /// Determine integration point weights and shape funcition derivatives from the element's geometry.
    virtual void CalculateGeometryData(Vector& rGaussWeights,
                                       Matrix& rNContainer,
                                       ShapeFunctionDerivativesArrayType& rDN_DX) const
    {
        const GeometryType& r_geometry = this->GetGeometry();

        RansCalculationUtilities().CalculateGeometryData(
            r_geometry, this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
    }

    void GetConvectionOperator(BoundedVector<double, TNumNodes>& rOutput,
                               const array_1d<double, 3>& rVector,
                               const Matrix& rShapeDerivatives) const
    {
        rOutput.clear();
        for (unsigned int i = 0; i < TNumNodes; ++i)
            for (unsigned int j = 0; j < TDim; j++)
            {
                rOutput[i] += rVector[j] * rShapeDerivatives(i, j);
            }
    }

    void CalculateGradient(BoundedMatrix<double, TDim, TDim>& rOutput,
                           const Variable<array_1d<double, 3>>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const
    {
        const GeometryType& r_geometry = this->GetGeometry();

        RansCalculationUtilities().CalculateGradient<TDim>(
            rOutput, r_geometry, rVariable, rShapeDerivatives, Step);
    }

    void CalculateGradient(array_1d<double, 3>& rOutput,
                           const Variable<double>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const
    {
        const GeometryType& r_geometry = this->GetGeometry();
        RansCalculationUtilities().CalculateGradient(
            rOutput, r_geometry, rVariable, rShapeDerivatives, Step);
    }

    void CalculateSymmetricGradientMatrix(BoundedMatrix<double, TDim, TDim>& rOutput,
                                          const Variable<array_1d<double, 3>>& rVariable,
                                          const BoundedMatrix<double, TDim, TDim>& rGradientMatrix,
                                          const Matrix& rShapeDerivatives,
                                          const int Step = 0) const
    {
        const double variable_divergence =
            this->GetDivergenceOperator(rVariable, rShapeDerivatives, Step);
        identity_matrix<double> identity(TDim);

        rOutput.clear();
        noalias(rOutput) = rGradientMatrix + trans(rGradientMatrix) -
                           (2.0 / 3.0) * variable_divergence * identity;
    }

    void AddLumpedMassMatrix(MatrixType& rMassMatrix, const double Mass)
    {
        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
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

    double GetDeltaTime(const ProcessInfo& rProcessInfo) const
    {
        if (this->Has(DELTA_TIME))
        {
            return this->GetValue(DELTA_TIME);
        }
        else
        {
            return rProcessInfo[DELTA_TIME];
        }
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

}; // Class StabilizedConvectionDiffusionReactionElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::istream& operator>>(
    std::istream& rIStream,
    StabilizedConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const StabilizedConvectionDiffusionReactionElement<TDim, TNumNodes, TConvectionDiffusionReactionData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

///@}

} // namespace Kratos.

#endif // KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ELEMENT_H_INCLUDED defined
