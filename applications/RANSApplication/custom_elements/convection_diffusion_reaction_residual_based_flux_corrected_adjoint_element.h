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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_ELEMENT)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_ELEMENT

// System includes
#include <array>

// External includes

// Project includes
#include "includes/element.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/time_discretization.h"

// Application includes
#include "custom_elements/convection_diffusion_reaction_stabilization_adjoint_utilities.h"
#include "custom_elements/convection_diffusion_reaction_stabilization_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <unsigned int TDim, unsigned int TNumNodes, class TConvectionDiffusionReactionAdjointData>
class ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Element
    KRATOS_CLASS_POINTER_DEFINITION(ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement);

    /// base type: an GeometricalObject that automatically has a unique number
    using BaseType = Element;

    /// definition of node type (default is: Node<3>)
    using NodeType = Node<3>;

    /**
     * Properties are used to store any parameters
     * related to the constitutive law
     */
    using PropertiesType = Properties;

    /// definition of the geometry type with given NodeType
    using GeometryType = Geometry<NodeType>;

    /// definition of nodes container type, redefined from GeometryType
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;

    using VectorType = Vector;

    using MatrixType = Matrix;

    using IndexType = std::size_t;

    using EquationIdVectorType = std::vector<std::size_t>;

    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    /// Type definition for integration methods
    using IntegrationMethod = GeometryData::IntegrationMethod;

    using ConvectionDiffusionReactionAdjointDataType = TConvectionDiffusionReactionAdjointData;

    using CurrentElementType =
        ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement<TDim, TNumNodes, TConvectionDiffusionReactionAdjointData>;

    // overriden methods
    using BaseType::CalculateFirstDerivativesLHS;
    using BaseType::CalculateLeftHandSide;
    using BaseType::CalculateSecondDerivativesLHS;
    using BaseType::CalculateSensitivityMatrix;

    ///@}

    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    explicit ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement(IndexType NewId = 0)
        : BaseType(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement(IndexType NewId,
                                                                        const NodesArrayType& ThisNodes)
        : BaseType(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement(IndexType NewId,
                                                                        GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement(
        IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement(
        ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement const& rOther)
        : BaseType(rOther)
    {
    }

    ~ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.

    ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement& operator=(
        ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement const& rOther)
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
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(
            NewId, Element::GetGeometry().Create(ThisNodes), pProperties);
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
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(NewId, pGeom, pProperties);
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
        KRATOS_TRY
        return Kratos::make_intrusive<CurrentElementType>(
            NewId, Element::GetGeometry().Create(ThisNodes), Element::pGetProperties());
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
        const Variable<double>& r_adjoint_Variable =
            TConvectionDiffusionReactionAdjointData::GetAdjointScalarVariable();
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
        const Variable<double>& r_adjoint_Variable =
            TConvectionDiffusionReactionAdjointData::GetAdjointScalarVariable();
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
        const Variable<double>& r_adjoint_variable =
            TConvectionDiffusionReactionAdjointData::GetAdjointScalarVariable();
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
            TConvectionDiffusionReactionAdjointData::GetAdjointScalarRateVariable();
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
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        BoundedMatrix<double, TNumNodes, TNumNodes> local_matrix;
        this->CalculateFirstDerivativesLHS(local_matrix, rCurrentProcessInfo);
        if (rLeftHandSideMatrix.size1() != local_matrix.size1() ||
            rLeftHandSideMatrix.size2() != local_matrix.size2())
            rLeftHandSideMatrix.resize(local_matrix.size1(), local_matrix.size2(), false);

        noalias(rLeftHandSideMatrix) = local_matrix;
    }

    void CalculateFirstDerivativesLHS(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
        const Variable<double>& r_derivative_variable =
            TConvectionDiffusionReactionAdjointData::GetScalarVariable();
        CalculateElementSteadyResidualScalarDerivatives(
            rLeftHandSideMatrix, r_derivative_variable, rCurrentProcessInfo);
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

        BoundedMatrix<double, TNumNodes, TNumNodes> output_matrix_scalar;
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes> output_matrix_vector;

        bool is_scalar = true;

        if (rVariable == TURBULENT_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
        {
            CalculateElementSteadyResidualScalarDerivatives(
                output_matrix_scalar, TURBULENT_ENERGY_DISSIPATION_RATE, rCurrentProcessInfo);
        }
        else if (rVariable == TURBULENT_KINETIC_ENERGY_PARTIAL_DERIVATIVE)
        {
            CalculateElementSteadyResidualScalarDerivatives(
                output_matrix_scalar, TURBULENT_KINETIC_ENERGY, rCurrentProcessInfo);
        }
        else if (rVariable == TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_PARTIAL_DERIVATIVE)
        {
            CalculateElementSteadyResidualScalarDerivatives(
                output_matrix_scalar,
                TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, rCurrentProcessInfo);
        }
        else if (rVariable == VELOCITY_PARTIAL_DERIVATIVE)
        {
            is_scalar = false;
            CalculateElementSteadyResidualVelocityDerivatives(
                output_matrix_vector, rCurrentProcessInfo);
        }

        if (is_scalar)
        {
            if (rOutput.size1() != TNumNodes || rOutput.size2() != TNumNodes)
            {
                rOutput.resize(TNumNodes, TNumNodes);
            }
            noalias(rOutput) = output_matrix_scalar;
        }
        else
        {
            if (rOutput.size1() != TNumNodes * TDim || rOutput.size2() != TNumNodes)
            {
                rOutput.resize(TNumNodes * TDim, TNumNodes);
            }
            noalias(rOutput) = output_matrix_vector;
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
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes> local_matrix;
        CalculateSensitivityMatrix(rSensitivityVariable, local_matrix, rCurrentProcessInfo);

        if (rOutput.size1() != TNumNodes * TDim || rOutput.size2() != TNumNodes)
        {
            rOutput.resize(TNumNodes * TDim, TNumNodes);
        }

        noalias(rOutput) = local_matrix;
    }

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if (rSensitivityVariable == SHAPE_SENSITIVITY)
        {
            CalculateElementSteadyResidualShapeDerivatives(rOutput, rCurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR << "Unsupported sensitivity variable.";
        }

        KRATOS_CATCH("")
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

        return 0;

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjoint"
                  "Element #"
               << this->Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjoi"
                    "ntElement #"
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

    void AddPrimalDampingMatrixGaussPointContributions(
        BoundedMatrix<double, TNumNodes, TNumNodes>& rOutput,
        const double EffectiveKinematicViscosity,
        const double ReactionTerm,
        const double Tau,
        const double Weight,
        const Vector& rVelocityConvectiveTerms,
        const Vector& rShapeFunctions,
        const Matrix& rShapeFunctionDerivatives) const
    {
        // calculating primal damping matrix
        for (IndexType a = 0; a < TNumNodes; ++a)
        {
            for (IndexType b = 0; b < TNumNodes; ++b)
            {
                double dNa_dNb = 0.0;
                for (IndexType i = 0; i < TDim; ++i)
                    dNa_dNb += rShapeFunctionDerivatives(a, i) *
                               rShapeFunctionDerivatives(b, i);

                double value = 0.0;

                value += rShapeFunctions[a] * rVelocityConvectiveTerms[b];
                value += rShapeFunctions[a] * ReactionTerm * rShapeFunctions[b];
                value += EffectiveKinematicViscosity * dNa_dNb;

                // // Adding SUPG stabilization terms
                value += Tau *
                         (rVelocityConvectiveTerms[a] + ReactionTerm * rShapeFunctions[a]) *
                         rVelocityConvectiveTerms[b];
                value += Tau *
                         (rVelocityConvectiveTerms[a] + ReactionTerm * rShapeFunctions[a]) *
                         ReactionTerm * rShapeFunctions[b];

                rOutput(a, b) += Weight * value;
            }
        }
    }

    double GetDeltaTime(const ProcessInfo& rProcessInfo) const
    {
        return -1.0 * rProcessInfo[DELTA_TIME];
    }

    void CalculateElementSteadyResidualScalarDerivatives(
        BoundedMatrix<double, TNumNodes, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        using adjoint_utilities =
            ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;

        noalias(rOutput) = ZeroMatrix(TNumNodes, TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        const GeometryType& r_geometry = this->GetGeometry();
        TConvectionDiffusionReactionAdjointData r_current_data(r_geometry);

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const double element_length = this->GetGeometry().Length();
        const double discrete_upwind_operator_coefficient =
            rCurrentProcessInfo[RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT];
        const double diagonal_positivity_preserving_coefficient =
            rCurrentProcessInfo[RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT];

        r_current_data.CalculateConstants(rCurrentProcessInfo);

        const Variable<double>& r_primal_variable =
            TConvectionDiffusionReactionAdjointData::GetScalarVariable();

        double scalar_multiplier = 0.0;

        array_1d<double, 3> scalar_variable_gradient;

        BoundedVector<double, TNumNodes> source_term_derivatives,
            reaction_term_derivatives, tau_derivatives, primal_values_vector,
            velocity_convective_terms, effective_kinematic_viscosity_derivatives,
            absolute_residual_unrelated_scalar_derivatives,
            absolute_residual_related_scalar_derivatives,
            scalar_multiplier_unrelated_scalar_derivatives,
            scalar_multiplier_related_scalar_derivatives;

        BoundedMatrix<double, TNumNodes, TDim> effective_velocity_derivatives;

        BoundedMatrix<double, TNumNodes, TNumNodes> primal_damping_matrix,
            discrete_diffusion_matrix, velocity_convective_term_derivatives,
            discrete_diffusion_matrix_residual_contribution_derivatives,
            positivity_preserving_matrix_residual_contribution_derivatives;

        BoundedVector<BoundedMatrix<double, TNumNodes, TNumNodes>, TNumNodes> primal_damping_matrix_derivatives;

        primal_damping_matrix.clear();
        scalar_multiplier_unrelated_scalar_derivatives.clear();
        scalar_multiplier_related_scalar_derivatives.clear();

        for (IndexType c = 0; c < TNumNodes; ++c)
        {
            primal_damping_matrix_derivatives[c].clear();
        }

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& r_shape_functions = row(shape_functions, g);
            const double weight = gauss_weights[g];

            r_current_data.CalculateGaussPointData(r_shape_functions, r_shape_derivatives);

            // primal calculations
            const array_1d<double, 3>& velocity = r_current_data.CalculateEffectiveVelocity(
                r_shape_functions, r_shape_derivatives);
            const double effective_kinematic_viscosity =
                r_current_data.CalculateEffectiveKinematicViscosity(
                    r_shape_functions, r_shape_derivatives);

            const double reaction = r_current_data.CalculateReactionTerm(
                r_shape_functions, r_shape_derivatives);
            const double source = r_current_data.CalculateSourceTerm(
                r_shape_functions, r_shape_derivatives);

            const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
                element_length, norm_2(velocity), reaction, effective_kinematic_viscosity,
                bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            this->CalculateGradient(scalar_variable_gradient, r_primal_variable,
                                    r_shape_derivatives);

            const double scalar_variable_value =
                this->EvaluateInPoint(r_primal_variable, r_shape_functions);

            const double velocity_dot_scalar_variable_gradient =
                this->GetDivergence(velocity, scalar_variable_gradient);

            // adjoint calculations
            r_current_data.CalculateSourceTermDerivatives(
                source_term_derivatives, rDerivativeVariable, r_shape_functions,
                r_shape_derivatives);

            r_current_data.CalculateReactionTermDerivatives(
                reaction_term_derivatives, rDerivativeVariable,
                r_shape_functions, r_shape_derivatives);

            r_current_data.CalculateEffectiveVelocityDerivatives(
                effective_velocity_derivatives, rDerivativeVariable,
                r_shape_functions, r_shape_derivatives);

            r_current_data.CalculateEffectiveKinematicViscosityDerivatives(
                effective_kinematic_viscosity_derivatives, rDerivativeVariable,
                r_shape_functions, r_shape_derivatives);

            adjoint_utilities::CalculateTauScalarDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity,
                element_length, reaction, effective_kinematic_viscosity_derivatives,
                reaction_term_derivatives);

            double residual = this->EvaluateInPoint(
                TConvectionDiffusionReactionAdjointData::GetScalarRelaxedRateVariable(),
                r_shape_functions);
            residual += velocity_dot_scalar_variable_gradient;
            residual += reaction * scalar_variable_value;
            residual -= source;

            adjoint_utilities::CalculateAbsoluteResidualUnrelatedScalarDerivatives(
                absolute_residual_unrelated_scalar_derivatives, residual,
                scalar_variable_value, effective_velocity_derivatives, scalar_variable_gradient,
                reaction_term_derivatives, source_term_derivatives);

            adjoint_utilities::CalculateAbsoluteResidualRelatedScalarDerivativeAdditionalTerms(
                absolute_residual_related_scalar_derivatives, residual, reaction,
                ZeroVector(TNumNodes), velocity, r_shape_functions, r_shape_derivatives);

            residual = std::abs(residual);
            if (scalar_variable_value > 0.0)
                scalar_multiplier += residual * tau / scalar_variable_value;

            adjoint_utilities::AddRFCBetaUnrelatedScalarDerivatives(
                scalar_multiplier_unrelated_scalar_derivatives,
                scalar_variable_value, residual, tau,
                absolute_residual_unrelated_scalar_derivatives, tau_derivatives);

            adjoint_utilities::AddRFCBetaRelatedScalarDerivativeAdditionalTerms(
                scalar_multiplier_related_scalar_derivatives,
                scalar_variable_value, residual, tau,
                absolute_residual_related_scalar_derivatives, r_shape_functions);

            // calculating primal damping matrix
            this->AddPrimalDampingMatrixGaussPointContributions(
                primal_damping_matrix, effective_kinematic_viscosity, reaction, tau, weight,
                velocity_convective_terms, r_shape_functions, r_shape_derivatives);

            noalias(velocity_convective_term_derivatives) =
                prod(effective_velocity_derivatives, trans(r_shape_derivatives));

            // calculating primal damping matrix derivatives
            for (IndexType c = 0; c < TNumNodes; ++c)
            {
                BoundedMatrix<double, TNumNodes, TNumNodes>& r_derivatives_matrix =
                    primal_damping_matrix_derivatives[c];

                const BoundedVector<double, TDim>& r_velocity_derivative =
                    row(effective_velocity_derivatives, c);

                for (IndexType a = 0; a < TNumNodes; ++a)
                {
                    for (IndexType b = 0; b < TNumNodes; ++b)
                    {
                        const double dNa_dNb = this->GetDivergence(
                            row(r_shape_derivatives, a), row(r_shape_derivatives, b));

                        double value = 0.0;

                        value += r_shape_functions[a] *
                                 velocity_convective_term_derivatives(c, b);
                        value += r_shape_functions[a] *
                                 reaction_term_derivatives[c] * r_shape_functions[b];
                        value += dNa_dNb * effective_kinematic_viscosity_derivatives[c];

                        // adding derivatives of SUPG
                        value += (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 tau_derivatives[c] * velocity_convective_terms[b];
                        value += (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 tau * velocity_convective_term_derivatives(c, b);
                        value += velocity_convective_term_derivatives(c, a) *
                                 tau * velocity_convective_terms[b];
                        value += reaction_term_derivatives[c] * r_shape_functions[a] *
                                 tau * velocity_convective_terms[b];

                        value += (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 tau_derivatives[c] * reaction * r_shape_functions[b];
                        value += (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 tau * reaction_term_derivatives[c] *
                                 r_shape_functions[b];
                        value += velocity_convective_term_derivatives(c, a) *
                                 tau * reaction * r_shape_functions[b];
                        value += reaction_term_derivatives[c] * r_shape_functions[a] *
                                 tau * reaction * r_shape_functions[b];

                        r_derivatives_matrix(a, b) += weight * value;
                    }
                }

                const double dU_ic_dPhi_i = this->GetDivergence(
                    scalar_variable_gradient, r_velocity_derivative);

                for (IndexType a = 0; a < TNumNodes; ++a)
                {
                    const double dNa_i_dPhi_i = this->GetDivergence(
                        scalar_variable_gradient, row(r_shape_derivatives, a));
                    double value = 0.0;

                    // add RHS derivative contributions
                    value += r_shape_functions[a] * source_term_derivatives[c];
                    value += tau_derivatives[c] *
                             (velocity_convective_terms[a] +
                              reaction * r_shape_functions[a]) *
                             source;
                    value += tau *
                             (velocity_convective_terms[a] +
                              reaction * r_shape_functions[a]) *
                             source_term_derivatives[c];
                    value +=
                        tau * (reaction_term_derivatives[c] * r_shape_functions[a]) * source;
                    value += tau * velocity_convective_term_derivatives(c, a) * source;

                    // add LHS derivative contributions
                    value -= r_shape_functions[a] * dU_ic_dPhi_i;
                    value -= r_shape_functions[a] * reaction_term_derivatives[c] * scalar_variable_value;
                    value -= effective_kinematic_viscosity_derivatives[c] * dNa_i_dPhi_i;
                    value -= tau_derivatives[c] *
                             (velocity_convective_terms[a] +
                              reaction * r_shape_functions[a]) *
                             velocity_dot_scalar_variable_gradient;
                    value -= tau *
                             (velocity_convective_terms[a] +
                              reaction * r_shape_functions[a]) *
                             dU_ic_dPhi_i;
                    value -= tau * reaction_term_derivatives[c] *
                             r_shape_functions[a] * velocity_dot_scalar_variable_gradient;
                    value -= tau * velocity_convective_term_derivatives(c, a) *
                             velocity_dot_scalar_variable_gradient;

                    value -= tau_derivatives[c] *
                             (velocity_convective_terms[a] +
                              reaction * r_shape_functions[a]) *
                             reaction * scalar_variable_value;
                    value -= tau *
                             (velocity_convective_terms[a] +
                              reaction * r_shape_functions[a]) *
                             reaction_term_derivatives[c] * scalar_variable_value;
                    value -= tau * reaction_term_derivatives[c] *
                             r_shape_functions[a] * reaction * scalar_variable_value;
                    value -= tau * velocity_convective_term_derivatives(c, a) *
                             reaction * scalar_variable_value;

                    rOutput(c, a) += value * weight;
                }
            }
        }

        // adding derivatives of discrete upwind operator
        this->GetPrimalValuesVector(primal_values_vector);
        double matrix_norm;
        ConvectionDiffusionReactionStabilizationUtilities::CalculateDiscreteUpwindOperator<TNumNodes>(
            matrix_norm, discrete_diffusion_matrix, primal_damping_matrix);

        const Vector& discrete_diffusion_matrix_residual_contribution =
            prod(discrete_diffusion_matrix, primal_values_vector);

        adjoint_utilities::template CalculateDiscreteUpwindOperatorResidualContributionDerivatives<TNumNodes>(
            discrete_diffusion_matrix_residual_contribution_derivatives, primal_values_vector,
            primal_damping_matrix, primal_damping_matrix_derivatives);

        noalias(rOutput) -= discrete_diffusion_matrix_residual_contribution_derivatives *
                            discrete_upwind_operator_coefficient * scalar_multiplier;

        // adding derivatives of positivity preserving matrix
        const double positivity_preserving_matrix_coefficient =
            ConvectionDiffusionReactionStabilizationUtilities::CalculatePositivityPreservingMatrix(
                primal_damping_matrix);

        adjoint_utilities::template CalculatePositivityPreservingMatrixResidualContributionDerivatives<TNumNodes>(
            positivity_preserving_matrix_residual_contribution_derivatives,
            positivity_preserving_matrix_coefficient, primal_values_vector,
            primal_damping_matrix, primal_damping_matrix_derivatives);

        noalias(rOutput) -=
            positivity_preserving_matrix_residual_contribution_derivatives *
            diagonal_positivity_preserving_coefficient * scalar_multiplier;

        for (IndexType c = 0; c < TNumNodes; ++c)
        {
            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                double value = 0.0;

                value -= scalar_multiplier_unrelated_scalar_derivatives[c] *
                         discrete_diffusion_matrix_residual_contribution[a] *
                         discrete_upwind_operator_coefficient;
                value -= scalar_multiplier_unrelated_scalar_derivatives[c] *
                         primal_values_vector[a] * positivity_preserving_matrix_coefficient *
                         diagonal_positivity_preserving_coefficient;

                rOutput(c, a) += value;
            }
        }

        // add additional stuff required from primal LHS
        if (r_primal_variable == rDerivativeVariable)
        {
            noalias(rOutput) -= trans(primal_damping_matrix);
            noalias(rOutput) -= discrete_diffusion_matrix *
                                (discrete_upwind_operator_coefficient * scalar_multiplier);
            noalias(rOutput) -=
                IdentityMatrix(TNumNodes) * positivity_preserving_matrix_coefficient *
                diagonal_positivity_preserving_coefficient * scalar_multiplier;

            for (IndexType c = 0; c < TNumNodes; ++c)
            {
                for (IndexType a = 0; a < TNumNodes; ++a)
                {
                    double value = 0.0;

                    value -= scalar_multiplier_related_scalar_derivatives[c] *
                             discrete_diffusion_matrix_residual_contribution[a] *
                             discrete_upwind_operator_coefficient;
                    value -= scalar_multiplier_related_scalar_derivatives[c] *
                             primal_values_vector[a] * positivity_preserving_matrix_coefficient *
                             diagonal_positivity_preserving_coefficient;

                    rOutput(c, a) += value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AddElementDynamicResidualContributionScalarDerivatives(
        BoundedMatrix<double, TNumNodes, TNumNodes>& rOutput,
        const Variable<double>& rDerivativeVariable,
        const ProcessInfo& rCurrentProcessInfo) const
    {
    }

    void CalculateElementSteadyResidualVelocityDerivatives(
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        using adjoint_utilities =
            ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;

        noalias(rOutput) = ZeroMatrix(TNumNodes * TDim, TNumNodes);

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        const GeometryType& r_geometry = this->GetGeometry();
        TConvectionDiffusionReactionAdjointData r_current_data(r_geometry);

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const double element_length = this->GetGeometry().Length();
        const double discrete_upwind_operator_coefficient =
            rCurrentProcessInfo[RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT];
        const double diagonal_positivity_preserving_coefficient =
            rCurrentProcessInfo[RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT];

        r_current_data.CalculateConstants(rCurrentProcessInfo);

        const Variable<double>& r_primal_variable =
            TConvectionDiffusionReactionAdjointData::GetScalarVariable();

        double scalar_multiplier = 0.0;

        array_1d<double, 3> scalar_variable_gradient;

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            scalar_variable_gradient_dot_shape_derivatives;
        BoundedVector<double, TNumNodes * TDim> velocity_dot_scalar_variable_gradient_derivatives;
        BoundedMatrix<double, TNumNodes, TNumNodes> primal_damping_matrix,
            discrete_diffusion_matrix;
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes> velocity_convective_term_derivatives,
            discrete_diffusion_matrix_residual_contribution_derivatives,
            positivity_preserving_matrix_residual_contribution_derivatives;
        BoundedMatrix<double, TNumNodes * TDim, TDim> effective_velocity_derivatives;
        BoundedMatrix<double, TNumNodes, TDim> source_term_derivatives,
            effective_kinematic_viscosity_derivatives, reaction_term_derivatives,
            tau_derivatives, effective_velocity_magnitude_derivatives,
            absolute_residual_derivatives, scalar_multiplier_derivatives;

        BoundedVector<BoundedMatrix<double, TNumNodes, TNumNodes>, TNumNodes * TDim> primal_damping_matrix_derivatives;

        primal_damping_matrix.clear();
        scalar_multiplier_derivatives.clear();
        for (IndexType c = 0; c < TNumNodes * TDim; ++c)
        {
            primal_damping_matrix_derivatives[c].clear();
        }

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& r_shape_functions = row(shape_functions, g);
            const double weight = gauss_weights[g];

            r_current_data.CalculateGaussPointData(r_shape_functions, r_shape_derivatives);

            const array_1d<double, 3>& r_velocity = r_current_data.CalculateEffectiveVelocity(
                r_shape_functions, r_shape_derivatives);
            const double velocity = norm_2(r_velocity);
            const double effective_kinematic_viscosity =
                r_current_data.CalculateEffectiveKinematicViscosity(
                    r_shape_functions, r_shape_derivatives);
            const double reaction = r_current_data.CalculateReactionTerm(
                r_shape_functions, r_shape_derivatives);

            const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
                element_length, velocity, reaction, effective_kinematic_viscosity,
                bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

            this->GetConvectionOperator(velocity_convective_terms, r_velocity,
                                        r_shape_derivatives);

            // calculating primal damping matrix
            this->AddPrimalDampingMatrixGaussPointContributions(
                primal_damping_matrix, effective_kinematic_viscosity, reaction, tau, weight,
                velocity_convective_terms, r_shape_functions, r_shape_derivatives);

            const double source = r_current_data.CalculateSourceTerm(
                r_shape_functions, r_shape_derivatives);
            this->CalculateGradient(scalar_variable_gradient, r_primal_variable,
                                    r_shape_derivatives);
            const double velocity_dot_scalar_variable_gradient =
                inner_prod(r_velocity, scalar_variable_gradient);
            const double scalar_variable_value =
                this->EvaluateInPoint(r_primal_variable, r_shape_functions);

            // calculating derivatives
            r_current_data.CalculateSourceTermDerivatives(
                source_term_derivatives, VELOCITY, r_shape_functions, r_shape_derivatives);

            r_current_data.CalculateReactionTermDerivatives(
                reaction_term_derivatives, VELOCITY, r_shape_functions, r_shape_derivatives);

            r_current_data.CalculateEffectiveKinematicViscosityDerivatives(
                effective_kinematic_viscosity_derivatives, VELOCITY,
                r_shape_functions, r_shape_derivatives);

            r_current_data.CalculateEffectiveVelocityDerivatives(
                effective_velocity_derivatives, VELOCITY, r_shape_functions, r_shape_derivatives);

            adjoint_utilities::CalculateEffectiveVelocityMagnitudeVelocityDerivative(
                effective_velocity_magnitude_derivatives, velocity, r_velocity,
                effective_velocity_derivatives);

            adjoint_utilities::CalculateTauVelocityDerivatives(
                tau_derivatives, tau, velocity, effective_kinematic_viscosity,
                element_length, reaction, effective_velocity_magnitude_derivatives,
                effective_kinematic_viscosity_derivatives, reaction_term_derivatives);

            noalias(velocity_convective_term_derivatives) =
                prod(effective_velocity_derivatives, trans(r_shape_derivatives));

            noalias(velocity_dot_scalar_variable_gradient_derivatives) =
                adjoint_utilities::template MatrixVectorProduct<TNumNodes * TDim>(
                    effective_velocity_derivatives, scalar_variable_gradient);

            noalias(scalar_variable_gradient_dot_shape_derivatives) =
                adjoint_utilities::template MatrixVectorProduct<TNumNodes>(
                    r_shape_derivatives, scalar_variable_gradient);

            double residual = this->EvaluateInPoint(
                TConvectionDiffusionReactionAdjointData::GetScalarRelaxedRateVariable(),
                r_shape_functions);
            residual += velocity_dot_scalar_variable_gradient;
            residual += reaction * scalar_variable_value;
            residual -= source;

            adjoint_utilities::CalculateAbsoluteResidualVelocityDerivativeTerms(
                absolute_residual_derivatives, residual, scalar_variable_value,
                scalar_variable_gradient, ZeroMatrix(TNumNodes, TDim),
                effective_velocity_derivatives, reaction_term_derivatives,
                source_term_derivatives);

            residual = std::abs(residual);
            if (scalar_variable_value > 0.0)
                scalar_multiplier += residual * tau / scalar_variable_value;

            adjoint_utilities::AddRFCBetaVelocityDerivative(
                scalar_multiplier_derivatives, scalar_variable_value, residual,
                tau, absolute_residual_derivatives, tau_derivatives);

            BoundedMatrix<double, TNumNodes, TNumNodes> dNa_dNb =
                prod(r_shape_derivatives, trans(r_shape_derivatives));

            for (IndexType c = 0; c < TNumNodes; ++c)
            {
                for (IndexType k = 0; k < TDim; ++k)
                {
                    const IndexType row = c * TDim + k;

                    // constructing primal damping matrix derivatives
                    BoundedMatrix<double, TNumNodes, TNumNodes>& r_derivatives_matrix =
                        primal_damping_matrix_derivatives[row];

                    for (IndexType a = 0; a < TNumNodes; ++a)
                    {
                        for (IndexType b = 0; b < TNumNodes; ++b)
                        {
                            double value = 0.0;

                            value += r_shape_functions[a] *
                                     velocity_convective_term_derivatives(row, b);
                            value += r_shape_functions[a] *
                                     reaction_term_derivatives(c, k) *
                                     r_shape_functions[b];
                            value += dNa_dNb(a, b) *
                                     effective_kinematic_viscosity_derivatives(c, k);

                            // adding derivatives of SUPG
                            value += (velocity_convective_terms[a] +
                                      reaction * r_shape_functions[a]) *
                                     tau_derivatives(c, k) *
                                     velocity_convective_terms[b];
                            value += (velocity_convective_terms[a] +
                                      reaction * r_shape_functions[a]) *
                                     tau * velocity_convective_term_derivatives(row, b);

                            value += velocity_convective_term_derivatives(row, a) *
                                     tau * velocity_convective_terms[b];
                            value += reaction_term_derivatives(c, k) *
                                     r_shape_functions[a] * tau *
                                     velocity_convective_terms[b];

                            value += (velocity_convective_terms[a] +
                                      reaction * r_shape_functions[a]) *
                                     tau_derivatives(c, k) * reaction *
                                     r_shape_functions[b];
                            value += (velocity_convective_terms[a] +
                                      reaction * r_shape_functions[a]) *
                                     tau * reaction_term_derivatives(c, k) *
                                     r_shape_functions[b];
                            value += velocity_convective_term_derivatives(row, a) *
                                     tau * reaction * r_shape_functions[b];
                            value += reaction_term_derivatives(c, k) *
                                     r_shape_functions[a] * tau * reaction *
                                     r_shape_functions[b];

                            r_derivatives_matrix(a, b) += weight * value;
                        }
                    }

                    // constructing residual derivatives
                    for (IndexType a = 0; a < TNumNodes; ++a)
                    {
                        double value = 0.0;

                        // adding contributions from RHS
                        value += r_shape_functions[a] * source_term_derivatives(c, k);
                        value += (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 tau_derivatives(c, k) * source;
                        value += (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 tau * source_term_derivatives(c, k);
                        value += (velocity_convective_term_derivatives(row, a) +
                                  reaction_term_derivatives(c, k) * r_shape_functions[a]) *
                                 tau * source;

                        // adding contributions from LHS
                        value -= r_shape_functions[a] *
                                 velocity_dot_scalar_variable_gradient_derivatives[row];
                        value -= r_shape_functions[a] *
                                 reaction_term_derivatives(c, k) * scalar_variable_value;
                        value -= effective_kinematic_viscosity_derivatives(c, k) *
                                 scalar_variable_gradient_dot_shape_derivatives[a];

                        value -= tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 velocity_dot_scalar_variable_gradient;
                        value -= tau *
                                 (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 velocity_dot_scalar_variable_gradient_derivatives[row];
                        value -= tau *
                                 (velocity_convective_term_derivatives(row, a) +
                                  reaction_term_derivatives(c, k) * r_shape_functions[a]) *
                                 velocity_dot_scalar_variable_gradient;

                        value -= tau_derivatives(c, k) *
                                 (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 reaction * scalar_variable_value;
                        value -= tau *
                                 (velocity_convective_terms[a] +
                                  reaction * r_shape_functions[a]) *
                                 reaction_term_derivatives(c, k) * scalar_variable_value;
                        value -= tau *
                                 (velocity_convective_term_derivatives(row, a) +
                                  reaction_term_derivatives(c, k) * r_shape_functions[a]) *
                                 reaction * scalar_variable_value;

                        rOutput(row, a) += value * weight;
                    }
                }
            }
        }

        // adding derivatives of discrete upwind operator
        BoundedVector<double, TNumNodes> primal_values_vector;
        this->GetPrimalValuesVector(primal_values_vector);
        double matrix_norm;
        ConvectionDiffusionReactionStabilizationUtilities::CalculateDiscreteUpwindOperator<TNumNodes>(
            matrix_norm, discrete_diffusion_matrix, primal_damping_matrix);

        const Vector& discrete_diffusion_matrix_residual_contribution =
            prod(discrete_diffusion_matrix, primal_values_vector);

        adjoint_utilities::template CalculateDiscreteUpwindOperatorResidualContributionDerivatives<TNumNodes * TDim>(
            discrete_diffusion_matrix_residual_contribution_derivatives, primal_values_vector,
            primal_damping_matrix, primal_damping_matrix_derivatives);

        noalias(rOutput) -= discrete_diffusion_matrix_residual_contribution_derivatives *
                            discrete_upwind_operator_coefficient * scalar_multiplier;

        // adding derivatives of positivity preserving matrix
        const double positivity_preserving_matrix_coefficient =
            ConvectionDiffusionReactionStabilizationUtilities::CalculatePositivityPreservingMatrix(
                primal_damping_matrix);

        adjoint_utilities::template CalculatePositivityPreservingMatrixResidualContributionDerivatives<TNumNodes * TDim>(
            positivity_preserving_matrix_residual_contribution_derivatives,
            positivity_preserving_matrix_coefficient, primal_values_vector,
            primal_damping_matrix, primal_damping_matrix_derivatives);

        noalias(rOutput) -=
            positivity_preserving_matrix_residual_contribution_derivatives *
            diagonal_positivity_preserving_coefficient * scalar_multiplier;

        for (IndexType c = 0; c < TNumNodes; ++c)
        {
            for (IndexType k = 0; k < TDim; ++k)
            {
                for (IndexType a = 0; a < TNumNodes; ++a)
                {
                    double value = 0.0;

                    value -= scalar_multiplier_derivatives(c, k) *
                             discrete_diffusion_matrix_residual_contribution[a] *
                             discrete_upwind_operator_coefficient;
                    value -= scalar_multiplier_derivatives(c, k) *
                             primal_values_vector[a] * positivity_preserving_matrix_coefficient *
                             diagonal_positivity_preserving_coefficient;

                    rOutput(c * TDim + k, a) += value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AddElementDynamicResidualContributionVelocityDerivatives(
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) const
    {
    }

    void CalculateElementSteadyResidualShapeDerivatives(
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        using adjoint_utilities =
            ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;

        rOutput.clear();

        const GeometryType& r_geometry = this->GetGeometry();

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const IndexType num_gauss_points = gauss_weights.size();

        ShapeParameter derivative;
        Geometry<Point>::JacobiansType J;
        r_geometry.Jacobian(J, this->GetIntegrationMethod());
        GeometricalSensitivityUtility::ShapeFunctionsGradientType DN_DX_derivatives;
        Geometry<Point>::ShapeFunctionsGradientsType DN_De;
        DN_De = r_geometry.ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const double element_length = this->GetGeometry().Length();

        TConvectionDiffusionReactionAdjointData r_current_data(r_geometry);
        r_current_data.CalculateConstants(rCurrentProcessInfo);

        BoundedVector<double, TNumNodes> velocity_convective_terms,
            velocity_convective_term_derivatives,
            shape_derivative_dot_scalar_variable_gradient,
            shape_derivative_dot_scalar_variable_gradient_derivative,
            shape_derivative_derivative_dot_scalar_variable_gradient;

        array_1d<double, 3> scalar_variable_gradient, scalar_variable_gradient_derivative;

        const Variable<double>& r_primal_variable =
            TConvectionDiffusionReactionAdjointData::GetScalarVariable();

        for (IndexType g = 0; g < num_gauss_points; ++g)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& r_shape_functions = row(shape_functions, g);
            const double weight = gauss_weights[g];

            r_current_data.CalculateGaussPointData(r_shape_functions, r_shape_derivatives);

            const array_1d<double, 3>& effective_velocity =
                r_current_data.CalculateEffectiveVelocity(r_shape_functions, r_shape_derivatives);
            const double effective_velocity_magnitude = norm_2(effective_velocity);
            const double effective_kinematic_viscosity =
                r_current_data.CalculateEffectiveKinematicViscosity(
                    r_shape_functions, r_shape_derivatives);

            const double reaction_term = r_current_data.CalculateReactionTerm(
                r_shape_functions, r_shape_derivatives);
            const double source_term = r_current_data.CalculateSourceTerm(
                r_shape_functions, r_shape_derivatives);

            const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
                element_length, effective_velocity_magnitude, reaction_term,
                effective_kinematic_viscosity, bossak_alpha, bossak_gamma,
                delta_time, dynamic_tau);

            this->CalculateGradient(scalar_variable_gradient, r_primal_variable,
                                    r_shape_derivatives);

            this->GetConvectionOperator(velocity_convective_terms,
                                        effective_velocity, r_shape_derivatives);

            const double effective_velocity_dot_scalar_variable_gradient = inner_prod(effective_velocity, scalar_variable_gradient);
            const double scalar_variable_value = this->EvaluateInPoint(r_primal_variable, r_shape_functions);

            noalias(shape_derivative_dot_scalar_variable_gradient) =
                adjoint_utilities::template MatrixVectorProduct<TNumNodes>(
                    r_shape_derivatives, scalar_variable_gradient);

            const Matrix& rJ = J[g];
            const Matrix& rDN_De = DN_De[g];
            const double inv_detJ = 1.0 / MathUtils<double>::DetMat(rJ);
            GeometricalSensitivityUtility geom_sensitivity(rJ, rDN_De);

            for (IndexType c = 0; c < TNumNodes; ++c)
            {
                const IndexType block_size = c * TDim;
                for (IndexType k = 0; k < TDim; ++k)
                {
                    derivative.NodeIndex = c;
                    derivative.Direction = k;

                    double detJ_derivative;
                    geom_sensitivity.CalculateSensitivity(
                        derivative, detJ_derivative, DN_DX_derivatives);
                    const double weight_derivative =
                        detJ_derivative * inv_detJ * gauss_weights[g];

                    const double source_term_derivative =
                        r_current_data.CalculateSourceTermShapeDerivatives(
                            derivative, r_shape_functions, r_shape_derivatives,
                            detJ_derivative, DN_DX_derivatives);

                    const double reaction_term_derivative =
                        r_current_data.CalculateReactionTermShapeDerivatives(
                            derivative, r_shape_functions, r_shape_derivatives,
                            detJ_derivative, DN_DX_derivatives);

                    const double effective_kinematic_viscosity_derivative =
                        r_current_data.CalculateEffectiveKinematicViscosityShapeDerivatives(
                            derivative, r_shape_functions, r_shape_derivatives,
                            detJ_derivative, DN_DX_derivatives);

                    const array_1d<double, 3>& effective_velocity_derivative =
                        r_current_data.CalculateEffectiveVelocityShapeDerivatives(
                            derivative, r_shape_functions, r_shape_derivatives,
                            detJ_derivative, DN_DX_derivatives);

                    const double effective_velocity_magnitude_derivative =
                        adjoint_utilities::CalculateEffectiveVelocityMagnitudeShapeDerivative(
                            effective_velocity_magnitude, effective_velocity,
                            effective_velocity_derivative);

                    const double tau_derivative = adjoint_utilities::CalculateTauShapeDerivatives(
                        tau, effective_velocity_magnitude,
                        effective_kinematic_viscosity, element_length,
                        reaction_term, effective_velocity_magnitude_derivative,
                        effective_kinematic_viscosity_derivative,
                        reaction_term_derivative, detJ_derivative);

                    this->CalculateGradient(scalar_variable_gradient_derivative,
                                            r_primal_variable, DN_DX_derivatives);

                    this->GetConvectionOperator(velocity_convective_term_derivatives,
                                                effective_velocity, DN_DX_derivatives);

                    const double effective_velocity_derivative_dot_scalar_variable_gradient = inner_prod(effective_velocity_derivative, scalar_variable_gradient);
                    const double effective_velocity_dot_scalar_variable_gradient_derivative = inner_prod(effective_velocity, scalar_variable_gradient_derivative);

                    noalias(shape_derivative_derivative_dot_scalar_variable_gradient) = adjoint_utilities::template MatrixVectorProduct<TNumNodes>(DN_DX_derivatives, scalar_variable_gradient);
                    noalias(shape_derivative_dot_scalar_variable_gradient_derivative) = adjoint_utilities::template MatrixVectorProduct<TNumNodes>(r_shape_derivatives, scalar_variable_gradient_derivative);

                    for (IndexType a = 0; a < TNumNodes; ++a)
                    {
                        double tau_operator = velocity_convective_terms[a] +
                                              reaction_term * r_shape_functions[a];
                        double tau_operator_derivative =
                            velocity_convective_term_derivatives[a] +
                            reaction_term_derivative * r_shape_functions[a];

                        double value = 0.0;

                        // adding RHS contributions
                        value += r_shape_functions[a] * source_term_derivative * weight;
                        value += r_shape_functions[a] * source_term * weight_derivative;

                        value += tau_derivative * tau_operator * source_term * weight;
                        value += tau * tau_operator_derivative * source_term * weight;
                        value += tau * tau_operator * source_term_derivative * weight;
                        value += tau * tau_operator * source_term * weight_derivative;

                        // adding LHS contributions
                        value -= r_shape_functions[a] * effective_velocity_dot_scalar_variable_gradient_derivative * weight;
                        value -= r_shape_functions[a] * effective_velocity_derivative_dot_scalar_variable_gradient * weight;
                        value -= r_shape_functions[a] * effective_velocity_dot_scalar_variable_gradient * weight_derivative;

                        value -= r_shape_functions[a] * reaction_term_derivative * scalar_variable_value * weight;
                        value -= r_shape_functions[a] * reaction_term * scalar_variable_value * weight_derivative;

                        value -= effective_kinematic_viscosity_derivative * shape_derivative_dot_scalar_variable_gradient[a] * weight;
                        value -= effective_kinematic_viscosity * shape_derivative_derivative_dot_scalar_variable_gradient[a] * weight;
                        value -= effective_kinematic_viscosity * shape_derivative_dot_scalar_variable_gradient_derivative[a] * weight;
                        value -= effective_kinematic_viscosity * shape_derivative_dot_scalar_variable_gradient[a] * weight_derivative;

                        value -= tau_derivative * tau_operator * effective_velocity_dot_scalar_variable_gradient * weight;
                        value -= tau * tau_operator_derivative * effective_velocity_dot_scalar_variable_gradient * weight;
                        value -= tau * tau_operator * effective_velocity_derivative_dot_scalar_variable_gradient * weight;
                        value -= tau * tau_operator * effective_velocity_dot_scalar_variable_gradient_derivative * weight;
                        value -= tau * tau_operator * effective_velocity_dot_scalar_variable_gradient * weight_derivative;

                        value -= tau_derivative * tau_operator * reaction_term * scalar_variable_value * weight;
                        value -= tau * tau_operator_derivative * reaction_term * scalar_variable_value * weight;
                        value -= tau * tau_operator * reaction_term_derivative * scalar_variable_value * weight;
                        value -= tau * tau_operator * reaction_term * scalar_variable_value * weight_derivative;

                        rOutput(block_size + k, a) += value;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AddElementDynamicResidualContributionShapeDerivatives(
        BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) const
    {
    }

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

    void CalculateGradient(array_1d<double, 3>& rOutput,
                           const Variable<double>& rVariable,
                           const Matrix& rShapeDerivatives,
                           const int Step = 0) const
    {
        RansCalculationUtilities::CalculateGradient(
            rOutput, this->GetGeometry(), rVariable, rShapeDerivatives, Step);
    }

    double GetDivergence(const Vector& rVector1, const Vector& rVector2) const
    {
        double result = 0.0;
        for (IndexType i = 0; i < TDim; ++i)
        {
            result += rVector1[i] * rVector2[i];
        }
        return result;
    }

    void GetPrimalValuesVector(BoundedVector<double, TNumNodes>& rOutput,
                               const int Step = 0) const
    {
        const GeometryType& r_geometry = this->GetGeometry();
        const Variable<double>& r_primal_variable =
            TConvectionDiffusionReactionAdjointData::GetScalarVariable();
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            rOutput[i] = r_geometry[i].FastGetSolutionStepValue(r_primal_variable, Step);
        }
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

    void AddPrimalSteadyTermScalarDerivatives(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                              const Variable<double>& rDerivativeVariable,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }

    void AddPrimalSteadyTermScalarRateDerivatives(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                                  const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }

    void AddPrimalDampingMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rPrimalDampingMatrix,
                                ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }

    void AddPrimalMassMatrix(BoundedMatrix<double, TNumNodes, TNumNodes>& rPrimalMassMatrix,
                             const double ScalingFactor,
                             const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }

    void AddMassTermScalarDerivatives(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                      const Variable<double>& rDerivativeVariable,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }

    void CalculateResidualShapeSensitivity(BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
                                           const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

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
    ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement<TDim, TNumNodes, TElementData>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TElementData>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointElement<TDim, TNumNodes, TElementData>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

} // namespace Kratos.
#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_ELEMENT defined
