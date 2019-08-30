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

// External includes

// Project includes
#include "includes/element.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_utilities/rans_variable_utils.h"
#include "includes/cfd_variables.h"
#include "rans_modelling_application_variables.h"
#include "stabilized_convection_diffusion_reaction_utilities.h"
#include "utilities/geometrical_sensitivity_utility.h"
#include "utilities/time_discretization.h"

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
 * @tparam TMonolithicAssemblyNodalDofSize Block size of the parent monolithic element
 * @tparam TMonolithicNodalEquationIndex   Equation index in the monolithic element
 */
template <unsigned int TDim, unsigned int TNumNodes, class TElementData, unsigned int TMonolithicAssemblyNodalDofSize = 1, unsigned int TMonolithicNodalEquationIndex = 0>
class StabilizedConvectionDiffusionReactionAdjointElement : public Element
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Element
    KRATOS_CLASS_POINTER_DEFINITION(StabilizedConvectionDiffusionReactionAdjointElement);

    constexpr static unsigned int TMonolithicAssemblyLocalSize =
        TNumNodes * TMonolithicAssemblyNodalDofSize;

    constexpr static unsigned int TVelPrBlockSize = TDim + 1;

    constexpr static unsigned int TVelPrLocalSize = TNumNodes * TVelPrBlockSize;

    constexpr static bool TMonolithicMatrixConstruction =
        (TMonolithicAssemblyLocalSize != TNumNodes);

    /// base type: an GeometricalObject that automatically has a unique number
    typedef Element BaseType;

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
    ///@}

    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId = 0)
        : Element(NewId)
    {
    }

    /**
     * Constructor using an array of nodes
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        const NodesArrayType& ThisNodes)
        : Element(NewId, ThisNodes)
    {
    }

    /**
     * Constructor using Geometry
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry)
    {
    }

    /**
     * Constructor using Properties
     */
    StabilizedConvectionDiffusionReactionAdjointElement(IndexType NewId,
                                                        GeometryType::Pointer pGeometry,
                                                        PropertiesType::Pointer pProperties)
        : Element(NewId, pGeometry, pProperties)
    {
    }

    /**
     * Copy Constructor
     */
    StabilizedConvectionDiffusionReactionAdjointElement(
        StabilizedConvectionDiffusionReactionAdjointElement const& rOther)
        : Element(rOther)
    {
    }

    ~StabilizedConvectionDiffusionReactionAdjointElement() override
    {
    }

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
        KRATOS_TRY
        return Kratos::make_intrusive<StabilizedConvectionDiffusionReactionAdjointElement>(
            NewId, GetGeometry().Create(ThisNodes), pGetProperties());
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
        if (rResult.size() != TMonolithicAssemblyLocalSize)
            rResult.resize(TMonolithicAssemblyLocalSize, false);

        const Variable<double>& r_adjoint_Variable = this->GetAdjointVariable();
        const IndexType dof_index =
            static_cast<IndexType>(rCurrentProcessInfo[this->GetPrimalVariable()]);

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
        {
            rResult[i_node * TMonolithicAssemblyNodalDofSize + dof_index] =
                this->GetGeometry()[i_node].GetDof(r_adjoint_Variable).EquationId();
        }
    }

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override
    {
        if (rElementalDofList.size() != TMonolithicAssemblyLocalSize)
            rElementalDofList.resize(TMonolithicAssemblyLocalSize);

        const Variable<double>& r_adjoint_Variable = this->GetAdjointVariable();
        const IndexType dof_index =
            static_cast<IndexType>(rCurrentProcessInfo[this->GetPrimalVariable()]);

        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
        {
            rElementalDofList[i_node * TMonolithicAssemblyNodalDofSize + dof_index] =
                this->GetGeometry()[i_node].pGetDof(r_adjoint_Variable);
        }
    }

    void GetValuesVector(VectorType& rValues, int Step = 0) override
    {
        if (rValues.size() != TMonolithicAssemblyLocalSize)
            rValues.resize(TMonolithicAssemblyLocalSize, false);

        const Variable<double>& r_adjoint_variable = this->GetAdjointVariable();

        const GeometryType& r_geometry = this->GetGeometry();
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const double r_value =
                r_geometry[i].FastGetSolutionStepValue(r_adjoint_variable, Step);
            rValues[i * TMonolithicAssemblyNodalDofSize + TMonolithicNodalEquationIndex] = r_value;
        }
    }

    void GetFirstDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        if (rValues.size() != TMonolithicAssemblyLocalSize)
            rValues.resize(TMonolithicAssemblyLocalSize, false);

        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            rValues[i * TMonolithicAssemblyNodalDofSize + TMonolithicNodalEquationIndex] = 0.0;
        }
    }

    void GetSecondDerivativesVector(VectorType& rValues, int Step = 0) override
    {
        if (rValues.size() != TMonolithicAssemblyLocalSize)
            rValues.resize(TMonolithicAssemblyLocalSize, false);

        const Variable<double>& r_adjoint_second_variable =
            this->GetAdjointSecondVariable();

        const GeometryType& r_geometry = this->GetGeometry();
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            const double r_value = r_geometry[i].FastGetSolutionStepValue(
                r_adjoint_second_variable, Step);
            rValues[i * TMonolithicAssemblyNodalDofSize + TMonolithicNodalEquationIndex] = r_value;
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
        if (!TMonolithicMatrixConstruction)
        {
            if (rLeftHandSideMatrix.size1() != TNumNodes ||
                rLeftHandSideMatrix.size2() != TNumNodes)
                rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);

            rLeftHandSideMatrix.clear();
        }
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
        const Variable<double>& r_derivative_variable = this->GetPrimalVariable();
        CalculateElementTotalResidualScalarDerivatives(
            rLeftHandSideMatrix, r_derivative_variable, rCurrentProcessInfo);
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
        if (!TMonolithicMatrixConstruction)
        {
            if (rLeftHandSideMatrix.size1() != TMonolithicAssemblyLocalSize ||
                rLeftHandSideMatrix.size2() != TMonolithicAssemblyLocalSize)
                rLeftHandSideMatrix.resize(TMonolithicAssemblyLocalSize,
                                           TMonolithicAssemblyLocalSize, false);
            rLeftHandSideMatrix.clear();
        }

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
                   Matrix& Output,
                   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (rVariable == RANS_VELOCITY_PRESSURE_PARTIAL_DERIVATIVE)
        {
            CalculateElementTotalResidualVelocityDerivatives(Output, rCurrentProcessInfo);
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
            this->CalculateElementTotalResidualShapeSensitivity(rOutput, rCurrentProcessInfo);
        }
        else
        {
            KRATOS_ERROR << "Sensitivity variable " << rSensitivityVariable
                         << " not supported." << std::endl;
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

        const unsigned int primal_dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[r_primal_variable]);
        if (TMonolithicNodalEquationIndex != primal_dof_index)
        {
            KRATOS_ERROR << this->Info() << "'s equation index and Processinfo equation index mismatch. Processinfo.EquationIndex != Element.EquationIndex [ "
                         << primal_dof_index
                         << " != " << TMonolithicNodalEquationIndex << " ]";
        }

        KRATOS_CHECK_VARIABLE_KEY(r_primal_variable);
        KRATOS_CHECK_VARIABLE_KEY(r_primal_relaxed_rate_variable);
        KRATOS_CHECK_VARIABLE_KEY(VELOCITY);
        KRATOS_CHECK_VARIABLE_KEY(BOSSAK_ALPHA);
        KRATOS_CHECK_VARIABLE_KEY(NEWMARK_GAMMA);
        KRATOS_CHECK_VARIABLE_KEY(DELTA_TIME);
        KRATOS_CHECK_VARIABLE_KEY(r_adjoint_variable);

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

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "StabilizedConvectionDiffusionReactionAdjointElement #" << Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StabilizedConvectionDiffusionReactionAdjointElement #" << Id();
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        pGetGeometry()->PrintData(rOStream);
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

        return RANS_SCALAR_1_ADJOINT_1;

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

        return RANS_SCALAR_1_ADJOINT_1;

        KRATOS_CATCH("");
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
    virtual void CalculateElementData(TElementData& rData,
                                      const Vector& rShapeFunctions,
                                      const Matrix& rShapeFunctionDerivatives,
                                      const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY;
        KRATOS_ERROR << "Attempting to call base "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
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
     * @param rCurrentData
     * @param rCurrentProcessInfo
     * @return double
     */
    virtual double CalculateEffectiveKinematicViscosity(const TElementData& rCurrentData,
                                                        const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateEffectiveKinematicViscosity "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates reaction coefficient
     *
     * This method calculates reaction coefficient (i.e. $s_\phi$).
     * This method is called for each gauss point.
     * This method should be implemented by the derrived class.
     *
     * @param rCurrentData
     * @param rCurrentProcessInfo
     * @return double
     */
    virtual double CalculateReactionTerm(const TElementData& rCurrentData,
                                         const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateReactionTerm "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

        KRATOS_CATCH("");
    }

    /**
     * @brief Calculates source term
     *
     * This method calculates source term (i.e. $f_\phi$).
     * This method is called for each gauss point.
     * This method should be implemented by the derrived class.
     *
     * @param rCurrentData
     * @param rCurrentProcessInfo
     * @return double
     */
    virtual double CalculateSourceTerm(const TElementData& rCurrentData,
                                       const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_TRY

        KRATOS_ERROR << "Calling base CalculateSourceTerm "
                        "method in "
                        "StabilizedConvectionDiffusionReactionAdjointElement "
                        "class. Please implement it in the derrived class.";

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
        Vector& rOutput,
        const Variable<double>& rDerivativeVariable,
        const TElementData& rCurrentData,
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
    virtual void CalculateReactionTermScalarDerivatives(Vector& rOutput,
                                                        const Variable<double>& rDerivativeVariable,
                                                        const TElementData& rCurrentData,
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
     * @param rOutput              Output vector containing partial derivatives for each node
     * @param rDerivativeVariable  Derivative variable (i.e. $w$)
     * @param rCurrentData         Data required to calculate partial derivatives
     * @param rCurrentProcessInfo  Current process info
     */
    virtual void CalculateSourceTermScalarDerivatives(Vector& rOutput,
                                                      const Variable<double>& rDerivativeVariable,
                                                      const TElementData& rCurrentData,
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
     * @param rOutput             Output vector containing partial derivatives for each node in rows, and each dimension in columns
     * @param rCurrentData        Data required to calculate partial derivatives
     * @param rCurrentProcessInfo Current process info
     */
    virtual void CalculateEffectiveKinematicViscosityVelocityDerivatives(
        Matrix& rOutput, const TElementData& rCurrentData, const ProcessInfo& rCurrentProcessInfo) const
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
     * @param rOutput             Output vector containing partial derivatives for each node in rows, and each dimension in columns
     * @param rCurrentData        Data required to calculate partial derivatives
     * @param rCurrentProcessInfo Current process info
     */
    virtual void CalculateReactionTermVelocityDerivatives(Matrix& rOutput,
                                                          const TElementData& rCurrentData,
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
     * @param rOutput             Output vector containing partial derivatives for each node in rows, and each dimension in columns
     * @param rCurrentData        Data required to calculate partial derivatives
     * @param rCurrentProcessInfo Current process info
     */
    virtual void CalculateSourceTermVelocityDerivatives(Matrix& rOutput,
                                                        const TElementData& rCurrentData,
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
     * @param rCurrentData         Data required to calculate partial derivatives
     * @param rShapeDerivative     Current derivative (i.e. $x^c_k$, where $c$ is the node, $k$ is the dimension)
     * @param detJ_deriv           Derivative of determinant of jacobian (i.e. $|J|$) w.r.t. $x^c_k$
     * @param rDN_Dx_deriv         Derivative of shape function gradients w.r.t. $x^c_k$
     * @param rCurrentProcessInfo  Current process info
     * @return double              Effective kinematic viscosity shape derivative w.r.t. $x^c_k$
     */
    virtual double CalculateEffectiveKinematicViscosityShapeSensitivity(
        const TElementData& rCurrentData,
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
     * @param rCurrentData         Data required to calculate partial derivatives
     * @param rShapeDerivative     Current derivative (i.e. $x^c_k$, where $c$ is the node, $k$ is the dimension)
     * @param detJ_deriv           Derivative of determinant of jacobian (i.e. $|J|$) w.r.t. $x^c_k$
     * @param rDN_Dx_deriv         Derivative of shape function gradients w.r.t. $x^c_k$
     * @param rCurrentProcessInfo  Current process info
     * @return double              Reaction coefficient shape derivative w.r.t. $x^c_k$
     */
    virtual double CalculateReactionTermShapeSensitivity(
        const TElementData& rCurrentData,
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
     * @param rCurrentData         Data required to calculate partial derivatives
     * @param rShapeDerivative     Current derivative (i.e. $x^c_k$, where $c$ is the node, $k$ is the dimension)
     * @param detJ_deriv           Derivative of determinant of jacobian (i.e. $|J|$) w.r.t. $x^c_k$
     * @param rDN_Dx_deriv         Derivative of shape function gradients w.r.t. $x^c_k$
     * @param rCurrentProcessInfo  Current process info
     * @return double              Source term shape derivative w.r.t. $x^c_k$
     */
    virtual double CalculateSourceTermShapeSensitivity(
        const TElementData& rCurrentData,
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
    void CalculateElementTotalResidualScalarDerivatives(Matrix& rResidualDerivatives,
                                                        const Variable<double>& rDerivativeVariable,
                                                        const ProcessInfo& rCurrentProcessInfo)
    {
        if (!TMonolithicMatrixConstruction)
        {
            if (rResidualDerivatives.size1() != TMonolithicAssemblyLocalSize ||
                rResidualDerivatives.size2() != TMonolithicAssemblyLocalSize)
                rResidualDerivatives.resize(TMonolithicAssemblyLocalSize,
                                            TMonolithicAssemblyLocalSize, false);
            rResidualDerivatives.clear();
        }

        AddPrimalSteadyTermScalarDerivatives(
            rResidualDerivatives, rDerivativeVariable, rCurrentProcessInfo);
        AddMassTermScalarDerivatives(rResidualDerivatives, rDerivativeVariable,
                                     rCurrentProcessInfo);
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

        RansCalculationUtilities().CalculateGeometryData(
            r_geometry, this->GetIntegrationMethod(), rGaussWeights, rNContainer, rDN_DX);
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
        return RansCalculationUtilities().CalculateGeometryParameterDerivatives(
            r_geometry, this->GetIntegrationMethod());
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
        return RansCalculationUtilities().EvaluateInPoint(
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
        return RansCalculationUtilities().EvaluateInPoint(
            this->GetGeometry(), rVariable, rShapeFunction, Step);
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
        for (unsigned int i = 0; i < TNumNodes; ++i)
            for (unsigned int j = 0; j < TDim; j++)
            {
                rOutput[i] += rVector[j] * rShapeDerivatives(i, j);
            }
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

        RansCalculationUtilities().CalculateGradient<TDim>(
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
        RansCalculationUtilities().CalculateGradient(
            rOutput, r_geometry, rVariable, rShapeDerivatives, Step);
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

    /**
     * @brief Calculates stabilization tau scalar derivatives
     *
     * \[
     *  -\tau_\phi^3\left[144\frac{\nu_\phi}{h^4_2}\left(\nu_{\phi,w}\right)^c + s_\phi\left(s_{\phi,w}\right)^c\right]
     * \]
     *
     * Where $w$ is the derivative variable
     *
     * @param rOutput                                        Scalar derivatives for each node w.r.t. $w$
     * @param Tau                                            Stabilization tau
     * @param EffectiveKinematicViscosity                    Effective kinematic viscosity $\nu_\phi$
     * @param Reaction                                       Reaction coefficient $s_\phi$
     * @param ElementLength                                  Element length $h_2$
     * @param rEffectiveKinematicViscosityScalarDerivatives  Scalar derivatives of effective kinematic viscosity $\left(\nu_{\phi,w}\right)$
     * @param rReactionScalarDerivatives                     Reaction scalar derivatives $\left(s_{\phi,w}\right)$
     */
    void CalculateStabilizationTauScalarDerivatives(Vector& rOutput,
                                                    const double Tau,
                                                    const double EffectiveKinematicViscosity,
                                                    const double Reaction,
                                                    const double ElementLength,
                                                    const Vector& rEffectiveKinematicViscosityScalarDerivatives,
                                                    const Vector& rReactionScalarDerivatives) const
    {
        noalias(rOutput) =
            (rEffectiveKinematicViscosityScalarDerivatives *
                 (144 * EffectiveKinematicViscosity / std::pow(ElementLength, 4)) +
             rReactionScalarDerivatives * (Reaction)) *
            (-1.0 * std::pow(Tau, 3));
    }

    void CalculateStabilizationTauVelocityDerivatives(
        Matrix& rOutput,
        const double Tau,
        const double EffectiveKinematicViscosity,
        const double Reaction,
        const double ElementLength,
        const array_1d<double, 3>& rVelocity,
        const Matrix& rContravariantMetricTensor,
        const Matrix& rEffectiveKinematicViscosityVelocityDerivatives,
        const Matrix& rReactionVelocityDerivatives,
        const Matrix& rElementLengthDerivatives,
        const Vector& rGaussShapeFunctions) const
    {
        Vector contravariant_metric_velocity(TDim);
        const Vector& velocity = RansCalculationUtilities().GetVector<TDim>(rVelocity);

        noalias(contravariant_metric_velocity) =
            prod(rContravariantMetricTensor, velocity) +
            prod(trans(rContravariantMetricTensor), velocity);

        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
            for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) = 0.5 * rGaussShapeFunctions[i_node] *
                                         contravariant_metric_velocity[i_dim];

        noalias(rOutput) +=
            rEffectiveKinematicViscosityVelocityDerivatives *
            (144.0 * EffectiveKinematicViscosity / std::pow(ElementLength, 4));
        noalias(rOutput) -= rElementLengthDerivatives *
                            (288.0 * std::pow(EffectiveKinematicViscosity, 2) /
                             std::pow(ElementLength, 5));
        noalias(rOutput) += rReactionVelocityDerivatives * (Reaction);
        noalias(rOutput) = rOutput * (-1.0 * std::pow(Tau, 3));
    }

    double CalculateStabilizationTauShapeSensitivity(const double Tau,
                                                     const double VelocityMagnitude,
                                                     const double ElementLength,
                                                     const double ElementLengthDeriv,
                                                     const double EffectiveKinematicViscosity,
                                                     const double EffectiveKinematicViscosityDeriv,
                                                     const double Reaction,
                                                     const double ReactionDeriv) const
    {
        double shape_sensitivity = 0.0;

        shape_sensitivity += 4.0 * std::pow(VelocityMagnitude, 2) *
                             ElementLengthDeriv / std::pow(ElementLength, 3);
        shape_sensitivity -= 144.0 * EffectiveKinematicViscosity *
                             EffectiveKinematicViscosityDeriv /
                             std ::pow(ElementLength, 4);
        shape_sensitivity += 288.0 * std::pow(EffectiveKinematicViscosity, 2) *
                             ElementLengthDeriv / std::pow(ElementLength, 5);
        shape_sensitivity -= Reaction * ReactionDeriv;

        shape_sensitivity *= std::pow(Tau, 3);

        return shape_sensitivity;
    }

    void CalculateVelocityMagnitudeVelocityDerivative(Matrix& rOutput,
                                                      const double VelocityMagnitude,
                                                      const array_1d<double, 3>& rVelocity,
                                                      const Vector& rGaussShapeFunctions) const
    {
        if (VelocityMagnitude <= std::numeric_limits<double>::epsilon())
        {
            rOutput.clear();
        }
        else
        {
            for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
                for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
                    rOutput(i_node, i_dim) =
                        rVelocity[i_dim] * rGaussShapeFunctions[i_node] / VelocityMagnitude;
        }
    }

    void CalculateElementLengthH2VelocityDerivative(Matrix& rOutput,
                                                    const double VelocityMagnitude,
                                                    const array_1d<double, 3>& rVelocity,
                                                    const Matrix& rVelocityMagnitudeVelocityDerivatives,
                                                    const Matrix& rContravariantMetricTensor,
                                                    const Vector& rGaussShapeFunctions) const
    {
        if (VelocityMagnitude <= std::numeric_limits<double>::epsilon())
        {
            rOutput.clear();
        }
        else
        {
            const Vector& velocity = RansCalculationUtilities().GetVector<TDim>(rVelocity);

            const double sqrt_u_e_u = std::sqrt(inner_prod(
                velocity, prod(rContravariantMetricTensor, velocity)));

            Vector contravariant_metric_velocity(TDim);
            noalias(contravariant_metric_velocity) =
                prod(rContravariantMetricTensor, velocity) +
                prod(trans(rContravariantMetricTensor), velocity);

            for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
                for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim)
                    rOutput(i_node, i_dim) = rGaussShapeFunctions[i_node] *
                                             contravariant_metric_velocity[i_dim];

            noalias(rOutput) =
                rOutput * (-1.0 * VelocityMagnitude / std::pow(sqrt_u_e_u, 3));
            noalias(rOutput) += (rVelocityMagnitudeVelocityDerivatives) * (2.0 / sqrt_u_e_u);
        }
    }

    double CalculateElementLengthH2ShapeSensitivity(const double VelocityMagnitude,
                                                    const array_1d<double, 3>& rVelocity,
                                                    const Matrix& rContravariantMetricTensor,
                                                    const Matrix& rContravariantMetricTensorShapeSensitivity)
    {
        if (VelocityMagnitude <= std::numeric_limits<double>::epsilon())
        {
            double sensitivity = 0.0;
            double element_length = 0.0;
            for (unsigned int i = 0; i < TDim; ++i)
                for (unsigned int j = 0; j < TDim; ++j)
                {
                    sensitivity += rContravariantMetricTensorShapeSensitivity(i, j);
                    element_length += rContravariantMetricTensor(i, j);
                }
            element_length = std::sqrt(1.0 / element_length) * 2.0;

            return sensitivity * std::pow(element_length, 3) * (-1.0 / 8.0);
        }
        else
        {
            const Vector& velocity = RansCalculationUtilities().GetVector<TDim>(rVelocity);

            const double u_e_u = std::pow(
                inner_prod(velocity, prod(rContravariantMetricTensor, velocity)), 1.5);

            return -VelocityMagnitude *
                   (inner_prod(velocity, prod(rContravariantMetricTensorShapeSensitivity,
                                              velocity))) /
                   u_e_u;
        }
    }

    void CalculateChiScalarDerivatives(Vector& rOutput,
                                       const double Chi,
                                       const double ElementLength,
                                       const double BossakAlpha,
                                       const double BossakGamma,
                                       const double DeltaTime,
                                       const double Reaction,
                                       const double DynamicTau,
                                       const Vector& rReactionScalarDerivatives)
    {
        const double reaction_tilde =
            Reaction + DynamicTau * (1 - BossakAlpha) / (BossakGamma * DeltaTime);

        CalculateAbsoluteScalarValueScalarDerivatives(
            rOutput, reaction_tilde, rReactionScalarDerivatives);
        noalias(rOutput) = rOutput * (-0.5 * std::pow(Chi, 2) * ElementLength);
    }

    void CalculateChiVelocityDerivatives(Matrix& rOutput,
                                         const double Chi,
                                         const double ElementLength,
                                         const double BossakAlpha,
                                         const double BossakGamma,
                                         const double DeltaTime,
                                         const double Reaction,
                                         const double DynamicTau,
                                         const Matrix& rReactionDerivatives,
                                         const Matrix& rVelocityMagnitudeDerivatives,
                                         const Matrix& rElementLengthDerivatives)
    {
        const double reaction_tilde =
            Reaction + DynamicTau * (1 - BossakAlpha) / (BossakGamma * DeltaTime);
        const double abs_reaction_tilde = std::abs(reaction_tilde);

        CalculateAbsoluteScalarValueVectorDerivatives(rOutput, reaction_tilde,
                                                      rReactionDerivatives);

        noalias(rOutput) = (rOutput * ElementLength + rElementLengthDerivatives * abs_reaction_tilde +
                            rVelocityMagnitudeDerivatives * 2.0) *
                           (-0.5 * std::pow(Chi, 2));
    }

    double CalculateChiShapeSensitivity(const double chi,
                                        const double reaction,
                                        const double reaction_deriv,
                                        const double element_length,
                                        const double element_length_deriv,
                                        const double bossak_alpha,
                                        const double bossak_gamma,
                                        const double delta_time,
                                        const double DynamicTau)
    {
        const double reaction_tilde =
            reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double abs_reaction_tilde = std::abs(reaction_tilde);

        return -0.5 * std::pow(chi, 2) *
               (abs_reaction_tilde * element_length_deriv +
                reaction_tilde * element_length * reaction_deriv /
                    (abs_reaction_tilde + std::numeric_limits<double>::epsilon()));
    }

    void CalculateAbsoluteScalarGradientScalarDerivative(Vector& rOutput,
                                                         const array_1d<double, 3> rScalarGradient,
                                                         const Matrix& rShapeFunctionDerivatives)
    {
        const double scalar_gradient_norm = norm_2(rScalarGradient);

        if (scalar_gradient_norm <= std::numeric_limits<double>::epsilon())
        {
            rOutput.clear();
        }
        else
        {
            for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
            {
                const Vector& shape_function_gradient =
                    row(rShapeFunctionDerivatives, i_node);
                rOutput[i_node] =
                    CalculateScalarProduct(shape_function_gradient, rScalarGradient) /
                    scalar_gradient_norm;
            }
        }
    }

    double CalculateAbsoluteScalarGradientShapeSensitivity(const array_1d<double, 3>& rScalarGradient,
                                                           const Matrix& rShapeFunctionDerivShapeSensitivity,
                                                           const Vector& rNodalScalarValues)
    {
        const double scalar_gradient_norm = norm_2(rScalarGradient);

        if (scalar_gradient_norm <= std::numeric_limits<double>::epsilon())
        {
            return 0.0;
        }
        else
        {
            Vector scalar_gradient_shape_sensitivity(TDim);
            noalias(scalar_gradient_shape_sensitivity) =
                prod(trans(rShapeFunctionDerivShapeSensitivity), rNodalScalarValues);

            const Vector& scalar_gradient =
                RansCalculationUtilities().GetVector<TDim>(rScalarGradient);
            return inner_prod(scalar_gradient_shape_sensitivity, scalar_gradient) / scalar_gradient_norm;
        }
    }

    void CalculateResidualScalarDerivative(Vector& rOutput,
                                           const double scalar_value,
                                           const double reaction,
                                           const array_1d<double, 3>& rVelocity,
                                           const Vector& rReactionScalarDerivatives,
                                           const Vector& rSourceScalarDerivatives,
                                           const Vector& rShapeFunctions,
                                           const Matrix& rShapeFunctionDerivatives,
                                           const Variable<double>& rDerivativeVariable)
    {
        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
        {
            const Vector& shape_function_gradient = row(rShapeFunctionDerivatives, i_node);
            double value = 0.0;

            value += scalar_value * rReactionScalarDerivatives[i_node];
            value -= rSourceScalarDerivatives[i_node];

            if (this->GetPrimalVariable() == rDerivativeVariable)
            {
                value += reaction * rShapeFunctions[i_node];
                value += CalculateScalarProduct(shape_function_gradient, rVelocity);
            }

            rOutput[i_node] = value;
        }
    }

    void CalculateResidualVelocityDerivative(Matrix& rOutput,
                                             const double primal_variable_value,
                                             const array_1d<double, 3>& rPrimalVariableGradient,
                                             const Matrix& rReactionDerivatives,
                                             const Matrix& rSourceDerivatives,
                                             const Vector& rGaussShapeFunctions)
    {
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) =
                    rGaussShapeFunctions[i_node] * rPrimalVariableGradient[i_dim];

        noalias(rOutput) = rOutput + rReactionDerivatives * primal_variable_value - rSourceDerivatives;
    }

    double CalculateResidualShapeSensitivity(const double residual,
                                             const array_1d<double, 3>& rVelocity,
                                             const Matrix& rShapeFunctionDerivShapeSensitivity,
                                             const double scalar_value,
                                             const Vector& rNodalScalarValues,
                                             const double reaction_deriv,
                                             const double source_deriv)
    {
        const double abs_residual = std::abs(residual);

        if (abs_residual <= std::numeric_limits<double>::epsilon())
        {
            return 0.0;
        }
        else
        {
            const Vector& r_velocity =
                RansCalculationUtilities().GetVector<TDim>(rVelocity);
            Vector primal_variable_gradient_shape_sensitivity(TDim);
            noalias(primal_variable_gradient_shape_sensitivity) =
                prod(trans(rShapeFunctionDerivShapeSensitivity), rNodalScalarValues);

            return residual *
                   (inner_prod(r_velocity, primal_variable_gradient_shape_sensitivity) +
                    reaction_deriv * scalar_value - source_deriv) /
                   abs_residual;
        }
    }

    void CalculatePositivityPreservationCoefficientScalarDerivatives(
        Vector& rOutput,
        const double chi,
        const double residual,
        const double scalar_gradient_norm,
        const double velocity_norm_square,
        const Vector& rChiScalarDerivatives,
        const Vector& rAbsoluteResidualScalarDerivatives,
        const Vector& rAbsoluteScalarGradientScalarDerivative,
        const Variable<double>& rDerivativeVariable)
    {
        const double abs_residual = std::abs(residual);

        if (scalar_gradient_norm <= std::numeric_limits<double>::epsilon() ||
            velocity_norm_square <= std::numeric_limits<double>::epsilon())
        {
            rOutput.clear();
        }
        else
        {
            noalias(rOutput) = rAbsoluteResidualScalarDerivatives *
                               (chi / (velocity_norm_square * scalar_gradient_norm));
            noalias(rOutput) +=
                rChiScalarDerivatives *
                (abs_residual / (velocity_norm_square * scalar_gradient_norm));

            if (this->GetPrimalVariable() == rDerivativeVariable)
                noalias(rOutput) -=
                    rAbsoluteScalarGradientScalarDerivative *
                    (chi * abs_residual /
                     (std::pow(scalar_gradient_norm, 2) * velocity_norm_square));
        }
    }

    void CalculatePositivityPreservationCoefficientVelocityDerivatives(
        Matrix& rOutput,
        const double absolute_residual,
        const double primal_variable_gradient_norm,
        const double velocity_magnitude,
        const double chi,
        const Matrix& rChiDerivatives,
        const Matrix& rAbsoluteResidualDerivatives,
        const Matrix& rVelocityMagnitudeDerivatives)
    {
        const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

        if (primal_variable_gradient_norm <= std::numeric_limits<double>::epsilon() ||
            velocity_magnitude_square <= std::numeric_limits<double>::epsilon())
        {
            rOutput.clear();
        }

        else
        {
            noalias(rOutput) =
                (rVelocityMagnitudeDerivatives * (-2.0 * chi / velocity_magnitude) + rChiDerivatives) *
                    (absolute_residual /
                     (velocity_magnitude_square * primal_variable_gradient_norm)) +
                rAbsoluteResidualDerivatives *
                    (chi / (primal_variable_gradient_norm * velocity_magnitude_square));
        }
    }

    double CalculatePositivityPreservationCoefficientShapeSensitivity(
        const double chi,
        const double chi_deriv,
        const double abs_residual,
        const double abs_residual_deriv,
        const double velocity_magnitude_square,
        const double scalar_gradient_norm,
        const double scalar_gradient_norm_deriv)
    {
        if (scalar_gradient_norm <= std::numeric_limits<double>::epsilon() ||
            velocity_magnitude_square <= std::numeric_limits<double>::epsilon())
        {
            return 0.0;
        }
        else
        {
            return chi_deriv * abs_residual / (scalar_gradient_norm * velocity_magnitude_square) +
                   chi * abs_residual_deriv / (scalar_gradient_norm * velocity_magnitude_square) -
                   chi * abs_residual * scalar_gradient_norm_deriv /
                       (std::pow(scalar_gradient_norm, 2) * velocity_magnitude_square);
        }
    }

    void CalculatePsiOneScalarDerivatives(Vector& rOutput,
                                          const double velocity_norm,
                                          const double reaction_tilde,
                                          const double tau,
                                          const Vector& rTauScalarDerivatives,
                                          const Vector& rAbsoluteReactionTildeScalarDerivatives)
    {
        const double absolute_reaction_tilde = std::abs(reaction_tilde);

        noalias(rOutput) = rTauScalarDerivatives * (velocity_norm * absolute_reaction_tilde);
        noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * velocity_norm);
    }

    void CalculatePsiOneVelocityDerivatives(Matrix& rOutput,
                                            const double velocity_norm,
                                            const double reaction_tilde,
                                            const double tau,
                                            const Matrix& rTauDerivatives,
                                            const Matrix& rAbsoluteReactionTildeDerivatives,
                                            const Matrix& rVelocityMagnitudeDerivatives)
    {
        noalias(rOutput) = rVelocityMagnitudeDerivatives +
                           rTauDerivatives * (velocity_norm * reaction_tilde) +
                           rVelocityMagnitudeDerivatives * (tau * reaction_tilde) +
                           rAbsoluteReactionTildeDerivatives * (tau * velocity_norm);
    }

    double CalculatePsiOneShapeSensitivity(const double tau,
                                           const double tau_deriv,
                                           const double velocity_magnitude,
                                           const double reaction,
                                           const double reaction_deriv,
                                           const double bossak_alpha,
                                           const double bossak_gamma,
                                           const double delta_time,
                                           const double DynamicTau)
    {
        const double reaction_dynamics =
            reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double abs_reaction_dynamics = std::abs(reaction_dynamics);

        return tau_deriv * velocity_magnitude * abs_reaction_dynamics +
               tau * velocity_magnitude * reaction_dynamics * reaction_deriv /
                   (abs_reaction_dynamics + std::numeric_limits<double>::epsilon());
    }

    void CalculatePsiTwoScalarDerivatives(Vector& rOutput,
                                          const double element_length,
                                          const double tau,
                                          const double reaction_tilde,
                                          const Vector& rTauScalarDerivatives,
                                          const Vector& rReactionTildeDerivatives,
                                          const Vector& rAbsoluteReactionTildeScalarDerivatives)
    {
        const double absolute_reaction_tilde = std::abs(reaction_tilde);

        noalias(rOutput) = rReactionTildeDerivatives;
        noalias(rOutput) +=
            rTauScalarDerivatives * (reaction_tilde * absolute_reaction_tilde);
        noalias(rOutput) += rReactionTildeDerivatives * (tau * absolute_reaction_tilde);
        noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * reaction_tilde);
        noalias(rOutput) = rOutput * (std::pow(element_length, 2) / 6.0);
    }

    void CalculatePsiTwoVelocityDerivatives(Matrix& rOutput,
                                            const double reaction_tilde,
                                            const double tau,
                                            const double element_length,
                                            const Matrix& rTauDerivatives,
                                            const Matrix& rReactionTildeDerivatives,
                                            const Matrix& rAbsoluteReactionTildeDerivatives,
                                            const Matrix& rElementLengthDerivatives)
    {
        const double abs_reaction_tilde = std::abs(reaction_tilde);

        noalias(rOutput) =
            (rReactionTildeDerivatives + rTauDerivatives * (reaction_tilde * abs_reaction_tilde) +
             rReactionTildeDerivatives * (tau * abs_reaction_tilde) +
             rAbsoluteReactionTildeDerivatives * (tau * reaction_tilde)) *
                std::pow(element_length, 2) / 6.0 +
            rElementLengthDerivatives *
                (element_length *
                 (reaction_tilde + tau * reaction_tilde * abs_reaction_tilde) / 3.0);
    }

    double CalculatePsiTwoShapeSensitivity(const double psi_two,
                                           const double element_length,
                                           const double element_length_deriv,
                                           const double reaction,
                                           const double reaction_deriv,
                                           const double tau,
                                           const double tau_deriv,
                                           const double bossak_alpha,
                                           const double bossak_gamma,
                                           const double delta_time,
                                           const double DynamicTau)
    {
        double shape_sensitivity = 0.0;

        const double reaction_dynamics =
            reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double abs_reaction_dynamics = std::abs(reaction_dynamics);

        shape_sensitivity += reaction_deriv;
        shape_sensitivity += tau_deriv * reaction_dynamics * abs_reaction_dynamics;
        shape_sensitivity += tau * reaction_deriv * abs_reaction_dynamics;
        shape_sensitivity +=
            tau * reaction_dynamics * reaction_dynamics * reaction_deriv /
            (abs_reaction_dynamics + std::numeric_limits<double>::epsilon());

        shape_sensitivity *= std::pow(element_length, 2) / 6.0;

        shape_sensitivity += 2.0 * psi_two * element_length_deriv / element_length;

        return shape_sensitivity;
    }

    void CalculateStreamLineDiffusionCoeffScalarDerivatives(
        Vector& rOutput,
        const double element_length,
        const double tau,
        const double velocity_norm,
        const double reaction_tilde,
        const double psi_one,
        const double psi_two,
        const Vector& rPsiOneScalarDerivatives,
        const Vector& rPsiTwoScalarDerivatives,
        const Vector& rTauScalarDerivatives,
        const Vector& rReactionTildeScalarDerivatives,
        const Vector& rEffectiveViscosityScalarDerivatives)
    {
        noalias(rOutput) = rPsiOneScalarDerivatives;
        noalias(rOutput) -= rTauScalarDerivatives * (velocity_norm * reaction_tilde);
        noalias(rOutput) -= rReactionTildeScalarDerivatives * (tau * velocity_norm);

        const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
        noalias(rOutput) = rOutput * (0.5 * element_length * (coeff) / (std::abs(coeff)) +
                                      std::numeric_limits<double>::epsilon());

        noalias(rOutput) += rPsiTwoScalarDerivatives;
        noalias(rOutput) -= rEffectiveViscosityScalarDerivatives;
        noalias(rOutput) -= rTauScalarDerivatives * std::pow(velocity_norm, 2);
    }

    void CalculateStreamLineDiffusionCoeffVelocityDerivatives(
        Matrix& rOutput,
        const double element_length,
        const double tau,
        const double velocity_norm,
        const double reaction_tilde,
        const double psi_one,
        const double psi_two,
        const Matrix& rVelocityMagnitudeDerivatives,
        const Matrix& rPsiOneDerivatives,
        const Matrix& rPsiTwoDerivatives,
        const Matrix& rTauDerivatives,
        const Matrix& rReactionTildeDerivatives,
        const Matrix& rEffectiveViscosityDerivatives,
        const Matrix& rElementLengthDerivatives)
    {
        noalias(rOutput) =
            (rPsiOneDerivatives - rTauDerivatives * (velocity_norm * reaction_tilde) -
             rVelocityMagnitudeDerivatives * (tau * reaction_tilde) -
             rReactionTildeDerivatives * (tau * velocity_norm));

        const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
        noalias(rOutput) = rOutput * (0.5 * element_length * (coeff) / (std::abs(coeff)) +
                                      std::numeric_limits<double>::epsilon());

        noalias(rOutput) += rElementLengthDerivatives * (0.5 * std::abs(coeff));

        noalias(rOutput) += rPsiTwoDerivatives;
        noalias(rOutput) -= rEffectiveViscosityDerivatives;
        noalias(rOutput) -= rTauDerivatives * std::pow(velocity_norm, 2);
        noalias(rOutput) -= rVelocityMagnitudeDerivatives * (2.0 * tau * velocity_norm);
    }

    double CalculateStreamLineDiffusionCoeffShapeSensitivity(
        const double psi_one,
        const double psi_one_deriv,
        const double tau,
        const double tau_deriv,
        const double velocity_magnitude,
        const double reaction,
        const double reaction_deriv,
        const double element_length,
        const double element_length_deriv,
        const double effective_kinematic_viscosity_deriv,
        const double psi_two_deriv,
        const double bossak_alpha,
        const double bossak_gamma,
        const double delta_time,
        const double DynamicTau)
    {
        const double reaction_dynamics =
            reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double coeff = psi_one - tau * velocity_magnitude * reaction_dynamics;
        const double abs_coeff = std::abs(coeff);
        double shape_sensitivity = 0.0;

        shape_sensitivity += psi_one_deriv - tau_deriv * velocity_magnitude * reaction_dynamics -
                             tau * velocity_magnitude * reaction_deriv;
        shape_sensitivity *= 0.5 * coeff * element_length /
                             (abs_coeff + std::numeric_limits<double>::epsilon());
        shape_sensitivity += 0.5 * abs_coeff * element_length_deriv;
        shape_sensitivity -= effective_kinematic_viscosity_deriv +
                             tau_deriv * std::pow(velocity_magnitude, 2);
        shape_sensitivity += psi_two_deriv;

        return shape_sensitivity;
    }

    void CalculateCrossWindDiffusionCoeffScalarDerivatives(
        Vector& rOutput,
        const double psi_one,
        const double element_length,
        const Vector& rPsiOneScalarDerivatives,
        const Vector& rPsiTwoScalarDerivatives,
        const Vector& rEffectiveKinematicViscosityScalarDerivatives)
    {
        noalias(rOutput) = rPsiOneScalarDerivatives *
                           (0.5 * psi_one * element_length / (std::abs(psi_one)) +
                            std::numeric_limits<double>::epsilon());
        noalias(rOutput) -= rEffectiveKinematicViscosityScalarDerivatives;
        noalias(rOutput) += rPsiTwoScalarDerivatives;
    }

    void CalculateCrossWindDiffusionCoeffVelocityDerivatives(
        Matrix& rOutput,
        const double psi_one,
        const double element_length,
        const Matrix& rPsiOneDerivatives,
        const Matrix& rPsiTwoDerivatives,
        const Matrix& rEffectiveKinematicViscosityDerivatives,
        const Matrix& rElementLengthDerivatives)
    {
        const double abs_psi_one = std::abs(psi_one);

        noalias(rOutput) =
            rPsiOneDerivatives * (0.5 * psi_one * element_length /
                                  (abs_psi_one + std::numeric_limits<double>::epsilon())) +
            rElementLengthDerivatives * (0.5 * abs_psi_one) -
            rEffectiveKinematicViscosityDerivatives + rPsiTwoDerivatives;
    }

    double CalculateCrossWindDiffusionCoeffShapeSensitivity(const double psi_one,
                                                            const double psi_one_deriv,
                                                            const double element_length,
                                                            const double element_length_deriv,
                                                            const double effective_kinematic_viscosity_deriv,
                                                            const double psi_two_deriv)
    {
        const double abs_psi_one = std::abs(psi_one);
        return 0.5 * psi_one * element_length * psi_one_deriv /
                   (abs_psi_one + std::numeric_limits<double>::epsilon()) +
               0.5 * abs_psi_one * element_length_deriv -
               effective_kinematic_viscosity_deriv + psi_two_deriv;
    }

    void CalculateAbsoluteScalarValueScalarDerivatives(Vector& rOutput,
                                                       const double scalar_value,
                                                       const Vector& rScalarValueDerivatives)
    {
        noalias(rOutput) = rScalarValueDerivatives *
                           (scalar_value / (std::abs(scalar_value) +
                                            std::numeric_limits<double>::epsilon()));
    }

    void CalculateAbsoluteScalarValueVectorDerivatives(Matrix& rOutput,
                                                       const double scalar_value,
                                                       const Matrix& rScalarValueDerivatives)
    {
        noalias(rOutput) = rScalarValueDerivatives *
                           (scalar_value / (std::abs(scalar_value) +
                                            std::numeric_limits<double>::epsilon()));
    }

    double CalculateScalarProduct(const Vector& rVector1, const array_1d<double, 3>& rVector2)
    {
        double result = 0.0;
        for (std::size_t i_dim = 0; i_dim < rVector1.size(); ++i_dim)
            result += rVector1[i_dim] * rVector2[i_dim];
        return result;
    }

    void AddPrimalSteadyTermScalarDerivatives(MatrixType& rLeftHandSideMatrix,
                                              const Variable<double>& rDerivativeVariable,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int TLocalMatrixSize = TMonolithicAssemblyLocalSize;

        if (rLeftHandSideMatrix.size1() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rLeftHandSideMatrix size mismatch. "
                            "rLeftHandSideMatrix.size1() != Expected size [ "
                         << rLeftHandSideMatrix.size1()
                         << " != " << TLocalMatrixSize << " ]";
        }

        if (rLeftHandSideMatrix.size2() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rLeftHandSideMatrix size mismatch. "
                            "rLeftHandSideMatrix.size2() != Expected size [ "
                         << rLeftHandSideMatrix.size2()
                         << " != " << TLocalMatrixSize << " ]";
        }

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const unsigned int equation_dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[primal_variable]);
        const unsigned int derivative_dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[rDerivativeVariable]);

        BoundedVector<double, TNumNodes> velocity_convective_terms, scalar_convective_terms;

        Matrix contravariant_metric_tensor(TDim, TDim);

        Vector effective_kinematic_viscosity_derivatives(TNumNodes),
            reaction_derivatives(TNumNodes), source_derivatives(TNumNodes),
            tau_derivatives(TNumNodes), s_derivatives(TNumNodes),
            chi_derivatives(TNumNodes), scalar_gradient_norm_derivative(TNumNodes),
            residual_derivatives(TNumNodes), absolute_residual_derivatives(TNumNodes),
            positivity_preserving_coeff_derivatives(TNumNodes),
            absolute_reaction_tilde_derivatives(TNumNodes),
            psi_one_derivatives(TNumNodes), psi_two_derivatives(TNumNodes),
            streamline_diffusion_coeff_derivatives(TNumNodes),
            crosswind_diffusion_coeff_derivatives(TNumNodes);

        array_1d<double, 3> scalar_gradient;

        TElementData current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

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

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityScalarDerivatives(
                effective_kinematic_viscosity_derivatives, rDerivativeVariable,
                current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            this->CalculateReactionTermScalarDerivatives(
                reaction_derivatives, rDerivativeVariable, current_data, rCurrentProcessInfo);

            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);
            this->CalculateSourceTermScalarDerivatives(
                source_derivatives, rDerivativeVariable, current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);
            this->CalculateStabilizationTauScalarDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction, element_length,
                effective_kinematic_viscosity_derivatives, reaction_derivatives);

            const double s = std::abs(reaction);
            CalculateAbsoluteScalarValueScalarDerivatives(
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
            if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
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
            }

            this->CalculateChiScalarDerivatives(
                chi_derivatives, chi, element_length, bossak_alpha, bossak_gamma,
                delta_time, reaction, dynamic_tau, reaction_derivatives);

            this->CalculateAbsoluteScalarGradientScalarDerivative(
                scalar_gradient_norm_derivative, scalar_gradient, r_shape_derivatives);

            this->CalculateResidualScalarDerivative(
                residual_derivatives, scalar_value, reaction, velocity,
                reaction_derivatives, source_derivatives, gauss_shape_functions,
                r_shape_derivatives, rDerivativeVariable);

            this->CalculateAbsoluteScalarValueScalarDerivatives(
                absolute_residual_derivatives, residual, residual_derivatives);

            this->CalculatePositivityPreservationCoefficientScalarDerivatives(
                positivity_preserving_coeff_derivatives, chi, residual,
                scalar_gradient_norm, velocity_magnitude_square,
                chi_derivatives, absolute_residual_derivatives,
                scalar_gradient_norm_derivative, rDerivativeVariable);

            const double reaction_tilde =
                reaction + dynamic_tau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
            this->CalculateAbsoluteScalarValueScalarDerivatives(
                absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

            const double psi_one =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
                    velocity_magnitude, tau, reaction_tilde);
            this->CalculatePsiOneScalarDerivatives(
                psi_one_derivatives, velocity_magnitude, reaction_tilde, tau,
                tau_derivatives, absolute_reaction_tilde_derivatives);

            const double psi_two =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
                    reaction_tilde, tau, element_length);
            this->CalculatePsiTwoScalarDerivatives(
                psi_two_derivatives, element_length, tau, reaction_tilde, tau_derivatives,
                reaction_derivatives, absolute_reaction_tilde_derivatives);

            this->CalculateStreamLineDiffusionCoeffScalarDerivatives(
                streamline_diffusion_coeff_derivatives, element_length, tau,
                velocity_magnitude, reaction_tilde, psi_one, psi_two,
                psi_one_derivatives, psi_two_derivatives, tau_derivatives,
                reaction_derivatives, effective_kinematic_viscosity_derivatives);

            this->CalculateCrossWindDiffusionCoeffScalarDerivatives(
                crosswind_diffusion_coeff_derivatives, psi_one, element_length, psi_one_derivatives,
                psi_two_derivatives, effective_kinematic_viscosity_derivatives);

            // calculating primal damping matrix scalar derivatives
            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    double dNa_dNc = 0.0;
                    for (unsigned int i = 0; i < TDim; i++)
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
                    rLeftHandSideMatrix(
                        c * TMonolithicAssemblyNodalDofSize + derivative_dof_index,
                        a * TMonolithicAssemblyNodalDofSize + equation_dof_index) +=
                        -1.0 * gauss_weights[g] * value;
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
                    rLeftHandSideMatrix(
                        c * TMonolithicAssemblyNodalDofSize + derivative_dof_index,
                        a * TMonolithicAssemblyNodalDofSize + equation_dof_index) +=
                        gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AddPrimalSteadyTermScalarRateDerivatives(MatrixType& rLeftHandSideMatrix,
                                                  const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int TLocalMatrixSize = TMonolithicAssemblyLocalSize;

        if (rLeftHandSideMatrix.size1() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rLeftHandSideMatrix size mismatch. "
                            "rLeftHandSideMatrix.size1() != Expected size [ "
                         << rLeftHandSideMatrix.size1()
                         << " != " << TLocalMatrixSize << " ]";
        }

        if (rLeftHandSideMatrix.size2() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rLeftHandSideMatrix size mismatch. "
                            "rLeftHandSideMatrix.size2() != Expected size [ "
                         << rLeftHandSideMatrix.size2()
                         << " != " << TLocalMatrixSize << " ]";
        }

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const unsigned int dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[primal_variable]);

        BoundedVector<double, TNumNodes> velocity_convective_terms, scalar_convective_terms;

        Matrix contravariant_metric_tensor(TDim, TDim);

        Vector positivity_preserving_coeff_derivatives(TNumNodes);

        array_1d<double, 3> scalar_gradient;

        TElementData current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

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

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);

            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);

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
            if (scalar_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                residual = relaxed_scalar_rate;
                residual += velocity_dot_scalar_gradient;
                residual += reaction * scalar_value;
                residual -= source;

                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);

                noalias(positivity_preserving_coeff_derivatives) =
                    gauss_shape_functions *
                    ((residual / (std::abs(residual) + std::numeric_limits<double>::epsilon())) *
                     chi / (scalar_gradient_norm * velocity_magnitude_square));
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
                    rLeftHandSideMatrix(c * TMonolithicAssemblyNodalDofSize + dof_index,
                                        a * TMonolithicAssemblyNodalDofSize + dof_index) +=
                        -1.0 * gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AddPrimalDampingMatrix(MatrixType& rPrimalDampingMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int TLocalMatrixSize = TMonolithicAssemblyLocalSize;

        if (rPrimalDampingMatrix.size1() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rPrimalDampingMatrix size mismatch. "
                            "rPrimalDampingMatrix.size1() != Expected size [ "
                         << rPrimalDampingMatrix.size1()
                         << " != " << TLocalMatrixSize << " ]";
        }

        if (rPrimalDampingMatrix.size2() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rPrimalDampingMatrix size mismatch. "
                            "rPrimalDampingMatrix.size2() != Expected size [ "
                         << rPrimalDampingMatrix.size2()
                         << " != " << TLocalMatrixSize << " ]";
        }

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const unsigned int num_gauss_points = gauss_weights.size();

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const unsigned int dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[primal_variable]);

        BoundedVector<double, TNumNodes> velocity_convective_terms;

        Matrix contravariant_metric_tensor(TDim, TDim);

        array_1d<double, 3> variable_gradient;

        TElementData r_current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];

            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);
            const double velocity_magnitude = norm_2(velocity);

            this->CalculateElementData(r_current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            this->CalculateGradient(variable_gradient, primal_variable, r_shape_derivatives);
            const double variable_gradient_norm = norm_2(variable_gradient);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(r_current_data, rCurrentProcessInfo);

            const double relaxed_variable_acceleration = this->EvaluateInPoint(
                this->GetPrimalRelaxedRateVariable(), gauss_shape_functions);

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

            const double velocity_dot_variable_gradient =
                inner_prod(velocity, variable_gradient);
            const double variable_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);

            if (variable_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
            {
                const double source =
                    this->CalculateSourceTerm(r_current_data, rCurrentProcessInfo);

                double residual = relaxed_variable_acceleration;
                residual += velocity_dot_variable_gradient;
                residual += reaction * variable_value;
                residual -= source;
                residual = std::abs(residual);
                residual /= variable_gradient_norm;

                double chi, k1, k2;
                StabilizedConvectionDiffusionReactionUtilities::CalculateCrossWindDiffusionParameters(
                    chi, k1, k2, velocity_magnitude, tau,
                    effective_kinematic_viscosity, reaction, bossak_alpha,
                    bossak_gamma, delta_time, element_length, dynamic_tau);

                stream_line_diffusion = residual * chi * k1 / velocity_magnitude_square;
                cross_wind_diffusion = residual * chi * k2 / velocity_magnitude_square;
            }

            const double s = std::abs(reaction);

            for (unsigned int a = 0; a < TNumNodes; a++)
            {
                for (unsigned int c = 0; c < TNumNodes; c++)
                {
                    double dNa_dNc = 0.0;
                    for (unsigned int i = 0; i < TDim; i++)
                        dNa_dNc += r_shape_derivatives(a, i) * r_shape_derivatives(c, i);

                    double value = 0.0;

                    value += gauss_shape_functions[a] * velocity_convective_terms[c];
                    value += gauss_shape_functions[a] * reaction *
                             gauss_shape_functions[c]; // * positive_values_list[c];
                    value += effective_kinematic_viscosity * dNa_dNc;

                    // Adding SUPG stabilization terms
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             velocity_convective_terms[c];
                    value += tau *
                             (velocity_convective_terms[a] + s * gauss_shape_functions[a]) *
                             reaction * gauss_shape_functions[c]; // * positive_values_list[c];

                    // Adding cross wind dissipation
                    value += cross_wind_diffusion * dNa_dNc * velocity_magnitude_square;
                    value -= cross_wind_diffusion * velocity_convective_terms[a] *
                             velocity_convective_terms[c];

                    // Adding stream line dissipation
                    value += stream_line_diffusion * velocity_convective_terms[a] *
                             velocity_convective_terms[c];

                    rPrimalDampingMatrix(c * TMonolithicAssemblyNodalDofSize + dof_index,
                                         a * TMonolithicAssemblyNodalDofSize + dof_index) +=
                        -1.0 * gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void AddLumpedMassMatrix(MatrixType& rMassMatrix, const double Mass, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int TLocalMatrixSize = TMonolithicAssemblyLocalSize;

        if (rMassMatrix.size1() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rMassMatrix size mismatch. "
                            "rMassMatrix.size1() != Expected size [ "
                         << rMassMatrix.size1() << " != " << TLocalMatrixSize << " ]";
        }

        if (rMassMatrix.size2() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rMassMatrix size mismatch. "
                            "rMassMatrix.size2() != Expected size [ "
                         << rMassMatrix.size2() << " != " << TLocalMatrixSize << " ]";
        }

        const unsigned int dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[this->GetPrimalVariable()]);

        for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
            rMassMatrix(iNode * TMonolithicAssemblyNodalDofSize + dof_index,
                        iNode * TMonolithicAssemblyNodalDofSize + dof_index) += Mass;

        KRATOS_CATCH("");
    }

    void AddPrimalMassMatrix(MatrixType& rMassMatrix,
                             const double ScalingFactor,
                             const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int TLocalMatrixSize = TMonolithicAssemblyLocalSize;

        if (rMassMatrix.size1() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rMassMatrix size mismatch. "
                            "rMassMatrix.size1() != Expected size [ "
                         << rMassMatrix.size1() << " != " << TLocalMatrixSize << " ]";
        }

        if (rMassMatrix.size2() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rMassMatrix size mismatch. "
                            "rMassMatrix.size2() != Expected size [ "
                         << rMassMatrix.size2() << " != " << TLocalMatrixSize << " ]";
        }

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();
        const unsigned int num_gauss_points = gauss_weights.size();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const unsigned int dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[this->GetPrimalVariable()]);

        TElementData current_data;

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
            this->AddLumpedMassMatrix(rMassMatrix, mass * ScalingFactor, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            BoundedVector<double, TNumNodes> velocity_convective_terms;
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            this->CalculateElementData(current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            const double s = std::abs(reaction);

            // Add mass stabilization terms
            for (unsigned int i = 0; i < TNumNodes; ++i)
                for (unsigned int j = 0; j < TNumNodes; ++j)
                    rMassMatrix(j * TMonolithicAssemblyNodalDofSize + dof_index,
                                i * TMonolithicAssemblyNodalDofSize + dof_index) +=
                        ScalingFactor * gauss_weights[g] * tau *
                        (velocity_convective_terms[i] + s * gauss_shape_functions[i]) *
                        gauss_shape_functions[j];
        }

        KRATOS_CATCH("");
    }

    void AddMassTermScalarDerivatives(MatrixType& rLeftHandSideMatrix,
                                      const Variable<double>& rDerivativeVariable,
                                      const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int TLocalMatrixSize = TMonolithicAssemblyLocalSize;

        if (rLeftHandSideMatrix.size1() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rLeftHandSideMatrix size mismatch. "
                            "rLeftHandSideMatrix.size1() != Expected size [ "
                         << rLeftHandSideMatrix.size1()
                         << " != " << TLocalMatrixSize << " ]";
        }

        if (rLeftHandSideMatrix.size2() != TLocalMatrixSize)
        {
            KRATOS_ERROR << "rLeftHandSideMatrix size mismatch. "
                            "rLeftHandSideMatrix.size2() != Expected size [ "
                         << rLeftHandSideMatrix.size2()
                         << " != " << TLocalMatrixSize << " ]";
        }

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const unsigned int equation_dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[this->GetPrimalVariable()]);
        const unsigned int derivative_dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[rDerivativeVariable]);

        Matrix contravariant_metric_tensor(TDim, TDim);

        BoundedVector<double, TNumNodes> velocity_convective_terms;

        TElementData current_data;

        Vector effective_kinematic_viscosity_derivatives(TNumNodes),
            reaction_derivatives(TNumNodes), tau_derivatives(TNumNodes),
            s_derivatives(TNumNodes);

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];

            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            this->CalculateElementData(current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityScalarDerivatives(
                effective_kinematic_viscosity_derivatives, rDerivativeVariable,
                current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            this->CalculateReactionTermScalarDerivatives(
                reaction_derivatives, rDerivativeVariable, current_data, rCurrentProcessInfo);

            const double relaxed_scalar_rate = this->EvaluateInPoint(
                this->GetPrimalRelaxedRateVariable(), gauss_shape_functions);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);
            this->CalculateStabilizationTauScalarDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction, element_length,
                effective_kinematic_viscosity_derivatives, reaction_derivatives);

            const double s = std::abs(reaction);
            this->CalculateAbsoluteScalarValueScalarDerivatives(
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

                    rLeftHandSideMatrix(
                        c * TMonolithicAssemblyNodalDofSize + derivative_dof_index,
                        a * TMonolithicAssemblyNodalDofSize + equation_dof_index) +=
                        -1.0 * gauss_weights[g] * value;
                }
            }
        }

        KRATOS_CATCH("");
    }

    void CalculateElementTotalResidualShapeSensitivity(Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int TLocalCoordsSize = TNumNodes * TDim;
        constexpr unsigned int TLocalMatrixSize = TMonolithicAssemblyLocalSize;

        if (!TMonolithicMatrixConstruction)
        {
            if (rOutput.size1() != TLocalCoordsSize || rOutput.size2() != TLocalMatrixSize)
                rOutput.resize(TLocalCoordsSize, TLocalMatrixSize, false);
            rOutput.clear();
        }

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

        RansCalculationUtilities rans_calculation_utilities;
        RansVariableUtils rans_variable_utils;

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const unsigned int dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[primal_variable]);

        const GeometryType::ShapeFunctionsGradientsType& r_dn_de =
            this->GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());

        Matrix contravariant_metric_tensor(TDim, TDim),
            parameter_derivatives_shape_derivs(TDim, TDim),
            contravariant_metric_tensor_deriv(TDim, TDim);

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

        for (unsigned int g = 0; g < num_gauss_points; g++)
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

            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

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

            rans_variable_utils.GetNodalArray(
                primal_variable_relaxed_rate_nodal_values, *this,
                this->GetPrimalRelaxedRateVariable());

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            double chi{0.0}, stream_line_diffusion_coeff{0.0},
                cross_wind_diffusion_coeff{0.0}, residual{0.0},
                positivity_preserving_coeff{0.0};

            if (primal_variable_gradient_norm > std::numeric_limits<double>::epsilon() &&
                velocity_magnitude_square > std::numeric_limits<double>::epsilon())
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

            rans_variable_utils.GetNodalArray(primal_variable_nodal_values,
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
                    deriv.NodeIndex = c;
                    deriv.Direction = k;

                    double detJ_deriv;
                    geom_sensitivity.CalculateSensitivity(deriv, detJ_deriv, DN_DX_deriv);
                    const double gauss_weight_deriv =
                        detJ_deriv * inv_detJ * gauss_weights[g];

                    rans_calculation_utilities.CalculateGeometryParameterDerivativesShapeSensitivity(
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
                        current_data, deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);
                    const double effective_kinematic_viscosity_deriv =
                        this->CalculateEffectiveKinematicViscosityShapeSensitivity(
                            current_data, deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);

                    this->GetConvectionOperator(
                        convective_primal_variable_gradient_terms_deriv,
                        primal_variable_gradient_deriv, r_shape_derivatives);
                    this->GetConvectionOperator(convective_deriv_primal_variable_gradient_terms,
                                                primal_variable_gradient, DN_DX_deriv);

                    const double element_length_deriv =
                        this->CalculateElementLengthH2ShapeSensitivity(
                            velocity_magnitude, velocity, contravariant_metric_tensor,
                            contravariant_metric_tensor_deriv);

                    const double tau_deriv = this->CalculateStabilizationTauShapeSensitivity(
                        tau, velocity_magnitude, element_length,
                        element_length_deriv, effective_kinematic_viscosity,
                        effective_kinematic_viscosity_deriv, reaction, reaction_deriv);

                    const double chi_deriv = this->CalculateChiShapeSensitivity(
                        chi, reaction, reaction_deriv, element_length, element_length_deriv,
                        bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

                    const double primal_variable_gradient_norm_deriv =
                        this->CalculateAbsoluteScalarGradientShapeSensitivity(
                            primal_variable_gradient, DN_DX_deriv, primal_variable_nodal_values);

                    const double source_deriv = this->CalculateSourceTermShapeSensitivity(
                        current_data, deriv, detJ_deriv, DN_DX_deriv, rCurrentProcessInfo);

                    const double residual_deriv = this->CalculateResidualShapeSensitivity(
                        residual, velocity, DN_DX_deriv, primal_variable_value,
                        primal_variable_nodal_values, reaction_deriv, source_deriv);

                    const double positivity_preserving_coeff_deriv =
                        CalculatePositivityPreservationCoefficientShapeSensitivity(
                            chi, chi_deriv, residual, residual_deriv,
                            velocity_magnitude_square, primal_variable_gradient_norm,
                            primal_variable_gradient_norm_deriv);

                    const double psi_one_deriv = CalculatePsiOneShapeSensitivity(
                        tau, tau_deriv, velocity_magnitude, reaction, reaction_deriv,
                        bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

                    const double psi_two_deriv = CalculatePsiTwoShapeSensitivity(
                        psi_two, element_length, element_length_deriv, reaction,
                        reaction_deriv, tau, tau_deriv, bossak_alpha,
                        bossak_gamma, delta_time, dynamic_tau);

                    const double stream_line_diffusion_coeff_deriv =
                        CalculateStreamLineDiffusionCoeffShapeSensitivity(
                            psi_one, psi_one_deriv, tau, tau_deriv, velocity_magnitude,
                            reaction, reaction_deriv, element_length, element_length_deriv,
                            effective_kinematic_viscosity_deriv, psi_two_deriv,
                            bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

                    const double cross_wind_diffusion_coeff_deriv =
                        CalculateCrossWindDiffusionCoeffShapeSensitivity(
                            psi_one, psi_one_deriv, element_length, element_length_deriv,
                            effective_kinematic_viscosity_deriv, psi_two_deriv);

                    const double s_deriv = reaction * reaction_deriv / s;

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

                        rOutput(block_size + k, a * TMonolithicAssemblyNodalDofSize +
                                                    dof_index) -= value;
                    }
                }
            }
        }

        KRATOS_CATCH("");
    }

    void CalculateElementTotalResidualVelocityDerivatives(Matrix& rResidualDerivatives,
                                                          const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        constexpr unsigned int TLocalMatrixSize = TMonolithicAssemblyLocalSize;
        const unsigned int local_vel_pr_size = std::max(TVelPrLocalSize, TLocalMatrixSize);

        if (!TMonolithicMatrixConstruction)
        {
            if (rResidualDerivatives.size1() != local_vel_pr_size ||
                rResidualDerivatives.size2() != TLocalMatrixSize)
                rResidualDerivatives.resize(local_vel_pr_size, TLocalMatrixSize, false);
            rResidualDerivatives.clear();
        }

        // Get Shape function data
        Vector gauss_weights;
        Matrix shape_functions;
        ShapeFunctionDerivativesArrayType shape_derivatives;
        this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
        const unsigned int num_gauss_points = gauss_weights.size();

        const ShapeFunctionDerivativesArrayType& r_parameter_derivatives =
            this->GetGeometryParameterDerivatives();

        const Variable<double>& primal_variable = this->GetPrimalVariable();

        const double delta_time = -1.0 * rCurrentProcessInfo[DELTA_TIME];
        const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
        const double bossak_gamma =
            TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
        const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
        const unsigned int equation_dof_index =
            static_cast<unsigned int>(rCurrentProcessInfo[primal_variable]);
        const unsigned int nodal_vel_pr_derivative_dof_size =
            std::max(TMonolithicAssemblyNodalDofSize, TVelPrBlockSize);

        BoundedVector<double, TNumNodes> primal_variable_gradient_convective_terms,
            velocity_convective_terms;

        Matrix effective_kinematic_viscosity_derivatives(TNumNodes, TDim),
            reaction_derivatives(TNumNodes, TDim),
            velocity_magnitude_derivatives(TNumNodes, TDim),
            element_length_derivatives(TNumNodes, TDim),
            tau_derivatives(TNumNodes, TDim), source_derivatives(TNumNodes, TDim),
            chi_derivatives(TNumNodes, TDim), residual_derivatives(TNumNodes, TDim),
            absolute_residual_derivatives(TNumNodes, TDim),
            positivity_preservation_coeff_derivatives(TNumNodes, TDim),
            absolute_reaction_tilde_derivatives(TNumNodes, TDim),
            psi_one_derivatives(TNumNodes, TDim), psi_two_derivatives(TNumNodes, TDim),
            stream_line_diffusion_coeff_derivatives(TNumNodes, TDim),
            cross_wind_diffusion_coeff_derivatives(TNumNodes, TDim),
            s_derivatives(TNumNodes, TDim), contravariant_metric_tensor(TDim, TDim);

        array_1d<double, 3> primal_variable_gradient;

        TElementData current_data;

        for (unsigned int g = 0; g < num_gauss_points; g++)
        {
            const Matrix& r_shape_derivatives = shape_derivatives[g];
            const Vector& gauss_shape_functions = row(shape_functions, g);

            const Matrix& r_parameter_derivatives_g = r_parameter_derivatives[g];
            noalias(contravariant_metric_tensor) =
                prod(trans(r_parameter_derivatives_g), r_parameter_derivatives_g);

            this->CalculateElementData(current_data, gauss_shape_functions,
                                       r_shape_derivatives, rCurrentProcessInfo);

            const double primal_variable_value =
                this->EvaluateInPoint(primal_variable, gauss_shape_functions);
            this->CalculateGradient(primal_variable_gradient, primal_variable,
                                    r_shape_derivatives);
            this->GetConvectionOperator(primal_variable_gradient_convective_terms,
                                        primal_variable_gradient, r_shape_derivatives);

            const double effective_kinematic_viscosity =
                this->CalculateEffectiveKinematicViscosity(current_data, rCurrentProcessInfo);
            this->CalculateEffectiveKinematicViscosityVelocityDerivatives(
                effective_kinematic_viscosity_derivatives, current_data, rCurrentProcessInfo);

            const double reaction =
                this->CalculateReactionTerm(current_data, rCurrentProcessInfo);
            this->CalculateReactionTermVelocityDerivatives(
                reaction_derivatives, current_data, rCurrentProcessInfo);

            const array_1d<double, 3>& velocity =
                this->EvaluateInPoint(VELOCITY, gauss_shape_functions);
            const double velocity_magnitude = norm_2(velocity);
            this->CalculateVelocityMagnitudeVelocityDerivative(
                velocity_magnitude_derivatives, velocity_magnitude, velocity,
                gauss_shape_functions);
            this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

            double tau, element_length;
            StabilizedConvectionDiffusionReactionUtilities::CalculateStabilizationTau(
                tau, element_length, velocity, contravariant_metric_tensor,
                reaction, effective_kinematic_viscosity, bossak_alpha,
                bossak_gamma, delta_time, dynamic_tau);

            this->CalculateElementLengthH2VelocityDerivative(
                element_length_derivatives, velocity_magnitude, velocity,
                velocity_magnitude_derivatives, contravariant_metric_tensor,
                gauss_shape_functions);

            this->CalculateStabilizationTauVelocityDerivatives(
                tau_derivatives, tau, effective_kinematic_viscosity, reaction,
                element_length, velocity, contravariant_metric_tensor,
                effective_kinematic_viscosity_derivatives, reaction_derivatives,
                element_length_derivatives, gauss_shape_functions);

            const double source =
                this->CalculateSourceTerm(current_data, rCurrentProcessInfo);
            this->CalculateSourceTermVelocityDerivatives(
                source_derivatives, current_data, rCurrentProcessInfo);

            const double velocity_dot_primal_variable_gradient =
                inner_prod(velocity, primal_variable_gradient);

            const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);
            const double primal_variable_gradient_norm = norm_2(primal_variable_gradient);

            double chi{0.0}, k1{0.0}, k2{0.0}, residual{0.0},
                positivity_preservation_coeff{0.0};

            const double primal_variable_relaxed_rate = this->EvaluateInPoint(
                this->GetPrimalRelaxedRateVariable(), gauss_shape_functions);

            if (velocity_magnitude_square > std::numeric_limits<double>::epsilon() &&
                primal_variable_gradient_norm > std::numeric_limits<double>::epsilon())
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
            }

            this->CalculateChiVelocityDerivatives(
                chi_derivatives, chi, element_length, bossak_alpha, bossak_gamma,
                delta_time, reaction, dynamic_tau, reaction_derivatives,
                velocity_magnitude_derivatives, element_length_derivatives);

            this->CalculateResidualVelocityDerivative(
                residual_derivatives, primal_variable_value, primal_variable_gradient,
                reaction_derivatives, source_derivatives, gauss_shape_functions);

            this->CalculateAbsoluteScalarValueVectorDerivatives(
                absolute_residual_derivatives, residual, residual_derivatives);

            this->CalculatePositivityPreservationCoefficientVelocityDerivatives(
                positivity_preservation_coeff_derivatives, std::abs(residual),
                primal_variable_gradient_norm, velocity_magnitude, chi, chi_derivatives,
                absolute_residual_derivatives, velocity_magnitude_derivatives);

            const double reaction_tilde =
                reaction + dynamic_tau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
            this->CalculateAbsoluteScalarValueVectorDerivatives(
                absolute_reaction_tilde_derivatives, reaction_tilde, reaction_derivatives);

            const double psi_one =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiOne(
                    velocity_magnitude, tau, reaction_tilde);
            this->CalculatePsiOneVelocityDerivatives(
                psi_one_derivatives, velocity_magnitude, reaction_tilde, tau, tau_derivatives,
                absolute_reaction_tilde_derivatives, velocity_magnitude_derivatives);

            const double psi_two =
                StabilizedConvectionDiffusionReactionUtilities::CalculatePsiTwo(
                    reaction_tilde, tau, element_length);
            this->CalculatePsiTwoVelocityDerivatives(
                psi_two_derivatives, reaction_tilde, tau, element_length,
                tau_derivatives, reaction_derivatives,
                absolute_reaction_tilde_derivatives, element_length_derivatives);

            this->CalculateStreamLineDiffusionCoeffVelocityDerivatives(
                stream_line_diffusion_coeff_derivatives, element_length, tau,
                velocity_magnitude, reaction_tilde, psi_one, psi_two,
                velocity_magnitude_derivatives, psi_one_derivatives,
                psi_two_derivatives, tau_derivatives, reaction_derivatives,
                effective_kinematic_viscosity_derivatives, element_length_derivatives);

            this->CalculateCrossWindDiffusionCoeffVelocityDerivatives(
                cross_wind_diffusion_coeff_derivatives, psi_one, element_length,
                psi_one_derivatives, psi_two_derivatives,
                effective_kinematic_viscosity_derivatives, element_length_derivatives);

            const double s = std::abs(reaction);
            this->CalculateAbsoluteScalarValueVectorDerivatives(
                s_derivatives, reaction, reaction_derivatives);

            for (unsigned int a = 0; a < TNumNodes; ++a)
            {
                for (unsigned int c = 0; c < TNumNodes; ++c)
                {
                    unsigned int block_size = c * nodal_vel_pr_derivative_dof_size;
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
                        rResidualDerivatives(block_size + k, a * TMonolithicAssemblyNodalDofSize +
                                                                 equation_dof_index) -=
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
                        rResidualDerivatives(block_size + k, a * TMonolithicAssemblyNodalDofSize +
                                                                 equation_dof_index) +=
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

                        rResidualDerivatives(block_size + k, a * TMonolithicAssemblyNodalDofSize +
                                                                 equation_dof_index) -=
                            value * gauss_weights[g];
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
template <unsigned int TDim, unsigned int TNumNodes, class TElementData, unsigned int TMonolithicAssemblyNodalDofSize = 1>
inline std::istream& operator>>(
    std::istream& rIStream,
    StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, TElementData, TMonolithicAssemblyNodalDofSize>& rThis);

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes, class TElementData, unsigned int TMonolithicAssemblyNodalDofSize = 1>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const StabilizedConvectionDiffusionReactionAdjointElement<TDim, TNumNodes, TElementData, TMonolithicAssemblyNodalDofSize>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}
///@}

} // namespace Kratos.
#endif // KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_ELEMENT defined
