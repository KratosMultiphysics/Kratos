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

#ifndef KRATOS_RANS_EVM_VMS_MONOLITHIC_ADJOINT_WALL_CONDITION_H
#define KRATOS_RANS_EVM_VMS_MONOLITHIC_ADJOINT_WALL_CONDITION_H

// System includes

// External includes

// Project includes
#include "includes/condition.h"

// Application includes

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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
 * @brief Epsilon Neumann wall condition
 *
 * This is a Neumann wall condition for $\epsilon$ equation in $k-\epsilon$
 * formulation of RANS based on eddy viscosity model formulation.
 *
 * @tparam TNumNodes Number of nodes in the wall condition
 */

template <unsigned int TDim>
class RansEvmVmsMonolithicAdjointWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansEvmVmsMonolithicAdjointWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RansEvmVmsMonolithicAdjointWallCondition);

    constexpr static unsigned int TNumNodes = TDim;

    constexpr static unsigned int TBlockSize = TDim + 1;

    constexpr static unsigned int TFluidLocalSize = TBlockSize * TNumNodes;

    constexpr static unsigned int TCoordLocalSize = TDim * TNumNodes;

    typedef Element::IndexType IndexType;

    typedef Element::SizeType SizeType;

    typedef Element::GeometryType GeometryType;

    typedef Element::PropertiesType PropertiesType;

    typedef Element::NodesArrayType NodesArrayType;

    typedef Element::VectorType VectorType;

    typedef std::array<double, TFluidLocalSize> ArrayType;

    typedef Element::MatrixType MatrixType;

    typedef Element::DofsVectorType DofsVectorType;

    typedef std::array<Dof<double>::Pointer, TFluidLocalSize> DofsArrayType;

    typedef Element::EquationIdVectorType EquationIdVectorType;

    typedef std::array<std::size_t, TFluidLocalSize> EquationIdArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    explicit RansEvmVmsMonolithicAdjointWallCondition(IndexType NewId = 0)
        : Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    RansEvmVmsMonolithicAdjointWallCondition(IndexType NewId, const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    RansEvmVmsMonolithicAdjointWallCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    RansEvmVmsMonolithicAdjointWallCondition(IndexType NewId,
                                             GeometryType::Pointer pGeometry,
                                             PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    RansEvmVmsMonolithicAdjointWallCondition(RansEvmVmsMonolithicAdjointWallCondition const& rOther)
        : Condition(rOther)
    {
    }

    /// Destructor.
    ~RansEvmVmsMonolithicAdjointWallCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    RansEvmVmsMonolithicAdjointWallCondition& operator=(RansEvmVmsMonolithicAdjointWallCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new RansEvmVmsMonolithicAdjointWallCondition object.
    /**
      @param NewId Index of the new condition
      @param ThisNodes An array containing the nodes of the new condition
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeom,
                              PropertiesType::Pointer pProperties) const override;

    /**
     * Clones the selected element variables, creating a new one
     * @param NewId the ID of the new element
     * @param ThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */

    Condition::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override;

    /// Return local contributions of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    /// Return a matrix of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see DampingMatrix
      */
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo);

    /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& rCurrentProcessInfo) override;

    /// Check that all data required by this condition is available and reasonable
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /// Provides the global indices for each one of this element's local rows.
    /** This determines the elemental equation ID vector for all elemental DOFs
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object (unused)
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdArray(EquationIdArrayType& rResult, ProcessInfo& rCurrentProcessInfo);

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofArray(DofsArrayType& rConditionDofList, ProcessInfo& rCurrentProcessInfo);

    void GetValuesVector(VectorType& rValues, int Step = 0) override;

    void GetValuesArray(ArrayType& rValues, int Step = 0);

    /// Returns TURBULENT_ENERGY_DISSIPATION_RATE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    void GetFirstDerivativesArray(ArrayType& rValues, int Step = 0);

    /// Returns TURBULENT_ENERGY_DISSIPATION_RATE_2 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

    void GetSecondDerivativesArray(ArrayType& rValues, int Step = 0);

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
                                      ProcessInfo& rCurrentProcessInfo) override;

    // Calculates epsilon wall condition derivative w.r.t. epsilon variable
    void CalculateFirstDerivativesLHS(BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rLeftHandSideMatrix,
                                      const ProcessInfo& rCurrentProcessInfo);

    void Calculate(const Variable<Matrix>& rVariable,
                   Matrix& rOutput,
                   const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateConditionResidualTurbulentKineticEnergyDerivatives(
        BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rOutput,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateConditionResidualTurbulentEnergyDissipationRateDerivatives(
        BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rOutput,
        const ProcessInfo& rCurrentProcessInfo);

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
                                       ProcessInfo& rCurrentProcessInfo) override;

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
                                    const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rSensitivityVariable,
                                    BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize>& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo);

    void CalculateResidualShapeSensitivity(BoundedMatrix<double, TCoordLocalSize, TFluidLocalSize>& rOutput,
                                           const ProcessInfo& rCurrentProcessInfo);

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
        buffer << "RansEvmVmsMonolithicAdjointWallCondition" << TNumNodes << "N";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RansEvmVmsMonolithicAdjointWallCondition";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
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
     * @brief Calculates scalar value for given gauss point
     *
     * @param rVariable      Scalar variable
     * @param rShapeFunction Gauss point shape functions
     * @param Step           Step
     * @return double        Gauss point scalar value
     */
    double EvaluateInPoint(const Variable<double>& rVariable,
                           const Vector& rShapeFunction,
                           const int Step = 0) const;

    array_1d<double, 3> EvaluateInPoint(const Variable<array_1d<double, 3>>& rVariable,
                                        const Vector& rShapeFunction,
                                        const int Step = 0) const;

    void CalculateNormal(array_1d<double, 3>& An);

    void CalculateUnitNormalShapeSensitivities(BoundedMatrix<double, TCoordLocalSize, TDim>& rOutput);

    void ApplyRansBasedWallLawFirstDerivatives(
        BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize>& rLocalMatrix,
        const ProcessInfo& rCurrentProcessInfo);

    void ApplyRansBasedWallLawTurbulentKineticEnergyDerivatives(
        BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rLocalMatrix,
        const ProcessInfo& rCurrentProcessInfo);

    void ApplyRansBasedWallLawTurbulentEnergyDissipationRateDerivatives(
        BoundedMatrix<double, TNumNodes, TFluidLocalSize>& rLocalMatrix,
        const ProcessInfo& rCurrentProcessInfo);

    void CalculateVelocityMagnitudeVelocityDerivative(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const double VelocityMagnitude,
        const array_1d<double, 3>& rVelocity,
        const Vector& rGaussShapeFunctions) const;

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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition);
    }

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    MatrixType GetJacobian(GeometryData::IntegrationMethod QuadratureOrder,
                           unsigned int IntegrationPointIndex) const;

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

}; // Class RansEvmVmsMonolithicAdjointWallCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                RansEvmVmsMonolithicAdjointWallCondition<TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEvmVmsMonolithicAdjointWallCondition<TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EVM_VMS_MONOLITHIC_ADJOINT_WALL_CONDITION_H
