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

#ifndef KRATOS_RANS_EVM_EPSILON_ADJOINT_WALL_CONDITION_H
#define KRATOS_RANS_EVM_EPSILON_ADJOINT_WALL_CONDITION_H

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

template <unsigned int TNumNodes, unsigned int TDim = TNumNodes>
class RansEvmEpsilonAdjointWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansEvmEpsilonAdjointWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RansEvmEpsilonAdjointWallCondition);

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

    typedef BoundedMatrix<double, TNumNodes, TNumNodes> BoundedMatrixNN;

    typedef BoundedVector<double, TNumNodes> BoundedVectorN;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    explicit RansEvmEpsilonAdjointWallCondition(IndexType NewId = 0)
        : Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    RansEvmEpsilonAdjointWallCondition(IndexType NewId, const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    RansEvmEpsilonAdjointWallCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    RansEvmEpsilonAdjointWallCondition(IndexType NewId,
                                       GeometryType::Pointer pGeometry,
                                       PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    RansEvmEpsilonAdjointWallCondition(RansEvmEpsilonAdjointWallCondition const& rOther)
        : Condition(rOther)
    {
    }

    /// Destructor.
    ~RansEvmEpsilonAdjointWallCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    RansEvmEpsilonAdjointWallCondition& operator=(RansEvmEpsilonAdjointWallCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new RansEvmEpsilonAdjointWallCondition object.
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

    void EquationIdArray(std::array<std::size_t, TNumNodes>& rResult,
                         ProcessInfo& rCurrentProcessInfo);

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofArray(std::array<Dof<double>::Pointer, TNumNodes>& rConditionDofList,
                     ProcessInfo& rCurrentProcessInfo);

    void GetValuesVector(VectorType& rValues, int Step = 0) override;

    void GetValuesArray(std::array<double, TNumNodes>& rValues, int Step = 0);

    /// Returns TURBULENT_ENERGY_DISSIPATION_RATE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    void GetFirstDerivativesArray(std::array<double, TNumNodes>& rValues, int Step = 0);

    /// Returns TURBULENT_ENERGY_DISSIPATION_RATE_2 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

    void GetSecondDerivativesArray(std::array<double, TNumNodes>& rValues, int Step = 0);

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
    void CalculateFirstDerivativesLHS(BoundedMatrix<double, TNumNodes, TNumNodes>& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo);

    void Calculate(const Variable<Matrix>& rVariable,
                   Matrix& rOutput,
                   const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateConditionResidualTurbulentKineticEnergyDerivatives(
        BoundedMatrixNN& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void CalculateConditionResidualTurbulentEnergyDissipationRateDerivatives(
        BoundedMatrixNN& rOutput, const ProcessInfo& rCurrentProcessInfo);

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
                                    BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo);

    void CalculateResidualShapeSensitivity(BoundedMatrix<double, TNumNodes * TDim, TNumNodes>& rOutput,
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
        buffer << "RansEvmEpsilonAdjointWallCondition" << TNumNodes << "N";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "RansEvmEpsilonAdjointWallCondition";
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

}; // Class RansEvmEpsilonAdjointWallCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                RansEvmEpsilonAdjointWallCondition<TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEvmEpsilonAdjointWallCondition<TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EVM_EPSILON_ADJOINT_WALL_CONDITION_H
