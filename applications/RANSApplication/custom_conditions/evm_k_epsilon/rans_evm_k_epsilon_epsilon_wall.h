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

#ifndef KRATOS_RANS_EVM_K_EPSILON_EPSILON_WALL_H
#define KRATOS_RANS_EVM_K_EPSILON_EPSILON_WALL_H

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
 * @tparam TDim dimension of the wall condition
 * @tparam TNumNodes Number of nodes in the wall condition
 */

template <unsigned int TDim, unsigned int TNumNodes = TDim>
class RansEvmKEpsilonEpsilonWall : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of RansEvmKEpsilonEpsilonWall
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(RansEvmKEpsilonEpsilonWall);

    using NodeType = Node<3>;
    using PropertiesType = Properties;
    using GeometryType = Geometry<NodeType>;
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;
    using VectorType = Vector;
    using MatrixType = Matrix;
    using IndexType = std::size_t;
    using EquationIdVectorType = std::vector<IndexType>;
    using DofsVectorType = std::vector<Dof<double>::Pointer>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
      @param NewId Index for the new condition
      */
    explicit RansEvmKEpsilonEpsilonWall(IndexType NewId = 0) : Condition(NewId)
    {
    }

    /// Constructor using an array of nodes
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    RansEvmKEpsilonEpsilonWall(IndexType NewId, const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    RansEvmKEpsilonEpsilonWall(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties
    /**
     @param NewId Index of the new element
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the element's properties
     */
    RansEvmKEpsilonEpsilonWall(IndexType NewId,
                               GeometryType::Pointer pGeometry,
                               PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    RansEvmKEpsilonEpsilonWall(RansEvmKEpsilonEpsilonWall const& rOther)
        : Condition(rOther)
    {
    }

    /// Destructor.
    ~RansEvmKEpsilonEpsilonWall() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator
    RansEvmKEpsilonEpsilonWall& operator=(RansEvmKEpsilonEpsilonWall const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new RansEvmKEpsilonEpsilonWall object.
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

    /// Return local right hand side of the correct size, filled with zeros (for compatibility with time schemes).
    /** The actual local contributions are computed in the Damping functions
      @see CalculateLocalVelocityContribution
      */
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix,
                                ProcessInfo& rCurrentProcessInfo) override;

    /// Calculate wall stress term for all nodes with SLIP == true
    /**
      @param rDampingMatrix Left-hand side matrix
      @param rRightHandSideVector Right-hand side vector
      @param rCurrentProcessInfo ProcessInfo instance (unused)
      */
    void CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
                                            VectorType& rRightHandSideVector,
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

    /// Returns a list of the element's Dofs
    /**
     * @param ElementalDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& ConditionDofList, ProcessInfo& CurrentProcessInfo) override;

    void GetValuesVector(VectorType& rValues, int Step = 0) override;

    /// Returns TURBULENT_ENERGY_DISSIPATION_RATE for each node
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

    /// Returns TURBULENT_ENERGY_DISSIPATION_RATE_2 0 for each node
    /**
     * @param Values Vector of nodal second derivatives
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

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
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

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

    void AddLocalVelocityContribution(MatrixType& rDampingMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo);
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

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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

}; // Class RansEvmKEpsilonEpsilonWall

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                RansEvmKEpsilonEpsilonWall<TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const RansEvmKEpsilonEpsilonWall<TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_RANS_EVM_K_EPSILON_EPSILON_WALL_H
