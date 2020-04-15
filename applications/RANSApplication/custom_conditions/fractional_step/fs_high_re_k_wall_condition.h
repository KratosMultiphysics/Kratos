//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#ifndef KRATOS_FS_HIGH_RE_K_WALL_CONDITION
#define KRATOS_FS_HIGH_RE_K_WALL_CONDITION

// System includes

// External includes

// Project includes
#include "includes/condition.h"
#include "includes/serializer.h"

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

template <unsigned int TDim, unsigned int TNumNodes = TDim>
class FSHighReKWallCondition : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FSHighReKWallCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(FSHighReKWallCondition);

    using BaseType = Condition;
    using NodeType = Node<3>;
    using PropertiesType = Properties;
    using GeometryType = Geometry<NodeType>;
    using PointType = Geometry<NodeType>::PointType;
    using NodesArrayType = Geometry<NodeType>::PointsArrayType;
    using GeometriesArrayType = Geometry<NodeType>::GeometriesArrayType;
    using ElementWeakPointerType = Element::WeakPointer;
    using ElementPointerType = Element::Pointer;
    using IndexType = std::size_t;
    using SizeType = std::size_t;
    using EquationIdVectorType = std::vector<std::size_t>;
    using DofsVectorType = std::vector<Dof<double>::Pointer>;
    using ShapeFunctionsType = Vector;
    using ShapeFunctionDerivativesType = Matrix;
    using ShapeFunctionDerivativesArrayType = GeometryType::ShapeFunctionsGradientsType;

    constexpr static IndexType VelocityLocalSize = TDim * TNumNodes;
    constexpr static IndexType PressureLocalSize = TNumNodes;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    /** Admits an Id as a parameter.
     @param NewId Index of the new condition
     */
    FSHighReKWallCondition(IndexType NewId = 0) : Condition(NewId)
    {
    }

    /// Constructor using an array of nodes.
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     */
    FSHighReKWallCondition(IndexType NewId, const NodesArrayType& ThisNodes)
        : Condition(NewId, ThisNodes)
    {
    }

    /// Constructor using Geometry.
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     */
    FSHighReKWallCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry)
    {
    }

    /// Constructor using Properties.
    /**
     @param NewId Index of the new condition
     @param pGeometry Pointer to a geometry object
     @param pProperties Pointer to the condition's properties
     */
    FSHighReKWallCondition(IndexType NewId,
                           GeometryType::Pointer pGeometry,
                           PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties)
    {
    }

    /// Copy constructor.
    FSHighReKWallCondition(FSHighReKWallCondition const& rOther)
        : Condition(rOther)
    {
    }

    /// Destructor.
    ~FSHighReKWallCondition() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Copy constructor.
    FSHighReKWallCondition& operator=(FSHighReKWallCondition const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /// Create a new FSHighReKWallCondition object.
    /**
     @param NewId Index of the new condition
     @param ThisNodes An array containing the nodes of the new condition
     @param pProperties Pointer to the condition's properties
     */
    Condition::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FSHighReKWallCondition>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    /// Create a new FSHighReKWallCondition object.
    /**
      @param NewId Index of the new condition
      @param pGeom A pointer to the condition's geometry
      @param pProperties Pointer to the element's properties
      */
    Condition::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeom,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<FSHighReKWallCondition>(NewId, pGeom, pProperties);
    }

    /// Find the condition's parent element.
    void Initialize() override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override;

    /// Calculate wall stress term for all nodes with SLIP set.
    /**
     @param rLeftHandSideMatrix Left-hand side matrix
     @param rRightHandSideVector Right-hand side vector
     @param rCurrentProcessInfo ProcessInfo instance
     */
    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              ProcessInfo& rCurrentProcessInfo) override;

    /// Check that all data required by this condition is available and reasonable.
    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    /// Provides the global indices for each one of this condition's local rows.
    /** This determines the equation ID vector for all DOFs.
     * @param rResult A vector containing the global Id of each row
     * @param rCurrentProcessInfo the current process info object
     */
    void EquationIdVector(EquationIdVectorType& rResult,
                          ProcessInfo& rCurrentProcessInfo) override;

    /// Returns a list of the condition's Dofs.
    /**
     * @param ConditionDofList the list of DOFs
     * @param rCurrentProcessInfo the current process info instance
     */
    void GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo) override;

    /// Returns VELOCITY_X, VELOCITY_Y, (VELOCITY_Z) for each node.
    /**
     * @param Values Vector of nodal unknowns
     * @param Step Get result from 'Step' steps back, 0 is current step. (Must be smaller than buffer size)
     */
    void GetValuesVector(Vector& Values, int Step = 0) override;

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

    /// Compute the wall stress and add corresponding terms to the system contributions.
    /**
     @param rLocalMatrix Local system matrix
     @param rLocalVector Local right hand side
     */
    void ApplyWallLaw(MatrixType& rLocalMatrix, VectorType& rLocalVector, ProcessInfo& rCurrentProcessInfo);

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

    double mWallHeight;

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

}; // Class FSHighReKWallCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::istream& operator>>(std::istream& rIStream,
                                FSHighReKWallCondition<TDim, TNumNodes>& rThis)
{
    return rIStream;
}

/// output stream function
template <unsigned int TDim, unsigned int TNumNodes>
inline std::ostream& operator<<(std::ostream& rOStream,
                                const FSHighReKWallCondition<TDim, TNumNodes>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_FS_HIGH_RE_K_WALL_CONDITION
