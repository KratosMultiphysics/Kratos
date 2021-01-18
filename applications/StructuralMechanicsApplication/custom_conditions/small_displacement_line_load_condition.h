// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#if !defined(KRATOS_SMALL_DISPLACEMENT_LINE_LOAD_CONDITION_H_INCLUDED )
#define  KRATOS_SMALL_DISPLACEMENT_LINE_LOAD_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_conditions/line_load_condition.h"

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
 * @class SmallDisplacementLineLoadCondition
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the responsible to add the contributions of the RHS and LHS of the line loads of the structure
 * @details It allows to consider different types of pressure and line loads
 * @tparam TDim The dimension of the condition
 * @author Vicente Mataix Ferrandiz
 */
template<std::size_t TDim>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SmallDisplacementLineLoadCondition
    : public LineLoadCondition<TDim>
{
public:
    ///@name Type Definitions
    ///@{

    /// We define the base class Base LoadCondition
    typedef LineLoadCondition<TDim> BaseType;

    /// Definition of the vector type
    typedef typename BaseLoadCondition::VectorType VectorType;
    
    /// Definition of the matrix type
    typedef typename BaseLoadCondition::MatrixType MatrixType;
    
    /// Definition of the index type
    typedef typename BaseLoadCondition::IndexType IndexType;

    /// Definition of the size type
    typedef typename BaseLoadCondition::SizeType SizeType;

    /// Definition of the node type
    typedef typename BaseLoadCondition::NodeType NodeType;

    /// Definition of the properties type
    typedef typename BaseLoadCondition::PropertiesType PropertiesType;

    /// Definition of the geometry type with given NodeType
    typedef typename BaseLoadCondition::GeometryType GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef typename BaseLoadCondition::NodesArrayType NodesArrayType;

    /// Counted pointer of SmallDisplacementLineLoadCondition
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION( SmallDisplacementLineLoadCondition );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor using an array of nodes
    SmallDisplacementLineLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        );

    // Constructor using an array of nodes with properties
    SmallDisplacementLineLoadCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        );

    /// Destructor.
    ~SmallDisplacementLineLoadCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer
     * @param NewId the ID of the new condition
     * @param pGeom the geometry to be employed
     * @param pProperties the properties assigned to the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * @brief Creates a new condition pointer and clones the previous condition data
     * @param NewId the ID of the new condition
     * @param ThisNodes the nodes of the new condition
     * @return a Pointer to the new condition
     */
    Condition::Pointer Clone (
        IndexType NewId,
        NodesArrayType const& ThisNodes
        ) const override;

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
        buffer << "Small displacement line load condition #" << this->Id();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SmallDisplacementLineLoadCondition #" << this->Id();
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

    /**
     * @brief This functions calculates both the RHS and the LHS
     * @param rLeftHandSideMatrix: The LHS
     * @param rRightHandSideVector: The RHS
     * @param rCurrentProcessInfo: The current process info instance
     * @param CalculateStiffnessMatrixFlag: The flag to set if compute the LHS
     * @param CalculateResidualVectorFlag: The flag to set if compute the RHS
     */
    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
        ) override;

    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{

    // A protected default constructor necessary for serialization
    SmallDisplacementLineLoadCondition() {};

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
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType );
    }

    void load( Serializer& rSerializer ) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class SmallDisplacementLineLoadCondition

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
template<std::size_t TDim>
inline std::istream& operator >> (std::istream& rIStream,
        SmallDisplacementLineLoadCondition<TDim>& rThis);
/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream,
        const SmallDisplacementLineLoadCondition<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_SMALL_DISPLACEMENT_LINE_LOAD_CONDITION_H_INCLUDED  defined


