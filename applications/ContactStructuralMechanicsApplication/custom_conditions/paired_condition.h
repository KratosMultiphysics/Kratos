// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_PAIRED_CONDITION_H_INCLUDED )
#define  KRATOS_PAIRED_CONDITION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/condition.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    typedef Point                                                PointType;
    typedef Node<3>                                               NodeType;
    typedef Geometry<NodeType>                                GeometryType;

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
 * @ingroup ContactStructuralMechanicsApplication
 * @class PairedCondition
 * @brief This is a base class for the conditions paired
 * @details This is a base class for the conditions paired, it is basically equal to the base condition, with a pointer to the paired geoemtry
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) PairedCondition
    : public Condition
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of PairedCondition
    KRATOS_CLASS_POINTER_DEFINITION( PairedCondition );

    typedef Condition                                                           BaseType;

    typedef BaseType::VectorType                                              VectorType;

    typedef BaseType::MatrixType                                              MatrixType;

    typedef BaseType::IndexType                                                IndexType;

    typedef BaseType::GeometryType::Pointer                          GeometryPointerType;

    typedef BaseType::NodesArrayType                                      NodesArrayType;

    typedef BaseType::PropertiesType::Pointer                      PropertiesPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    PairedCondition()
        : Condition(),
          mpPairedGeometry(nullptr)
    {}

    // Constructor 1
    PairedCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry
        ) :Condition(NewId, pGeometry),
           mpPairedGeometry(nullptr)
    {}

    // Constructor 2
    PairedCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
        ) :Condition( NewId, pGeometry, pProperties ),
           mpPairedGeometry(nullptr)
    {}

    // Constructor 3
    PairedCondition(
        IndexType NewId,
        GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        GeometryType::Pointer pPairedGeometry
        )
        :Condition( NewId, pGeometry, pProperties ),
         mpPairedGeometry(pPairedGeometry)
    {}

    ///Copy constructor
    PairedCondition( PairedCondition const& rOther){}

    /// Destructor.
    ~PairedCondition() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
    * Called at the beginning of each solution step
    */
    void Initialize() override;

    /**
     * Creates a new element pointer from an arry of nodes
     * @param NewId the ID of the new element
     * @param rThisNodes the nodes of the new element
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * Creates a new element pointer from an existing geometry
     * @param NewId the ID of the new element
     * @param pGeom the  geometry taken to create the condition
     * @param pProperties the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    /**
     * Creates a new element pointer from an existing geometry
     * @param NewId the ID of the new element
     * @param pGeom the  geometry taken to create the condition
     * @param pProperties the properties assigned to the new element
     * @param pPairedGeom the paired geometry
     * @return a Pointer to the new element
     */
    virtual Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties,
        GeometryType::Pointer pPairedGeom
        ) const;

    ///@}
    ///@name Access
    ///@{

    GeometryType::Pointer pGetPairedGeometry()
    {
        return mpPairedGeometry;
    }

    const GeometryType::Pointer pGetPairedGeometry() const
    {
        return mpPairedGeometry;
    }

    GeometryType& GetPairedGeometry()
    {
        return *mpPairedGeometry;
    }

    GeometryType const& GetPairedGeometry() const
    {
        return *mpPairedGeometry;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

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

    GeometryType::Pointer mpPairedGeometry; // The geometry of the pair "condition"

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

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Condition );
        rSerializer.save("PairedGeometry", mpPairedGeometry);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Condition );
        rSerializer.load("PairedGeometry", mpPairedGeometry);
    }

    ///@}

}; // Class PairedCondition

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}// namespace Kratos.

#endif // KRATOS_PAIRED_CONDITION_H_INCLUDED  defined
