//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//
//

#if !defined(KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED )
#define  KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED

// System includes
#include <atomic>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "containers/flags.h"
#include "geometries/geometry.h"

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
 * @class GeometricalObject
 * @ingroup KratosCore
 * @brief This defines the geometrical object, base definition of the element and condition entities
 * @details Derives from IndexedObject, so it has an ID, and from Flags
 * @author Pooyan Dadvand
*/
class GeometricalObject : public IndexedObject, public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObject
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeometricalObject);

    /// Definition of the node type
    typedef Node <3> NodeType;

    /// The geometry type definition
    typedef Geometry<NodeType> GeometryType;

    /// Defines the index type
    typedef std::size_t IndexType;

    /// Defines the result type
    typedef std::size_t result_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit GeometricalObject(IndexType NewId = 0)
        : IndexedObject(NewId),
          Flags(),
          mpGeometry(),
          mReferenceCounter(0)
    {}

    /// Default constructor.
    GeometricalObject(IndexType NewId, GeometryType::Pointer pGeometry)
        : IndexedObject(NewId),
          Flags(),
          mpGeometry(pGeometry),
          mReferenceCounter(0)
    {}

    /// Destructor.
    ~GeometricalObject() override {}

    /// Copy constructor.
    GeometricalObject(GeometricalObject const& rOther)
        : IndexedObject(rOther.Id()),
          Flags(rOther),
          mpGeometry(rOther.mpGeometry),
          mReferenceCounter(0)
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GeometricalObject& operator=(GeometricalObject const& rOther)
    {
        IndexedObject::operator=(rOther);
        Flags::operator =(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Sets the pointer to the geometry
     * @param pGeometry The pointer of the geometry
     */
    virtual void SetGeometry(GeometryType::Pointer pGeometry)
    {
        mpGeometry = pGeometry;
    }

    /**
     * @brief Returns the pointer to the geometry
     * @return The pointer of the geometry
     */
    GeometryType::Pointer pGetGeometry()
    {
        return mpGeometry;
    }

    /**
     * @brief Returns the pointer to the geometry (const version)
     * @return The pointer of the geometry
     */
    const GeometryType::Pointer pGetGeometry() const
    {
        return mpGeometry;
    }

    /**
     * @brief Returns the reference of the geometry
     * @return The reference of the geometry
     */
    GeometryType& GetGeometry()
    {
        return *mpGeometry;
    }

    /**
     * @brief Returns the reference of the geometry (const version)
     * @return The reference of the geometry
     */
    GeometryType const& GetGeometry() const
    {
        return *mpGeometry;
    }

    /**
     * @brief Returns the flags of the object
     * @return The  flags of the object
     */
    Flags& GetFlags()
    {
        return *this;
    }

    /**
     * @brief Returns the flags of the object (const version)
     * @return The  flags of the object
     */
    Flags const& GetFlags() const
    {
        return *this;
    }

    /**
     * @brief Sets the flags of the object
     * @param rThisFlags The flags to be set
     */
    void SetFlags(Flags const& rThisFlags)
    {
        Flags::operator=(rThisFlags);
    }

    ///@}
    ///@name Data
    ///@{

    /**
     * Access Data:
     */
    DataValueContainer& Data()
    {
        return pGetGeometry()->GetData();
    }

    DataValueContainer const& GetData() const
    {
        return GetGeometry().GetData();
    }

    void SetData(DataValueContainer const& rThisData)
    {
        return GetGeometry().SetData(rThisData);
    }

    /**
     * Check if the Data exists with Has(..) methods:
     */
    template<class TDataType> bool Has(const Variable<TDataType>& rThisVariable) const
    {
        return GetData().Has(rThisVariable);
    }

    template<class TAdaptorType> bool Has(
        const VariableComponent<TAdaptorType>& rThisVariable) const
    {
        return GetData().Has(rThisVariable);
    }

    /**
     * Set Data with SetValue and the Variable to set:
     */
    template<class TVariableType> void SetValue(
        const TVariableType& rThisVariable,
        typename TVariableType::Type const& rValue)
    {
        Data().SetValue(rThisVariable, rValue);
    }

    /**
     * Get Data with GetValue and the Variable to get:
     */
    template<class TVariableType> typename TVariableType::Type& GetValue(
        const TVariableType& rThisVariable)
    {
        return Data().GetValue(rThisVariable);
    }

    template<class TVariableType> typename TVariableType::Type const& GetValue(
        const TVariableType& rThisVariable) const
    {
        return GetData().GetValue(rThisVariable);
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief Checks if two GeometricalObject have the same type
     * @return True if the objects are the same type, false otherwise
     */
    inline static bool HasSameType(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return (typeid(rLHS) == typeid(rRHS));
    }

    /**
     * @brief Checks if two GeometricalObject have the same type (pointer version)
     * @return True if the objects are the same type, false otherwise
     */
    inline static bool HasSameType(const GeometricalObject * rLHS, const GeometricalObject* rRHS) {
        return GeometricalObject::HasSameType(*rLHS, *rRHS);
    }

    /**
     * @brief Checks if two GeometricalObject have the same geometry type
     * @return True if the geometries are the same type, false otherwise
     */
    inline static bool HasSameGeometryType(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return (rLHS.GetGeometry().GetGeometryType() == rRHS.GetGeometry().GetGeometryType());
    }

    /**
     * @brief Checks if two GeometricalObject have the same geometry type (pointer version)
     * @return True if the geometries are the same type, false otherwise
     */
    inline static bool HasSameGeometryType(const GeometricalObject* rLHS, const GeometricalObject* rRHS) {
        return GeometricalObject::HasSameGeometryType(*rLHS, *rRHS);
    }

    /**
     * @brief Checks if two GeometricalObject are the same
     * @return True if the object is the same, false otherwise
     */
    inline static bool IsSame(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return GeometricalObject::HasSameType(rLHS, rRHS) && GeometricalObject::HasSameGeometryType(rLHS, rRHS);
    }

    /**
     * @brief Checks if two GeometricalObject are the same (pointer version)
     * @return True if the object is the same, false otherwise
     */
    inline static bool IsSame(const GeometricalObject* rLHS, const GeometricalObject* rRHS) {
        return GeometricalObject::HasSameType(*rLHS, *rRHS) && GeometricalObject::HasSameGeometryType(*rLHS, *rRHS);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Geometrical object # "
               << Id();
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    //*********************************************
    //public API of intrusive_ptr
    unsigned int use_count() const noexcept
    {
        return mReferenceCounter;
    }
    //*********************************************

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

    GeometryType::Pointer mpGeometry; /// Pointer to the entity geometry

    ///@}
    ///@name Private Operators
    ///@{

    //*********************************************
    //this block is needed for refcounting
    mutable std::atomic<int> mReferenceCounter;

    friend void intrusive_ptr_add_ref(const GeometricalObject* x)
    {
        x->mReferenceCounter.fetch_add(1, std::memory_order_relaxed);
    }

    friend void intrusive_ptr_release(const GeometricalObject* x)
    {
        if (x->mReferenceCounter.fetch_sub(1, std::memory_order_release) == 1) {
        std::atomic_thread_fence(std::memory_order_acquire);
        delete x;
        }
    }
    //*********************************************


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, IndexedObject );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("Geometry",mpGeometry);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
        rSerializer.load("Geometry",mpGeometry);
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

}; // Class GeometricalObject

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GeometricalObject& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeometricalObject& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GEOMETRICAL_OBJECT_H_INCLUDED  defined


