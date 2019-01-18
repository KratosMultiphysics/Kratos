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

// External includes

// Project includes
#include "includes/define.h"
#include "includes/node.h"
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

/// Short class definition.
/** Detail class definition.
*/
class GeometricalObject : public IndexedObject
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometricalObject
    KRATOS_CLASS_POINTER_DEFINITION(GeometricalObject);

    typedef IndexedObject BaseType;

    typedef Node < 3 > NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef std::size_t IndexType;

    typedef std::size_t result_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    explicit GeometricalObject(IndexType NewId = 0) : BaseType(NewId),
        mpGeometry()
    {}

    /// Default constructor.
    GeometricalObject(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId),
        mpGeometry(pGeometry)
    {}

    /// Destructor.
    ~GeometricalObject() override {}

    /// Copy constructor.
    GeometricalObject(GeometricalObject const& rOther) : BaseType(rOther.Id()),
        mpGeometry(rOther.mpGeometry)
    {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GeometricalObject& operator=(GeometricalObject const& rOther)
    {
        BaseType::operator=(rOther);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    GeometryType::Pointer pGetGeometry()
    {
        return mpGeometry;
    }

    const GeometryType::Pointer pGetGeometry() const
    {
        return mpGeometry;
    }

    GeometryType& GetGeometry()
    {
        return *mpGeometry;
    }

    GeometryType const& GetGeometry() const
    {
        return *mpGeometry;
    }

    ///@}
    ///@name Inquiry
    ///@{

    inline static bool HasSameType(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return (typeid(rLHS) == typeid(rRHS));
    }

    inline static bool HasSameType(const GeometricalObject * rLHS, const GeometricalObject* rRHS) {
        return GeometricalObject::HasSameType(*rLHS, *rRHS);
    }

    inline static bool HasSameGeometryType(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return (rLHS.GetGeometry().GetGeometryType() == rRHS.GetGeometry().GetGeometryType());
    }

    inline static bool HasSameGeometryType(const GeometricalObject* rLHS, const GeometricalObject* rRHS) {
        return GeometricalObject::HasSameGeometryType(*rLHS, *rRHS);
    }

    inline static bool IsSame(const GeometricalObject& rLHS, const GeometricalObject& rRHS) {
        return GeometricalObject::HasSameType(rLHS, rRHS) && GeometricalObject::HasSameGeometryType(rLHS, rRHS);
    }

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
        buffer << "geometrical object # "
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

    /**
     * pointer to the condition geometry
     */
    GeometryType::Pointer mpGeometry;

    ///@}
    ///@name Private Operators
    ///@{


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
        rSerializer.save("Geometry",mpGeometry);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, IndexedObject );
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


