//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Tobias Teschemacher
//
//


#if !defined(KRATOS_GEOMETRY_MODEL_H_INCLUDED )
#define  KRATOS_GEOMETRY_MODEL_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>
#include <unordered_set>


// External includes


// Project includes
#include "includes/define.h"
#include "containers/pointer_vector_set.h"
#include "utilities/indexed_object.h"
#include "geometries/geometry.h"
#include "containers/flags.h"
#include "containers/data_value_container.h"


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

/// Mesh is the second level of abstraction in the data structure which hold Nodes, Elements and Conditions and their Properties.
/**
 * Mesh is the second level of abstraction in the data structure which hold Nodes, Elements and Conditions and their Properties.
 * In other words, Mesh is a complete pack of all type of entities without any additional data associated with them.
 * So a set of Elements and Conditions with their Nodes and Properties can be grouped together as a Mesh and send to
 * procedures like mesh refinement, material optimization, mesh movement or any other procedure which works on entities
 * without needing additional data for their processes.
*/
template<class TGeometryType>
class GeometryModel : public DataValueContainer, public Flags
{
public:


    ///@name Type Definitions
    ///@{

    /// Pointer definition of Mesh
    KRATOS_CLASS_POINTER_DEFINITION(GeometryModel);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef TGeometryType GeometryType;

    typedef GeometryModel<TGeometryType> GeometryModelType;

    /// Geometries container. A vector set of Geometries with their Pointers as key.
    //typedef std::unordered_set<typename GeometryType::Pointer> GeometriesContainerType;
    //typedef std::unordered_set<typename GeometryModelType::Pointer> GeometriesContainerType;
    typedef PointerVectorSet<GeometryType> GeometriesContainerType;
        //typename GeometryType::Pointer,
        //std::less<typename GeometryType::Pointer>,
        //std::equal_to<typename GeometryType::Pointer>,
        //typename GeometryType::Pointer,
        //std::vector< typename GeometryType::Pointer >
        //> 

    /** Iterator over the Geometries. This iterator is an indirect
    iterator over Geometries::Pointer which turn back a reference to
    Geometry by * operator and not a pointer for more convenient
    usage. */
    typedef typename GeometriesContainerType::iterator GeometryIterator;

    /** Const iterator over the Geometries. This iterator is an indirect
    iterator over Geometries::Pointer which turn back a reference to
    Geometry by * operator and not a pointer for more convenient
    usage. */
    typedef typename GeometriesContainerType::const_iterator GeometryConstantIterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeometryModel()
        : Flags()
        , mGeometries(GeometriesContainerType())
    {}

    /// Copy constructor.
    GeometryModel(GeometryModel const& rOther)
        : Flags(rOther)
        , mGeometries(rOther.mGeometries)
    {}

    /// Components constructor.
    GeometryModel(GeometriesContainerType NewGeometries)
        : Flags()
        , mGeometries(NewGeometries)
    {}


    /// Destructor.
    ~GeometryModel() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    GeometryModel Clone()
    {
        typename GeometriesContainerType::Pointer p_geometries(new GeometriesContainerType(*mGeometries));

        return GeometryModel(p_geometries);
    }

    void Clear()
    {
        Flags::Clear();
        DataValueContainer::Clear();
        mGeometries.clear();
    }

    ///@}
    ///@name Informations
    ///@{

    ///@}
    ///@name Geometries
    ///@{

    SizeType NumberOfGeometries() const
    {
        return mGeometries.size();
    }

    /** Inserts a condition in the mesh.
    */
    void AddGeometry(typename GeometryType::Pointer pNewGeometry)
    {
        //mGeometries->insert(pNewGeometry);
    }

    ///** Returns the Geometry::Pointer  corresponding to it's identifier */
    //typename GeometryType::Pointer pGetGeometry(IndexType GeometryId)
    //{
    //    auto i = mGeometries->find(GeometryId);
    //    KRATOS_ERROR_IF(i == mGeometries->end()) << " geometry index not found: " << GeometryId << ".";
    //    return *i.base();
    //}

    ///** Returns a reference geometry corresponding to it's identifier */
    //GeometryType& GetGeometry(IndexType GeometryId)
    //{
    //    auto i = mGeometries->find(GeometryId);
    //    KRATOS_ERROR_IF(i == mGeometries->end()) << " geometry index not found: " << GeometryId << ".";
    //    return *i;
    //}

    //const GeometryType& GetGeometry(IndexType GeometryId) const
    //{
    //    auto i = mGeometries->find(GeometryId);
    //    KRATOS_ERROR_IF(i == mGeometries->end()) << " geometry index not found: " << GeometryId << ".";
    //    return *i;
    //}

    ///** Remove the condition with given Id from mesh.
    //*/
    //void RemoveGeometry(IndexType GeometryId)
    //{
    //    mGeometries->erase(GeometryId);
    //}

    ///** Remove given condition from mesh.
    //*/
    //void RemoveGeometry(GeometryType& ThisGeometry)
    //{
    //    mGeometries->erase(ThisGeometry.Id());
    //}

    /** Remove given condition from mesh.
    */
    void RemoveGeometry(typename GeometryType::Pointer pThisGeometry)
    {
        //mGeometries->erase(pThisGeometry);
    }

    //GeometryIterator GeometriesBegin()
    //{
    //    return mGeometries.begin();
    //}

    //GeometryConstantIterator GeometriesBegin() const
    //{
    //    return mGeometries.begin();
    //}

    //GeometryIterator GeometriesEnd()
    //{
    //    return mGeometries.end();
    //}

    //GeometryConstantIterator GeometriesEnd() const
    //{
    //    return mGeometries.end();
    //}

    GeometriesContainerType& Geometries()
    {
        return mGeometries;
    }

    //const GeometriesContainerType& Geometries() const
    //{
    //    return *mGeometries;
    //}

    void SetGeometries(GeometriesContainerType OtherGeometries)
    {
        mGeometries = OtherGeometries;
    }

    bool HasGeometry() const
    {
        return true;// (mGeometries->find(GeometryId) != mGeometries->end());
    }

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
        return "GeometryModel";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "    Number of Geometries  : " << mGeometries.size() << std::endl;
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream, std::string const& PrefixString) const
    {
        rOStream << PrefixString << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& PrefixString ) const
    {
        rOStream << PrefixString << "    Number of Geometries  : " << mGeometries.size() << std::endl;
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

    GeometriesContainerType mGeometries;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DataValueContainer );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        //rSerializer.save("Geometries",mGeometries);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DataValueContainer );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
        //rSerializer.load("Geometries",mGeometries);
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

    /// Assignment operator.
    GeometryModel& operator=(const GeometryModel& rOther)
    {
        Flags::operator =(rOther);
        mGeometries = rOther.mGeometries;
    }

    ///@}

}; // Class GeometryModel

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TGeometryType>
inline std::istream& operator >> (
    std::istream& rIStream,
    GeometryModel<TGeometryType>& rThis);

/// output stream function
template<class TGeometryType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const GeometryModel<TGeometryType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_MODEL_H_INCLUDED  defined


