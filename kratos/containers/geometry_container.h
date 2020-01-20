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
//
//


#if !defined(KRATOS_GEOMETRY_CONTAINER_H_INCLUDED )
#define  KRATOS_GEOMETRY_CONTAINER_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "containers/pointer_hash_map_set.h"
#include "containers/data_value_container.h"


namespace Kratos
{

///@name Kratos Classes
///@{

template<class TGeometryType>
class GeometryContainer
{
    class GetGeometryId : public std::unary_function<const TGeometryType* const, std::size_t>
    {
    public:
        std::size_t const& operator()(const TGeometryType& rGeometry) const
        {
            return rGeometry.Id();
        }
    };

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryContainer
    KRATOS_CLASS_POINTER_DEFINITION(GeometryContainer);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    /* Geometry Hash Map Container.
    *  Hash of Id are keys to corresponding intrusive pointer */
    typedef PointerHashMapSet<
        TGeometryType,
        std::hash<std::size_t>,
        GetGeometryId,
        typename TGeometryType::Pointer
    > GeometriesContainerType;

    /// Geometry Iterator
    typedef typename GeometriesContainerType::iterator GeometryIterator;

    /// Const Geometry Iterator
    typedef typename GeometriesContainerType::const_iterator GeometryConstantIterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor
    GeometryContainer()
        : mGeometries()
    {}

    /// Copy Constructor
    GeometryContainer(GeometryContainer const& rOther)
        : mGeometries(rOther.mGeometries)
    {}

    /// Components Constructor
    GeometryContainer(
        GeometriesContainerType& NewGeometries)
        : mGeometries(NewGeometries)
    {}

    /// Destructor
    ~GeometryContainer() = default;

    ///@}
    ///@name Operations
    ///@{

    GeometryContainer Clone()
    {
        typename GeometriesContainerType::Pointer p_geometries(new GeometriesContainerType(*mGeometries));

        return GeometryContainer(p_geometries);
    }

    void Clear()
    {
        mGeometries.clear();
    }

    ///@}
    ///@name Geometries
    ///@{

    /// Return number of geometries stored inside this geometry container
    SizeType NumberOfGeometries() const
    {
        return mGeometries.size();
    }

    ///@}
    ///@name Add Functions
    ///@{

    /// Adds a geometry to the geometry container.
    GeometryIterator AddGeometry(typename TGeometryType::Pointer pNewGeometry)
    {
        auto i = mGeometries.find(pNewGeometry->Id());
        if(i == mGeometries.end())
            return mGeometries.insert(pNewGeometry);
        else
        {
            KRATOS_ERROR << "Geometry with Id: " << pNewGeometry->Id()
                << " exists already.";
        }
    }

    ///@}
    ///@name Get Functions
    ///@{

    /// Returns the Geometry::Pointer corresponding to its Id
    typename TGeometryType::Pointer pGetGeometry(IndexType GeometryId)
    {
        auto i = mGeometries.find(GeometryId);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryId << ".";
        return (i.base()->second);
    }

    /// Returns the const Geometry::Pointer corresponding to its Id
    const typename TGeometryType::Pointer pGetGeometry(IndexType GeometryId) const
    {
        auto i = mGeometries.find(GeometryId);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryId << ".";
        return (i.base()->second);
    }

    /// Returns the Geometry::Pointer corresponding to its name
    typename TGeometryType::Pointer pGetGeometry(std::string GeometryName)
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);
        auto i = mGeometries.find(hash_index);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryName << ".";
        return (i.base()->second);
    }

    /// Returns the Geometry::Pointer corresponding to its name
    const typename TGeometryType::Pointer pGetGeometry(std::string GeometryName) const
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);
        auto i = mGeometries.find(hash_index);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryName << ".";
        return (i.base()->second);
    }

    /// Returns a reference geometry corresponding to the id
    TGeometryType& GetGeometry(IndexType GeometryId)
    {
        auto i = mGeometries.find(GeometryId);
        KRATOS_ERROR_IF(i == mGeometries.end()) << " geometry index not found: " << GeometryId << ".";
        return *i;
    }

    /// Returns a const reference geometry corresponding to the id
    const TGeometryType& GetGeometry(IndexType GeometryId) const
    {
        auto i = mGeometries.find(GeometryId);
        KRATOS_ERROR_IF(i == mGeometries.end()) << " geometry index not found: " << GeometryId << ".";
        return *i;
    }

    /// Returns a reference geometry corresponding to the name
    TGeometryType& GetGeometry(std::string GeometryName)
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);
        auto i = mGeometries.find(hash_index);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryName << ".";
        return *i;
    }

    /// Returns a const reference geometry corresponding to the name
    const TGeometryType& GetGeometry(std::string GeometryName) const
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);
        auto i = mGeometries.find(hash_index);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryName << ".";
        return *i;
    }

    ///@}
    ///@name Remove Functions
    ///@{

    /// Remove the geometry with given Id from geometry container
    void RemoveGeometry(IndexType GeometryId)
    {
        mGeometries.erase(GeometryId);
    }

    /// Remove the geometry with given name from geometry container
    void RemoveGeometry(std::string GeometryName)
    {
        auto index = TGeometryType::GenerateId(GeometryName);

        mGeometries.erase(index);
    }

    ///@}
    ///@name Search Functions
    ///@{

    bool HasGeometry(IndexType GeometryId) const
    {
        return (mGeometries.find(GeometryId) != mGeometries.end());
    }

    bool HasGeometry(std::string GeometryName) const
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);

        return (mGeometries.find(hash_index) != mGeometries.end());
    }

    ///@}
    ///@name Iterators
    ///@{

    GeometryIterator GeometriesBegin()
    {
        return mGeometries.begin();
    }

    GeometryConstantIterator GeometriesBegin() const
    {
        return mGeometries.begin();
    }

    GeometryIterator GeometriesEnd()
    {
        return mGeometries.end();
    }

    GeometryConstantIterator GeometriesEnd() const
    {
        return mGeometries.end();
    }

    ///@}
    ///@name Container Functions
    ///@{

    GeometriesContainerType& Geometries()
    {
        return mGeometries;
    }

    const GeometriesContainerType& Geometries() const
    {
        return mGeometries;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Return information
    std::string Info() const
    {
        return "GeometryContainer";
    }

    /// Print information about this object
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Number of Geometries: " << mGeometries.size() << std::endl;
    }

    /// Print information about this object
    virtual void PrintInfo(std::ostream& rOStream, std::string const& PrefixString) const
    {
        rOStream << PrefixString << Info();
    }

    /// Print object's data
    virtual void PrintData(std::ostream& rOStream, std::string const& PrefixString ) const
    {
        rOStream << PrefixString << "Number of Geometries: " << mGeometries.size() << std::endl;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// Geometry Container
    GeometriesContainerType mGeometries;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
        rSerializer.save("Geometries", mGeometries);
    }

    void load(Serializer& rSerializer)
    {
        rSerializer.load("Geometries", mGeometries);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    GeometryContainer& operator=(const GeometryContainer& rOther)
    {
        Flags::operator =(rOther);
        mGeometries = rOther.mGeometries;
    }

    ///@}

}; // Class GeometryContainer

///@}
///@name Input and output
///@{

/// input stream function
template<class TGeometryType>
inline std::istream& operator >> (std::istream& rIStream,
                                  GeometryContainer<TGeometryType>& rThis);

/// output stream function
template<class TGeometryType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeometryContainer<TGeometryType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_CONTAINER_H_INCLUDED  defined


