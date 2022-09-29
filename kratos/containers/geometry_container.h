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

    typedef typename TGeometryType::Pointer GeometryPointerType;


    /// Geometry Hash Map Container.
    // Stores with hash of Ids to corresponding geometries.
    typedef PointerHashMapSet<
        TGeometryType,
        std::hash<std::size_t>,
        GetGeometryId,
        GeometryPointerType
        > GeometriesMapType;

    /// Geometry Iterator
    typedef typename GeometriesMapType::iterator GeometryIterator;

    /// Const Geometry Iterator
    typedef typename GeometriesMapType::const_iterator GeometryConstantIterator;

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
        GeometriesMapType& NewGeometries)
        : mGeometries(NewGeometries)
    {}

    /// Destructor
    ~GeometryContainer() = default;

    ///@}
    ///@name Operations
    ///@{

    GeometryContainer Clone()
    {
        typename GeometriesMapType::Pointer p_geometries(new GeometriesMapType(*mGeometries));

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
    GeometryIterator AddGeometry(GeometryPointerType pNewGeometry)
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
    GeometryPointerType pGetGeometry(IndexType GeometryId)
    {
        auto i = mGeometries.find(GeometryId);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryId << ".";
        return (i.base()->second);
    }

    /// Returns the const Geometry::Pointer corresponding to its Id
    const GeometryPointerType pGetGeometry(IndexType GeometryId) const
    {
        auto i = mGeometries.find(GeometryId);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryId << ".";
        return (i.base()->second);
    }

    /// Returns the Geometry::Pointer corresponding to its name
    GeometryPointerType pGetGeometry(std::string GeometryName)
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);
        auto i = mGeometries.find(hash_index);
        KRATOS_ERROR_IF(i == mGeometries.end())
            << " geometry index not found: " << GeometryName << ".";
        return (i.base()->second);
    }

    /// Returns the Geometry::Pointer corresponding to its name
    const GeometryPointerType pGetGeometry(std::string GeometryName) const
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
        return *pGetGeometry(GeometryId);
    }

    /// Returns a const reference geometry corresponding to the id
    const TGeometryType& GetGeometry(IndexType GeometryId) const
    {
        return *pGetGeometry(GeometryId);
    }

    /// Returns a reference geometry corresponding to the name
    TGeometryType& GetGeometry(std::string GeometryName)
    {
        return *pGetGeometry(GeometryName);
    }

    /// Returns a const reference geometry corresponding to the name
    const TGeometryType& GetGeometry(std::string GeometryName) const
    {
        return *pGetGeometry(GeometryName);
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

    GeometriesMapType& Geometries()
    {
        return mGeometries;
    }

    const GeometriesMapType& Geometries() const
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
    GeometriesMapType mGeometries;

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
        mGeometries = rOther.mGeometries;
        return *this;
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


