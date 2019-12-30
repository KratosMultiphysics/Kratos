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
        : mpGeometries(new GeometriesContainerType())
    {}

    /// Copy Constructor
    GeometryContainer(GeometryContainer const& rOther)
        : mpGeometries(rOther.mpGeometries)
    {}

    /// Components Constructor
    GeometryContainer(
        typename GeometriesContainerType::Pointer NewGeometries)
        : mpGeometries(NewGeometries)
    {}

    /// Destructor
    ~GeometryContainer() override
    {}

    ///@}
    ///@name Operations
    ///@{

    GeometryContainer Clone()
    {
        typename GeometriesContainerType::Pointer p_geometries(new GeometriesContainerType(*mpGeometries));

        return GeometryContainer(p_geometries);
    }

    void Clear()
    {
        mpGeometries->clear();
    }

    ///@}
    ///@name Geometries
    ///@{

    /// Return number of geometries stored inside this geometry container
    SizeType NumberOfGeometries() const
    {
        return mpGeometries->size();
    }

    ///@}
    ///@name Add Functions
    ///@{

    /// Adds a geometry to the geometry container.
    GeometryIterator AddGeometry(typename TGeometryType::Pointer pNewGeometry)
    {
        auto i = mpGeometries->find(pNewGeometry->Id());
        if(i == mpGeometries->end())
            return mpGeometries->insert(pNewGeometry);
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
        auto i = mpGeometries->find(GeometryId);
        KRATOS_ERROR_IF(i == mpGeometries->end())
            << " geometry index not found: " << GeometryId << ".";
        return (i.base()->second);
    }

    /// Returns the Geometry::Pointer corresponding to its name
    typename TGeometryType::Pointer pGetGeometry(std::string GeometryName)
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);
        auto i = mpGeometries->find(hash_index);
        KRATOS_ERROR_IF(i == mpGeometries->end())
            << " geometry index not found: " << GeometryName << ".";
        return (i.base()->second);
    }

    /// Returns a reference geometry corresponding to the id
    TGeometryType& GetGeometry(IndexType GeometryId)
    {
        auto i = mpGeometries->find(GeometryId);
        KRATOS_ERROR_IF(i == mpGeometries->end()) << " geometry index not found: " << GeometryId << ".";
        return *i;
    }

    /// Returns a const reference geometry corresponding to the name
    const TGeometryType& GetGeometry(std::string GeometryName) const
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);
        auto i = mpGeometries->find(hash_index);
        KRATOS_ERROR_IF(i == mpGeometries->end())
            << " geometry index not found: " << GeometryName << ".";
        return *i;
    }

    ///@}
    ///@name Remove Functions
    ///@{

    /// Remove the geometry with given Id from geometry container
    void RemoveGeometry(IndexType GeometryId)
    {
        mpGeometries->erase(GeometryId);
    }

    /// Remove the geometry with given name from geometry container
    void RemoveGeometry(std::string GeometryName)
    {
        auto index = TGeometryType::GenerateId(GeometryName);

        mpGeometries->erase(index);
    }

    /// Remove given geometry from geometry container
    void RemoveGeometry(TGeometryType& ThisGeometry)
    {
        mpGeometries->erase(ThisGeometry.Id());
    }

    /// Remove given geometry from geometry container
    void RemoveGeometry(typename TGeometryType::Pointer pThisGeometry)
    {
        mpGeometries->erase(pThisGeometry->Id());
    }

    ///@}
    ///@name Search Functions
    ///@{

    bool HasGeometry(IndexType GeometryId) const
    {
        return (mpGeometries->find(GeometryId) != mpGeometries->end());
    }

    bool HasGeometry(std::string GeometryName) const
    {
        auto hash_index = TGeometryType::GenerateId(GeometryName);

        return (mpGeometries->find(hash_index) != mpGeometries->end());
    }

    ///@}
    ///@name Iterators
    ///@{

    GeometryIterator GeometriesBegin()
    {
        return mpGeometries->begin();
    }

    GeometryConstantIterator GeometriesBegin() const
    {
        return mpGeometries->begin();
    }

    GeometryIterator GeometriesEnd()
    {
        return mpGeometries->end();
    }

    GeometryConstantIterator GeometriesEnd() const
    {
        return mpGeometries->end();
    }

    ///@}
    ///@name Container Functions
    ///@{

    GeometriesContainerType& Geometries()
    {
        return *mpGeometries;
    }

    const GeometriesContainerType& Geometries() const
    {
        return *mpGeometries;
    }

    typename GeometriesContainerType::Pointer pGeometries()
    {
        return mpGeometries;
    }

    void SetGeometries(typename GeometriesContainerType::Pointer pOtherGeometries)
    {
        mpGeometries = pOtherGeometries;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Return information
    std::string Info() const override
    {
        return "GeometryContainer";
    }

    /// Print information about this object
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Number of Geometries: " << mpGeometries->size() << std::endl;
    }

    /// Print information about this object
    virtual void PrintInfo(std::ostream& rOStream, std::string const& PrefixString) const
    {
        rOStream << PrefixString << Info();
    }

    /// Print object's data
    virtual void PrintData(std::ostream& rOStream, std::string const& PrefixString ) const
    {
        rOStream << PrefixString << "Number of Geometries: " << mpGeometries->size() << std::endl;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// Geometry Container
    typename GeometriesContainerType::Pointer mpGeometries;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        rSerializer.save("pGeometries", mpGeometries);
    }

    void load(Serializer& rSerializer) override
    {
        rSerializer.load("pGeometries", mpGeometries);
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    GeometryContainer& operator=(const GeometryContainer& rOther)
    {
        Flags::operator =(rOther);
        mpGeometries = rOther.mpGeometries;
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


