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
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/pointer_vector_set.h"
#include "containers/pointer_vector_map.h"
#include "utilities/indexed_object.h"
#include "geometries/geometry.h"
#include "containers/flags.h"
#include "containers/data_value_container.h"


namespace Kratos
{

///@name Kratos Classes
///@{

template<class TNodeType>
class GeometryContainer : public DataValueContainer, public Flags
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryContainer
    KRATOS_CLASS_POINTER_DEFINITION(GeometryContainer);

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    typedef TNodeType NodeType;

    typedef Geometry<NodeType> GeometryType;

    typedef GeometryContainer<TNodeType> GeometryContainerType;

    /// Element container. A vector set of Elements with their Id's as key.
    typedef PointerVectorSet< GeometryType,
                            IndexedObject,
                            std::less<typename IndexedObject::result_type>,
                            std::equal_to<typename IndexedObject::result_type>,
                            typename GeometryType::Pointer,
                            std::vector< typename GeometryType::Pointer >
                            > GeometriesContainerType;

    /// Geometry Iterator
    typedef typename GeometriesContainerType::iterator GeometryIterator;

    /// Const Geometry Iterator
    typedef typename GeometriesContainerType::const_iterator GeometryConstantIterator;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GeometryContainer() : Flags()
        , mpGeometries(new GeometriesContainerType())
    {}

    /// Copy constructor.
    GeometryContainer(GeometryContainer const& rOther)
        : Flags(rOther)
        , mpGeometries(rOther.mpGeometries)
    {}

    /// Components constructor.
    GeometryContainer(typename NodesContainerType::Pointer NewNodes,
         typename GeometriesContainerType::Pointer NewGeometries)
        : Flags()
        , mpGeometries(NewGeometries)
    {}

    /// Destructor.
    ~GeometryContainer() override
    {}

    ///@}
    ///@name Operations
    ///@{

    GeometryContainer Clone()
    {
        typename GeometryContainerType::Pointer p_geometries(new GeometryContainerType(*mpGeometries));

        return GeometryContainer(p_geometries);
    }

    void Clear()
    {
        Flags::Clear();
        DataValueContainer::Clear();
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
    ///@name Add / Get / Remove Functions
    ///@{

    /// Adds a geometry to the geometry container.
    void AddGeometry(typename GeometryType::Pointer pNewGeometry)
    {
        mpGeometries->push_back(pNewGeometry);
    }

    /// Returns the Geometry::Pointer corresponding to it's identifier
    typename GeometryType::Pointer pGetGeometry(IndexType GeometryId)
    {
        // HOW?
        //auto i = mpGeometries->find(GeometryId);
        KRATOS_ERROR_IF(i == mpGeometries->end()) << " geometry index not found: " << GeometryId << ".";
        return *i.base();
    }

    /// Returns a reference geometry corresponding to it's identifier
    GeometryType& GetGeometry(IndexType GeometryId)
    {
        // HOW?
        //auto i = mpGeometries->find(GeometryId);
        KRATOS_ERROR_IF(i == mpGeometries->end()) << " geometry index not found: " << GeometryId << ".";
        return *i;
    }

    /// Returns a const reference geometry corresponding to it's identifier
    const GeometryType& GetGeometry(IndexType GeometryId) const
    {
        // HOW?
        //auto i = mpGeometries->find(GeometryId);
        KRATOS_ERROR_IF(i == mpGeometries->end()) << " geometry index not found: " << GeometryId << ".";
        return *i;
    }

    /// Remove the geometry with given Id from geometry container
    void RemoveGeometry(IndexType GeometryId)
    {
        // HOW?
        //mpGeometries->erase(GeometryId);
    }

    /// Remove given element from geometry container
    void RemoveGeometry(GeometryType& ThisGeometry)
    {
        // HOW?
        //mpGeometries->erase(ThisGeometry.Id());
    }

    /// Remove given element from geometry container
    void RemoveGeometry(typename GeometryType::Pointer pThisGeometry)
    {
        // HOW?
        //mpGeometries->erase(ThisGeometry->Id());
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

    typename GeometriesContainerType::ContainerType& GeometriesArray()
    {
        return mpGeometries->GetContainer();
    }

    ///@}
    ///@name Container Functions
    ///@{

    GeometriesContainerType& Geometries(GeometryData::KratosGeometryType GeometryType)
    {
        return *mpGeometries->where(...);
    }

    ///@}
    ///@name Container Functions
    ///@{

    GeometriesContainerType& Geometries(GeometryData::KratosGeometryType GeometryType)
    {
        return *mpGeometries->where(...);
    }

    ///@}
    ///@name Search Functions
    ///@{

    bool HasGeometry(IndexType GeometryId) const
    {
        // HOW?
        //return (mpGeometries->find(GeometryId) != mpGeometries->end());
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GeometryContainer";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << "Number of Geometries: " << mpGeometries->size() << std::endl;
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream, std::string const& PrefixString) const
    {
        rOStream << PrefixString << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream, std::string const& PrefixString ) const
    {
        rOStream << PrefixString << "Number of Geometries: " << mpGeometries->size() << std::endl;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    /// Geometry Container
    typename GeometryContainerType::Pointer mpGeometries;

    /// Model Tolerance
    double mModelTolerance;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DataValueContainer );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags );
        rSerializer.save("pGeometries", mpGeometries);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DataValueContainer );
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags );
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
template<class TNodeType>
inline std::istream& operator >> (std::istream& rIStream,
                                  GeometryContainer<TNodeType>& rThis);

/// output stream function
template<class TNodeType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeometryContainer<TNodeType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_CONTAINER_H_INCLUDED  defined


