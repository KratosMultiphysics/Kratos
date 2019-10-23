//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//                   Philip Bucher
//                   Pooyan Dadvand
//
//

#if !defined(KRATOS_GEOMETRY_DIMENSION_H_INCLUDED )
#define  KRATOS_GEOMETRY_DIMENSION_H_INCLUDED

// System includes

// External includes

// Project includes

namespace Kratos
{

///@name Kratos Classes
///@{

/** GeometryDimension
    Has the values for geometrical dimension,
    working space dimension and local space dimension.

    @see Geometry
    @see Element
*/

class GeometryDimension
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeometryDimension
    KRATOS_CLASS_POINTER_DEFINITION( GeometryDimension );

    typedef std::size_t IndexType;
    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** This constructor gives a complete set of dimensional
    arguments neede for geometries.

    @param ThisDimension Dimension of this geometry.

    @param ThisWorkingSpaceDimension Working space dimension. for
    example a triangle 3d is a 2 dimensional shape but can be used
    in 3 dimensional space.

    @param ThisLocalSpaceDimension Local space dimension.
    for example a triangle is a 2 dimensional shape but
    can have 3 dimensional area coordinates l1, l2, l3.
    */
    GeometryDimension( SizeType ThisDimension,
                  SizeType ThisWorkingSpaceDimension,
                  SizeType ThisLocalSpaceDimension )
        : mDimension( ThisDimension )
        , mWorkingSpaceDimension( ThisWorkingSpaceDimension )
        , mLocalSpaceDimension( ThisLocalSpaceDimension )
    {
    }

    /** Copy constructor.
    Construct this geometry data as a copy of given geometry data.
    */
    GeometryDimension( const GeometryDimension& rOther )
        : mDimension( rOther.mDimension )
        , mWorkingSpaceDimension( rOther.mWorkingSpaceDimension )
        , mLocalSpaceDimension( rOther.mLocalSpaceDimension )
    {
    }

    /// Destructor. Do nothing!!!
    virtual ~GeometryDimension() {}

    ///@}
    ///@name Operators
    ///@{

    GeometryDimension& operator=( const GeometryDimension& rOther )
    {
        mDimension = rOther.mDimension;
        mWorkingSpaceDimension = rOther.mWorkingSpaceDimension;
        mLocalSpaceDimension = rOther.mLocalSpaceDimension;

        return *this;
    }

    ///@}
    ///@name Informations
    ///@{

    /** Dimension of the geometry for example a triangle2d is a 2
    dimensional shape. A line is 1 dimensional shape and
    a point 0 dimensional.

    @return SizeType, dimension of this geometry.
    @see WorkingSpaceDimension()
    @see LocalSpaceDimension()
    */
    inline SizeType Dimension() const
    {
        return mDimension;
    }

    /** Working space dimension. for example a triangle as part
    of a surface is a 2 dimensional shape but can be used in
    3 dimensional space.

    @return SizeType, working space dimension of this geometry.
    @see Dimension()
    @see LocalSpaceDimension()
    */
    inline SizeType WorkingSpaceDimension() const
    {
        return mWorkingSpaceDimension;
    }

    /** Local space dimension. for example a triangle is a 2
    dimensional shape but can have 3 dimensional area
    coordinates l1, l2, l3.

    @return SizeType, local space dimension of this geometry.
    @see Dimension()
    @see WorkingSpaceDimension()
    */
    inline  SizeType LocalSpaceDimension() const
    {
        return mLocalSpaceDimension;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "geometry dimension";
    }

    /// Print information about this object.
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "geometry dimension";
    }

    /// Print object's data.
    virtual void PrintData( std::ostream& rOStream ) const
    {
        rOStream << "    Dimension               : " << mDimension << std::endl;
        rOStream << "    Working space dimension : " << mWorkingSpaceDimension << std::endl;
        rOStream << "    Local space dimension   : " << mLocalSpaceDimension;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    SizeType mDimension;

    SizeType mWorkingSpaceDimension;

    SizeType mLocalSpaceDimension;

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        rSerializer.save( "Dimension", mDimension );
        rSerializer.save( "WorkingSpaceDimension", mWorkingSpaceDimension );
        rSerializer.save( "LocalSpaceDimension", mLocalSpaceDimension );
    }

    virtual void load( Serializer& rSerializer )
    {
        rSerializer.load( "Dimension", mDimension );
        rSerializer.load( "WorkingSpaceDimension", mWorkingSpaceDimension );
        rSerializer.load( "LocalSpaceDimension", mLocalSpaceDimension );
    }

    // Private default constructor for serialization
    GeometryDimension()
    {
    }

    ///@}

}; // Class GeometryDimension

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> ( std::istream& rIStream,
                                   GeometryDimension& rThis );

/// output stream function
inline std::ostream& operator << ( std::ostream& rOStream,
                                   const GeometryDimension& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_GEOMETRY_DIMENSION_H_INCLUDED  defined


