//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//

#if !defined(KRATOS_POINT_H_INCLUDED)
#define KRATOS_POINT_H_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <sstream>
#include <cstddef>

// External includes

// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"
#include "includes/serializer.h"
#include "includes/kratos_components.h"

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

/// Point class.
/** Point class. Stores coordinates of a point and have some basic
    operations defined. 

@see Geometry
@see Node
@see IntegrationPoint
*/
class Point : public array_1d<double, 3>
{
    static constexpr int mDimension = 3;

public:
    ///@name Constants
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Point
    KRATOS_CLASS_POINTER_DEFINITION(Point);

    typedef array_1d<double, mDimension> BaseType;

    typedef BaseType CoordinatesArrayType;

    typedef std::size_t SizeType;

    typedef std::size_t IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Point() : BaseType(mDimension)
    {
        SetAllCoordinates();
    }

    /// 3d constructor.
    Point(double NewX, double NewY = 0, double NewZ = 0) : BaseType(mDimension)
    {
        this->operator()(0) = NewX;
        this->operator()(1) = NewY;
        this->operator()(2) = NewZ;
    }

    /** Copy constructor. Initialize this point with the coordinates
    of given point.*/
    Point(Point const &rOtherPoint)
        : BaseType(rOtherPoint) {}

    /** Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array. */
    Point(CoordinatesArrayType const &rOtherCoordinates)
        : BaseType(rOtherCoordinates) {}

    /** Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array. */
    template <class TVectorType>
    Point(vector_expression<TVectorType> const &rOtherCoordinates)
        : BaseType(rOtherCoordinates) {}

    /** Constructor using coordinates stored in given std::vector. Initialize
    this point with the coordinates in the array. */
    Point(std::vector<double> const &rOtherCoordinates) : BaseType(mDimension)
    {
        SizeType size = rOtherCoordinates.size();
        size = (mDimension < size) ? mDimension : size;
        for (IndexType i = 0; i < size; i++)
            this->operator[](i) = rOtherCoordinates[i];
    }

    /// Destructor.
    virtual ~Point() {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Point &operator=(const Point &rOther)
    {
        CoordinatesArrayType::operator=(rOther);
        return *this;
    }

    bool operator==(const Point &rOther)
    {
        return std::equal(this->begin(), this->end(), rOther.begin());
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    static constexpr IndexType Dimension()
    {
        return 3;
    }

    /** Returns X coordinate */
    double X() const
    {
        return this->operator[](0);
    }

    /** Returns Y coordinate */
    double Y() const
    {
        return this->operator[](1);
    }

    /** Returns Z coordinate */
    double Z() const
    {
        return this->operator[](2);
    }

    double &X()
    {
        return this->operator[](0);
    }

    /** Returns Y coordinate */
    double &Y()
    {
        return this->operator[](1);
    }

    /** Returns Z coordinate */
    double &Z()
    {
        return this->operator[](2);
    }

    /** This is an access method to point's coordinate by indices. For example this
    function return x, y and z coordinate whith 1, 2 and 3 as input
    respectively.
    @note: This method is deprecated: use X(), Y(), Z() or the base class' operator[] method instead.
    */
    KRATOS_DEPRECATED double Coordinate(IndexType CoordinateIndex) const
    {
        KRATOS_DEBUG_ERROR_IF((CoordinateIndex > 3)||(CoordinateIndex == 0)) << "Coordinate index = " << CoordinateIndex << " is out of range [1..3]";
        return this->operator[](CoordinateIndex - 1);
    }

    /** This is an access method to get a reference to point's coordinate by
    indices. For example this function return references to x, y and z coordinate whith 1, 2
    and 3 as input respectively.
    @note: This method is deprecated: use X(), Y(), Z() or the base class' operator[] method instead.
    */
    KRATOS_DEPRECATED double &Coordinate(IndexType CoordinateIndex)
    {
        KRATOS_DEBUG_ERROR_IF((CoordinateIndex > 3)||(CoordinateIndex == 0)) << "Coordinate index = " << CoordinateIndex << " is out of range [1..3]";
        return this->operator[](CoordinateIndex - 1);
    }

    CoordinatesArrayType const &Coordinates() const
    {
        return *this;
    }

    CoordinatesArrayType &Coordinates()
    {
        return *this;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Point";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << this->Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
        rOStream << "(" << this->operator[](0)
                        << this->operator[](1)  
                        << this->operator[](2) 
                 << ")";
    }

    ///@}

  private:
    ///@name Private Operations
    ///@{

    void SetAllCoordinates(double const &Value = double())
    {
        for (IndexType i = 0; i < mDimension; i++)
            this->operator()(i) = Value;
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
        rSerializer.save_base("BaseClass", *static_cast<const array_1d<double, mDimension> *>(this));
    }

    virtual void load(Serializer &rSerializer)
    {
        rSerializer.load_base("BaseClass", *static_cast<array_1d<double, mDimension> *>(this));
    }

    ///@}

}; // Class Point

///@}

template class KRATOS_API(KRATOS_CORE) KratosComponents<Point>;

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                Point &rThis){
                                    return rIStream;
                                }

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const Point &rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_POINT_H_INCLUDED  defined
