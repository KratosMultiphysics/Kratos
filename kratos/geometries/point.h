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
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Point
    KRATOS_CLASS_POINTER_DEFINITION(Point);

    typedef array_1d<double, Dimension> BaseType;

    typedef BaseType CoordinatesArrayType;

    typedef typename std::size_t SizeType;

    typedef typename std::size_t IndexType;

    ///@}
    ///@name Constants
    ///@{

    constexpr int Dimension = 3;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Point() : BaseType(Dimension)
    {
        SetAllCoordinates();
    }

    /// 3d constructor.
    Point(double NewX, double NewY = 0, double NewZ = 0) : BaseType(Dimension)
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
    Point(std::vector<double> const &rOtherCoordinates) : BaseType(Dimension)
    {
        SizeType size = rOtherCoordinates.size();
        size = (Dimension < size) ? Dimension : size;
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

    /// Assignment operator.
    template <SizeType TOtherDimension>
    Point &operator=(const Point<TOtherDimension> &rOther)
    {
        KRATOS_TRY_LEVEL_4
        IndexType size = (Dimension < TOtherDimension) ? Dimension : TOtherDimension;
        IndexType i;

        for (i = 0; i < size; i++)
            this->operator[](i) = rOther[i];

        for (i = size; i < Dimension; i++)
            this->operator[](i) = double();

        return *this;

        KRATOS_CATCH_LEVEL_4(*this)
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    static IndexType Dimension()
    {
        return Dimension;
    }

    /** Returns X coordinate */
    double X() const
    {
        KRATOS_TRY_LEVEL_4
        return this->operator[](0);
        KRATOS_CATCH_LEVEL_4(*this)
    }

    /** Returns Y coordinate */
    double Y() const
    {
        KRATOS_TRY_LEVEL_4
        return this->operator[](1);
        KRATOS_CATCH_LEVEL_4(*this)
    }

    /** Returns Z coordinate */
    double Z() const
    {
        KRATOS_TRY_LEVEL_4
        return this->operator[](2);
        KRATOS_CATCH_LEVEL_4(*this)
    }

    double &X()
    {
        KRATOS_TRY_LEVEL_4
        return this->operator[](0);
        KRATOS_CATCH_LEVEL_4(*this)
    }

    /** Returns Y coordinate */
    double &Y()
    {
        KRATOS_TRY_LEVEL_4
        return this->operator[](1);
        KRATOS_CATCH_LEVEL_4(*this)
    }

    /** Returns Z coordinate */
    double &Z()
    {
        KRATOS_TRY_LEVEL_4
        return this->operator[](2);
        KRATOS_CATCH_LEVEL_4(*this)
    }

    /** This is an access method to point's coordinate by indices. For example this
    function return x, y and z coordinate whith 1, 2 and 3 as input
    respectively.
    */
    double Coordinate(IndexType CoordinateIndex) const
    {
        KRATOS_TRY_LEVEL_4
        return this->operator[](CoordinateIndex - 1);
        KRATOS_CATCH_LEVEL_4(*this)
    }

    /** This is an access method to get a reference to point's coordinate by
    indices. For example this function return references to x, y and z coordinate whith 1, 2
    and 3 as input respectively.
    */
    double &Coordinate(IndexType CoordinateIndex)
    {
        KRATOS_TRY_LEVEL_4
        return this->operator[](CoordinateIndex - 1);
        KRATOS_CATCH_LEVEL_4(*this)
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
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << Dimension << " dimensional point";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << Dimension << " dimensional point";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
        if (!Dimension)
            return;

        rOStream << "(" << this->operator[](0);

        for (IndexType i = 1; i < Dimension; i++)
            rOStream << " , " << this->operator[](i);
        rOStream << ")";
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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void SetAllCoordinates(double const &Value = double())
    {
        KRATOS_TRY_LEVEL_4
        for (IndexType i = 0; i < Dimension; i++)
            this->operator()(i) = Value;
        KRATOS_CATCH_LEVEL_4(*this)
    }

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer &rSerializer) const
    {
        rSerializer.save_base("BaseClass", *static_cast<const array_1d<double, Dimension> *>(this));
        //rSerializer.save_base("BaseData",*dynamic_cast<const array_1d<double, Dimension>*>(this));
    }

    virtual void load(Serializer &rSerializer)
    {
        rSerializer.load_base("BaseClass", *static_cast<array_1d<double, Dimension> *>(this));
        //	  rSerializer.load_base("BaseData",*dynamic_cast<array_1d<double, Dimension>*>(this));
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

}; // Class Point

///@}

template class KRATOS_API(KRATOS_CORE) KratosComponents<Point<3, double>>;

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
template <std::size_t Dimension, class double>
inline std::istream &operator>>(std::istream &rIStream,
                                Point<Dimension, double> &rThis);

/// output stream function
template <std::size_t Dimension, class double>
inline std::ostream &operator<<(std::ostream &rOStream,
                                const Point<Dimension, double> &rThis)
{
    rThis.PrintInfo(rOStream);
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_POINT_H_INCLUDED  defined
