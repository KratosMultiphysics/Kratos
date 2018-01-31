//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Janosch Stascheit
//                   Felix Nagel
//  contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//

// Project includes
#include "geometries/point.h"

namespace Kratos {

/// Default constructor.
Point::Point() : BaseType(mDimension)
{
    SetAllCoordinates();
}

/// 3d constructor.
Point::Point(double NewX, double NewY, double NewZ) : BaseType(mDimension)
{
    this->operator()(0) = NewX;
    this->operator()(1) = NewY;
    this->operator()(2) = NewZ;
}

/** Copy constructor. Initialize this point with the coordinates
of given point.*/
Point::Point(Point const &rOtherPoint) : BaseType(rOtherPoint) {}

/** Constructor using coordinates stored in given array. Initialize
this point with the coordinates in the array. */
Point::Point(CoordinatesArrayType const &rOtherCoordinates) : BaseType(rOtherCoordinates) {}

/** Constructor using coordinates stored in given std::vector. Initialize
this point with the coordinates in the array. */
Point::Point(std::vector<double> const &rOtherCoordinates) : BaseType(mDimension)
{
    SizeType size = rOtherCoordinates.size();
    size = (mDimension < size) ? mDimension : size;
    for (IndexType i = 0; i < size; i++)
        this->operator[](i) = rOtherCoordinates[i];
}

/// Destructor.
Point::~Point() {}

/// Assignment operator.
Point &Point::operator=(const Point &rOther)
{
    CoordinatesArrayType::operator=(rOther);
    return *this;
}

bool Point::operator==(const Point &rOther)
{
    return std::equal(this->begin(), this->end(), rOther.begin());
}

constexpr Point::IndexType Point::Dimension()
{
    return 3;
}

/** Returns X coordinate */
double Point::X() const
{
    return this->operator[](0);
}

/** Returns Y coordinate */
double Point::Y() const
{
    return this->operator[](1);
}

/** Returns Z coordinate */
double Point::Z() const
{
    return this->operator[](2);
}

double & Point::X()
{
    return this->operator[](0);
}

/** Returns Y coordinate */
double & Point::Y()
{
    return this->operator[](1);
}

/** Returns Z coordinate */
double & Point::Z()
{
    return this->operator[](2);
}

Point::CoordinatesArrayType const & Point::Coordinates() const
{
    return *this;
}

Point::CoordinatesArrayType & Point::Coordinates()
{
    return *this;
}

/// Turn back information as a string.
std::string Point::Info() const
{
    return "Point";
}

/// Print information about this object.
void Point::PrintInfo(std::ostream &rOStream) const
{
    rOStream << this->Info();
}

/// Print object's data.
void Point::PrintData(std::ostream &rOStream) const
{
    rOStream << "(" << this->operator[](0)
                    << this->operator[](1)  
                    << this->operator[](2) 
             << ")";
}

void Point::SetAllCoordinates(double const &Value)
{
    for (IndexType i = 0; i < mDimension; i++)
        this->operator()(i) = Value;
}

void Point::save(Serializer &rSerializer) const
{
    rSerializer.save_base("BaseClass", *static_cast<const array_1d<double, mDimension> *>(this));
}

void Point::load(Serializer &rSerializer)
{
    rSerializer.load_base("BaseClass", *static_cast<array_1d<double, mDimension> *>(this));
}

// Explicit instantiation of the component
template class KratosComponents<Point>;

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

} // namespace Kratos.
