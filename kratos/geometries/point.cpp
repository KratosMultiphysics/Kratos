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
//  Contributors:    Hoang Giang Bui
//                   Josep Maria Carbonell
//

// System includes

// External includes

// Project includes
#include "geometries/point.h"

namespace Kratos
{

Point::Point() : BaseType()
{
    SetAllCoordinates();
}

/***********************************************************************************/
/***********************************************************************************/

Point::Point(
    const double NewX,
    const double NewY,
    const double NewZ
    ) : BaseType()
{
    this->operator()(0) = NewX;
    this->operator()(1) = NewY;
    this->operator()(2) = NewZ;
}

/***********************************************************************************/
/***********************************************************************************/

Point::Point(const Point& rOtherPoint)
    : BaseType(rOtherPoint)
{

}

/***********************************************************************************/
/***********************************************************************************/

Point::Point(const CoordinatesArrayType& rOtherCoordinates)
    : BaseType(rOtherCoordinates)
{

}

/***********************************************************************************/
/***********************************************************************************/

Point::Point(const std::vector<double>& rOtherCoordinates)
    : BaseType()
{
    SizeType size = rOtherCoordinates.size();
    size = (mDimension < size) ? mDimension : size;
    for (IndexType i = 0; i < size; i++)
        this->operator[](i) = rOtherCoordinates[i];
}

/***********************************************************************************/
/***********************************************************************************/

Point& Point::operator=(const Point &rOther)
{
    CoordinatesArrayType::operator=(rOther);
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

bool Point::operator==(const Point &rOther) const
{
    return std::equal(this->begin(), this->end(), rOther.begin());
}

/***********************************************************************************/
/***********************************************************************************/

double Point::SquaredDistance(const Point& rOtherPoint) const
{
    const array_1d<double, 3> diff_vector = this->Coordinates() - rOtherPoint.Coordinates();
    return (std::pow(diff_vector[0], 2) + std::pow(diff_vector[1], 2) + std::pow(diff_vector[2], 2));
}

/***********************************************************************************/
/***********************************************************************************/

double Point::Distance(const Point& rOtherPoint) const
{
    return norm_2(this->Coordinates() - rOtherPoint.Coordinates());
}

/***********************************************************************************/
/***********************************************************************************/

double Point::X() const
{
    return this->operator[](0);
}

/***********************************************************************************/
/***********************************************************************************/

double Point::Y() const
{
    return this->operator[](1);
}

/***********************************************************************************/
/***********************************************************************************/

double Point::Z() const
{
    return this->operator[](2);
}

/***********************************************************************************/
/***********************************************************************************/

double& Point::X()
{
    return this->operator[](0);
}

/***********************************************************************************/
/***********************************************************************************/

double& Point::Y()
{
    return this->operator[](1);
}

/***********************************************************************************/
/***********************************************************************************/

double& Point::Z()
{
    return this->operator[](2);
}

/***********************************************************************************/
/***********************************************************************************/

const Point::CoordinatesArrayType& Point::Coordinates() const
{
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

Point::CoordinatesArrayType& Point::Coordinates()
{
    return *this;
}

/***********************************************************************************/
/***********************************************************************************/

std::string Point::Info() const
{
    return "Point";
}

/***********************************************************************************/
/***********************************************************************************/

void Point::PrintInfo(std::ostream &rOStream) const
{
    rOStream << this->Info();
}

/***********************************************************************************/
/***********************************************************************************/

void Point::PrintData(std::ostream &rOStream) const
{
    rOStream << " ("  << this->operator[](0)
             << ", " << this->operator[](1)
             << ", " << this->operator[](2)
             << ")";
}

/***********************************************************************************/
/***********************************************************************************/

void Point::SetAllCoordinates(const double Value)
{
    for (IndexType i = 0; i < mDimension; i++)
        this->operator()(i) = Value;
}

/***********************************************************************************/
/***********************************************************************************/

void Point::save(Serializer &rSerializer) const
{
    rSerializer.save_base("BaseClass", *static_cast<const array_1d<double, mDimension> *>(this));
}

/***********************************************************************************/
/***********************************************************************************/

void Point::load(Serializer &rSerializer)
{
    rSerializer.load_base("BaseClass", *static_cast<array_1d<double, mDimension> *>(this));
}

} // namespace Kratos.
