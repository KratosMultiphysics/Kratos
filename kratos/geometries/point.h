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
class KRATOS_API(KRATOS_CORE) Point : public array_1d<double, 3>
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
    Point();

    /// 3d constructor.
    Point(double NewX, double NewY = 0, double NewZ = 0);

    /** Copy constructor. Initialize this point with the coordinates
    of given point.*/
    Point(Point const &rOtherPoint);

    /** Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array. */
    Point(CoordinatesArrayType const &rOtherCoordinates);

    /** Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array. */
    template <class TVectorType>
    Point(vector_expression<TVectorType> const &rOtherCoordinates) : BaseType(rOtherCoordinates) {}

    /** Constructor using coordinates stored in given std::vector. Initialize
    this point with the coordinates in the array. */
    Point(std::vector<double> const &rOtherCoordinates);

    /// Destructor.
    virtual ~Point();

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Point &operator=(const Point &rOther);
    bool operator==(const Point &rOther);

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    static constexpr IndexType Dimension();

    /** Returns X coordinate */
    double X() const;

    /** Returns Y coordinate */
    double Y() const;

    /** Returns Z coordinate */
    double Z() const;

    /** Returns X coordinate */
    double &X();

    /** Returns Y coordinate */
    double &Y();

    /** Returns Z coordinate */
    double &Z();

    CoordinatesArrayType const &Coordinates() const;
    CoordinatesArrayType &Coordinates();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const;

    ///@}

  private:
    ///@name Private Operations
    ///@{

    void SetAllCoordinates(double const &Value = double());

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer &rSerializer) const;
    virtual void load(Serializer &rSerializer);

    ///@}

}; // Class Point

///@}

// Explicit instantiation declaration
extern template class KRATOS_API(KRATOS_CORE) KratosComponents<Point>;

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_POINT_H_INCLUDED  defined
