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

#pragma once

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/serializer.h"

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

/**
 * @class Point
 * @ingroup KratosCore
 * @brief Point class.
 * @details Point class. Stores coordinates of a point and have some basic operations defined.
 * @see Geometry
 * @see Node
 * @see IntegrationPoint
 * @author Riccardo Rossi
*/
class Point : public array_1d<double, 3>
{
    /// Dimension definition
    static constexpr std::size_t mDimension = 3;

public:
    ///@name Constants
    ///@{

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Point
    KRATOS_CLASS_POINTER_DEFINITION(Point);

    /// The base type definition
    using BaseType = array_1d<double, mDimension>;

    /// The coordinates array type definition
    using CoordinatesArrayType = BaseType;

    /// The size type definition
    using SizeType = std::size_t;

    /// The index type definition
    using IndexType = std::size_t;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Point();

    /**
     * @brief 3d constructor.
     * @param NewX New X coordinate
     * @param NewY New Y coordinate
     * @param NewZ New Z coordinate
     */
    Point(
        const double NewX,
        const double NewY = 0,
        const double NewZ = 0
        );

    /**
     * @brief Copy constructor. Initialize this point with the coordinates
    of given point.
     * @param rOtherPoint The other point to copy
    */
    Point(const Point& rOtherPoint);

    /**
     * @brief Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array.
     * @param rOtherCoordinates The coordinates to be used
     */
    explicit Point(const CoordinatesArrayType& rOtherCoordinates);

    /**
     * @brief Constructor using coordinates stored in given array. Initialize
    this point with the coordinates in the array. */
    template <class TVectorType>
    explicit Point(const vector_expression<TVectorType>& rOtherCoordinates)
        : BaseType(rOtherCoordinates) {}

    /**
     * @brief Constructor using coordinates stored in given std::vector. Initialize
    this point with the coordinates in the array.
     * @param rOtherCoordinates The coordinates to be used
     */
    explicit Point(const std::vector<double>& rOtherCoordinates);

    /// Destructor.
    virtual ~Point() {}

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Point &operator=(const Point &rOther);

    /// == operator.
    bool operator==(const Point &rOther) const;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method computes the distance between this point and another one (squared)
     * @details In order to avoid square root operation if faster computation is needed
     * @param rOtherPoint The other point to compute the distance
     * @return The squared distance between this and another point
     */
    double SquaredDistance(const Point& rOtherPoint) const;

    /**
     * @brief This method computes the distance between this point and another one
     * @details Using norm_2 to take benefic of the SIMD optimization of the library
     * @param rOtherPoint The other point to compute the distance
     * @return The distance between this and another point
     */
    double Distance(const Point& rOtherPoint) const;

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Returns the dimension of the point
     * @details It is always 3
     * @return The dimension of the point
     */
    static constexpr IndexType Dimension()
    {
        return mDimension;
    }

    /**
     * @brief Returns X coordinate
     * @details const version
     * @return X coordinate
     */
    double X() const;

    /**
     * @brief Returns Y coordinate
     * @details const version
     * @return Z coordinate
     */
    double Y() const;

    /**
     * @brief Returns Z coordinate
     * @details const version
     * @return Z coordinate
     */
    double Z() const;

    /**
     * @brief Returns X coordinate
     * @details reference version
     * @return X coordinate
     */
    double& X();

    /**
     * @brief Returns Y coordinate
     * @details reference version
     * @return Y coordinate
     */
    double& Y();

    /**
     * @brief Returns Z coordinate
     * @details reference version
     * @return Z coordinate
     */
    double& Z();

    /**
     * @brief Returns the coordinates
     * @details const version
     * @return The coordinates
     */
    const CoordinatesArrayType& Coordinates() const;

    /**
     * @brief Returns the coordinates
     * @details reference version
     * @return The coordinates
     */
    CoordinatesArrayType& Coordinates();

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

    /**
     * @brief Set all coordinates
     * @param Value The new value
     */
    void SetAllCoordinates(const double Value = double());

    ///@}
    ///@name Serialization
    ///@{

    /// Friend class for serialization
    friend class Serializer;

    /**
     * @brief Serialization save function
     * @param rSerializer The serializer
     */
    virtual void save(Serializer& rSerializer) const;

    /**
     * @brief Serialization load function
     * @param rSerializer The serializer
     */
    virtual void load(Serializer& rSerializer);

    ///@}
}; // Class Point

///@}

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
