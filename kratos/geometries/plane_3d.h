//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//  contributors:    Ruben Zorrilla
//

#if !defined(KRATOS_PLANE_3D_H_INCLUDED)
#define  KRATOS_PLANE_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/point.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"

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
*/
class Plane3D
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Point2D
    KRATOS_CLASS_POINTER_DEFINITION(Plane3D);

    ///@}
    ///@name Life Cycle
    ///@{

    Plane3D() = delete;

    Plane3D(array_1d<double, 3> const &rNormal, double DistanceToOrigin) : mD(DistanceToOrigin), mNormal(rNormal) {}

    Plane3D(array_1d<double, 3> const &rNormal, const Point &rReferencePoint)
    {
        // Compute the unit normal
        mNormal = rNormal;
        auto normal_length = norm_2(mNormal);
        KRATOS_DEBUG_CHECK_GREATER(normal_length, std::numeric_limits<double>::epsilon());
        mNormal /= normal_length;
        // Compute the distance to origin
        mD = -inner_prod(mNormal,rReferencePoint);
    }   

    Plane3D(const Point &Point1, const Point &Point2, const Point &Point3)
    {
        // Compute the unit normal
        array_1d<double,3> v_1 = Point2 - Point1;
        array_1d<double,3> v_2 = Point3 - Point1;
        MathUtils<double>::CrossProduct(mNormal, v_1, v_2);
        auto normal_length = norm_2(mNormal);
        KRATOS_DEBUG_CHECK_GREATER(normal_length, std::numeric_limits<double>::epsilon());
        mNormal /= normal_length;
        // Compute the distance to origin
        mD = -inner_prod(mNormal, Point1);
    }

    /// Destructor. Do nothing!!!
    ~Plane3D() {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /** Return the plane normal vector
    @return Array containing the plane normal
    */
    array_1d<double,3> const& GetNormal()
    { 
        return mNormal; 
    }

    /** Return the plane distance value
    @return Plane distance value
    */
    double GetDistanceToOrigin()
    { 
        return mD; 
    }

    /** Calculates the plane signed distance value
    @return Plane signed distance value
    */
    double CalculateSignedDistance(Point const &rPoint)
    {
        return inner_prod(mNormal, rPoint) + mD;
    }

    /** Turn back information as a string.
    @return String contains information about this geometry.
    @see PrintData()
    @see PrintInfo()
    */
    std::string Info() const
    {
        return "a 3D plane auxiliar class";
    }

    /** Print information about this object.
    @param rOStream Stream to print into it.
    @see PrintData()
    @see Info()
    */
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "a 3D plane auxiliar class";
    }

    /** Print geometry's data into given stream. Prints it's points
    by the order they stored in the geometry and then center
    point of geometry.
    @param rOStream Stream to print into it.
    @see PrintInfo()
    @see Info()
    */
    void PrintData(std::ostream& rOStream) const
    {
        rOStream << "a 3D plane auxiliar class with:\n";
        rOStream << "\t- distance to origin: " << mD << "\n";
        rOStream << "\t- normal: (" << mNormal[0] << " , " << mNormal[1] << " , " << mNormal[2] << ")";
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:

    ///@name Member Variables
    ///@{

    double                       mD;
    array_1d<double,3>      mNormal;
   
    ///@}

}; // Class Plane3D

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  Plane3D& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Plane3D& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_PLANE_3D_H_INCLUDED  defined 


