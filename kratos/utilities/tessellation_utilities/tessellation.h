//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

#if !defined(KRATOS_TESSELLATION_H_INCLUDED )
#define  KRATOS_TESSELLATION_H_INCLUDED

namespace Kratos {

template <class TContainerPointType>
class Tessellation
{
public:
    typedef Geometry<typename TContainerPointType::value_type> GeometryType;
    typedef typename GeometryType::CoordinatesArrayType CoordinatesArrayType;

    ///@name Constructors
    ///@{

    /// Conctructor for tessellation
    Tessellation()
    {
    }

    /// Copy Constructor
    Tessellation(Tessellation<TContainerPointType> const& rOther)
    {
    }

    /// Assignment Operator
    Tessellation<TContainerPointType>& operator=(const Tessellation<TContainerPointType>& rOther)
    {
        return *this;
    }

    ///@}
    ///@name Constructors
    ///@{

    virtual void Tessellate(
        const GeometryType& rGeometry,
        const double Tolerance,
        const int NumberOfGuessesPerInterval = 1)
    {
        KRATOS_ERROR << "Calling Tessellate from base Tessellation class."
            << std::endl;
    }

    /* @brief This method returns closest point within the tesselation.
     */
    virtual void GetClosestPoint(
        const CoordinatesArrayType& rPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointGlobalCoordinates,
        CoordinatesArrayType& rClosestPointLocalCoordinates) const
    {
        KRATOS_ERROR << "Calling GetClosestPoint from base Tessellation class."
            << std::endl;
    }

};

} // namespace Tessellation

#endif // KRATOS_TESSELLATION_H_INCLUDED defined
