// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 

#if !defined(POINT_BELONG_DEFINED )
#define  POINT_BELONG_DEFINED

// System includes

// External includes

// Project includes
#include "geometries/point.h"

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
    
#if !defined(POINT_BELONGS)
#define POINT_BELONGS
    
    enum PointBelongs 
    {
        Master       = 0, 
        Slave        = 1, 
        Intersection = 2
    };
    
    enum PointBelongsLine2D2N 
    {
        MasterLine2D2N0      = 0, 
        MasterLine2D2N1      = 1, 
        SlaveLine2D2N0       = 2, 
        SlaveLine2D2N1       = 3, 
        IntersectionLine2D2N = 4
    };
    
    enum PointBelongsTriangle3D3N 
    {
        MasterTriangle3D3N0          = 0, 
        MasterTriangle3D3N1          = 1, 
        MasterTriangle3D3N2          = 2, 
        SlaveTriangle3D3N0           = 3, 
        SlaveTriangle3D3N1           = 4, 
        SlaveTriangle3D3N2           = 5, 
        IntersectionTriangle3D3N     = 6,
        IntersectionTriangle3D3N0101 = 10106,
        IntersectionTriangle3D3N1001 = 10016,
        IntersectionTriangle3D3N1201 = 10216,
        IntersectionTriangle3D3N2101 = 10126,
        IntersectionTriangle3D3N0201 = 10206,
        IntersectionTriangle3D3N2001 = 10026,
        IntersectionTriangle3D3N0110 = 1106,
        IntersectionTriangle3D3N1010 = 1016,
        IntersectionTriangle3D3N1210 = 1216,
        IntersectionTriangle3D3N2110 = 1126,
        IntersectionTriangle3D3N0210 = 1206,
        IntersectionTriangle3D3N2010 = 1026,
        IntersectionTriangle3D3N0112 = 21106,
        IntersectionTriangle3D3N1012 = 21016,
        IntersectionTriangle3D3N1212 = 21216,
        IntersectionTriangle3D3N2112 = 21126,
        IntersectionTriangle3D3N0212 = 21206,
        IntersectionTriangle3D3N2012 = 21026,
        IntersectionTriangle3D3N0121 = 12106,
        IntersectionTriangle3D3N1021 = 12016,
        IntersectionTriangle3D3N1221 = 12216,
        IntersectionTriangle3D3N2121 = 12126,
        IntersectionTriangle3D3N0221 = 12206,
        IntersectionTriangle3D3N2021 = 12026,
        IntersectionTriangle3D3N0102 = 20106,
        IntersectionTriangle3D3N1002 = 20016,
        IntersectionTriangle3D3N1202 = 20216,
        IntersectionTriangle3D3N2102 = 20126,
        IntersectionTriangle3D3N0202 = 20206,
        IntersectionTriangle3D3N2002 = 20026,
        IntersectionTriangle3D3N0120 = 2106,
        IntersectionTriangle3D3N1020 = 2016,
        IntersectionTriangle3D3N1220 = 2216,
        IntersectionTriangle3D3N2120 = 2126,
        IntersectionTriangle3D3N0220 = 2206,
        IntersectionTriangle3D3N2020 = 2026
    };
    
    enum PointBelongsQuadrilateral3D4N 
    {
        MasterQuadrilateral3D4N0          = 0, 
        MasterQuadrilateral3D4N1          = 1, 
        MasterQuadrilateral3D4N2          = 2, 
        MasterQuadrilateral3D4N3          = 3, 
        SlaveQuadrilateral3D4N0           = 4, 
        SlaveQuadrilateral3D4N1           = 5, 
        SlaveQuadrilateral3D4N2           = 6, 
        SlaveQuadrilateral3D4N3           = 7, 
        IntersectionQuadrilateral3D4N     = 8,
        IntersectionQuadrilateral3D4N0101 = 10108,
        IntersectionQuadrilateral3D4N1001 = 10018,
        IntersectionQuadrilateral3D4N1201 = 10218,
        IntersectionQuadrilateral3D4N2101 = 10128,
        IntersectionQuadrilateral3D4N2301 = 10328,
        IntersectionQuadrilateral3D4N3201 = 10238,
        IntersectionQuadrilateral3D4N3001 = 10038,
        IntersectionQuadrilateral3D4N0301 = 10308,
        IntersectionQuadrilateral3D4N0110 = 1108,
        IntersectionQuadrilateral3D4N1010 = 1018,
        IntersectionQuadrilateral3D4N1210 = 1218,
        IntersectionQuadrilateral3D4N2110 = 1128,
        IntersectionQuadrilateral3D4N2310 = 1328,
        IntersectionQuadrilateral3D4N3210 = 1238,
        IntersectionQuadrilateral3D4N3010 = 1038,
        IntersectionQuadrilateral3D4N0310 = 1308,
        IntersectionQuadrilateral3D4N0112 = 21108,
        IntersectionQuadrilateral3D4N1012 = 21018,
        IntersectionQuadrilateral3D4N1212 = 21218,
        IntersectionQuadrilateral3D4N2112 = 21128,
        IntersectionQuadrilateral3D4N2312 = 21328,
        IntersectionQuadrilateral3D4N3212 = 21238,
        IntersectionQuadrilateral3D4N3012 = 21038,
        IntersectionQuadrilateral3D4N0312 = 21308,
        IntersectionQuadrilateral3D4N0121 = 12108,
        IntersectionQuadrilateral3D4N1021 = 12018,
        IntersectionQuadrilateral3D4N1221 = 12218,
        IntersectionQuadrilateral3D4N2121 = 12128,
        IntersectionQuadrilateral3D4N2321 = 12328,
        IntersectionQuadrilateral3D4N3221 = 12238,
        IntersectionQuadrilateral3D4N3021 = 12038,
        IntersectionQuadrilateral3D4N0321 = 12308,
        IntersectionQuadrilateral3D4N0123 = 32108,
        IntersectionQuadrilateral3D4N1023 = 32018,
        IntersectionQuadrilateral3D4N1223 = 32218,
        IntersectionQuadrilateral3D4N2123 = 32128,
        IntersectionQuadrilateral3D4N2323 = 32328,
        IntersectionQuadrilateral3D4N3223 = 32238,
        IntersectionQuadrilateral3D4N3023 = 32038,
        IntersectionQuadrilateral3D4N0323 = 32308,
        IntersectionQuadrilateral3D4N0132 = 23108,
        IntersectionQuadrilateral3D4N1032 = 23018,
        IntersectionQuadrilateral3D4N1232 = 23218,
        IntersectionQuadrilateral3D4N2132 = 23128,
        IntersectionQuadrilateral3D4N2332 = 23328,
        IntersectionQuadrilateral3D4N3232 = 23238,
        IntersectionQuadrilateral3D4N3032 = 23038,
        IntersectionQuadrilateral3D4N0332 = 23308,
        IntersectionQuadrilateral3D4N0130 = 3108,
        IntersectionQuadrilateral3D4N1030 = 3018,
        IntersectionQuadrilateral3D4N1230 = 3218,
        IntersectionQuadrilateral3D4N2130 = 3128,
        IntersectionQuadrilateral3D4N2330 = 3328,
        IntersectionQuadrilateral3D4N3230 = 3238,
        IntersectionQuadrilateral3D4N3030 = 3038,
        IntersectionQuadrilateral3D4N0330 = 3308,
        IntersectionQuadrilateral3D4N0103 = 30108,
        IntersectionQuadrilateral3D4N1003 = 30018,
        IntersectionQuadrilateral3D4N1203 = 30218,
        IntersectionQuadrilateral3D4N2103 = 30128,
        IntersectionQuadrilateral3D4N2303 = 30328,
        IntersectionQuadrilateral3D4N3203 = 30238,
        IntersectionQuadrilateral3D4N3003 = 30038,
        IntersectionQuadrilateral3D4N0303 = 30308
    };
    
#endif
    
///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/** @brief Custom Point container to be used by the mapper
 */

template<unsigned int TNumNodes>
class PointBelong : public Point<3>
{
public:
    ///@name Type Definitions
    ///@{
    
    typedef typename std::conditional<TNumNodes == 2, PointBelongsLine2D2N, typename std::conditional<TNumNodes == 3, PointBelongsTriangle3D3N, PointBelongsQuadrilateral3D4N>::type>::type BelongType;
    
    /// Counted pointer of PointBelong
    KRATOS_CLASS_POINTER_DEFINITION( PointBelong );

    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructors
    PointBelong():
        Point<3>()
    {}

    PointBelong(const array_1d<double, 3> Coords):
        Point<3>(Coords)
    {}
    
    PointBelong(const array_1d<double, 3> Coords, const BelongType& ThisBelongs):
        Point<3>(Coords),
        mBelongs(ThisBelongs)
    {}
    
    /// Destructor.
    ~PointBelong() override= default;
    
    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This method allows to set where the point belongs
     */
    void SetBelong(BelongType ThisBelongs)
    {
        mBelongs = ThisBelongs;
    }
    
    /**
     * This method recovers where the point belongs
     */
    BelongType GetBelong() const
    {
        return mBelongs;
    }
    
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

    BelongType mBelongs; // To know if the point belongs to the master/slave/intersection (just 3D) side          

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class PointBelong 

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // POINT_BELONG_DEFINED  defined
