// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(POINT_BELONG_DEFINED )
#define  POINT_BELONG_DEFINED

// System includes
#include <iostream>
#include <vector>
#include "boost/smart_ptr.hpp"

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"

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
    enum PointBelongs {Master = 0, Slave = 1, Intersection = 2};
    enum PointBelongsLine2D2N {MasterLine2D2N0 = 0, MasterLine2D2N1 = 1, SlaveLine2D2N0 = 2, SlaveLine2D2N1 = 3, IntersectionLine2D2N = 4};
    enum PointBelongsTriangle3D3N {MasterTriangle3D3N0 = 0, MasterTriangle3D3N1 = 1, MasterTriangle3D3N2 = 2, SlaveTriangle3D3N0 = 3, SlaveTriangle3D3N1 = 4, SlaveTriangle3D3N2 = 5, IntersectionTriangle3D3N = 6};
    enum PointBelongsQuadrilateral3D4N {MasterQuadrilateral3D4N0 = 0, MasterQuadrilateral3D4N1 = 1, MasterQuadrilateral3D4N2 = 2, MasterQuadrilateral3D4N3 = 3, SlaveQuadrilateral3D4N0 = 4, SlaveQuadrilateral3D4N1 = 5, SlaveQuadrilateral3D4N2 = 6, SlaveQuadrilateral3D4N3 = 7, IntersectionQuadrilateral3D4N = 8};
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
class PointBelong : public PointType
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
        PointType()
    {}

    PointBelong(const array_1d<double, 3> Coords):
        PointType(Coords)
    {}
    
    PointBelong(const array_1d<double, 3> Coords, const BelongType& ThisBelongs):
        PointType(Coords),
        mBelongs(ThisBelongs)
    {}
    
    /// Destructor.
    ~PointBelong() override{}
    
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
