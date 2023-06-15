//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/geometry.h"


namespace Kratos {

///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Kratos Classes
///@{

/// Tools to calculate the nearest point in different geometries
/** These tools are generic enough to be used in different contexts while used in the geometries
*/
class KRATOS_API(KRATOS_CORE) NearestPointUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestPointUtilities
    KRATOS_CLASS_POINTER_DEFINITION(NearestPointUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NearestPointUtilities() = delete;

    /// Copy constructor.
    NearestPointUtilities(NearestPointUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NearestPointUtilities& operator=(NearestPointUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{


        template<class TPointType, class TGeometryType>
        static Point LineNearestPoint(TPointType const& ThePoint, TGeometryType const& TheLine){
            KRATOS_DEBUG_ERROR_IF_NOT(TheLine.size() == 2) << "This function only accepts Line2D2 as input" << std::endl;
            Point result;
            const Point line_projected_point = GeometricalProjectionUtilities::FastProjectOnLine(TheLine, ThePoint, result);
            array_1d<double,3> projected_local;
            if(TheLine.IsInside(result.Coordinates(), projected_local))
                return result;
            
            double distance1 = norm_2(TheLine[0] - result);
            double distance2 = norm_2(TheLine[1] - result);

            result = (distance1 < distance2) ? TheLine[0] : TheLine[1];
            return result;
        }

    //  Dviding the plane of the triangle in 7 zones and find the nearest reflecting those zones
    //
    //             \       /
    //              \  6  /
    //               \   /
    //                \ /
    //                 /\ 
    //                /  \ 
    //          2    /    \     3
    //              /   1  \ 
    //             /        \ 
    //    _______ /__________\_____________ 
    //           /            \ 
    //     5    /       4      \    7
    //         /                \ 
    //        /                  \ 

        template<class TPointType, class TGeometryType>
        static Point TriangleNearestPoint(TPointType const& ThePoint, TGeometryType const& TheTriangle){
            constexpr double Tolerance = 1e-12;
            array_1d<double, 3> result;
            array_1d<double,3> local_coordinates = ZeroVector(3);
            const Point center = TheTriangle.Center();
            const array_1d<double, 3> normal = TheTriangle.UnitNormal(local_coordinates);
            double distance = 0.00;
            const Point point_projected = GeometricalProjectionUtilities::FastProject( center, ThePoint, normal, distance);
            TheTriangle.PointLocalCoordinates(local_coordinates, point_projected);

            if(local_coordinates[0] < -Tolerance) { // case 2,5,6
                if(local_coordinates[1] < -Tolerance) { // case 5
                    result = TheTriangle[0];
                }
                else if ((local_coordinates[0] + local_coordinates[1]) > (1.0+Tolerance)) { // case 6
                    result = TheTriangle[2];
                }
                else {
                    result = GetLineNearestPoint(ThePoint, Line3D2<Node<3>>(TheTriangle.pGetPoint(0), TheTriangle.pGetPoint(2)));
                }
            }
            else if(local_coordinates[1] < -Tolerance) { // case 4,7 (case 5 is already covered in previous if)
                if ((local_coordinates[0] + local_coordinates[1]) > (1.0+Tolerance)) { // case 7
                    result = TheTriangle[1];
                }
                else { // case 4
                    result = GetLineNearestPoint(ThePoint, Line3D2<Node<3>>(TheTriangle.pGetPoint(0), TheTriangle.pGetPoint(1)));
                }
            }
            else if ((local_coordinates[0] + local_coordinates[1]) > (1.0+Tolerance)) { // case 3
                result = GetLineNearestPoint(ThePoint, Line3D2<Node<3>>(TheTriangle.pGetPoint(1), TheTriangle.pGetPoint(2)));
            } 
            else {  // inside
                result = point_projected;
            }
                
            return result;
        }


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


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}

}; // Class NearestPointUtilities

///@}

///@name Type Definitions
///@{

///@}

///@} addtogroup block

} // namespace Kratos
