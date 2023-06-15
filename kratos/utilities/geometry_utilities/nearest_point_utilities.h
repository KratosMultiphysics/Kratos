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
