//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (based on work of Pablo Becker)
//


// System includes


// External includes


// Project includes
#include "point_locator.h"


namespace Kratos
{

    bool PointLocator::Find(const Point& rThePoint)
    {
        KRATOS_ERROR_IF(mIsInitalized) << "This instance can only be used for one Point" << std::endl;

        // note that this cannot be omp bcs breaking is not allowed in omp
        for (const auto& r_elem : mrModelPart.Elements())
        {
            // if found
            //     break
        }
    }



    void PointLocator::InterpolateValue(const Variable<double>& rVariable, double& rValue)
    {

    }



}  // namespace Kratos.


