//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Veronika Singer
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/particle_based_conditions/mpm_particle_base_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
    double MPMParticleBaseLoadCondition::GetPointLoadIntegrationWeight()
    {
        return 1.0;
    }

    /**
     * Shape function values in given point. This method calculate the shape function
     * vector in given point.
     *
     * @param rPoint point which shape function values have to
     * be calculated in it.
     *
     * @return Vector of double which is shape function vector \f$ N \f$ in given point.
     *
    */
    Vector& MPMParticleBaseLoadCondition::MPMShapeFunctionPointValues( Vector& rResult, const array_1d<double,3>& rPoint )
    {
        KRATOS_TRY

        rResult = MPMParticleBaseCondition::MPMShapeFunctionPointValues(rResult, rPoint);

        // Additional check to eliminate loss of point load quantity
        GeometryType& r_geometry = GetGeometry();
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();

        double denominator = 1.0;
        for ( unsigned int i = 0; i < number_of_nodes; i++ )
        {
            if (r_geometry[i].FastGetSolutionStepValue(NODAL_MASS, 0) <= std::numeric_limits<double>::epsilon()){
                denominator -= rResult[i];
                rResult[i] = 0;
            }
        }

        rResult = rResult/denominator;

        return rResult;

        KRATOS_CATCH( "" )
    }



} // Namespace Kratos
