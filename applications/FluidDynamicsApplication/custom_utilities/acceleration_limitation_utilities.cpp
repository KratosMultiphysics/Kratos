//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "utilities/openmp_utils.h"
#include "utilities/variable_utils.h"

// Application includes
#include "acceleration_limitation_utilities.h"

namespace Kratos
{

    /* Public functions *******************************************************/

    /**
     * @brief Iteration over all nodes in ModelPart to find and reduce excessive accelerations.
     * 
     * @param rModelPart model part to be controled
     */

    void AccelerationLimitationUtilities::LimitAccelerationAtNodes(ModelPart &rModelPart) {
        // Sum the reactions in the model part of interest.
        // Note that the reactions are assumed to be already computed.

        double dt = rModelPart.GetProcessInfo().GetValue(DELTA_TIME);
        double time = rModelPart.GetProcessInfo().GetValue(TIME);

        // iterating over all elements
        #pragma omp parallel for schedule(dynamic) shared(dt, time)
        for (int i = 0; i < static_cast<int>(rModelPart.Elements().size()); ++i){
            
            auto it_elem = rModelPart.ElementsBegin() + i;
            // iterating over all nodes of the current element
            for (unsigned int nnode = 0; nnode < it_elem->GetGeometry().size(); nnode++){

                // retrieving velocities of current and last time step
                array_1d<double, 3> v  = it_elem->GetGeometry()[nnode].FastGetSolutionStepValue( Kratos::VELOCITY, 0 );
                array_1d<double, 3> vn = it_elem->GetGeometry()[nnode].FastGetSolutionStepValue( Kratos::VELOCITY, 1 );

                array_1d<double, 3> delta_v = v - vn;
                double norm_delta_v = sqrt( delta_v[0]*delta_v[0] + delta_v[1]*delta_v[1] + delta_v[2]*delta_v[2] );

                // a multiple of the 
                double alpha = norm_delta_v / ( dt * mMaximalAccelaration * 9.81 );

                if ( alpha > 1.0){
                    // setting a new and "reasonable" velocity by scaling
                    v = vn + ( 1 / alpha ) * delta_v;
                    it_elem->GetGeometry()[nnode].FastGetSolutionStepValue( Kratos::VELOCITY, 0 ) = v;
                    // std::cout << "AccelerationLimitationUtilities : High accelaration was detected and reduced" << std::endl;

                    std::ofstream myfile;
                    myfile.open ("accelerationLog.txt", std::ios::out | std::ios::app);
                    if ( myfile.is_open() ){
                        myfile << "time = " << time << " s   |   node = " << it_elem->GetGeometry()[nnode].Id() << "   :  ";
                        myfile << " High accelaration was detected and reduced \n";
                        myfile.close();
                    }
                }
            }       
        }
    }

    /**
     * @brief Set a new value for the maximal acceleration
     * The value is given as a multiple of the gravitational acceleration (1g = 9.81 m/s²) 
     *  
     * @param newMaxAcc new maximal acceleration
     */

    void AccelerationLimitationUtilities::SetLimitAsMultipleOfGravitionalAcceleration( double& newMaxAcc ){

        mMaximalAccelaration = newMaxAcc;
        std::cout << "Maximal acceleration is limited to " << newMaxAcc * 9.81 << " m/s2" << std::endl;

    }

    /**
     * @brief Execute the utility based on already given data
     * 
     */
    void AccelerationLimitationUtilities::Execute(){

        LimitAccelerationAtNodes( *mpModelPart );
    }


    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const AccelerationLimitationUtilities& rThis ) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
