//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// Project includes 
#include "includes/define.h"
#include "includes/variables.h" 
#include "includes/model_part.h"


namespace Kratos
{
    // This class contains utilities which are used for the User Vehicle model
    class UvecUtilities 
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(UvecUtilities);

        UvecUtilities() = default;


        /// <summary>
        /// Sets the load on the first condition that has a non-zero value for the specified variable.
        /// </summary>
        /// <param name="rModelPart"></param>
        /// <param name="rVariable"></param>
        /// <param name="Precision"></param>
        static void SetLoadOnCondition(ModelPart& rModelPart,const Variable<array_1d<double, 3>>& rVariable, double Precision)
        {
            const array_1d<double, 3>& global_value = rModelPart.GetValue(rVariable);

            for (auto& r_condition : rModelPart.Conditions())
            {
                const array_1d<double, 3>& cond_value = r_condition.GetValue(rVariable);
                
                // Check if any component is greater than precision
                if (std::abs(cond_value[0]) >= Precision ||
                    std::abs(cond_value[1]) >= Precision ||
                    std::abs(cond_value[2]) >= Precision)
                {
                    r_condition.SetValue(rVariable, global_value);

                    break; // Exit the loop after setting the value
                }
            }
        }
        
        /// <summary>
        /// Gets the value of a moving condition variable by summing the values of all active conditions in the model part. This assumes there is only one active condition with non zero values
        /// </summary>
        /// <param name="rModelPart"></param>
        /// <param name="rVariable"></param>
        /// <returns></returns>
        static array_1d<double, 3> GetMovingConditionVariable(ModelPart& rModelPart, const Variable<array_1d<double, 3>>& rVariable)
        {
            array_1d<double, 3> result = ZeroVector(3);
            for (const auto& r_condition : rModelPart.Conditions())
            {
                if (r_condition.IsActive())
                {
                    const array_1d<double, 3>& cond_value = r_condition.GetValue(rVariable);
                    result += cond_value;
                }
            }
            return result;
        }

    
    protected:


    private:


    };

}  // namespace Kratos.