//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(SCALAR_CO_SOLVING_UTILITIES_H_INCLUDED)
#define SCALAR_CO_SOLVING_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes

namespace Kratos
{
class ScalarCoSolvingUtilities
{
public:
    ScalarCoSolvingUtilities() : mEchoLevel(0)
    {
    }

    ScalarCoSolvingUtilities(int EchoLevel) : mEchoLevel(EchoLevel)
    {
    }

    bool IsSlipConditionsUsed(const ModelPart::NodesContainerType& rNodes, const Flags& rFlag)
    {
        const int number_of_nodes = rNodes.size();
        unsigned int number_of_slip_condition_nodes = 0;
#pragma omp parallel for reduction(+ : number_of_slip_condition_nodes)
        for (int i_node = 0; i_node < number_of_nodes; ++i_node)
        {
            const ModelPart::NodeType& r_node = *(rNodes.begin() + i_node);
            if (r_node.Is(rFlag))
            {
                const array_1d<double, 3>& r_velocity =
                    r_node.FastGetSolutionStepValue(VELOCITY);
                const double velocity_norm = norm_2(r_velocity);

                if (velocity_norm > std::numeric_limits<double>::epsilon())
                {
                    number_of_slip_condition_nodes++;
                }
            }
        }

        if (number_of_slip_condition_nodes > 0)
        {
            KRATOS_INFO_IF("IsSlipConditionsUsed", mEchoLevel > 0)
                << "Slip condition is being used by "
                << number_of_slip_condition_nodes << " nodes.\n";
            return true;
        }
        else
        {
            return false;
        }
    }

private:
    int mEchoLevel;
};
} // namespace Kratos

#endif