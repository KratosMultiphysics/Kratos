// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//
#include "fix_model_part_above_phreatic_line_process.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

FixModelPartAbovePhreaticLineProcess::FixModelPartAbovePhreaticLineProcess(ModelPart& rModelPart, const Parameters& rProcessSettings)
    :  ApplyConstantPhreaticMultiLinePressureProcess(rModelPart, rProcessSettings)

{
    mMoveMeshFlag = rProcessSettings["move_mesh_flag"].GetBool();
}

int FixModelPartAbovePhreaticLineProcess::FindIndex(const Node& rNode) const
{
    array_1d<double, 3> coords = rNode.Coordinates();
    if (!mMoveMeshFlag) {
        std::transform(coords.begin(),
            coords.end(), rNode.GetSolutionStepValue(TOTAL_DISPLACEMENT).begin(),
            coords.begin(), std::plus<>());
    }
    const auto  number_of_coordinates = static_cast<int>(HorizontalDirectionCoordinates().size());
    for (int index = 0; index < number_of_coordinates; ++index) {
        if (HorizontalDirectionCoordinates()[index] >= coords[HorizontalDirection()]) {
            return index == 0 ? index : index - 1;
        }
    }

    return number_of_coordinates - 2;
}


double FixModelPartAbovePhreaticLineProcess::CalculateDistanceToPhreaticLine(const Node& rNode) const
{
    // find nodes in horizontalDirectionCoordinates
    const int firstPointIndex  = FindIndex(rNode);
    const int secondPointIndex = firstPointIndex + 1;

    array_1d<double, 2> y;
    y[0] = GravityDirectionCoordinates()[firstPointIndex] ;
    y[1] = GravityDirectionCoordinates()[secondPointIndex];

    const double slope = (y[1] - y[0]) / (HorizontalDirectionCoordinates()[secondPointIndex] -
                                          HorizontalDirectionCoordinates()[firstPointIndex]);

    auto adjustedPosition = rNode.Coordinates();

    if (!mMoveMeshFlag) {
        adjustedPosition = rNode.Coordinates() + rNode.GetSolutionStepValue(TOTAL_DISPLACEMENT);
    }

    const double height = slope * (adjustedPosition[HorizontalDirection()] -
                                   HorizontalDirectionCoordinates()[firstPointIndex]) + y[0];
    const double distance = height - adjustedPosition[GravityDirection()] ;
    return distance;
}
void FixModelPartAbovePhreaticLineProcess::ExecuteInitialize()
{
}
void FixModelPartAbovePhreaticLineProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY
    // create new model part
    KRATOS_INFO("Fixing Model Part Above (SOLUTION STEP)")<< mrModelPart.Name() <<std::endl;
    const Variable<double>& var = KratosComponents<Variable<double>>::Get(VariableName());

    // loop over nodes
    block_for_each(mrModelPart.Nodes(), [&var, this](Node& rNode) {
            const double distance = CalculateDistanceToPhreaticLine(rNode);
            if (distance <= 1.0e-6) {
                rNode.Free(var);
                rNode.SetValue(var, 0.0);
                rNode.Fix(var);
                KRATOS_INFO("Fixing Node:")<< rNode.Id() << " Distance: " << distance <<
                    "TOTAL_DISPLACEMENT:" << rNode.GetSolutionStepValue(TOTAL_DISPLACEMENT) <<
                    "DISPLACEMENT" << rNode.GetSolutionStepValue(DISPLACEMENT) <<  std::endl;
                }
            else {
                rNode.Free(var);
            }
        });
    KRATOS_CATCH("")
}

FixModelPartAbovePhreaticLineProcess::~FixModelPartAbovePhreaticLineProcess() = default;

} // namespace Kratos
