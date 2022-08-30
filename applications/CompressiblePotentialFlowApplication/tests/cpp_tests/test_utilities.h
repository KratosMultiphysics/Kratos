//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Inigo Lopez and Marc Núñez
//

#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
    // forward-declaring
    class ModelPart;
    class Element;

namespace PotentialFlowTestUtilities
{

template <int NumNodes>
void AssignPotentialsToNormalElement(Element& rElement, const std::array<double, NumNodes> rPotential);

template <int NumNodes>
void AssignPotentialsToWakeElement(Element& rElement, const array_1d<double, NumNodes>& rDistances, const std::array<double, 2*NumNodes>& rPotential);

template <int NumNodes>
BoundedVector<double,NumNodes> AssignDistancesToElement();

template <int NumNodes>
BoundedVector<double,NumNodes> GetDistanceValues();

void ComputeElementalSensitivitiesMatrixRow(ModelPart& rModelPart, double delta, unsigned int row, Matrix& rLHS_original, Vector& rRHS_original, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical);

template <int NumNodes>
void ComputeElementalSensitivities(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, NumNodes> rPotential);

template <int NumNodes>
void ComputeWakeElementalSensitivities(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 2*NumNodes> rPotential);

} // namespace PotentialFlowTestUtilities
} // namespace Kratos