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


#include "tests/cpp_tests/test_utilities.h"
#include "includes/model_part.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos {
namespace PotentialFlowTestUtilities {

template <int NumNodes>
void AssignPotentialsToNormalElement(Element& rElement, const std::array<double, NumNodes> rPotential)
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
}

template <int NumNodes>
void AssignPotentialsToWakeElement(Element& rElement, const array_1d<double, NumNodes>& rDistances, const std::array<double, 2*NumNodes>& rPotential)
{
    for (unsigned int i = 0; i < NumNodes; i++){
        if (rDistances(i) > 0.0)
            rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i];
        else
            rElement.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = rPotential[i];
    }
    for (unsigned int i = 0; i < NumNodes; i++){
        if (rDistances(i) < 0.0)
            rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) = rPotential[i+NumNodes];
        else
            rElement.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) = rPotential[i+NumNodes];
    }
}

template <>
BoundedVector<double,3> GetDistanceValues<3>()
{
    BoundedVector<double,3> distances;
    distances(0) = 1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    return distances;
}

template <>
BoundedVector<double,4> GetDistanceValues<4>()
{
    BoundedVector<double,4> distances;
    distances(0) = -1.0;
    distances(1) = -1.0;
    distances(2) = -1.0;
    distances(3) = 1.0;
    return distances;
}

template <int NumNodes>
BoundedVector<double,NumNodes> AssignDistancesToElement()
{
    return GetDistanceValues<NumNodes>();
}

void ComputeElementalSensitivitiesMatrixRow(ModelPart& rModelPart, double delta, unsigned int row, Matrix& rLHS_original, Vector& rRHS_original, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical){
    Element::Pointer p_element = rModelPart.pGetElement(1);
    const unsigned int number_of_nodes = p_element->GetGeometry().size();

    // Compute pinged LHS and RHS
    Vector RHS_pinged = ZeroVector(number_of_nodes);
    Matrix LHS_pinged = ZeroMatrix(number_of_nodes, number_of_nodes);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS_pinged, RHS_pinged, r_current_process_info);

    for(unsigned int k = 0; k < rLHS_original.size2(); k++){
        // Compute the finite difference estimate of the sensitivity
        rLHS_finite_diference( k, row) = -(RHS_pinged(k)-rRHS_original(k)) / delta;
        // Compute the average of the original and pinged analytic sensitivities
        rLHS_analytical( k, row) = 0.5 * (rLHS_original(k,row) + LHS_pinged(k,row));
    }
}

template <int NumNodes>
void ComputeElementalSensitivities(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, NumNodes> rPotential){
    Element::Pointer p_element = rModelPart.pGetElement(1);

    AssignPotentialsToNormalElement<NumNodes>(*p_element, rPotential);

    // Compute original RHS and LHS
    Vector RHS_original = ZeroVector(NumNodes);
    Matrix LHS_original = ZeroMatrix(NumNodes, NumNodes);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < NumNodes; i++){
        // Pinging
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

        ComputeElementalSensitivitiesMatrixRow(rModelPart, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

        // Unpinging
        p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
    }
}

template <int NumNodes>
void ComputeWakeElementalSensitivities(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 2*NumNodes> rPotential){
    Element::Pointer p_element = rModelPart.pGetElement(1);

    BoundedVector<double, NumNodes> distances = AssignDistancesToElement<NumNodes>();
    p_element->GetValue(WAKE_ELEMENTAL_DISTANCES) = distances;
    p_element->GetValue(WAKE) = true;

    AssignPotentialsToWakeElement<NumNodes>(*p_element, distances, rPotential);

    // Compute original RHS and LHS
    Vector RHS_original = ZeroVector(2*NumNodes);
    Matrix LHS_original = ZeroMatrix(2*NumNodes, 2*NumNodes);
    const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
    p_element->CalculateLocalSystem(LHS_original, RHS_original, r_current_process_info);

    double delta = 1e-3;
    for(unsigned int i = 0; i < 2*NumNodes; i++){
        if(i < NumNodes){
            // Pinging
            if (distances(i) > 0.0)
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;
            else
                p_element->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;

            ComputeElementalSensitivitiesMatrixRow(rModelPart, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

            // Unpinging
            if (distances(i) > 0.0)
                p_element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
            else
                p_element->GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
        }
        else{
            // Pinging
            if (distances(i-NumNodes) > 0.0)
                p_element->GetGeometry()[i-NumNodes].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) += delta;
            else
                p_element->GetGeometry()[i-NumNodes].FastGetSolutionStepValue(VELOCITY_POTENTIAL) += delta;

            ComputeElementalSensitivitiesMatrixRow(rModelPart, delta, i, LHS_original, RHS_original, rLHS_finite_diference, rLHS_analytical);

            // Unpinging
            if (distances(i-NumNodes) > 0.0)
                p_element->GetGeometry()[i-NumNodes].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL) -= delta;
            else
                p_element->GetGeometry()[i-NumNodes].FastGetSolutionStepValue(VELOCITY_POTENTIAL) -= delta;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template instantiation

// 2D
template void AssignPotentialsToNormalElement<3>(Element& rElement, const std::array<double, 3> rPotential);
template void AssignPotentialsToWakeElement<3>(Element& rElement, const array_1d<double, 3>& rDistances, const std::array<double, 6>& rPotential);
template BoundedVector<double,3> AssignDistancesToElement<3>();
template void ComputeElementalSensitivities<3>(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical,  const std::array<double, 3> rPotential);
template void ComputeWakeElementalSensitivities<3>(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 6> rPotential);

// 3D
template void AssignPotentialsToNormalElement<4>(Element& rElement, const std::array<double, 4> rPotential);
template void AssignPotentialsToWakeElement<4>(Element& rElement, const array_1d<double, 4>& rDistances, const std::array<double, 8>& rPotential);
template BoundedVector<double,4> AssignDistancesToElement<4>();
template void ComputeElementalSensitivities<4>(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical,  const std::array<double, 4> rPotential);
template void ComputeWakeElementalSensitivities<4>(ModelPart& rModelPart, Matrix& rLHS_finite_diference, Matrix& rLHS_analytical, const std::array<double, 8> rPotential);


} // namespace PotentialFlowTestUtilities
} // namespace Kratos
