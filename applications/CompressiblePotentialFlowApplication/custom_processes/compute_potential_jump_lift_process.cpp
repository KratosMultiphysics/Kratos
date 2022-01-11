//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Marc Nunez
//


#include "compute_potential_jump_lift_process.h"
#include "compressible_potential_flow_application_variables.h"
// #include "includes/cfd_variables.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
// Constructor for ComputePotentialJumpLiftProcess Process
template <unsigned int Dim, unsigned int NumNodes>
ComputePotentialJumpLiftProcess<Dim, NumNodes>::ComputePotentialJumpLiftProcess(ModelPart& rModelPart,
                Vector& rResultantForce
                ):
        Process(),
        mrModelPart(rModelPart),
        mrResultantForce(rResultantForce)
    {
    }

template <unsigned int Dim, unsigned int NumNodes>
void ComputePotentialJumpLiftProcess<Dim, NumNodes>::Execute()
{
    KRATOS_TRY;

    mrResultantForce = ZeroVector(3);

    //Declaring auxilary variables needed to use with omp
    double potential_integral = 0.0;
    // double drag_integral = 0.0;
    double area = 0.0;
    double max_x = 0.0;

    for (auto& r_cond : mrModelPart.Conditions()) {
        double x_center = r_cond.GetGeometry().Center().X();
        max_x = std::max(x_center, max_x);
    }

    for (auto& r_cond : mrModelPart.Conditions()) {

        array_1d<double, 3> distances;
        bool not_zero = true;
        for (std::size_t i =0; i<r_cond.GetGeometry().size(); i++) {
            distances[i] = r_cond.GetGeometry()[i].GetValue(WAKE_DISTANCE);
            if (std::abs(distances[i]) < 1e-6) {
                not_zero = false;
            }
        }

        const bool is_wake = PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(distances);
        const bool is_max = std::abs( max_x- r_cond.GetGeometry().Center().X()) < 1e-3;

        if (is_wake && not_zero && is_max) {
            Geometry<Node<3>>::PointsArrayType points_array;
            points_array.push_back(Node<3>::Pointer(new Node<3>(1, r_cond.GetGeometry().Points()[0].Y(),r_cond.GetGeometry().Points()[0].Z(), 0.0)));
            points_array.push_back(Node<3>::Pointer(new Node<3>(2, r_cond.GetGeometry().Points()[1].Y(),r_cond.GetGeometry().Points()[1].Z(), 0.0)));
            points_array.push_back(Node<3>::Pointer(new Node<3>(3, r_cond.GetGeometry().Points()[2].Y(),r_cond.GetGeometry().Points()[2].Z(), 0.0)));

            Triangle2D3<Node<3>> triangle_geom(points_array);
            const auto p_geom = triangle_geom.Create(triangle_geom);
            ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(p_geom, Vector(distances));


            Matrix rInterfacePositiveSideShapeFunctionsValues;
            ModifiedShapeFunctions::ShapeFunctionsGradientsType rInterfacePositiveSideShapeFunctionsGradientsValues;
            Vector rInterfacePositiveSideWeightsValues;

            pModifiedShFunc ->ComputeInterfacePositiveSideShapeFunctionsAndGradientsValues(rInterfacePositiveSideShapeFunctionsValues,
                                                                        rInterfacePositiveSideShapeFunctionsGradientsValues,
                                                                        rInterfacePositiveSideWeightsValues,
                                                                        GeometryData::GI_GAUSS_1);
            area += rInterfacePositiveSideWeightsValues(0);

            double potential_jump = 0.0;
            for (std::size_t i =0; i<r_cond.GetGeometry().size(); i++) {
                double potential = r_cond.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
                double aux_potential = r_cond.GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
                potential_jump += std::abs(potential-aux_potential)*rInterfacePositiveSideShapeFunctionsValues(0, i);
            }

            potential_integral += potential_jump*rInterfacePositiveSideWeightsValues(0);

        }

    }

    // Storing final result
    mrResultantForce[0] = potential_integral;
    mrResultantForce[1] = 0.0;
    mrResultantForce[2] = area;

    KRATOS_CATCH("");
}

template<>
ModifiedShapeFunctions::Pointer ComputePotentialJumpLiftProcess<2, 3>::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
        return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template<>
ModifiedShapeFunctions::Pointer ComputePotentialJumpLiftProcess<3, 4>::pGetModifiedShapeFunctions(const GeomPointerType pGeometry, const Vector& rDistances) const {
        return Kratos::make_unique<Triangle2D3ModifiedShapeFunctions>(pGeometry, rDistances);
}

template class ComputePotentialJumpLiftProcess<2, 3>;
template class ComputePotentialJumpLiftProcess<3, 4>;
}// Namespace Kratos
