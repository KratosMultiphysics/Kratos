//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/parallel_distance_calculator.h"

namespace Kratos
{

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::CalculateDistances(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar,
    const unsigned int MaxLevels,
    const double MaxDistance,
    Flags Options)
{
    KRATOS_TRY

    Check(rModelPart, rDistanceVar, rAreaVar);

    ResetVariables(rModelPart,rDistanceVar, MaxDistance);

    CalculateExactDistancesOnDividedElements(rModelPart, rDistanceVar, rAreaVar, MaxDistance, Options);

    ExtendDistancesByLayer(rModelPart, rDistanceVar, rAreaVar, MaxLevels, MaxDistance);

    AssignDistanceSign(rModelPart, rDistanceVar, rAreaVar, MaxDistance);

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::CalculateInterfacePreservingDistances(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar,
    const unsigned int MaxLevels,
    const double MaxDistance)
{
    KRATOS_TRY

    Check(rModelPart, rDistanceVar, rAreaVar);

    ResetVariables(rModelPart,rDistanceVar, MaxDistance);

    AbsDistancesOnDividedElements(rModelPart, rDistanceVar, rAreaVar, MaxDistance);

    ExtendDistancesByLayer(rModelPart, rDistanceVar, rAreaVar, MaxLevels, MaxDistance);

    AssignDistanceSign(rModelPart, rDistanceVar, rAreaVar, MaxDistance);

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::CalculateDistancesLagrangianSurface(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar,
    const unsigned int MaxLevels,
    const double MaxDistance)
{
    KRATOS_TRY

    //check that variables needed are in the model part
    const bool is_distributed = rModelPart.IsDistributed();
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rDistanceVar)) << "Distance variable is not in the model part" << std::endl;
    KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(rAreaVar)) << "Area Variable is not in the model part" << std::endl;
    if(is_distributed)
        KRATOS_ERROR_IF_NOT(rModelPart.HasNodalSolutionStepVariable(PARTITION_INDEX)) << "PARTITION_INDEX variable is not in the model part" << std::endl;

    // set to zero the distance
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        double& r_area = rNode.FastGetSolutionStepValue(rAreaVar);
        double& r_distance = rNode.FastGetSolutionStepValue(rDistanceVar);
        rNode.GetValue(rDistanceVar) = r_distance;
        if (rNode.IsNot(VISITED)) {
            r_area = 0.0;
            r_distance = 0.0;
        } else {
            r_area = 1.0;
        }
    });

    // Set the TLS container
    array_1d<double,TDim+1> visited, N;
    BoundedMatrix <double, TDim+1,TDim> DN_DX;
    typedef std::tuple<array_1d<double,TDim+1>, array_1d<double,TDim+1>, BoundedMatrix<double, TDim+1, TDim>> TLSType;
    TLSType tls_container = std::make_tuple(visited, N, DN_DX);

    // Extend the distances layer by layer up to a maximum level of layers
    for(unsigned int level=0; level<MaxLevels; level++)
    {
        //loop on active elements and advance the distance computation
        block_for_each(rModelPart.Elements(), tls_container, [&](Element& rElement, TLSType& rTLSContainer){
            auto& r_geom = rElement.GetGeometry();
            auto& r_visited = std::get<0>(rTLSContainer);
            auto& r_N = std::get<1>(rTLSContainer);
            auto& r_DN_DX = std::get<2>(rTLSContainer);
            for (unsigned int j=0; j<TDim+1; j++) {
                r_visited[j] = r_geom[j].Is(VISITED);
            }
            if (IsActive(r_visited)) {
                double volume;
                GeometryUtils::CalculateGeometryData(r_geom, r_DN_DX, r_N, volume);
                AddDistanceToNodes(rDistanceVar, rAreaVar, r_geom, r_DN_DX, volume);
            }
        });

        //mpi sync variables
        if(is_distributed)
        {
            block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
                if (rNode.Is(VISITED)) {
                    double& r_distance = rNode.FastGetSolutionStepValue(rDistanceVar);
                    rNode.GetValue(rDistanceVar) = r_distance;
                    r_distance = 0.0;
                } else {
                    rNode.GetValue(rDistanceVar) = 0.0;
                }
            });

            rModelPart.GetCommunicator().AssembleCurrentData(rAreaVar);
            rModelPart.GetCommunicator().AssembleCurrentData(rDistanceVar);

            block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
                rNode.FastGetSolutionStepValue(rDistanceVar) += rNode.GetValue(rDistanceVar);
            });

            rModelPart.GetCommunicator().GetDataCommunicator().Barrier();
        }

        //finalize the computation of the distance
        block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
            double& r_area = rNode.FastGetSolutionStepValue(rAreaVar);
            if(r_area > 1e-20 && rNode.IsNot(VISITED)) //this implies that node was computed at the current level and not before
            {
                double& r_distance = rNode.FastGetSolutionStepValue(rDistanceVar);
                r_distance /= r_area;
                rNode.Set(VISITED, true);
            }
        });
    }

    //assign the sign to the distance function according to the original distribution. Set to max for nodes that were not calculated
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        const double area = rNode.FastGetSolutionStepValue(rAreaVar);
        double& r_dist = rNode.FastGetSolutionStepValue(rDistanceVar);
        if(r_dist > MaxDistance || area < 1e-20) {
            r_dist = MaxDistance;
        }

        // if (rNode.Is(FLUID)) {
        //     r_dist = -std::abs(r_dist);
        // } else {
        //     r_dist = std::abs(r_dist);
        // }
    });

    KRATOS_CATCH("")
}

template<unsigned int TDim>
double ParallelDistanceCalculator<TDim>::FindMaximumEdgeSize(ModelPart& rModelPart)
{
    KRATOS_TRY

    double h_max = 0.0;

    //TODO: This can be done with the parallel utils reductions
    for(ModelPart::ElementsContainerType::iterator it=rModelPart.ElementsBegin(); it!=rModelPart.ElementsEnd(); it++)
    {
        Geometry<NodeType >&geom = it->GetGeometry();

        double h = 0.0;

        for(unsigned int i=0; i<TDim+1; i++)
        {

            double xc = geom[i].X();
            double yc = geom[i].Y();
            double zc = geom[i].Z();
            for(unsigned int j=i+1; j<TDim+1; j++)
            {
                double x = geom[j].X();
                double y = geom[j].Y();
                double z = geom[j].Z();
                double l = (x - xc)*(x - xc);
                l += (y - yc)*(y - yc);
                l += (z - zc)*(z - zc);

                if (l > h) h = l;
            }
        }

        h = sqrt(h);

        if(h > h_max) h_max = h;

    }

    h_max = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(h_max);

    return h_max;

    KRATOS_CATCH("");
}

template<unsigned int TDim>
bool ParallelDistanceCalculator<TDim>::IsDivided(const array_1d<double,TDim+1>& rDistance)
{
    unsigned int positive = 0;
    unsigned int negative = 0;

    for(unsigned int i=0; i<TDim+1; i++)
    {
        if(rDistance[i] >= 0)
            positive++;
        else
            negative++;
    }

    bool is_divided = false;
    if(positive > 0 && negative>0)
        is_divided = true;

    return is_divided;
}

template<unsigned int TDim>
bool ParallelDistanceCalculator<TDim>::IsActive(const array_1d<double,TDim+1>& rVisited)
{
    unsigned int positive = 0;

    for(unsigned int i=0; i<TDim+1; i++)
        if(rVisited[i] > 0.9999999999) //node was considered
            positive++;

    bool is_active = false;
    if(positive == TDim)
        is_active = true;

    return is_active;
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::AddDistanceToNodes(
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar,
    Geometry<NodeType>& rGeometry,
    const BoundedMatrix<double,TDim+1,TDim>& rDN_DX,
    const double& Volume)
{
    unsigned int unknown_node_index = 0;
    array_1d<double,TDim> d;
    double nodal_vol = Volume/static_cast<double>(TDim+1);
    double avg_dist = 0.0;

    //compute discriminant and find the index of the unknown node
    noalias(d) = ZeroVector(TDim);
    for (unsigned int iii = 0; iii < TDim + 1; iii++)
    {
        if (rGeometry[iii].Is(VISITED)) //identyfing the unknown node
        {
            const double distance = rGeometry[iii].FastGetSolutionStepValue(rDistanceVar);
            avg_dist += distance;
            for (unsigned int jjj = 0; jjj < TDim; jjj++)
                d[jjj] += rDN_DX(iii, jjj) * distance;
        }
        else
            unknown_node_index = iii;
    }
    avg_dist /= static_cast<double>(TDim);

    //finalizing computation of discriminant
    double c = -1.0;
    double a = 0.0;
    double b = 0.0;
    for (unsigned int jjj = 0; jjj < TDim; jjj++)
    {
        a += rDN_DX(unknown_node_index, jjj) * rDN_DX(unknown_node_index, jjj);
        b += d[jjj] * rDN_DX(unknown_node_index, jjj);
        c += d[jjj] * d[jjj];
    }
    b *= 2.0;

    //here we require (a*x^2 + b*x + c)^2 to be minimum (x represents the unknown distance)
    //this implies setting to zero
    //(a*x^2 + b*x + c)*(2ax+b) = 0
    double distance;

    double discriminant = b * b - 4.0 * a*c;

    if (discriminant < 0.0) //here we solve (2ax+b) = 0
    {
//                  double numerator = 0.0;
//                  double denominator = 0.0;
//                  for(unsigned int i=0; i<TDim+1; i++)
//                  {
//                      for (unsigned int jjj = 0; jjj < TDim; jjj++)
//                      {
//                          if(i != unknown_node_index)
//                            numerator += rDN_DX(unknown_node_index, jjj) * rDN_DX(i, jjj);
//                          else
//                            denominator += rDN_DX(unknown_node_index, jjj)*rDN_DX(unknown_node_index, jjj);
//                      }
//                  }
//                  distance = - numerator/denominator;
//
//                  KRATOS_WATCH(rGeometry[unknown_node_index].Id());


// 		KRATOS_WATCH(discriminant);
        distance = -b / (2.0*a); //avg_dist ; //
    }
    else //in this case we solve (a*x^2 + b*x + c)=0
    {
        //(accurate) computation of the distance
        //requires the solution of a*x^2+b*x+c=0
        double q, root1, root2;
        double sqrt_det = sqrt(discriminant);
        if (a != 0.0)
        {
            if (b > 0) q = -0.5 * (b + sqrt_det);
            else q = -0.5 * (b - sqrt_det);
            root1 = q / a;
            distance = root1;
            if(std::abs(q) > 0.0) {
                root2 = c / q;
                if (root2 > root1) distance = root2;
            }
        }
        else   //in this case we have a linear equation
        {
            distance = -c / b;
        }
    }

    if(distance < 0.0)
        distance = 1e-15;

    rGeometry[unknown_node_index].SetLock();
    rGeometry[unknown_node_index].FastGetSolutionStepValue(rDistanceVar) += distance*nodal_vol;
    rGeometry[unknown_node_index].FastGetSolutionStepValue(rAreaVar) += nodal_vol;
    rGeometry[unknown_node_index].UnSetLock();
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::Check(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar)
{
    KRATOS_TRY

    //check that variables needed are in the model part
    const int node_size = rModelPart.Nodes().size();
    KRATOS_ERROR_IF(node_size && !(rModelPart.NodesBegin()->SolutionStepsDataHas(rDistanceVar))) << "Distance variable is not in the model part." << std::endl;
    KRATOS_ERROR_IF(node_size && !(rModelPart.NodesBegin()->SolutionStepsDataHas(rAreaVar))) << "Area variable is not in the model part." << std::endl;
    if (rModelPart.IsDistributed()) {
        KRATOS_ERROR_IF(node_size && !(rModelPart.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX))) << "Area variable is not in the model part." << std::endl;
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::ResetVariables(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const double MaxDistance)
{
    KRATOS_TRY

    //reset the variables needed
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode) {
        //it->FastGetSolutionStepValue(rAreaVar) = 0.0;
        double& r_dist = rNode.FastGetSolutionStepValue(rDistanceVar);
        rNode.SetValue(rDistanceVar, r_dist); //here we copy the distance function to the fixed database

        rNode.Set(FLUID, r_dist < 0.0);

        r_dist = MaxDistance;

        rNode.Set(VISITED,false);
    });

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::CalculateExactDistancesOnDividedElements(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar,
    const double MaxDistance,
    Flags Options)
{
    KRATOS_TRY

    //identify the list of elements divided by the original distance distribution and recompute an "exact" distance
    //attempting to mantain the original position of the free surface
    //note that the backup value is used in calculating the position of the free surface and the divided elements
    array_1d<double,TDim+1> dist;
    block_for_each(rModelPart.Elements(), dist, [&](Element& rElement, array_1d<double,TDim+1>& rDist){
        auto& element_geometry = rElement.GetGeometry();
        // Set distances vector from non-historical database
        for (unsigned int j = 0; j < TDim + 1; j++) {
            rDist[j] = element_geometry[j].GetValue(rDistanceVar);
        }

        bool is_divided = IsDivided(rDist);
        if (is_divided) {
            if (Options.Is(CALCULATE_EXACT_DISTANCES_TO_PLANE)) {
                GeometryUtils::CalculateExactDistancesToPlane(element_geometry, rDist);
            } else {
                GeometryUtils::CalculateTetrahedraDistances(element_geometry, rDist);
            }

            // loop over nodes and apply the new distances.
            for (unsigned int i_node = 0; i_node < TDim+1; i_node++) {
                double& r_distance = element_geometry[i_node].GetSolutionStepValue(rDistanceVar);
                const double new_distance = rDist[i_node];

                element_geometry[i_node].SetLock();
                if (std::abs(r_distance) > std::abs(new_distance)) {
                    r_distance = new_distance;
                }
                element_geometry[i_node].Set(VISITED, true);
                element_geometry[i_node].UnSetLock();
            }
        }
    });

    //mpi sync variables
    rModelPart.GetCommunicator().AssembleCurrentData(rAreaVar);
    rModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);
    rModelPart.GetCommunicator().SynchronizeCurrentDataToMin(rDistanceVar);

    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        if(rNode.IsNot(VISITED)) {
            rNode.FastGetSolutionStepValue(rAreaVar) = 0.0;
            rNode.FastGetSolutionStepValue(rDistanceVar) = 0.0;
        } else {
            rNode.GetSolutionStepValue(rAreaVar) = 1.00; // This is not correct
        }
    });

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::AbsDistancesOnDividedElements(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar,
    const double MaxDistance)
{
    KRATOS_TRY

    //identify the list of elements divided by the original distance distribution and recompute an "exact" distance
    //attempting to mantain the original position of the free surface
    //note that the backup value is used in calculating the position of the free surface and the divided elements
    array_1d<double, TDim+1> dist;
    block_for_each(rModelPart.Elements(), dist, [&](Element& rElement, array_1d<double,TDim+1>& rDist){
        // Set distances vector from non-historical database
        auto& element_geometry = rElement.GetGeometry();
        for (unsigned int j = 0; j < TDim + 1; j++) {
            rDist[j] = element_geometry[j].GetValue(rDistanceVar);
        }

        // Check intersection
        bool is_divided = IsDivided(dist);
        if (is_divided) {
            // loop over nodes and apply the new distances.
            for (unsigned int i_node = 0; i_node < TDim + 1; i_node++) {
                element_geometry[i_node].SetLock();
                element_geometry[i_node].GetSolutionStepValue(rDistanceVar) = std::abs(rDist[i_node]);
                element_geometry[i_node].Set(VISITED, true);
                element_geometry[i_node].UnSetLock();
            }
        }
    });

    //mpi sync variables
    rModelPart.GetCommunicator().AssembleCurrentData(rAreaVar);
    rModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);
    rModelPart.GetCommunicator().SynchronizeCurrentDataToMin(rDistanceVar);

    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        if (rNode.IsNot(VISITED)) {
            rNode.FastGetSolutionStepValue(rAreaVar) = 0.0;
            rNode.FastGetSolutionStepValue(rDistanceVar) = 0.0;
        } else {
            rNode.FastGetSolutionStepValue(rAreaVar) = 1.0; // This is not correct
        }
    });

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::ExtendDistancesByLayer(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar,
    const unsigned int MaxLevels,
    const double MaxDistance)
{
    KRATOS_TRY

    // Set the TLS container
    array_1d<double,TDim+1> visited, N;
    BoundedMatrix <double, TDim+1,TDim> DN_DX;
    typedef std::tuple<array_1d<double,TDim+1>, array_1d<double,TDim+1>, BoundedMatrix<double, TDim+1, TDim>> TLSType;
    TLSType tls_container = std::make_tuple(visited, N, DN_DX);

    //now extend the distances layer by layer up to a maximum level of layers
    for(unsigned int level=0; level<MaxLevels; level++)
    {
        //loop on active elements and advance the distance computation
        block_for_each(rModelPart.Elements(), tls_container, [&](Element& rElement, TLSType& rTLSContainer){
            auto& r_geom = rElement.GetGeometry();
            auto& r_visited = std::get<0>(rTLSContainer);
            auto& r_N = std::get<1>(rTLSContainer);
            auto& r_DN_DX = std::get<2>(rTLSContainer);
            for (unsigned int j=0; j<TDim+1; j++) {
                r_visited[j] = r_geom[j].Is(VISITED);
            }
            if (IsActive(r_visited)) {
                double volume;
                GeometryUtils::CalculateGeometryData(r_geom, r_DN_DX, r_N, volume);
                AddDistanceToNodes(rDistanceVar, rAreaVar, r_geom, r_DN_DX, volume);
            }
        });

        //mpi sync variables
        if(rModelPart.IsDistributed())
        {
            block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
                if (rNode.Is(VISITED)) {
                    double& r_distance = rNode.FastGetSolutionStepValue(rDistanceVar);
                    rNode.GetValue(rDistanceVar) = r_distance;
                    r_distance = 0.0;
                } else {
                    rNode.GetValue(rDistanceVar) = 0.0;
                }
            });

            rModelPart.GetCommunicator().AssembleCurrentData(rAreaVar);
            rModelPart.GetCommunicator().AssembleCurrentData(rDistanceVar);

            block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
                rNode.FastGetSolutionStepValue(rDistanceVar) += rNode.GetValue(rDistanceVar);
            });

            rModelPart.GetCommunicator().GetDataCommunicator().Barrier();
        }

        //finalize the computation of the distance
        block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
            const double area = rNode.FastGetSolutionStepValue(rAreaVar);
            if (area > 1e-20 && rNode.IsNot(VISITED)) { //this implies that node was computed at the current level and not before
                double& r_distance = rNode.FastGetSolutionStepValue(rDistanceVar);
                r_distance /= area;
                rNode.Set(VISITED, true);
            }
        });
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculator<TDim>::AssignDistanceSign(
    ModelPart& rModelPart,
    const Variable<double>& rDistanceVar,
    const Variable<double>& rAreaVar,
    const double MaxDistance)
{
    KRATOS_TRY

    //assign the sign to the distance function according to the original distribution. Set to max for nodes that were not calculated
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        const double area = rNode.FastGetSolutionStepValue(rAreaVar);
        double& r_dist = rNode.FastGetSolutionStepValue(rDistanceVar);

        KRATOS_ERROR_IF(r_dist < 0.0) << "IMPOSSIBLE negative distance found !!" << std::endl;
        if (r_dist > MaxDistance || area < 1e-20) {
            r_dist = MaxDistance;
        }

        if(rNode.Is(FLUID)) {
            r_dist = -std::abs(r_dist);
        } else {
            r_dist = std::abs(r_dist);
        }
    });

    KRATOS_CATCH("")
}

}  // namespace Kratos.
