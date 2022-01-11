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
#include "containers/model.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"
#include "processes/parallel_distance_calculation_process.h"

namespace Kratos
{

template<unsigned int TDim>
ParallelDistanceCalculationProcess<TDim>::ParallelDistanceCalculationProcess(
    Model& rModel,
    Parameters Settings)
    : ParallelDistanceCalculationProcess<TDim>(rModel.GetModelPart(Settings["model_part_name"].GetString()),Settings)
{}

template<unsigned int TDim>
ParallelDistanceCalculationProcess<TDim>::ParallelDistanceCalculationProcess(
    ModelPart& rModelPart,
    Parameters Settings)
    : mrModelPart(rModelPart)
{
    Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());
    mpDistanceVar = &KratosComponents< Variable<double> >::Get(Settings["distance_variable"].GetString());
    mpAreaVar = &KratosComponents< Variable<double> >::Get(Settings["nodal_area_variable"].GetString());

    mMaxLevels = Settings["max_levels"].GetInt();
    mMaxDistance = Settings["max_distance"].GetDouble();
    mCalculateExactDistancesToPlane = Settings["calculate_exact_distances_to_plane"].GetBool();
}

template<unsigned int TDim>
const Parameters ParallelDistanceCalculationProcess<TDim>::GetDefaultParameters() const
{
    return Parameters(R"({
        "model_part_name" : "PLEASE_SPECIFY_MODEL_PART",
        "distance_variable" : "DISTANCE",
        "nodal_area_variable" : "NODAL_AREA",
        "max_levels" : 25,
        "max_distance" : 1.0,
        "calculate_exact_distances_to_plane" : false
    })");
}

template<unsigned int TDim>
double ParallelDistanceCalculationProcess<TDim>::FindMaximumEdgeSize()
{
    KRATOS_TRY

    double h_max = 0.0;

    for(ModelPart::ElementsContainerType::iterator it=mrModelPart.ElementsBegin(); it!=mrModelPart.ElementsEnd(); it++)
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

    h_max = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(h_max);

    return h_max;

    KRATOS_CATCH("");
}

template<unsigned int TDim>
bool ParallelDistanceCalculationProcess<TDim>::IsDivided(const array_1d<double,TDim+1>& rDistance)
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
bool ParallelDistanceCalculationProcess<TDim>::IsActive(const array_1d<double,TDim+1>& rVisited)
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
void ParallelDistanceCalculationProcess<TDim>::AddDistanceToNodes(
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
            const double distance = rGeometry[iii].FastGetSolutionStepValue(*mpDistanceVar);
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
    rGeometry[unknown_node_index].FastGetSolutionStepValue(*mpDistanceVar) += distance*nodal_vol;
    rGeometry[unknown_node_index].FastGetSolutionStepValue(*mpAreaVar) += nodal_vol;
    rGeometry[unknown_node_index].UnSetLock();
}

template<unsigned int TDim>
int ParallelDistanceCalculationProcess<TDim>::Check()
{
    KRATOS_TRY

    //check that variables needed are in the model part
    const int node_size = mrModelPart.Nodes().size();
    KRATOS_ERROR_IF(node_size && !(mrModelPart.NodesBegin()->SolutionStepsDataHas(*mpDistanceVar))) << "Distance variable is not in the model part." << std::endl;
    KRATOS_ERROR_IF(node_size && !(mrModelPart.NodesBegin()->SolutionStepsDataHas(*mpAreaVar))) << "Area variable is not in the model part." << std::endl;
    if (mrModelPart.IsDistributed()) {
        KRATOS_ERROR_IF(node_size && !(mrModelPart.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX))) << "Area variable is not in the model part." << std::endl;
    }

    return 0;

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculationProcess<TDim>::Execute()
{
    ResetVariables();

    CalculateExactDistancesOnDividedElements();

    ExtendDistancesByLayer();

    AssignDistanceSign();
}

template<unsigned int TDim>
void ParallelDistanceCalculationProcess<TDim>::ResetVariables()
{
    KRATOS_TRY

    //reset the variables needed
    block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode) {
        //it->FastGetSolutionStepValue(*mpAreaVar) = 0.0;
        double& r_dist = rNode.FastGetSolutionStepValue(*mpDistanceVar);
        rNode.SetValue(*mpDistanceVar, r_dist); //here we copy the distance function to the fixed database

        rNode.Set(FLUID, r_dist < 0.0);

        r_dist = mMaxDistance;

        rNode.Set(VISITED,false);
    });

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculationProcess<TDim>::CalculateExactDistancesOnDividedElements()
{
    KRATOS_TRY

    //identify the list of elements divided by the original distance distribution and recompute an "exact" distance
    //attempting to mantain the original position of the free surface
    //note that the backup value is used in calculating the position of the free surface and the divided elements
    array_1d<double,TDim+1> dist;
    block_for_each(mrModelPart.Elements(), dist, [&](Element& rElement, array_1d<double,TDim+1>& rDist){
        auto& element_geometry = rElement.GetGeometry();
        // Set distances vector from non-historical database
        for (unsigned int j = 0; j < TDim + 1; j++) {
            rDist[j] = element_geometry[j].GetValue(*mpDistanceVar);
        }

        bool is_divided = IsDivided(rDist);
        if (is_divided) {
            if (mCalculateExactDistancesToPlane) {
                GeometryUtils::CalculateExactDistancesToPlane(element_geometry, rDist);
            } else {
                GeometryUtils::CalculateTetrahedraDistances(element_geometry, rDist);
            }

            // loop over nodes and apply the new distances.
            for (unsigned int i_node = 0; i_node < TDim+1; i_node++) {
                double& r_distance = element_geometry[i_node].GetSolutionStepValue(*mpDistanceVar);
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
    mrModelPart.GetCommunicator().AssembleCurrentData(*mpAreaVar);
    mrModelPart.GetCommunicator().SynchronizeOrNodalFlags(VISITED);
    mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(*mpDistanceVar);

    block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
        if(rNode.IsNot(VISITED)) {
            rNode.FastGetSolutionStepValue(*mpAreaVar) = 0.0;
            rNode.FastGetSolutionStepValue(*mpDistanceVar) = 0.0;
        } else {
            rNode.GetSolutionStepValue(*mpAreaVar) = 1.00; // This is not correct
        }
    });

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculationProcess<TDim>::ExtendDistancesByLayer()
{
    KRATOS_TRY

    // Set the TLS container
    array_1d<double,TDim+1> visited, N;
    BoundedMatrix <double, TDim+1,TDim> DN_DX;
    typedef std::tuple<array_1d<double,TDim+1>, array_1d<double,TDim+1>, BoundedMatrix<double, TDim+1, TDim>> TLSType;
    TLSType tls_container = std::make_tuple(visited, N, DN_DX);

    //now extend the distances layer by layer up to a maximum level of layers
    for(unsigned int level = 0; level < mMaxLevels; level++)
    {
        //loop on active elements and advance the distance computation
        block_for_each(mrModelPart.Elements(), tls_container, [&](Element& rElement, TLSType& rTLSContainer){
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
                AddDistanceToNodes(r_geom, r_DN_DX, volume);
            }
        });

        //mpi sync variables
        if(mrModelPart.IsDistributed())
        {
            block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
                if (rNode.Is(VISITED)) {
                    double& r_distance = rNode.FastGetSolutionStepValue(*mpDistanceVar);
                    rNode.GetValue(*mpDistanceVar) = r_distance;
                    r_distance = 0.0;
                } else {
                    rNode.GetValue(*mpDistanceVar) = 0.0;
                }
            });

            mrModelPart.GetCommunicator().AssembleCurrentData(*mpAreaVar);
            mrModelPart.GetCommunicator().AssembleCurrentData(*mpDistanceVar);

            block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
                rNode.FastGetSolutionStepValue(*mpDistanceVar) += rNode.GetValue(*mpDistanceVar);
            });

            mrModelPart.GetCommunicator().GetDataCommunicator().Barrier();
        }

        //finalize the computation of the distance
        block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
            const double area = rNode.FastGetSolutionStepValue(*mpAreaVar);
            if (area > 1e-20 && rNode.IsNot(VISITED)) { //this implies that node was computed at the current level and not before
                double& r_distance = rNode.FastGetSolutionStepValue(*mpDistanceVar);
                r_distance /= area;
                rNode.Set(VISITED, true);
            }
        });
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim>
void ParallelDistanceCalculationProcess<TDim>::AssignDistanceSign()
{
    KRATOS_TRY

    //assign the sign to the distance function according to the original distribution. Set to max for nodes that were not calculated
    block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
        const double area = rNode.FastGetSolutionStepValue(*mpAreaVar);
        double& r_dist = rNode.FastGetSolutionStepValue(*mpDistanceVar);

        KRATOS_ERROR_IF(r_dist < 0.0) << "IMPOSSIBLE negative distance found !!" << std::endl;
        if (r_dist > mMaxDistance || area < 1e-20) {
            r_dist = mMaxDistance;
        }

        if(rNode.Is(FLUID)) {
            r_dist = -std::abs(r_dist);
        } else {
            r_dist = std::abs(r_dist);
        }
    });

    KRATOS_CATCH("")
}

template class ParallelDistanceCalculationProcess<2>;
template class ParallelDistanceCalculationProcess<3>;

}  // namespace Kratos.
