//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED )
#define  KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
template< unsigned int TDim>
class ParallelDistanceCalculator
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    KRATOS_DEFINE_LOCAL_FLAG(CALCULATE_EXACT_DISTANCES_TO_PLANE);

    /// Pointer definition of ParallelDistanceCalculator
    KRATOS_CLASS_POINTER_DEFINITION(ParallelDistanceCalculator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParallelDistanceCalculator() {};

    /// Destructor.
    virtual ~ParallelDistanceCalculator() {};

    ///Function to calculate a signed distance function suitable for calculations using the Level Set Method
    ///the function assumes given a "signed distance" distributions and recomputes the distances
    ///respecting as accurately as possible the position of the zero of the original distributions
    ///@param rModelPart is the ModelPart on which we will operate
    ///@param rDistanceVar is the Variable that we will use in calculating the distance
    ///@param rAreaVar is the Variable that we will use for L2 projections
    ///@param max_levels is the number of maximum "layers" of element that will be used in the calculation of the distances
    ///@param max_distance distances will not be computed after reaching this limit
    void CalculateDistances(ModelPart& rModelPart,
                            const Variable<double>& rDistanceVar,
                            const Variable<double>& rAreaVar,
                            const unsigned int max_levels,
                            const double max_distance,
                            Flags Options = CALCULATE_EXACT_DISTANCES_TO_PLANE.AsFalse())
    {
        KRATOS_TRY

		Check(rModelPart, rDistanceVar, rAreaVar);

		ResetVariables(rModelPart,rDistanceVar, max_distance);

		CalculateExactDistancesOnDividedElements(rModelPart, rDistanceVar, rAreaVar, max_distance, Options);

        ExtendDistancesByLayer(rModelPart, rDistanceVar, rAreaVar, max_levels, max_distance);

		AssignDistanceSign(rModelPart, rDistanceVar, rAreaVar, max_distance);

		KRATOS_CATCH("")
    }


    ///Function to calculate a signed distance function suitable for calculations using the Level Set Method
	///The difference of this function with previous one is the fact that it wont recalculate the exact distance
	///in divided elements in order to preserve the current distance.
    ///the function assumes given a "signed distance" distributions and recomputes the distances
    ///respecting as accurately as possible the position of the zero of the original distributions
    ///@param rModelPart is the ModelPart on which we will operate
    ///@param rDistanceVar is the Variable that we will use in calculating the distance
    ///@param rAreaVar is the Variable that we will use for L2 projections
    ///@param max_levels is the number of maximum "layers" of element that will be used in the calculation of the distances
    ///@param max_distance distances will not be computed after reaching this limit
    void CalculateInterfacePreservingDistances(ModelPart& rModelPart,
                            const Variable<double>& rDistanceVar,
                            const Variable<double>& rAreaVar,
                            const unsigned int max_levels,
                            const double max_distance)
    {
        KRATOS_TRY

		Check(rModelPart, rDistanceVar, rAreaVar);

		ResetVariables(rModelPart,rDistanceVar, max_distance);

		AbsDistancesOnDividedElements(rModelPart, rDistanceVar, rAreaVar, max_distance);

        ExtendDistancesByLayer(rModelPart, rDistanceVar, rAreaVar, max_levels, max_distance);

		AssignDistanceSign(rModelPart, rDistanceVar, rAreaVar, max_distance);

		KRATOS_CATCH("")
    }


    /// A simplified version of CalculateDistances to be used when the rDistanceVar == 0 surface is described by a set of nodes
    /**
     * @param rModelPart is the ModelPart on which we will operate
     * @param rDistanceVar is the Variable that we will use in calculating the distance
     * @param rAreaVar is the Variable that we will use for L2 projections
     * @param max_levels is the number of maximum "layers" of element that will be used in the calculation of the distances
     * @param max_distance distances will not be computed after reaching this limit
     * @see ParallelDistanceCalculator::CalculateDistances
     */
    void CalculateDistancesLagrangianSurface(
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

    //**********************************************************************************
    //**********************************************************************************
    double FindMaximumEdgeSize(ModelPart& r_model_part)
    {
        KRATOS_TRY

        double h_max = 0.0;

        for(ModelPart::ElementsContainerType::iterator it=r_model_part.ElementsBegin(); it!=r_model_part.ElementsEnd(); it++)
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

        h_max = r_model_part.GetCommunicator().GetDataCommunicator().MaxAll(h_max);

        return h_max;

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "ParallelDistanceCalculator" << TDim << "D";
        return buffer.str();
    };

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ParallelDistanceCalculator" << TDim << "D";
    };

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const {};


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{

    //*******************************************************************
    bool IsDivided(array_1d<double,TDim+1>& dist)
    {
        unsigned int positive = 0;
        unsigned int negative = 0;

        for(unsigned int i=0; i<TDim+1; i++)
        {
            if(dist[i] >= 0)
                positive++;
            else
                negative++;
        }

        bool is_divided = false;
        if(positive > 0 && negative>0)
            is_divided = true;

        return is_divided;
    }

    //*******************************************************************
    bool IsActive(array_1d<double,TDim+1>& visited)
    {
        unsigned int positive = 0;

        for(unsigned int i=0; i<TDim+1; i++)
            if(visited[i] > 0.9999999999) //node was considered
                positive++;

        bool is_active = false;
        if(positive == TDim)
            is_active = true;

        return is_active;
    }

    //*******************************************************************
    void ComputeExactDistances(const BoundedMatrix <double, TDim+1,TDim>& DN_DX,
                               const double& Area,
                               Geometry<NodeType >& geom,
                               const array_1d<double,TDim+1>& distances,
                               array_1d<double,TDim+1>& exact_dist
                              )
    {
        array_1d<double,TDim> grad_d;
        array_1d<double,3> coord_on_0 = ZeroVector(3);
        array_1d<double,3> temp;

        //compute the gradient of the distance and normalize it
        noalias(grad_d) = prod(trans(DN_DX),distances);
        double norm = norm_2(grad_d);
        grad_d /= norm;

        //find one division point on one edge
        for(unsigned int i = 1; i<TDim+1; i++)
        {
            if(distances[0]*distances[i]<=0.0) //if the edge is divided
            {
                double delta_d = fabs(distances[i]) + fabs(distances[0]);

                if(delta_d>1e-20)
                {
                    double Ni = fabs(distances[0]) / delta_d;
                    double N0 = fabs(distances[i]) / delta_d;

                    noalias(coord_on_0) = N0 * geom[0].Coordinates();
                    noalias(coord_on_0) += Ni * geom[i].Coordinates();
                }
                else
                    noalias(coord_on_0) = geom[0].Coordinates();

                break;

            }
        }


        //now calculate the distance of all the nodes from the elemental free surface
        for(unsigned int i = 0; i<TDim+1; i++)
        {
            noalias(temp) = geom[i].Coordinates();
            noalias(temp) -= coord_on_0 ;

            double real_distance = 0.0;
            for(unsigned int k=0; k<TDim; k++)
                real_distance += temp[k]*grad_d[k];
            real_distance = fabs(real_distance);

            exact_dist[i] = real_distance;
        }
    }

    //*******************************************************************
    void AddDistanceToNodesNew(const Variable<double>& rDistanceVar,
                            const Variable<double>& rAreaVar,
                            Geometry<NodeType >& geom,
                            const BoundedMatrix <double, TDim+1,TDim>& DN_DX,
                            const double& Volume
                           )
    {
        unsigned int unknown_node_index = 0;
        array_1d<double,TDim> d;
        double nodal_vol = Volume/static_cast<double>(TDim+1);
        double avg_dist = 0.0;


		Matrix coord_a(3,3);
		int row = 0;
		int reference_node_index;

        //compute discriminant and find the index of the unknown node
        noalias(d) = ZeroVector(TDim);
        for (unsigned int iii = 0; iii < TDim + 1; iii++)
        {
            if (geom[iii].Is(VISITED)) //identyfing the known node
            {
				reference_node_index = iii;
				for(int i_coord = 0 ; i_coord < 3 ; i_coord++)
					coord_a(row,i_coord) = geom[iii].Coordinates()[i_coord];


                d[row] = geom[iii].FastGetSolutionStepValue(rDistanceVar);
                avg_dist += d[row];
				row++;
            }
            else
                unknown_node_index = iii;
        }
        avg_dist /= static_cast<double>(TDim);

		Matrix inverse_a(3,3);
		double det_a;
		MathUtils<double>::InvertMatrix3(coord_a,inverse_a,det_a);
		array_1d<double,TDim> x;  // normal to the surface
		noalias(x) = prod(inverse_a,d);
		double norm_x = norm_2(x);
		x /= norm_x;
		array_1d<double,TDim> v = geom[unknown_node_index].Coordinates() - geom[reference_node_index].Coordinates();

		double distance = inner_prod(x,v);
		distance += geom[reference_node_index].FastGetSolutionStepValue(rDistanceVar);
		//KRATOS_WATCH(coord_a)
		//KRATOS_WATCH(distance)

        geom[unknown_node_index].SetLock();
        geom[unknown_node_index].FastGetSolutionStepValue(rDistanceVar) += distance*nodal_vol;
        geom[unknown_node_index].FastGetSolutionStepValue(rAreaVar) += nodal_vol;
        geom[unknown_node_index].UnSetLock();

        //GeometryUtils::CalculateTetrahedraDistances(element_geometry, dist);

     }


    //*******************************************************************
    void AddDistanceToNodes(const Variable<double>& rDistanceVar,
                            const Variable<double>& rAreaVar,
                            Geometry<NodeType >& geom,
                            const BoundedMatrix <double, TDim+1,TDim>& DN_DX,
                            const double& Volume
                           )
    {
        unsigned int unknown_node_index = 0;
        array_1d<double,TDim> d;
        double nodal_vol = Volume/static_cast<double>(TDim+1);
        double avg_dist = 0.0;

        //compute discriminant and find the index of the unknown node
        noalias(d) = ZeroVector(TDim);
        for (unsigned int iii = 0; iii < TDim + 1; iii++)
        {
            if (geom[iii].Is(VISITED)) //identyfing the unknown node
            {
                const double distance = geom[iii].FastGetSolutionStepValue(rDistanceVar);
                avg_dist += distance;
                for (unsigned int jjj = 0; jjj < TDim; jjj++)
                    d[jjj] += DN_DX(iii, jjj) * distance;
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
            a += DN_DX(unknown_node_index, jjj) * DN_DX(unknown_node_index, jjj);
            b += d[jjj] * DN_DX(unknown_node_index, jjj);
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
//                            numerator += DN_DX(unknown_node_index, jjj) * DN_DX(i, jjj);
//                          else
//                            denominator += DN_DX(unknown_node_index, jjj)*DN_DX(unknown_node_index, jjj);
//                      }
//                  }
//                  distance = - numerator/denominator;
//
//                  KRATOS_WATCH(geom[unknown_node_index].Id());


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
                root2 = c / q;
                if (root1 > root2) distance = root1;
                else distance = root2;
            }
            else   //in this case we have a linear equation
            {
                distance = -c / b;
            }
        }

        if(distance < 0.0)
            distance = 1e-15;

        geom[unknown_node_index].SetLock();
        geom[unknown_node_index].FastGetSolutionStepValue(rDistanceVar) += distance*nodal_vol;
        geom[unknown_node_index].FastGetSolutionStepValue(rAreaVar) += nodal_vol;
        geom[unknown_node_index].UnSetLock();
    }




    ///@}

    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void Check(
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

    void ResetVariables(
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

    void CalculateExactDistancesOnDividedElements(
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

    void AbsDistancesOnDividedElements(
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

	void ExtendDistancesByLayer(
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

    void AssignDistanceSign(
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

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ParallelDistanceCalculator<TDim>& operator=(ParallelDistanceCalculator<TDim> const& rOther) {};

    /// Copy constructor.
    ParallelDistanceCalculator(ParallelDistanceCalculator<TDim> const& rOther) {};


    ///@}

}; // Class ParallelDistanceCalculator

///@}

///@name Type Definitions
///@{

template< unsigned int TDim>
const Kratos::Flags ParallelDistanceCalculator<TDim>::CALCULATE_EXACT_DISTANCES_TO_PLANE(Kratos::Flags::Create(0));

///@}
///@name Input and output
///@{


/// input stream function
template<unsigned int TDim>
inline std::istream& operator >> (std::istream& rIStream,
                                  ParallelDistanceCalculator<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<unsigned int TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ParallelDistanceCalculator<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED  defined
