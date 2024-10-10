//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors: Riccardo Rossi
//                Joaquín Irazábal
//
//

#ifndef _EDGEBASED_LEVELSET_AUXILIARY_UTILS_H
#define	_EDGEBASED_LEVELSET_AUXILIARY_UTILS_H
/* System includes */


/* External includes */


/* Project includes */
//#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "utilities/body_distance_calculation_utils.h"


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/** Short class definition.
Detail class definition.

*/
template< unsigned int TDim>
class EdgebasedLevelsetAuxiliaryUtils
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */

    /** Destructor.
    */

    /*@} */
    /**@name Operators
    */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */

    void CalculateDistances(
        ModelPart& r_model_part,
        Variable<double>& rDistanceVar,
        const double max_distance
    )
    {
        if constexpr (TDim == 2)
            CalculateDistances2D(r_model_part, rDistanceVar, max_distance);
        else if constexpr (TDim == 3)
            CalculateDistances3D(r_model_part, rDistanceVar, max_distance);

     //   r_model_part.GetCommunicator().SynchronizeCurrentDataToMin(rDistanceVar);
    }

    //***********************************************************************
    //***********************************************************************
    void CalculateDistances2D(
        ModelPart& r_model_part,
        Variable<double>& rDistanceVar,
        const double max_distance
    )
    {
        KRATOS_TRY

        double large_distance = 1e6;
        double tol = 1.0/large_distance;

        //copy rDistance var to GetValue database
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            it->GetValue(rDistanceVar) = it->FastGetSolutionStepValue(rDistanceVar);
            it->FastGetSolutionStepValue(rDistanceVar) = large_distance;
            it->GetValue(IS_VISITED) = 0;
        }

        BoundedMatrix<double,TDim+1,TDim> DN_DX;
        array_1d<double,TDim+1> N, distances;
        array_1d<double,TDim> grad_d;
        array_1d<double,3> coord_on_0(3,0.0);
        array_1d<double,3> temp;
        BoundedMatrix<double, TDim, TDim> free_surface_points;

        //fill the list of the first elements to be solved for the "distance"
        for(ElementsArrayType::iterator it =  r_model_part.ElementsBegin(); it !=r_model_part.ElementsEnd(); it++)
        {
            Element::GeometryType& geom = it->GetGeometry();

            unsigned int n_positive = 0;
            unsigned int n_negative = 0;
            for(unsigned int kk = 0; kk<TDim+1 ; kk++)
            {
                const double dist = geom[kk].GetValue(rDistanceVar);
                distances[kk] = dist;
                if(dist < 0)
                    n_negative += 1;
                else
                    n_positive += 1;
            }

            if(n_negative > 0 && n_positive > 0) //ELEMENT IS CROSSED BY THE INTERFACE!
            {
                double Volume;
                GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

                //compute the gradient of the distance and normalize it
                noalias(grad_d) = prod(trans(DN_DX),distances);
                double norm = norm_2(grad_d);
                grad_d /= norm;

                //find one division point on one edge
                for(unsigned int i = 1; i<TDim+1; i++)
                {
                    if(distances[0]*distances[i]<=0) //if the edge is divided
                    {
                        double delta_d = fabs(distances[i]) + fabs(distances[0]);

                        if(delta_d>1e-15)
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

                    double& dist_i = geom[i].FastGetSolutionStepValue(rDistanceVar);
                    if(real_distance < dist_i)
                        dist_i = real_distance;

                    //select the nodes for the computation of the distance
                    geom[i].GetValue(IS_VISITED) = 1;
                }

            }



        }

        //loop on all nodes to treat correctly the case of the surface coinciding with a node
        //copy rDistance var to GetValue database
        array_1d<double,3> aux;
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(fabs(it->GetValue(rDistanceVar)) < tol)
            {
                it->FastGetSolutionStepValue(rDistanceVar) = 0.0;
                it->GetValue(IS_VISITED) = 1;
            }
        }

        //compute the distance using the element based approach
        BodyDistanceCalculationUtils util;
        util.CalculateDistances<TDim>(r_model_part.Elements(),rDistanceVar,max_distance);

        //finally change the sign to the distance as needed
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(it->GetValue(rDistanceVar) < 0)
                it->FastGetSolutionStepValue(rDistanceVar) = -it->FastGetSolutionStepValue(rDistanceVar);
        }

        //check if something is wrong
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(fabs(it->FastGetSolutionStepValue(rDistanceVar)) == large_distance)
            {
                KRATOS_WATCH("error in the calculation of the distance for node")
                KRATOS_WATCH(it->Id());

            }
        }

        KRATOS_CATCH("")
    }

    void CalculateDistances3D(
        ModelPart& r_model_part,
        Variable<double>& rDistanceVar,
        const double max_distance
    )
    {
        KRATOS_TRY

        const double large_distance = 1e9;

        std::vector<array_1d<double,TDim + 1> > distances;
        std::vector<Element::Pointer> interface_elements;

         //copy rDistance var to GetValue database
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            it->GetValue(rDistanceVar) = it->FastGetSolutionStepValue(rDistanceVar);
            it->GetValue(IS_VISITED) = 0;
        }

        // Make a list of all elements crossed by interface and get their nodal distances
        for(ElementsArrayType::iterator i_element =  r_model_part.ElementsBegin(); i_element !=r_model_part.ElementsEnd(); i_element++)
        {
            // Get a reference to the geometry
            Element::GeometryType& element_geometry = i_element->GetGeometry();
            array_1d<double,TDim + 1> element_distance;

            bool positive = false;
            bool negative = false;

            // loop over nodes to see if they have positive or negative distance.
            for(unsigned int i_node = 0; i_node < element_geometry.size() ; i_node++)
            {
                const double distance = element_geometry[i_node].GetSolutionStepValue(rDistanceVar);
                element_distance[i_node] = distance;
                negative |= (distance < 0);
                positive |= (distance >= 0);
            }

            if(negative && positive) //ELEMENT IS CROSSED BY THE INTERFACE!
            {
                interface_elements.push_back(*(i_element.base()));
                distances.push_back(element_distance);
            }
        }

        // Reset the distance to a large value maintaining the sign
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            double& distance = it->FastGetSolutionStepValue(rDistanceVar);

            if(distance >= 0.00)
                distance = large_distance;
            else
                distance = -large_distance;
        }

        // Now calculating the new distance for interface elements
        int size = static_cast<int>(interface_elements.size());
        for(int i = 0 ; i < size ; i++)
        {
            // Get a reference to the geometry
            Element::GeometryType& element_geometry = interface_elements[i]->GetGeometry();

            GeometryUtils::CalculateTetrahedraDistances(element_geometry, distances[i]);

            // loop over nodes and apply the new distances.
            for(unsigned int i_node = 0; i_node < element_geometry.size() ; i_node++)
            {
                double& distance = element_geometry[i_node].GetSolutionStepValue(rDistanceVar);
                double new_distance = distances[i][i_node];

                if(fabs(distance) > fabs(new_distance))
                {
                    if(distance >= 0)
                        distance = new_distance;
                    else
                        distance = new_distance;
                }
                element_geometry[i_node].GetValue(IS_VISITED) = 1;
            }
       }
        //compute the distance using the element based approach
        BodyDistanceCalculationUtils util;
        util.CalculateDistances<TDim>(r_model_part.Elements(),rDistanceVar,max_distance);

        //finally change the sign to the distance as needed
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(it->GetValue(rDistanceVar) < 0)
                it->FastGetSolutionStepValue(rDistanceVar) = -it->FastGetSolutionStepValue(rDistanceVar);
        }

        //check if something is wrong
        for(ModelPart::NodesContainerType::iterator it =  r_model_part.NodesBegin(); it !=r_model_part.NodesEnd(); it++)
        {
            if(fabs(it->FastGetSolutionStepValue(rDistanceVar)) == large_distance)
            {
                KRATOS_WATCH("error in the calculation of the distance for node")
                KRATOS_WATCH(it->Id());

            }
        }

        KRATOS_CATCH("")
    }

    //**********************************************************************************+
    //**********************************************************************************+
    double FindMaximumEdgeSize(ModelPart& r_model_part)
    {
        KRATOS_TRY

        ModelPart::NodesContainerType& rNodes = r_model_part.Nodes();

        double h_max = 0.0;

        for (ModelPart::NodesContainerType::iterator in = rNodes.begin(); in != rNodes.end(); in++)
        {
            double xc = in->X();
            double yc = in->Y();
            double zc = in->Z();

            double h = 0.0;
            for (GlobalPointersVector< Node >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                    i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
            {
                double x = i->X();
                double y = i->Y();
                double z = i->Z();
                double l = (x - xc)*(x - xc);
                l += (y - yc)*(y - yc);
                l += (z - zc)*(z - zc);

                if (l > h) h = l;
            }
            h = sqrt(h);

            if(h > h_max) h_max = h;

        }

        return h_max;

        KRATOS_CATCH("");
    }

    /*@} */
    /**@name Acces */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Acces */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/


#endif	/* _EDGEBASED_LEVELSET_AUXILIARY_UTILS_H */

