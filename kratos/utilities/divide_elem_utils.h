//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    
//                    
//



#if !defined(KRATOS_DIVIDE_ELEM_UTILS )
#define  KRATOS_DIVIDE_ELEM_UTILS


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"


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

///Various mathematical utilitiy functions.
/**
* Defines several utility functions
 */
class DivideElemUtils
{
public:

// 		BoundedMatrix<double,3,2>& DN_DX
    static void DivideElement_2D( Element::GeometryType& geom,
// 					const double toll,
// 					double is_divided,
                                  BoundedMatrix<double,4,2>& aux_gp,
                                  array_1d<double,4>& A_on_agp,
                                  BoundedMatrix<double,4,3>&  N_on_agp,
                                  array_1d<double,4>& dist_on_agp
                                )
    {
        KRATOS_TRY
        //clearing auxiliar vectors

        aux_gp = ZeroMatrix(4,2);
        A_on_agp = ZeroVector(4);
        N_on_agp = ZeroMatrix(4,3);
        dist_on_agp = ZeroVector(4);

        double x0 = geom[0].X();
        double y0 = geom[0].Y();
        double x1 = geom[1].X();
        double y1 = geom[1].Y();
        double x2 = geom[2].X();
        double y2 = geom[2].Y();

        double Area = ((x1-x0)*(y2-y0)-(x2-x0)*(y1-y0))*0.5;
// 				KRATOS_WATCH(Area);


        double toll =  0.1*sqrt(Area * 2.30940108);//0.1*(h in a equailateral triangle of given area)
        array_1d<double,3> dist = ZeroVector(3);
        for (unsigned int i = 0; i < dist.size(); i++)
        {
            dist[i] = geom[i].GetSolutionStepValue(DISTANCE);
            if(dist[i] < 0.0 && dist[i] > -toll)
                dist[i] = 0.0;
        }



        double sign0 = dist[1] * dist[2];
        double sign1 = dist[2] * dist[0];
        double sign2 = dist[0] * dist[1];

        array_1d<double,3> abs_dist = ZeroVector(3);
        for (unsigned int i = 0;  i< dist.size(); i++)
            abs_dist[i] = fabs(dist[i]);

        //Calculating the 3 VIRTUAL NODES:
        //if the edge intersects the free surface: proportional calculation
        //otherwise mean point of the edge
        BoundedMatrix<double,3,2> virt_nodes = ZeroMatrix(3,2);
// 				noalias(virt_nodes) = ZeroMatrix(3,2);
        if(sign0 < 0)
        {
            virt_nodes(0,0) = geom[1].X() + (geom[2].X() - geom[1].X()) * abs_dist[1]/(abs_dist[1]+abs_dist[2]);
            virt_nodes(0,1) = geom[1].Y() + (geom[2].Y() - geom[1].Y()) * abs_dist[1]/(abs_dist[1]+abs_dist[2]);
        }
        else
        {
            virt_nodes(0,0) = (geom[1].X() + geom[2].X())*0.5;
            virt_nodes(0,1) = (geom[1].Y() + geom[2].Y())*0.5;
        }
        if(sign1 < 0)
        {
            virt_nodes(1,0) = geom[2].X() + (geom[0].X() - geom[2].X()) * abs_dist[2]/(abs_dist[2]+abs_dist[0]);
            virt_nodes(1,1) = geom[2].Y() + (geom[0].Y() - geom[2].Y()) * abs_dist[2]/(abs_dist[2]+abs_dist[0]);
        }
        else
        {
            virt_nodes(1,0) = (geom[2].X() + geom[0].X())*0.5;
            virt_nodes(1,1) = (geom[2].Y() + geom[0].Y())*0.5;
        }
        if(sign2 < 0)
        {
            virt_nodes(2,0) = geom[0].X() + (geom[1].X() - geom[0].X()) * abs_dist[0]/(abs_dist[0]+abs_dist[1]);
            virt_nodes(2,1) = geom[0].Y() + (geom[1].Y() - geom[0].Y()) * abs_dist[0]/(abs_dist[0]+abs_dist[1]);
        }
        else
        {
            virt_nodes(2,0) = (geom[0].X() + geom[1].X())*0.5;
            virt_nodes(2,1) = (geom[0].Y() + geom[1].Y())*0.5;
        }


        double x3 = virt_nodes(0,0);
        double y3 = virt_nodes(0,1);
        double x4 = virt_nodes(1,0);
        double y4 = virt_nodes(1,1);
        double x5 = virt_nodes(2,0);
        double y5 = virt_nodes(2,1);
// 				KRATOS_WATCH(virt_nodes);

        if(dist[0] == 0.0)
        {
            //Calculating coordinates of the Auxiliary Gauss Points of the 4 sub elements
            aux_gp(0,0) = (x0 + x5 + x3)*0.3333333333;
            aux_gp(0,1) = (y0 + y5 + y3)*0.3333333333;
            aux_gp(1,0) = (x5 + x1 + x3)*0.3333333333;
            aux_gp(1,1) = (y5 + y1 + y3)*0.3333333333;
            aux_gp(2,0) = (x4 + x3 + x2)*0.3333333333;
            aux_gp(2,1) = (y4 + y3 + y2)*0.3333333333;
            aux_gp(3,0) = (x0 + x3 + x4)*0.3333333333;
            aux_gp(3,1) = (y0 + y3 + y4)*0.3333333333;
            // 				KRATOS_WATCH(aux_gp);

            //Calculating the area of the virtual elements
            A_on_agp[0] = ((x5 - x0)*(y3 - y0)-(x3 - x0)*(y5 - y0))*0.5;
            A_on_agp[1] = ((x1 - x5)*(y3 - y5)-(x3 - x5)*(y1 - y5))*0.5;
            A_on_agp[2] = ((x3 - x4)*(y2 - y4)-(x2 - x4)*(y3 - y4))*0.5;
            A_on_agp[3] = ((x3 - x0)*(y4 - y0)-(x4 - x0)*(y3 - y0))*0.5;
            for (unsigned int i = 0; i < A_on_agp.size(); i++)
            {
                if(A_on_agp[i] <= 0.0)
                {
                    KRATOS_WATCH(geom[0].Id())
                    KRATOS_WATCH(geom[1].Id())
                    KRATOS_WATCH(geom[2].Id())
                    KRATOS_WATCH("NODE ZERO ZERO DISTANCE")
                    KRATOS_WATCH(A_on_agp)
                    KRATOS_WATCH("NEGATIVE AREAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                }
            }

        }
        else if (dist[1] == 0.0)
        {
            //Calculating coordinates of the Auxiliary Gauss Points of the 4 sub elements
            aux_gp(0,0) = (x0 + x5 + x4)*0.3333333333;
            aux_gp(0,1) = (y0 + y5 + y4)*0.3333333333;
            aux_gp(1,0) = (x5 + x1 + x4)*0.3333333333;
            aux_gp(1,1) = (y5 + y1 + y4)*0.3333333333;
            aux_gp(2,0) = (x4 + x3 + x2)*0.3333333333;
            aux_gp(2,1) = (y4 + y3 + y2)*0.3333333333;
            aux_gp(3,0) = (x1 + x3 + x4)*0.3333333333;
            aux_gp(3,1) = (y1 + y3 + y4)*0.3333333333;
            // 				KRATOS_WATCH(aux_gp);

            //Calculating the area of the virtual elements
            A_on_agp[0] = ((x5 - x0)*(y4 - y0)-(x4 - x0)*(y5 - y0))*0.5;
            A_on_agp[1] = ((x1 - x5)*(y4 - y5)-(x4 - x5)*(y1 - y5))*0.5;
            A_on_agp[2] = ((x3 - x4)*(y2 - y4)-(x2 - x4)*(y3 - y4))*0.5;
            A_on_agp[3] = ((x3 - x1)*(y4 - y1)-(x4 - x1)*(y3 - y1))*0.5;
            for (unsigned int i = 0; i < A_on_agp.size(); i++)
            {
                if(A_on_agp[i] <= 0.0)
                {
                    KRATOS_WATCH(geom[0].Id())
                    KRATOS_WATCH(geom[1].Id())
                    KRATOS_WATCH(geom[2].Id())
                    KRATOS_WATCH("NODE 1 ZERO DISTANCE")
                    KRATOS_WATCH(A_on_agp)
                    KRATOS_WATCH("NEGATIVE AREAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                }
            }
        }
        else if (dist[2] == 0.0)
        {
            //Calculating coordinates of the Auxiliary Gauss Points of the 4 sub elements
            aux_gp(0,0) = (x0 + x5 + x4)*0.3333333333;
            aux_gp(0,1) = (y0 + y5 + y4)*0.3333333333;
            aux_gp(1,0) = (x5 + x1 + x3)*0.3333333333;
            aux_gp(1,1) = (y5 + y1 + y3)*0.3333333333;
            aux_gp(2,0) = (x4 + x5 + x2)*0.3333333333;
            aux_gp(2,1) = (y4 + y5 + y2)*0.3333333333;
            aux_gp(3,0) = (x5 + x3 + x2)*0.3333333333;
            aux_gp(3,1) = (y5 + y3 + y2)*0.3333333333;
            // 				KRATOS_WATCH(aux_gp);

            //Calculating the area of the virtual elements
            A_on_agp[0] = ((x5 - x0)*(y4 - y0)-(x4 - x0)*(y5 - y0))*0.5;
            A_on_agp[1] = ((x1 - x5)*(y3 - y5)-(x3 - x5)*(y1 - y5))*0.5;
            A_on_agp[2] = ((x5 - x4)*(y2 - y4)-(x2 - x4)*(y5 - y4))*0.5;
            A_on_agp[3] = ((x3 - x5)*(y2 - y5)-(x2 - x5)*(y3 - y5))*0.5;
            for (unsigned int i = 0; i < A_on_agp.size(); i++)
            {
                if(A_on_agp[i] <= 0.0)
                {
                    KRATOS_WATCH(geom[0].Id())
                    KRATOS_WATCH(geom[1].Id())
                    KRATOS_WATCH(geom[2].Id())
                    KRATOS_WATCH("NODE 2 ZERO DISTANCE")
                    KRATOS_WATCH(A_on_agp)
                    KRATOS_WATCH("NEGATIVE AREAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                }
            }
        }
        else
        {
            //Calculating coordinates of the Auxiliary Gauss Points of the 4 sub elements
            aux_gp(0,0) = (x0 + x5 + x4)*0.3333333333;
            aux_gp(0,1) = (y0 + y5 + y4)*0.3333333333;
            aux_gp(1,0) = (x5 + x1 + x3)*0.3333333333;
            aux_gp(1,1) = (y5 + y1 + y3)*0.3333333333;
            aux_gp(2,0) = (x4 + x3 + x2)*0.3333333333;
            aux_gp(2,1) = (y4 + y3 + y2)*0.3333333333;
            aux_gp(3,0) = (x5 + x3 + x4)*0.3333333333;
            aux_gp(3,1) = (y5 + y3 + y4)*0.3333333333;
// 				KRATOS_WATCH(aux_gp);

            //Calculating the area of the virtual elements
            A_on_agp[0] = ((x5 - x0)*(y4 - y0)-(x4 - x0)*(y5 - y0))*0.5;
            A_on_agp[1] = ((x1 - x5)*(y3 - y5)-(x3 - x5)*(y1 - y5))*0.5;
            A_on_agp[2] = ((x3 - x4)*(y2 - y4)-(x2 - x4)*(y3 - y4))*0.5;
            A_on_agp[3] = ((x3 - x5)*(y4 - y5)-(x4 - x5)*(y3 - y5))*0.5;
// 			       KRATOS_WATCH(A_on_agp)
            for (unsigned int i = 0; i < A_on_agp.size(); i++)
            {
                if(A_on_agp[i] <= 0.0)
                {
                    KRATOS_WATCH(geom[0].Id())
                    KRATOS_WATCH(geom[1].Id())
                    KRATOS_WATCH(geom[2].Id())
                    KRATOS_WATCH("NO ZERO DISTANCE")
                    KRATOS_WATCH(A_on_agp)
                    KRATOS_WATCH("NEGATIVE AREAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                }
            }
        }




// 				KRATOS_WATCH(A_on_agp);
        //Evaluating the element shape functions on the auxiliary gauss points

        N_on_agp(0,0) = ((x1 - aux_gp(0,0)) * (y2 - aux_gp(0,1))-(x2 - aux_gp(0,0)) * (y1 - aux_gp(0,1)))*0.5/Area;
        N_on_agp(0,1) = ((aux_gp(0,0) - x0) * (y2 - y0)-(x2 - x0) * (aux_gp(0,1) - y0))*0.5/Area;
        N_on_agp(0,2) = ((x1 - x0) * (aux_gp(0,1) - y0)-(aux_gp(0,0) - x0) * (y1 - y0))*0.5/Area;

        N_on_agp(1,0) = ((x1 - aux_gp(1,0)) * (y2 - aux_gp(1,1))-(x2 - aux_gp(1,0)) * (y1 - aux_gp(1,1)))*0.5/Area;
        N_on_agp(1,1) = ((aux_gp(1,0) - x0) * (y2 - y0)-(x2 - x0) * (aux_gp(1,1) - y0))*0.5/Area;
        N_on_agp(1,2) = ((x1 - x0) * (aux_gp(1,1) - y0)-(aux_gp(1,0) - x0) * (y1 - y0))*0.5/Area;

        N_on_agp(2,0) = ((x1 - aux_gp(2,0)) * (y2 - aux_gp(2,1))-(x2 - aux_gp(2,0)) * (y1 - aux_gp(2,1)))*0.5/Area;
        N_on_agp(2,1) = ((aux_gp(2,0) - x0) * (y2 - y0)-(x2 - x0) * (aux_gp(2,1) - y0))*0.5/Area;
        N_on_agp(2,2) = ((x1 - x0) * (aux_gp(2,1) - y0)-(aux_gp(2,0) - x0) * (y1 - y0))*0.5/Area;

        N_on_agp(3,0) = ((x1 - aux_gp(3,0)) * (y2 - aux_gp(3,1))-(x2 - aux_gp(3,0)) * (y1 - aux_gp(3,1)))*0.5/Area;
        N_on_agp(3,1) = ((aux_gp(3,0) - x0) * (y2 - y0)-(x2 - x0) * (aux_gp(3,1) - y0))*0.5/Area;
        N_on_agp(3,2) = ((x1 - x0) * (aux_gp(3,1) - y0)-(aux_gp(3,0) - x0) * (y1 - y0))*0.5/Area;

// 				KRATOS_WATCH(N_on_agp);

        //Calculating DISTANCE on the aGP
        for (unsigned int i = 0; i < dist_on_agp.size(); i++)
        {
            for (unsigned j = 0; j < dist.size(); j++)
            {
                dist_on_agp[i] += N_on_agp(i,j) * dist[j];
            }
        }
// 				KRATOS_WATCH(dist);
// 				KRATOS_WATCH(dist_on_agp);




        KRATOS_CATCH("")
    }

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

#endif /* KRATOS_DIVIDE_ELEM_UTILS  defined */

