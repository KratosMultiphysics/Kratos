/*
==============================================================================
KratosPFEMApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


/* *********************************************************
*
*   Last Modified by:    $Author: antonia $
*   Date:                $Date: 2009-10-19 12:12:16 $
*   Revision:            $Revision: 0.1 $
*
* ***********************************************************/

#if !defined(KRATOS_EROSION_UTILS )
#define  KRATOS_EROSION_UTILS


/* System includes */
#include <string>
#include <iostream>
#include <stdlib.h>
#include <algorithm>

/* External includes */


/* Project includes */
#include "utilities/math_utils.h"
#include "includes/model_part.h"
#include "includes/define.h"
// #include "utilities/projection.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/deprecated_variables.h"
#include "PFEM_application.h"
#include "includes/deprecated_variables.h"

//Database includes
#include "spatial_containers/spatial_containers.h"

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

  \URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

	\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

	  \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

		\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


			\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

			  \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

				\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

				  \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


*/
template< std::size_t TDim>
class ErosionUtils
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

    bool CheckErosionableNodes(
        ModelPart& rOrigin_ModelPart ,
        ModelPart& rDestination_ModelPart,
        const Vector& rCriticalVel, //TO BE INSERTED
        const double rCriticalEnergy,
        const double fluid_nu,
        const double fluid_rho
    )
    {

        KRATOS_TRY
// KRATOS_WATCH(	"Line 170")
// 			//properties to be used in the generation
// 			Properties::Pointer properties = rDestination_ModelPart.GetMesh().pGetProperties(1);
//
// 			//defintions for spatial search
        typedef Node<3> PointType;
        typedef Node<3>::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer>           PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double>               DistanceVector;
        typedef std::vector<double>::iterator     DistanceIterator;

        //*************
        // Bucket types
        typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
// 			   typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
// 			   typedef BinsDynamic< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
        //*************
        // DynamicBins;

        typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;
// 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
// 			   typedef Tree< StaticBins > tree; 		     	//Binstree;
// 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
// 			   typedef typename KdtreeBins::Partitions SubPartitions;
// 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
        /*
        			   typedef Bins< TDim, PointType, stdPointVector> stdBins;
        			   typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/
// KRATOS_WATCH(	"Line 199")

        for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
                node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
        {
            Node<3>::Pointer pnode = *(node_it.base());
            pnode->GetValue(IS_VISITED) = 0.0;
        }
        //********************************************************************
        //Build an auxiliary ModelPart with all the PFEM free surface nodes
        //********************************************************************
        PointVector list_of_erosionable_nodes;

        for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
                node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
        {

            //PointType::Pointer pnode(new PointType(*node_it));
            Node<3>::Pointer pnode = *(node_it.base());

            if(pnode->FastGetSolutionStepValue(IS_EROSIONABLE)==1.0 && pnode->FastGetSolutionStepValue(IS_STRUCTURE)!=1.0)
            {
                //putting the surface nodes of the destination_model part in an auxiliary list
                list_of_erosionable_nodes.push_back( pnode );

                WeakPointerVector< Node<3> >& neighb_nodes = pnode->GetValue(NEIGHBOUR_NODES);

                for( WeakPointerVector< Node<3> >::iterator j = neighb_nodes.begin(); j != neighb_nodes.end(); j++)
                {


                    if(j->GetValue(IS_VISITED)==0.0 && j->FastGetSolutionStepValue(IS_STRUCTURE)!=1.0 && j->FastGetSolutionStepValue(IS_FREE_SURFACE)!=1.0)
                    {
                        list_of_erosionable_nodes.push_back( Node<3>::Pointer( *(j.base() ) ) );


                    }
                }
            }
        }

        //********************************************************************
        //Calculate velocity of the free surface nodes by direct variable interpolation
        //********************************************************************
        //create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;

        tree nodes_tree(list_of_erosionable_nodes.begin(),list_of_erosionable_nodes.end(),bucket_size);
        //work arrays
        Node<3> work_point(0,0.0,0.0,0.0);
        unsigned int MaximumNumberOfResults = 10000;
        PointVector Results(MaximumNumberOfResults);
        DistanceVector ResultsDistances(MaximumNumberOfResults);

        array_1d<double,TDim+1> N; //Shape functions vector//

        double critical_vel2 = 0.0;
        for (unsigned int i=0; i< TDim; i++)
            critical_vel2 +=  rCriticalVel[i]*rCriticalVel[i];
        double fluidified_density = 1000.0;
        double fluidified_viscosity = 0.000001;
        int rCalculateDam = 0;

        //loop over all of the elements in the "old" list to perform the interpolation
        for( ModelPart::ElementsContainerType::iterator el_it = rOrigin_ModelPart.ElementsBegin();
                el_it != rOrigin_ModelPart.ElementsEnd(); el_it++)
        {
            Geometry<Node<3> >&geom = el_it->GetGeometry();
            //check if the element is a fluid element or an element of the extrapolation domain
            bool is_fluid = false;
            unsigned int fluid_nodes_num = 0;
// 				double fluid_nu = 0.0;
// 				double fluid_rho = 0.0;

            for (unsigned int i = 0; i<TDim+1 ; i++)
            {
                if (geom[i].FastGetSolutionStepValue(DISTANCE) <= 0)
                {
                    fluid_nodes_num ++;
                }
            }
            if(fluid_nodes_num >= TDim )
            {
                is_fluid = true;

                //find the center and "radius" of the element
                double xc, yc, zc, radius;
                CalculateCenterAndSearchRadius( geom,xc,yc,zc, radius, N);

                //find all of the new nodes within the radius
                int number_of_points_in_radius;
                work_point.X() = xc;
                work_point.Y() = yc;

                //look between the new nodes which of them is inside the radius of the circumscribed cyrcle
                number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                                             ResultsDistances.begin(),  MaximumNumberOfResults);
                //check if inside
                for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
                {

                    bool is_inside = false;
                    //once we are sure the node in inside the circle we have to see if it is inside the triangle i.e. if all the Origin element shape functions are >1
                    double is_visited = (*it_found)->GetValue(IS_VISITED);
                    is_inside = CalculatePosition(geom,	(*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);

                    //if the node falls inside the element interpolate
                    if(is_inside == true && is_visited != 1.0 && is_fluid == true)
                    {


                        array_1d<double,3> InterpolatedVelocity;

                        //Interpolating all the rVariables of the rOrigin_ModelPart to get their nodal value in the rDestination_ModelPart
                        Interpolate(  el_it,  N, *it_found , InterpolatedVelocity , VELOCITY );

                        double fluid_vel2 = 0.0;
                        // 						    //Interpolated vel (not node real velocity)
                        for (unsigned int i=0; i< TDim; i++)
                        {
                            fluid_vel2 +=  InterpolatedVelocity[i] * InterpolatedVelocity[i];
                        }

                        //********************************************************************
                        //Check if the velocity of PFEM free surface node is >= v_critical (Shield)
                        //********************************************************************

                        double dt = rDestination_ModelPart.GetProcessInfo()[DELTA_TIME];
                        double coeff = 0.25 * fluid_nu * dt * fluid_rho * fluid_vel2;

                        if(fluid_vel2 > critical_vel2 || (*it_found)->FastGetSolutionStepValue(FRICTION_COEFFICIENT) > 0.00001 )
                        {


// 							//Starting the pfem solution
                            if(rCalculateDam == 0)
                            {
                                rCalculateDam = 1;
// 							    KRATOS_WATCH("Solving PFEM stucture too ::::::::::::::::::::::::::::::::::::::::::::")
                            }
                            (*it_found)->FastGetSolutionStepValue(FRICTION_COEFFICIENT) += coeff;


                        }

                        double volumetric_parameter = (*it_found)->FastGetSolutionStepValue(NODAL_H);
                        for(unsigned int i = 0; i<TDim-1; i++)
                            volumetric_parameter *= (*it_found)->FastGetSolutionStepValue(NODAL_H);


                        if( (*it_found)->FastGetSolutionStepValue(FRICTION_COEFFICIENT)>=(rCriticalEnergy * volumetric_parameter))
                        {
                            (*it_found)->FastGetSolutionStepValue(VISCOSITY) = fluidified_viscosity;
                            (*it_found)->FastGetSolutionStepValue(DENSITY) = fluidified_density;

// 							//In order not to allow the friction coefficient to increase too much
                            (*it_found)->FastGetSolutionStepValue(FRICTION_COEFFICIENT) = 100*(rCriticalEnergy * volumetric_parameter);
                        }
                        //}
                    }
                }
            }
        }


        //CHANGE DENSITY AND VISCOSITY IF ALL THE NODES OF AN ELEMENT ARE FREE SURFACE
        for( ModelPart::ElementsContainerType::iterator el_it = rDestination_ModelPart.ElementsBegin();
                el_it != rDestination_ModelPart.ElementsEnd(); el_it++)
        {
            Geometry< Node<3> >& geom = el_it->GetGeometry();
            unsigned int n_freesurf_nodes = 0;
            for(unsigned int i =0; i<TDim+1; i++)
            {
                if(geom[i].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
                    n_freesurf_nodes ++;
            }
            if(n_freesurf_nodes == TDim+1)
            {
                for(unsigned int i =0; i<TDim+1; i++)
                {
                    geom[i].FastGetSolutionStepValue(VISCOSITY) = fluidified_viscosity;
                    geom[i].FastGetSolutionStepValue(DENSITY) = fluidified_density;

                }
            }
        }
        //CHANGE DENSITY AND VISCOSITY IF A NODE HAS MORE THAN 3 NEIGHBOURS THAT ARE FREE SURFACE
        for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
                node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
        {
            double count = 0.0;
            //PointType::Pointer pnode(new PointType(*node_it));
            Node<3>::Pointer pnode = *(node_it.base());
            WeakPointerVector< Node<3> >& neighb_nodes = pnode->GetValue(NEIGHBOUR_NODES);
            for( WeakPointerVector< Node<3> >::iterator j = neighb_nodes.begin(); j != neighb_nodes.end(); j++)
            {
                if( j->FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
                    count ++;
            }
            if (count >= (TDim + 2))
            {
                (pnode)->FastGetSolutionStepValue(VISCOSITY) = fluidified_viscosity;
                (pnode)->FastGetSolutionStepValue(DENSITY) = fluidified_density;
            }
        }

        if(rCalculateDam == 1)
            return true;

        return false;

        KRATOS_CATCH("")

// KRATOS_WATCH("check erosion finished+++++++++++++++++++++++++++++++++++++++++++++++++")
    }






    void SetErosionableNodes(	ModelPart& rModelPart)
    {
        KRATOS_TRY

// 			for(ModelPart::NodesContainerType::iterator node_it = rModelPart.NodesBegin();
// 						node_it != rModelPart.NodesEnd(); ++node_it)
// 			{
// 				//PointType::Pointer pnode(new PointType(*node_it));
//  				Node<3>::Pointer pnode = *(node_it.base());
// 				if(pnode->FastGetSolutionStepValue(IS_EROSIONABLE)==1.0)
// 				{
// 					WeakPointerVector< Node<3> >& neighb_nodes = pnode->GetValue(NEIGHBOUR_NODES);
// 					for( WeakPointerVector< Node<3> >::iterator j = neighb_nodes.begin(); j != neighb_nodes.end(); j++)
// 					{
// 					     if(j->FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
// 					     {
// 						  j->FastGetSolutionStepValue(IS_EROSIONABLE) = 1.0;
// 					     }
//
// 					}
// 				}
// 			}
        KRATOS_CATCH("")
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
    ///@name Static Member rVariables
    ///@{


    ///@}
    ///@name Member rVariables
    ///@{
    inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
            double& xc, double& yc, double& zc, double& R, array_1d<double,3>& N
                                              )
    {
        double x0 = geom[0].X();
        double  y0 = geom[0].Y();
        double x1 = geom[1].X();
        double  y1 = geom[1].Y();
        double x2 = geom[2].X();
        double  y2 = geom[2].Y();


        xc = 0.3333333333333333333*(x0+x1+x2);
        yc = 0.3333333333333333333*(y0+y1+y2);
        zc = 0.0;

        double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
        double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
        double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);

        R = R1;
        if(R2 > R) R = R2;
        if(R3 > R) R = R3;

        R = 1.01 * sqrt(R);
    }
    //***************************************
    //***************************************
    inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
            double& xc, double& yc, double& zc, double& R, array_1d<double,4>& N

                                              )
    {
        double x0 = geom[0].X();
        double  y0 = geom[0].Y();
        double  z0 = geom[0].Z();
        double x1 = geom[1].X();
        double  y1 = geom[1].Y();
        double  z1 = geom[1].Z();
        double x2 = geom[2].X();
        double  y2 = geom[2].Y();
        double  z2 = geom[2].Z();
        double x3 = geom[3].X();
        double  y3 = geom[3].Y();
        double  z3 = geom[3].Z();


        xc = 0.25*(x0+x1+x2+x3);
        yc = 0.25*(y0+y1+y2+y3);
        zc = 0.25*(z0+z1+z2+z3);

        double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0) + (zc-z0)*(zc-z0);
        double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1) + (zc-z1)*(zc-z1);
        double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2) + (zc-z2)*(zc-z2);
        double R4 = (xc-x3)*(xc-x3) + (yc-y3)*(yc-y3) + (zc-z3)*(zc-z3);

        R = R1;
        if(R2 > R) R = R2;
        if(R3 > R) R = R3;
        if(R4 > R) R = R4;

        R = sqrt(R);
    }
    //***************************************
    //***************************************
    inline double CalculateVol(	const double x0, const double y0,
                                const double x1, const double y1,
                                const double x2, const double y2
                              )
    {
        return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
    }
    //***************************************
    //***************************************
    inline double CalculateVol(	const double x0, const double y0, const double z0,
                                const double x1, const double y1, const double z1,
                                const double x2, const double y2, const double z2,
                                const double x3, const double y3, const double z3
                              )
    {
        double x10 = x1 - x0;
        double y10 = y1 - y0;
        double z10 = z1 - z0;

        double x20 = x2 - x0;
        double y20 = y2 - y0;
        double z20 = z2 - z0;

        double x30 = x3 - x0;
        double y30 = y3 - y0;
        double z30 = z3 - z0;

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
        return  detJ*0.1666666666666666666667;

        //return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
    }
    //***************************************
    //***************************************
    inline bool CalculatePosition(	Geometry<Node<3> >&geom,
                                    const double xc, const double yc, const double zc,
                                    array_1d<double,3>& N
                                 )
    {
        double x0 = geom[0].X();
        double  y0 = geom[0].Y();
        double x1 = geom[1].X();
        double  y1 = geom[1].Y();
        double x2 = geom[2].X();
        double  y2 = geom[2].Y();

        double area = CalculateVol(x0,y0,x1,y1,x2,y2);
        double inv_area = 0.0;
        if(area == 0.0)
        {
            KRATOS_THROW_ERROR(std::logic_error,"element with zero area found","");
        }
        else
        {
            inv_area = 1.0 / area;
        }


        N[0] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
        N[1] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;
        N[2] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;


        if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <=1.0 && N[1]<= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
            return true;

        return false;
    }

    //***************************************
    //***************************************

    inline bool CalculatePosition(	Geometry<Node<3> >&geom,
                                    const double xc, const double yc, const double zc,
                                    array_1d<double,4>& N
                                 )
    {

        double x0 = geom[0].X();
        double  y0 = geom[0].Y();
        double  z0 = geom[0].Z();
        double x1 = geom[1].X();
        double  y1 = geom[1].Y();
        double  z1 = geom[1].Z();
        double x2 = geom[2].X();
        double  y2 = geom[2].Y();
        double  z2 = geom[2].Z();
        double x3 = geom[3].X();
        double  y3 = geom[3].Y();
        double  z3 = geom[3].Z();

        double vol = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

        double inv_vol = 0.0;
        if(vol < 0.0000000000001)
        {
            KRATOS_THROW_ERROR(std::logic_error,"element with zero vol found","");
        }
        else
        {
            inv_vol = 1.0 / vol;
        }

        N[0] = CalculateVol(x1,y1,z1,x3,y3,z3,x2,y2,z2,xc,yc,zc) * inv_vol;
        N[1] = CalculateVol(x3,y3,z3,x0,y0,z0,x2,y2,z2,xc,yc,zc) * inv_vol;
        N[2] = CalculateVol(x3,y3,z3,x1,y1,z1,x0,y0,z0,xc,yc,zc) * inv_vol;
        N[3] = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) * inv_vol;


        if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >=0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <=1.0)			//if the xc yc zc is inside the tetrahedron return true
            return true;

        return false;
    }


    static inline void Calculate_DN_DX(Geometry<Node<3> >&geom,
                                       boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX
                                      )
    {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();

        DN_DX(0,0) = -y20 + y10      ;
        DN_DX(0,1) = x20 - x10;
        DN_DX(1,0) =  y20	   ;
        DN_DX(1,1) = -x20     ;
        DN_DX(2,0) = -y10	   ;
        DN_DX(2,1) = x10	;
    }

    static inline void Calculate_DN_DX(Geometry<Node<3> >&geom,
                                       boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX
                                      )
    {
        double x10 = geom[1].X() - geom[0].X();
        double y10 = geom[1].Y() - geom[0].Y();
        double z10 = geom[1].Z() - geom[0].Z();

        double x20 = geom[2].X() - geom[0].X();
        double y20 = geom[2].Y() - geom[0].Y();
        double z20 = geom[2].Z() - geom[0].Z();

        double x30 = geom[3].X() - geom[0].X();
        double y30 = geom[3].Y() - geom[0].Y();
        double z30 = geom[3].Z() - geom[0].Z();

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

        DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
        DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
        DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
        DN_DX(1,0) = y20 * z30 - y30 * z20;
        DN_DX(1,1) = z20 * x30 - x20 * z30;
        DN_DX(1,2) = x20 * y30 - y20 * x30;
        DN_DX(2,0) = -y10 * z30 + z10 * y30;
        DN_DX(2,1) = x10 * z30 - z10 * x30;
        DN_DX(2,2) = -x10 * y30 + y10 * x30;
        DN_DX(3,0) = y10 * z20 - z10 * y20;
        DN_DX(3,1) = -x10 * z20 + z10 * x20;
        DN_DX(3,2) = x10 * y20 - y10 * x20;

        DN_DX /= detJ;	   ;
    }

    inline void CalculateNodalVelocityMatrix(	Geometry<Node<3> >&geom,
            boost::numeric::ublas::bounded_matrix<double,3,2>& vel
                                            )
    {

        //Nodal velocity
        const array_1d<double,3> vel0 = geom[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3> vel1 = geom[1].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3> vel2 = geom[2].FastGetSolutionStepValue(VELOCITY);


        for(unsigned int j=0; j<TDim; j++)
        {
            vel(0,j) = vel0[j];
            vel(1,j) = vel1[j];
            vel(2,j) = vel2[j];
        }
    }

    inline void CalculateNodalVelocityMatrix(	Geometry<Node<3> >&geom,
            boost::numeric::ublas::bounded_matrix<double,4,3>& vel
                                            )
    {

        //Nodal velocity
        const array_1d<double,3> vel0 = geom[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3> vel1 = geom[1].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3> vel2 = geom[2].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double,3> vel3 = geom[3].FastGetSolutionStepValue(VELOCITY);


        for(unsigned int j=0; j<TDim; j++)
        {
            vel(0,j) = vel0[j];
            vel(1,j) = vel1[j];
            vel(2,j) = vel2[j];
            vel(3,j) = vel3[j];
        }
    }




//el_it		     	Element iterator
//N			Shape functions
//step_data_size
//pnode			pointer to the node
    //total model part
    //projecting an array1D 2Dversion
    void Interpolate(
        ModelPart::ElementsContainerType::iterator el_it,
        const array_1d<double,3>& N,
        Node<3>::Pointer pnode,
        array_1d<double,3>& Interpolated_velocity,
        Variable<array_1d<double,3> >& rVariable)
    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();
// 			KRATOS_WATCH(geom[0].Id());
// 			KRATOS_WATCH(geom[1].Id());
// 			KRATOS_WATCH(geom[2].Id());
// 			unsigned int buffer_size = pnode->GetBufferSize();

// 			for(unsigned int step = 0; step<buffer_size; step++)
// 			{
        //getting the data of the solution step
// 				array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rVariable , step);
        //Reference or no reference???//CANCELLA
        const array_1d<double,3> node0_data = geom[0].FastGetSolutionStepValue(rVariable);
// KRATOS_WATCH(node0_data);
        const array_1d<double,3> node1_data = geom[1].FastGetSolutionStepValue(rVariable);
// KRATOS_WATCH(node1_data);
        const array_1d<double,3> node2_data = geom[2].FastGetSolutionStepValue(rVariable);
// KRATOS_WATCH(node2_data);


        //copying this data in the position of the vector we are interested in
        for(unsigned int j= 0; j< TDim; j++)
        {
// 					step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
            Interpolated_velocity[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];

        }
// KRATOS_WATCH(Interpolated_velocity);

// 			}
        pnode->GetValue(IS_VISITED) = 1.0;
    }

    //projecting an array1D 3Dversion
    void Interpolate(
        ModelPart::ElementsContainerType::iterator el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        array_1d<double,3>& Interpolated_velocity,
        Variable<array_1d<double,3> >& rVariable)
    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();

// 			unsigned int buffer_size = pnode->GetBufferSize();

// 			for(unsigned int step = 0; step<buffer_size; step++)
// 			{
        //getting the data of the solution step
        /*				array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rVariable , step);*/
        //Reference or no reference???//CANCELLA
        const array_1d<double,3> node0_data = geom[0].FastGetSolutionStepValue(rVariable);
        const array_1d<double,3> node1_data = geom[1].FastGetSolutionStepValue(rVariable);
        const array_1d<double,3> node2_data = geom[2].FastGetSolutionStepValue(rVariable);
        const array_1d<double,3> node3_data = geom[3].FastGetSolutionStepValue(rVariable);

        //copying this data in the position of the vector we are interested in
        for(unsigned int j= 0; j< TDim; j++)
        {
            Interpolated_velocity[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j] + N[3]*node3_data[j];
        }
// 			}
        pnode->GetValue(IS_VISITED) = 1.0;
    }


// 		//scalar
// 		void Interpolate(
// 				ModelPart::ElementsContainerType::iterator el_it,
// 				const array_1d<double,3>& N,
//       				Node<3>::Pointer pnode,
// 				Variable<double>& rVariable)
// 		{
// 			//Geometry element of the rOrigin_ModelPart
// 			Geometry< Node<3> >& geom = el_it->GetGeometry();
//
// 			unsigned int buffer_size = pnode->GetBufferSize();
// 		//facendo un loop sugli step temporali step_data come salva i dati al passo anteriore? Cioś dove passiamo l'informazione ai nodi???
// 			for(unsigned int step = 0; step<buffer_size; step++)
// 			{
// 				//getting the data of the solution step
// 				double& step_data = (pnode)->FastGetSolutionStepValue(rVariable , step);
// 				//Reference or no reference???//CANCELLA
// 				double& node0_data = geom[0].FastGetSolutionStepValue(rVariable , step);
// 				double& node1_data = geom[1].FastGetSolutionStepValue(rVariable , step);
// 				double& node2_data = geom[2].FastGetSolutionStepValue(rVariable , step);
//
// 				//copying this data in the position of the vector we are interested in
//
// 				step_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data;
//
// 			}
// 		}


// 		inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it , Variable<array_1d<double,3> >& rVariable)
// 		{
// 			array_1d<double, 3>& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
//
// 			noalias(Aux_var) = ZeroVector(3);
//
// 		}


// 		inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it,  Variable<double>& rVariable)
// 		{
// 			double& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);
//
// 			Aux_var = 0.0;
//
// 		}



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

    //ErosionUtils(void);

    //ErosionUtils(ErosionUtils& rSource);


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_NORMAL_TO_WALL_CALCULATIO_UTILS  defined */

