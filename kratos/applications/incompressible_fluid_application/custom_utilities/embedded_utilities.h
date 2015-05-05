/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
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


//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pavel $
//   Date:                $Date: 2009-01-15 14:50:24 $
//   Revision:            $Revision: 1.12 $
//
//

#if !defined(KRATOS_EMBEDDED_UTILITIES_INCLUDED )
#define  KRATOS_EMBEDDED_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "utilities/math_utils.h"
#include "geometries/tetrahedra_3d_4.h"
#include "incompressible_fluid_application.h"

namespace Kratos
{
class EmbeddedUtils
{
public:
    typedef Node<3> NodeType;    
    ////////////////////////////////////////////////////////////////////////////////////////
    void CreateIntersConditions(ModelPart& model_part, ModelPart& interface_conditions_model_part)
    {
        KRATOS_TRY
	KRATOS_WATCH("CREATING THE INTERSECTION CONDITIONS")
        interface_conditions_model_part.Conditions().clear();
        interface_conditions_model_part.Elements().clear();
        interface_conditions_model_part.Nodes().clear();

        //reset the IS_INTERFACE flag - if the distance is negative (i.e. the fluid node lies inside of the solid object and is fictitious - set 1 otherwise 0
        for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin() ;  in != model_part.NodesEnd() ; ++in)
        {
            if (in->FastGetSolutionStepValue(DISTANCE)<0.0)
                in->FastGetSolutionStepValue(IS_INTERFACE)=1.0;
            else
                in->FastGetSolutionStepValue(IS_INTERFACE)=0.0;
        }
	//////////////////////////////////////////////////////////////////////////////////////////////////

       
        

        for(ModelPart::ElementsContainerType::iterator im = model_part.ElementsBegin() ;  im != model_part.ElementsEnd() ; ++im)
        {
            //intersection Points - at most we shall consider 4 points
            std::vector<array_1d<double,3> > IntersectionPoints;
	    array_1d<double,3> IntersectionVel(3, 0.0);
            //std::vector<array_1d<double,3> > IntersectionVel;

            IntersectionPoints.reserve(8);
            //IntersectionVel.reserve(4);

            array_1d<double,3> Point;
            //identify
            int intersection_count=0;
            //iterating over edges (01 02 03 12 13 23)
            for (int i=0; i<3; i++)
            {
                for (int j=i+1; j<4; j++)
                {
                    //std::cout<<"edge ij "<<i<<j<<std::endl;
                    double d0=im->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
                    double d1=im->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);

                    //if the product of distances of two nodes is negative - the edge is crossed
                    if (d0*d1<0.0)
                    {
                        //std::cout<<"Intersected edge "<<i<<j<<std::endl;

                        double x0=im->GetGeometry()[i].X();
                        double x1=im->GetGeometry()[j].X();
                        double k=0.0;
                        double b=0.0;
                        double x_inters=0.0;
                        double y_inters=0.0;
                        double z_inters=0.0;
                        if (x1!=x0)
                        {
                            k=(d1-d0)/(x1-x0);
                            b=d0-k*x0;
                            x_inters=-b/k;
                        }
                        else
                            x_inters=x0;

                        double y0=im->GetGeometry()[i].Y();
                        double y1=im->GetGeometry()[j].Y();
                        if (y1!=y0)
                        {
                            k=(d1-d0)/(y1-y0);
                            b=d0-k*y0;
                            y_inters=-b/k;
                        }
                        else
                            y_inters=y0;

                        double z0=im->GetGeometry()[i].Z();
                        double z1=im->GetGeometry()[j].Z();
                        if (z1!=z0)
                        {
                            k=(d1-d0)/(z1-z0);
                            b=d0-k*z0;
                            z_inters=-b/k;
                        }
                        else
                            z_inters=z0;

                        Point[0]=x_inters;
                        Point[1]=y_inters;
                        Point[2]=z_inters;

			//we assume that the velcoity is unique over the intersection.. so just one EMBEDDED_VELOCITY per element
			
			array_1d<double,3> vel;
			vel=im->GetValue(EMBEDDED_VELOCITY);
			//KRATOS_WATCH("THE CUT ELEMENT HAS A VELOCITY OF...")
			//KRATOS_WATCH(im->GetValue(EMBEDDED_VELOCITY))
                        //KRATOS_WATCH(Point)
			

                        IntersectionPoints.push_back(Point);
			IntersectionVel=vel;
			//IntersectionVel.push_back(vel);		
                        intersection_count++; 			
				
                    }

                }
            }
	    //KRATOS_WATCH(intersection_count)
	    IntersectionPoints.resize(intersection_count);
            //if the element is intersected by the embedded skin, create condition
            //if (intersection_count!=0)
	    if (intersection_count==3 || intersection_count==4)
            {
                //the size of array should represent actual intersection number
                //IntersectionPoints.resize(intersection_count);
		//IntersectionVel.resize(intersection_count);
                //KRATOS_WATCH(IntersectionPoints.size())


                //for now we assume zero velocity of the structure TO BE completed later...
                //array_1d<double,3> ZeroVel=ZeroVector(3);
                //for (unsigned int i=0; i<IntersectionVel.size(); i++)
                 //   IntersectionVel.push_back(ZeroVel);

                //////////////////////////////////////////////////////////////////////
                Geometry< Node<3> >::Pointer geom = im->pGetGeometry();
                Properties::Pointer properties = model_part.GetMesh().pGetProperties(1);

                int id=interface_conditions_model_part.Conditions().size()+1;

                if (IntersectionPoints.size()==3)
                {
                    Condition::Pointer p_condition(new ProjDirichletCond3D(id, geom,properties, IntersectionPoints[0], IntersectionPoints[1], IntersectionPoints[2], IntersectionVel));
                    interface_conditions_model_part.Conditions().push_back(p_condition);
                }
                else if (IntersectionPoints.size()==4)
                {
                    Condition::Pointer p_condition(new ProjDirichletCond3D(id, geom,properties, IntersectionPoints[0], IntersectionPoints[1], IntersectionPoints[2], IntersectionPoints[3], IntersectionVel));

                    interface_conditions_model_part.Conditions().push_back(p_condition);
                }
	    }
	    else if (intersection_count>4)
	    {
		//if there is a strangely large number of intersections we just treat this case as if there are just three...
		Geometry< Node<3> >::Pointer geom = im->pGetGeometry();
                Properties::Pointer properties = model_part.GetMesh().pGetProperties(1);             

		////////////////////////////////////////////////////////////////
		array_1d<double,3> AuxPoint=ZeroVector(3);
                array_1d<double,3> AuxPoint2=ZeroVector(3);
		array_1d<double,3> AuxPoint3=ZeroVector(3);
		array_1d<double,4> N;
                //array_1d<double,10> which_edge;
		Vector which_edge=ZeroVector(10);
		which_edge.resize(IntersectionPoints.size());

                for (unsigned int i=0; i<IntersectionPoints.size(); i++)
                {
                    which_edge[i]=100;
                }
                for (unsigned int kk=0; kk<IntersectionPoints.size(); kk++)
                {
                    AuxPoint=IntersectionPoints[kk];
                    CalculateN_at_Point(im->GetGeometry(), AuxPoint[0], AuxPoint[1], AuxPoint[2],  N);
                    for (int i=0; i<4; i++)
                    {
                        if (N[i]<0.00000000000001)
                            which_edge[kk]=i;
                    }
                }
                int which_edge_first;
                AuxPoint=IntersectionPoints[0];
                which_edge_first=which_edge[0];

                for (unsigned int kk=1; kk<IntersectionPoints.size(); kk++)
                {
                    if (which_edge[kk]!=which_edge_first)
                    {
                        AuxPoint2=IntersectionPoints[kk];
                    }
                }
		int which_edge_second=which_edge[1];
              
		for (unsigned int kk=2; kk<IntersectionPoints.size(); kk++)
                {
                    if (which_edge[kk]!=which_edge_first && which_edge[kk]!=which_edge_second)
                    {
                        AuxPoint3=IntersectionPoints[kk];
                    }
                }
                IntersectionPoints.resize(3);
                IntersectionPoints[0]=AuxPoint;
                IntersectionPoints[1]=AuxPoint2;
                IntersectionPoints[2]=AuxPoint3;

 		int id=interface_conditions_model_part.Conditions().size()+1;
	    	Condition::Pointer p_condition(new ProjDirichletCond3D(id, geom,properties, IntersectionPoints[0], IntersectionPoints[1], IntersectionPoints[2], IntersectionVel));
                    
		interface_conditions_model_part.Conditions().push_back(p_condition);


	    }	
            else if (intersection_count!=0)
		{
		KRATOS_WATCH("Strange number of intersections")
		KRATOS_WATCH(intersection_count)     
		//KRATOS_THROW_ERROR(std::logic_error,  "Strange number of intersections - neither 3 nor 4 - check  CreateIntersectionConditions function" , "");                         

		}
                   

            
            //end loop over elements
        }
        KRATOS_WATCH(interface_conditions_model_part)
        KRATOS_CATCH("")


    }


    ///////////////////////////////////////////////////////////////////////////
    ////////	SUBDOMAIN DISABLING			//////////////////
    ///////////////////////////////////////////////////////////////////////////
    void DisableSubdomain(ModelPart& full_model_part, ModelPart& reduced_model_part)
    {
        KRATOS_TRY
  KRATOS_THROW_ERROR(std::logic_error,  "USE THE PROCESS INSTEAD.... " , "");
       /*
        reduced_model_part.Conditions().clear();
        reduced_model_part.Elements().clear();
        reduced_model_part.Nodes().clear();

        reduced_model_part.Conditions().reserve(full_model_part.Conditions().size());
        reduced_model_part.Elements().reserve(full_model_part.Elements().size());
        reduced_model_part.Nodes().reserve(full_model_part.Nodes().size());

        int n_int=0;
        int n_disabled=0;

        for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ; in != full_model_part.NodesEnd() ; ++in)
        {
            in->FastGetSolutionStepValue(DISABLE)=false;
        }
        for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ; im != full_model_part.ElementsEnd() ; ++im)
        {
            n_int=im->GetGeometry()[0].FastGetSolutionStepValue(IS_INTERFACE);
            n_int+=im->GetGeometry()[1].FastGetSolutionStepValue(IS_INTERFACE);
            n_int+=im->GetGeometry()[2].FastGetSolutionStepValue(IS_INTERFACE);
	    n_int+=im->GetGeometry()[3].FastGetSolutionStepValue(IS_INTERFACE);

            if(n_int==4)
            {
                im->GetGeometry()[0].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
                im->GetGeometry()[1].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
                im->GetGeometry()[2].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
		im->GetGeometry()[3].FastGetSolutionStepValue(MATERIAL_VARIABLE)=true;
            }

            if (n_int<4)
            {
                reduced_model_part.Elements().push_back(*(im.base()));
                for (int i=0; i<3; i++)
                {
                    im->GetGeometry()[i].FastGetSolutionStepValue(DISABLE)=true;
                }

            }

            if (n_int>4)
                KRATOS_THROW_ERROR(std::logic_error,  "Number of DISABLE flags cant exceed number of the element nodes.... " , "");
            
        }
        for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ;
                in != full_model_part.NodesEnd() ; ++in)
        {
            n_disabled=in->FastGetSolutionStepValue(DISABLE);
            if (n_disabled==1.0)
            {
                reduced_model_part.Nodes().push_back(*(in.base()));
            }
        }

*/

        KRATOS_CATCH("")
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////
    void ApplyProjDirichlet(ModelPart& full_model_part)
    {
        KRATOS_TRY
	KRATOS_THROW_ERROR(std::logic_error,  "USE THE PROCESS INSTEAD.... " , "");
	/*
        unsigned int n_old_int;

        //first we remove the Dirichlet conditions from the nodes that were defining the interface in the previous step:
        for(ModelPart::NodesContainerType::iterator in = full_model_part.NodesBegin() ;
                in != full_model_part.NodesEnd() ; ++in)
        {
            n_old_int=in->FastGetSolutionStepValue(IS_INTERFACE,1);

            //make the velocity free at the nodes that are not the "Dangerosu ones"
            if (n_old_int>0.0 && in->FastGetSolutionStepValue(IS_INTERFACE)!=100.0)
            {
                //KRATOS_WATCH("OLD INTERFACEEEEEEEEEEEEEEEE!!!!!!!!!!!!!!!!!!!!!!1")
                in->FastGetSolutionStepValue(DISABLE)=false;
                in->Free(VELOCITY_X);
                in->Free(VELOCITY_Y);
                in->Free(VELOCITY_Z);

                in->Free(AUX_VEL_X);
                in->Free(AUX_VEL_Y);
                in->Free(AUX_VEL_Z);

                in->Free(PRESSURE);
            }
            //fix the IS_INt to 1 additionally at the nodes that are too close to the interface
            if (in->FastGetSolutionStepValue(IS_INTERFACE)==100.0)
            {
                in->FastGetSolutionStepValue(IS_INTERFACE)=1.0;
                KRATOS_WATCH("BAD NODE IS")
                KRATOS_WATCH(in->GetId())
            }

        }

        for(ModelPart::ElementsContainerType::iterator im = full_model_part.ElementsBegin() ;
                im != full_model_part.ElementsEnd() ; ++im)
        {
            double n_int=0.0;
            //iterate over the of nodes in the element (interface element)
            for (unsigned int i=0; i<im->GetGeometry().size(); i++)
            {
                n_int=im->GetGeometry()[i].FastGetSolutionStepValue(IS_INTERFACE);
                // if the node is lying on fictitious side - apply the project Dirichlet condition
                if (n_int==1.0)
                    //if (((ic->GetGeometry()[i]).GetDof(AUX_VEL_X)).IsFixed())
                {
                    im->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)=im->GetGeometry()[i].FastGetSolutionStepValue(AUX_VEL);
                    im->GetGeometry()[i].Fix(VELOCITY_X);
                    im->GetGeometry()[i].Fix(VELOCITY_Y);
                    im->GetGeometry()[i].Fix(VELOCITY_Z);
                }
            }

        }


	*/

        KRATOS_CATCH("")
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //													   //
    //			AUXILIARY FUNCTIONS								   //
    //													   //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
double CalculateVol(	const double x0, const double y0, const double z0,
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

    double detJ = (x10 * y20 * z30 - x10 * y30 * z20) + (y10 * z20 * x30 - y10 * x20 * z30) + (z10 * x20 * y30 - z10 * y20 * x30);
    return  detJ*0.1666666666666666666667;

}
    /////////////////////////////////////////////////////////////
inline void CalculateN_at_Point(Element::GeometryType& geom, const double xc, const double yc, const double zc, array_1d<double,4>& N_at_c)
{
    double x0=geom[0].X();
    double x1=geom[1].X();
    double x2=geom[2].X();
    double x3=geom[3].X();

    double y0=geom[0].Y();
    double y1=geom[1].Y();
    double y2=geom[2].Y();
    double y3=geom[3].Y();

    double z0=geom[0].Z();
    double z1=geom[1].Z();
    double z2=geom[2].Z();
    double z3=geom[3].Z();

    double vol = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

    double inv_vol = 0.0;
    if(vol < 0.000000000000001)
    {
        KRATOS_THROW_ERROR(std::logic_error,"element with zero vol found","");
    }
    else
    {
        inv_vol = 1.0 / vol;
    }

    N_at_c[0] = CalculateVol(x1,y1,z1,x3,y3,z3,x2,y2,z2,xc,yc,zc) * inv_vol;
    N_at_c[1] = CalculateVol(x3,y3,z3,x0,y0,z0,x2,y2,z2,xc,yc,zc) * inv_vol;
    N_at_c[2] = CalculateVol(x3,y3,z3,x1,y1,z1,x0,y0,z0,xc,yc,zc) * inv_vol;
    N_at_c[3] = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) * inv_vol;

    //if the intersection is close one of the tetrahedra's vertices ->send this message
    /*
    if (  (N_at_c[0]<0.01 && N_at_c[1]<0.01 && N_at_c[2]<0.01) ||
            (N_at_c[0]<0.01 && N_at_c[2]<0.01 && N_at_c[3]<0.01) ||
            (N_at_c[2]<0.01 && N_at_c[1]<0.01 && N_at_c[3]<0.01) ||
            (N_at_c[0]<0.01 && N_at_c[1]<0.01 && N_at_c[3]<0.01)     )
        KRATOS_WATCH("Dangerous VERTICES!!! Intersection is very close to the node")
        //KRATOS_THROW_ERROR(std::logic_error,  "Too close to the node is the INTERSECTION!!!! " , "")
	*/

    }
/*
    inline void CalculateN_at_Point(Element::GeometryType& geom, const double xc, const double yc, const double zc, array_1d<double,4>& N_at_c)
    {

        //first we calculate the volume of the whole tetrahedron
	double Tot_Vol=GeometryUtils::CalculateVolume3D(geom);


	//xc, 1 , 2,  3
	double x1c = geom[1].X() - xc;
        double y1c = geom[1].Y() - yc;
        double z1c = geom[1].Z() - zc;

        double x2c = geom[2].X() - xc;
        double y2c = geom[2].Y() - yc;
        double z2c = geom[2].Z() - zc;

        double x3c = geom[3].X() - xc;
        double y3c = geom[3].Y() - yc;
        double z3c = geom[3].Z() - zc;

        double detJ = x1c * y2c * z3c - x1c * y3c * z2c + y1c * z2c * x3c - y1c * x2c * z3c + z1c * x2c * y3c - z1c * y2c * x3c;
	double Vol0= detJ*0.1666666666666666666667;


	//////////////////        
        // xc, 0, 2, 3
        x1c = geom[0].X() - xc ;
        y1c = geom[0].Y() - yc ;
        z1c = geom[0].Z() - zc ;

        x2c = geom[2].X() - xc;
        y2c = geom[2].Y() - yc;
        z2c = geom[2].Z() - zc;

        x3c = geom[3].X() - xc;
        y3c = geom[3].Y() - yc;
        z3c = geom[3].Z() - zc;
	detJ = x1c * y2c * z3c - x1c * y3c * z2c + y1c * z2c * x3c - y1c * x2c * z3c + z1c * x2c * y3c - z1c * y2c * x3c;
	double Vol1= detJ*0.1666666666666666666667;

	//////////////////        
        // xc, 0, 1, 3
	x1c = geom[0].X() - xc ;
        y1c = geom[0].Y() - yc ;
        z1c = geom[0].Z() - zc ;

        x2c = geom[1].X() - xc;
        y2c = geom[1].Y() - yc;
        z2c = geom[1].Z() - zc;

        x3c = geom[3].X() - xc;
        y3c = geom[3].Y() - yc;
        z3c = geom[3].Z() - zc;
	detJ = x1c * y2c * z3c - x1c * y3c * z2c + y1c * z2c * x3c - y1c * x2c * z3c + z1c * x2c * y3c - z1c * y2c * x3c;
	double Vol2= detJ*0.1666666666666666666667;


	//////////////////        
        // xc, 0, 1, 2
	x1c = geom[0].X() - xc ;
        y1c = geom[0].Y() - yc ;
        z1c = geom[0].Z() - zc ;

        x1c = geom[1].X() - xc;
        y1c = geom[1].Y() - yc;
        z1c = geom[1].Z() - zc;

        x2c = geom[2].X() - xc;
        y2c = geom[2].Y() - yc;
        z2c = geom[2].Z() - zc;
	detJ = x1c * y2c * z3c - x1c * y3c * z2c + y1c * z2c * x3c - y1c * x2c * z3c + z1c * x2c * y3c - z1c * y2c * x3c;
	double Vol3= detJ*0.1666666666666666666667;
	
        
        if (Tot_Vol<0.00000000000000000001)
            KRATOS_THROW_ERROR(std::logic_error,  "Your element Proj DIrichlet Cond has a zero volume!!!! " , "");
        //and now we fill in the array of shape functions values:
        // 1 0 2
        N_at_c[0]=fabs(Vol0/Tot_Vol);
        N_at_c[1]=fabs(Vol1/Tot_Vol);
        N_at_c[2]=fabs(Vol2/Tot_Vol);
	N_at_c[3]=fabs(Vol3/Tot_Vol);
        if (  (N_at_c[0]<0.05 && N_at_c[1]<0.05) || (N_at_c[0]<0.05 && N_at_c[2]<0.05) || (N_at_c[2]<0.05 && N_at_c[1]<0.05))
            KRATOS_WATCH("Dangerous VERTICES!!!")
            //KRATOS_THROW_ERROR(std::logic_error,  "Too close to the node is the INTERSECTION!!!! " , "")
	KRATOS_WATCH(N_at_c)
	
        }
    //////////////////////////////////////////////////////////////////////////////////////////
   
*/

private:




};

}  // namespace Kratos.

#endif // KRATOS_EMBEDDED_UTILITIES_INCLUDED  defined 


