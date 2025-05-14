//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti.
//

#if !defined(KRATOS_AULF_UTILITIES_INCLUDED )
#define  KRATOS_AULF_UTILITIES_INCLUDED






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
#include "geometries/tetrahedra_3d_4.h"
//#include "pfem_2_application_variables.h"
#include "boost/smart_ptr.hpp"
#include "includes/cfd_variables.h"
#include "includes/deprecated_variables.h"

#include "processes/process.h"
#include "includes/condition.h"
#include "processes/find_nodal_neighbours_process.h"
#include "pfem_melting_application.h"
//#include "incompressible_fluid_application.h"



namespace Kratos
{

  template<std::size_t TDim> class Pfem2Utils
    {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(Pfem2Utils<TDim>);




//
//class Pfem2Utils
//{
//public:
//gggggggggg
 //   typedef Node<3> NodeType;
    //**********************************************************************************************
    //**********************************************************************************************

    //functions to apply the boundary conditions

    //**********************************************************************************************
    //**********************************************************************************************
    void MarkExcessivelyCloseNodes(ModelPart& ThisModelPart, const double admissible_distance_factor)
    {
    
         
        KRATOS_TRY;
        KRATOS_WATCH("ENTERD Mark close nodes")
        double fact2 = admissible_distance_factor*admissible_distance_factor;


        for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
        {
            if(in->FastGetSolutionStepValue(IS_INTERFACE) == 0) //if it is not a wall node i can erase
            {
                double hnode2 = in->FastGetSolutionStepValue(NODAL_H);
                hnode2 *= hnode2; //take the square

                //loop on neighbours and erase if they are too close
                for( GlobalPointersVector< Node >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                        i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
                {
                    if( bool(i->Is(TO_ERASE)) == false) //we can erase the current node only if the neighb is not to be erased
                    {
                        double dx = i->X() - in->X();
                        double dy = i->Y() - in->Y();
                        double dz = i->Z() - in->Z();

                        double dist2 = dx*dx + dy*dy + dz*dz;

                        if(dist2 < fact2 *  hnode2)
                            in->Set(TO_ERASE, true);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }
    
    
    //**********************************************************************************************
    //**********************************************************************************************

        void MarkNodesCloseToWall(ModelPart& ThisModelPart, double alpha_shape)
        //void MarkExcessivelyCloseNodes(ModelPart& ThisModelPart, const double admissible_distance_factor)
    {
   

        if (TDim==2)
        {

            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                int n_str=0;

                //counting number on nodes at the wall
                Geometry< Node >& geom = i->GetGeometry();
                n_str = int(geom[0].FastGetSolutionStepValue(IS_SOLID));
                n_str+= int(geom[1].FastGetSolutionStepValue(IS_SOLID));
                n_str+= int(geom[2].FastGetSolutionStepValue(IS_SOLID));
                //if two walls are at the wall, we check if the third node is close to it or not by passing the alpha-shape
                if (n_str==2)
                {
                    //KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl; 	
                    //if alpha shape tells to preserve
                    if (AlphaShape(alpha_shape, geom)==false)
                    {
                    	//KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl; 
                        for (int i=0; i<3; i++)
                        {
                            //if thats not the wall node, remove it
                            if (geom[i].FastGetSolutionStepValue(IS_SOLID)==0.0)
                            {
                                //KRATOS_ERROR << "AddTimeIntegratedLHS is not implemented." << std::endl; 
                            	 KRATOS_WATCH("AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII") 
                                geom[i].Set(TO_ERASE,true);
                                //KRATOS_WATCH("NODE CLOSE TO THE WALL - WILL BE ERASED!!!!")
                            }
                        }
                    }
                }

            }
        }
        if (TDim==3)
        {
            //KRATOS_WATCH("Inside mark nodes close to wall process")
            for(ModelPart::ElementsContainerType::iterator i = ThisModelPart.ElementsBegin();
                    i!=ThisModelPart.ElementsEnd(); i++)
            {
                int n_str=0;
                int n_fl=0;
                int n_lag=0;
                int n_interf=0;
                //counting number on nodes at the wall
                Geometry< Node >& geom = i->GetGeometry();
                for (unsigned int iii=0; iii<geom.size(); iii++)
                {
                    n_lag += int(geom[iii].FastGetSolutionStepValue(IS_LAGRANGIAN_INLET));
                    n_str += int(geom[iii].FastGetSolutionStepValue(IS_STRUCTURE));
                    n_fl += int(geom[iii].FastGetSolutionStepValue(IS_FLUID));
                    n_interf += int(geom[iii].FastGetSolutionStepValue(IS_INTERFACE));
                }


                //if three nodes are at the wall, we check if the fourth node is close to it or not by passing the alpha-shape
                //if (geom.size()==4.0 && n_str==3.0 && n_fl==4.0)
                if (geom.size()==4 && n_interf==3 && n_lag==0)
                {
                    //if alpha shape tells to preserve
                    if (AlphaShape3D(alpha_shape, geom)==false)
                    {
                        for (unsigned int iii=0; iii<geom.size(); iii++)
                        {
                            //if thats not the wall node, remove it
                            if (geom[iii].FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
                            {
                                geom[iii].Set(TO_ERASE,true);
                                KRATOS_WATCH("NODE CLOSE TO THE WALL - WILL BE ERASED!!!!")
                            }
                        }
                    }
                }


            }
            
          
        }
    }

    bool AlphaShape(double alpha_param, Geometry<Node >& pgeom)
    //bool AlphaShape(double alpha_param, Triangle2D<Node<3> >& pgeom)
    {
        KRATOS_TRY
        BoundedMatrix<double,2,2> J; //local jacobian
        BoundedMatrix<double,2,2> Jinv; //local jacobian
        static array_1d<double,2> c; //center pos
        static array_1d<double,2> rhs; //center pos

        double x0 = pgeom[0].X();
        double x1 = pgeom[1].X();
        double x2 = pgeom[2].X();

        double y0 = pgeom[0].Y();
        double y1 = pgeom[1].Y();
        double y2 = pgeom[2].Y();

        J(0,0)=2.0*(x1-x0);
        J(0,1)=2.0*(y1-y0);
        J(1,0)=2.0*(x2-x0);
        J(1,1)=2.0*(y2-y0);


        double detJ = J(0,0)*J(1,1)-J(0,1)*J(1,0);

        Jinv(0,0) =  J(1,1);
        Jinv(0,1) = -J(0,1);
        Jinv(1,0) = -J(1,0);
        Jinv(1,1) =  J(0,0);

        BoundedMatrix<double,2,2> check;


        if(detJ < 1e-12)
        {
            std::cout << "detJ = " << detJ << std::endl;
            ////mark as boundary
            //pgeom[0].GetSolutionStepValue(IS_INTERFACE) = 1;
            //pgeom[1].GetSolutionStepValue(IS_INTERFACE) = 1;
            //pgeom[2].GetSolutionStepValue(IS_INTERFACE) = 1;
            return false;
        }

        else
        {

            double x0_2 = x0*x0 + y0*y0;
            double x1_2 = x1*x1 + y1*y1;
            double x2_2 = x2*x2 + y2*y2;

            //finalizing the calculation of the inverted matrix
            //std::cout<<"MATR INV"<<MatrixInverse(msJ)<<std::endl;
            Jinv /= detJ;
            //calculating the RHS
            rhs[0] = (x1_2 - x0_2);
            rhs[1] = (x2_2 - x0_2);

            //calculate position of the center
            noalias(c) = prod(Jinv,rhs);

            double radius = sqrt(pow(c[0]-x0,2)+pow(c[1]-y0,2));

            //calculate average h
            double h;
            h =  pgeom[0].FastGetSolutionStepValue(NODAL_H);
            h += pgeom[1].FastGetSolutionStepValue(NODAL_H);
            h += pgeom[2].FastGetSolutionStepValue(NODAL_H);
            h *= 0.333333333;
            if (radius < h*alpha_param)
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        KRATOS_CATCH("")
    }
    bool AlphaShape3D( double alpha_param, Geometry<Node >& geom	)
    {
        KRATOS_TRY

        BoundedMatrix<double,3,3> J; //local jacobian
        BoundedMatrix<double,3,3> Jinv; //local jacobian
        array_1d<double,3> Rhs; //rhs
        array_1d<double,3> xc;
        double radius=0.0;

        const double x0 = geom[0].X();
        const double y0 = geom[0].Y();
        const double z0 = geom[0].Z();
        const double x1 = geom[1].X();
        const double y1 = geom[1].Y();
        const double z1 = geom[1].Z();
        const double x2 = geom[2].X();
        const double y2 = geom[2].Y();
        const double z2 = geom[2].Z();
        const double x3 = geom[3].X();
        const double y3 = geom[3].Y();
        const double z3 = geom[3].Z();

        //calculation of the jacobian
        J(0,0) = x1-x0;
        J(0,1) = y1-y0;
        J(0,2) = z1-z0;
        J(1,0) = x2-x0;
        J(1,1) = y2-y0;
        J(1,2) = z2-z0;
        J(2,0) = x3-x0;
        J(2,1) = y3-y0;
        J(2,2) = z3-z0;
// 			KRATOS_WATCH("33333333333333333333333");

        //inverse of the jacobian
        //first column
        Jinv(0,0) = J(1,1)*J(2,2) - J(1,2)*J(2,1);
        Jinv(1,0) = -J(1,0)*J(2,2) + J(1,2)*J(2,0);
        Jinv(2,0) = J(1,0)*J(2,1) - J(1,1)*J(2,0);
        //second column
        Jinv(0,1) = -J(0,1)*J(2,2) + J(0,2)*J(2,1);
        Jinv(1,1) = J(0,0)*J(2,2) - J(0,2)*J(2,0);
        Jinv(2,1) = -J(0,0)*J(2,1) + J(0,1)*J(2,0);
        //third column
        Jinv(0,2) = J(0,1)*J(1,2) - J(0,2)*J(1,1);
        Jinv(1,2) = -J(0,0)*J(1,2) + J(0,2)*J(1,0);
        Jinv(2,2) = J(0,0)*J(1,1) - J(0,1)*J(1,0);
        //calculation of determinant (of the input matrix)

// 			KRATOS_WATCH("44444444444444444444444444");
        double detJ = J(0,0)*Jinv(0,0)
                      + J(0,1)*Jinv(1,0)
                      + J(0,2)*Jinv(2,0);

        //volume = detJ * 0.16666666667;
// 			KRATOS_WATCH("55555555555555555555555");


        double x0_2 = x0*x0 + y0*y0 + z0*z0;
        double x1_2 = x1*x1 + y1*y1 + z1*z1;
        double x2_2 = x2*x2 + y2*y2 + z2*z2;
        double x3_2 = x3*x3 + y3*y3 + z3*z3;

        //finalizing the calculation of the inverted matrix
        Jinv /= detJ;

        //calculating the RHS
        //calculating the RHS
        Rhs[0] = 0.5*(x1_2 - x0_2);
        Rhs[1] = 0.5*(x2_2 - x0_2);
        Rhs[2] = 0.5*(x3_2 - x0_2);

        //calculate position of the center
        noalias(xc) = prod(Jinv,Rhs);
        //calculate radius
        radius = pow(xc[0] - x0,2);
        radius		  += pow(xc[1] - y0,2);
        radius		  += pow(xc[2] - z0,2);
        radius = sqrt(radius);

        double h;
        h =  geom[0].FastGetSolutionStepValue(NODAL_H);
        h += geom[1].FastGetSolutionStepValue(NODAL_H);
        h += geom[2].FastGetSolutionStepValue(NODAL_H);
        h += geom[3].FastGetSolutionStepValue(NODAL_H);
        h *= 0.250;

        if (radius < h*alpha_param)
        {
            return true;
        }
        else
        {
            return false;
        }

        KRATOS_CATCH("")
    }


    //**********************************************************************************************
    //**********************************************************************************************
    //imposes the velocity that corresponds to a
    void MoveLonelyNodes(ModelPart& ThisModelPart)
    {
        KRATOS_TRY;
        //KRATOS_WATCH("MOVING LONELY NODES")
        double Dt = ThisModelPart.GetProcessInfo()[DELTA_TIME];

        array_1d<double,3> DeltaDisp, acc;

        //const array_1d<double,3> body_force = ThisModelPart.ElementsBegin()->GetProperties()[BODY_FORCE];
        for(ModelPart::NodeIterator i = ThisModelPart.NodesBegin() ;
                i != ThisModelPart.NodesEnd() ; ++i)
        {
            if(
                (i)->FastGetSolutionStepValue(IS_INTERFACE) == 0 && //if it is not a wall node
                (i)->GetValue(NEIGHBOUR_ELEMENTS).size() == 0 &&//and it is lonely
                ((i)->GetDof(VELOCITY_X).IsFixed() == false || (i)->GetDof(VELOCITY_Y).IsFixed() == false || (i)->GetDof(VELOCITY_Z).IsFixed() == false) //and not the node of the wall, where the 0-displ is prescribed

            )
            {
                //i->Set(TO_ERASE,true);
                //set to zero the pressure
                (i)->FastGetSolutionStepValue(PRESSURE) = 0;

                const array_1d<double,3>& old_vel = (i)->FastGetSolutionStepValue(VELOCITY,1);
                array_1d<double,3>& vel = (i)->FastGetSolutionStepValue(VELOCITY);
                //array_1d<double,3>& acc = (i)->FastGetSolutionStepValue(ACCELERATION);

                noalias(acc) =  (i)->FastGetSolutionStepValue(BODY_FORCE);

                noalias(vel) = old_vel;
                noalias(vel) += Dt * acc ;

                //calculate displacements
/*                noalias(DeltaDisp) = Dt * vel;

                array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
                noalias(disp) = i->FastGetSolutionStepValue(DISPLACEMENT,1);
                noalias(disp) += DeltaDisp;

*/

            }

        }

        KRATOS_CATCH("")
    }
private:

    //aux vars
    static BoundedMatrix<double,3,3> msJ; //local jacobian
    static BoundedMatrix<double,3,3> msJinv; //inverse jacobian
    static array_1d<double,3> msc; //center pos
    static array_1d<double,3> ms_rhs; //center pos


};

}  // namespace Kratos.

#endif // KRATOS_AULF_UTILITIES_INCLUDED  defined
