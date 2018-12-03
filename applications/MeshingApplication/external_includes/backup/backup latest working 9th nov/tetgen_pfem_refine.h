//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-15 14:50:34 $
//   Revision:            $Revision: 1.8 $
//


 

#if !defined(KRATOS_TETGEN_PFEM_MODELER_H_INCLUDED )
#define  KRATOS_TETGEN_PFEM_MODELER_H_INCLUDED
  


// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>
#include <boost/timer.hpp>



#include "tetgen.h" // Defined tetgenio, tetrahedralize().

// Project includes
#include "includes/define.h"
#include "utilities/geometry_utilities.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"


#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "meshing_application.h"
#include "processes/node_erase_process.h"
#include "includes/deprecated_variables.h"

#include "spatial_containers/spatial_containers.h"
//#include "containers/bucket.h"
//#include "containers/kd_tree.h"
//#include "external_includes/trigen_refine.h"
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
	class TetGenPfemModeler  
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TetGenPfemModeler
		KRATOS_CLASS_POINTER_DEFINITION(TetGenPfemModeler);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
	    TetGenPfemModeler() :
		mJ(ZeroMatrix(3,3)), //local jacobian
		mJinv(ZeroMatrix(3,3)), //inverse jacobian
		mc(ZeroVector(3)), //dimension = number of nodes
		mRhs(ZeroVector(3)){} //dimension = number of nodes

		/// Destructor.
		virtual ~TetGenPfemModeler()= default;


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{


		//*******************************************************************************************
		//*******************************************************************************************
		void ReGenerateMesh(
			ModelPart& ThisModelPart , 
			Element const& rReferenceElement, 
			Condition const& rReferenceBoundaryCondition,
			NodeEraseProcess& node_erase, bool rem_nodes = true, bool add_nodes=true,
			double alpha_param = 1.4, double h_factor=0.5)
		{

			KRATOS_TRY
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_FREE_SURFACE)==false )
				KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_FREE_SURFACE---- variable!!!!!! ERROR","");
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_STRUCTURE)==false )
				KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_STRUCTURE---- variable!!!!!! ERROR","");
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_BOUNDARY)==false )
				KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_BOUNDARY---- variable!!!!!! ERROR","");
			if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_FLUID)==false )
				KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_FLUID---- variable!!!!!! ERROR","");

			KRATOS_WATCH(" ENTERED TETGENMESHSUITE PFEM of Meshing Application")
			
			//clearing elements
                        
                        int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();

			ThisModelPart.Elements().clear();
// 			ThisModelPart.Conditions().clear();
                        

                        
//20151124 ajarauta -> TO FINISH!!! INSERT NODE AT CENTER OF FACE IF TRIANGLE IS TOO BIG
			int id = (ThisModelPart.Nodes().end() - 1)->Id() + 1;
			
                        
                        //n_struc==ic->GetGeometry().size())
                        
                        
                        
                        double mesh_size = 0.0;
                        for(ModelPart::ConditionsContainerType::iterator ic = ThisModelPart.ConditionsBegin() ;
			    ic != ThisModelPart.ConditionsEnd() ; ic++)
			    {
			      if (ic->GetGeometry().size()==3)
			      {
                                  
                                unsigned int n_st = ic->GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);
				n_st		 += ic->GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);
				n_st		 += ic->GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);
                                
				unsigned int n_tp = ic->GetGeometry()[0].FastGetSolutionStepValue(TRIPLE_POINT);
				n_tp		 += ic->GetGeometry()[1].FastGetSolutionStepValue(TRIPLE_POINT);
				n_tp		 += ic->GetGeometry()[2].FastGetSolutionStepValue(TRIPLE_POINT);
                                
                                unsigned int n_fs = ic->GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE);
				n_fs		 += ic->GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE);
				n_fs		 += ic->GetGeometry()[2].FastGetSolutionStepValue(IS_FREE_SURFACE);
                                
				/////// 1
				if ((n_st) == 3 && (n_tp == 0))
				{
                                    mesh_size  = ic->GetGeometry()[0].FastGetSolutionStepValue(NODAL_H);
                                    mesh_size += ic->GetGeometry()[1].FastGetSolutionStepValue(NODAL_H);
                                    mesh_size += ic->GetGeometry()[2].FastGetSolutionStepValue(NODAL_H);
                                    mesh_size *= 0.33333333333;
                                  
                                    int ii = 0;
                                    int jj = 1;
                                    int kk = 2;
                                    // Node with index kk is the one that is not on the contact line
                                    double xi=ic->GetGeometry()[ii].X();
                                    double yi=ic->GetGeometry()[ii].Y();
                                    double zi=ic->GetGeometry()[ii].Z();
                                    double xj=ic->GetGeometry()[jj].X();
                                    double yj=ic->GetGeometry()[jj].Y();
                                    double zj=ic->GetGeometry()[jj].Z();
                                    double xk=ic->GetGeometry()[kk].X();
                                    double yk=ic->GetGeometry()[kk].Y();
                                    double zk=ic->GetGeometry()[kk].Z();
                                    
                                    array_1d<double,3> rij = Vector3D(xi,yi,zi,xj,yj,zj);
                                    array_1d<double,3> rik = Vector3D(xi,yi,zi,xk,yk,zk);
                                    array_1d<double,3> rjk = Vector3D(xj,yj,zj,xk,yk,zk);
                                    
                                    double dist_ij = Norm3D(rij);
                                    double dist_ik = Norm3D(rik);
                                    double dist_jk = Norm3D(rjk);

                                    //if the distance between IS_FREE_SURFACE node and TRIPLE_POINT node is too long, we add a node in the middle of the face
                                    Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();
                                
                                    double factor = 3.0;	
                                    if ((dist_ij > factor*mesh_size) || (dist_ik > factor*mesh_size) || (dist_jk > factor*mesh_size))
                                    {
                                        id++;
                                        double x = 0.3333333333333*(xi+xj+xk);
                                        double y = 0.3333333333333*(yi+yj+yk);
                                        double z = 0.3333333333333*(zi+zj+zk);
                            
                                        Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

                                        //putting the new node also in an auxiliary list
                                        KRATOS_WATCH("adding nodes to list")										

                                        std::cout << "new node id = " << pnode->Id() << std::endl;
                                        //generating the dofs
                                        for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
                                        {
                                        Node<3>::DofType& rDof = *iii;
                                        Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
                                        (p_new_dof)->FreeDof();
                                        }
                                        Geometry<Node<3> >& geom = ic->GetGeometry();

                                        InterpolateOnFace(geom, ii, jj, kk, step_data_size, pnode);
                                        const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
                                        pnode->X0() = pnode->X() - disp[0];
                                        pnode->Y0() = pnode->Y() - disp[1];
                                        pnode->Z0() = pnode->Z() - disp[2];
                                        KRATOS_WATCH("Added node at the FACE with all IS_STRUCTURE, no TRIPLE_POINT")
                                        
//                                         n_struc==ic->GetGeometry().size()

                    // 			
                                    }
                                    
                                 
                                    if ((dist_ik < mesh_size * 0.333 ) || (dist_jk < mesh_size * 0.333) || (dist_ij < mesh_size * 0.333))
                                    {
//                                         ic->GetGeometry()[kk].Set(TO_ERASE,true);
                                        ic->GetGeometry()[ii].Set(TO_ERASE,true);
//                                          ic->GetGeometry()[jj].Set(TO_ERASE,true);
                                        KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF IS_STRUCTURE IS_STRUCTURE ELAF kkkkk")
                                    }
                                
//                                     if ((dist_ik < mesh_size *  0.67) || (dist_ij < mesh_size *  0.67))
//                                     {
// //                                          ic->GetGeometry()[kk].Set(TO_ERASE,true);
//                                         ic->GetGeometry()[ii].Set(TO_ERASE,true);
// //                                          ic->GetGeometry()[jj].Set(TO_ERASE,true);
//                                         KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF IS_STRUCTURE IS_STRUCTURE ELAF iiii")
//                                     }
//                                     
//                                     if ((dist_jk < mesh_size *  0.67 ) || (dist_ij < mesh_size * 0.67))
//                                     {
// //                                             ic->GetGeometry()[kk].Set(TO_ERASE,true);
// //                                          ic->GetGeometry()[ii].Set(TO_ERASE,true);
//                                         ic->GetGeometry()[jj].Set(TO_ERASE,true);
//                                         KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF IS_STRUCTURE IS_STRUCTURE ELAF jjjjj")
//                                     }
                                        
                                    
//                                     
 				}
                                /////// 2
                                
                                if (((n_st) == 3 && (n_tp == 2)))
				{
                                        mesh_size  = ic->GetGeometry()[0].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size += ic->GetGeometry()[1].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size += ic->GetGeometry()[2].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size *= 0.33333333333;
                                        
                                        int ii = 0;
                                        int jj = 1;
                                        int kk = 2;
                                        if (ic->GetGeometry()[0].FastGetSolutionStepValue(TRIPLE_POINT) == 0.0)
                                        {
                                            ii = 1;
                                            jj = 2;
                                            kk = 0;
                                        }
                                        if (ic->GetGeometry()[1].FastGetSolutionStepValue(TRIPLE_POINT) == 0.0)
                                        {
                                            ii = 0;
                                            jj = 2;
                                            kk = 1;
                                        }				    
                                        double xi=ic->GetGeometry()[ii].X();
                                        double yi=ic->GetGeometry()[ii].Y();
                                        double zi=ic->GetGeometry()[ii].Z();
                                        double xj=ic->GetGeometry()[jj].X();
                                        double yj=ic->GetGeometry()[jj].Y();
                                        double zj=ic->GetGeometry()[jj].Z();
                                        double xk=ic->GetGeometry()[kk].X();
                                        double yk=ic->GetGeometry()[kk].Y();
                                        double zk=ic->GetGeometry()[kk].Z();
                                        
                                        array_1d<double,3> rij = Vector3D(xi,yi,zi,xj,yj,zj);
                                        array_1d<double,3> rik = Vector3D(xi,yi,zi,xk,yk,zk);
                                        array_1d<double,3> rjk = Vector3D(xj,yj,zj,xk,yk,zk);
                                        
                                        double dist_ij = Norm3D(rij);
                                        double dist_ik = Norm3D(rik);
                                        double dist_jk = Norm3D(rjk);

                                        //if the distance between IS_FREE_SURFACE node and TRIPLE_POINT node is too long, we add a node in the middle of the face
                                        Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();
                                    
                                        double factor = 2.0;	
                                        if ((dist_ij > factor*mesh_size) || (dist_ik > factor*mesh_size) || (dist_jk > factor*mesh_size))
                                        {
                                            id++;
                                            double x = 0.3333333333333*(xi+xj+xk);
                                            double y = 0.3333333333333*(yi+yj+yk);
                                            double z = 0.3333333333333*(zi+zj+zk);
                                
                                            Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

                                            //putting the new node also in an auxiliary list
                                            KRATOS_WATCH("adding nodes to list")										

                                            std::cout << "new node id = " << pnode->Id() << std::endl;
                                            //generating the dofs
                                            for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
                                            {
                                            Node<3>::DofType& rDof = *iii;
                                            Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
                                            (p_new_dof)->FreeDof();
                                            }
                                            Geometry<Node<3> >& geom = ic->GetGeometry();

                                            InterpolateOnFace(geom, ii, jj, kk, step_data_size, pnode);
                                            const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
                                            pnode->X0() = pnode->X() - disp[0];
                                            pnode->Y0() = pnode->Y() - disp[1];
                                            pnode->Z0() = pnode->Z() - disp[2];
                                            KRATOS_WATCH("Added node at the EDGE with all IS_STRUCTURE, two TRIPLE_POINT")
                                    }
                                    
//                                         if ((dist_ik < mesh_size * 0.5 ) || (dist_jk < mesh_size * 0.5))
//                                         {
//                                         ic->GetGeometry()[kk].Set(TO_ERASE,true);
//                                         KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF two_TP two_TP ELAF kkkkk")
//                                         }
                                
                                        if ((dist_ik < mesh_size * 0.5 ) || (dist_jk < mesh_size * 0.5) || (dist_ij < mesh_size * 0.5))
                                        {
                                            if ( ic->GetGeometry()[ii].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
                                            {
                                                ic->GetGeometry()[ii].Set(TO_ERASE,true);
                                                KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF two_TP  ELAF iiii")
                                            }
                                                
                                            else if ( ic->GetGeometry()[jj].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
                                            {
                                                ic->GetGeometry()[jj].Set(TO_ERASE,true);
                                                KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF two_TP  ELAF jjjj")
                                            }
                                                
                                            else if ( ic->GetGeometry()[kk].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
                                            {
                                                ic->GetGeometry()[kk].Set(TO_ERASE,true);
                                                KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF two_TP  ELAF kkkk")
                                            }
                                        }




                                }
                                
                                ////// 3
                                
                                if ((n_st) == 3 && (n_tp == 1))
				{
                                    mesh_size  = ic->GetGeometry()[0].FastGetSolutionStepValue(NODAL_H);
                                    mesh_size += ic->GetGeometry()[1].FastGetSolutionStepValue(NODAL_H);
                                    mesh_size += ic->GetGeometry()[2].FastGetSolutionStepValue(NODAL_H);
                                    mesh_size *= 0.33333333333;
                                    
                                    int ii = 0;
				    int jj = 1;
				    int kk = 2;
				    if (ic->GetGeometry()[1].FastGetSolutionStepValue(TRIPLE_POINT) == 1.0)
				    {
				        ii = 1;
					jj = 0;
					kk = 2;
				    }
				    if (ic->GetGeometry()[2].FastGetSolutionStepValue(TRIPLE_POINT) == 1.0)
				    {
				        ii = 2;
					jj = 0;
					kk = 1;
				    }				    
				    // Node with index ii is the TRIPLE_POINT
                                    double xi=ic->GetGeometry()[ii].X();
                                    double yi=ic->GetGeometry()[ii].Y();
                                    double zi=ic->GetGeometry()[ii].Z();
                                    double xj=ic->GetGeometry()[jj].X();
                                    double yj=ic->GetGeometry()[jj].Y();
                                    double zj=ic->GetGeometry()[jj].Z();
                                    double xk=ic->GetGeometry()[kk].X();
                                    double yk=ic->GetGeometry()[kk].Y();
                                    double zk=ic->GetGeometry()[kk].Z();
                                    
                                    array_1d<double,3> rij = Vector3D(xi,yi,zi,xj,yj,zj);
                                    array_1d<double,3> rik = Vector3D(xi,yi,zi,xk,yk,zk);
                                    array_1d<double,3> rjk = Vector3D(xj,yj,zj,xk,yk,zk);
                                    
                                    double dist_ij = Norm3D(rij);
                                    double dist_ik = Norm3D(rik);
                                    double dist_jk = Norm3D(rjk);

                                    //if the distance between IS_FREE_SURFACE node and TRIPLE_POINT node is too long, we add a node in the middle of the face
                                    Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();
                                
                                    double factor = 2.1;	
                                    if ((dist_ij > factor*mesh_size) || (dist_ik > factor*mesh_size) || (dist_jk > factor*mesh_size))
                                    {
                                        id++;
                                        double x = 0.3333333333333*(xi+xj+xk);
                                        double y = 0.3333333333333*(yi+yj+yk);
                                        double z = 0.3333333333333*(zi+zj+zk);
                            
                                        Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

                                        //putting the new node also in an auxiliary list
                                        KRATOS_WATCH("adding nodes to list")										

                                        std::cout << "new node id = " << pnode->Id() << std::endl;
                                        //generating the dofs
                                        for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
                                        {
                                        Node<3>::DofType& rDof = *iii;
                                        Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
                                        (p_new_dof)->FreeDof();
                                        }
                                        Geometry<Node<3> >& geom = ic->GetGeometry();

                                        InterpolateOnFace(geom, ii, jj, kk, step_data_size, pnode);
                                        const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
                                        pnode->X0() = pnode->X() - disp[0];
                                        pnode->Y0() = pnode->Y() - disp[1];
                                        pnode->Z0() = pnode->Z() - disp[2];
                                        KRATOS_WATCH("Added node at the FACE with all IS_STRUCTURE, one TRIPLE_POINT")
                                        
//                                         double factor_2 = 0.75;			
//                                         if ( (dist_ki < factor_2*nodal_h) || (dist_kj < factor_2*nodal_h) )
//                                             ic->GetGeometry()[kk].Set(TO_ERASE,true);
                                        
				    }
				    
				    
//                                     if ((dist_ik < mesh_size * 0.5 ) || (dist_jk < mesh_size * 0.5))
//                                     {
//                                         ic->GetGeometry()[kk].Set(TO_ERASE,true);
// //                                         ic->GetGeometry()[jj].Set(TO_ERASE,true);
// //                                          ic->GetGeometry()[jj].Set(TO_ERASE,true);
//                                         KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF one_tp one_tp ELAF kkkkk")
//                                     }
//                                     
//                                     if ((dist_jk < mesh_size * 0.5 ) || (dist_ij < mesh_size * 0.5))
//                                     {
// //                                             ic->GetGeometry()[kk].Set(TO_ERASE,true);
// //                                          ic->GetGeometry()[ii].Set(TO_ERASE,true);
//                                         ic->GetGeometry()[jj].Set(TO_ERASE,true);
//                                         KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF one_tp one_tp ELAF jjjjj")
//                                     }


                                    if ((dist_ik < mesh_size * 0.5 ) || (dist_jk < mesh_size * 0.5) || (dist_ij < mesh_size * 0.5))
                                    {
                                        if ( ic->GetGeometry()[ii].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
                                        {
                                            ic->GetGeometry()[ii].Set(TO_ERASE,true);
                                            KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF one_tp one_tp ELAF iiii")
                                        }
                                            
//                                         else if ( ic->GetGeometry()[jj].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
//                                         {
//                                             ic->GetGeometry()[jj].Set(TO_ERASE,true);
//                                             KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF one_tp one_tp ELAF jjjj")
//                                         }
                                            
                                        else if ( ic->GetGeometry()[kk].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
                                        {
                                            ic->GetGeometry()[kk].Set(TO_ERASE,true);
                                            KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF one_tp one_tp ELAF kkkk")
                                        }
                                        }

				    
                                }
                                /////// 4
                                
                                if ((n_fs == 1) && (n_tp == 2) || (n_fs == 2) && (n_tp == 1) )
                                {
                                        mesh_size  = ic->GetGeometry()[0].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size += ic->GetGeometry()[1].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size += ic->GetGeometry()[2].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size *= 0.33333333333;
                                        int ii = 0;
                                        int jj = 1;
                                        int kk = 2;
                                        if (ic->GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
                                        {
                                            ii = 1;
                                            jj = 2;
                                            kk = 0;
                                        }
                                        if (ic->GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
                                        {
                                            ii = 0;
                                            jj = 2;
                                            kk = 1;
                                        }
                                        // Node with index kk is the one that is not on the contact line
                                        double xi=ic->GetGeometry()[ii].X();
                                        double yi=ic->GetGeometry()[ii].Y();
                                        double zi=ic->GetGeometry()[ii].Z();
                                        double xj=ic->GetGeometry()[jj].X();
                                        double yj=ic->GetGeometry()[jj].Y();
                                        double zj=ic->GetGeometry()[jj].Z();
                                        double xk=ic->GetGeometry()[kk].X();
                                        double yk=ic->GetGeometry()[kk].Y();
                                        double zk=ic->GetGeometry()[kk].Z();
                                        
                                        
                                        array_1d<double,3> rij = Vector3D(xi,yi,zi,xj,yj,zj);
                                        array_1d<double,3> rik = Vector3D(xi,yi,zi,xk,yk,zk);
                                        array_1d<double,3> rjk = Vector3D(xj,yj,zj,xk,yk,zk);
                                        
                                        double dist_ij = Norm3D(rij);
                                        double dist_ik = Norm3D(rik);
                                        double dist_jk = Norm3D(rjk);
                                        
                                        //if the distance between IS_FREE_SURFACE node and TRIPLE_POINT node is too long, we add a node in the middle of the face
                                        Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();
                                    
                                        double factor = 1.5;	
                                        if ((dist_ik > factor*mesh_size) || (dist_jk > factor*mesh_size) || (dist_ij > factor*mesh_size))
                                        {
                                        id++;
                                        double x = 0.3333333333333*(xi+xj+xk);
                                        double y = 0.3333333333333*(yi+yj+yk);
                                        double z = 0.3333333333333*(zi+zj+zk);
                                
                                        Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

                                        //putting the new node also in an auxiliary list
                                        KRATOS_WATCH("adding nodes to list")										

                                        std::cout << "new node id = " << pnode->Id() << std::endl;
                                        //generating the dofs
                                        for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
                                        {
                                            Node<3>::DofType& rDof = *iii;
                                            Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
                                            (p_new_dof)->FreeDof();
                                        }
                                        Geometry<Node<3> >& geom = ic->GetGeometry();

    // 				      InterpolateOnEdge(geom, ii, kk, step_data_size, pnode);
                                        InterpolateOnFace(geom, ii, jj, kk, step_data_size, pnode);
                                        const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
                                        pnode->X0() = pnode->X() - disp[0];
                                        pnode->Y0() = pnode->Y() - disp[1];
                                        pnode->Z0() = pnode->Z() - disp[2];
                                        KRATOS_WATCH("Added node at the FACE with IS_FREE_SURFACE")
                                        }
                                        

                                        
//                                         if ((dist_ik < mesh_size * 0.5 ) || (dist_jk < mesh_size * 0.5 ) || (dist_ij < mesh_size * 0.5))
//                                             {
//                                                 ic->GetGeometry()[kk].Set(TO_ERASE,true);
// //                                              ic->GetGeometry()[ii].Set(TO_ERASE,true);
// //                                                 ic->GetGeometry()[jj].Set(TO_ERASE,true);
//                                                 KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF FACE FACE ELAF kkkkkk")
// //                                                 KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF FACE FACE ELAF jjjjj")
//                                             }

                                        if ((dist_ik < mesh_size * 0.5 ) || (dist_jk < mesh_size * 0.5) || (dist_ij < mesh_size * 0.5))
                                        {
                                            if ( ic->GetGeometry()[ii].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
                                            {
                                                ic->GetGeometry()[ii].Set(TO_ERASE,true);
                                                KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF one_tp one_tp ELAF iiii")
                                            }
                                                
                                            else if ( ic->GetGeometry()[jj].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
                                            {
                                                ic->GetGeometry()[jj].Set(TO_ERASE,true);
                                                KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF one_tp one_tp ELAF jjjj")
                                            }
                                                
                                            else if ( ic->GetGeometry()[kk].FastGetSolutionStepValue(TRIPLE_POINT) < 1e-15)
                                            {
                                                ic->GetGeometry()[kk].Set(TO_ERASE,true);
                                                KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF one_tp one_tp ELAF kkkk")
                                            }
                                        }
                                            
                                        
				}
                                
                                
                         /////// 5
                                
                                if ((n_fs == 3) && (n_tp == 0) )
                                {
                                        mesh_size  = ic->GetGeometry()[0].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size += ic->GetGeometry()[1].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size += ic->GetGeometry()[2].FastGetSolutionStepValue(NODAL_H);
                                        mesh_size *= 0.33333333333;
                                        int ii = 0;
                                        int jj = 1;
                                        int kk = 2;
                                        if (ic->GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
                                        {
                                            ii = 1;
                                            jj = 2;
                                            kk = 0;
                                        }
                                        if (ic->GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE) == 1.0)
                                        {
                                            ii = 0;
                                            jj = 2;
                                            kk = 1;
                                        }
                                        // Node with index kk is the one that is not on the contact line
                                        double xi=ic->GetGeometry()[ii].X();
                                        double yi=ic->GetGeometry()[ii].Y();
                                        double zi=ic->GetGeometry()[ii].Z();
                                        double xj=ic->GetGeometry()[jj].X();
                                        double yj=ic->GetGeometry()[jj].Y();
                                        double zj=ic->GetGeometry()[jj].Z();
                                        double xk=ic->GetGeometry()[kk].X();
                                        double yk=ic->GetGeometry()[kk].Y();
                                        double zk=ic->GetGeometry()[kk].Z();
                                        
                                        
                                        array_1d<double,3> rij = Vector3D(xi,yi,zi,xj,yj,zj);
                                        array_1d<double,3> rik = Vector3D(xi,yi,zi,xk,yk,zk);
                                        array_1d<double,3> rjk = Vector3D(xj,yj,zj,xk,yk,zk);
                                        
                                        double dist_ij = Norm3D(rij);
                                        double dist_ik = Norm3D(rik);
                                        double dist_jk = Norm3D(rjk);
                                        
                                        //if the distance between IS_FREE_SURFACE node and TRIPLE_POINT node is too long, we add a node in the middle of the face
                                        Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();
                                    
                                        double factor = 2.5;	
                                        if ((dist_ik > factor*mesh_size) || (dist_jk > factor*mesh_size) || (dist_ij > factor*mesh_size))
                                        {
                                        id++;
                                        double x = 0.3333333333333*(xi+xj+xk);
                                        double y = 0.3333333333333*(yi+yj+yk);
                                        double z = 0.3333333333333*(zi+zj+zk);
                                
                                        Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

                                        //putting the new node also in an auxiliary list
                                        KRATOS_WATCH("adding nodes to list")										

                                        std::cout << "new node id = " << pnode->Id() << std::endl;
                                        //generating the dofs
                                        for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
                                        {
                                            Node<3>::DofType& rDof = *iii;
                                            Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
                                            (p_new_dof)->FreeDof();
                                        }
                                        Geometry<Node<3> >& geom = ic->GetGeometry();

    // 				      InterpolateOnEdge(geom, ii, kk, step_data_size, pnode);
                                        InterpolateOnFace(geom, ii, jj, kk, step_data_size, pnode);
                                        const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
                                        pnode->X0() = pnode->X() - disp[0];
                                        pnode->Y0() = pnode->Y() - disp[1];
                                        pnode->Z0() = pnode->Z() - disp[2];
                                        KRATOS_WATCH("Added node at the FACE with IS_FREE_SURFACE")
                                        }
                                        
                                        if ((dist_ik < mesh_size * 0.35 ) || (dist_jk < mesh_size * 0.35))
                                            {
                                            ic->GetGeometry()[kk].Set(TO_ERASE,true);
//                                             ic->GetGeometry()[ii].Set(TO_ERASE,true);
//                                             ic->GetGeometry()[jj].Set(TO_ERASE,true);
                                            KRATOS_WATCH("we find some points to be deleted ELAF ELAF ELAF ELAF FACE ELAF kkkkk")
                                            }
                                    

                                            
                                            ///////////////////
                                        
				}
                                
                                ///////
			    }
                            }
                        
  
			//20161129 ajarauta - DELETE IS_STRUCTURE NODE IF IT GETS TOO CLOSE TO AN EDGE FORMED BY TWO TRIPLE_POINT NODES
// 			for(ModelPart::ConditionsContainerType::iterator ic = ThisModelPart.ConditionsBegin() ; ic != ThisModelPart.ConditionsEnd() ; ic++)
// 			{
// 			    if (ic->GetGeometry().size()==3)
// 			    {
// 				unsigned int n_trip = ic->GetGeometry()[0].FastGetSolutionStepValue(TRIPLE_POINT);
// 				n_trip    += ic->GetGeometry()[1].FastGetSolutionStepValue(TRIPLE_POINT);
// 				n_trip    += ic->GetGeometry()[2].FastGetSolutionStepValue(TRIPLE_POINT);
// 				
// 				unsigned int n_struc = ic->GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);			
// 				n_struc    += ic->GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);
// 				n_struc    += ic->GetGeometry()[2].FastGetSolutionStepValue(IS_STRUCTURE);				
// 								
// 				if ((n_struc==ic->GetGeometry().size()) && (n_trip == 2))
// 				{
// 				    int ii = 0;
// 				    int jj = 1;
// 				    int kk = 2;
// 				    if (ic->GetGeometry()[0].FastGetSolutionStepValue(TRIPLE_POINT) == 0.0)
// 				    {
// 				        ii = 1;
// 					jj = 2;
// 					kk = 0;
// 				    }
// 				    if (ic->GetGeometry()[1].FastGetSolutionStepValue(TRIPLE_POINT) == 0.0)
// 				    {
// 				        ii = 0;
// 					jj = 2;
// 					kk = 1;
// 				    }
// 				    // Node with index kk is the one that is not on the contact line
// 				    double xi=ic->GetGeometry()[ii].X();
// 				    double yi=ic->GetGeometry()[ii].Y();
// 				    double zi=ic->GetGeometry()[ii].Z();
// 				    double xj=ic->GetGeometry()[jj].X();
// 				    double yj=ic->GetGeometry()[jj].Y();
// 				    double zj=ic->GetGeometry()[jj].Z();
// 				    double xk=ic->GetGeometry()[kk].X();
// 				    double yk=ic->GetGeometry()[kk].Y();
// 				    double zk=ic->GetGeometry()[kk].Z();
// 				    
// 				    array_1d<double,3> rki = Vector3D(xk,yk,zk,xi,yi,zi);
// 				    array_1d<double,3> rkj = Vector3D(xk,yk,zk,xj,yj,zj);
// 				    
// 				    double dist_ki = Norm3D(rki);
// 				    double dist_kj = Norm3D(rkj);
// 				    
// 				    double nodal_h = ic->GetGeometry()[0].FastGetSolutionStepValue(NODAL_H);
// 				    nodal_h       += ic->GetGeometry()[1].FastGetSolutionStepValue(NODAL_H);
// 				    nodal_h       += ic->GetGeometry()[2].FastGetSolutionStepValue(NODAL_H);
// 				    nodal_h       *= 0.33333333333333333;
// 
// 				    double factor = 0.75;			
// 				    if ( (dist_ki < factor*nodal_h) || (dist_kj < factor*nodal_h) )
// 					ic->GetGeometry()[kk].Set(TO_ERASE,true);
// 				}
// 			    }
// 			}
			
			//20161129 ajarauta -> Find those interior nodes close to the contact line and IS_STRUCTURE nodes to avoid deletion of interface elements
// 			for(ModelPart::NodesContainerType::iterator im = ThisModelPart.NodesBegin() ; im != ThisModelPart.NodesEnd() ; im++)
// 			{
// 			    if ((im->FastGetSolutionStepValue(IS_BOUNDARY) == 0.0) && (im->Z() < 0.00004) && (im->Z() > 0.00002))
// 			        im->GetValue(ERASE_FLAG) = 1.0;
// 			}
// 			
// 			ThisModelPart.Conditions().clear();

                        

			boost::timer auxiliary;
			////////////////////////////////////////////////////////////
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointPointerType;
			//typedef PointerVector<PointType>           PointVector;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef PointVector::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;

// 			int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();

			// bucket types
			//typedef Bucket<3, PointType, ModelPart::NodesContainerType, PointPointerType, PointIterator, DistanceIterator > BucketType;
			//typedef Bins< 3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > StaticBins;
			// bucket types
			typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
		
				
			//*************
			// DynamicBins;	
			typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;
			//typedef Tree< StaticBins > Bin; 			     //Binstree;
			unsigned int bucket_size = 20;
			
			//performing the interpolation - all of the nodes in this list will be preserved
			unsigned int max_results = 50;
			//PointerVector<PointType> res(max_results);
			//NodeIterator res(max_results);
			PointVector res(max_results);
			DistanceVector res_distances(max_results);
			Node<3> work_point(0,0.0,0.0,0.0);
			KRATOS_WATCH(h_factor)
 			//if the remove_node switch is activated, we check if the nodes got too close
			if (rem_nodes==true)
			{
 			
				PointVector list_of_nodes;
				list_of_nodes.reserve(ThisModelPart.Nodes().size());
				for(ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; i_node++)
				{
						(list_of_nodes).push_back(*(i_node.base()));
						//(list_of_nodes).push_back(i_node.base());
				}

				kd_tree  nodes_tree1(list_of_nodes.begin(),list_of_nodes.end(), bucket_size);
				//std::cout<<nodes_tree2<<std::endl;				
			
				unsigned int n_points_in_radius;			
				//radius means the distance, closer than which no node shall be allowd. if closer -> mark for erasing
				double radius;

				for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin();
					in != ThisModelPart.NodesEnd(); in++)
					{
					radius=h_factor*in->FastGetSolutionStepValue(NODAL_H);
					
					work_point[0]=in->X();
					work_point[1]=in->Y();
					work_point[2]=in->Z();
				
					n_points_in_radius = nodes_tree1.SearchInRadius(work_point, radius, res.begin(),res_distances.begin(), max_results);
						if (n_points_in_radius>1)
						{
						if ((in->FastGetSolutionStepValue(IS_BOUNDARY)==0.0 && in->FastGetSolutionStepValue(IS_STRUCTURE)==0.0 && in->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)==0.0) || (in->FastGetSolutionStepValue(IS_BOUNDARY)==0.0 && in->FastGetSolutionStepValue(IS_STRUCTURE)==0.0))
							{
							in->Set(TO_ERASE,true);
							//below is just for the bladder example
							
							}
						}
					
					}				
			
				node_erase.Execute();
			}			
			/////////////////////////////////////////////////////////////////
			/////// 	ALPHA SHAPE		/////////////////////////
			/////////////////////////////////////////////////////////////////
			

			tetgenio in, out, in2, outnew;

			// All indices start from 1.
 			in.firstnumber = 1;

			in.numberofpoints = ThisModelPart.Nodes().size();
			in.pointlist = new REAL[in.numberofpoints * 3];

			//writing the point coordinates in a vector
			ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();

			//reorder node Ids
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
                                (nodes_begin + i)->SetId(i+1);
//				(nodes_begin + i)->Id() = i+1;
			}

			//give the corrdinates to the mesher
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				int base = i*3;
				in.pointlist[base] = (nodes_begin + i)->X();
				in.pointlist[base+1] = (nodes_begin + i)->Y();
				in.pointlist[base+2] = (nodes_begin + i)->Z();
			}
			char tetgen_options[] = "SQJ";


			tetrahedralize(tetgen_options, &in, &out); //with option to remove slivers

			
			double first_part_time = auxiliary.elapsed();
			std::cout << "mesh generation time = " << first_part_time << std::endl;

			//generate Kratos Tetrahedra3D4
			int el_number = out.numberoftetrahedra;

			boost::timer alpha_shape_time;
			//calculate r , center for all of the elements
			std::vector<int> preserved_list(el_number);
			array_1d<double,3> x1,x2,x3,x4, xc;
			int point_base;
			int number_of_preserved_elems = 0;
			for(int el = 0; el< el_number; el++)
			{
				int base = el * 4;

				//coordinates
				point_base = (out.tetrahedronlist[base] - 1)*3;
				x1[0] = out.pointlist[point_base]; 
				x1[1] = out.pointlist[point_base+1]; 
				x1[2] = out.pointlist[point_base+2];

				point_base = (out.tetrahedronlist[base+1] - 1)*3;
				x2[0] = out.pointlist[point_base]; 
				x2[1] = out.pointlist[point_base+1]; 
				x2[2] = out.pointlist[point_base+2];

				point_base = (out.tetrahedronlist[base+2] - 1)*3;
				x3[0] = out.pointlist[point_base]; 
				x3[1] = out.pointlist[point_base+1]; 
				x3[2] = out.pointlist[point_base+2];

				point_base = (out.tetrahedronlist[base+3] - 1)*3;
				x4[0] = out.pointlist[point_base]; 
				x4[1] = out.pointlist[point_base+1]; 
				x4[2] = out.pointlist[point_base+2];

				//calculate the geometrical data 
				//degenerate elements are given a very high radius
				double geometrical_hmin, geometrical_hmax;
				double radius;
				double vol;
				CalculateElementData(x1,x2,x3,x4,vol,xc,radius,geometrical_hmin,geometrical_hmax);

				//calculate the prescribed h
				double prescribed_h = (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(NODAL_H);
				prescribed_h += (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(NODAL_H);
				prescribed_h += (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(NODAL_H);
				prescribed_h += (nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(NODAL_H);
				prescribed_h *= 0.25;

				//check the number of nodes on the wall
				int nb = int( (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(IS_STRUCTURE) );
				nb += int( (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(IS_STRUCTURE) );
				nb += int( (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(IS_STRUCTURE) );
				nb += int((nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(IS_STRUCTURE) );

				//check the number of nodes of bo			node_erase.Execute();undary
				int nfs = int( (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(IS_FREE_SURFACE) );
				nfs += int( (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(IS_FREE_SURFACE) );
				nfs += int( (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(IS_FREE_SURFACE) );
				nfs += int((nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(IS_FREE_SURFACE) );
				
				//check the number of nodes of boundary
				int nfluid = int( (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(IS_FLUID) );
				nfluid += int( (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(IS_FLUID) );
				nfluid += int( (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(IS_FLUID) );
				nfluid += int((nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(IS_FLUID) );
				
				int n_lag = int( (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) );
				n_lag += int( (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) );
				n_lag += int( (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) );
				n_lag += int((nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) );
                                
                                int ntp = int( (nodes_begin + out.tetrahedronlist[base]-1)->FastGetSolutionStepValue(TRIPLE_POINT) );
				ntp += int( (nodes_begin + out.tetrahedronlist[base+1]-1)->FastGetSolutionStepValue(TRIPLE_POINT) );
				ntp += int( (nodes_begin + out.tetrahedronlist[base+2]-1)->FastGetSolutionStepValue(TRIPLE_POINT) );
				ntp += int((nodes_begin + out.tetrahedronlist[base+3]-1)->FastGetSolutionStepValue(TRIPLE_POINT) );
                                
                                
                                if (vol < 1E-40) // slivers
				{
					preserved_list[el] = false;
					KRATOS_WATCH("SLIVER FOUND!!!!")
				}
				else
				{
//                                     if ((ntp >=3 && nb >=3) || (ntp >=2 && nfs >=2)) //at least node node at contact line
// 				    {
// 					preserved_list[el] = false;
//                                     }
				    if(nb == 4) // 4 nodes on the wall
					preserved_list[el] = false;
// 				    else if (ntp >=1) //at least node node at contact line
// 				    {
// 					preserved_list[el] = true;
// 					number_of_preserved_elems += 1;
// 				    }
				    else
				    {
				      if (nfs != 0) //close to the free surface or external
				      {
					if( radius  < prescribed_h * alpha_param ) //the nodes are not all on the boundary
					{
// 					    if(nb == 1 && n_lag == 3)// && nb != 0 && nfv == 3) //20161217 ajarauta: this is for Lap_Young example with "Step" initial shape
// 						preserved_list[el] = false;
// 					    else if(nb == 2 && n_lag == 2)// && nb != 0 && nfv == 3)
// 						preserved_list[el] = false;
// 					    else if(nb == 3 && n_lag == 1)// && nb != 0 && nfv == 3)
// 						preserved_list[el] = false;
// 					    else if(nb == 3 && n_lag == 3)// && (nb == 1 || nb == 2))
// 						preserved_list[el] = false;	
// 					    else if(nint == 2 && n_lag == 3 && (nb == 1 || nb == 2))
// 						preserved_list[el] = false;						    
// 					    else if(nint == 3 && n_lag == 2 && (nb == 1 || nb == 2))
// 						preserved_list[el] = false;					    
// 					    else if(nint == 3 && n_lag == 2)// && nb != 0)
// 						preserved_list[el] = false;
// 					    else if(nint == 2 && n_lag == 3)// && nb != 0)
// 						preserved_list[el] = false;				    
// 					    else
// 					    {
						preserved_list[el] = true; //preserve!!
						number_of_preserved_elems += 1;
// 					    }
					}			
					else
					{
						preserved_list[el] = false;					
					}
				      }
				      else if (nfluid < 3.9)
				      {
					if( radius  < prescribed_h * alpha_param ) //the nodes are not all on the boundary
					{
						preserved_list[el] = true; //preserve!!
						number_of_preserved_elems += 1;
					}			
					else
					{
						preserved_list[el] = false;					
					}				      
				      }
				      else
				      {
					preserved_list[el] = true;
					number_of_preserved_elems += 1;
				      }
				    }
				}

				
// 				//cases:
// 				//4 nodes on the wall - elminate
// 				// at least one node of boundary OR at least one node NOT of fluid --> pass alpha shape
// 				/*
// 				if(nb == 4) // 4 nodes on the wall
// 					preserved_list[el] = false;
// 				else if (nboundary != 0 || nfluid != 4) //close to the free surface or external
// 				{
// 					if( radius  < prescribed_h * alpha_param && //alpha shape says to preserve
// 						nb!=4) //the nodes are not all on the boundary
// 					{
// 						preserved_list[el] = true; //preserve!!
// 						number_of_preserved_elems += 1;
// 					}					
// 				}
// 				else
// 				{
// 					preserved_list[el] = true;
// 					number_of_preserved_elems += 1;
// 				}
// 				*/
// 				if(nb == 4) // 4 nodes on the wall
// 				{
// 					preserved_list[el] = false;
// 					
// 				}
// 				else 
// 				{
// 					//if (nfs != 0 || nfluid != 4) //close to the free surface or external
// 					if (n_lag >= 2 && nb==3 ) //close to the free surface or external
// 					{
// 						if( radius  < prescribed_h * 1.5*alpha_param ) //alpha shape says to preserve
// 							 //the nodes are not all on the boundary
// 						{
// 							preserved_list[el] = true; //preserve!!
// 							number_of_preserved_elems += 1;
// 						}
// 						else
// 						{
// 							preserved_list[el] = false;
// 							
// 						}					
// 					}
// 					else if (nfs != 0 ) //close to the free surface or external
// 					{
// 						if( radius  < prescribed_h * alpha_param ) //alpha shape says to preserve
// 							 //the nodes are not all on the boundary
// 						{
// 							preserved_list[el] = true; //preserve!!
// 							number_of_preserved_elems += 1;
// 						}
// 						else
// 						{
// 							preserved_list[el] = false;
// 							
// 						}
// 					}
// 
// 					
// 					else if (nfluid < 3.9) //close to the free surface or external
// 					{
// 						if( radius  < prescribed_h * alpha_param ) //alpha shape says to preserve
// 							 //the nodes are not all on the boundary
// 						{
// 							preserved_list[el] = true; //preserve!!
// 							number_of_preserved_elems += 1;
// 						}
// 						else
// 						{
// 							preserved_list[el] = false;
// 							
// 						}
// 					}
// 					
// 					else //internal elements should be preserved as much as possible not to create holes
// 					{
// 						if( radius  < prescribed_h * alpha_param * 5.0 ) 
// 						{
// //std::cout << "element not deleted" <<std::endl;
// 							preserved_list[el] = true; //preserve!!
// 							number_of_preserved_elems += 1;
// 						}
// 						else
// 						{
// //std::cout << "sliver removed" << std::endl;
// 							preserved_list[el] = false;
// 							
// 						}
// //						preserved_list[el] = true;
// //						number_of_preserved_elems += 1;
// 					}
// 				}

			}
			std::cout << "time for passing alpha shape" << alpha_shape_time.elapsed() << std::endl;
			
			//freeing unnecessary memory
			in.deinitialize();
			in.initialize();
			
			//setting the new new ids in a aux vector			
//			tetgenio in2;

  			in2.firstnumber = 1;
			in2.numberofpoints = ThisModelPart.Nodes().size();
			in2.pointlist = new REAL[in2.numberofpoints * 3];

			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				int base = i*3;
				in2.pointlist[base] = (nodes_begin + i)->X();
				in2.pointlist[base+1] = (nodes_begin + i)->Y();
				in2.pointlist[base+2] = (nodes_begin + i)->Z();
			}
			std::cout << "qui" << std::endl;
			in2.numberoftetrahedra = number_of_preserved_elems;
			in2.tetrahedronlist = new int[in2.numberoftetrahedra * 4];
			in2.tetrahedronvolumelist = new double[in2.numberoftetrahedra];

			int counter = 0;
			for(int el = 0; el< el_number; el++)
			{
				if( preserved_list[el] ) 
				{
					//saving the compact element list
					int new_base = counter*4;
					int old_base = el*4;
					in2.tetrahedronlist[new_base] = out.tetrahedronlist[old_base];
					in2.tetrahedronlist[new_base+1] = out.tetrahedronlist[old_base+1];
					in2.tetrahedronlist[new_base+2] = out.tetrahedronlist[old_base+2];
					in2.tetrahedronlist[new_base+3] = out.tetrahedronlist[old_base+3];

					//calculate the prescribed h
					double prescribed_h = (nodes_begin + out.tetrahedronlist[old_base]-1)->FastGetSolutionStepValue(NODAL_H);
					prescribed_h += (nodes_begin + out.tetrahedronlist[old_base+1]-1)->FastGetSolutionStepValue(NODAL_H);
					prescribed_h += (nodes_begin + out.tetrahedronlist[old_base+2]-1)->FastGetSolutionStepValue(NODAL_H);
					prescribed_h += (nodes_begin + out.tetrahedronlist[old_base+3]-1)->FastGetSolutionStepValue(NODAL_H);
					prescribed_h *= 0.25; 
					//if h is the height of a perfect tetrahedra, the edge size is edge = sqrt(3/2) h
					//filling in the list of "IDEAL" tetrahedron volumes=1/12 * (edge)^3 * sqrt(2)~0.11785* h^3=
					//0.2165063509*h^3
					

					//chapuza to be checked -  I just didnt want so much of refinement,,,, 
					in2.tetrahedronvolumelist[counter] = 2000.0*0.217*prescribed_h*prescribed_h*prescribed_h;


					//in2.tetrahedronvolumelist[counter] = 0.217*prescribed_h*prescribed_h*prescribed_h;
					//in2.tetrahedronvolumelist[counter] = 0.0004;
					//KRATOS_WATCH(in2.tetrahedronvolumelist[counter])
					counter += 1;
				}
			
			}

			//NOW WE SHALL IDENTIFY FLYING NODES
			/*
			std::vector<double> aux_before_ref(ThisModelPart.Nodes().size());
			
			for (unsigned int i=0; i<aux_before_ref.size();i++)
			{
			aux_before_ref[i]=0.0;
			}
			
			for(int el = 0; el< el_number; el++)
			{
			int old_base=el*4;
			if (preserved_list[el]==true)
				{
			
				//we add a non-zero value for the node of the element that has a non-zero value
				aux_before_ref[out.tetrahedronlist[old_base]-1]+=1;
				aux_before_ref[out.tetrahedronlist[old_base+1]-1]+=1;
				aux_before_ref[out.tetrahedronlist[old_base+2]-1]+=1;
				aux_before_ref[out.tetrahedronlist[old_base+3]-1]+=1;		
				}
			}		
			*/
			//creating a new mesh
			boost::timer mesh_recreation_time;

			//freeing unnecessary memory
			out.deinitialize();
			out.initialize();
			

			//HERE WE ADD THE VOLUME CONSTRAINT  -the desired volume of the equidistant tertrahedra
			//based upon average nodal_h (that is the  "a" switch
			//char regeneration_options[] = "rQJYq1.8anS";
			if (add_nodes==true)
				{
				char mesh_regen_opts[] = "rQJYYqnS";
				tetrahedralize(mesh_regen_opts, &in2, &outnew);
				KRATOS_WATCH("Adaptive remeshing executed")
				}
			else 
				{
				char mesh_regen_opts[] = "rQJYnS";
				tetrahedralize(mesh_regen_opts, &in2, &outnew);
				KRATOS_WATCH("Non-Adaptive remeshing executed")
				}

			//q - creates quality mesh, with the default radius-edge ratio set to 2.0

 			std::cout << "mesh recreation time" << mesh_recreation_time.elapsed() << std::endl;


			//PAVEL

			//putting the new nodes in a spatial container
			//PointerVector< Element >& neighb

			/*
			typedef Node<3> PointType;
			typedef Node<3>::Pointer PointPointerType;
			typedef std::vector<PointType::Pointer>           PointVector;
			typedef PointVector::iterator PointIterator;
			typedef std::vector<double>               DistanceVector;
			typedef std::vector<double>::iterator     DistanceIterator;
			int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();
			// bucket types
			typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
			typedef Bins< 3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > StaticBins;
			// DynamicBins;	
			typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;
			//typedef Tree< StaticBins > Bin; 			     //Binstree;
			*/
 			PointVector list_of_new_nodes;
			
			Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();

			int n_points_before_refinement = in2.numberofpoints;
			//if the refinement was performed, we need to add it to the model part.
			if (outnew.numberofpoints>n_points_before_refinement)
			{
				for(int i = n_points_before_refinement; i<outnew.numberofpoints;i++)
				{
					int id=i+1;
					int base = i*3;
					double& x= outnew.pointlist[base];
					double& y= outnew.pointlist[base+1];
					double& z= outnew.pointlist[base+2];

					Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

					//putting the new node also in an auxiliary list
					//KRATOS_WATCH("adding nodes to list")					
					list_of_new_nodes.push_back( pnode );
					
					//std::cout << "new node id = " << pnode->Id() << std::endl;
					//generating the dofs
					for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
					{
						Node<3>::DofType& rDof = *iii;
						Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
						
						(p_new_dof)->FreeDof();
					}
					
				}
			}	
			
			std::cout << "During refinement we added " << outnew.numberofpoints-n_points_before_refinement<< "nodes " <<std::endl;

						
			bucket_size = 50;
			
			//performing the interpolation - all of the nodes in this list will be preserved
			max_results = 800;
			PointVector results(max_results);
			DistanceVector results_distances(max_results);
			array_1d<double,4> N;
			//int data_size = ThisModelPart.GetNodalSolutionStepTotalDataSize();

			
			//double* work_array;
			//Node<3> work_point(0,0.0,0.0,0.0);

			

			//NOW WE SHALL IDENTIFY LONELY BUT NOT FLYING NODES (compare with  the flying nodes identification above (after first remeshing step))			
			/*
			std::vector<double> aux_after_ref(outnew.numberofpoints);
			
			for (unsigned int i=0; i<aux_after_ref.size();i++)
			{
				aux_after_ref[i]=0.0;
			}
			int el_number_ref= outnew.numberoftetrahedra;
						
			for(int el = 0; el< el_number_ref; el++)
			{
			int base=el*4;
			
				
				//we add a non-zero value for the node of the element that has a non-zero value
				aux_after_ref[outnew.tetrahedronlist[base]-1]+=1;
				aux_after_ref[outnew.tetrahedronlist[base+1]-1]+=1;
				aux_after_ref[outnew.tetrahedronlist[base+2]-1]+=1;
				aux_after_ref[outnew.tetrahedronlist[base+3]-1]+=1;		
				
			}		
			*/
 
			if(outnew.numberofpoints-n_points_before_refinement > 0) //if we added points
			{
				kd_tree  nodes_tree2(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);
				//std::cout<<nodes_tree2<<std::endl;				
				nodes_begin = ThisModelPart.NodesBegin();

				for(int el = 0; el< in2.numberoftetrahedra; el++)
				//for(unsigned int el = 0; el< outnew.numberoftetrahedra; el++)
				{	
					
					int base = el * 4;
					//coordinates
					
					point_base = (in2.tetrahedronlist[base] - 1)*3;
					x1[0] = in2.pointlist[point_base]; 
					x1[1] = in2.pointlist[point_base+1]; 
					x1[2] = in2.pointlist[point_base+2];

					point_base = (in2.tetrahedronlist[base+1] - 1)*3;
					x2[0] = in2.pointlist[point_base]; 
					x2[1] = in2.pointlist[point_base+1]; 
					x2[2] = in2.pointlist[point_base+2];

					point_base = (in2.tetrahedronlist[base+2] - 1)*3;
					x3[0] = in2.pointlist[point_base]; 
					x3[1] = in2.pointlist[point_base+1]; 
					x3[2] = in2.pointlist[point_base+2];

					point_base = (in2.tetrahedronlist[base+3] - 1)*3;
					x4[0] = in2.pointlist[point_base]; 
					x4[1] = in2.pointlist[point_base+1]; 
					x4[2] = in2.pointlist[point_base+2];
					
					//calculate the geometrical data 
					//degenerate elements are given a very high radius
					double geometrical_hmin, geometrical_hmax;
					double radius;
					double vol;

					array_1d<double,3> xc;

					//it calculates the data necessary to perfom the search, like the element center etc., search radius
					//and writes it to radius, xc etc etc
				
					CalculateElementData(x1,x2,x3,x4,vol,xc,radius,geometrical_hmin,geometrical_hmax);
					
					//find all of the nodes in a radius centered in xc
				
					std::size_t number_of_points_in_radius;
					work_point.X() = xc[0]; work_point.Y() = xc[1]; work_point.Z() = xc[2];

					number_of_points_in_radius = nodes_tree2.SearchInRadius(work_point, radius, results.begin(),results_distances.begin(), max_results);
				
							
					//for each of the nodes in radius find if it is inside the element or not
					//if inside interpolate		
				

					Tetrahedra3D4<Node<3> > geom(
						*( (nodes_begin +  in2.tetrahedronlist[base]-1).base() 	), 
						*( (nodes_begin +  in2.tetrahedronlist[base+1]-1).base() 	), 
						*( (nodes_begin +  in2.tetrahedronlist[base+2]-1).base() 	), 
						*( (nodes_begin +  in2.tetrahedronlist[base+3]-1).base() 	) 
						);	
					
									
					//KRATOS_WATCH(results.size())
					for(auto it=results.begin(); it!=results.begin() + number_of_points_in_radius; it++)
	 				{
						bool is_inside=false; 
														
						is_inside = CalculatePosition(x1[0],x1[1],x1[2], 
										x2[0],x2[1],x2[2],
										x3[0],x3[1],x3[2],
										x4[0],x4[1],x4[2],
										(*it)->X(),(*it)->Y(), (*it)->Z(),N);

						if(is_inside == true)
							{	
								
								Interpolate(  geom,  N, step_data_size, *(it) );

							}
											
	 				}		
				}
			}
			
			ThisModelPart.Elements().clear();
			ThisModelPart.Conditions().clear();
			
			//set the coordinates to the original value
			for(auto & list_of_new_node : list_of_new_nodes)
			{				
				const array_1d<double,3>& disp = list_of_new_node->FastGetSolutionStepValue(DISPLACEMENT);
				list_of_new_node->X0() = list_of_new_node->X() - disp[0];
				list_of_new_node->Y0() = list_of_new_node->Y() - disp[1];
				list_of_new_node->Z0() = list_of_new_node->Z() - disp[2];	
			}
			//cleaning unnecessary data
			in2.deinitialize();
			in2.initialize();

			//***********************************************************************************
			//***********************************************************************************
			boost::timer adding_elems;
			//add preserved elements to the kratos
			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
			nodes_begin = ThisModelPart.NodesBegin();
			(ThisModelPart.Elements()).reserve(outnew.numberoftetrahedra);
			
			for(int iii = 0; iii< outnew.numberoftetrahedra; iii++)
			{
				int id = iii + 1;
				int base = iii * 4;
				Tetrahedra3D4<Node<3> > geom(
					*( (nodes_begin +  outnew.tetrahedronlist[base]-1).base() 		), 
					*( (nodes_begin +  outnew.tetrahedronlist[base+1]-1).base() 	), 
					*( (nodes_begin +  outnew.tetrahedronlist[base+2]-1).base() 	), 
					*( (nodes_begin +  outnew.tetrahedronlist[base+3]-1).base() 	) 
					);


#ifdef _DEBUG
ModelPart::NodesContainerType& ModelNodes = ThisModelPart.Nodes();
				if( *(ModelNodes).find( outnew.tetrahedronlist[base]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node","");
				if( *(ModelNodes).find( outnew.tetrahedronlist[base+1]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node","");
				if( *(ModelNodes).find( outnew.tetrahedronlist[base+2]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node","");
				if( *(ModelNodes).find( outnew.tetrahedronlist[base+3]).base() == *(ThisModelPart.Nodes().end()).base() ) 
					KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node","");
#endif

				Element::Pointer p_element = rReferenceElement.Create(id, geom, properties);
				(ThisModelPart.Elements()).push_back(p_element);

			}
			std::cout << "time for adding elems" << adding_elems.elapsed() << std::endl;;
			ThisModelPart.Elements().Sort();
	
			boost::timer adding_neighb;
//			//filling the neighbour list 
			ModelPart::ElementsContainerType::const_iterator el_begin = ThisModelPart.ElementsBegin();
			for(ModelPart::ElementsContainerType::const_iterator iii = ThisModelPart.ElementsBegin();
				iii != ThisModelPart.ElementsEnd(); iii++)
			{
				//Geometry< Node<3> >& geom = iii->GetGeometry();
				int base = ( iii->Id() - 1 )*4;

				(iii->GetValue(NEIGHBOUR_ELEMENTS)).resize(4);
				WeakPointerVector< Element >& neighb = iii->GetValue(NEIGHBOUR_ELEMENTS);

				for(int i = 0; i<4; i++)
				{
					int index = outnew.neighborlist[base+i];
					if(index > 0)
						neighb(i) = *((el_begin + index-1).base());
					else
						neighb(i) = Element::WeakPointer();
				}
			}
			std::cout << "time for adding neigbours" << adding_neighb.elapsed() << std::endl;;

						
			
		

			//***********************************************************************************
			//***********************************************************************************
			//mark boundary nodes 
			//reset the boundary flag
			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
			{
				in->FastGetSolutionStepValue(IS_BOUNDARY) = 0;
			}


			//***********************************************************************************
			//***********************************************************************************
			boost::timer adding_faces;

			(ThisModelPart.Conditions()).reserve(outnew.numberoftrifaces   );

			//creating the faces
			for(ModelPart::ElementsContainerType::const_iterator iii = ThisModelPart.ElementsBegin();
				iii != ThisModelPart.ElementsEnd(); iii++)
			{
				
				int base = ( iii->Id() - 1 )*4;
				
				//create boundary faces and mark the boundary nodes 
				//each face is opposite to the corresponding node number so
				// 0 ----- 1 2 3
				// 1 ----- 0 3 2
				// 2 ----- 0 1 3
				// 3 ----- 0 2 1

				//node 1

				//if( neighb(0).expired()  );
				if( outnew.neighborlist[base] == -1)
				{
					CreateBoundaryFace(1, 2, 3, ThisModelPart,   0, *(iii.base()), rReferenceBoundaryCondition, properties );
				}
				//if(neighb(1).expired() );
				if( outnew.neighborlist[base+1] == -1)
				{
					CreateBoundaryFace(0,3,2, ThisModelPart,   1, *(iii.base()), rReferenceBoundaryCondition, properties );
				}
				if( outnew.neighborlist[base+2] == -1)
				//if(neighb(2).expired() );
				{
					CreateBoundaryFace(0,1,3, ThisModelPart,   2, *(iii.base()), rReferenceBoundaryCondition, properties );
				}
				if( outnew.neighborlist[base+3] == -1)
				//if(neighb(3).expired() );
				{
					CreateBoundaryFace(0,2,1, ThisModelPart,   3, *(iii.base()), rReferenceBoundaryCondition, properties );
				}

			}
			outnew.deinitialize();
			outnew.initialize();
			std::cout << "time for adding faces" << adding_faces.elapsed() << std::endl;;


			
			//here we remove lonely nodes that are not teh flying nodes, but are the lonely nodes inside water vol
			/*
			for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
			{
			//if this is a node added during refinement, but isnt contained in any element
			if (aux_after_ref[in->GetId()]==0 && in->GetId()>aux_before_ref.size() && in->GetId()<aux_after_ref.size())
				{
				in->Set(TO_ERASE,true);
				KRATOS_WATCH("This is that ugly ugly lonely interior node a666a. IT SHALL BE TERRRRMINATED!!!!!!")
				KRATOS_WATCH(in->GetId())
				}
			//if this was an interior node after first step and became single due to derefinement - erase it			
			if (aux_after_ref[in->GetId()]==0 && in->GetId()<aux_before_ref.size() && aux_before_ref[in->GetId()]!=0)
				{
				in->Set(TO_ERASE,true);
				KRATOS_WATCH("This is that ugly ugly lonely interior node. IT SHALL BE TERRRRMINATED!!!!!!")
				}
			}
			*/


			double second_part_time = auxiliary.elapsed();
			std::cout << "second part time = " << second_part_time - first_part_time << std::endl;


			KRATOS_CATCH("")
		}

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
		virtual std::string Info() const{return "";}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const{}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const{}


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
		boost::numeric::ublas::bounded_matrix<double,3,3> mJ; //local jacobian
		boost::numeric::ublas::bounded_matrix<double,3,3> mJinv; //inverse jacobian
		array_1d<double,3> mc; //center pos
		array_1d<double,3> mRhs; //center pos
		
		array_1d<double,3> Vector3D(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1)
		{
		    array_1d<double,3> r01;
		    r01[0] = x1 - x0;
		    r01[1] = y1 - y0;
		    r01[2] = z1 - z0;
		    return r01;
		}
		
		double Norm3D(const array_1d<double,3>& a)
		{
		    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
		}
		
		double DotProduct3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
		{
		    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
		}
		
		double Angle2vecs3D(const array_1d<double,3>& a, const array_1d<double,3>& b)
		{
		    double norm_a = Norm3D(a);
		    double norm_b = Norm3D(b);
		    double temp = 0.0;
		    if (norm_a*norm_b == 0.0)
		        temp = 0.0;
		    else
		        temp = DotProduct3D(a,b)/(norm_a*norm_b);
		    return acos(temp);
		}
		
		array_1d<double,3> NormalizeVec3D(array_1d<double,3>& input)
		{
		    double norm = Norm3D(input);
		    if (norm != 0.0)
		    {
		        input[0] /= norm;
		        input[1] /= norm;
		        input[2] /= norm;
		    }
		    return input;
		}		
		
		double DistPointLine3D(const double x0, const double y0, const double z0, const double x1, const double y1, const double z1, const double x2, const double y2, const double z2)
		{
		    // This function returns the distance of point 2 from the line formed by points 0 and 1
		    array_1d<double,3> r01 = Vector3D(x0,y0,z0,x1,y1,z1); //vector director (d) of the line L
		    array_1d<double,3> r02 = Vector3D(x0,y0,z0,x2,y2,z2);
		    double proj = (DotProduct3D(r02,r01))/Norm3D(r01);
		    
		    array_1d<double,3> d_vec = NormalizeVec3D(r01);
		    
		    array_1d<double,3> x = ZeroVector(3);
		    x[0] = r02[0] - proj*d_vec[0];
		    x[1] = r02[1] - proj*d_vec[1];
		    x[2] = r02[2] - proj*d_vec[2];
		    
		    return Norm3D(x);
		}




		void CreateBoundaryFace(const int& i1, const int& i2, const int& i3, ModelPart& ThisModelPart, const int& outer_node_id, Element::Pointer origin_element, Condition const& rReferenceBoundaryCondition, Properties::Pointer properties)
		{
			KRATOS_TRY

			Geometry<Node<3> >& geom = origin_element->GetGeometry();
			//mark the nodes as free surface
			geom[i1].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
			geom[i2].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
			geom[i3].FastGetSolutionStepValue(IS_BOUNDARY) = 1;

			//generate a face condition
			Condition::NodesArrayType temp;
			temp.reserve(3);
			temp.push_back(geom(i1)); 
			temp.push_back(geom(i2));
			temp.push_back(geom(i3));
			Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D3< Node<3> >(temp) );
			//Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Triangle3D< Node<3> >(temp) );
			int id = (origin_element->Id()-1)*4;
			Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, temp, properties);	

			//assigning the neighbour node
			(p_cond->GetValue(NEIGHBOUR_NODES)).clear();
			(p_cond->GetValue(NEIGHBOUR_NODES)).push_back( Node<3>::WeakPointer( geom(outer_node_id) ) );
			(p_cond->GetValue(NEIGHBOUR_ELEMENTS)).clear();
			(p_cond->GetValue(NEIGHBOUR_ELEMENTS)).push_back( Element::WeakPointer( origin_element ) );
			ThisModelPart.Conditions().push_back(p_cond);
			KRATOS_CATCH("")

		}
		
		//////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////////////////////
		void Interpolate( Tetrahedra3D4<Node<3> >& geom, const array_1d<double,4>& N, 
				  unsigned int step_data_size,
      				Node<3>::Pointer pnode)
		{
			unsigned int buffer_size = pnode->GetBufferSize();


			for(unsigned int step = 0; step<buffer_size; step++)
			{	
				
						//getting the data of the solution step
				double* step_data = (pnode)->SolutionStepData().Data(step);
				
											
				double* node0_data = geom[0].SolutionStepData().Data(step);
				double* node1_data = geom[1].SolutionStepData().Data(step);
				double* node2_data = geom[2].SolutionStepData().Data(step);
				double* node3_data = geom[3].SolutionStepData().Data(step);
				
				//copying this data in the position of the vector we are interested in
				for(unsigned int j= 0; j<step_data_size; j++)
				{ 

					step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j] + N[3]*node3_data[j];
								

				}						
			}
			//now we assure that the flag variables are set coorect!! since we add nodes inside of the fluid volume only
			//we manually reset the IS_BOUNDARY, IS_FLUID, IS_STRUCTURE, IS_FREE_SURFACE values in a right way
			//not to have values, like 0.33 0.66 resulting if we would have been interpolating them in the same way 		
			//as the normal variables, like Velocity etc		
			
                        pnode->FastGetSolutionStepValue(IS_BOUNDARY)=0.0;
			pnode->FastGetSolutionStepValue(IS_STRUCTURE)=0.0;
			pnode->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)=0.0;
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=0.0;
			pnode->FastGetSolutionStepValue(IS_FLUID)=1.0;
			pnode->FastGetSolutionStepValue(MEAN_CURVATURE_3D)=0.0;
			pnode->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
			pnode->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
		}
		
		
		void InterpolateOnFace( Geometry<Node<3> >& geom, int point1, int point2, int point3, unsigned int step_data_size, Node<3>::Pointer pnode)
		{
		    unsigned int buffer_size = pnode->GetBufferSize();
		    KRATOS_WATCH(buffer_size)

		    for(unsigned int step = 0; step<buffer_size; step++)
		    {	
			//getting the data of the solution step
			double* step_data = (pnode)->SolutionStepData().Data(step);
					    
			if (point1>2 or point1<0 or point2>2 or point2<0 or point3>2 or point3<0 or (point1==point2) or (point1==point3) or (point3==point2))
			{
			    KRATOS_WATCH(point1)
			    KRATOS_WATCH(point2)
			    KRATOS_WATCH(point3)
// 			    KRATOS_ERROR(std::logic_error,"THE EDGE POINTS ARE INVALID ", "");
			}
												    
			double* node0_data = geom[point1].SolutionStepData().Data(step);
			double* node1_data = geom[point2].SolutionStepData().Data(step);
			double* node2_data = geom[point3].SolutionStepData().Data(step);						
					    
			//copying this data in the position of the vector we are interested in
			for(unsigned int j= 0; j<step_data_size; j++)
			{ 
			    step_data[j] = 0.3333333333333*(node0_data[j] + node1_data[j] + node2_data[j]);
			}						
		    }
		    //now we assure that the flag variables are set coorect!! since we add nodes inside of the fluid volume only
		    //we manually reset the IS_BOUNDARY, IS_FLUID, IS_STRUCTURE, IS_FREE_SURFACE values in a right way
		    //not to have values, like 0.33 0.66 resulting if we would have been interpolating them in the same way 		
		    //as the normal variables, like Velocity etc			

		    pnode->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)=0.0;
		    pnode->FastGetSolutionStepValue(IS_BOUNDARY)=1.0;
		    pnode->FastGetSolutionStepValue(IS_FLUID)=1.0;
// 		    pnode->FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT_X)=0.0;
		    pnode->FastGetSolutionStepValue(MEAN_CURVATURE_3D) =0.0;
		    pnode->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
		    pnode->FastGetSolutionStepValue(NODAL_LENGTH) = 0.0;
		    pnode->FastGetSolutionStepValue(NORMAL_TRIPLE_POINT) = ZeroVector(3);
		    pnode->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
				    
		    pnode->Set(TO_ERASE,false);;

		    if (pnode->FastGetSolutionStepValue(IS_INTERFACE)>0.2)
			pnode->FastGetSolutionStepValue(IS_INTERFACE)=1.0;
		    else 
			pnode->FastGetSolutionStepValue(IS_INTERFACE)=0.0;
		    if (pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)>0.2)
		    {
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=1.0;
			pnode->FastGetSolutionStepValue(NORMAL_EQUILIBRIUM_X) = 0.0;
			pnode->FastGetSolutionStepValue(NORMAL_EQUILIBRIUM_Y) = 0.0;
			pnode->FastGetSolutionStepValue(NORMAL_EQUILIBRIUM_Z) = 0.0;
		    }
		    else 
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=0.0;
		    if (pnode->FastGetSolutionStepValue(FLAG_VARIABLE)>0.2)
			pnode->FastGetSolutionStepValue(FLAG_VARIABLE)=1.0;
		    else 
			pnode->FastGetSolutionStepValue(FLAG_VARIABLE)=0.0;	
// 		    if (pnode->FastGetSolutionStepValue(TRIPLE_POINT)>0.7)
// 			pnode->FastGetSolutionStepValue(TRIPLE_POINT)=1.0;
// 		    else 
// 		    {
// 			pnode->FastGetSolutionStepValue(TRIPLE_POINT)	= 0.0;
// 			pnode->FastGetSolutionStepValue(NODAL_LENGTH)	= 0.0;
// 			pnode->FastGetSolutionStepValue(NORMAL_TP_X)	= 0.0;
// 			pnode->FastGetSolutionStepValue(NORMAL_TP_Y)	= 0.0;
// 			pnode->FastGetSolutionStepValue(NORMAL_TP_Z)	= 0.0;
// 		    }
		    if (pnode->FastGetSolutionStepValue(IS_STRUCTURE)>0.7)
		    {
			pnode->FastGetSolutionStepValue(IS_STRUCTURE)=1.0;
			pnode->FastGetSolutionStepValue(IS_INTERFACE)=0.0;
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=0.0;
			pnode->FastGetSolutionStepValue(MEAN_CURVATURE_3D) =0.0;
			pnode->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
			pnode->FastGetSolutionStepValue(NODAL_LENGTH) = 0.0;
			pnode->FastGetSolutionStepValue(NORMAL_TRIPLE_POINT) = ZeroVector(3);
			pnode->FastGetSolutionStepValue(CONTACT_ANGLE) = 0.0;
		    }
		    else 
			pnode->FastGetSolutionStepValue(IS_STRUCTURE)=0.0;	
		    
		}
		  
		void InterpolateOnEdge( Geometry<Node<3> >& geom, int point1, int point2, unsigned int step_data_size, Node<3>::Pointer pnode)
		{
		    unsigned int buffer_size = pnode->GetBufferSize();

		    for(unsigned int step = 0; step<buffer_size; step++)
		    {	
			//getting the data of the solution step
			double* step_data = (pnode)->SolutionStepData().Data(step);
					    
			if (point1>2 or point1<0 or point2>2 or point2<0 or (point1==point2))
			{
			    KRATOS_WATCH(point1)
			    KRATOS_WATCH(point2)
// 			    KRATOS_ERROR(std::logic_error,"THE EDGE POINTS ARE INVALID ", "");
			}
												    
			double* node0_data = geom[point1].SolutionStepData().Data(step);
			double* node1_data = geom[point2].SolutionStepData().Data(step);						
					    
			//copying this data in the position of the vector we are interested in
			for(unsigned int j= 0; j<step_data_size; j++)
			{ 
			    step_data[j] = 0.5*(node0_data[j] + node1_data[j]);
			}						
		    }
		    //now we assure that the flag variables are set coorect!! since we add nodes inside of the fluid volume only
		    //we manually reset the IS_BOUNDARY, IS_FLUID, IS_STRUCTURE, IS_FREE_SURFACE values in a right way
		    //not to have values, like 0.33 0.66 resulting if we would have been interpolating them in the same way 		
		    //as the normal variables, like Velocity etc		
				  
		    //pnode->FastGetSolutionStepValue(IS_BOUNDARY)=1.0;
		    //pnode->FastGetSolutionStepValue(FLAG_VARIABLE)=1.0;			

		    pnode->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)=0.0;
		    pnode->FastGetSolutionStepValue(IS_BOUNDARY)=1.0;
		    pnode->FastGetSolutionStepValue(FLAG_VARIABLE) = 1.0;
		    pnode->FastGetSolutionStepValue(IS_FLUID)=1.0;
// 		    pnode->FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT_X)=0.0;
		    pnode->FastGetSolutionStepValue(MEAN_CURVATURE_3D)=0.0;
		    pnode->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
		    pnode->FastGetSolutionStepValue(DISPLACEMENT_Z) = 0.0;
				    
		    pnode->Set(TO_ERASE,false);;

		    if (pnode->FastGetSolutionStepValue(IS_INTERFACE)>0.4)
		    {
			pnode->FastGetSolutionStepValue(IS_INTERFACE)  = 1.0;
// 			pnode->FastGetSolutionStepValue(FLAG_VARIABLE) = 1.0;
		    }
		    else
		    {
			pnode->FastGetSolutionStepValue(IS_INTERFACE)  = 0.0;
// 			pnode->FastGetSolutionStepValue(FLAG_VARIABLE) = 0.0;
		    }
		    if (pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)>0.4)
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=1.0;
		    else 
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=0.0;
		    if (pnode->FastGetSolutionStepValue(IS_STRUCTURE) > 0.9)
			pnode->FastGetSolutionStepValue(IS_STRUCTURE) = 1.0;
		    else
			pnode->FastGetSolutionStepValue(IS_STRUCTURE) = 0.0;
		    if (pnode->FastGetSolutionStepValue(TRIPLE_POINT) > 0.9)
			pnode->FastGetSolutionStepValue(TRIPLE_POINT) = 1.0;
		    else
		    {
			pnode->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
		    }
		}		  
		
		
		
		
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
		
		inline bool CalculatePosition(	const double x0, const double y0, const double z0,
						const double x1, const double y1, const double z1,
   						const double x2, const double y2, const double z2,
						const double x3, const double y3, const double z3,
						const double xc, const double yc, const double zc,
						array_1d<double,4>& N		
					  )
		{ 
			
					
			double vol = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

			double inv_vol = 0.0;
			if(vol < 1e-26)
			  {
				//KRATOS_THROW_ERROR(std::logic_error,"element with zero vol found","");
				KRATOS_WATCH("WARNING!!!!!!!!!! YOU HAVE A SLIVER HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
				return false;
			  }
			else
			  {
				inv_vol = 1.0 / vol;
				N[0] = CalculateVol(x1,y1,z1,x3,y3,z3,x2,y2,z2,xc,yc,zc) * inv_vol;
				N[1] = CalculateVol(x3,y3,z3,x0,y0,z0,x2,y2,z2,xc,yc,zc) * inv_vol;
				N[2] = CalculateVol(x3,y3,z3,x1,y1,z1,x0,y0,z0,xc,yc,zc) * inv_vol;
				N[3] = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) * inv_vol;

						
				if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >=0.0 &&   N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <=1.0)
				{
				  //if the xc yc zc is inside the tetrahedron return true
				return true;
				}
				else 
				  return false;
			  }
			
			  
		}	
		//***************************************
		//***************************************
		inline void CalculateCenterAndSearchRadius(const double x0, const double y0, const double z0,
						const double x1, const double y1, const double z1,
    						const double x2, const double y2, const double z2,
    						const double x3, const double y3, const double z3,
	  					double& xc, double& yc, double& zc, double& R
					       )
		{
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
		//*****************************************************************************************
		//*****************************************************************************************
		void CalculateElementData(const array_1d<double,3>& c1,
			const array_1d<double,3>& c2,
			const array_1d<double,3>& c3,
			const array_1d<double,3>& c4,
			double& volume,
			array_1d<double,3>& xc,
			double& radius,
			double& hmin,
			double& hmax
			)
		{
			KRATOS_TRY
			const double x0 = c1[0]; const double y0 = c1[1]; const double z0 = c1[2];
			const double x1 = c2[0]; const double y1 = c2[1]; const double z1 = c2[2];
			const double x2 = c3[0]; const double y2 = c3[1]; const double z2 = c3[2];
			const double x3 = c4[0]; const double y3 = c4[1]; const double z3 = c4[2];

// 			KRATOS_WATCH("111111111111111111111");
			//calculate min side lenght 
			//(use xc as a auxiliary vector) !!!!!!!!!!!!
			double aux;
			noalias(xc) = c4; noalias(xc)-=c1; aux = inner_prod(xc,xc); hmin = aux; hmax = aux;
			noalias(xc) = c3; noalias(xc)-=c1; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c2; noalias(xc)-=c1; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c4; noalias(xc)-=c2; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c3; noalias(xc)-=c2; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			noalias(xc) = c4; noalias(xc)-=c3; aux = inner_prod(xc,xc); if(aux<hmin) hmin = aux; else if(aux>hmax) hmax = aux;
			hmin = sqrt(hmin);
			hmax = sqrt(hmax);

			//calculation of the jacobian
			mJ(0,0) = x1-x0; mJ(0,1) = y1-y0; mJ(0,2) = z1-z0;
			mJ(1,0) = x2-x0; mJ(1,1) = y2-y0; mJ(1,2) = z2-z0;
			mJ(2,0) = x3-x0; mJ(2,1) = y3-y0; mJ(2,2) = z3-z0;
// 			KRATOS_WATCH("33333333333333333333333");

			//inverse of the jacobian
			//first column
			mJinv(0,0) = mJ(1,1)*mJ(2,2) - mJ(1,2)*mJ(2,1);
			mJinv(1,0) = -mJ(1,0)*mJ(2,2) + mJ(1,2)*mJ(2,0);
			mJinv(2,0) = mJ(1,0)*mJ(2,1) - mJ(1,1)*mJ(2,0);		
			//second column
			mJinv(0,1) = -mJ(0,1)*mJ(2,2) + mJ(0,2)*mJ(2,1);
			mJinv(1,1) = mJ(0,0)*mJ(2,2) - mJ(0,2)*mJ(2,0);
			mJinv(2,1) = -mJ(0,0)*mJ(2,1) + mJ(0,1)*mJ(2,0);
			//third column
			mJinv(0,2) = mJ(0,1)*mJ(1,2) - mJ(0,2)*mJ(1,1);
			mJinv(1,2) = -mJ(0,0)*mJ(1,2) + mJ(0,2)*mJ(1,0);
			mJinv(2,2) = mJ(0,0)*mJ(1,1) - mJ(0,1)*mJ(1,0);
			//calculation of determinant (of the input matrix)

// 			KRATOS_WATCH("44444444444444444444444444");
			double detJ = mJ(0,0)*mJinv(0,0) 
				+ mJ(0,1)*mJinv(1,0)
				+ mJ(0,2)*mJinv(2,0);

			volume = detJ * 0.16666666667;
// 			KRATOS_WATCH("55555555555555555555555");

			if(volume < 1e-3 * hmax*hmax*hmax)  //this is a sliver and we will remove it
			{
 			//KRATOS_WATCH("666666666666666666666666666");
				//very bad element // ser a very high center for it to be eliminated
				xc[0]=0.0; xc[1]=0.0; xc[2]=0.0;
				radius = 10000000;
// 			KRATOS_WATCH("777777777777777777777777");
			}
			else
			{

// 			KRATOS_WATCH("888888888888888888888888");
				double x0_2 = x0*x0 + y0*y0 + z0*z0;
				double x1_2 = x1*x1 + y1*y1 + z1*z1; 
				double x2_2 = x2*x2 + y2*y2 + z2*z2; 
				double x3_2 = x3*x3 + y3*y3 + z3*z3; 

				//finalizing the calculation of the inverted matrix
				mJinv /= detJ;

				//calculating the RHS
				//calculating the RHS
				mRhs[0] = 0.5*(x1_2 - x0_2);
				mRhs[1] = 0.5*(x2_2 - x0_2);
				mRhs[2] = 0.5*(x3_2 - x0_2);

				//calculate position of the center
				noalias(xc) = prod(mJinv,mRhs);

				//calculate radius
				radius = pow(xc[0] - x0,2);
				radius		  += pow(xc[1] - y0,2);
				radius		  += pow(xc[2] - z0,2);
				radius = sqrt(radius);
// 			KRATOS_WATCH("999999999999999999999999999");
			}
			KRATOS_CATCH("")
		}

		template< class T, std::size_t dim >
				class PointDistance{
					public:
						double operator()( T const& p1, T const& p2 ){
							double dist = 0.0;
							for( std::size_t i = 0 ; i < dim ; i++){
								double tmp = p1[i] - p2[i];
								dist += tmp*tmp;
							}
							return sqrt(dist);
						}
				};		
		
		template< class T, std::size_t dim >
			class DistanceCalculator{
				public:
					double operator()( T const& p1, T const& p2 ){
						double dist = 0.0;
						for( std::size_t i = 0 ; i < dim ; i++){
							double tmp = p1[i] - p2[i];
							dist += tmp*tmp;
						}
						return dist;
					}
			};

		///@} 
		///@name Private Operators
		///@{ 

		///@} 
		///@name Private Operations
		///@{ 


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
		TetGenPfemModeler& operator=(TetGenPfemModeler const& rOther);


		///@}    

	}; // Class TetGenPfemModeler 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		TetGenPfemModeler& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const TetGenPfemModeler& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_TETGEN_PFEM_MODELER_H_INCLUDED  defined 
