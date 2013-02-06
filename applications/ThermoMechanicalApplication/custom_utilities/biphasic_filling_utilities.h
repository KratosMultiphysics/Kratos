/*
==============================================================================
KratosIncompressibleFluidApplication 
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
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: kkamran $
//   Date:                $Date: 2008-10-13 06:58:23 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_BIPHASIC_FILLING_UTILITIES_INCLUDED )
#define  KRATOS_BIPHASIC_FILLING_UTILITIES_INCLUDED



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
//#include "geometries/tetrahedra_3d_4.h"
#include "geometries/point.h"
#include "thermo_mechanical_application.h"
// #include "custom_conditions/environment_contact.h"
//#include "includes/variables.h"



namespace Kratos
{
 	
class BiphasicFillingUtilities
 {
  public:

		//**********************************************************************************************
		//**********************************************************************************************
		double CreateAutoExitAssignAirSmagorinsky(ModelPart& ThisModelPart, double y_wall, double C_Smagorinsky)
		{			
			KRATOS_TRY;
			int node_size = ThisModelPart.Nodes().size();
			for (int ii = 0; ii < node_size; ii++)
			 {
                 ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
				 double str_flag = it->GetValue(IS_STRUCTURE);
				 double slip_flag = it->GetSolutionStepValue(IS_SLIP);

				 if (str_flag == 0.0 && slip_flag>=10.0)
					 return 1.0;
			 }
			// if there is no dry node
		    double is_dry_node = 0.0;
#pragma omp parallel for firstprivate(node_size)
			for (int ii = 0; ii < node_size; ii++)
			 {
                ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
				double dist = it->FastGetSolutionStepValue(DISTANCE);
				double slip_flag = it->GetSolutionStepValue(IS_SLIP);
				double str_flag = it->GetValue(IS_STRUCTURE);

				if(dist > 0.0){
					is_dry_node = 1.0;
					//slip_flag=10.0 refers to the boundary nodes with well-defined normal
					if(slip_flag==10.0 && str_flag!=0.0){
						it->SetValue(IS_STRUCTURE,0.0); 
						it->SetValue(Y_WALL,y_wall); 
					}
				}
			 }

			//assign smagorinsky at air element
			AirSmagorinskey(ThisModelPart, C_Smagorinsky);

		    return is_dry_node;
		
			KRATOS_CATCH("")
		}
	//**********************************************************************************************
/*for node in fluid_model_part.Nodes:
	slip_flag = node.GetSolutionStepValue(IS_SLIP)
	nd_dist = node.GetSolutionStepValue(DISTANCE)
        if((slip_flag == 20.0 or slip_flag == 30.0 )):# 
	  if(nd_dist< 0.0):
	    node.SetValue(IS_STRUCTURE,1.0)	  
	    node.SetValue(Y_WALL,y_wall_val)
	  else:
	    node.SetValue(IS_STRUCTURE,0.0)	  
	    node.SetValue(Y_WALL,y_wall_val*y_wall_fac)	    

	if(slip_flag == 10.0):
	  if(nd_dist< 0.0):
	    node.SetValue(Y_WALL,y_wall_val)	      
	    node.SetValue(IS_STRUCTURE,1.0)	  
	  else:
	    node.SetValue(Y_WALL,y_wall_val*y_wall_fac)
	    node.SetValue(IS_STRUCTURE,0.0)*/
	//**********************************************************************************************
    void AssignSmoothBoundaryAirExit(ModelPart& ThisModelPart, bool air_exit_flag, const double y_wall_val, const double  y_wall_fac)
		{			
		  KRATOS_TRY;
			int node_size = ThisModelPart.Nodes().size();
			double is_str = 1.0;
			if(air_exit_flag) is_str = 0.0;

#pragma omp parallel for firstprivate(node_size)
			for (int ii = 0; ii < node_size; ii++)
			 {
                 ModelPart::NodesContainerType::iterator it = ThisModelPart.NodesBegin() + ii;
				double dist = it->FastGetSolutionStepValue(DISTANCE);
				double slip_flag = it->GetSolutionStepValue(IS_SLIP);

				if(slip_flag == 20.0 || slip_flag == 30.0 )//edges(20) and corners(30) are automatic air exits till they are wetten
				  if(dist<0.0){
					    it->SetValue(IS_STRUCTURE,1.0);	  
					    it->SetValue(Y_WALL,y_wall_val);}
				  else{
						it->SetValue(IS_STRUCTURE,0.0);	  
						it->SetValue(Y_WALL,y_wall_val*y_wall_fac);}	
				else if(slip_flag == 10.0)//smooth boundaries(10), if dry, can be air exit or not
				  {
				  if(dist<0.0){
					    it->SetValue(IS_STRUCTURE,1.0);	  
					    it->SetValue(Y_WALL,y_wall_val);}
				  else{
						it->SetValue(IS_STRUCTURE,is_str);	  
						it->SetValue(Y_WALL,y_wall_val*y_wall_fac);}	
				  }
			}
		  KRATOS_CATCH("")
		}

  private:
 	
	void AirSmagorinskey(ModelPart& ThisModelPart, double C_Smagorinsky)
	{
        int elem_size = ThisModelPart.Elements().size();

        #pragma omp parallel for firstprivate(elem_size)
        for(int ii = 0; ii<elem_size; ii++)
         {
            PointerVector< Element>::iterator iel=ThisModelPart.ElementsBegin()+ii;
			double dist_sign = 1.0;
            Geometry< Node<3> >& geom = iel->GetGeometry();
            for(unsigned int i =0; i<geom.size(); i++)
            {
				double dist = geom[i].FastGetSolutionStepValue(DISTANCE);
				if(dist_sign*dist < 0.0){
					dist_sign = -1.0;
					break;}
			}
			// to be sure to not apply to the cutted elements and just to all air elements
			if(dist_sign == 1.0)
				iel->SetValue(C_SMAGORINSKY, C_Smagorinsky );

		}
	}

 };

}  // namespace Kratos.

#endif // KRATOS_BIPHASIC_FILLING_UTILITIES_INCLUDED  defined 


