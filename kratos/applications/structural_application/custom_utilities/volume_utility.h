/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
*   Last Modified by:    $Author: virginia $
*   Date:                $Date: 2011-02-24 12:14:58 $
*   Revision:            $Revision: 1.12 $
*
* ***********************************************************/

#if !defined(KRATOS_VOLUME_UTILITY_H_INCLUDED )
#define  KRATOS_VOLUME_UTILITY_H_INCLUDED

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"
// #include "solving_strategies/strategies/solving_strategy.h"
// #include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
// #include "solving_strategies/convergencecriterias/convergence_criteria.h"


#include "structural_application.h"

//default builder and solver
// #include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

namespace Kratos
{
    /**
     * auxiliary functions for calculating volume of a 3D structure
     */
    class VolumeUtility
    {
        public:
            
            /**
             * Type Definitions 
             */
//             typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;
            
            /** 
             * Counted pointer of ContactUtility 
             */
            KRATOS_CLASS_POINTER_DEFINITION( VolumeUtility );
            
            typedef std::size_t IndexType;
            
            typedef ModelPart::ConditionsContainerType ConditionsArrayType;
            
            typedef Geometry<Node<3> > GeometryType;
            
            typedef Properties PropertiesType;

            /**
             * Life Cycle 
             */
            
            /**
             * Constructor.
             */
            VolumeUtility( int echo_level )
            {
                mEchoLevel = echo_level;
            }
            
            /**
             * Destructor.
             */
            virtual ~VolumeUtility() {}
            
            /**
             * Operations
             */
            
          
          /**
             * This function calculates the volume of a closed surface 
             * @param mr_model_part the model part
             */
      // VM !!!!!!!!!  
 //double CalculateVolume(ModelPart& model_part, Vector& center,double area,double volume)  // input: elements, output: center, area, volume
	    double CalculateVolume(ModelPart& model_part)  // input: elements, output: center, area, volume
	    // using & above to replace the initial value for center, area and volume
	    {
		Vector center;
		noalias(center) = ZeroVector(3);
                double area=0.0;
		double volume=0.0;
		int fail=0.0;
		//int idx;
		//ModelPart::ElementsContainerType& rElements = model_part.Elements();  // to introduce elements container
		//ModelPart::NodesContainerType& rNodes = model_part.Nodes(); // to introduce Nodes container
		//Elem* elm;  // classe de elementos
		//Element::GeometryType& geom = iel->GetGeometry();  // takes the element geometry
		//unsigned int number_of_nodes = geom.PointsNumber(); // gets the total number of nodes of a element
		//KRATOS_WATCH (number_of_nodes);
		//int inode; // "i" counter
		vector<double> ce(3),ne(3); //ce = center, ne= normal
		noalias(ce) = ZeroVector(3);
		noalias(ne) = ZeroVector(3);
		//double ce[3],ne[3]; // vectors
		//double v1[3],v2[3]; // vectors
		vector<double> v1 (3), v2(3) ; //auxiliary vectors
		noalias(v1) = ZeroVector(3);
		noalias(v2) = ZeroVector(3);
		//const vector<double> coords(3); // coordenadas dos nós
		double athird=1.0/3.0;
		//double area=0.0;
		//double volume=0.0;
		//vector<double> center (3);
		//if(!ElemsPlex.Count()) return fail; // 
		unsigned int number_of_total_element_in_the_model_part = model_part.ElementsEnd() - model_part.ElementsBegin();
		//KRATOS_WATCH (number_of_total_element_in_the_model_part);

// 		if(number_of_total_element_in_the_model_part==0) 
// 		  {
// 		    return fail;
// 		  } // if # elem = 0, fail!
		
		//center=ZeroVector(3);
		//vzero(center);
		//idx=0;
		
		//while ( ( elm = (Elem*)ElemsPlex.next(idx)) != NULL)
		// CENTER
// 		  Element::GeometryType& geom = iel->GetGeometry();  // takes the element geometry
// 		  unsigned int number_of_nodes = geom.PointsNumber(); // gets the total number of nodes of a element
// 		  std::vector< Point<3, double > > Puntos; /// vector of points for the element
// 		  Puntos.resize( number_of_nodes);

		  for( ModelPart::ElementsContainerType::iterator iel = model_part.ElementsBegin();  // loop over elements
		      iel != model_part.ElementsEnd();
		      iel++)	

	    	  {
		    //unsigned int number_of_total_element_in_the_model_part = model_part.ElementsEnd() - model_part.ElementsBegin(); // number of elements
		  Element::GeometryType& geom = iel->GetGeometry();  // takes the element geometry
		  unsigned int number_of_nodes = geom.PointsNumber(); // gets the total number of nodes of a element
		  std::vector< Point<3, double > > Puntos; /// vector of points for the element
		  Puntos.resize( number_of_nodes);

		    for (unsigned int point = 0; point < number_of_nodes; point ++ )
		      {
			Puntos[point][0] = geom[point].X();
			Puntos[point][1] = geom[point].Y();
			Puntos[point][2] = geom[point].Z();
		    
		      }
		    //noalias(ce)=Puntos[0];
		    //vcopy(coords[0],ce);
		    //noalias(ce)+=Puntos[1];   
		    //vadd(ce,coords[1],ce);
		    //noalias(ce)+=Puntos[2];
		    //vadd(ce,coords[2],ce);
		    noalias(ce)=(Puntos[0]+Puntos[1]+Puntos[2]);
		    ce*=athird;
		    //vscale(ce,athird);
		    noalias(center)+=ce;
		    //vadd(center,ce,center);
		    //noalias(x) = Puntos[0]+Puntos[1]; // sum of coordinates of 2 poimts
		  }
		//noalias(center)*=(1.0/number_of_total_element_in_the_model_part);
		//KRATOS_WATCH (number_of_total_element_in_the_model_part);
		//KRATOS_WATCH(center); 
		center *= (1.0/number_of_total_element_in_the_model_part);
		//KRATOS_WATCH("center after division"); 
		//KRATOS_WATCH(center); 
		//vscale(center,1.0/((double)ElemsPlex.Count()));
		//idx=0;
		
		//  NORM
		for( ModelPart::ElementsContainerType::iterator iel = model_part.ElementsBegin();
		      iel != model_part.ElementsEnd();
		      iel++)
		//while ( ( elm = (Elem*)ElemsPlex.next(idx)) != NULL) 
		  {
		    //unsigned int number_of_total_element_in_the_model_part = model_part.ElementEnd() - model_part.ElementsBegin(); // number of elements
		    Element::GeometryType& geom = iel->GetGeometry();  // takes the element geometry
		    unsigned int number_of_nodes = geom.PointsNumber(); // gets the total number of nodes of a element
		    std::vector< Point<3, double > > Puntos; /// vector of points for the element
		    Puntos.resize( number_of_nodes);

		    for (unsigned int point = 0; point < number_of_nodes; point ++ )
		      {
			Puntos[point][0] = geom[point].X();
			Puntos[point][1] = geom[point].Y();
			Puntos[point][2] = geom[point].Z();
		    
		      }
		      noalias(ce)=Puntos[0];
		      //vcopy(coords[0],ce);
		      noalias(ce)+=Puntos[1];   
		      //vadd(ce,coords[1],ce);
		      noalias(ce)+=Puntos[2];
		      //vadd(ce,coords[2],ce);
		      ce*=athird;
		      //vscale(ce,athird);
		      noalias(center)-=ce;
		      //noalias(center)-=ce;
		      //vsub(ce,center,ce);//numerically better translating to center
		      //noalias(v1)=Puntos[1];
		      //noalias(v1)-=Puntos[0];
                      noalias(v1)=Puntos[1]-Puntos[0];
		      //vsub(coords[1],coords[0],v1);
		      noalias(v2)=Puntos[2]-Puntos[1];
		      //noalias(v2)=Puntos[2];
		      //noalias(v2)-=Puntos[1];
		      //vsub(coords[2],coords[1],v2);
		      //CalcVol:CrossProduct(ne,v1,v2); // Uncomment when using the function Cross Product
		      // TRIAL CROSS PRODUCT
		      ne[0] = v1[1]*v2[2] - v1[2]*v2[1];
		      ne[1] = v1[2]*v2[0] - v1[0]*v2[2];
		      ne[2] = v1[0]*v2[1] - v1[1]*v2[0];
		      // END TRIAL
		      //vcross(v1,v2,ne); KK
		      double Norm=MathUtils<double>::Norm3(ne); // norm of ne is computed
		      //MathUtils<double>::Norm3(ne); // Norm of vector size 3 is computed
		      area+=Norm;  // TO CHECK!!!!!!!!!!
		      //area+=norm_1(ne);  // TO CHECK!!!!!!!!!!
		      //area+=vlength(ne);
		      //KRATOS_WATCH(volume); 
		      double Dot=MathUtils<double>::Dot3(ce, ne);  // inner product (dot/scalar product)
		      volume+=Dot;
		      //KRATOS_WATCH("after dot product"); 
		      //KRATOS_WATCH(volume); 
		      //volume+=inner_prod (ce, ne); // inner product (dot/scalar product)
		      //volume+=vdot(ce,ne);
		      //area*=0.5;  // used for triangular elements
		      //KRATOS_WATCH(area);
		      //volume/=-6.0; // used for triangular elements
		      //KRATOS_WATCH(volume); 
		  }
		      area*=0.5;  // used for triangular elements
		      //KRATOS_WATCH(area);
		      //volume/=-6.0;
		      volume/=6.0; // used for triangular elements  // before with a "-" above
		//KRATOS_WATCH(volume); 

		//Element::GeometryType& geom = iel->GetGeometry();  // takes the element geometry
		//elm=(Elem*)(ElemsPlex)(ElemsPlex.first());
		for( ModelPart::ElementsContainerType::iterator iel = model_part.ElementsBegin();  // loop over elements
		      iel != model_part.ElementsEnd();
		      iel++)	

		    {
		    //unsigned int number_of_total_element_in_the_model_part = model_part.ElementsEnd() - model_part.ElementsBegin(); // number of elements
		    Element::GeometryType& geom = iel->GetGeometry();  // takes the element geometry
		    unsigned int number_of_nodes = geom.PointsNumber(); // gets the total number of nodes of a element
		    std::vector< Point<3, double > > Puntos; /// vector of points for the element
		    Puntos.resize( number_of_nodes);
		    //KRATOS_WATCH(number_of_nodes);

		    if (number_of_nodes==3)  // TRIANGLE
		    //if (Element::GeometryType& geom = iel->GetGeometry()==Triangle)
		    //if(elm->GiveElemType()==Triangle)
		      {
		      //KRATOS_WATCH("TRIANGLE");
		      //area*=0.5;
		      //KRATOS_WATCH(area);
		      //volume/=-6.0;
		      //KRATOS_WATCH(volume);
		      }
 		     if (number_of_nodes==4)  // QUADRILATERAL
//  				//else if(Element::GeometryType& geom = iel->GetGeometry()==Quadrilateral)
// 				//else if(elm->GiveElemType()==Quadrilateral)
 			{
 		  	KRATOS_WATCH("QUADRILATERAL, please change the formula for volume");
			//volume/=-3.0;
 			} 
		    else 
			{
			fail=1.0;
			}
		      }
		   //return fail;
		   return volume;
		   //KRATOS_WATCH (volume);
	        
     //int mEchoLevel;
     //ModelPart& model_part;    
      }


// VM !!!!!!!!!!
      
            
        private:
       

    


            int mEchoLevel;
            //double k_contact; //VM
            //double k_contact_t;  //VM
            
    };//class VolumeUtility
}  /* namespace Kratos.*/

#endif /* KRATOS_VOLUME_UTILITY_H_INCLUDED  defined */
