/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Last Modified by:    $Author: pablo $
//   Date:                $Date: 2011-09-20 12:51:34 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_POINT_LOCATION_INCLUDED )
#define  KRATOS_POINT_LOCATION_INCLUDED



// System includes
#include <string>
#include <iostream> 
using namespace std;


// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h" 


namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
	typedef  ModelPart::NodesContainerType NodesContainerType;
	typedef  ModelPart::ElementsContainerType ElementsContainerType;
	
  
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
  class PointLocation
	//: public
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of CalculateNodalAreaProcess
     // KRATOS_CLASS_POINTER_DEFINITION(PointLocation);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
		  
		  
		///// Default constructor.
		PointLocation(ModelPart& model_part)  : m_model_part(model_part) {}      

		/// Destructor.
//		virtual ~UbicacionDepoint() 
	//	{}
	 
				
    /// this function only returns whether an element has been found or not. to be called AFTER execute!!!!
    bool found()  
	{
				KRATOS_TRY 
   return misinside; 
   		KRATOS_CATCH ("");
	} 
 

     /// MAIN Function to be called from python, the others below are subroutines. 2D version
      int Find2D( double mxpoint, double mypoint, Vector& Avector)   //note that except mxpoint and mypoint all variables are passed by reference  ( & ). this way we'll be able to modify their values and return data
	{     

			KRATOS_TRY
			
			
		//ModelPart::ElementsContainerType& rElements = r_model_part.Elements();
		melem = 0;		
		array_1d<double,3> N;      // nodal coordinates
		misinside = false ;     //we start with 0, if this value doesnt change then no element has been found ---> we're out of the boundaries of the model
		N[0] = N[1] = N[2] = 0.0 ;
				
		
		for(ModelPart::ElementsContainerType::iterator i = m_model_part.ElementsBegin(); 
						 i!=m_model_part.ElementsEnd(); i++)
			{	
				melem = melem + 1;

				Geometry< Node<3> >& geom = i->GetGeometry();  // current element's coordinates
					
				misinside = IsInside ( geom, mxpoint, mypoint, N);  //we check if it's inside
 
                //cout << endl << "quepaso?" << N[0] << endl << endl  << mA2 << endl << mA3 << endl << endl ;
 
                if (misinside==true) break ;       //if it's inside this element then there's no need to keep searching
			}		
			
			
		if (misinside == false) melem = 0 ; //if it has found no element after the loop we set melem to zero, this will tell the .py file that there's no elem

        mA1 = N[0]; mA2 = N[1] ; mA3 = N[2];          //we save in two different formats, might be useful depending on the case
        mA[0] = N[0]; mA[1] = N[1] ; mA[2] = N[2];	
         		
		// at the end we have the element and the area coordinates of our point
		
		Avector[0] = N[0];
		Avector[1] = N[1];
		Avector[2] = N[2];
		
		return melem; //note: It's NOT THE ELEMENT ID!! it's the position (current-begin)
		
		KRATOS_CATCH ("");	
		
	}
 
      
     /// MAIN Function to be called from python, the others below are subroutines. 3D
      int Find3D( double mxpoint, double mypoint, double mzpoint, Vector& Avector)   //note that except mxpoint and mypoint all variables are passed by reference  ( & ). this way we'll be able to modify their values and return data
	{     

			KRATOS_TRY
			
			
		//ModelPart::ElementsContainerType& rElements = r_model_part.Elements();
		melem = 0;		
		array_1d<double,4> N;      // nodal coordinates
		misinside = false ;     //we start with 0, if this value doesnt change then no element has been found ---> we're out of the boundaries of the model
		N[0] = N[1] = N[2] = N[3] = 0.0 ;
				
		
		for(ModelPart::ElementsContainerType::iterator i = m_model_part.ElementsBegin(); 
						 i!=m_model_part.ElementsEnd(); i++)
			{	
				melem = i->Id();

				Geometry< Node<3> >& geom = i->GetGeometry();  // current element's coordinates
					
				misinside = IsInside ( geom, mxpoint, mypoint, mzpoint, N);  //we check if it's inside
 
                if (misinside==true) break ;       //if it's inside this element then there's no need to keep searching
			}		
			
			
		if (misinside == false) melem = 0 ; //if it has found no element after the loop we set melem to zero, this will tell the .py file that there's no elem

        mA1 = N[0]; mA2 = N[1] ; mA3 = N[2] ; mA4 = N[3] ;	
        mA[0] = N[0]; mA[1] = N[1] ; mA[2] = N[2] ; mA[3] = N[3] ;	
         		
		// at the end we have the element and the area coordinates of our point
		
		Avector[0] = N[0];
		Avector[1] = N[1];
		Avector[2] = N[2];
		Avector[3] = N[3];
		
		return melem; //note: It's THE ELEMENT ID!! (it's NOT the position (current-begin))
		
		KRATOS_CATCH ("");	
		
	}

		//***************************************   
		//***************************************
    bool IsInside (  Geometry<Node<3> >&geom, double xpoint, double ypoint, array_1d<double,3>& N)     //this function "draws" circles around the nodes (using r= max distance between nodes) to see if our point is inside the element
    {
		double x1=geom[0].X(); double x2=geom[1].X(); double x3=geom[2].X();
		double y1=geom[0].Y(); double y2=geom[1].Y(); double y3=geom[2].Y();
		
		double edge1= (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
		double edge2= (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3);      //creating the 3 sides of the triangle
		double edge3= (x1-x3)*(x1-x3)+(y1-y3)*(y1-y3);
		
		//cout << endl << "los lados" << " x1 "<< x1 << " x2 "<< x2 << " x3 " << x3 << " y1 " << y1 << " y2 " << y2 << " y3 " << y3 << "lado 1 "  << lado1 << "   lado 2 " << lado2  << "    lado 3" << lado3 << endl << endl ;
		
		double radius;
		radius = edge1;
		if (edge2>radius) radius=edge2    ;                       //i pick the longest
	    if (edge3>radius) radius=edge3    ;
		
		bool isinside=false;   	
		bool isincircles;
		isincircles=true;     
		
		double delta;
		delta= ( (xpoint-x1)*(xpoint-x1)  ) + ( (ypoint-y1)*(ypoint-y1) ) ;  // checking if the point is inside the circle surronding the first node of the triangle
		delta *= 0.9 ;
		if (delta>radius)  isincircles=false;      // if it is outside, then this is not the element we're looking for
		
		delta=0.0;
		delta= ( (xpoint-x2)*(xpoint-x2) ) + ( (ypoint-y2)*(ypoint-y2) ) ;
		delta *= 0.9 ;
		if (delta>radius)  isincircles=false;
		 
		delta=0.0; 
		delta= ( (xpoint-x3)*(xpoint-x3) ) + ( (ypoint-y3)*(ypoint-y3) ) ;        //checked the third node
		delta *= 0.9 ;
		if (delta>radius)   isincircles=false;
		
 		
        if (isincircles)    // this means there's a possibily that the point is inside the triangle boundaries, to be sure we have to calculate the area coordinates.
        {
			isinside= CalculatePosition(geom, xpoint, ypoint, N);      //it's the same as the original CalculatePosition from projection.h, but the zc coord has been removed
		}
            //now we know for sure. we return our answer. the area coordinates are returned in N     
        return isinside;	
	} //end IsInside


		//***************************************   
		//***************************************

    
    bool IsInside (  Geometry<Node<3> >&geom, double xpoint, double ypoint, double zpoint, array_1d<double,4>& N)      //  3D   : this function "draws" spheres around the nodes (using r= max distance between nodes) to see if our point is inside the element. this is the first check. if passed then calculates the weight nodal funcions of the points to make sure it's inside the element.
    {

		array_1d<double, 3 > node_coord;             // 
        array_1d<double, 3 > neigh_coord;  
		double radius = 0.0; //actually its ^2. but no need to find the real one since we're comparing squares
		double dist_node_neigh = 0.0; 
		for(unsigned int i = 0; i < geom.size() ; i++) {//edge i
			node_coord[0] = geom[i].X();
			node_coord[1] = geom[i].Y();
			node_coord[2] = geom[i].Z();
			for(unsigned int j = 0; j < geom.size() ; j++) {//edge i
				neigh_coord[0] = geom[j].X();
				neigh_coord[1] = geom[j].Y();
				neigh_coord[2] = geom[j].Z();
				dist_node_neigh = ( pow((node_coord[0]- neigh_coord[0]),2) + pow((node_coord[1]- neigh_coord[1]),2) + pow((node_coord[2]- neigh_coord[2]),2) ) ; // square distance between node and neighbour
				if ((dist_node_neigh>radius) && (i!=j)) radius=dist_node_neigh; //saving the longest found up to now
			}//closing j loop
		} //closing i loop
		
		bool isinside=false;   	
		bool isincircles;
		isincircles=true;     
		double delta = 0.0;
		
		for (unsigned int i=0; i!=4; ++i) { //cheching the 4 nodes of the element
		delta= ( pow((geom[i].X()-xpoint ),2) + pow((geom[i].Y()-ypoint ),2) + pow((geom[i].Z()-zpoint ),2) ) ;  // checking if the point is inside the circle surronding the first node of the triangle
		delta *= 0.95 ;
		if (delta>radius)  isincircles=false;      // if it is outside, then this is not the element we're looking for
		}


        if (isincircles)    // this means there's a possibily that the point is inside the triangle boundaries, to be sure we have to calculate the area coordinates.
        {
			isinside= CalculatePosition(geom, xpoint, ypoint, zpoint, N);      //it's the same as the original CalculatePosition from projection.h
		}
            //now we know for sure. we return our answer. the area coordinates are returned in N
         
            
        return isinside;

		
		
	} //end IsInside
    
    
		//***************************************   for 3D
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
        
    	//*************************************** //FOR 2D
		//***************************************
		inline double CalculateVol(	const double x0, const double y0,
						const double x1, const double y1,
    						const double x2, const double y2 )
		{
			return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
		}
		
		//***************************************
		//***************************************
		
		
		
    	inline bool CalculatePosition(	Geometry<Node<3> >&geom, //for triangles (2D) (2 double  + N(3) )
						const double xc, const double yc,      
						array_1d<double,3>& N	)
		{
			 double x0 = geom[0].X();double  y0 = geom[0].Y(); 
			 double x1 = geom[1].X();double  y1 = geom[1].Y(); 
			 double x2 = geom[2].X();double  y2 = geom[2].Y(); 
			 
			double area = CalculateVol(x0,y0,x1,y1,x2,y2);
			double inv_area = 0.0;
			if(area == 0.0)
			  {

// 				KRATOS_ERROR(std::logic_error,"element with zero area found","");
				//The interpolated node will not be inside an elemente with zero area
				return false;
				
			  }
			else
			  {
				inv_area = 1.0 / area;
			  }
			
			  
			  N[0] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
			  N[1] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;
			  N[2] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;
		
			
			if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <=1.0 && N[1]<= 1.0 && N[2] <= 1.0) //if (xc,yc) is inside the triangle return true
				return true;
		 	
			return false;
		} //end CalculatePosition	
		
		
		//***************************************
		//***************************************		
		
		
				inline bool CalculatePosition(	Geometry<Node<3> >&geom,      //for tetraedra (3D) ( 3 doubles + N(4) )
						const double xc, const double yc, const double zc,
						array_1d<double,4>& N		
					  )
		{ 
			
			 double x0 = geom[0].X();double  y0 = geom[0].Y();double  z0 = geom[0].Z();
			 double x1 = geom[1].X();double  y1 = geom[1].Y();double  z1 = geom[1].Z();
			 double x2 = geom[2].X();double  y2 = geom[2].Y();double  z2 = geom[2].Z();
			 double x3 = geom[3].X();double  y3 = geom[3].Y();double  z3 = geom[3].Z();	
	
			double vol = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

			double inv_vol = 0.0;
			if(vol < 0.0000000000001)
			  {

// 				KRATOS_ERROR(std::logic_error,"element with zero vol found","");
				//The interpolated node will not be inside an elemente with zero volume
				return false;
				KRATOS_WATCH("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
			  }
			else
			  {
				inv_vol = 1.0 / vol;
			  }

			  N[0] = CalculateVol(x1,y1,z1,x3,y3,z3,x2,y2,z2,xc,yc,zc) * inv_vol;
			  N[1] = CalculateVol(x3,y3,z3,x0,y0,z0,x2,y2,z2,xc,yc,zc) * inv_vol;
			  N[2] = CalculateVol(x3,y3,z3,x1,y1,z1,x0,y0,z0,xc,yc,zc) * inv_vol;
	   		  N[3] = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) * inv_vol;

						
			if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >=0.0 &&
			   N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <=1.0)
			//if the xc yc zc is inside the tetrahedron return true
				return true;
			
			return false;
		}
		
		
		//***************************************
		//***************************************		
		
		///FUNCTIONS TO BE CALLED TO INTERPOLATE DATA FROM THE MODELPART
		
		//default point, one scalar variable
		double ReturnDefaultPointData_scalar(Variable<double>& ThisVariable)
		{
		ModelPart::ElementsContainerType::iterator it_elem = m_model_part.Elements().find(melem);				
		Geometry< Node<3> >& geom = it_elem->GetGeometry();
		unsigned int number_of_nodes = it_elem->GetGeometry().size(); //3 for 2d case, 4 for 3d case (tetraedra)
		double out_scalar_temp	= 0.0; //temp scalar to add te contribution of the nodes
		for (unsigned int i=0 ; i != number_of_nodes; ++i) //looping the 3/4 nodes of the element
			out_scalar_temp += (geom[i].FastGetSolutionStepValue(ThisVariable , 0)) * mA[i]; //adding the contribution of the i node (interpolating)
		return out_scalar_temp;	
		}

		//default point, vector(3) variable		
		array_1d<double,3> ReturnDefaultPointData_vector(Variable<array_1d<double,3> >& ThisVariable) //Vector& out_vector,
		{
		ModelPart::ElementsContainerType::iterator it_elem = m_model_part.Elements().find(melem);	
		Geometry< Node<3> >& geom = it_elem->GetGeometry();				
		unsigned int number_of_nodes = it_elem->GetGeometry().size(); //3 for 2d case, 4 for 3d case (tetraedra)
		array_1d<double,3> out_vector_temp;  
		for (unsigned int index = 0 ; index !=3 ; ++index) out_vector_temp[index] = 0.0; //initializing the temp vector
		for (unsigned int i=0 ; i != number_of_nodes; ++i) 
			out_vector_temp += (geom[i].FastGetSolutionStepValue(ThisVariable , 0)) * mA[i]; //adding the contribution of the i node (interpolating
		//out_vector=out_vector_temp;			
		return out_vector_temp;	
		}			
		
		//custom point, one scalar variable		
		double ReturnCustomPointData_scalar(int element, Vector& Avector , Variable<double>& ThisVariable)
		{
		ModelPart::ElementsContainerType::iterator it_elem = m_model_part.Elements().find(element);	
		Geometry< Node<3> >& geom = it_elem->GetGeometry();
		unsigned int number_of_nodes = it_elem->GetGeometry().size(); //3 for 2d case, 4 for 3d case (tetraedra)				
		double out_scalar_temp	= 0.0;
		for (unsigned int i=0 ; i != number_of_nodes; ++i) 
			out_scalar_temp += (geom[i].FastGetSolutionStepValue(ThisVariable , 0)) * Avector[i]; //adding the contribution of the i node (interpolating)
		return out_scalar_temp;					
		}
		
		//custom point, vector (3) variable		
		array_1d<double,3> ReturnCustomPointData_vector(int element, Vector& Avector , Variable<array_1d<double,3> >& ThisVariable)
		{
		ModelPart::ElementsContainerType::iterator it_elem = m_model_part.Elements().find(element);	
		Geometry< Node<3> >& geom = it_elem->GetGeometry();
		unsigned int number_of_nodes = it_elem->GetGeometry().size(); //3 for 2d case, 4 for 3d case (tetraedra)			
		array_1d<double,3> out_vector_temp;  
		for (unsigned int index = 0 ; index !=3 ; ++index) out_vector_temp[index] = 0.0; //initializing the temp vector
		for (unsigned int i=0 ; i != number_of_nodes; ++i) 
			out_vector_temp += (geom[i].FastGetSolutionStepValue(ThisVariable , 0)) * Avector[i]; //adding the contribution of the i node (interpolating
		return out_vector_temp;					
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


            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
    //protected:
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
	bool misinside;
	int melem;
	array_1d<double,4> mA;
	double mA1;
	double mA2;
	double mA3;
	double mA4;
	ModelPart& m_model_part;
	//unsigned int mdomain_size;
        
        
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
    //  UbicacionDepoint& operator=(UbicacionDepoint const& rOther);

      /// Copy constructor.

        
      ///@}    
        
    }; // Class UbicacionDepoint

  
}  // namespace Kratos.

#endif // KRATOS_POINT_LOCATION_INCLUDED  defined 


