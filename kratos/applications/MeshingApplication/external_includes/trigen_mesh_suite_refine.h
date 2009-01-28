//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-03-11 10:29:26 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_TRIGEN_MODELER_H_INCLUDED )
#define  KRATOS_TRIGEN_MODELER_H_INCLUDED

 

// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>


#include "triangle.h" 

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "meshing_application.h"




namespace Kratos
{
	extern "C" {
		void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
		//void trifree();
	}


			
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
	class TriGenModeler  
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TriGenModeler
		KRATOS_CLASS_POINTER_DEFINITION(TriGenModeler);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		TriGenModeler() :
		    mJ(ZeroMatrix(2,2)), //local jacobian
		    mJinv(ZeroMatrix(2,2)), //inverse jacobian
		    mC(ZeroVector(2)), //dimension = number of nodes
		    mRhs(ZeroVector(2)){} //dimension = number of nodes

		/// Destructor.
		virtual ~TriGenModeler(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{


		//*******************************************************************************************
		//*******************************************************************************************
		void RefineMesh(
			ModelPart& ThisModelPart , 
			Element const& rReferenceElement, 
			Condition const& rReferenceBoundaryCondition,
			double my_alpha = 1.4)
		{

			KRATOS_TRY

			struct triangulateio in;
			struct triangulateio out;
			struct triangulateio vorout;
			
			//sizing array of coordinates
			in.numberofpoints = ThisModelPart.Nodes().size();
			in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

			//set the array of attributes 
			in.numberofpointattributes = 0; 
			in.pointattributelist = (REAL *) NULL; //in my opinion this is not needed

			//reserving memory for the element list
			in.numberoftriangles =  ThisModelPart.Elements().size();
			in.trianglelist = (int *) malloc(in.numberoftriangles * 3 * sizeof(int) ];

			//reserving area for volume constraints
			in.trianglearealist = new REAL[in.numberoftetrahedra];			

			in.pointmarkerlist = (int *) NULL;
			in.regionlist = (REAL *) NULL; //in my opinion not needed

			out.pointattributelist = (REAL *) NULL; //in my opinion not needed
			out.trianglelist = (int *) NULL;          
			out.triangleattributelist = (REAL *) NULL; //in my opinion not needed
			out.neighborlist = (int *) NULL; //in my opinion not needed
			//*******************************************************************

			//writing the points coordinates in a vector
			ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				int base = i*2;
				//int base = ((nodes_begin + i)->Id()   -  1 ) * 2;

				//from now on it is consecutive
				(nodes_begin + i)->Id() = i+1;

				in.pointlist[base] = (nodes_begin + i)->X();
				in.pointlist[base+1] = (nodes_begin + i)->Y();

				//providing to the mesher the nodal attributea
				...
			}

			//provide the list of triangles and the area desired to the mesher
			unsigned int i = 0;
			for(ModelPart::ElementsContainerType iel =  ThisModelPart.ElementsBegin(); iel!=  ThisModelPart.ElementsEnd(); iel++)
			{
				Geometry< Node<3> >& geom = iel->GetGeometry(); 
				int base = i*3;
				in.trianglelist[base]   = geom[0]->Id();
				in.trianglelist[base+1] = geom[1]->Id();
				in.trianglelist[base+2] = geom[2]->Id();

				double h =  geom[0].FastGetSolutionStepValue(NODAL_H);
				h +=  geom[1].FastGetSolutionStepValue(NODAL_H);
				h +=  geom[2].FastGetSolutionStepValue(NODAL_H);
				h *= 0.33333333333333;

				//a given triangle (equilatera) of edge h has area A = sqrt(3)/4* h*h = 0.43301*h*h
				double desired_area = 0.43301 * h*h;

				//set the desired area in the list
				in.trianglearealist[i] = desired_area;
				
				//update counter
				i = i + 1;
			}	

			//erase elements and conditions
			ThisModelPart.Elements().clear();
			ThisModelPart.Conditions().clear();

			//perform the meshing - it also finds the neighbours
			triangulate("rNPBn", &in, &out, &vorout); //version with no elemental neighbours

			//generate the new node list
			unsigned int original_size = ThisModelPart.Nodes().size();
			if( out.numberofpoints > original_size;  )
			{
				for(unsigned int i = original_size; i < out.numberofpoints; i++)
				{
					int id = original_size+1
					int base = i*3;
					double& x = out.pointlist[base];		
					double& y = out.pointlist[base+1];
					double& z = out.pointlist[base+2];
					
					//create new node with the x,y,z corresponding and the attribute list given
					ThisModelPart.AddNode(id,x,y,z, attributelist);
				}
			}

			//properties to be used in the generation of the elements
			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);

			//generate Kratos triangle
			int el_number = out.numberoftriangles;

				
			//note that the node list can not be changed starting from here
			ModelPart::NodesContainerType& ModelNodes = ThisModelPart.Nodes();

			for(int el = 0; el< el_number; el++)
			{
				int id = el + 1;
				int base = el * 3;

				//generating the geometry
				Geometry<Node<3> > temp;
				temp.push_back( *((ModelNodes).find( out.trianglelist[base]).base()	) );
				temp.push_back( *((ModelNodes).find( out.trianglelist[base+1]).base()	) );
				temp.push_back( *((ModelNodes).find( out.trianglelist[base+2]).base()	) );

				//generating a new element
				Element::Pointer p_element = rReferenceElement.Create(id, temp, properties);
				ThisModelPart.Elements().push_back(p_element);
			}
KRATOS_WATCH("trigenmodeler qua")
	
			ThisModelPart.Elements().Unique();
			
			for(ModelPart::ElementsContainerType::iterator iii = ThisModelPart.ElementsBegin();	iii != ThisModelPart.ElementsEnd(); iii++)
			{
				int base = ( iii->Id() - 1 )*3;
				
				ModelPart::ElementsContainerType::iterator el_neighb;
				/*each face is opposite to the corresponding node number so
				 0 ----- 1 2 
				 1 ----- 2 0
				 2 ----- 0 1 
				*/

				////finding boundaries and creating the "skin"
				//
				//********************************************************************
				//first face
				el_neighb = (ThisModelPart.Elements()).find( out.neighborlist[base] );
				if( el_neighb == elements_end )
				{
					//Generate condition
					Condition::NodesArrayType temp;
					temp.reserve(2);
					temp.push_back(iii->GetGeometry()(1)); 
					temp.push_back(iii->GetGeometry()(2));
				
					Geometry< Node<3> >::Pointer cond = 
						Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp) );
					int id = (iii->Id()-1)*3;
					
					Condition::Pointer p_cond = 
						Condition::Pointer(new Condition(id, cond, properties) );	

					ThisModelPart.Conditions().push_back(p_cond);

				}
				
				//********************************************************************
				//second face
				el_neighb = (ThisModelPart.Elements()).find( out.neighborlist[base+1] );
				//if( el != ThisModelPart.Elements().end() )
				if( el_neighb == elements_end )
				{
					//Generate condition
					Condition::NodesArrayType temp;
					temp.reserve(2);
					temp.push_back(iii->GetGeometry()(2)); 
					temp.push_back(iii->GetGeometry()(0));
					
					Geometry< Node<3> >::Pointer cond = 
						Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp) );
					int id = (iii->Id()-1)*3+1;
					//
					Condition::Pointer p_cond = 
						Condition::Pointer(new Condition(id, cond, properties) );
					
					ThisModelPart.Conditions().push_back(p_cond);
									

				}

				//********************************************************************
				//third face
				el_neighb = (ThisModelPart.Elements()).find( out.neighborlist[base+2] );
				if( el_neighb == elements_end )
				{
//					Generate condition
					Condition::NodesArrayType temp;
					temp.reserve(2);
					temp.push_back(iii->GetGeometry()(0)); 
					temp.push_back(iii->GetGeometry()(1));
					Geometry< Node<3> >::Pointer cond = 
						Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp) );
					int id = (iii->Id()-1)*3+2;
					
					Condition::Pointer p_cond = 
						Condition::Pointer(new Condition(id, cond, properties) );
					
					ThisModelPart.Conditions().push_back(p_cond);
					

				}

			
			}

KRATOS_WATCH("free surface identified")

			free( in.pointlist );
			free( in.numberoftriangles );
			free( in.trianglelist );
			free( in.trianglearealist );

/*			free( in.pointmarkerlist);
			free( in.pointmarkerlist );
			free( in.regionlist );*/
	
			free( out.pointattributelist );
			free( out.trianglelist );
// 			free( out.triangleattributelist );
			free( out.neighborlist );
KRATOS_WATCH("everything freed")

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
		 boost::numeric::ublas::bounded_matrix<double,2,2> mJ; //local jacobian
		 boost::numeric::ublas::bounded_matrix<double,2,2> mJinv; //inverse jacobian
		 array_1d<double,2> mC; //center pos
		 array_1d<double,2> mRhs; //center pos



		///@} 
		///@name Private Operators
		///@{ 
		//returns false if it should be removed
		bool AlphaShape(double alpha_param, Geometry<Node<3> >& pgeom)
			//bool AlphaShape(double alpha_param, Triangle2D<Node<3> >& pgeom)
		{
			KRATOS_TRY
			
			
			double x0 = pgeom[0].X();
			double x1 = pgeom[1].X();
			double x2 = pgeom[2].X();
			
			double y0 = pgeom[0].Y();
			double y1 = pgeom[1].Y();
			double y2 = pgeom[2].Y();
						
			mJ(0,0)=2.0*(x1-x0);	mJ(0,1)=2.0*(y1-y0);
			mJ(1,0)=2.0*(x2-x0);	mJ(1,1)=2.0*(y2-y0);
			
			
			double detJ = mJ(0,0)*mJ(1,1)-mJ(0,1)*mJ(1,0);
						
			mJinv(0,0) =  mJ(1,1); mJinv(0,1) = -mJ(0,1);
			mJinv(1,0) = -mJ(1,0); mJinv(1,1) =  mJ(0,0);
		
			bounded_matrix<double,2,2> check;
		
			
			if(detJ < 1e-12) 
			{
				//std::cout << "detJ = " << detJ << std::endl;
				////mark as boundary
				pgeom[0].GetSolutionStepValue(IS_BOUNDARY) = 1;
				pgeom[1].GetSolutionStepValue(IS_BOUNDARY) = 1;
				pgeom[2].GetSolutionStepValue(IS_BOUNDARY) = 1;
				return false;
			}
			
			else
			{

				double x0_2 = x0*x0 + y0*y0;
				double x1_2 = x1*x1 + y1*y1; 
				double x2_2 = x2*x2 + y2*y2; 

				//finalizing the calculation of the inverted matrix
				//std::cout<<"MATR INV"<<MatrixInverse(mJ)<<std::endl;
				mJinv /= detJ;
				//calculating the RHS
				mRhs[0] = (x1_2 - x0_2);
				mRhs[1] = (x2_2 - x0_2);

				//calculate position of the center
				noalias(mC) = prod(mJinv,mRhs);

				double radius = sqrt(pow(mC[0]-x0,2)+pow(mC[1]-y0,2));

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
		TriGenModeler& operator=(TriGenModeler const& rOther);


		///@}    

	}; // Class TriGenModeler 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		TriGenModeler& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const TriGenModeler& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_TRIGEN_MODELER_H_INCLUDED  defined 


