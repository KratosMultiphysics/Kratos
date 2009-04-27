//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2008-02-25 19:05:54 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_TRIGEN_CDT_H_INCLUDED )
#define  KRATOS_TRIGEN_CDT_H_INCLUDED

 

// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>

#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h" 
#endif

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
	class TriGenCDT  
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of TriGenCDT
		KRATOS_CLASS_POINTER_DEFINITION(TriGenCDT);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		TriGenCDT() {} //

		/// Destructor.
		virtual ~TriGenCDT(){}


		///@}
		///@name Operators 
		///@{


		///@}
		///@name Operations
		///@{


		//*******************************************************************************************
		//*******************************************************************************************
		void GenerateCDT(
			ModelPart& ThisModelPart , 
			Element const& rReferenceElement
			)
		{

			KRATOS_TRY


			struct triangulateio in;
			struct triangulateio mid;
			struct triangulateio out;
			struct triangulateio vorout;
			
			clean_triangulateio(in);
			clean_triangulateio(mid);
			clean_triangulateio(out);
			clean_triangulateio(vorout);
			
			//*********************************************************************
			//input mesh			
			in.numberofpoints = ThisModelPart.Nodes().size();
			in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
			//writing the points coordinates in a vector
			ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();
			for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
			{
				int base = i*2;
				//int base = ((nodes_begin + i)->Id()   -  1 ) * 2;

				//from now on it is consecutive
                                (nodes_begin + i)->SetId(i+1);
//				(nodes_begin + i)->Id() = i+1;

				in.pointlist[base] = (nodes_begin + i)->X();
				in.pointlist[base+1] = (nodes_begin + i)->Y();
			}			
			
			in.numberoftriangles = ThisModelPart.Elements().size();
			in.trianglelist = (int*) malloc(in.numberoftriangles * 3 * sizeof(int));
			//copying the elements in the input file
			ModelPart::ElementsContainerType::iterator elem_begin = ThisModelPart.ElementsBegin();
			for(unsigned int el = 0; el<ThisModelPart.Elements().size(); el++)
			{
				int base = el * 3;
				
				Geometry<Node<3> >& geom = (elem_begin+el)->GetGeometry();
				

				in.trianglelist[base] = geom[0].Id();
				in.trianglelist[base+1] = geom[1].Id();
				in.trianglelist[base+2] = geom[2].Id();
			}
			
			//clearing elements
			ThisModelPart.Elements().clear();
			
			KRATOS_WATCH(in.numberofsegments);
			KRATOS_WATCH(in.numberofpoints);
			KRATOS_WATCH(in.numberoftriangles);
			KRATOS_WATCH(in.numberofholes);  
			
			//read and regenerate the existing mesh ... to generate the boundaries
			char options[] = "rcEBQ";
			triangulate(options, &in, &mid, &vorout);
			
			//free the memory used in the first step
			free( in.pointlist );
			free( in.trianglelist );
			
			//uses the boundary list generated at the previous step to generate the "skin"
			mid.numberoftriangles = 0;
			free(mid.trianglelist);
			char regeneration_option[] = "pBQY";
			triangulate(regeneration_option, &mid, &out, &vorout);
								
			KRATOS_WATCH(out.numberofsegments);
			KRATOS_WATCH(out.numberofpoints);
			KRATOS_WATCH(out.numberoftriangles);
			KRATOS_WATCH(out.numberofholes); 

			//free the memory used in the intermediate step
			free( mid.pointlist );
			free( mid.segmentlist );
			
			//*******************************************************************
			//properties to be used in the generation
			Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);

			//generate Kratos triangle
			int el_number = out.numberoftriangles;
			
			//note that the node list can not be changed starting from here
			ModelPart::NodesContainerType& ModelNodes = ThisModelPart.Nodes();
			
			//generate kratos elements (conditions are not touched)
			for(int el = 0; el< el_number; el++)
			{
				int id = el + 1;
				int base = el * 3;

				Geometry<Node<3> > temp;
				temp.push_back( *((ModelNodes).find( out.trianglelist[base]).base()	) );
				temp.push_back( *((ModelNodes).find( out.trianglelist[base+1]).base()	) );
				temp.push_back( *((ModelNodes).find( out.trianglelist[base+2]).base()	) );
				
				Element::Pointer p_element = rReferenceElement.Create(id, temp, properties);
				ThisModelPart.Elements().push_back(p_element);
			}
			

/*			free( in.pointmarkerlist);
			free( in.pointmarkerlist );
			free( in.regionlist );*/
	


			free( out.pointattributelist );
			free( out.trianglelist );
// 			free( out.triangleattributelist );
//			free( out.neighborlist );


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
		void clean_triangulateio( triangulateio& tr ){

			tr.pointlist                  = (REAL*) NULL;
			tr.pointattributelist         = (REAL*) NULL;
			tr.pointmarkerlist            = (int*) NULL;
			tr.numberofpoints             = 0;
			tr.numberofpointattributes    = 0;
			tr.trianglelist               = (int*) NULL;
			tr.triangleattributelist      = (REAL*) NULL;
			tr.trianglearealist           = (REAL*) NULL;
			tr.neighborlist               = (int*) NULL;
			tr.numberoftriangles          = 0;
			tr.numberofcorners            = 3;
			tr.numberoftriangleattributes = 0;
			tr.segmentlist                = (int*) NULL;
			tr.segmentmarkerlist          = (int*) NULL;
			tr.numberofsegments           = 0;
			tr.holelist                   = (REAL*) NULL;
			tr.numberofholes              = 0;
			tr.regionlist                 = (REAL*) NULL;
			tr.numberofregions            = 0;
			tr.edgelist                   = (int*) NULL;
			tr.edgemarkerlist             = (int*) NULL;
			tr.normlist                   = (REAL*) NULL;
			tr.numberofedges              = 0;

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
		TriGenCDT& operator=(TriGenCDT const& rOther);


		///@}    

	}; // Class TriGenCDT 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		TriGenCDT& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const TriGenCDT& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_TRIGEN_CDT_H_INCLUDED  defined 


