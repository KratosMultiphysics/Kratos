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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_MODELER_H_INCLUDED )
#define  KRATOS_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"


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
  class Modeler
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Modeler
      KRATOS_CLASS_POINTER_DEFINITION(Modeler);

	  typedef std::size_t SizeType;
	  typedef std::size_t IndexType;
	  	  
	  //typedef ModelPart::GeometricalDataContainerType GeometricalDataContainerType;

	  //typedef GeometricalDataContainerType::GeometricalDataType GeometricalDataType;

	  //typedef ModelPart::GeometryType GeometryType;

	  //typedef ModelPart::GeometriesContainerType GeometriesContainerType;


      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Modeler(){}

      /// Destructor.
	  virtual ~Modeler(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{

	  
	  /*	  void CollapsePoints(ModelPart& rThisModelPart, double Tolerance)
	  {
		  double distance;
		  Point3D<Node<3> >::Pointer p_founded_point;

		  GeometriesContainerType& geometries = rThisModelPart.Geometries();
		  GeometriesContainerType::PointsArrayType& r_points_array = geometries.Points();

		  if(geometries.NumberOfGeometries() == 0)
			  return;

		  for(GeometriesContainerType::GeometryIterator i_geometry = geometries.GeometriesBegin() ; 
			  i_geometry != geometries.GeometriesEnd() ; i_geometry++)
		  {
		  // At this moment a brute force search is used
			  for(GeometryType::iterator i_point = i_geometry->begin() ; i_point != i_geometry->end() ; i_point++)
			  {
				  bool founded = false;
				  for(GeometriesContainerType::PointIterator j_point = r_points_array.begin() ; 
					  j_point != r_points_array.end() ; j_point++)
				  {
					distance = j_point->Distance(*i_point); 
					if(distance < Tolerance)
					{
						founded = true;
						p_founded_point = *(j_point.base());
						break;
					}
				  }
				  if(founded)
				  {
					  *(i_point.base()) = p_founded_point->pGetPoint(0);
				  }
				  else
				  {
					  r_points_array.push_back(Point3D<Node<3> >(*(i_point.base())));
				  }
			  }
		  }

	  }
	  */


//	  virtual void GenerateMesh(ModelPart& ThisModelPart, Element const& rReferenceElement)
//	  {
//		  KRATOS_ERROR(std::logic_error, "This modeler CAN NOT be used for mesh generation.", "");
//	  }
//

	  virtual void GenerateMesh(ModelPart& ThisModelPart, Element const& rReferenceElement, Condition const& rReferenceBoundaryCondition)
	  {
		  KRATOS_ERROR(std::logic_error, "This modeler CAN NOT be used for mesh generation.", "");
	  }

	  virtual void GenerateNodes(ModelPart& ThisModelPart)
	  {
		  KRATOS_ERROR(std::logic_error, "This modeler CAN NOT be used for node generation.", "");
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
      virtual std::string Info() const
	  {
		  return "Modeler";
	  }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	  {
		  rOStream << Info();
	  }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
	  {
	  }
      
            
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
      Modeler& operator=(Modeler const& rOther);

      /// Copy constructor.
      Modeler(Modeler const& rOther);

        
      ///@}    
        
    }; // Class Modeler 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream, 
				    Modeler& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Modeler& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_MODELER_H_INCLUDED  defined 


