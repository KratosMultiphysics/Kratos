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
//   Last Modified by:    $Author: janosch $
//   Date:                $Date: 2008-04-29 11:37:37 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_POINT_H_INCLUDED )
#define  KRATOS_POINT_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <sstream>
#include <cstddef>


// Project includes
#include "includes/define.h"
#include "containers/array_1d.h"


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
  
  /// Point class. 
  /** Point class. Stores coordinates of a point and have some basic
      operations defined.  Point class has two template parameter:

      - TDimension which define the dimension of the point.  TDataType
      - which specifies the point coordinate value type. This type by
        default is double.

	@see Geometry
	@see Node
	@see IntegrationPoint
  */
	template<std::size_t TDimension, class TDataType = double>
    class Point : public array_1d<TDataType, TDimension>
    {
      public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of Point
      KRATOS_CLASS_POINTER_DEFINITION(Point);
      
      typedef array_1d<TDataType, TDimension> BaseType;
  
      typedef Point<TDimension, TDataType> Type;

      typedef BaseType CoordinatesArrayType;

	  typedef typename std::size_t SizeType;

	  typedef typename std::size_t IndexType;

      ///@}
      ///@name Constants 
      ///@{ 

      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      Point() : BaseType(TDimension)
      {
	KRATOS_TRY_LEVEL_4
	SetAllCoordinates();
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /// 1d constructor.
      Point(TDataType const& NewX) : BaseType(TDimension)
      {
	KRATOS_TRY_LEVEL_4
	SetAllCoordinates();
	this->operator()(0) = NewX;
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /// 2d constructor.
      Point(TDataType const& NewX, TDataType const& NewY) : BaseType(TDimension)
      {
	KRATOS_TRY_LEVEL_4
	SetAllCoordinates();
	this->operator()(0) = NewX;
	this->operator()(1) = NewY;
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /// 3d constructor.
      Point(TDataType const& NewX, TDataType const& NewY, TDataType const& NewZ) : BaseType(TDimension)
      {
	KRATOS_TRY_LEVEL_4
	SetAllCoordinates();
	this->operator()(0) = NewX;
	this->operator()(1) = NewY;
	this->operator()(2) = NewZ;
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /** Copy constructor. Initialize this point with the coordinates
	  of given point.*/
      Point(Point const& rOtherPoint)
	: BaseType(rOtherPoint){}
      
      /** Copy constructor from a point with different dimension. Initialize this point with the coordinates
	  of given point.*/
      template<SizeType TOtherDimension>
      Point(Point<TOtherDimension> const& rOtherPoint) : BaseType(TDimension)
      {
	KRATOS_TRY_LEVEL_4
	IndexType size = (TDimension < TOtherDimension) ? TDimension : TOtherDimension;
	IndexType i;

	for(i = 0 ; i < size ; i++)
	  this->operator[](i)=rOtherPoint[i];

	for(i = size ; i < TDimension ; i++)
	  this->operator[](i)= TDataType();

	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /** Constructor using coordinates stored in given array. Initialize
	  this point with the coordinates in the array. */
      Point(CoordinatesArrayType const& rOtherCoordinates) 
	: BaseType(rOtherCoordinates){}
      
        /** Constructor using coordinates stored in given array. Initialize
        this point with the coordinates in the array. */
      template<class TVectorType>
      Point(vector_expression<TVectorType> const&  rOtherCoordinates)
	: BaseType(rOtherCoordinates){}
        
        /** Constructor using coordinates stored in given std::vector. Initialize
        this point with the coordinates in the array. */
      Point(std::vector<TDataType> const&  rOtherCoordinates) : BaseType(TDimension)
      {
	KRATOS_TRY_LEVEL_4
	SizeType size = rOtherCoordinates.size();
	size = (TDimension < size) ? TDimension : size;
	for(IndexType i = 0 ; i < size ; i++)
	  this->operator[](i)=rOtherCoordinates[i];
	KRATOS_CATCH_LEVEL_4(*this)
      }
        
      /// Destructor.
      virtual ~Point(){}
      

      ///@}
      ///@name Operators 
      ///@{
      
      /// Assignment operator.
      Point& operator=(const Point& rOther)
      {
	CoordinatesArrayType::operator=(rOther);
	return *this;
      }
      
      bool operator==(const Point& rOther)
      {
	return std::equal(this->begin(), this->end(), rOther.begin());
      }
      
      /// Assignment operator.
      template<SizeType TOtherDimension>
      Point& operator=(const Point<TOtherDimension>& rOther)
      {
	KRATOS_TRY_LEVEL_4
	IndexType size = (TDimension < TOtherDimension) ? TDimension : TOtherDimension;
	IndexType i;

	for(i = 0 ; i < size ; i++)
	  this->operator[](i)=rOther[i];

	for(i = size ; i < TDimension ; i++)
	  this->operator[](i)= TDataType();

	return *this;

	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      static IndexType Dimension() {return TDimension;}

      /** Returns X coordinate */
      TDataType X() const
      {
	KRATOS_TRY_LEVEL_4
	return this->operator[](0);
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /** Returns Y coordinate */
      TDataType Y() const
      {
	KRATOS_TRY_LEVEL_4
	return this->operator[](1);
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /** Returns Z coordinate */
      TDataType Z() const
      {
	KRATOS_TRY_LEVEL_4
	return this->operator[](2);
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      TDataType& X()
      {
	KRATOS_TRY_LEVEL_4
	return this->operator[](0);
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /** Returns Y coordinate */
      TDataType& Y()
      {
	KRATOS_TRY_LEVEL_4
	return this->operator[](1);
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /** Returns Z coordinate */
      TDataType& Z()
      {
	KRATOS_TRY_LEVEL_4
	return this->operator[](2);
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      
      
      /** This is an access method to point's coordinate by indices. For example this
	  function return x, y and z coordinate whith 1, 2 and 3 as input
	  respectively. 
	  */
      TDataType  Coordinate(IndexType CoordinateIndex) const
      {
	KRATOS_TRY_LEVEL_4
	return this->operator[](CoordinateIndex - 1);
	KRATOS_CATCH_LEVEL_4(*this)
      }
      
      /** This is an access method to get a reference to point's coordinate by
	  indices. For example this function return references to x, y and z coordinate whith 1, 2
	  and 3 as input respectively. 
	  */
      TDataType& Coordinate(IndexType CoordinateIndex)
      {
	KRATOS_TRY_LEVEL_4
	return this->operator[](CoordinateIndex - 1);
	KRATOS_CATCH_LEVEL_4(*this)
      }

      CoordinatesArrayType const& Coordinates() const
      {
	return *this;
      }
      
      CoordinatesArrayType& Coordinates()
      {
	return *this;
      }
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	  std::stringstream buffer;
	  buffer << TDimension << " dimensional point";
	  return buffer.str();
      }
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	  rOStream << TDimension << " dimensional point";
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
      {
	if(!TDimension)
	  return;

	  rOStream << "("  << this->operator[](0);
	
	  for(IndexType i = 1 ; i < TDimension  ; i++)
	    rOStream << " , " << this->operator[](i);
	  rOStream << ")";
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

      void SetAllCoordinates(TDataType const& Value = TDataType())
      {
	KRATOS_TRY_LEVEL_4
	for(IndexType i = 0 ; i < TDimension ; i++)
	this->operator()(i) = Value;
	KRATOS_CATCH_LEVEL_4(*this)
      }
        
        
//       friend class boost::serialization::access;
      
//       template<class TArchive>
// 	  void serialize(TArchive & ThisArchive, const unsigned int ThisVersion)
// 	  {
//  	      ThisArchive & boost::serialization::base_object<array_1d<TDataType, TDimension> >(*this); 
// 	  }
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
        
      ///@}    
        
    }; // Class Point 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<std::size_t TDimension, class TDataType>
  inline std::istream& operator >> (std::istream& rIStream, 
				    Point<TDimension, TDataType>& rThis);

  /// output stream function
  template<std::size_t TDimension, class TDataType>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const Point<TDimension, TDataType>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_POINT_H_INCLUDED  defined 


